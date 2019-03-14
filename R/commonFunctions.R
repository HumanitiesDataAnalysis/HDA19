readfile <- function(filename) {
  read_lines(filename) %>%
    tibble(text = .) %>%
    filter(text != "") %>%
    # Add a counter for the paragraph number
    mutate(paragraph = 1:n()) %>%
    unnest_tokens(word, text) %>%
    mutate(filename = filename %>%
      str_replace(".*/", "") %>%
      str_replace(".txt", ""))
}

read_all_SOTUs <- function() {
  raw <- list.files("../data/SOTUS/", full.names = TRUE) %>% 
    map_dfr(readfile) %>% 
    mutate(year = as.numeric(filename))
  
  meta <- read_csv("../data/SOTUmetadata.csv") %>% distinct(year, president, party, sotu_type)
  inner_join(raw, meta)
}


summarize_pmi <- function(data, token, count = rep(1, n())) {
  token <- enquo(token)
  count <- enquo(count)
  groupings <- groups(data)
  data %>%
    group_by(!!!groupings, !!token) %>%
    summarize(count = sum(!!count)) %>%
    ungroup() %>%
    mutate(total_words = sum(count)) %>%
    group_by(!!token) %>%
    mutate(word_share = sum(count) / total_words) %>%
    group_by(!!!groupings) %>%
    mutate(p_share = count / sum(count)) %>%
    mutate(pmi = log(p_share / word_share)) %>%
    select(!!!groupings, !!token, pmi)
}

#' Add TF-IDF summary based on groupings.
#'
#' @description The TF-IDF function in tidytext requires an explicit 'doc'
#' parameter; this applies it on the existing dataset groups.
#'
#' @param frame A grouped data from
#' @param word The unquoted variable name storing terms for frequency
#' @param count The unquoted variable storing the count
#'
#' @return A data_frame with a column tfidf added.
#' @export
summarize_tf_idf <- function(data, word, count = rep(1, n())) {
  token <- enquo(word)
  count <- enquo(count)
  groupings <- groups(data)
  n_docs <- data %>% n_groups()
  data <- data %>%
    group_by(!!token, add = TRUE) %>%
    summarize(.count := sum(!!count))
  data %>%
    distinct(!!!groupings, !!token) %>%
    group_by(!!token) %>%
    summarize(idf = log(n_docs / n())) %>%
    inner_join(data) %>%
    group_by(!!!groupings) %>%
    mutate(doc_total = sum(.count)) %>%
    group_by(!!token, add = TRUE) %>%
    summarize(tf = sum(.count) / doc_total[1], idf = idf[1], tfidf = tf * idf) %>%
    ungroup()
}


#' Summarize the log-likelihood ratio across a grouped data frame
#'
#' @param data A data frame
#' @param token The column indicating a token.
#' @param count The column indicating wordcount data.
#'
#' @return A dataframe with the supplied grouping and a log-likelihood for each token in that grouping.
#' Strongly positive numbers are over-represented; strongly negative numbers are under-represented.
#' Also included are unadjusted p-values assuming a single degree of freedom. It would be wise to correct this
#' for the number of words using Sedak's correction or some other method; in practice, the size of this correction
#' is highly responsive.
#' @export
#'
summarize_llr <- function(data, token, count = rep(1, n())) {
  token <- enquo(token)
  count <- enquo(count)
  groupings <- groups(data)

  data %>%
    group_by(!!token, add = TRUE) %>%
    # Numeric to fix some integer overflow problems in very large sets.
    summarize(count = sum(!!count) %>% as.numeric()) %>%
    group_by(!!token) %>%
    # Total number of words in the set.
    mutate(grandtot = sum(count)) %>%
    group_by(!!!groupings) %>%
    # Three levels of counts: for this group, for all the other groups together *except* this one (count.y),
    # and for all groups together (grandtot)
    mutate(count.x = count, count.y = grandtot - count) %>%
    # We need not just the token-level counts, but the counts summed across all words.
    mutate(c = sum(count.x), d = sum(count.y), totalWords = c + d) %>%
    # A phony correction to allow logarithms; if a word does not appear in the 
    # 
    mutate(count.y = ifelse(count.y == 0, 1, count.y)) %>%
    # Expected counts based on uniform rates.
    mutate(exp.x = c * grandtot / totalWords, exp.y = d * grandtot / totalWords) %>%
    # The Dunning log-likelihood based on mutual information.
    mutate(score = 2 * (
      (count * log(
        count / exp.x
      )) + count.y * log(count.y / exp.y))) %>%
    mutate(p = pchisq(abs(score), df = 1)) %>%
    # Use negative signage to indicate *under-represented* words.
    mutate(score = ifelse((count.y - exp.y) > 0, -score, score)) %>%
    # Throw away some of the on-the-way calculations.
    select(!!!groupings, !!token, dunning_llr = score, dunning_p = p) %>%
    ungroup()
}

add_chunks <- function(data, count = NULL, chunk_length = 2000) {
  count <- enquo(count)
  group_names <- group_vars(data)
  groupings <- groups(data)
  # If there are no counts, then the previous ordering is preserved exactly 
  if (rlang::quo_is_null(count)) {
    data %>% 
      mutate(.count = 1) %>%
      mutate(.tally = cumsum(.count), chunk = 1 + .tally %/% chunk_length) %>%
      select(-.tally, -.count)
  } else {
    d = data %>%
      summarize(.count = sum(!!count)) %>%
      mutate(.tally = cumsum(.count), chunk = 1 + .tally %/% chunk_length) %>%
      select(-.tally, -.count)
    data %>% inner_join(d, by = group_names)
  }
}

# Utility function: change 'summarize_tf' to 'add_tf', where the original data is maintained.
lift_to_merger = function(f) {
  function(data, ...) {
    summary_version = f(data, ...)
    # Early return when no merge is needed.
    cat(group_vars(data))
    if (length(group_vars(data) == 0)) {
      summary_version
    } else {
      data %>% left_join(summary_version)
    }
  }
}

summarize_mann_whitney <- function(data, token, ..., count = rep(1, n())) {
  token <- enquo(token)
  count <- enquo(count)
  
  # Mann Whitney needs initial groupings as well as a higher-level
  # document representation at the end. For example, you might 
  # group by federalist paper, but then compare authors of 
  # federalists papers as the final document.

  doc = enquos(...)
  
  groupings <- groups(data)
  
  d = data %>% group_by(!!token, add = TRUE) %>% summarize(.count = sum(!!count))
  
  long_long_long = d %>%
    group_by(!!token) %>% ungroup %>% complete(word, nesting(title, author), fill = list(.count = 0)) %>%
    mutate(rank = rank(.count/sum(.count))) %>%
    group_by(!!token, !!!doc) %>%
    summarize(U = sum(rank), n = n()) %>%
    group_by(!!token) %>% mutate(n = sum(n), UU = sum(U))
  long_long_long %>% mutate(p = psignrank(U, n)) %>% arrange(-p)
    group_by(word, !!!doc) %>% mutate(rank = rank(share))
    filter(word == "upon") -> e
  
  rows = 1:nrow(e)
  share = e$share
  map_dbl(rows, ~wilcox.test(unlist(share[.x]), unlist(share[-.x]))$p.value * 100)
  
}
