
getTTOL <- function(sp1, sp2) {
  # search timetree via html
  # TODO: avoid line numbers, control for failure
  getVal <- function(strng) {
    stop <- regexpr(" Mya", strng) - 1
    as.numeric(substr(strng, start=26, stop=stop))
  }
  url <- "http://www.timetree.org/search/pairwise/"
  sp1 <- gsub("\\s+", "%20", sp1)
  sp2 <- gsub("\\s+", "%20", sp2)
  qry <- paste0(url, sp1, '/', sp2)
  res <- readLines(qry)
  unlink(qry)
  mean_ttol <- getVal(res[[195]])
  median_ttol <- getVal(res[[199]])
  data.frame(mean_ttol, median_ttol)
}