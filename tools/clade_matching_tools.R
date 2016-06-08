getMtchScrs <- function(kids_vctr, kids_lst) {
  # Calculate match scores for a vector of kids against
  # a list of vectors
  # Match scores = 2, if both vectors share the exact names
  .mtch <- function(qry) {
    mtch_scr <- (sum(kids_vctr %in% qry)/length(kids_vctr)) +
      (sum(qry %in% kids_vctr)/length(qry))
    mtch_scr
  }
  sapply(kids_lst, .mtch)
}