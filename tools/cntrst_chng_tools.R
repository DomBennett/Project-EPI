lookUpRsqs <- function(nms) {
  # Look up mean R2s for characters from rsq_obj
  .srch1 <- function(i) {
    res <- NA
    tst <- rsq_obj[['cmbs']][,i] %in% nms
    if(sum(tst) == 2) {
      res <- i
    }
    res
  }
  .getrsq <- function(i, nm) {
    res <- NA
    if(nm %in% rsq_obj[['cmbs']][,i]) {
      res <- rsq_obj[['rsq']][i]
    }
    res
  }
  .srch2 <- function(nm) {
    res <- sapply(srch, .getrsq, nm=nm)
    res <- res[!is.na(res)]
    mean(res)
  }
  srch <- sapply(1:ncol(rsq_obj[['cmbs']]), .srch1)
  srch <- srch[!is.na(srch)]
  sapply(nms, .srch2)
}