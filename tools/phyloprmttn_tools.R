meanPD <- function(nms, trees, ntrees) {
  dst_mtrx <- calcDstMtrx(trees[[1]], ids=nms)
  pd <- mean(dst_mtrx)
  if(ntrees > 1) {
    for(i in 2:ntrees) {
      dst_mtrx <- calcDstMtrx(trees[[i]], ids=nms)
      tmp_pd <- mean(dst_mtrx)
      pd <- (tmp_pd + pd)/2
    }
  }
  pd
}