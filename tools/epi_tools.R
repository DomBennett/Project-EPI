
normalise <- function(xs) {
  (xs - min(xs, na.rm=TRUE)) /
    ( max(xs, na.rm=TRUE) - min(xs, na.rm=TRUE) )
}


EPIChecker <- function (metrics, cut) {
  hist(metrics$time)
  hist(metrics$change)
  hist(metrics$success)
  plot(time ~ epi, data=metrics)
  abline (lm (time ~ epi, data=metrics), col = "red")
  plot(change ~ epi, data=metrics)
  abline (lm (change ~ epi, data=metrics), col = "red")
  plot(success ~ epi, data=metrics)
  abline (lm (success ~ epi, data=metrics), col = "red")
  plot(epi_nc ~ epi, data=metrics)
  abline (lm (epi_nc ~ epi, data=metrics), col = "red")
  hist (metrics$epi, xlab = 'EPI', ylab = NULL, main = NULL,
        col = 'cornflowerblue')
  cutoff <- quantile (metrics$epi, probs = cut, na.rm=TRUE)
  abline (v = cutoff, col = "red", lwd = 2)
  text (labels = '<-- Living fossils', x=cutoff, y=nrow(metrics)/100, cex = 0.8)
  hist (metrics$epi_nc, xlab = 'EPI (No change)', ylab = NULL, main = NULL,
        col = 'cornflowerblue')
  cutoff <- quantile (metrics$epi_nc, probs = cut, na.rm=TRUE)
  abline (v = cutoff, col = "red", lwd = 2)
  text (labels = '<-- Living fossils', x=cutoff, y=nrow(metrics)/100, cex = 0.8)
}

plotEPI <- function(epi, n=20) {
  ordrd_top <- order(epi[['pepi']])[1:n]
  val <- c(epi[['pepi']][ordrd_top],
           epi[['epi']][ordrd_top],
           epi[['ed']][ordrd_top])
  type <- rep(c('pEPI', 'EPI', 'ED'), each=length(ordrd_top))
  nms <- epi[['nm']][ordrd_top]
  nm <- rep(nms, 3)
  nm <- factor(nm, levels=nms[length(nms):1])
  p_data <- data.frame(nm, val, type)
  p <- ggplot(data=p_data, aes(x=nm, y=val, colour=type)) +
    geom_point() +
    coord_flip() +
    scale_y_reverse() +
    xlab("") + ylab("")
  p
}
