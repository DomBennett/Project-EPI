
genDataframe <- function(cnddts, node_obj) {
  # generate clade data frame for clades of interest
  # ensure sisters are accounted for in order to test autocorrelation
  cnddts <- cnddts[!duplicated(cnddts)]
  cnddts <- unique(unlist(sapply(cnddts,
                                 function(x) c(x, node_obj[[x]][['sstr']]))))
  nms <- tmsplts <- cntrst_ns <- cntrst_chngs <-
    chngs <- eds <- ns <- rep(NA, length(cnddts))
  # loop through candidate IDs and extract info
  for(i in 1:length(cnddts)) {
    tmsplt <- node_obj[[cnddts[i]]][["tmsplt"]]
    # if no tmsplt, use sstrs
    if(is.null(tmsplt)) {
      sstrs <- node_obj[[cnddts[i]]][['sstr']]
      if(length(sstrs) > 0) {
        sstr_tmsplts <- sapply(sstrs, function(x) {
          res <- NA
          if(!is.null(node_obj[[x]][['tmsplt']])) {
            res <- node_obj[[x]][['tmsplt']]
          }
          res})
        tmsplt <- mean(sstr_tmsplts, na.rm=TRUE)
      }
    }
    cntrst_n <- node_obj[[cnddts[i]]][["cntrst_n"]]
    if(!is.null(tmsplt) & !is.null(cntrst_n)) {
      tmsplts[i] <- tmsplt
      cntrst_ns[i] <- cntrst_n
      ns[i] <- length(node_obj[[cnddts[i]]][['kids']])
      nms[i] <- node_obj[[cnddts[i]]][["nm"]][["scientific name"]]
    }
    if(!is.null(node_obj[[cnddts[i]]][['ed']])) {
      eds[i] <- node_obj[[cnddts[i]]][['ed']]
    }
    if(!is.null(node_obj[[cnddts[i]]][["cntrst_chng"]])) {
      cntrst_chngs[i] <- node_obj[[cnddts[i]]][["cntrst_chng"]]
      chngs[i] <- node_obj[[cnddts[i]]][["chng"]]
    }
  }
  bool <- !is.na(cntrst_ns)
  cntrst_ns <- cntrst_ns[bool]
  tmsplts <- tmsplts[bool]
  cnddts <- cnddts[bool]
  nms <- nms[bool]
  cntrst_chngs <- cntrst_chngs[bool]
  chngs <- chngs[bool]
  eds <- eds[bool]
  ns <- ns[bool]
  cld_data <- data.frame(nm=nms, txid=cnddts, tmsplt=tmsplts, cntrst_n=cntrst_ns,
                         chng=chngs, ed=eds, cntrst_chng=cntrst_chngs, n=ns,
                         stringsAsFactors=FALSE)
  cld_data
}


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

plotEPI <- function(epi, n=100) {
  pepi_z <- (cld_data[['pepi']] - mean(cld_data[['pepi']], na.rm=TRUE)) /
    sd(cld_data[['pepi']], na.rm=TRUE)
  epi_z <- (cld_data[['epi']] - mean(cld_data[['epi']], na.rm=TRUE)) /
    sd(cld_data[['epi']], na.rm=TRUE)
  log_ed <- log(cld_data[['ed']])
  ed_z <- (mean(log_ed, na.rm=TRUE) - log_ed) /
    sd(log_ed, na.rm=TRUE)
  ordrd_top <- order(cld_data[['pepi']])[1:n]
  val <- c(pepi_z[ordrd_top],
           epi_z[ordrd_top],
           ed_z[ordrd_top])
  type <- rep(c('pEPI', 'EPI', 'ED'), each=length(ordrd_top))
  nms <- cld_data[['nm']][ordrd_top]
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
