plotEPI <- function(cld_data, n=100, order_by='pepi') {
  pepi_z <- (cld_data[['pepi']] - mean(cld_data[['pepi']], na.rm=TRUE)) /
    sd(cld_data[['pepi']], na.rm=TRUE)
  epi_z <- (cld_data[['epi']] - mean(cld_data[['epi']], na.rm=TRUE)) /
    sd(cld_data[['epi']], na.rm=TRUE)
  log_ed <- log(cld_data[['ed']])
  ed_z <- (mean(log_ed, na.rm=TRUE) - log_ed) /
    sd(log_ed, na.rm=TRUE)
  ordrd_top <- order(cld_data[[order_by]])[1:n]
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

#http://docs.ggplot2.org/current/geom_smooth.html
binomial_smooth <- function(...) {
  geom_smooth(method = "glm", method.args = list(family = "binomial"), ...)
}

lookUp <- function(ident, slt) {
  pull <- cld_data[[slt]] %in% ident
  cld_data[pull, ]
}

getWikiByQuantile <- function(indx) {
  qts <- quantile(cld_data[[indx]], na.rm=TRUE)
  res <- rep(NA, length(qts) - 1)
  for(i in 2:length(qts)) {
    pull <- cld_data[[indx]] > qts[i-1]
    pull <- which(pull & cld_data[[indx]] < qts[i])
    res[i-1] <- sum(cld_data$wiki[pull], na.rm=TRUE)
  }
  names(res) <- c('0-25%', '25-50%', '50-75%', '75-100%')
  res
}