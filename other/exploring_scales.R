# Change and success are both contrasted measures
#  but they may still be on different scales

library(ggplot2)

combHisto <- function(cntrst_c, cntrst_s) {
  cntrst_c <- data.frame(val=cntrst_c)
  cntrst_c$type <- 'change'
  cntrst_s <- data.frame(val=cntrst_s)
  cntrst_s$type <- 'success'
  p_data <- rbind(cntrst_s, cntrst_c)
  p_data$type <- as.factor(p_data$type)
  p <- ggplot(p_data, aes(val, fill=type, colour=type)) + geom_density(alpha=0.75)
  p
}

# normal
changes <- rnorm(n=100, mean=1000, sd=100)
success <- rnorm(n=100, mean=10000, sd=1000)
cntrst_c <- changes[seq(1, 99, 2)]/changes[seq(2, 100, 2)]
cntrst_s <- success[seq(1, 99, 2)]/success[seq(2, 100, 2)]
combHisto(cntrst_c, cntrst_s)
# this is no problem if both distributions are normal

# varying distributions
changes <- rnorm(n=100, mean=1000, sd=100)
success <- rexp(100, rate = 0.01)  # more realistic
cntrst_c <- changes[seq(1, 99, 2)]/changes[seq(2, 100, 2)]
cntrst_s <- success[seq(1, 99, 2)]/success[seq(2, 100, 2)]
combHisto(cntrst_c, cntrst_s)
# problems arise if the distributions differ

# distribution matters, scale doesn't


chars1a <- sample(1:3, 100, replace = TRUE)
chars1b <- sample(1:3, 100, replace = TRUE)
chars1 <- chars1a/chars1b
chars2a <- sample(1:10, 100, replace = TRUE)
chars2b <- sample(1:10, 100, replace = TRUE)
chars2 <- chars2a/chars2b
combHisto(chars1, chars2)

mean(chars2)
mean(chars1)
wts <- rep(1/c(3, 10), each=100)
weighted.mean(c(chars1, chars2), w=wts)


# TODO:
# -- update change step to calculate rsqs
# -- create function set to look up rsqs for cntrst change
# -- recaluclate a weighted mean change score


char_nms <- sample(colnames(chars), 3)

lookUpRsqs <- function(char_nms) {
  for(i in 1:ncol(cmbs)) {
    
  }
  
}

nchars <- c(3, 3, 5, 5, 10, 10)
rsqs <- c(.1, .9, .1, .9, .1, .9)
wts <- (1/nchars)/rsqs



for(i in columns) {
  tmp_columns <- columns[-i]
  for(j in 2:ncol) {
    tmp_columns
  }
}