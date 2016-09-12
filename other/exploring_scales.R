# Change and success are both contrasted measures
#  but they may still be on different scales

library(ggplot2)

# normal
changes <- rnorm(n=100, mean=1000, sd=100)
success <- rnorm(n=100, mean=10000, sd=1000)
cntrst_c <- changes[seq(1, 99, 2)]/changes[seq(2, 100, 2)]
cntrst_s <- success[seq(1, 99, 2)]/success[seq(2, 100, 2)]

p_data <- data.frame(val=c(cntrst_c, cntrst_s),
                     type=rep(c('change', 'success'), each=100))
p_data <- data.frame(val=c(cntrst_c),
                     type=factor(rep(c('change'), each=100)))
p <- ggplot(p_data, aes(val, fill=factor(type))) + geom_histogram()
p
# this is no problem if both distributions are normal

# varying distributions
changes <- rnorm(n=100, mean=1000, sd=100)
success <- rexp(100, rate = 0.01)  # more realistic
cntrst_c <- changes[seq(1, 99, 2)]/changes[seq(2, 100, 2)]
cntrst_s <- success[seq(1, 99, 2)]/success[seq(2, 100, 2)]

p_data <- data.frame(val=c(cntrst_c, cntrst_s),
                     type=rep(c('change', 'success'), each=100))
p <- ggplot(p_data, aes(val, fill=factor(type))) + geom_histogram()
p
# problems arise if the distributions differ

# distribution matters, scale doesn't