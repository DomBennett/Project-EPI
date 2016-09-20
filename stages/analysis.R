# ANALYSIS STAGE

# LIBS
library(ggplot2)
source(file.path("tools", "analysis_tools.R"))

# DIRS
input_dir <- "9_epi"
output_dir <- 'results'
input_file <- file.path(input_dir, "res.RData")

# LOAD
load(input_file)

# GENERAL PLOTS
# top living fossils and their indexes
plotEPI(cld_data, n=50)
# time-success plot
ggplot(data=cld_data, aes(x=log(tmsplt), y=log(cntrst_n), colour=pepi)) +
  geom_point() + theme_bw()
# hists of key stats
ggplot(data=cld_data, aes(epi)) + geom_density() + theme_bw()
ggplot(data=cld_data, aes(pepi)) + geom_density() + theme_bw()
# both look normal

# BASIC STATS
qntls <- quantile(cld_data[['pepi']], na.rm=TRUE)
cld_data$qntl <- NA
cld_data$qntl[cld_data[['pepi']] < qntls[2]] <- '0-25%'
cld_data$qntl[cld_data[['pepi']] > qntls[2] &
                cld_data[['pepi']] < qntls[4]] <- '25-75%'
cld_data$qntl[cld_data[['pepi']] > qntls[4]] <- '75-100%'
cld_data$qntl <- factor(cld_data$qntl, levels=c('0-25%', '25-75%', '75-100%'))
min(cld_data[['pepi']], na.rm=TRUE)
mean(cld_data[['pepi']], na.rm=TRUE)
max(cld_data[['pepi']], na.rm=TRUE)

# WIKI
sum(!is.na(cld_data$wiki))
p <- ggplot(cld_data, aes(x=pepi, y=as.numeric(wiki))) + binomial_smooth()
p + theme_bw() + xlab('pEPI') + ylab("'Living Fossil' in Wikipedia")

# WHICH LARGE CLADES HAVE LOW PEPI?
pull <- !is.na(cld_data[['pepi']]) &
  cld_data[['n']] > 100
#sum(pull)
#cld_data[['nm']][pull]
#cld_data[['txid']][pull]
plt_data <- cld_data[pull, c('nm', 'pepi', 'qntl')]
plt_data <- plt_data[order(plt_data[['pepi']], decreasing=FALSE), ]
plt_data[['nm']] <- factor(plt_data[['nm']], levels=plt_data[['nm']])
p <- ggplot(plt_data, aes(x=nm, y=pepi, colour=qntl)) +
  geom_point() + coord_flip() + xlab('') + ylab('pEPI')
p


# WHAT ABOUT US?
pull <- cld_data[['txid']] %in% c('9606', '207598', '9604', '314295', '9526')
cld_data[pull, ]  # apes are performing the worst


# COMPARE EPI, PEPI AND ED
pull <- !is.na(cld_data[['ed']]) & !is.na(cld_data[['epi']]) &
  cld_data[['epi']] != 0
plotEPI(cld_data[pull,], n=50)
# which have the biggest uncertainty?
m1 <- lm(epi~pepi, data=cld_data)

cld_data[['pepi']] < -10
sum(cld_data[['pepi']] < -10, na.rm=TRUE)
order(abs(m1$residuals))
cld_data[['nm']][]

pepi_rnkngs <- order(cld_data[['pepi']][pull], decreasing=FALSE)
epi_rnkngs <- order(cld_data[['epi']][pull], decreasing=FALSE)
ed_rnkngs <- order(cld_data[['ed']][pull], decreasing=TRUE)
cor(pepi_rnkngs, epi_rnkngs)
pepi_rnkngs[1:10]
epi_rnkngs[1:10]

# ED~PEPI
pull <- cld_data[['n']] == 1  # spp only
ggplot(data=cld_data[pull, ], aes(ed)) + geom_density() + theme_bw()
ggplot(data=cld_data[pull, ], aes(log(ed))) + geom_density() + theme_bw()
# logging makes ed normal
cor1 <- cor.test(log(cld_data[['ed']][pull]), cld_data[['pepi']][pull])
# strong correlation, -0.4
m1 <- lm(pepi~log(ed), data=cld_data[pull,])
summary(m1)  # strong relationship
ggplot(data=cld_data[pull,], aes(x=pepi, y=log(ed), colour=factor(txnmcgrp))) +
  geom_point() + stat_smooth(method='lm') + theme_bw()

# EPI and PEPI
cor.test(cld_data$pepi, cld_data$epi)  # correlate strongly, R = ~0.76
# this is what you would expect though, how much does change by itself correlate with pepi
pull <- !is.na(cld_data$cntrst_chng)
cor.test(cld_data$pepi[pull], cld_data$cntrst_chng[pull])
# log because not normal
pull <- !is.na(cld_data$cntrst_chng) & cld_data$cntrst_chng != 0
cor.test(cld_data$pepi[pull], log(cld_data$cntrst_chng[pull]))
ggplot(cld_data[pull,], aes(x=pepi, log(cntrst_chng))) + geom_point() + stat_smooth(method='lm')
m1 <- lm(pepi~log(cntrst_chng), data=cld_data[pull, ])
summary(m1)
plot(m1)  # weak, postive relationship. lower pepi, lower change

# MEASURABLES
pull <- !is.na(cld_data[['cntrst_n']]) &
  !is.na(cld_data[['cntrst_chng']]) &
  cld_data[['cntrst_chng']] != 0
ggplot(cld_data[pull,], aes(y=log(cntrst_n), x=log(tmsplt), colour=log(cntrst_chng))) + geom_point()
# WHATS THE RELATIONSHIP BETWEEN SUCCESS AND N?
lg_n <- log(cld_data$cntrst_n)
p_data <- data.frame(n=abs(lg_n - mean(lg_n)),
                     t=log(cld_data$tmsplt),
                     grp=cld_data$txnmcgrp)
ggplot(p_data, aes(y=n, x=t)) + geom_point() + stat_smooth(formula=y~exp(x), method='lm')


# OTHER
spp <- cld_data[['n']] == 1
unique(cld_data[spp, 'cntrst_chng'])  # error in change causing no change for species
hist(cld_data[spp, 'cntrst_chng'])
