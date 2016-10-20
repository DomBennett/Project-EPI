# ANALYSIS STAGE

# LIBS
library(ggplot2)
library(nlme)
library(gridExtra)
source(file.path("tools", "analysis_tools.R"))

# DIRS
input_dir <- "9_epi"
results_path <- 'results'
input_file <- file.path(input_dir, "res.RData")

# LOAD
load(input_file)

# GENERAL STATS
# pepi
sum(!is.na(cld_data$pepi))
sort(table(cld_data$txnmcgrp[!is.na(cld_data$pepi)], useNA="always"))
quantile(cld_data$pepi, na.rm=TRUE)
ggplot(cld_data, aes(pepi)) + geom_density()
# top mammals
mmls <- cld_data[cld_data$txnmcgrp == "mammals", ]
ordd <- order(mmls$pepi)
mmls[ordd[1:10], c('scinm', 'pepi')]
# top birds
brds <- cld_data[cld_data$txnmcgrp == "birds", ]
ordd <- order(brds$pepi)
brds[ordd[1:10], c('scinm', 'pepi')]
# top plants
plnts <- cld_data[cld_data$txnmcgrp == "plants", ]
ordd <- order(plnts$pepi)
plnts[ordd[1:10], c('scinm', 'pepi')]
# top arthropods
arthrpds <- cld_data[cld_data$txnmcgrp == "arthropods", ]
ordd <- order(arthrpds$pepi)
arthrpds[ordd[1:10], c('scinm', 'pepi')]
# top bony fish
bfs <- cld_data[cld_data$txnmcgrp == "bony_fish", ]
ordd <- order(bfs$pepi)
bfs[ordd[1:10], c('scinm', 'pepi')]
# top metazoa
bfs <- cld_data[cld_data$txnmcgrp != "plants", ]
ordd <- order(bfs$pepi)
bfs[ordd[1:10], c('scinm', 'pepi')]

# epi
sum(!is.na(cld_data$epi))
sort(table(cld_data$txnmcgrp[!is.na(cld_data$epi)], useNA="always"))
quantile(cld_data$epi, na.rm=TRUE)
ggplot(cld_data, aes(epi)) + geom_density()
# top mammals
mmls <- cld_data[cld_data$txnmcgrp == "mammals", ]
ordd <- order(mmls$epi)
mmls[ordd[1:10], c('scinm', 'epi')]
# top birds
brds <- cld_data[cld_data$txnmcgrp == "birds", ]
ordd <- order(brds$epi)
brds[ordd[1:10], c('scinm', 'epi')]

# WHAT ABOUT US?
pull <- cld_data[['txid']] %in% c('9606', '207598', '9604', '314295', '9526')
cld_data[pull, ]  # apes are performing the worst

# CORRELATION BETWEEN EPI AND PEPI AND ED
cor.test(cld_data$epi, cld_data$pepi, method='pearson')
cor.test(cld_data$epi, cld_data$pepi, method="spearman")
cor.test(cld_data$pepi, cld_data$cntrst_chng, method="spearman")
cor.test(log(cld_data$ed), cld_data$pepi, method='pearson')
cor.test(log(cld_data$ed), cld_data$epi, method="pearson")

# VARIANCE FROM EACH VARIABLE
pull <- !is.na(cld_data$epi) & !is.na(cld_data$pepi)
epi_data <- cld_data[pull, ]
sd(epi_data$cntrst_n)/mean(epi_data$cntrst_n)
sd(epi_data$tmsplt, na.rm=TRUE)/mean(epi_data$tmsplt, na.rm=TRUE)
sd(epi_data$cntrst_chng, na.rm=TRUE)/mean(epi_data$cntrst_chng, na.rm=TRUE)

# CORRELATION BETWEEN MEASURABLES AND EPIs
cor.test(cld_data$cntrst_n, cld_data$pepi, method='spearman')
cor.test(cld_data$tmsplt, cld_data$pepi, method='spearman')
cor.test(epi_data$cntrst_n, epi_data$epi, method='spearman')
cor.test(epi_data$cntrst_chng, epi_data$epi, method='spearman')
cor.test(epi_data$tmsplt, epi_data$epi, method='spearman')

# WHICH MOVED THE MOST?
pull <- !is.na(cld_data$epi) & !is.na(cld_data$pepi)
pull <- pull & cld_data$n > 2
epi_data <- cld_data[pull, ]
epi_data$z_epi <- (epi_data$epi - mean(epi_data$epi))/
  sd(epi_data$epi)
epi_data$z_pepi <- (epi_data$pepi - mean(epi_data$pepi))/
  sd(epi_data$pepi)
epi_data$resids <- epi_data$z_epi - epi_data$z_pepi
# higher pepi than epi
ordd <- order(epi_data$resids)
epi_data[ordd[1:10], c('scinm', 'epi', 'z_epi', 'pepi', 'z_pepi', 'resids')]
# higher epi than pepi
ordd <- order(epi_data$resids, decreasing=TRUE)
epi_data[ordd[1:10], c('scinm', 'cntrst_chng', 'epi', 'z_epi', 'pepi', 'z_pepi', 'resids')]
ggplot(epi_data, aes(resids, log(n))) +geom_point()

# WIKI
sum(!is.na(cld_data$wiki))
sum(cld_data$wiki, na.rm=TRUE)
# mentions
getWikiByQuantile('epi')
getWikiByQuantile('pepi')
getWikiByQuantile('ed')
# model
pull <- cld_data[['n']] == 1 & !is.na(cld_data[['ed']]) &
  !is.na(cld_data[['pepi']]) & !is.na(cld_data[['epi']])
m1 <- glm(wiki~1,data=cld_data[pull,],family=binomial())
summary(m1)
m2 <- glm(wiki~log(ed),data=cld_data[pull,],family=binomial())
summary(m2)
m3 <- glm(wiki~pepi,data=cld_data[pull,],family=binomial())
summary(m3)
m4 <- glm(wiki~epi,data=cld_data[pull,],family=binomial())
summary(m4)
exp(coef(m4))
anova(m1, m2, test="Chisq")
anova(m1, m3, test="Chisq")
anova(m1, m4, test="Chisq")
anova(m2, m3, test="Chisq")
anova(m2, m4, test="Chisq")
anova(m2, m3, test="Chisq")  # ed, epi and pepi have equal explanatory power

# TABLE 1
pull <- c(which(cld_data$scinm == 'Amphibia'), 
          which(cld_data$txnmcgrp == 'amphibia')[1:5],
          which(cld_data$scinm == 'Arthropoda'), 
          which(cld_data$txnmcgrp == 'arthropods')[1:5],
          which(cld_data$scinm == 'Aves'), 
          which(cld_data$txnmcgrp == 'birds')[1:5],
          which(cld_data$scinm == 'Embryophyta'), 
          which(cld_data$txnmcgrp == 'plants')[1:5],
          which(cld_data$scinm == 'Lepidosauria'), 
          which(cld_data$txnmcgrp == 'lepidosaurs')[1:5],
          which(cld_data$scinm == 'Mammalia'), 
          which(cld_data$txnmcgrp == 'mammals')[1:5],
          which(cld_data$scinm == 'Metazoa'), 
          which(cld_data$txnmcgrp == 'metazoa')[1:5],
          which(cld_data$scinm == 'Vertebrata'), 
          which(cld_data$txnmcgrp == 'vertebrates')[1:5],
          which(cld_data$scinm == 'Teleostei'), 
          which(cld_data$txnmcgrp == 'bony_fish')[1:5])
write.csv(cld_data[pull, ], file=file.path(results_path, 'table1.csv'))

# PLOTS
# correlation
p_epipepi <- ggplot(cld_data, aes(epi, pepi)) + geom_point(alpha=.5) +
  stat_smooth(method='lm') + theme_bw() + xlab('EPI') + ylab('pEPI') + ggtitle('A')
p_epied <- ggplot(cld_data, aes(epi, log(ed))) + geom_point(alpha=.5) +
  stat_smooth(method='lm') + theme_bw() + xlab('EPI') + ylab('log(ED)') + ggtitle('B')
p_pepied <- ggplot(cld_data, aes(pepi, log(ed))) + geom_point(alpha=.5) +
  stat_smooth(method='lm') + theme_bw() + xlab('pEPI') + ylab('log(ED)') + ggtitle('C')
tiff(file.path(results_path, 'figure_2.tiff'), width=18, height=9, units="cm",
     res=1200)
grid.arrange(p_epipepi, p_epied, p_pepied, ncol=3)
dev.off()
# data dist.
p_pepi <- ggplot(cld_data, aes(log(tmsplt), log(cntrst_n), colour=pepi)) +
  geom_point(alpha=.5) + theme_bw() + xlab('log(Time since split)') +
  ylab('log(Contrasted N)') + ggtitle('A')
p_epi <- ggplot(cld_data, aes(log(tmsplt), log(cntrst_chng), colour=epi)) +
  geom_point(alpha=.5) + theme_bw() + xlab('log(Time since split)') +
  ylab('log(Contrasted Change Score)') + ggtitle('B')
tiff(file.path(results_path, 'figure_3.tiff'), width=18, height=9, units="cm",
     res=1200)
grid.arrange(p_pepi, p_epi, ncol=2)
dev.off()
# wiki
p_EPI <- ggplot(cld_data, aes(x=epi, y=as.numeric(wiki))) +
  binomial_smooth()
p_EPI <- p_EPI + theme_bw() + xlab('EPI') + ylab("PROBABILITY") +
  ylim(0, 1) + ggtitle("A")
p_pEPI <- ggplot(cld_data, aes(x=pepi, y=as.numeric(wiki))) +
  binomial_smooth()
p_pEPI <- p_pEPI + theme_bw() + xlab('pEPI') + ylab("") + ylim(0, 1) + ggtitle("B")
p_ED <- ggplot(cld_data[cld_data[['n']] == 1, ], aes(x=log(ed), y=as.numeric(wiki))) +
  binomial_smooth()
p_ED <- p_ED + theme_bw() + xlab('ED') + ylab("")  + ylim(0, 1) + ggtitle("C")
tiff(file.path(results_path, 'figure_4.tiff'), width=18, height=9, units="cm",
     res=1200)
grid.arrange(p_EPI, p_pEPI, p_ED, ncol=3)
dev.off()


# GENERAL PLOTS
# top living fossils and their indexes
plotEPI(cld_data, n=50, order_by='pepi')
plotEPI(cld_data, n=50, order_by='epi')
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

# Correlations between the measurables
cor.test(log(cld_data$cntrst_n), log(cld_data$tmsplt))
pull <- cld_data$txnmcgrp == 'mammals'
cor.test(log(cld_data$cntrst_n[pull]), log(cld_data$tmsplt[pull]))
pull <- cld_data$txnmcgrp == 'birds'
cor.test(log(cld_data$cntrst_n[pull]), log(cld_data$tmsplt[pull]))
pull <- cld_data$txnmcgrp == 'plants'
cor.test(log(cld_data$cntrst_n[pull]), log(cld_data$tmsplt[pull]))

ggplot(cld_data[pull,], aes(log(tmsplt), log(cntrst_n), colour=txnmcgrp)) + geom_point() +
  stat_smooth(method="lm")
ggplot(cld_data, aes(log(cntrst_n))) + geom_density()






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
cor.test(cld_data$pepi, cld_data$epi)  # correlate strongly, R = ~0.84
# this is what you would expect though, how much does change by itself correlate with pepi
pull <- !is.na(cld_data$cntrst_chng)
cor.test(cld_data$pepi[pull], cld_data$cntrst_chng[pull])
# weak but significant correlation between chng and pepi, -0.11
# log because not normal
hist(cld_data$cntrst_chng)
pull <- !is.na(cld_data$cntrst_chng) & cld_data$cntrst_chng != 0
hist(log(cld_data$cntrst_chng[pull]))
cor.test(cld_data$pepi[pull], log(cld_data$cntrst_chng[pull]))
# no longer significant, -0.01
ggplot(cld_data[pull,], aes(x=pepi, log(cntrst_chng), colour=txnmcgrp)) +
  geom_point() + stat_smooth(method='lm')
ggplot(cld_data[pull,], aes(x=log(cntrst_n), log(cntrst_chng), colour=txnmcgrp)) +
  geom_point() + stat_smooth(method='lm')
pull_mmls <- cld_data[['txnmcgrp']] == 'mammal' & pull
m_mml <- lm(pepi~log(cntrst_chng), data=cld_data[pull_mmls, ])
summary(m_mml)
pull_brds <- cld_data[['txnmcgrp']] == 'bird' & pull
m_brd <- lm(pepi~log(cntrst_chng), data=cld_data[pull_brds, ])
summary(m_brd)

# CHNG
pull <- !is.na(cld_data$cntrst_chng) & cld_data$cntrst_chng < 20
ggplot(cld_data[pull,], aes(log(cntrst_chng), color=txnmcgrp, fill=txnmcgrp)) +
  geom_density(alpha=0.75) + theme_bw() + xlab('CONTRAST CHANGE')
pull <- !is.na(cld_data$cntrst_chng) & cld_data$n == 1
ggplot(cld_data[pull,], aes(log(cntrst_chng), color=txnmcgrp, fill=txnmcgrp)) +
  geom_density(alpha=0.75) + theme_bw() + xlab('CONTRAST CHANGE')
pull <- !is.na(cld_data$cntrst_chng) & cld_data$n > 1
ggplot(cld_data[pull,], aes(log(cntrst_chng), color=txnmcgrp, fill=txnmcgrp)) +
  geom_density(alpha=0.75) + theme_bw() + xlab('CONTRAST CHANGE')
# bimodal distributions for mammals and birds, not to do with species
# is there a correlation between n and chng?
ggplot(cld_data[pull,], aes(y=log(cntrst_chng), x=log(n), colour=txnmcgrp)) +
  geom_point(alpha=0.75) + stat_smooth(method='lm')
# not really....
# it's probably simply to do with difference sets of characters used for contrasts
# TODO: look up sisters and show that the difference is normal
cntrst_chng_diff <- rep(NA, nrow(cld_data))
for(i in 1:nrow(cld_data)) {
  chng1 <- cld_data[i, 'cntrst_chng']
  chng2 <- cld_data[["cntrst_chng"]][cld_data[i, 'tmsplt'] == cld_data[['tmsplt']]]
  cntrst_chng_diff[i] <- abs(log(chng1) - log(mean(chng2, na.rm=TRUE)))
}
hist(log(cntrst_chng_diff))


# it must be to do with differing sets of characters

# work out if bimodality is due to sisters
above_below <- rep(NA, nrow(cld_data))
for(i in 1:nrow(cld_data)) {
  chng1 <- cld_data[i, 'cntrst_chng']
  chng2 <- cld_data[["cntrst_chng"]][cld_data[i, 'tmsplt'] == cld_data[['tmsplt']]]
  above_below[i] <- log(chng1) > log(mean(chng2, na.rm=TRUE))
}

# N
pull <- !is.na(cld_data$cntrst_n)
ggplot(cld_data[pull,], aes(log(cntrst_n), color=txnmcgrp, fill=txnmcgrp)) +
  geom_density(alpha=0.75) + theme_bw() + xlab('CONTRAST N')

# N + CHNG
ggplot(cld_data, aes(x=log(cntrst_chng), y=log(cntrst_n))) + geom_point(alpha=0.75)
ggplot(cld_data, aes(x=log(cntrst_chng), y=log(tmsplt))) + geom_point(alpha=0.75)


# MEASURABLES
pull <- !is.na(cld_data[['cntrst_n']]) &
  !is.na(cld_data[['cntrst_chng']]) &
  cld_data[['cntrst_chng']] != 0
ggplot(cld_data[pull,], aes(y=log(cntrst_n), x=log(tmsplt), colour=log(cntrst_chng))) + geom_point()
ggplot(cld_data[pull,], aes(y=log(cntrst_n + cntrst_chng), x=log(tmsplt), colour=epi)) + geom_point()
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
