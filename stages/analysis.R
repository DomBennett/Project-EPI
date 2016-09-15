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
plotEPI(epi, n=50)
# time-success plot
ggplot(data=cld_data, aes(x=log(tmsplt), y=log(cntrst_n), colour=pepi)) +
  geom_point() + theme_bw()
# hists of key stats
ggplot(data=cld_data, aes(epi)) + geom_density() + theme_bw()
ggplot(data=cld_data, aes(pepi)) + geom_density() + theme_bw()
# both look normal

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
cor.test(cld_data$pepi, cld_data$epi)  # correlate strongly, R = ~0.7
# this is what you would expect though, how much does change by itself correlate with pepi
pull <- !is.na(cld_data$cntrst_chng)
cor.test(cld_data$pepi[pull], cld_data$cntrst_chng[pull])
# TODO: workout why change has NA, Nan and 1