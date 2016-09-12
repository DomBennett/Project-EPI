# CALCULATE EPIS FROM NODE_OBJ
# TODO: add ranks


# START
cat(paste0('\nStage `epi` started at [', Sys.time(), ']\n'))

# PARAMETERS
cat("Loading parameters ...\n")
source('parameters.R')

# LIBS
cat("Loading libraries ...\n")
library(ggplot2)
source(file.path("tools", "epi_tools.R"))
source(file.path('tools', 'node_obj_tools.R'))

# DIRS
cc_dir <- "8_cntrst_chng"
tt_dir <- "7_timetree"
output_dir <- '9_epi'
if (!file.exists(output_dir)) {
  dir.create(output_dir)
}
output_file <- file.path(output_dir, "res.RData")
cc_file <- file.path(cc_dir, "res.RData")
tt_file <- file.path(tt_dir, "res.RData")

# CLADE DATA FRAME
cat("Generating clade dataframe.... ")
load(cc_file)
cnddts <- ls(node_obj)
cld_data <- genDataframe(cnddts, node_obj)
load(tt_file)
cld_data <- rbind(cld_data, genDataframe(cnddts, tt_obj))
rm(tt_obj, node_obj)
cat("Done.\n")

# CLEAN-UP
# some data points don't make sense, reasons as yet unknown
cld_data <- cld_data[cld_data[['cntrst_n']] != Inf, ]  # two instances of this

# CALCULATE EPIS
cat("Calculating EPIs .... ")
cld_data$epi <- log((cld_data[['cntrst_chng']] + cld_data[['cntrst_n']]) / cld_data[['tmsplt']])
cld_data$pepi <- log(cld_data[['cntrst_n']] / cld_data[['tmsplt']])
cat("Done.\n")

# PLOTS
plotEPI(epi, n=25)




pdf(file.path(output_dir, "res.pdf"))
ggplot(data=p_data, aes(x=log(tmsplt), y=log(cntrst_n), colour=pepi)) +
  geom_point() + theme_bw()
ggplot(data=epi, aes(x=pepi, y=log(ed))) +
  geom_point() + theme_bw()
ggplot(data=epi, aes(x=pepi, y=chng)) +
  geom_point() + theme_bw()
dev.off()

# EPI and PEPI
cor.test(cld_data$pepi, cld_data$epi)  # correlate strongly, R = ~0.7
# this is what you would expect though, how much does change by itself correlate with pepi
pull <- !is.na(cld_data$cntrst_chng)
cor.test(cld_data$pepi[pull], cld_data$cntrst_chng[pull])
# TODO: workout why change has NA, Nan and 1

# ADDING TAXONOMIC INFO
txids <- ls(node_obj)
epi$txnmcgrp <- NA
mmls <- getGrpTxids(txids, "mammals")
epi$txnmcgrp[epi$txid %in% mmls] <- 'mammal'
brds <- getGrpTxids(txids, "birds")
epi$txnmcgrp[epi$txid %in% brds] <- 'bird'
lpdsrs <- getGrpTxids(txids, "lepidosaurs")
epi$txnmcgrp[epi$txid %in% lpdsrs] <- 'lepidosaur'

# OUTPUT
cat("Outputting ... ")
top50 <- epi[order(epi[['pepi']])[1:50], ]
write.csv(top50, file.path(output_dir, "top50.csv"), row.names=FALSE)
save(node_obj, epi, file=output_file)
cat('Done.\n')

# END
cat(paste0('\nStage `epi` finished at [', Sys.time(), ']\n'))