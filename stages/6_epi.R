
#TODO: Generate a random distribution of phylogenies to account for polytomies

# START
cat(paste0('\nStage `epi` started at [', Sys.time(), ']\n'))

# PARAMETERS
cat("Loading parameters ...\n")
source('parameters.R')

# LIBS
cat("Loading libraries ...\n")
library(ape)
library(plyr)
source(file.path("tools", "epi_tools.R"))

# DIRS
input.dir <- "1_measurables"
output.dir <- "2_epi"
if(!file.exists(output.dir)) {
  dir.create(output.dir)
}

# INPUT
cat("Importing data ...\n")
load(file.path(input.dir, paste0(stdy_grp, '.RData')))

# PROCESS
cat("Calculating EPI ...\n")
metrics <- calcMetricsPhylo(phylos[[1]])
time <- metrics$contrast.ed1
success <- metrics$contrast.n
success <- log(metrics$contrast.n)
change <- metrics$contrast.change
metrics$success <- success
metrics$time <- time
metrics$change <- change
metrics$epi <- ((change + success)/2) - time
metrics$epi_nc <- (success - time)
bool <- metrics$time.split < 50  # ignore less than 50MY
metrics$epi_nc[bool] <- NA
metrics$epi[bool] <- NA
pdf(file.path('figures', paste0(stdy_grp, '_epicheck.pdf')))
EPIChecker(metrics, 0.05)
dev.off()
#metrics$clade_label[order(metrics$epi_nc)[1:50]]
#metrics$clade_label[order(metrics$epi)[1:25]]

# OUTPUT
cat("Outputting ...\n")
save(metrics, file = file.path(output.dir, paste0(stdy_grp, ".RData")))

# END
cat(paste0('\nStage `epi` finished at [', Sys.time(), ']\n'))