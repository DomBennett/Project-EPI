
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
metrics <- calcMetrics(phylo)
time <- log (metrics$contrast.ed)
time <- time / max (time, na.rm=TRUE)
success <- log (metrics$contrast.n)
success <- success / max (success, na.rm=TRUE)
change <- log (metrics$contrast.change)  # +1?
change <- change/max (change, na.rm=TRUE)
metrics$epi <- ((change + success)/2) - time
pdf(file.path('figures', paste0(stdy_grp, '_epicheck.pdf')))
EPIChecker(time, change, success, 0.05)
dev.off()
metrics$success <- success
metrics$time <- time
metrics$change <- change


# OUTPUT
cat("Outputting ...\n")
save(metrics, file = file.path(output.dir, paste0(stdy_grp, ".RData")))

# END
cat(paste0('\nStage `epi` finished at [', Sys.time(), ']\n'))