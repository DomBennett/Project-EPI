# Generates an augmented phylo object with change, time and success values appended

# START
cat(paste0('\nStage `measurables` started at [', Sys.time(), ']\n'))

# PARAMETERS
cat("Loading parameters ...\n")
source('parameters.R')

# LIBS
cat("Loading functions ...\n")
library(ape)
library(plyr)
library (foreach)
library (doMC)
source(file.path("tools", "measurables_tools.R"))

# DIRS
input.dir <- "0_data"
output.dir <- "1_measurables"
if(!file.exists(output.dir)) {
  dir.create(output.dir)
}

# INPUT
cat("Importing data ...\n")
load(file.path(input.dir, paste0(stdy_grp, '.RData')))
phylos <- data[["phylos"]]
chars <- data[["chars"]]
rm(data)
# # test with mangeable size
# phylo <- drop.tip(phylo, sample(phylo$tip.label, length(phylo$tip.label) - 1000))
# chars <- chars[phylo$tip.label, sample(1:length(chars), size=100)]

# PROCESS
cat("Reducing character matrix ...\n")
chars <- reduceChrctrMtrx(chars, pcut=1)
counter <- foreach(i=1:length(phylos)) %dopar% {
  cat("Estimating ancestral node states ...\n")
  reconstruction.obj <- parsimonyReconstruction(chars, phylos[i])
  cat("Calculating edge changes ...\n")
  phylos[i] <- calcChange(phylos[i], reconstruction.obj)
  cat("Calculating fair proportion ...\n")
  phylos[i] <- calcTime(phylos[i])
  cat("Calculating N descendants ...\n")
  phylos[i] <- calcSuccess(phylos[i])
}

# OUTPUT
cat("Outputting ...\n")
save(phylo, file = file.path(output.dir, paste0(stdy_grp, ".RData")))

# END
cat(paste0('\nStage `mesurables` finished at [', Sys.time(), ']\n'))