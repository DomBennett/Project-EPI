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
phylo <- data[["phylo"]]
chars <- data[["chars"]]
rm(data)
# # test with mangeable size
# phylo <- drop.tip(phylo, sample(phylo$tip.label, length(phylo$tip.label) - 1000))
# chars <- chars[phylo$tip.label, sample(1:length(chars), size=100)]

# PROCESS
cat("Reducing character matrix ...\n")
chars <- reduceChrctrMtrx(chars, pcut=0.95)
cat("Estimating ancestral node states ...\n")
reconstruction.obj <- parsimonyReconstruction(chars, phylo)
cat("Calculating edge changes ...\n")
phylo <- calcChange(phylo, reconstruction.obj, parallel=parallel)
cat("Calculating fair proportion ...\n")
phylo <- calcTime(phylo)
cat("Calculating N descendants ...\n")
phylo <- calcSuccess(phylo)

# OUTPUT
cat("Outputting ...\n")
save(phylo, file = file.path(output.dir, paste0(stdy_grp, ".RData")))

# END
cat(paste0('\nStage `mesurables` finished at [', Sys.time(), ']\n'))