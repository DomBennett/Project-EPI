## No copyright, no warranty
## Dominic John Bennett
## 07/02/2014
## TODO: Generate a random distribution of phylogenies to account for polytomies

## Libraries
cat ("Loading libraries ...\n")
source (file.path ("functions", "LFITools.R"))
source (file.path ("functions", "taxaResolve.R"))
source (file.path ("functions", "EcoDataTools.R"))

## Directories
dir <- "1_measures"

## Import data (not automated...)
cat ("Importing data ...\n")
load (file.path(dir, "PhyloMeasures.RData"))

## Calculate LFI.data at different time points
cat ("Calculating LFI measures ...\n")
EDs <- calcFairProportion (phylo)
measures <- calcLFIMeasures (phylo, EDs)

## Output LFI.data
cat ("Outputting ...\n")
write.csv (x = measures, file = file.path (dir, "measures.csv"), row.names = FALSE)