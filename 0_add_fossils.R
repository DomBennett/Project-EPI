##Libraries
library(stringr)
library(plyr)
library(RJSONIO)
source("functions/TaxaResolve.R")
source("functions/EcoDataTools.R")
source("functions/LFITools.R")

cleanNames <- function (names) {
  # Remove name accronyms for name matching
  clean <- function (name) {
    name <- gsub ("\\s+cf\\.", "", name)
    name <- gsub ("\\s+sp\\.", "", name)
    name <- gsub ("\\s+ex\\.", "", name)
    name <- gsub ("\\s+gr\\.", "", name)
    name <- gsub ("\\s+aff\\.", "", name)
    name <- gsub ("\\s+n\\.", "", name)
    name <- gsub ("\\s+gen\\.", "", name)
    name <- gsub ("\\s+\\?", "", name)
    #name <- sub ("\\s+indet\\.", "", name)
  }
  mdply (.data = data.frame (name = names), .fun = clean)[ ,2]
}

## Input
output.dir <- "0_data"
input.dir <- file.path(output.dir, "raw")
phylo <- read.tree (file.path (input.dir, "bininda.txt"))
if (!is.binary.tree (phylo)) {
  phylo <- multi2di (phylo)
}
phylo$tip.label <- sub ("_", " ", phylo$tip.label)

## Process
cat ("Retrieving fossil records ... \n")
records <- palaeoPull("Mammalia", limit = 'all')
cat ("Cleaning names ... \n")
records$name <- cleanNames (records$name)
cat ('Labelling phylogeny ... \n')
phylo <- labelNodes (phylo)
cat ('Adding node ages ...\n')
phylo <- addNodeAges (phylo)
cat ('Adding fossils to phylogeny ...\n')
new.phylo <- addFossilsToPhylogeny (phylo, records)
cat ('Outputting ...\n')
write.tree (new.phylo, file.path (output.dir, "mammalia_w_fossils.tre"))