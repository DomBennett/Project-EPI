## No copyright, no warranty
## Dominic John Bennett
## 07/02/2014
## TODO: Generate a random distribution of phylogenies to account for polytomies

## Libraries
cat ("Loading libraries ...\n")
source (file.path ("functions", "LFITools.R"))
source (file.path ("functions", "EcoDataTools.R"))

## Directories
input.dir <- "0_data"
output.dir <- "1_measures"
if (!file.exists (output.dir)) {
  dir.create (output.dir)
}

## Import data (not automated...)
cat ("Importing data ...\n")
load (file.path(input.dir, "mammal.RData")) # or ladybird
phylo <- data[["phylo"]]
chars <- data[["chars"]]
rm (data)
# start with mangeable size
#phylo <- drop.tip (phylo, sample (phylo$tip.label, length (phylo$tip.label) - 1000))
#chars <- chars[phylo$tip.label,]

## Calculate Ancestral states
cat ("Estimating node states ...\n")
reconstruction.obj <- parsimonyReconstruction (chars, phylo)
cat ("Calculating edge changes ...\n")
phylo <- calcEdgeChanges (phylo, reconstruction.obj)

cat ("Outputting ...\n")
save (phylo, file = file.path (output.dir, "PhyloMeasures.RData"))