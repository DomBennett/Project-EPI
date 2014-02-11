## No copyright, no warranty
## Dominic John Bennett
## Automated calculation of LFI from Dryad data
## 07/02/2014
## TODO: Generate a random distribution of phylogenies to account for polytomies

## Libraries
source (file.path ("Functions", "LFITools.R"))
source (file.path ("Functions", "EcoDataTools.R"))

## Directories
input.dir <- "0_data"
output.dir <- "1_LFI"
if (!file.exists (output.dir)) {
  dir.create (output.dir)
}

## Import data (not automated...)
data <- read.nexus.data (file.path (input.dir, "ladybird_matrix.nex"))
data <- lapply (data, function (x) x[2254:2365])
chars <- matrix (unlist (data), nrow = length (data), byrow = TRUE)
rownames (chars) <- names (data)
phylo <- read.nexus (file.path (input.dir, "ladybird_tree.nex"))[[1]]
if (!all (c (phylo$tip.label %in% names (data), names (data) %in% phylo$tip.label))) {
  stop ("All labels do not match between data and phylogeny.")
}
# start with mangeable size
#phylo <- drop.tip (phylo, sample (phylo$tip.label, length (phylo$tip.label) - 10))
#chars <- chars[phylo$tip.label, 1:10]
if (!is.ultrametric (phylo)) {
  phylo <- multi2di (phylo)
}

## Calculate LFI
reconstruction.obj <- parsimonyReconstruction (chars, phylo)
phylo <- calcBranchChanges (phylo, reconstruction.obj)
plot(phylo)
edgelabels(text = phylo$edge.changes)