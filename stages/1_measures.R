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
phylo <- read.tree (file.path("0_data/raw/mammalia_w_fossils.tre"))

## Calculate EPI.data at different time points
cat ("Calculating EPI measures ...\n")
phylo$tip.label <- sub ("_", " ", phylo$tip.label)
#phylo <- drop.tip (phylo, phylo$tip.label[!phylo$tip.label %in% sample (phylo$tip.label, 1000)])
phylo$edge.changes <- phylo$edge.length
phylo <- addNodeAges (phylo)
phylo.age <- max (phylo$node.age)
#phylo <- labelNodes (phylo)
phylo$node.label <- paste0 ('id', 1:(length (phylo$tip.label) + phylo$Nnode))
n.time.points <- 2
increment <- phylo.age/(n.time.points)
# getTimeslice does not work with 0.0
time.point <- seq (from = 0.001, to = (phylo.age - increment), by = increment - 0.001)
res <- list ()
for (i in 1:n.time.points) {
  past.phylo <- getTimeslice (phylo, time.point[i])
  EDs <- calcFairProportion (past.phylo)
  measures <- calcLFIMeasures (past.phylo, EDs)
  res <- c (res, list (measures))
}

## Output LFI.data
cat ("Outputting ...\n")
save (res, file = file.path (dir, "measures.RData"))