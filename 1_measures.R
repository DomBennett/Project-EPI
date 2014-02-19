## No copyright, no warranty
## Dominic John Bennett
## 07/02/2014
## TODO: Generate a random distribution of phylogenies to account for polytomies

## Libraries
source (file.path ("Functions", "LFITools.R"))
source (file.path ("Functions", "EcoDataTools.R"))

## Directories
input.dir <- "0_data"
output.dir <- "1_measures"
if (!file.exists (output.dir)) {
  dir.create (output.dir)
}

## Import data (not automated...)
# Ladybird
#data <- read.nexus.data (file.path (input.dir, "ladybird_matrix.nex"))
#data <- lapply (data, function (x) x[2254:2365])
#chars <- matrix (unlist (data), nrow = length (data), byrow = TRUE)
#rownames (chars) <- names (data)
#chars[chars == "?"] <- rep (NA, sum (chars == "?"))
#phylo <- read.nexus (file.path (input.dir, "ladybird_tree.nex"))[[1]]
# Mammal (http://esapubs.org/archive/ecol/e090/184/)
data <- read.delim (file.path (input.dir, "panTHERIA.txt"), na.strings = -999,
                                        stringsAsFactors = FALSE)
chars <- data[ ,- c (1:5,36:55)]
for (i in 1:ncol (chars)) {
  temp.chars <- chars[!is.na (chars[ , i]),i]
  if (length (unique (temp.chars)) > 10) {
    chars[!is.na (chars[ , i]),i] <- cut(temp.chars, 10)
  }
}
rownames (chars) <- data[ ,"MSW93_Binomial"]
phylo <- read.tree (file.path (input.dir, "bininda.txt"))
phylo$tip.label <- sub ("_", " ", phylo$tip.label)
chars <- chars[rownames (chars) %in% phylo$tip.label, ]
if (!any (phylo$tip.label %in% rownames (chars))) {
  phylo <- drop.tip (phylo,
                     tip = phylo$tip.label[!phylo$tip.label %in% rownames (chars)])
}
if (!all (c (phylo$tip.label %in% rownames (chars),
             rownames (chars) %in% phylo$tip.label))) {
  stop ("All labels do not match between data and phylogeny.")
}
# start with mangeable size
#phylo <- drop.tip (phylo, sample (phylo$tip.label, length (phylo$tip.label) - 10))
#chars <- chars[phylo$tip.label,]
if (!is.binary.tree (phylo)) {
  phylo <- multi2di (phylo)
}

## Calculate LFI.data
reconstruction.obj <- parsimonyReconstruction (chars, phylo)
phylo <- calcEdgeChanges (phylo, reconstruction.obj)
#plotEdgeChanges (phylo, by.char = TRUE)
lfi.data <- calcLFI (phylo)

## Output LFI.data
write.csv (x = lfi.data, file = file.path (output.dir, "measures.csv"), row.names = FALSE)
save (phylo, file = file.path (output.dir, "PhyloMeasures.RData"))