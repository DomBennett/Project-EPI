## No copyright, no warranty
## Dominic John Bennett
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

## Calculate LFI
reconstruction.obj <- parsimonyReconstruction (chars, phylo)
phylo <- calcBranchChanges (phylo, reconstruction.obj, plot.test = FALSE)
plot (phylo)
edgelabels (text = phylo$edge.changes)
lfi.data <- calcLFI (phylo)

## Testing Metrics

## Metric 1 (rescale to 0-1)
# removing big clades
#lfi.data <- lfi.data[lfi.data$n < 40, ]
time <- lfi.data$time / max (lfi.data$time)
change <- lfi.data$s.edge.change + lfi.data$d.edge.change
change <- change / max (change)
performance <- lfi.data$d.edge.length
performance <- performance / max (performance)
lfi <- lfiChecker (time, change, performance)

## Exploring LFI
living.fossils <- lfi.data[lfi > cutoff, ]
living.fossils <- living.fossils[order (living.fossils$lfi, decreasing = TRUE), ]
plot (phylo)
nodelabels (text = round (living.fossils$lfi, digits = 3), node = living.fossils$node)