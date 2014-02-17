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
plot (phylo)
edgelabels (text = phylo$edge.changes)
lfi.data <- calcLFI (phylo)

## Testing Metrics
# function
lfiChecker <- function (time, change, ep, model = 1) {
  hist(time)
  hist(change)
  hist(ep)
  if (model == 1) {
    lfi <- time / change * ep
  } else if (model == 2) {
    lfi <- time / (change + ep)
  } else {
    stop ("Invalid model number")
  }
  plot(time ~ lfi)
  abline (lm (time ~ lfi))
  plot(change ~ lfi)
  abline (lm (change ~ lfi))
  plot(ep ~ lfi)
  abline (lm (ep ~ lfi))
  lfi
}
## Metric 1
#time <- lfi.data$time
#change <- lfi.data$s.edge.change + lfi.data$d.edge.change
#ep <- lfi.data$d.edge.length
#lfiChecker (time, change + 1, ep + 1, model = 1)

## Metric 1 (logged)
#time <- lfi.data$time
#change <- log (lfi.data$s.edge.change + lfi.data$d.edge.change + 1)
#ep <- log (lfi.data$d.edge.length + 1)
#lfiChecker (time, change, ep, model = 1) # produces NaNs because of 0s

## Metric 1 (sqrt)
#time <- lfi.data$time
#change <- sqrt (lfi.data$s.edge.change + lfi.data$d.edge.change)
#ep <- sqrt (lfi.data$d.edge.length)
#lfiChecker (time, change, ep, model = 1)

## Metric 1 (reset to 1 as min)
time <- lfi.data$time
change <- lfi.data$s.edge.change + lfi.data$d.edge.change
change <- change + (1 - min (change))
ep <- lfi.data$d.edge.length
ep <- ep  + (1 - min (ep))
lfi <- lfiChecker (time, change, ep, model = 1)

## Exploring LFI
hist (lfi)
lfi.data$lfi <- lfi
cutoff <- quantile (lfi, probs = 0.95)
abline (v = cutoff, col = "red")
living.fossils <- lfi.data[lfi > cutoff, ]
living.fossils <- living.fossils[order (living.fossils$lfi), ]
plot (phylo)
nodelables (text = living.fossils$lfi, node = living.fossils$node)
