## No copyright, no warranty
## Dominic John Bennett
## 25/02/2014

## Directories
input.dir <- "1_measures/mammal_all"
output.dir <- "2_metrics"
if (!file.exists (output.dir)) {
  dir.create (output.dir)
}

## Libraries
cat ("Loading libraries ...\n")
require (calibrate)
require (ggplot2)
require (scales)
source (file.path ("functions", "LFITools.R"))
source (file.path ("functions", "EcoDataTools.R"))

## Sanity check
load (file.path (input.dir, "PhyloMeasures.RData"))
#plotEdgeChanges (phylo,by.char=TRUE)
#phylo$edge.change.obj <- NULL

## Reading in measures
measures <- read.csv (file.path (input.dir, "measures.csv"), stringsAsFactors = FALSE)
measures <- measures[ ,names (measures) != 'sd.ED']
cat (paste0 ("[", sum (complete.cases (measures)), "] of [", nrow (measures),
             "] with complete data ... \n"))
max(measures$n)
colSums(is.na(measures))
measures <- na.omit (measures)

## Creating metrics
## Metric 1 (rescale to 0-1)
# removing big clades
#measures <- measures[measures$n < 40, ]
#measures <- measures[measures$s.edge.change > 0.01,]
#time <- log (measures$mean.ED)
#time <- time / max (time)
#change <- log (measures$mean.change)
#change <- change/max (change)
#performance <- log (measures$n)
#performance <- performance / max (performance)
#lfi <- lfiChecker (time, change, performance, cut = 0.99)

time <- log (measures$contrast.ED)
time <- time / max (time)
performance <- log (measures$contrast.n)
performance <- performance / max (performance)
change <- log (measures$contrast.change)
change <- change/max (change)
lfi <- lfiChecker (time, change, performance, cut = 0.99)

## Exploring LFI
#cutoff <- quantile (lfi, probs = 0.99)
measures$lfi <- lfi
measures$time <- time
measures$performance <- performance
measures$change <- change
measures <- measures[order (lfi, decreasing = TRUE), ]
measures[1:10, -2]

## Export
write.csv (x = measures, file = file.path (output.dir, "livingfossils.csv"),
           row.names = FALSE)

## PCA
pca.res <- prcomp (measures[ ,names (measures) %in% c ("time", "performance",
                                                                   "change")],
                   scale. = TRUE, center = TRUE)
pca.x <- as.data.frame(pca.res$x)
pca.rot <- as.data.frame (pca.res$rotation)
prop.var <- round (sapply (pca.res$sdev^2,
                           function (x) Reduce('+', x)/sum (pca.res$sdev^2)), 3)
comparisons <- list (c ("PC1", "PC2"), c ("PC2", "PC3"), c ("PC1", "PC3"))
for (comp in comparisons) {
  p <- ggplot (pca.x, aes_string (x = comp[1], y = comp[2], size = 2)) +
    geom_point (aes (colour = measures$lfi)) +
    scale_colour_gradient2 (low = "red", high = "blue") +
    xlab (paste0 (comp[1], " (", prop.var[1], ")")) +
    ylab (paste0 (comp[2], " (", prop.var[2], ")"))
  print (p)
  plot(x = pca.rot[ ,comp[1]], y =  pca.rot[ ,comp[2]], xlab = comp[1],
       ylab = comp[2], cex = 0.5, pch = 19)
  text (x = pca.rot[ ,comp[1]], y =  pca.rot[ ,comp[2]],
        rownames (pca.rot[comp[1]]), adj = 1)
  prop.var <- prop.var[-c(1,2)]
}