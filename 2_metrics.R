## No copyright, no warranty
## Dominic John Bennett
## 25/02/2014

## Directories
input.dir <- "1_measures/ladybird"
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
#load (file.path (input.dir, "PhyloMeasures.RData"))
#plotEdgeChanges (phylo,by.char=TRUE)
#phylo$edge.change.obj <- NULL

## Reading in measures
measures <- read.csv (file.path (input.dir, "measures.csv"))
cat (paste0 ("[", sum (complete.cases (measures)), "] of [", nrow (measures),
             "] with complete data ... \n"))
measures <- na.omit (measures)

## Creating metrics
## Metric 1 (rescale to 0-1)
# removing big clades
#measures <- measures[measures$n < 40, ]
time <- measures$time / max (measures$time)
change <- measures$s.edge.change + measures$d.edge.change
change <- change / max (change)
performance <- measures$d.edge.length
performance <- performance / max (performance)
lfi <- lfiChecker (time, change, performance, cut = 0.99)

## Exploring LFI
cutoff <- quantile (lfi, probs = 0.99)
measures$lfi <- lfi
measures <- measures[order (measures$lfi, decreasing = TRUE), ]
living.fossils <- measures[measures$lfi > cutoff, ]
#plot (phylo)
#nodelabels (text = round (living.fossils$lfi, digits = 3), node = living.fossils$node)

## Export
write.csv (x = living.fossils, file = file.path (output.dir, "livingfossils.csv"),
           row.names = FALSE)

## PCA
pca.res <- prcomp (measures[ ,!names (measures) %in% c("node", "clade", "lfi")],
                   scale. = TRUE, center = TRUE)
pca.x <- as.data.frame(pca.res$x)
pca.rot <- as.data.frame (pca.res$rotation)
prop.var <- round (sapply (pca.res$sdev^2,
                            function (x) Reduce('+', x)/sum (pca.res$sdev^2)), 3)
comparisons <- list (c ("PC1", "PC2"), c ("PC3", "PC4"), c ("PC5", "PC6"))
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