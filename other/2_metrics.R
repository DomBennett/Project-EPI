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
# proportion of missing values
#sum (is.na (phylo$edge.changes))*100/length (phylo$edge.length)
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
epi <- EPIChecker (time, change, performance, cut = 0.01)

## Exploring epi
measures$epi <- epi
# proportion between -0.5 and 0.5
#sum(measures$epi < 0.5 & measures$epi > -0.5)/nrow(measures)
# proportion below 1%
#cutoff <- quantile (epi, probs = 0.01)
#sum(measures$epi < cutoff)
measures$time <- time
measures$performance <- performance
measures$change <- change
measures <- measures[order (epi, decreasing = FALSE), ]
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
  EPI <- measures$epi
  p <- ggplot (pca.x, aes_string (x = comp[1], y = comp[2])) +
    geom_point (aes (colour = EPI, size = 2)) +
    scale_colour_gradient2 (low = "red", high = "blue") +
    xlab (paste0 (comp[1], " (", prop.var[1], ")")) +
    ylab (paste0 (comp[2], " (", prop.var[2], ")")) +
    guides (size = FALSE)
  print (p)
  plot(x = pca.rot[ ,comp[1]], y =  pca.rot[ ,comp[2]], xlab = comp[1],
       ylab = comp[2], cex = 0.5, pch = 19)
  text (x = pca.rot[ ,comp[1]], y =  pca.rot[ ,comp[2]],
        rownames (pca.rot[comp[1]]), adj = 1)
  prop.var <- prop.var[-c(1,2)]
}


## Report plot
# rodents node == 4520
# grepl ('Oryctolagus', measures[rodents, 2]) # no rabbits
# vesper bats are 7770
vesper <- which (measures$node == 7770)
rodents <- which (measures$node == 4520)
homo <- which (measures$clade == 'Homo_sapiens')
aardvark <- which (measures$clade == "Orycteropus_afer")
comp <- c ('PC1', 'PC2')
EPI <- measures$epi
prop.var <- round (sapply (pca.res$sdev^2,
                           function (x) Reduce('+', x)/sum (pca.res$sdev^2)), 3)
pdf ('PCA_EPI.pdf')
p <- ggplot (pca.x, aes_string (x = comp[1], y = comp[2])) +
  geom_point (aes (colour = EPI, size = 2)) +
  scale_colour_gradient2 (low = "red", high = "blue") +
  xlab (paste0 (comp[1], " (", prop.var[1], ")")) +
  ylab (paste0 (comp[2], " (", prop.var[2], ")")) +
  guides (size = FALSE) +
  annotate("text", label = 'Mo', x = pca.x[1,1], y = pca.x[1,2]+0.22, size = 3) +
  annotate("text", label = 'Ap', x = pca.x[2,1], y = pca.x[2,2]+0.22, size = 3) +
  annotate("text", label = 'Ro', x = pca.x[rodents,1], y = pca.x[rodents,2]+0.22, size = 3) +
  annotate("text", label = 'Ve', x = pca.x[vesper,1], y = pca.x[vesper,2]+0.22, size = 3) +
  annotate("text", label = 'Ho', x = pca.x[homo,1], y = pca.x[homo,2]+0.22, size = 3) +
  annotate("text", label = 'Th', x = pca.x[nrow (measures),1], y = pca.x[nrow (measures),2]+0.22, size = 3) +
  annotate("text", label = 'Or', x = pca.x[aardvark,1], y = pca.x[aardvark,2]+0.22, size = 3) +
  theme_bw()
print (p)
dev.off()
  