# Run this script to wrangle 0_data/raw files for pipleine

# LIBS
library(ape)
library(plyr)
source (file.path('tools', "ReadNexusData.R"))

# FUNCTIONS
calcScope <- function (r.phylo.names, f.phylo) {
  calcEachScope <- function (i) {
    r.phylo.dists <- abs (vcv (r.phylo, corr = TRUE) - 1)
    r.dists.i <- r.phylo.dists[i, -i]
    f.dists.i <- f.phylo.dists[i, -i]
    r.nn <- names (which (r.dists.i == min (r.dists.i)))[1]
    length (f.dists.i[f.dists.i < f.dists.i[names (f.dists.i) == r.nn]]) + 1
  }
  r.phylo <- drop.tip (f.phylo, tip = f.phylo$tip.label[!f.phylo$tip.label %in% r.phylo.names])
  f.phylo.dists <- abs (vcv (f.phylo, corr = TRUE) - 1)
  r.phylo$scope <- mdply (.data = data.frame (i = 1: length (r.phylo$tip.label)),
                          .fun = calcEachScope)[ ,1]
  r.phylo
}

matchPhyChar <- function (phylo, chars) {
  chars <- chars[rownames (chars) %in% phylo$tip.label, ]
  if (!all (phylo$tip.label %in% rownames (chars))) {
    phylo <- calcScope (rownames (chars), phylo)
  }
  if (!all (c (phylo$tip.label %in% rownames (chars),
               rownames (chars) %in% phylo$tip.label))) {
    stop ("All labels do not match between data and phylogeny.")
  }
  return (list (phylo = phylo, chars = chars))
}


# DIRS
output.dir <- "0_data"
input.dir <- file.path(output.dir, "raw")

# DATA -- Ladybird
data <- read.nexus.data (file.path (input.dir, "ladybird_matrix.nex"))
data <- lapply (data, function (x) x[2254:2365])
chars <- matrix (unlist (data), nrow = length (data), byrow = TRUE)
rownames (chars) <- names (data)
rownames (chars) <- sub (" ", "_", rownames (chars))
chars[chars == "?"] <- rep (NA, sum (chars == "?"))
phylo <- read.nexus (file.path (input.dir, "ladybird_tree.nex"))[[1]]
data <- matchPhyChar (phylo, chars)
save (data, file = file.path (output.dir, "ladybird.RData"))

# DATA -- Mammal (http://esapubs.org/archive/ecol/e090/184/)
data <- read.delim (file.path (input.dir, "panTHERIA.txt"), na.strings = -999,
                    stringsAsFactors = FALSE)
pantheria <- data[ ,- c (1:5,36:55)]
for (i in 1:ncol (pantheria)) {
  temp.chars <- pantheria[!is.na (pantheria[ , i]),i]
  if (length (unique (temp.chars)) > 10) {
    pantheria[!is.na (pantheria[ , i]),i] <- cut(temp.chars, 10)
  }
}
rownames (pantheria) <- data[ ,"MSW93_Binomial"]
rownames (pantheria) <- sub (" ", "_", rownames (pantheria))
file <- "morpho_matrix_forR.nex"
oleary <- readNexusData (file.path (input.dir, file))
chars <- merge (oleary, pantheria, by = 0, all = TRUE)
rownames (chars) <- c (rownames (oleary) [!rownames (oleary) %in% rownames (pantheria)], rownames (pantheria))
phylo <- read.tree (file.path (input.dir, "bininda.txt"))
clades <- MoreTreeTools::getClades(phylo)
children <- list()
children[clades$clade.node] <- gsub("_", " ", clades$clade.children)
#clade_labels <- paste0('n', 1:(length(phylo$tip.label) + phylo$Nnode))
clade_labels <- MoreTreeTools::getNodeLabels(phylo, all=FALSE, cache=TRUE,
                                             parent="Mammalia")
clade_labels <- paste0(clade_labels, '_', 1:phylo$Nnode)  # prevent duplicates
clade_labels <- c(gsub("_", " ", phylo$tip.label), clade_labels)
phylo$clade_labels <- clade_labels
chars <- chars[rownames (chars) %in% phylo$tip.label, ]
# missing data
#sum(is.na (chars))*100/(nrow (chars) * ncol (chars))
data <- list (phylos = list(phylo), chars = chars, children=children)
save (data, file = file.path (output.dir, "mammal.RData"))

# DATA -- Birds
file <- "X1228_Morphology Matrix_morphobank.nex"
livezy <- readNexusData(file.path (input.dir, file))
phylos <- read.tree(file.path(input.dir, 'jetz_aves.tre'))
for(i in 1:length(phylos)) {
  clades <- MoreTreeTools::getClades(phylos[[i]])
  children <- list()
  children[clades$clade.node] <- gsub("_", " ", clades$clade.children)
  clade_labels <- MoreTreeTools::getNodeLabels(phylos[[i]], all=FALSE, cache=TRUE,
                                               parent="Aves")
  clade_labels <- paste0(clade_labels, '_', 1:phylos[[i]]$Nnode)  # prevent duplicates
  clade_labels <- c(gsub("_", " ", phylos[[i]]$tip.label), clade_labels)
  phylos[[i]]$clade_labels <- clade_labels
}
# character data is for whole groups
# use character matching to assign the same value to memebers of the same group
# NAs for missing taxa
# use last phylo in loop
livezy_mod <- matrix(data=NA, nrow=length(phylos[i]$tip.label), ncol=ncol(livezy))
rownames(livezy_mod) <- phylos[i]$tip.label
for(i in 1:nrow(livezy)) {
  mtchs <- which(grepl(rownames(livezy)[i], phylos[i]$tip.label))
  if(length(mtchs) > 1) {
    for(j in mtchs) {
      livezy_mod[phylos[i]$tip.label[j], ] <- livezy[i, ]
    }
  }
}
data <- list (phylo = phylos, chars = livezy_mod, children=children)
save (data, file = file.path (output.dir, "bird.RData"))
