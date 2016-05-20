# DATA WRANGLING FOR 1_CHNG

# LIBS
library(ape)
source (file.path('tools', "data_tools.R"))

# DIRS
output.dir <- "0_data"
input.dir <- file.path(output.dir, "raw")

# DATA -- Mammal ()
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
p <- phylos[[1]]
livezy_mod <- matrix(data=NA, nrow=length(p$tip.label), ncol=ncol(livezy))
rownames(livezy_mod) <- p$tip.label
for(i in 1:nrow(livezy)) {
  mtchs <- which(grepl(rownames(livezy)[i], p$tip.label))
  if(length(mtchs) > 1) {
    for(j in mtchs) {
      livezy_mod[p$tip.label[j], ] <- livezy[i, ]
    }
  }
}
data <- list (phylos = phylos, chars = livezy_mod, children=children)
save (data, file = file.path (output.dir, "bird.RData"))