# DATA WRANGLING FOR 1_CHNG

# START
cat(paste0('\nStage `data` started at [', Sys.time(), ']\n'))

# FUNCTIONS
library(ape)
source(file.path('tools', "data_tools.R"))

# DIRS
tree_dir <- file.path('0_data', "trees")
chars_dir <- file.path('0_data', "chars")

# MAMMALS
cat('Working on mammals ....\n')
data <- read.delim(file.path(chars_dir, "panTHERIA.txt"), na.strings = -999,
                    stringsAsFactors = FALSE)
pantheria <- data[ ,- c(1:5,36:55)]
for(i in 1:ncol(pantheria)) {
  temp.chars <- pantheria[!is.na(pantheria[ , i]),i]
  if(length(unique(temp.chars)) > 10) {
    pantheria[!is.na(pantheria[ , i]),i] <- cut(temp.chars, 10)
  }
}
rownames(pantheria) <- data[ ,"MSW93_Binomial"]
rownames(pantheria) <- sub(" ", "_", rownames(pantheria))
file <- "morpho_matrix_forR.nex"
oleary <- readNexusData(file.path(chars_dir, file))
chars <- merge(oleary, pantheria, by = 0, all = TRUE)
rownames(chars) <- c(rownames(oleary) [!rownames(oleary) %in% rownames(pantheria)], rownames(pantheria))
tree <- read.tree(file.path(tree_dir, "bininda_mammalia.tre"))
chars <- chars[rownames(chars) %in% tree$tip.label, ]
clades_phylo <- MoreTreeTools::getClades(tree)
data <- list(tree=tree, chars=chars, clades_phylo=clades_phylo)
save(data, file=file.path(chars_dir, "mammal.RData"))
prep <- signif(mean(colSums(!is.na(chars)))/length(tree$tip.label), 3)
cat('Done. Found [', ncol(chars), '] characters each on average representing [', 
    prep, '%] of all tips\n', sep="")

# BIRDS
cat('Working on birds ....\n')
file <- "X1228_Morphology Matrix_morphobank.nex"
livezy <- readNexusData(file.path(chars_dir, file))
trees <- read.tree(file.path(tree_dir, 'jetz_aves.tre'))
tree <- consensus(trees)  # strict consensus tree
tree$edge.length <- rep(1, nrow(tree$edge))
clades_phylo <- MoreTreeTools::getClades(tree)
# character data is for whole groups
# use character matching to assign the same value to memebers of the same group
# NAs for missing taxa
# use last tree in loop
livezy_mod <- matrix(data=NA, nrow=length(tree$tip.label), ncol=ncol(livezy))
rownames(livezy_mod) <- tree$tip.label
for(i in 1:nrow(livezy)) {
  mtchs <- which(grepl(rownames(livezy)[i], tree$tip.label))
  if(length(mtchs) > 1) {
    for(j in mtchs) {
      livezy_mod[tree$tip.label[j], ] <- livezy[i, ]
    }
  }
}
data <- list(tree=tree, chars=livezy_mod, clades_phylo=clades_phylo)
save(data, file = file.path(chars_dir, "bird.RData"))
prep <- signif(mean(colSums(!is.na(livezy_mod)))/length(tree$tip.label), 3)
cat('Done. Found [', ncol(livezy_mod), '] characters each on average representing [', 
    prep, '%] of all tips\n', sep="")

# END
cat(paste0('\nStage `data` finished at [', Sys.time(), ']\n'))