# ESTIMATE CHANGE PER NODE
# TODO: where do NaNs originate?

# START
cat(paste0('\nStage `change` started at [', Sys.time(), ']\n'))

# PARAMETERS
cat("Loading parameters ...\n")
source('parameters.R')

# LIBS
cat("Loading functions ...\n")
library(ape)
source(file.path("tools", "chng_tools.R"))

# DIRS
input_dir <- file.path("1_wrngl")
output_dir <- "2_chng"
if(!file.exists(output_dir)) {
  dir.create(output_dir)
}

# INPUT
phychr_files <- list.files(input_dir, pattern=".RData")

# LOOP THROUGH TREE FILES
cat("Looping through all published trees and character sets ....\n")
ttl_cc <- 0
for(phychr_file in phychr_files) {
  cat('    Reading in [', phychr_file, '] ....', sep="")
  load(file.path(input_dir, phychr_file))
  tree <- data[["tree"]]
  chars <- data[["chars"]]
  rm(data)
  nspp <- sum(rowSums(!is.na(chars)) > 0)
  cat("Done, found character matrix of [", ncol(chars),
      "] characters and [", nspp, '/', length(tree$tip.label),
      '] tips\n', sep="")
  
  cat('Preparing character matrix....')
  # do any have fewer than 10 species?
  chars <- chars[ ,colSums(!is.na(chars)) > 10]
  # do any have only 1 unique character?
  nuchars <- apply(chars, 2, function(x) length(unique(x)))
  pull <- nuchars > 1
  chars <- chars[ ,pull]
  nspp <- sum(rowSums(!is.na(chars)) > 0)
  cat('Done, found [', ncol(chars), '] characters for [', 
      nspp, '/', length(tree$tip.label), '] tips\n', sep="")
  
  cat("    Estimating ancestral node states with parsimony .... ")
  # requires categorical variables
  reconstruction_obj <- parsimonyReconstruction(chars, tree)
  cat('Done\n')
  
  cat("    Calculating n. changes .... ")
  changes <- calcChange(tree, reconstruction_obj, parallel=TRUE)
  cat('Done\n')
  
  cat("    Calculating clades_change obj .... ")
  clades_phylo <- MoreTreeTools::getClades(tree)
  changes_by_node <- changesByNode(nds=clades_phylo[['clade.node']],
                                   changes, parallel=TRUE)
  clades_phylo$changes <- changes_by_node
  clades_change <- clades_phylo
  
  cat("    Outputting ... ")
  save(clades_change, file=file.path(output_dir, phychr_file))
  cat("Done.\n")
}

# END
cat(paste0('\nStage `change` finished at [', Sys.time(), ']\n'))