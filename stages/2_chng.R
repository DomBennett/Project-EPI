# ESTIMATE CHANGE PER NODE

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
cat("Looping through all published trees and character sets .... ")
ttl_cc <- 0
for(phychr_file in phychr_files) {
  cat('    Reading in [', phychr_file, '] ....', sep="")
  load(file.path(input_dir, phychr_file))
  tree <- data[["tree"]]
  chars <- data[["chars"]]
  clades_phylo <- data[['clades_phylo']]
  rm(data)
  cat("Done, found character matrix of [", nrow(chars),
      "] characters and [", ncol(chars), "] species\n")
  
  cat("    Reducing character matrix .... ")
  chars <- reduceChrctrMtrx(chars)
  nspp <- sum(colSums(!is.na(chars)) > 0)
  cat('Done, found [', ncol(chars), '] characters for [', 
      nspp, '/', length(tree$tip.label), '] tips\n', sep="")
  
  cat("    Estimating ancestral node states with parsimony .... ")
  reconstruction_obj <- parsimonyReconstruction(chars, tree)
  cat('Done\n')
  
  cat("    Calculating edge changes .... ")
  echanges <- calcChange(tree, reconstruction_obj, parallel=TRUE)
  cat('Done\n')
  
  cat("    Calculate mean change per clade .... ")
  clades_phylo <- calcMeanCladeChange(tree, clades_phylo, echanges)
  nds_wchrs <- sum(clades_phylo[['chng']] > 0, na.rm=TRUE)
  cat("Done, [", nds_wchrs, "] nodes with mean change.\n", sep="")
  
  cat("    Outputting ... ")
  save(clades_phylo, file=file.path(output_dir, phychr_file))
  cat("Done.\n")
}

# END
cat(paste0('\nStage `change` finished at [', Sys.time(), ']\n'))