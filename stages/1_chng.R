# ESTIMATE CHANGE PER NODE

# START
cat(paste0('\nStage `chng` started at [', Sys.time(), ']\n'))

# PARAMETERS
cat("Loading parameters ...\n")
source('parameters.R')

# LIBS
cat("Loading functions ...\n")
library(ape)
source(file.path("tools", "chng_tools.R"))

# DIRS
input_dir <- file.path('0_data', "chars")
output_dir <- "1_change"
if(!file.exists(output_dir)) {
  dir.create(output_dir)
}

# INPUT
phychr_files <- list.files(input_dir, pattern=".RData")

# LOOP THROUGH TREE FILES
cat("Looping through all published trees and character sets ....\n")
ttl_cc <- 0
for(phychr_file in phychr_files) {
  cat('    Reading in [', phychr_file, '] ....\n', sep="")
  load(file.path(input_dir, phychr_file))
  tree <- data[["tree"]]
  chars <- data[["chars"]]
  clades_phylo <- data[['clades_phylo']]
  rm(data)
  cat("    Done.\n")
  
  cat("    Reducing character matrix ....\n")
  chars <- reduceChrctrMtrx(chars, pcut=1)
  prep <- signif(mean(colSums(!is.na(chars)))/length(tree$tip.label), 3)
  cat('    Done. Found [', ncol(chars), '] characters each on average representing [', 
      prep, '%] of all tips\n', sep="")
  
  cat("    Estimating ancestral node states with parsimony ....\n")
  reconstruction_obj <- parsimonyReconstruction(chars, tree)
  cat('    Done\n')
  
  cat("    Calculating edge changes ....\n")
  echanges <- calcChange(tree, reconstruction_obj)
  cat('Done\n')
  
  cat("    Calculate mean change per clade ....\n")
  clades_phylo <- calcMeanCladeChange(tree, clades_phylo, echanges)
  nds_wchrs <- sum(clades_phylo[['chng']] > 0, na.rm=TRUE)
  cat("Done. [", nds_wchrs, "] nodes with mean change.\n", sep="")
  
  cat("    Outputting ...\n")
  save(clades_phylo, file=file.path(output_dir, phychr_file))
  cat("    Done.\n")
}

# END
cat(paste0('\nStage `chng` finished at [', Sys.time(), ']\n'))