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
phychr_files <- phychr_files[length(phychr_files):1]

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
  
  cat('    Preparing character matrix....')
  # do any have fewer than 10 species?
  chars <- chars[ ,colSums(!is.na(chars)) > 10]
  # convert non-numeric to numbers
  tmp_chars <- matrix(NA, nrow=nrow(chars), ncol=ncol(chars))
  rownames(tmp_chars) <- rownames(chars)
  colnames(tmp_chars) <- colnames(chars)
  for(i in 1:ncol(chars)) {
    if(any(tolower(chars[ ,i]) %in% letters)) {
      tmp <- tolower(chars[ ,i])
      chars[ ,i] <- match(tmp, letters)
      }
    tmp_chars[,i] <- as.numeric(chars[,i])
  }
  chars <- tmp_chars
  rm(tmp_chars)
  nspp <- sum(rowSums(!is.na(chars)) > 0)
  # do any have only 1 unique character?
  nstates <- apply(chars, 2, function(x) length(unique(x)))
  pull <- nstates > 1
  chars <- chars[ ,pull]
  nstates <- nstates[pull]
  cat('Done, found [', ncol(chars), '] characters for [', 
      nspp, '/', length(tree$tip.label), '] tips\n', sep="")
  
  cat("    Estimating ancestral node states with parsimony .... ")
  # requires categorical variables
  reconstruction_obj <- parsimonyReconstruction(chars, tree)
  char_labels <- sapply(reconstruction_obj, function(x) x[['character.name']])
  nstates <- nstates[char_labels]
  cat('Done\n')
  
  cat("    Calculating all pairwise R^2s .... ")
  # non-normal, categorical data, using spearman
  rsqs <- suppressWarnings(cor(chars, use="pairwise.complete.obs", method='spearman'))
  rsqs <- abs(rsqs)
  cat('Done\n')
  
  cat("    Calculating change scores .... ")
  changes <- calcChange(tree, reconstruction_obj, parallel=TRUE)
  cat('Done\n')
  
  cat("    Calculating clades_change obj .... ")
  clades_phylo <- MoreTreeTools::getClades(tree)
  changes_by_clade <- changesByClade(nds=clades_phylo[['clade.node']],
                                    changes, tree)
  # add labels and nstates
  for(i in 1:length(changes_by_clade)) {
    tmp <- data.frame(chng=changes_by_clade[[i]], nstate=nstates)
    changes_by_clade[[i]] <- tmp[!is.na(tmp[['chng']]), ]
  }
  clades_phylo$changes <- changes_by_clade
  clades_change <- clades_phylo
  cat('Done\n')
  
  cat("    Outputting ... ")
  save(clades_change, rsqs, file=file.path(output_dir, phychr_file))
  cat("Done.\n")
}

# END
cat(paste0('\nStage `change` finished at [', Sys.time(), ']\n'))