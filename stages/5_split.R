# SPLIT NODE OBJ

# START
cat(paste0('\nStage `split` started at [', Sys.time(), ']\n'))

# FUNCTIONS
library(treeman)
source(file.path('tools', 'node_obj_tools.R'))

# DIRS
outdir <- '5_split'
if(!file.exists(outdir)) {
  dir.create(outdir)
}
tree_dir <- file.path("0_data", "trees")
input_file <- file.path("4_cntrst_n", "res.RData")

# INPUT
load(input_file)
tree_files <- list.files(tree_dir, pattern='.RData')

# SPLIT NODE_OBJ
# split node_obj into one with phylo estiamtes and one without
cat("Splitting node_obj to save memory ....\n")
txids <- ls(node_obj)
ph_txids <- NULL
for(tree_file in tree_files) {
  grp <- sub("\\.RData", "", tree_file)
  ph_txids <- c(ph_txids, getGrpTxids(txids, grp=grp))
}
cat(".... extracting phylo IDs\n")
ph_obj <- new.env()
for(txid in ph_txids) {
  ph_obj[[txid]] <- node_obj[[txid]]
}
cat(".... outputting\n")
rm(list=ph_txids, envir=node_obj)
tt_obj <- node_obj
save(tt_obj, file=file.path(outdir, 'tt_obj.RData'))
save(ph_obj, file=file.path(outdir, 'ph_obj.RData'))
cat('Done.\n')

# END
cat(paste0('\nStage `split` finished at [', Sys.time(), ']\n'))