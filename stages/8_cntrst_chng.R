# CALCULATE CONTRAST CHANGE

# START
cat(paste0('\nStage `contrast change` started at [', Sys.time(), ']\n'))

# PARAMETERS
cat("Loading parameters ...\n")
source('parameters.R')

# LIBS
cat("Loading libraries ...\n")
source(file.path('tools', 'node_obj_tools.R'))
source(file.path('tools', 'clade_matching_tools.R'))

# DIRS
chng_dir <- "2_chng"
tt_dir <- "6_timetree"
output_dir <- '7_cntrst_chng'
if (!file.exists(output_dir)) {
  dir.create(output_dir)
}
output_file <- file.path(output_dir, "res.RData")
input_file <- file.path(tt_dir, "res.RData")
load(input_file)

# CHANGE
chng_fls <- list.files(chng_dir)
for(chng_fl in chng_fls) {
  # INPUT
  grp <- sub("\\.RData", "", chng_fl)
  cat('    Working on [', grp, '] ....\n', sep="")
  txids <- ls(node_obj)
  txids <- getGrpTxids(txids, grp=grp)
  spp <- getSppTxids(txids)
  txids_wchng <- NULL
  load(file.path(chng_dir, chng_fl))
  
  # MATCH SPECIES NAMES
  cat("    Matching change estimates for species.... ")
  cp_spp <- unlist(lapply(clades_phylo[['clade.children']],
                          function(x) length(x) == 1))
  cp_spp <- which(cp_spp)
  cp_nms <- unlist(lapply(clades_phylo[['clade.children']][cp_spp],
                          function(x) gsub("_", " ", x[[1]])))
  for(i in 1:length(spp)) {
    txid <- spp[i]
    nm <- node_obj[[txid]][['nm']][['scientific name']]
    if(nm %in% cp_nms) {
      j <- cp_spp[[which(cp_nms == nm)]]
      node_obj[[txid]][['chng']] <- clades_phylo[['chng']][[j]]
      txids_wchng <- c(txids_wchng, txid)
    }
  }
  cat("Done.\n")
  
  # MATCH KIDS NAMES
  cat("    Matching change estimates for internal nodes.... ")
  nids <- txids[!txids %in% spp]
  cp_nids <- 1:length(clades_phylo[['clade.children']])
  cp_nids <- cp_nids[!cp_nids %in% cp_spp]
  cp_kids <- lapply(cp_nids, function(x) gsub("_", " ", clades_phylo[['clade.children']][[x]]))
  for(i in 1:length(nids)) {
    txid <- nids[i]
    kids <- node_obj[[txid]][['kids']]
    kids <- unlist(lapply(kids, function(x) node_obj[[x]][['nm']][['scientific name']]))
    mtch_scrs <- getMtchScrs(kids, cp_kids)
    if(max(mtch_scrs) > 1) {
      bst_mtch <- which.max(mtch_scrs)[1]
      chng <- clades_phylo[['chng']][[cp_nids[bst_mtch]]]
      node_obj[[txid]][["chng"]] <- chng
      txids_wchng <- c(txids_wchng, txid)
    }
  }
  cat("Done.\n")
  
  # CALCULATE CONTRASTS
  cat("    Calculating contrast change.... ")
  for(txid in txids_wchng) {
    sstrs <- node_obj[[txid]][['sstr']]
    sstr_chng <- NULL
    for(sstr in sstrs) {
      sstr_chng <- c(sstr_chng, node_obj[[txid]][['chng']])
    }
    sstr_chng <- mean(sstr_chng, na.rm=TRUE)
    node_obj[[txid]][['cntrst_chng']] <-
      node_obj[[txid]][['chng']] / sstr_chng
  }
  cat("Done.\n")
}

# OUTPUT
cat("Outputting ...\n")
save(node_obj, cnddts, file=output_file)

# END
cat(paste0('\nStage `contrast change` finished at [', Sys.time(), ']\n'))