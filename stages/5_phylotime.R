# GET NODE TIMINGS WITH PUBLISHED PHYLOGENIES
# TODO: use analysis groups to reduce looping

# START
cat(paste0('\nStage `phylotime` started at [', Sys.time(), ']\n'))

# PARAMETERS
source('parameters.R')

# FUNCTIONS
library(treeman)
source(file.path('tools', 'phylotime_tools.R'))
source(file.path('tools', 'node_obj_tools.R'))

# DIRS
if(!file.exists('5_phylotime')) {
  dir.create('5_phylotime')
}
tree_dir <- file.path("0_data", "trees")
input_file <- file.path("4_cntrst_n", "res.RData")
output_file <- file.path('5_phylotime', 'res.RData')
dst_mtrx_dir <- file.path('1_wrngl')

# INPUT
load(input_file)
tree_files <- list.files(tree_dir)

# LOOP THROUGH TREE FILES
cat("Looping through all published trees ....\n")
ttl_cc <- 0
top_cnddts <- cnddts
for(tree_file in tree_files) {
  # INPUT
  grp <- sub("\\.tre", "", tree_file)
  cat('    Reading in [', grp, '] ....\n', sep="")
  txids <- ls(node_obj)
  txids <- getGrpTxids(txids, grp=grp)
  cnddts <- c(txids, cnddts)
  spp <- getSppTxids(txids)
  tree <- readTree(file.path(tree_dir, tree_file))
  map_obj <- list()
  cat("    Done.\n")
  
  # DISTANCE MATRIX
  cat("    Distance matrix .... ")
  # run this here while the tree is read in
  dst_mtrx <- calcDstMtrx(tree, tree['all'], .parallel=parallel)
  save(dst_mtrx, file=file.path(dst_mtrx_dir, paste0(grp, ".RData")))
  cat("Done.\n")
  
  # MATCH TIPS TO NMS
  cat("    Matching tree tips to named species in node_obj .... ")
  cc <- 0
  mtch_tps <- gsub("_", " ", tree['tips'])
  tps <- tree['tips']
  for(sp in spp) {
    bool <- node_obj[[sp]][['nm']] %in% mtch_tps
    if(any(bool)) {
      map_obj[[sp]][['nid']] <- tps[mtch_tps == node_obj[[sp]][['nm']][bool][1]]
      cc <- cc + 1
    }
  }
  cat("    Done, matched [", cc, "] nodes to tips in tree.\n", sep="")
  
  # MATCH NODES OF INTEREST TO NODES IN TREE
  cat("    Matching nodes in node_obj to tree ....")
  cc <- 0
  nids <- tree['nds']
  tree_kids <- getNdsKids(tree, ids=nids)
  for(txid in txids) {
    nd_kids <- node_obj[[txid]][['kids']]
    if(nd_kids[1] != "none") {
      nd_tips <- vector(length=length(nd_kids))
      for(j in 1:length(nd_kids)) {
        if(!is.null(map_obj[[nd_kids[j]]])) {
          nd_tips[j] <- map_obj[[nd_kids[j]]][['nid']]
        }
      }
      nd_tips <- nd_tips[nd_tips != FALSE]
      if(length(nd_tips) == 0) {
        next
      }
      if(length(nd_tips) == 1) {
        map_obj[[txid]][["nid"]] <- nd_tips
        cc <- cc + 1
      } else {
        mtch_scrs <- vector(length=length(tree_kids))
        for(j in 1:length(tree_kids)) {
          mtch_scrs[j] <- (sum(nd_tips %in% tree_kids[[j]])/length(nd_tips)) +
            (sum(tree_kids[[j]] %in% nd_tips)/length(tree_kids[[j]]))
        }
        if(max(mtch_scrs) > 1) {
          bst_mtch <- which.max(mtch_scrs)[1]
          map_obj[[txid]][["nid"]] <- nids[bst_mtch]
          cc <- cc + 1
        }
      }
    }
  }
  cat(" Done, matched [", cc, "] nodes to tree.\n", sep="")
  
  # FIND PD, ED, PE AND AGE FOR EVERY NODE
  cat("    Adding tree timings to node_obj .... ")
  ed_vals <- calcFrPrp(tree, tids=tree['tips'], .parallel=TRUE)
  cc <- 0
  mtxids <- names(map_obj)
  for(txid in mtxids) {
    nid <- map_obj[[txid]][['nid']]
    if(nid == "n1") {  # this is the root
      next
    }
    sid <- getNdSstr(tree, id=nid)[1]
    age <- getNdAge(tree, id=nid)
    kids <- getNdKids(tree, nid)
    ed <- mean(ed_vals[which(tps %in% kids)])
    kids <- getNdKids(tree, sid)
    sstr_ed <- mean(ed_vals[which(tps %in% kids)])
    spn <- getNdSlt(tree, slt_nm="spn", id=nid)
    sstr_spn <- getNdSlt(tree, slt_nm="spn", id=sid)
    pd <- getNdSlt(tree, slt_nm="pd", id=nid)
    sstr_pd <- getNdSlt(tree, slt_nm="pd", id=sid)
    node_obj[[txid]][['ed']] <- ed
    node_obj[[txid]][['cntrst_ed']] <- ed/sstr_ed
    node_obj[[txid]][['pe']] <- spn
    node_obj[[txid]][['cntrst_pe']] <- (spn + 1)/(sstr_spn + 1)
    node_obj[[txid]][['pd']] <- pd
    node_obj[[txid]][['cntrst_pd1']] <- pd/sstr_pd
    node_obj[[txid]][['cntrst_pd2']] <- pd - sstr_pd
    node_obj[[txid]][['age']] <- age
    node_obj[[txid]][['tmsplt']] <- spn + age
    cc <- cc + 1
    ttl_cc <- ttl_cc + 1
  }
  cat("Done, got timings for [", cc, "] nodes from tree.\n", sep="")
}
cat("Done. Got timings for [", ttl_cc, "] nodes from all trees.\n", sep="")

# OUTPUT
cat('Saving .... ')
save(node_obj, top_cnddts, cnddts, file=output_file)
cat('Done.\n')

# END
cat(paste0('\nStage `phylotime` finished at [', Sys.time(), ']\n'))