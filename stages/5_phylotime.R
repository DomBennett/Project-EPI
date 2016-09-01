# GET NODE TIMINGS WITH PUBLISHED PHYLOGENIES

# START
cat(paste0('\nStage `phylotime` started at [', Sys.time(), ']\n'))

# PARAMETERS
source('parameters.R')

# FUNCTIONS
library(treeman)
source(file.path('tools', 'phylotime_tools.R'))
source(file.path('tools', 'node_obj_tools.R'))
source(file.path('tools', 'clade_matching_tools.R'))
source(file.path('tools', 'i_tools.R'))

# DIRS
outdir <- '5_phylotime'
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
ph_obj <- new.env()
for(txid in ph_txids) {
  ph_obj[[txid]] <- node_obj[[txid]]
}
rm(list=ph_txids, envir=node_obj)
save(node_obj, file=file.path(outdir, 'node_obj.RData'))
rm(node_obj)
cat('Done.\n')

# LOOP THROUGH TREE FILES
cat("Looping through all published trees ....\n")
for(tree_file in tree_files) {
  # INPUT
  grp <- sub("\\.RData", "", tree_file)
  cat('    Working on [', grp, '] ....\n', sep="")
  load(file.path(outdir, tree_file))
  spp <- getSppTxids(txids)
  ntrees <- getNtrees(tree_file)
  
  # MAIN LOOP
  cat("    Main loop....\n")
  for(i in 1:ntrees) {
    iPrnt(i, ntrees)
    tree <- getTree(i, tree_file)
    
    # MATCH TIPS TO NMS
    map_obj <- list()
    mtch_tps <- gsub("_", " ", tree['tips'])
    tps <- tree['tips']
    for(sp in spp) {
      bool <- prt_obj[[sp]][['nm']] %in% mtch_tps
      if(any(bool)) {
        map_obj[[sp]][['nid']] <- tps[mtch_tps == prt_obj[[sp]][['nm']][bool][1]]
      }
    }
    
    # MATCH NODES OF INTEREST TO NODES IN TREE
    nids <- tree['nds']
    tree_kids <- getNdsKids(tree, ids=nids)
    for(txid in txids) {
      nd_kids <- prt_obj[[txid]][['kids']]
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
        } else {
          mtch_scrs <- getMtchScrs(nd_kids, nd_tips)
          if(max(mtch_scrs) > 1) {
            bst_mtch <- which.max(mtch_scrs)[1]
            map_obj[[txid]][["nid"]] <- nids[bst_mtch]
          }
        }
      }
    }
    
    # FIND PD, ED, PE AND AGE FOR EVERY NODE
    ed_vals <- calcFrPrp2(tree, tids=tree['tips'])
    tree_age <- tree['age']
    mtxids <- names(map_obj)
    for(txid in mtxids) {
      nid <- map_obj[[txid]][['nid']]
      if(nid == tree['root']) {  # this is the root
        next
      }
      sid <- getNdSstr(tree, id=nid)[1]
      age <- getNdAge(tree, id=nid, tree_age)
      kids <- getNdKids(tree, nid)
      ed <- mean(ed_vals[which(tps %in% kids)])
      kids <- getNdKids(tree, sid)
      sstr_ed <- mean(ed_vals[which(tps %in% kids)])
      spn <- getNdSlt(tree, slt_nm="spn", id=nid)
      sstr_spn <- getNdSlt(tree, slt_nm="spn", id=sid)
      pd <- getNdPD(tree, id=nid)
      sstr_pd <- getNdSlt(tree, slt_nm="pd", id=sid)
      assgnWMean(val=ed, nm="ed")
      assgnWMean(val=ed/sstr_ed, nm="cntrst_ed")
      assgnWMean(val=spn, nm="pe")
      assgnWMean(val=(spn + 1)/(sstr_spn + 1), nm="cntrst_pe")
      assgnWMean(val=pd, nm="pd")
      assgnWMean(val=pd/sstr_pd, nm="cntrst_pd1")
      assgnWMean(val=pd - sstr_pd, nm="cntrst_pd2")
      assgnWMean(val=age, nm="age")
      assgnWMean(val=spn + age, nm="tmsplt")
    }
  }
}
cat("Done.\n")

# OUTPUT
cat('Saving .... ')
save(ph_obj, file=file.path(outdir, 'ph_obj.RData'))
cat('Done.\n')

# END
cat(paste0('\nStage `phylotime` finished at [', Sys.time(), ']\n'))