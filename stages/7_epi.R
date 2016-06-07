# CALCULATE EPIS FROM NODE_OBJ
# TODO: separate this and make a separate stage for contrast change

# START
cat(paste0('\nStage `epi` started at [', Sys.time(), ']\n'))

# PARAMETERS
cat("Loading parameters ...\n")
source('parameters.R')

# LIBS
cat("Loading libraries ...\n")
source(file.path("tools", "epi_tools.R"))
source(file.path('tools', 'node_obj_tools.R'))

# DIRS
chng_dir <- "2_chng"
tt_dir <- "6_timetree"
output_dir <- '7_epi'
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
    }
  }
  cat("Done.\n")
  
  # MATCH KIDS NAMES
  cat("    Matching change estimates for internal nodes.... ")
  nids <- txids[!txids %in% spp]
  cp_nids <- 1:length(clades_phylo[['clade.children']])
  cp_nids <- cp_nids[!cp_nids %in% cp_spp]
  for(i in 1:length(nids)) {
    txid <- nids[i]
    kids <- node_obj[[txid]][['kids']]
    kids <- unlist(lapply(kids, function(x) node_obj[[x]][['nm']][['scientific name']]))
    mtch_scrs <- vector(length=length(cp_nids))
    for(j in 1:length(cp_nids)) {
      cp_kids <- clades_phylo[['clade.children']][[cp_nids[j]]]
      cp_kids <- gsub("_", " ", cp_kids)
      mtch_scrs[i] <- (sum(kids %in% cp_kids)/length(kids)) +
        (sum(cp_kids %in% kids)/length(cp_kids))
    }
    if(max(mtch_scrs) > 1) {
      bst_mtch <- which.max(mtch_scrs)[1]
      chng <- clades_phylo[['chng']][[cp_nids[bst_mtch]]]
      node_obj[[txid]][["chng"]] <- chng
    }
  }
  cat("Done.\n")
}

# GENERATE EPI DATAFRAME
cat("Generating EPI dataframe.... ")
txids <- ls(node_obj)
nms <- tmsplts <- cntrst_ns <- chngs <- rep(NA, length(txids))
for(i in 1:length(txids)) {
  tmsplt <- node_obj[[txids[i]]][["tmsplt"]]
  cntrst_n <- node_obj[[txids[i]]][["cntrst_n"]]
  if(!is.null(tmsplt) & !is.null(cntrst_n)) {
    tmsplts[i] <- tmsplt
    cntrst_ns[i] <- cntrst_n
    nms[i] <- node_obj[[txids[i]]][["nm"]][["scientific name"]]
  }
  if(!is.null(node_obj[[txids[i]]][["chng"]])) {
    chngs[i] <- node_obj[[txids[i]]][["chng"]]
  }
}
bool <- !is.na(cntrst_ns)
cntrst_ns <- cntrst_ns[bool]
tmsplts <- tmsplts[bool]
txids <- txids[bool]
nms <- nms[bool]
chngs <- chngs[bool]
epi <- data.frame(nm=nms, txid=txids, tmsplt=tmsplts, cntrst_n=cntrst_ns,
                  chng=chngs, stringsAsFactors=FALSE)
cat("Done.\n")

# CALCULATE EPIS
cat("Calculating EPIs .... ")
epi$epi <- (log(epi[['cntrst_chng']]) + log(epi[['cntrst_n']])) /
  log(epi[['tmsplt']])
epi$pepi <- log(epi[['cntrst_n']]) / log(epi[['tmsplt']])
cat("Done.\n")

# OUTPUT
cat("Outputting ...\n")
save(node_obj, epi, file=output_file)

# END
cat(paste0('\nStage `epi` finished at [', Sys.time(), ']\n'))