# CALCULATE CONTRAST CHANGE

# START
cat(paste0('\nStage `contrast change` started at [', Sys.time(), ']\n'))

# PARAMETERS
cat("Loading parameters ...\n")
source('parameters.R')
min_nchars <- 4

# LIBS
cat("Loading libraries ...\n")
source(file.path('tools', 'node_obj_tools.R'))
source(file.path('tools', 'clade_matching_tools.R'))

# DIRS
chng_dir <- "2_chng"
ph_dir <- "6_phylotime"
output_dir <- '8_cntrst_chng'
if (!file.exists(output_dir)) {
  dir.create(output_dir)
}
output_file <- file.path(output_dir, "res.RData")
input_file <- file.path(ph_dir, "ph_obj.RData")
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
  rsqs <- abs(rsqs)
  
  # MATCH
  cat("    Matching change estimates nodes in node_obj.... ")
  cc_nids <- 1:length(clades_change[['clade.children']])
  cc_kids <- lapply(cc_nids, function(x) gsub("_", " ",
                                              clades_change[['clade.children']][[x]]))
  cc_spp <- which(sapply(cc_kids, length) == 1)
  cc_spp_nms <- unlist(cc_kids[cc_spp])
  mtch_indx <- rep(NA, length(cc_nids))
  for(i in 1:length(txids)) {
    txid <- txids[i]
    kids <- node_obj[[txid]][['kids']]
    if(length(kids) == 1 && kids == 'none') {
      nm <- node_obj[[txid]][['nm']][['scientific name']]
      if(nm %in% cc_spp_nms) {
        mtch_indx[i] <- cc_spp[cc_spp_nms == nm]
      }
    } else {
      kids <- unlist(lapply(kids, function(x) node_obj[[x]][['nm']][['scientific name']]))
      mtch_scrs <- getMtchScrs(kids, cc_kids)
      if(max(mtch_scrs) > 1) {
        mtch_indx[i] <- which.max(mtch_scrs)[1]
      }
    }
  }
  cat("Done. Found [", sum(!is.na(mtch_indx)), "/", length(mtch_indx),
      '] of clade change estimates in node_obj.\n')
  
  # CALCULATE CNTRST CHNG
  cat("    Calculating contrast change.... ")
  cntr <- 0
  txids <- txids[!is.na(mtch_indx)]
  mtch_indx <- mtch_indx[!is.na(mtch_indx)]
  for(i in 1:length(txids)) {
    txid <- txids[i]
    sstr <- node_obj[[txid]][['sstr']]
    if(!any(sstr %in% txids)) {
      next
    }
    sstr_i <- which(txids %in% sstr)
    chngs <- clades_change[['changes']][[mtch_indx[i]]]
    if(nrow(chngs) < min_nchars) {
      next
    }
    chngs$sstr <- NA
    char_labels <- rownames(chngs)
    if(length(sstr_i) == 1) {
      sstr_chngs <- clades_change[['changes']][[mtch_indx[sstr_i]]]
      chngs$sstr <- sstr_chngs[['chng']][match(char_labels, rownames(sstr_chngs))]
    } else {
      # unpack
      sstr_chngs_list <- lapply(sstr_i,
                                function(x) clades_change[['changes']][[mtch_indx[x]]])
      # create list of vectors
      sstr_chngs <- vector("list", length=length(sstr_chngs_list))
      for(j in 1:length(sstr_chngs_list)) {
        tmp_chng_data <- sstr_chngs_list[[j]]
        sstr_chngs[[j]] <-  tmp_chng_data[['chng']][match(char_labels, rownames(tmp_chng_data))]
      }
      # create vector of means
      sstr_chng <- rep(NA, length(char_labels))
      for(j in 1:length(char_labels)) {
        sstr_vals <- sapply(sstr_chngs, function(x) x[j])
        if(all(is.na(sstr_vals))) {
          sstr_chng[j] <- NA
        } else {
          sstr_chng[j]<- max(sstr_vals, na.rm=TRUE)
        }
      }
      chngs$sstr <- sstr_chng
    }
    chngs <- chngs[!is.na(chngs$sstr), ]
    if(nrow(chngs) < min_nchars) {
      next
    }
    chngs$rsq <- rowMeans(rsqs[rownames(chngs), ], na.rm=TRUE)
    chngs$rsq[is.na(chngs$rsq)] <- 0
    chngs$cntrst_chng <- chngs$chng/chngs$sstr
    # weighted mean cntrst chng
    wts <- (1/chngs$nstate)*(1-chngs$rsq)
    cntrst_chng <- weighted.mean(chngs$cntrst_chng, w=wts, na.rm=TRUE)
    cntr <- cntr + 1
    node_obj[[txid]][['cntrst_chng']] <- list('cntrst_chng'=cntrst_chng,
                                              'chng_data'=chngs)
  }
  cat("Done. Calculated change for [", cntr, '/', length(txids),
      '] clades.\n')
}

# OUTPUT
cat("Outputting ...\n")
save(node_obj, file=output_file)

# END
cat(paste0('\nStage `contrast change` finished at [', Sys.time(), ']\n'))