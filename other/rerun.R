# RE-RUN SCRIPT FOR EDITING RESULTS OF STAGES WITHOUT RERUNNING ENTIRE STAGE

# START
cat(paste0('\nStage `contrast N` started at [', Sys.time(), ']\n'))

# PARAMETERS
source('parameters.R')

# FUNCTIONS
source(file.path('tools', 'node_obj_tools.R'))

# DIRS
output_dir <- "4_cntrst_n"
if(!file.exists(output_dir)) {
  dir.create(output_dir)
}

# INPUT
load(file=file.path(output_dir, "res_6.RData"))

# IDENTIFY AVOID IDS
cat("Identify and remove IDs to avoid ....")
txids <- ls(node_obj)
nms <- sapply(txids, function(x) node_obj[[x]][['nm']][['scientific name']])
avoid_bool <- rep(FALSE, length(nms))
for(pttrn in ignore_pttrns) {
  avoid_bool <- avoid_bool & grepl(pttrn, nms, ignore.case=TRUE)
}
txids <- txids[!avoid_bool]
avoid_ids <- txids[avoid_bool]
cat("Done.\n")

# REMOVE AVOID IDS
cat("Remove avoid IDs from sstr and kids ....")
for(avoid_id in avoid_ids) {
  rm(list=avoid_id, envir=node_obj)
}
for(txid in txids) {
  kids <- node_obj[[txid]][['kids']]
  sstrs <- node_obj[[txid]][['sstr']]
  node_obj[[txid]][['kids']] <- kids[!kids %in% avoid_ids]
  node_obj[[txid]][['sstr']] <- sstrs[!sstrs %in% avoid_ids]
}
cat("Done.\n")

# ADD CONTRAST N ESTIMATES
cat("Calculating contrast N ....")
txids <- ls(node_obj)
for(txid in txids) {
  sstrs <- node_obj[[txid]][['sstr']]
  if(!is.null(sstrs) && length(sstrs) > 0) {
    n <- length(node_obj[[txid]][['kids']])
    if(n > 0) {
      # calculate based on the maximum possible sstr_n
      sstr_n <- max(sapply(sstrs, function(x) length(node_obj[[x]][['kids']])))
      node_obj[[txid]][['sstr_n']] <- sstr_n
      node_obj[[txid]][['cntrst_n']] <- n/sstr_n
    }
  }
}
save(node_obj, file=file.path(output_dir, "res_7.RData"))
cat('Done.\n')

# IDENTIFYING LIKELY LFS BASED ON CNTRST_N AND PARENT SIZE
# this avoids searching timetree nodes that are likely to be young
cat("Identifying candidates ....")
txids <- ls(node_obj)
cnddts <- vector(length=length(txids))
for(i in 1:length(txids)) {
  cn <- node_obj[[txids[i]]][['cntrst_n']]
  prid <- node_obj[[txids[i]]][['prid']]
  sstrs <- node_obj[[txids[i]]][['sstr']]
  if(prid %in% txids) {
    prnt_n <- length(node_obj[[prid]][['kids']])
  }
  bool <- (!is.null(cn) && cn < max_cntrst_n) &
    (prnt_n > min_prnt_n | !prid %in% txids) &
    (length(sstrs) < max_nsstrs)
  if(bool) {
    cnddts[i] <- TRUE
  }
}
cnddts <- txids[cnddts]
cat("Done. [", length(cnddts), "] candidates nodes.\n", sep="")
save(node_obj, cnddts, file=file.path(output_dir, "res.RData"))

# END
cat(paste0('\nStage `contrast N` finished at [', Sys.time(), ']\n'))