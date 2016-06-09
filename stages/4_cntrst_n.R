# PARSE NODE OBJ AND CALCULATE CONTRST N

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
input_file <- file.path('3_node_obj', 'res.RData')

# INPUT
load(input_file)

# REMOVE SUBSPECIES
cat("Removing subspecies ....")
txids <- ls(node_obj)
for(txid in txids) {
  rank <- node_obj[[txid]][['rank']]
  if(rank == "subspecies") {
    rm(list=txid, envir=node_obj)
  }
}
save(node_obj, file=file.path(output_dir, "res_1.RData"))
cat("Done. [", length(node_obj), "] nodes.\n", sep="")

# PARSE NAMES
cat("Parsing names ....")
blcklst <- NULL
txids <- ls(node_obj)
for(txid in txids) {
  # drop names with ignore pattern
  bool <- vector(length=length(ignore_pttrns))
  for(i in 1:length(ignore_pttrns)) {
    bool[i] <- any(grepl(ignore_pttrns[i], node_obj[[txid]][['nm']],
                         ignore.case=TRUE))
  }
  if(any(bool)) {
    blcklst <- c(blcklst, txid)
    rm(list=txid, envir=node_obj)
  }
}
unfnshd <- TRUE
while(unfnshd) {
  unfnshd <- FALSE
  # drop all descending names
  txids <- ls(node_obj)
  for(txid in txids) {
    prid <- node_obj[[txid]][['prid']]
    if(prid %in% blcklst) {
      unfnshd <- TRUE
      blcklst <- c(blcklst, txid)
      rm(list=txid, envir=node_obj)
    }
  }
}
save(node_obj, file=file.path(output_dir, "res_2.RData"))
cat("Done. [", length(node_obj), "] nodes.\n", sep="")

# ADD POST-NODE ID SLOT
cat("Adding post-node IDs to node object ....")
txids <- ls(node_obj)
for(txid in txids) {
  prid <- node_obj[[txid]][['prid']]
  if(prid %in% txids) {
    node_obj[[prid]][['ptid']] <-
      c(node_obj[[prid]][['ptid']], txid)
  }
}
save(node_obj, file=file.path(output_dir, "res_3.RData"))
cat('Done.\n') 

# REMOVE SINGLETONS
cat("Removing singleton nodes ....")
txids <- ls(node_obj)
for(txid in txids) {
  ptids <- node_obj[[txid]][['ptid']]
  if(length(ptids) == 1) {
    prid <- node_obj[[txid]][['prid']]
    node_obj[[ptids]][['prid']] <- prid
    if(prid %in% txids) {
      prid_ptids <- node_obj[[prid]][['ptid']]
      prid_ptids <- prid_ptids[prid_ptids != txid]
      prid_ptids <- c(prid_ptids, ptids)
      node_obj[[prid]][['ptid']] <- prid_ptids
    }
    rm(list=txid, envir=node_obj)
  }
}
save(node_obj, file=file.path(output_dir, "res_4.RData"))
cat("Done. [", length(node_obj), "] furcating nodes.\n", sep="")

# ADD KIDS TO NODES
cat("Adding kids to node object ....")
txids <- ls(node_obj)
spp <- which(unlist(lapply(txids, function(x) node_obj[[x]][['rank']] == "species")))
for(i in spp) {
  kid <- txids[i]
  node_obj[[kid]][['kids']] <- 'none'
  prid <- node_obj[[kid]][['prid']]
  addKid(prid, kid)
}
save(node_obj, file=file.path(output_dir, "res_5.RData"))
cat("Done. [", length(node_obj), "] nodes with more than 1 kids.\n", sep="")

# ADD SISTER ID SLOT
cat("Adding sister IDs to node object ....")
for(txid in txids) {
  prid <- node_obj[[txid]][["prid"]]
  if(prid %in% txids) {
    cnddts <- node_obj[[prid]][["ptid"]]
    node_obj[[txid]][['sstr']] <- cnddts[cnddts != txid]
  }
}
save(node_obj, file=file.path(output_dir, "res_6.RData"))
cat('Done.\n')

# ADD CONTRAST N ESTIMATES
cat("Calculating contrast N ....")
for(txid in txids) {
  sstrs <- node_obj[[txid]][['sstr']]
  if(!is.null(sstrs) && length(sstrs) > 0) {
    n <- length(node_obj[[txid]][['kids']])
    if(n > 0) {
      sstr_n <- 0
      # calculate based on the maximum possible sstr_n
      for(j in 1:length(sstrs)) {
        sstr_n <- sstr_n + length(node_obj[[sstrs[j]]][['kids']])
      }
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
    (length(sstrs) > max_nsstrs)
  if(bool) {
    cnddts[i] <- TRUE
  }
}
cnddts <- txids[cnddts]
cat("Done. [", length(cnddts), "] candidates nodes.\n", sep="")
save(node_obj, cnddts, file=file.path(output_dir, "res.RData"))

# TOP-10
cat("And the top 100 NCBI contrast N nodes are....\n")
cntrst_ns <- lapply(cnddts, function(x) node_obj[[x]][['cntrst_n']])
bool <- unlist(lapply(cntrst_ns, function(x) is.null(x)))
cntrst_ns <- unlist(cntrst_ns)
nms <- unlist(lapply(cnddts[!bool], function(x) node_obj[[x]][['nm']][['scientific name']]))
rnks <- unlist(lapply(cnddts[!bool], function(x) node_obj[[x]][['rank']]))
ordrd <- order(cntrst_ns)[1:100]
cc <- 1
for(i in ordrd) {
  nm <- paste0(nms[i], ' (', rnks[i], ')')
  cn <- cntrst_ns[i]
  spcr1 <- paste0(rep(' ', 3 - nchar(cc)), collapse="")
  spcr2 <- paste0(rep(' ', 40 - nchar(nm)), collapse="")
  cat(spcr1, cc, ' | ', nm, spcr2, signif(cn, 3), "\n")
  cc <- cc + 1
}

# END
cat(paste0('\nStage `contrast N` finished at [', Sys.time(), ']\n'))