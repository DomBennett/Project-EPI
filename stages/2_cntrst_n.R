# PARSE NODE OBJ AND CALCULATE CONTRST N

# START
cat(paste0('\nStage `cntrst n` started at [', Sys.time(), ']\n'))

# PARAMETERS
source('parameters.R')

# FUNCTIONS
source(file.path('tools', 'node_obj_tools.R'))

# DIRS
if(!file.exists("2_cntrst_n")) {
  dir.create("2_cntrst_n")
}
output_file <- file.path("2_cntrst_n", "res.RData")
input_file <- file.path('1_node_obj', 'res.RData')

# INPUT
load(input_file)

# REMOVE SUBSPECIES
cat("Removing subspecies ....\n")
txids <- ls(node_obj)
for(txid in txids) {
  rank <- node_obj[[txid]][['rank']]
  if(rank == "subspecies") {
    rm(list=txid, envir=node_obj)
  }
}
cat("Done. [", length(node_obj), "] nodes.\n", sep="")

# PARSE NAMES
cat("Parsing names ....\n")
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
cat("Done. [", length(node_obj), "] nodes.\n", sep="")

# ADD POST-NODE ID SLOT
cat("Adding post-node IDs to node object ....\n")
txids <- ls(node_obj)
for(txid in txids) {
  prid <- node_obj[[txid]][['prid']]
  if(prid %in% txids) {
    node_obj[[prid]][['ptid']] <-
      c(node_obj[[prid]][['ptid']], txid)
  }
}
cat('Done.\n') 

# REMOVE SINGLETONS
cat("Removing singleton nodes ....\n")
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
cat("Done. [", length(node_obj), "] furcating nodes.\n", sep="")

# ADD KIDS TO NODES
cat("Adding kids to node object ....\n")
txids <- ls(node_obj)
spp <- which(unlist(lapply(txids, function(x) node_obj[[x]][['rank']] == "species")))
for(i in spp) {
  kid <- txids[i]
  node_obj[[kid]][['kids']] <- 'none'
  prid <- node_obj[[kid]][['prid']]
  addKid(prid, kid)
}
cat("Done. [", length(node_obj), "] nodes with more than 1 kids.\n", sep="")

# ADD SISTER ID SLOT
cat("Adding sister IDs to node object ....\n")
for(txid in txids) {
  prid <- node_obj[[txid]][["prid"]]
  if(prid %in% txids) {
    cnddts <- node_obj[[prid]][["ptid"]]
    node_obj[[txid]][['sstr']] <- cnddts[cnddts != txid]
  }
}
cat('Done.\n')

# ADD CONTRAST N ESTIMATES
cat("Calculating contrast N ....\n")
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
cat('Done.\n')

# TOP-10
cat("And the top 100 NCBI contrast N nodes are....\n")
cntrst_ns <- lapply(txids, function(x) node_obj[[x]][['cntrst_n']])
bool <- unlist(lapply(cntrst_ns, function(x) is.null(x)))
cntrst_ns <- unlist(cntrst_ns)
nms <- unlist(lapply(txids[!bool], function(x) node_obj[[x]][['nm']][['scientific name']]))
rnks <- unlist(lapply(txids[!bool], function(x) node_obj[[x]][['rank']]))
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

# OUTPUT
cat('Saving ....\n')
save(node_obj, file=output_file)
cat('Done.\n')

# END
cat(paste0('\nStage `cntrst n` finished at [', Sys.time(), ']\n'))