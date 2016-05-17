# PARSE NCBI TAXONOMY DMP
# ftp://ftp.ncbi.nih.gov/pub/taxonomy/ (downloaded: 16/05/2016)

# START
cat(paste0('\nStage `NCBI taxonomy` started at [', Sys.time(), ']\n'))

# FUNCTIONS
source(file.path('tools', 'ncbi_taxonomy_tools.R'))

# PARAMETERS
source('parameters.R')

# DIRS
if(!file.exists("ncbi_cache")) {
  dir.create("ncbi_cache")
}
if(!file.exists("x_ncbi_taxonomy")) {
  dir.create("x_ncbi_taxonomy")
}
output_file <- file.path("x_ncbi_taxonomy", "res.RData")
input_dir <- file.path('0_data', 'raw')

# INPUT
cat("Reading in nodes and names dump ....\n")
nodes_lines <- readLines(con=file.path(input_dir, 'nodes.dmp'), n=-1L)
nodes_lines <- strsplit(nodes_lines, split="\\|")
names_lines <- readLines(con=file.path(input_dir, 'names.dmp'), n=-1L)
names_lines <- strsplit(names_lines, split="\\|")
cat("Done.\n")

# CONVERT TO NODE-OBJ
cat("Parsing and generating node object ....\n")
node_obj <- vector("list", length=length(nodes_lines))
for(i in 1:length(node_obj)) {
  division_code <- as.numeric(nodes_lines[[i]][5])
  if(division_code %in% division_codes) {
    txid <- as.numeric(nodes_lines[[i]][1])
    prid <- as.numeric(nodes_lines[[i]][2])
    rank <- gsub('\t', '', nodes_lines[[i]][3])
    node_obj[[i]] <- list(txid=txid, prid=prid, ptid=NULL,
                          rank=rank, kids=NULL)
  }
}
fnds <- unlist(lapply(node_obj, function(x) !is.null(x)))
node_obj <- node_obj[fnds]
sbspp <- unlist(lapply(node_obj, function(x) x[['rank']] == "subspecies"))
node_obj <- node_obj[!sbspp]
txids <- unlist(lapply(node_obj, function(x) x[['txid']]))
cat("Done. Identified [", length(node_obj), "] valid taxonomic nodes.\n", sep="")

# PARSE NAMES
cat("Parsing names ....\n")
names_lines <- dropNamesLines(names_lines, txids, parallel=TRUE)
to_drp <- vector(length=length(txids))
nms_obj <- vector('list', length=length(txids))
names(nms_obj) <- txids
for(i in 1:length(names_lines)) {
  txid <- as.numeric(names_lines[[i]][1])
  if(txid %in% txids) {
    j <- which(txids == txid)
    nm <- gsub('\t', '', names_lines[[i]][2])
    to_drp[j] <- grepl("unclassified", nm)[[1]]
    nm_typ <- gsub('\t', '', names_lines[[i]][4])
    nms <- nms_obj[[which(txids == txid)]]
    nms_obj[[j]] <- c(nms, nm)
    names(nms_obj[[j]]) <- c(names(nms), nm_typ)
  }
}
nms_obj <- nms_obj[!to_drp]
node_obj <- node_obj[!to_drp]
txids <- txids[!to_drp]
cat("Done. [", length(node_obj), "] nodes with valid taxonomic names.\n", sep="")

# ADD KIDS TO NODES
cat("Adding kids to node object ....\n")
spp <- which(unlist(lapply(node_obj, function(x) x[['rank']] == "species")))
for(i in spp) {
  kid <- txids[i]
  node_obj[[i]][['kids']] <- 'none'
  prid <- node_obj[[i]][['prid']]
  node_obj <- addKid(node_obj, prid, kid)
}
cat('Done.\n')

# ADD PRE-NODE ID SLOT
cat("Adding pre-node IDs to node object ....\n")
for(i in 1:length(node_obj)) {
  prid <- node_obj[[i]][['prid']]
  prid_i <- which(txids == prid)
  if(length(prid_i) > 0) {
    node_obj[[prid_i]][['ptid']] <-
      c(node_obj[[prid_i]][['ptid']], node_obj[[i]][['txid']])
  }
}
cat('Done.\n')

# ADD SISTER ID SLOT
cat("Adding sister IDs to node object ....\n")
for(i in 1:length(node_obj)) {
  txid <- node_obj[[i]][['txid']]
  prid <- node_obj[[i]][['prid']]
  prid_i <- which(txids == prid)
  if(length(prid_i) > 0) {
    ptids <- node_obj[[prid_i]][['ptid']]
    sstr <- ptids[ptids != txid]
    node_obj[[i]][['sstr']] <- sstr
  }
}
cat('Done.\n')

# ADD CONTRAST N ESTIMATES
cat("Calculating contrast N ....\n")
for(i in 1:length(node_obj)) {
  sstrs <- node_obj[[i]][['sstr']]
  if(!is.null(sstrs) && length(sstrs) > 0) {
    cns <- vector(length=length(sstrs))
    names(cns) <- sstrs
    n <- length(node_obj[[i]][['kids']])
    if(n > 0) {
      for(j in 1:length(sstrs)) {
        si <- which(txids == sstrs[j])
        sstr_n <- length(node_obj[[si]][['kids']])
        cns[j] <- n/sstr_n
      }
      node_obj[[i]][['cntrst_n']] <- cns
    }
  }
}
cat('Done.\n')

# IGNORE BOOL
cat("Determining which nodes to ignore in later analysis ....\n")
ignore_bool <- rep(FALSE, length(node_obj))
for(i in 1:length(node_obj)) {
  cns <- node_obj[[i]][["cntrst_n"]]
  bool <- any(cns > contrst_n_min | cns < 1/contrst_n_min)
  if(is.na(bool) || !bool) {
    ignore_bool[i] <- TRUE
  }
}
cat('Done. Ignoring [', sum(ignore_bool),'] and keeping [',
    sum(!ignore_bool),'].\n', sep="")

# TOP-10
cat("And the top 100 NCBI contrast N nodes are....\n")
tmp_nodes <- node_obj[!ignore_bool]
tmp_nms <- nms_obj[!ignore_bool]
cntrst_ns <- unlist(lapply(tmp_nodes, function(x) min(x[['cntrst_n']])))
nms <- unlist(lapply(tmp_nms, function(x) x[['scientific name']]))
ordrd <- order(cntrst_ns)[1:100]
for(i in ordrd) {
  nm <- nms[i]
  cn <- cntrst_ns[i]
  spcr <- paste0(rep(' ', 33 - nchar(nm)), collapse="")
  cat(nm, spcr, signif(cn, 3), "\n")
}

# OUTPUT
save(node_obj, nms_obj, ignore_bool, file=output_file)

# END
cat(paste0('\nStage `NCBI taxonomy` finished at [', Sys.time(), ']\n'))