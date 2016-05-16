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
cat("Reading in nodes dump ....\n")
lines <- readLines(con=file.path(input_dir, 'nodes.dmp'), n=-1L)
lines <- strsplit(lines, split="\\|")
cat("Done.\n")

# CONVERT TO NODE-OBJ
cat("Parsing and generating node object ....\n")
node_obj <- vector("list", length=length(lines))
for(i in 1:length(node_obj)) {
  division_code <- as.numeric(lines[[i]][5])
  if(division_code %in% division_codes) {
    txid <- as.numeric(lines[[i]][1])
    prid <- as.numeric(lines[[i]][2])
    rank <- gsub('\t', '', lines[[i]][3])
    node_obj[[i]] <- list(txid=txid, prid=prid, ptid=NULL,
                          rank=rank, kids=NULL)
  }
}
fnds <- unlist(lapply(node_obj, function(x) !is.null(x)))
node_obj <- node_obj[fnds]
sbspp <- unlist(lapply(node_obj, function(x) x[['rank']] == "subspecies"))
node_obj <- node_obj[!sbspp]
cat("Done. Identified [", length(node_obj), "] valid taxonomic nodes.\n", sep="")

# RETRIEVE NAMES
cat("Searching names via NCBI ....\n")
for(i in 1:length(node_obj)) {
  txid <- node_obj[[i]][['txid']]
  nm <- fetchName(txid)
  node_obj[[i]][['name']] <- nm
}
cat('Done.\n')

# REMOVE IGNORE NAMES
cat("Dropping nodes with invalid names ....\n")
nms <- unlist(lapply(node_obj, function(x) x[['name']]))
to_drp <- grepl("unclassified", nms)
node_obj <- node_obj[!to_drp]
txids <- unlist(lapply(node_obj, function(x) x[['txid']]))
cat('Done. Dropped [', sum(to_drp),']\n', sep='')

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
  if(!is.null(sstrs) && length(sstrs) > 1) {
    cns <- vector(length=length(sstrs))
    names(cns) <- sstrs
    n <- length(length(node_obj[[i]][['kids']]))
    for(j in 1:length(sstrs)) {
      si <- which(txids == sstrs[j])
      sstr_n <- length(node_obj[[si]][['kids']])
      cns[j] <- n/sstr_n
    }
    node_obj[[i]][['cntrst_n']] <- cns
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
cat("And the top 10 NCBI contrast N nodes are....")
tmp <- node_obj[!ignore_bool]
cntrst_ns <- unlist(lapply(tmp, function(x) min(x[['cntrst_n']])))
nms <- unlist(lapply(tmp, function(x) x[['name']]))
print(nms[order(cntrst_ns)][1:10])

# OUTPUT
save(node_obj, ignore_bool, file=output_file)

# END
cat(paste0('\nStage `NCBI taxonomy` finished at [', Sys.time(), ']\n'))