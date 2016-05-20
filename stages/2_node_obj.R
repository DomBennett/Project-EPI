# PARSE NCBI TAXONOMY DMP AND GENERATE NODE_OBJ
# ftp://ftp.ncbi.nih.gov/pub/taxonomy/ (downloaded: 16/05/2016)

# START
cat(paste0('\nStage `node_obj` started at [', Sys.time(), ']\n'))

# PARAMETERS
source('parameters.R')

# DIRS
if(!file.exists("1_node_obj")) {
  dir.create("1_node_obj")
}
output_file <- file.path("1_node_obj", "res.RData")
input_dir <- file.path('0_data', 'ncbi_taxonomy')

# INPUT
cat("Reading in nodes and names dump ....\n")
nodes_lines <- readLines(con=file.path(input_dir, 'nodes.dmp'), n=-1L)
nodes_lines <- strsplit(nodes_lines, split="\\|")
names_lines <- readLines(con=file.path(input_dir, 'names.dmp'), n=-1L)
names_lines <- strsplit(names_lines, split="\\|")
cat("Done.\n")

# GENERATE NODE OBJECT
cat("Generating node object ....\n")
node_obj <- new.env()
for(i in 1:length(nodes_lines)) {
  division_code <- as.numeric(nodes_lines[[i]][5])
  if(division_code %in% division_codes) {
    txid <- gsub('\t', '', nodes_lines[[i]][1])
    prid <- gsub('\t', '', nodes_lines[[i]][2])
    rank <- gsub('\t', '', nodes_lines[[i]][3])
    node_obj[[txid]] <- list("prid"=prid, "rank"=rank,
                             'dcode'=division_code)
  }
}
cat("Done.\n")

# PARSE NAMES
cat("Parsing names ....\n")
txids <- ls(node_obj)
for(i in 1:length(names_lines)) {
  txid <- gsub("\t", "", names_lines[[i]][1])
  if(txid %in% txids) {
    nm <- gsub('\t', '', names_lines[[i]][2])
    nm_typ <- gsub('\t', '', names_lines[[i]][4])
    nms <- node_obj[[txid]][['nm']]
    node_obj[[txid]][['nm']] <- c(nms, nm)
    names(node_obj[[txid]][['nm']]) <- c(names(nms), nm_typ)
  }
}
cat("Done.\n")

# OUTPUT
cat('Saving ....\n')
save(node_obj, file=output_file)
cat('Done.\n')

# END
cat(paste0('\nStage `node_obj` finished at [', Sys.time(), ']\n'))