# GET NODE TIMINGS WITH TIMETREE

# START
cat(paste0('\nStage `timetree` started at [', Sys.time(), ']\n'))

# PARAMETERS
source('parameters.R')

# FUNCTIONS
source(file.path('tools', 'timetree_tools.R'))

# DIRS
if(!file.exists("timetree_cache")) {
  dir.create("timetree_cache")
}
if(!file.exists('4_timetree')) {
  dir.create('4_timetree')
}
input_file <- file.path("3_phylotime", "res.RData")
output_file <- file.path('4_timetree', 'res.RData')

# INPUT
load(input_file)

# LOOK UP TIMETREE DIVERGENCES
cat('Searching Time Tree ....\n')
cc <- 0
for(txid in cnddts) {
  if('tmsplt' %in% names(node_obj[[txid]])) {
    cc <- cc + 1
    next
  }
  # get time since split
  tmsplt <- getTmsplt(txid)
  if(is.na(tmsplt['mean_ttol'])) {
    next
  }
  # assign
  node_obj[[txid]][['tmsplt']] <- tmsplt[['mean_ttol']]
  cc <- cc + 1
}
cat("Done. Time data now available for [", cc, "] nodes.\n", sep="")

# OUTPUT
cat('Saving ....\n')
save(node_obj, cnddts, file=output_file)
cat('Done.\n')

# END
cat(paste0('\nStage `timetree` finished at [', Sys.time(), ']\n'))