# GET NODE TIMINGS WITH TIMETREE

# START
cat(paste0('\nStage `timetree` started at [', Sys.time(), ']\n'))

# PARAMETERS
source('parameters.R')

# FUNCTIONS
source(file.path('tools', 'timetree_tools.R'))
source(file.path('tools', 'i_tools.R'))

# DIRS
if(!file.exists('7_timetree')) {
  dir.create('7_timetree')
}
input_file <- file.path("5_split", "tt_obj.RData")
output_file <- file.path('7_timetree', 'res.RData')

# INPUT
cat('Loading input.... ')
load(input_file)
node_obj <- tt_obj
rm(tt_obj)
cat('Done.\n')

# LOOK UP TIMETREE DIVERGENCES
cat('Searching Time Tree ....\n')
txids <- ls(node_obj)
cnddts <- cnddts[cnddts %in% txids]
cc <- 0
cat('    [', length(cnddts), '] txids to search for ....\n')
for(i in 1:length(cnddts)) {
  txid <- cnddts[i]
  iPrnt(i, length(cnddts))
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
  # count
  cc <- cc + 1
}
cat("Done. Time data now available for [", cc, "] nodes.\n", sep="")

# OUTPUT
cat('Saving ....\n')
save(node_obj, cnddts, file=output_file)
cat('Done.\n')

# END
cat(paste0('\nStage `timetree` finished at [', Sys.time(), ']\n'))