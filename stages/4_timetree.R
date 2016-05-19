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

# LOOK UP SPN
cat('Searching Time Tree for pendant edges ....\n')
to_drp <- vector(length=length(cnddts))
for(i in 1:length(cnddts)) {
  txid <- cnddts[i]
  if('pe' %in% names(node_obj[[txid]])) {
    next
  }
  # get time since split
  sstrs <- node_obj[[txid]][['sstr']]
  tmsplt <- getDivergence(txid, sstrs[1])
  if(length(sids) > 1) {
    for(j in 2:length(sids)) {
      tmp <- getDivergence(txid, sids[j])
      bool <- tmsplt[['mean_ttol']] > tmp[['mean_ttol']]
      if(!is.na(bool) && bool) {
        tmsplt <- tmp
      }
    }
  }
  if(is.na(tmsplt['mean_ttol']) || tmsplt['mean_ttol'] < 50) {
    to_drp[i] <- TRUE
    next
  }
  # get age
  ptids <- node_obj[[i]][['ptid']]
  cmbs <- combn(ptids, 2)
  age <- getDivergence(cmbs[1,1], cmbs[2,1])
  if(ncol(cmbs) > 1) {
    for(j in 2:length(sids)) {
      tmp <- getDivergence(cmbs[1,j], cmbs[2,j])
      bool <- age[['mean_ttol']] > tmp[['mean_ttol']]
      if(!is.na(bool) && bool) {
        age <- tmp
      }
    }
  }
  node_obj[[i]][['age']] <- age
  node_obj[[i]][['tmsplt']] <- tmsplt
  node_obj[[i]][['pe']] <- tmsplt - age
}
cnddts <- cnddts[!to_drp]
cat("Done. Time data now available for [", length(cnddts), "] nodes.\n", sep="")

# OUTPUT
cat('Saving ....\n')
save(node_obj, cnddts, file=output_file)
cat('Done.\n')

# END
cat(paste0('\nStage `timetree` finished at [', Sys.time(), ']\n'))