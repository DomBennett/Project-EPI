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
cat('Searching Time Tree for pendant edges ....\n')
cc <- 0
for(txid in cnddts) {
  if('pe' %in% names(node_obj[[txid]])) {
    cc <- cc + 1
    next
  }
  print(node_obj[[txid]][['nm']][['scientific name']])
  # get time since split
  tmsplt <- getTmsplt(txid)
  if(is.na(tmsplt['mean_ttol'])) {
    next
  }
  # get ages
  if(node_obj[[txid]][['rank']] == 'species') {
    age <- 0
  } else {
    age <- getAge(txid)
  }
  if(is.na(age['mean_ttol']) || age['mean_ttol'] > tmsplt['mean_ttol']) {
    next
  }
  sstrs <- node_obj[[txid]][['sstr']]
  if(length(sstrs) > 1) {
    next
  }
  sstr_age <- NA
  for(sstr in sstrs) {
    tmp <- getAge(sstr)
    # always go for the greatest time
    bool <- is.na(sstr_age) || sstr_age[['mean_ttol']] > tmp[['mean_ttol']]
    if(bool) {
      sstr_age <- tmp
    }
  }
  if(is.na(sstr_age[['mean_ttol']])) {
    next
  }
  # get PEs
  pe <- tmsplt[['mean_ttol']] - age[['mean_ttol']]
  sstr_pe <- tmsplt[['mean_ttol']] - sstr_age[['mean_ttol']]
  cntrst_pe <- (pe + 1)/(sstr_pe + 1)
  # assign
  node_obj[[txid]][['age']] <- age
  node_obj[[txid]][['pe']] <- pe
  node_obj[[txid]][['cntrst_pe']] <- cntrst_pe
  cc <- cc + 1
}
cat("Done. Time data now available for [", cc, "] nodes.\n", sep="")

# OUTPUT
cat('Saving ....\n')
save(node_obj, cnddts, file=output_file)
cat('Done.\n')

# END
cat(paste0('\nStage `timetree` finished at [', Sys.time(), ']\n'))