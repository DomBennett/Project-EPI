# Search for individual ids on timetree again

# LIBS
source(file.path('tools', 'timetree_tools.R'))
source(file.path('tools', 'url_tools.R'))

# PARAMS
wt <- c(1, 3, 10, 60, 120, 600, 3600, 7200)  # wait times
max_tt <- 10
try_again <- TRUE

# DIRS
output_file <- file.path('7_timetree', 'res.RData')
load(output_file)

# SEARCH
txid <- '7736'
txid <- '9255'
tmsplt <- getTmsplt(txid)

# CORRECT
node_obj[[txid]][['tmsplt']] <- tmsplt
save(node_obj, cnddts, file=output_file)
