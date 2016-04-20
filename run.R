source('parameters.R')
cat('Running all stages for [', stdy_grp, '] ....\n', sep='')
source(file.path('stages', '1_measurables.R'),local=TRUE)
source(file.path('stages', '2_epi.R'), local=TRUE)
