# SET PARAMETERS
parallel <- TRUE
if(parallel) {
  library(doMC)
  ncps <- detectCores()
  registerDoMC(cores=ncps)
}
stdy_grp <- 'mammal'
