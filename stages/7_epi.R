# CALCULATE EPIS FROM NODE_OBJ

# START
cat(paste0('\nStage `epi` started at [', Sys.time(), ']\n'))

# PARAMETERS
cat("Loading parameters ...\n")
source('parameters.R')

# LIBS
cat("Loading libraries ...\n")
source(file.path("tools", "epi_tools.R"))

# DIRS
output_dir <- '7_epi'
if (!file.exists(output_dir)) {
  dir.create(output_dir)
}
input_file <- file.path("6_timetree", "res.RData")
output_file <- file.path(output_dir, "res.RData")

# INPUT
load(input_file)

# GENERATE EPI DATAFRAME
nms <- tmsplts <- cntrst_ns <- rep(NA, length(txids))
for(i in 1:length(txids)) {
  tmsplt <- node_obj[[txids[i]]][["tmsplt"]]
  cntrst_n <- node_obj[[txids[i]]][["cntrst_n"]]
  if(!is.null(tmsplt) & !is.null(cntrst_n)) {
    tmsplts[i] <- tmsplt
    cntrst_ns[i] <- cntrst_n
    nms[i] <- node_obj[[txids[i]]][["nm"]][["scientific name"]]
  }
}
bool <- !is.na(cntrst_ns)
cntrst_ns <- cntrst_ns[bool]
tmsplts <- tmsplts[bool]
txids <- txids[bool]
nms <- nms[bool]
epi <- data.frame(nm=nms, txid=txids, tmsplt=tmsplts, cntrst_n=cntrst_ns,
                  stringsAsFactors=FALSE)

# CALCULATE EPIS
epi$nt_indx <- epi[['cntrst_n']]/epi[['tmsplt']]
epi$nt_indx_lg <- log(epi[["nt_indx"]])

# VISUALISE

# OUTPUT
cat("Outputting ...\n")
save(node_obj, epi, file=output_file)

# END
cat(paste0('\nStage `epi` finished at [', Sys.time(), ']\n'))