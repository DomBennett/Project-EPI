# CALCULATE EPIS FROM NODE_OBJ

# START
cat(paste0('\nStage `epi` started at [', Sys.time(), ']\n'))

# PARAMETERS
cat("Loading parameters ...\n")
source('parameters.R')

# LIBS
cat("Loading libraries ...\n")
source(file.path("tools", "epi_tools.R"))
source(file.path('tools', 'node_obj_tools.R'))
source(file.path('tools', 'url_tools.R'))
source(file.path('tools', 'wiki_tools.R'))
source(file.path('tools', 'i_tools.R'))

# DIRS
cc_dir <- "8_cntrst_chng"
tt_dir <- "archive/7_timetree_by_name"
output_dir <- '9_epi'
if (!file.exists(output_dir)) {
  dir.create(output_dir)
}
output_file <- file.path(output_dir, "res.RData")
cc_file <- file.path(cc_dir, "res.RData")
tt_file <- file.path(tt_dir, "res.RData")

# CLADE DATA FRAME
cat("Generating clade dataframe.... ")
load(cc_file)
cnddts <- ls(node_obj)
cld_data <- genDataframe(cnddts, node_obj)
chngobj <- getCntrstChngObj(node_obj)
load(tt_file)
cld_data <- rbind(cld_data, genDataframe(cnddts, node_obj))
cld_data$cntrst_chng <- NA
indxs <- match(names(chngobj), cld_data$txid)
cld_data$cntrst_chng[indxs] <- sapply(chngobj, mean)
rm(node_obj)
cat("Done.\n")

# CLEAN-UP
# some data points don't make sense, reasons as yet unknown
cld_data <- cld_data[cld_data[['cntrst_n']] != Inf, ]  # two instances of this

# CALCULATE EPIS
cat("Calculating EPIs .... ")
cld_data$epi <- log((cld_data[['cntrst_chng']] + cld_data[['cntrst_n']]) / cld_data[['tmsplt']])
cld_data$pepi <- log(cld_data[['cntrst_n']] / cld_data[['tmsplt']])
#cld_data <- cld_data[!is.na(cld_data$epi), ]
cat("Done.\n")

# POST-FILTERING
# remove non-taxonomic entries
# mycorrhizal specimen is its own taxonomic clade, although practical this cannot be the true taxonomy
pull <- cld_data[['nm']] != "vouchered mycorrhizae (Pezizomycotina)"
# any names with numbers indicate species defined from genomic sequencing, true taxonomy cannot be known
pull <- pull & !grepl("[0-9]", cld_data[['nm']])
# environmental samples
pull <- pull & !grepl("environmental", cld_data[['nm']])
cld_data <- cld_data[pull, ]

# ADDING TAXONOMIC INFO
cat("Adding taxonomic information .... ")
txids <- cld_data[['txid']]
cld_data$txnmcgrp <- NA
grps <- c('metazoa', 'vertebrates', "mammals",
          "birds", "bony_fish", "lepidosaurs",
          "amphibia", 'plants', 'arthropods')
for(grp in grps) {
  grpids <- getGrpTxids(txids, grp)
  cld_data$txnmcgrp[cld_data$txid %in% grpids] <- grp
}
cat("Done.\n")

# WIKI SEARCH
cat("Searching wikipedia for living fossil mentions .... ")
cld_data$wiki <- NA
for(i in 1:nrow(cld_data)) {
  iPrnt(i, nrow(cld_data))
  txid <- cld_data[['txid']][i]
  # return TRUE if "living fossil" in wiki article
  # return FALSE if not
  # return NA if no wiki article found
  cld_data[['wiki']][i] <- getWikiLFMention(txid)
}
cat("Done.\n")

# OUTPUT
cat("Outputting ... ")
top100 <- cld_data[order(cld_data[['pepi']])[1:100], ]
write.csv(top100, file.path(output_dir, "top100.csv"), row.names=FALSE)
#save(node_obj, cld_data, file=output_file)
save(cld_data, chngobj, file=output_file)
cat('Done.\n')

# END
cat(paste0('\nStage `epi` finished at [', Sys.time(), ']\n'))