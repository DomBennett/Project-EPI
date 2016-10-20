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
tt_dir <- "7_timetree"
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
# MANUAL EDITS
# birds
node_obj[['8782']][['tmsplt']] <- 236.5
node_obj[['8782']][['sstr']] <- '1294634'
node_obj[['8782']][['tmsplt']] <- 236.5
sstr_n <- 26  # crocodylia
n <- length(node_obj[['8782']][['kids']])
node_obj[['8782']][['cntrst_n']] <- n/sstr_n
# mammals
node_obj[['40674']][['tmsplt']] <- 311.9
# gen cld_data
cc_ids <- ls(node_obj)
cld_data <- genDataframe(cc_ids, node_obj)
cchngs <- getCntrstChng(cc_ids, node_obj)
nchars <- getNChars(cc_ids, node_obj)
load(tt_file)
# MANUAL EDITS (to get estimates for key taxa)
cnddts <- c(cnddts, '8292', '8504')
node_obj[['8292']][['tmsplt']] <- 351.8
node_obj[['8504']][['tmsplt']] <- 279.7
cld_data <- rbind(cld_data, genDataframe(cnddts, node_obj))
rm(node_obj)
cat("Done.\n")

# CONTRAST CHANGE
cld_data$cntrst_chng <- NA
indxs <- match(cc_ids, cld_data$txid)
pull <- nchars > 4
cld_data$cntrst_chng[indxs[pull]] <- cchngs[pull]
plot(y=log(cchngs[pull]), x=log(nchars[pull]))
abline(lm(log(cchngs[pull]+1)~log(nchars[pull])))
tests <- cor.test(cchngs[pull], nchars[pull], method="spearman")
if(tests$p.value < 0.05 & tests$estimate > 0.1) {
  warning('nchars and contrast changes correlate')
}

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
pull <- cld_data[['scinm']] != "vouchered mycorrhizae (Pezizomycotina)"
# any names with numbers indicate species defined from genomic sequencing, true taxonomy cannot be known
pull <- pull & !grepl("[0-9]", cld_data[['scinm']])
# environmental samples
pull <- pull & !grepl("environmental", cld_data[['scinm']])
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
  #n <- length(node_obj[[anlyss_grps[[grp]]]][['kids']])
  cld_data$txnmcgrp[cld_data$txid %in% grpids] <- grp
  #cat('.... [', n, '] species for [', grp, ']\n')
}
# remove anything that isn't in these groups, I don't trust the species counts
cld_data <- cld_data[!is.na(cld_data$txnmcgrp), ]
# lepidoptera has more species than coleoptera:
# nleps <- length(node_obj[['7088']][['kids']])
# ncols <- length(node_obj[['7041']][['kids']])
#cat('.... nleps: [', nleps, '], ncols [', ncols, ']\n')
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
cld_data <- cld_data[order(cld_data[['pepi']]), ]
write.csv(cld_data, file.path(output_dir, "living_fossils.csv"), row.names=FALSE)
#save(node_obj, cld_data, file=output_file)
save(cld_data, file=output_file)
cat('Done.\n')

# END
cat(paste0('\nStage `epi` finished at [', Sys.time(), ']\n'))