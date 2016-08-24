# DOWNLOAD IUCN DATA ON LIVING FOSSILS AND NULL SPECIES

# START
cat(paste0('\nStage `iucn download` started at [', Sys.time(), ']\n'))

# FUNCTIONS
source(file.path('tools', 'hbbt_tools.R'))
source(file.path('tools', 'iucn_dwnld_tools.R'))
source(file.path('tools', 'node_obj_tools.R'))

# PARAMETERS
source('parameters.R')
token <- getToken()
min_nchar <- 3  # minimum number of characters in a meaningful word

# DIRS
output_dir <- '9_iucn_dwnld'
if (!file.exists(output_dir)) {
  dir.create(output_dir)
}
input_file <- file.path("8_epi", "res.RData")

# INPUT
load(input_file)
epi <- epi[!duplicated(epi$txid), ]

# LOOP
cds <- hbbts <- vector("list", length=nrow(epi))
cat('    Searching IUCN for habitats ....')
for(i in 1:nrow(epi)) {
  nms <- getKidNms(epi[i,'txid'])
  res_h <- res_c <- NULL
  for(nm in nms) {
    hs <- getIUCNHbbts(nm, token)
    if(class(hs) == "list" && length(hs[['result']]) > 0) {
      ss <- sapply(hs[['result']], function(x) x[['suitability']])
      cs <- sapply(hs[['result']], function(x) x[['code']])
      hs <- sapply(hs[['result']], function(x) x[['habitat']])
      res_h <- c(res_h, hs[ss == "Suitable"])
      res_c <- c(res_c, cs[ss == "Suitable"])
    }
  }
  if(length(res_h) > 0) {
    hbbts[[i]] <- res_h
    cds[[i]] <- res_c
  }
}
whbbts <- which(sapply(hbbts, length) > 0)
cat("Done, found data for [", length(whbbts), "/",
    length(hbbts), "].\n", sep="")

# IUCN HABITAT TYPES
htypes_of_interest <- c("forest", "subterranean",
                        "wetlands", "tundra",
                        "boreal", "temperate",
                        "rocky", "tropical")
ht_tests <- vector("list", length=length(htypes_of_interest))
names(ht_tests) <- htypes_of_interest
ht_pull <- scrs <- NULL
for(ht in htypes_of_interest) {
  vals <- as.numeric(sapply(hbbts[whbbts], function(x) any(grepl(ht, tolower(x)))))
  test <- wilcox.test(epi$pepi[whbbts][vals == 1],
                      epi$pepi[whbbts][vals == 0])
  means <- tapply(epi$pepi[whbbts], factor(vals), mean)
  ht_tests[[ht]] <- list(test, means)
  if(test$p.value < 0.05) {
    scrs <- c(scrs, vals)
    ht_pull <- c(ht_pull, TRUE)
  } else {
    ht_pull <- c(ht_pull, FALSE)
  }
}
hts <- factor(rep(htypes_of_interest[ht_pull], each=length(whbbts)))
hbbts_scrs <- data.frame(scr=scrs, type=hts, pepi=epi$pepi[whbbts])
p <- ggBinomial(hbbts_scrs)
p

# IUCN HABITAT TYPES BY CODE
cds_of_interest <- c("7.1", "7.2")  # subterranean
cds_of_interest <- as.character(seq(1.1, 1.9, .1))  # forests
cd_tests <- vector("list", length=length(cds_of_interest))
names(cd_tests) <- cds_of_interest
cd_pull <- scrs <- NULL
for(cd in cds_of_interest) {
  vals <- as.numeric(sapply(cds[whbbts], function(x) cd %in% x))
  test <- wilcox.test(epi$pepi[whbbts][vals == 1],
                      epi$pepi[whbbts][vals == 0])
  means <- tapply(epi$pepi[whbbts], factor(vals), mean)
  cd_tests[[cd]] <- list(test, means)
  if(test$p.value < 0.05) {
    scrs <- c(scrs, vals)
    cd_pull <- c(cd_pull, TRUE)
  } else {
    cd_pull <- c(cd_pull, FALSE)
  }
}
cds_scrs <- data.frame(scr=scrs, type=
                         factor(rep(cds_of_interest[cd_pull],
                                    each=length(whbbts))),
                       pepi=epi$pepi[whbbts])
p <- ggBinomial(cds_scrs)
p