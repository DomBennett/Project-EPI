# DOWNLOAD IUCN DATA ON LIVING FOSSILS AND NULL SPECIES

# START
cat(paste0('\nStage `iucn download` started at [', Sys.time(), ']\n'))

# FUNCTIONS
source(file.path('tools', 'i_tools.R'))
source(file.path('tools', 'iucn_dwnld_tools.R'))
source(file.path('tools', 'node_obj_tools.R'))

# PARAMETERS
source('parameters.R')
token <- getToken()

# DIRS
output_dir <- '9_iucn_dwnld'
if (!file.exists(output_dir)) {
  dir.create(output_dir)
}
input_file <- file.path("8_epi", "res.RData")

# INPUT
load(input_file)
epi <- epi[!duplicated(epi$txid), ]

grp <- "mammals"
cat('    Working on [', grp, '] ....\n', sep="")
txids <- ls(node_obj)
txids <- getGrpTxids(txids, grp=grp)
epi <- epi[epi[['txid']] %in% txids, ]

epi$nhbbts <- NA
for(i in 1:nrow(epi)) {
  nms <- getKidNms(epi[i,'txid'])
  all_n <- 0
  cc <- 0
  for(nm in nms) {
    hbbts <- getIUCNHbbts(nm, token)
    if(class(hbbts) == "list" && length(hbbts[['result']]) > 0) {
      cc <- cc + 1
      nhbbts <- sum(unlist(lapply(hbbts[['result']],
                                  function(x) x[['suitability']] == "Suitable")))
      all_n <- (all_n + nhbbts)/cc
    }
  }
  epi[['nhbbts']][i] <- all_n
}

# 

epi$rnkng <- order(epi$pepi)
plot(epi$nhbbts~epi$pepi)
hist(epi$nhbbts)
m1 <- lm(log(epi$nhbbts + 1)~epi$pepi)
summary(m1)
plot(m1)

plot(log(epi$nhbbts) ~ epi$pepi)
cor.test(log(epi$nhbbts + 1), epi$pepi)
cor.test(epi$nhbbts, epi$epi)

bool <- !is.na(epi$nhbbts) & epi$nhbbts > 20
epi[bool, ]

epi[epi$rnkng[1:10],]


