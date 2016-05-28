# DOWNLOAD IUCN DATA ON LIVING FOSSILS AND NULL SPECIES

# START
cat(paste0('\nStage `iucn` started at [', Sys.time(), ']\n'))

# FUNCTIONS
source(file.path('tools', 'iucn_tools.R'))
source(file.path('tools', 'node_obj_tools.R'))

# PARAMETERS
source('parameters.R')
token <- getToken()

# DIRS
output_dir <- '7_iucn'
if (!file.exists(output_dir)) {
  dir.create(output_dir)
}
input_file <- file.path("6_epi", "res.RData")

# INPUT
load(input_file)

# LOOP THROUGH ANALYSIS GROUP AND DOWNLOAD
cat('Looping through each analysis group ....\n')
res <- list()
grps <- names(anlyss_grps)
for(grp in grps) {
  # GET GROUP IDS
  cat('    Working on [', grp, '] ....\n', sep="")
  txids <- ls(node_obj)
  txids <- getGrpTxids(txids, grp=grp)
  spp <- getSppTxids(txids)
  output_file <- file.path('7_iucn', paste0(grp, ".RData"))
  
  # GET LIVING FOSSILS
  cat("    Finding [", nlfs, "] top most living fossil clades ....\n",
      sep="")
  epi <- epi[order(epi[['nt_indx_lg']], decreasing=FALSE), ]
  lf_txids <- epi[['txid']][epi[['txid']] %in% txids][1:nlfs]
  lf_data <- vector("list", length=length(lf_txids))
  names(lf_data) <- lf_txids
  cat("    Done.\n")
  
  # SEARCH IUCN FOR LIVING FOSSILS
  cat('    Searching IUCN for living fossil data ....\n')
  for(i in 1:length(lf_data)) {
    txid <- names(lf_data)[i]
    nms <- getKidNms(txid)
    res <- list()
    for(j in 1:length(nms)) {
      #cat('Searching [', nms[j], '] ....\n', sep="")
      cate <- getIUCNCat(nms[j], token)
      res_nm <- list()
      if(class(cate) == "list" && length(cate[['result']]) > 0) {
        res_nm[["cate"]] <- cate[['result']][[1]][['category']]
      }
      nrrtv <- getIUCNNrrtv(nms[j], token)
      if(class(nrrtv) == "list" && length(nrrtv[['result']]) > 0) {
        res_nm[["dscrptn"]] <- nrrtv[['result']][[1]][['habitat']]
      }
      hbbts <- getIUCNHbbts(nms[j], token)
      if(class(hbbts) == "list" && length(hbbts[['result']]) > 0) {
        nhbbts <- sum(unlist(lapply(hbbts[['result']],
                                    function(x) x[['suitability']] == "Suitable")))
        res_nm[['nhbbts']] <- nhbbts
      }
      cntrs <- getIUCNCntrs(nms[j], token)
      if(class(cntrs) == "list" && length(cntrs[['result']]) > 0) {
        ncntrs <- sum(unlist(lapply(cntrs[['result']],
                                    function(x) x[['origin']] == "Native")))
        res_nm[['ncntrs']] <- ncntrs
      }
      if(length(res_nm) > 0) {
        res[[nms[j]]] <- res_nm
      }
    }
    if(length(res) > 0) {
      lf_data[[txid]] <- res
    }
  }
  lf_data <- lf_data[unlist(lapply(lf_data, function(x) length(x[[1]]) > 0))]
  cat("    Done. Found data for [", length(lf_data), "] ....\n", sep="")
  
  # SEARCH IUCN FOR NULL CATES
  cat('    Searching IUCN for null category data ....\n')
  null_cate <- list()
  ntrys <- 0
  for(itrtn in 1:nitrtns) {
    cat('        Iteration [', itrtn, '] ....\n', sep="")
    null_cate[[itrtn]] <- list()
    cc <- 0
    while(cc < length(lf_data)) {
      sp <- sample(spp, size=1)
      nm <- node_obj[[sp]][['nm']][['scientific name']]
      cate <- getIUCNCat(nm, token)
      if(class(cate) == "list" && length(cate[['result']]) > 0) {
        null_cate[[itrtn]][[nm]] <- cate[['result']][[1]][['category']]
        cc <- cc + 1
        ntrys <- 0
      } else {
        ntrys <- ntrys + 1
      }
      if(ntrys > 100) {
        cat("    ----- can't find category data ----")
        break
      }
    }
  }
  cat("    Done.\n")
  
  # SEARCH IUCN FOR NULL HABITAT AND ECOLOGY DESCRIPTIONS
  cat('    Searching IUCN for null description data ....\n')
  null_dscrptn <- list()
  ntrys <- 0
  for(itrtn in 1:nitrtns) {
    cat('        Iteration [', itrtn, '] ....\n', sep="")
    null_dscrptn[[itrtn]] <- list()
    cc <- 0
    while(cc < length(lf_data)) {
      sp <- sample(spp, size=1)
      nm <- node_obj[[sp]][['nm']][['scientific name']]
      nrrtv <- getIUCNNrrtv(nm, token)
      if(class(nrrtv) == "list" && length(nrrtv[['result']]) > 0 &&
         !is.null(nrrtv[['result']][[1]][['habitat']])) {
        null_dscrptn[[itrtn]][[nm]] <- nrrtv[['result']][[1]][['habitat']]
        cc <- cc + 1
        ntrys <- 0
      } else {
        ntrys <- ntrys + 1
      }
      if(ntrys > 100) {
        cat("    ----- can't find description data ----")
        break
      }
    }
  }
  cat("    Done.\n")
  
  # SEARCH IUCN FOR NULL HABITAT
  cat('    Searching IUCN for null habitat data ....\n')
  null_nhbbts <- list()
  ntrys <- 0
  for(itrtn in 1:nitrtns) {
    cat('        Iteration [', itrtn, '] ....\n', sep="")
    null_nhbbts[[itrtn]] <- list()
    cc <- 0
    while(cc < length(lf_data)) {
      sp <- sample(spp, size=1)
      nm <- node_obj[[sp]][['nm']][['scientific name']]
      hbbts <- getIUCNHbbts(nm, token)
      if(class(hbbts) == "list" && length(hbbts[['result']]) > 0) {
        nhbbts <- sum(unlist(lapply(hbbts[['result']],
                                    function(x) x[['suitability']] == "Suitable")))
        null_nhbbts[[itrtn]][[nm]] <- nhbbts
        cc <- cc + 1
        ntrys <- 0
      } else {
        ntrys <- ntrys + 1
      }
      if(ntrys > 100) {
        cat("    ----- can't find habitats data ----")
        break
      }
    }
  }
  cat("    Done.\n")
  
  # SEARCH IUCN FOR NULL COUNTRY DATA
  cat('    Searching IUCN for null country data ....\n')
  null_ncntrs <- list()
  ntrys <- 0
  for(itrtn in 1:nitrtns) {
    cat('        Iteration [', itrtn, '] ....\n', sep="")
    null_ncntrs[[itrtn]] <- list()
    cc <- 0
    while(cc < length(lf_data)) {
      sp <- sample(spp, size=1)
      nm <- node_obj[[sp]][['nm']][['scientific name']]
      cntrs <- getIUCNCntrs(nm, token)
      if(class(cntrs) == "list" && length(cntrs[['result']]) > 0) {
        ncntrs <- sum(unlist(lapply(cntrs[['result']],
                                    function(x) x[['origin']] == "Native")))
        null_ncntrs[[itrtn]][[nm]] <- ncntrs
        cc <- cc + 1
        ntrys <- 0
      } else {
        ntrys <- ntrys + 1
      }
      if(ntrys > 100) {
        cat("    ----- can't find country data ----")
        break
      }
    }
  }
  cat("    Done.\n")
  
  # OUTPUT
  cat('    Saving ....\n')
  save(lf_data, null_cate, null_dscrptn, null_nhbbts,
       null_ncntrs, file=output_file)
  cat('    Done.\n')
}

# END
cat(paste0('\nStage `iucn` finished at [', Sys.time(), ']\n'))