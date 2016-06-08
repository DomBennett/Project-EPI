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
output_dir <- '8_iucn_dwnld'
if (!file.exists(output_dir)) {
  dir.create(output_dir)
}
input_file <- file.path("7_epi", "res.RData")

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
  output_file <- file.path(output_dir, paste0(grp, ".RData"))
  
  # GET LIVING FOSSILS
  cat("    Finding [", nlfs, "] top most living fossil clades ....",
      sep="")
  epi <- epi[order(epi[['nt_indx_lg']], decreasing=FALSE), ]
  lf_txids <- epi[['txid']][epi[['txid']] %in% txids][1:nlfs]
  lf_txids <- lf_txids[!is.na(lf_txids)]
  if(length(lf_txids) < 10) {
    cat("    \nToo few living fossils for this group!\n")
    next
  }
  cat("Done.\n")
  
  # SEARCH IUCN FOR LIVING FOSSIL CATES
  cat('    Searching IUCN for living fossil categories ....')
  ignr <- NULL  # ignore all DD species
  lf_cate <- list()
  for(txid in lf_txids) {
    nms <- getKidNms(txid)
    res <- list()
    for(nm in nms) {
      #cat('Searching [', nms[j], '] ....\n', sep="")
      cate <- getIUCNCat(nm, token)
      if(class(cate) == "list" && length(cate[['result']]) > 0) {
        cate <- cate[['result']][[1]][['category']]
        if(cate == "DD") {
          ignr <- c(ignr, nm)
        } else {
          res[[nm]] <- cate
        }
      }
    }
    if(length(res) > 0) {
      lf_cate[[txid]] <- res
    }
  }
  if(length(lf_cate) < 10) {
    cat("    \nToo little living fossil IUCN data for this group!\n")
    next
  }
  cat("Done, found data for [", length(lf_cate), "].\n", sep="")
  
  # SEARCH IUCN FOR LIVING FOSSIL NARRATIVES
  lf_dscrptn <- list()
  cat('    Searching IUCN for living fossil descriptions ....')
  for(txid in lf_txids) {
    nms <- getKidNms(txid)
    nms <- nms[!nms %in% ignr]
    res <- list()
    for(nm in nms) {
      nrrtv <- getIUCNNrrtv(nm, token)
      if(class(nrrtv) == "list" && length(nrrtv[['result']]) > 0 &&
         !is.null(nrrtv[['result']][[1]][['habitat']])) {
        nrrtv <- nrrtv[['result']][[1]][['habitat']]
        res[[nm]] <- gsub("<.*?>", "", nrrtv)  # remove html tags
      }
    }
    if(length(res) > 0) {
      lf_dscrptn[[txid]] <- res
    }
  }
  cat("Done, found data for [", length(lf_dscrptn), "].\n", sep="")
  
  # SEARCH IUCN FOR LIVING FOSSIL HABITATS
  lf_nhbbt <- list()
  cat('    Searching IUCN for living fossil habitats ....')
  for(txid in lf_txids) {
    nms <- getKidNms(txid)
    nms <- nms[!nms %in% ignr]
    res <- list()
    for(nm in nms) {
      hbbts <- getIUCNHbbts(nm, token)
      if(class(hbbts) == "list" && length(hbbts[['result']]) > 0) {
        nhbbts <- sum(unlist(lapply(hbbts[['result']],
                                    function(x) x[['suitability']] == "Suitable")))
        res[[nm]] <- nhbbts
      }
    }
    if(length(res) > 0) {
      lf_nhbbt[[txid]] <- res
    }
  }
  cat("Done, found data for [", length(lf_nhbbt), "].\n", sep="")
  
  # SEARCH IUCN FOR LIVING FOSSIL COUNTRIES
  lf_ncntr <- list()
  cat('    Searching IUCN for living fossil countries ....')
  for(txid in lf_txids) {
    nms <- getKidNms(txid)
    nms <- nms[!nms %in% ignr]
    res <- list()
    for(nm in nms) {
      cntrs <- getIUCNCntrs(nm, token)
      if(class(cntrs) == "list" && length(cntrs[['result']]) > 0) {
        ncntrs <- sum(unlist(lapply(cntrs[['result']],
                                    function(x) x[['origin']] == "Native")))
        res[[nm]] <- ncntrs
      }
    }
    if(length(res) > 0) {
      lf_ncntr[[txid]] <- res
    }
  }
  cat("Done, found data for [", length(lf_ncntr), "].\n", sep="")
  
  # SEARCH IUCN FOR NULL CATES
  cat('    Searching IUCN for null categories (iterating) ....\n')
  null_cate <- list()
  ntrys <- 0
  cat('    ')
  for(itrtn in 1:nitrtns) {
    iPrnt(itrtn, nitrtns)
    null_cate[[itrtn]] <- list()
    tmp_spp <- spp
    cc <- 0
    while(cc < length(lf_cate)) {
      sp <- sample(tmp_spp, size=1)
      nm <- node_obj[[sp]][['nm']][['scientific name']]
      cate <- getIUCNCat(nm, token)
      sccss <- FALSE
      if(class(cate) == "list" && length(cate[['result']]) > 0) {
        cate <- cate[['result']][[1]][['category']]
        sccss <- TRUE
      }
      if(sccss && cate != "DD") {
        null_cate[[itrtn]][[nm]] <- cate
        tmp_spp <- tmp_spp[tmp_spp != sp]
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
  cat("Done.\n")
  
  # SEARCH IUCN FOR NULL HABITAT AND ECOLOGY DESCRIPTIONS
  cat('    Searching IUCN for null descriptions (iterating) ....\n')
  null_dscrptn <- list()
  ntrys <- 0
  cat('    ')
  for(itrtn in 1:nitrtns) {
    iPrnt(itrtn, nitrtns)
    null_dscrptn[[itrtn]] <- list()
    tmp_spp <- spp
    cc <- 0
    while(cc < length(lf_dscrptn)) {
      sp <- sample(tmp_spp, size=1)
      nm <- node_obj[[sp]][['nm']][['scientific name']]
      nrrtv <- getIUCNNrrtv(nm, token)
      cate <- getIUCNCat(nm, token)
      sccss <- FALSE
      if(class(cate) == "list" && length(cate[['result']]) > 0) {
        cate <- cate[['result']][[1]][['category']]
        sccss <- cate != "DD"
      }
      if(sccss && class(nrrtv) == "list" && length(nrrtv[['result']]) > 0 &&
         !is.null(nrrtv[['result']][[1]][['habitat']])) {
        null_dscrptn[[itrtn]][[nm]] <- nrrtv[['result']][[1]][['habitat']]
        cc <- cc + 1
        tmp_spp <- tmp_spp[tmp_spp != sp]
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
  cat("       Done.\n")
  
  # SEARCH IUCN FOR NULL HABITAT
  cat('    Searching IUCN for null habitats (iterating) ....\n')
  null_nhbbt <- list()
  ntrys <- 0
  cat('    ')
  for(itrtn in 1:nitrtns) {
    iPrnt(itrtn, nitrtns)
    null_nhbbt[[itrtn]] <- list()
    cc <- 0
    tmp_spp <- spp
    while(cc < length(lf_nhbbt)) {
      sp <- sample(tmp_spp, size=1)
      nm <- node_obj[[sp]][['nm']][['scientific name']]
      hbbts <- getIUCNHbbts(nm, token)
      cate <- getIUCNCat(nm, token)
      sccss <- FALSE
      if(class(cate) == "list" && length(cate[['result']]) > 0) {
        cate <- cate[['result']][[1]][['category']]
        sccss <- cate != "DD"
      }
      if(sccss && class(hbbts) == "list" && length(hbbts[['result']]) > 0) {
        nhbbts <- sum(unlist(lapply(hbbts[['result']],
                                    function(x) x[['suitability']] == "Suitable")))
        null_nhbbt[[itrtn]][[nm]] <- nhbbts
        cc <- cc + 1
        tmp_spp <- tmp_spp[tmp_spp != sp]
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
  cat("Done.\n")
  
  # SEARCH IUCN FOR NULL COUNTRY DATA
  cat('    Searching IUCN for null countries (iterating) ....\n')
  null_ncntr <- list()
  ntrys <- 0
  cat('    ')
  for(itrtn in 1:nitrtns) {
    iPrnt(itrtn, nitrtns)
    null_ncntr[[itrtn]] <- list()
    cc <- 0
    tmp_spp <- spp
    while(cc < length(lf_ncntr)) {
      sp <- sample(tmp_spp, size=1)
      nm <- node_obj[[sp]][['nm']][['scientific name']]
      cntrs <- getIUCNCntrs(nm, token)
      cate <- getIUCNCat(nm, token)
      sccss <- FALSE
      if(class(cate) == "list" && length(cate[['result']]) > 0) {
        cate <- cate[['result']][[1]][['category']]
        sccss <- cate != "DD"
      }
      if(sccss && class(cntrs) == "list" && length(cntrs[['result']]) > 0) {
        ncntrs <- sum(unlist(lapply(cntrs[['result']],
                                    function(x) x[['origin']] == "Native")))
        null_ncntr[[itrtn]][[nm]] <- ncntrs
        cc <- cc + 1
        tmp_spp <- tmp_spp[tmp_spp != sp]
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
  cat('    Saving ....')
  save(lf_cate, lf_dscrptn, lf_nhbbt, lf_ncntr,
       null_cate, null_dscrptn, null_nhbbt, null_ncntr,
       file=output_file)
  cat('Done.\n')
}

# END
cat(paste0('\nStage `iucn download` finished at [', Sys.time(), ']\n'))