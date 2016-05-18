
getTTOL <- function(sp1, sp2) {
  # search timetree via html
  # TODO: avoid line numbers, control for failure
  getVal <- function(strng) {
    stop <- regexpr(" Mya", strng) - 1
    as.numeric(substr(strng, start=26, stop=stop))
  }
  fl <- file.path('timetree_cache', paste0(sp1, "_", sp2, ".RData"))
  if(file.exists(fl)) {
    load(fl)
  } else {
    Sys.sleep(4)  # prevent overloading timetree
    url <- "http://www.timetree.org/search/pairwise/"
    sp1 <- gsub("\\s+", "%20", sp1)
    sp2 <- gsub("\\s+", "%20", sp2)
    qry <- paste0(url, sp1, '/', sp2)
    res <- readLines(qry)
    unlink(qry)
    mean_ttol <- getVal(res[[195]])
    median_ttol <- getVal(res[[199]])
    res <- data.frame(mean_ttol, median_ttol)
    save(res, file=fl)
  }
  res
}

findCommonName <- function(ids) {
  # Return the most likely name to be found in timetree
  # (Based on number of name entries in genbank)
  is <- which(txids %in% ids)
  nnms <- unlist(lapply(nms_obj[is], function(x) length(x)))
  nnms <- which.max(nnms)
  if(length(nnms) > 1) {
    nnms <- sample(nnms, 1)
  }
  i <- is[nnms]
  nms_obj[[i]][['scientific name']]
}

getDivergence <- function(id1, id2) {
  i <- which(txids == id1)
  j <- which(txids == id2)
  kids_1 <- node_obj[[i]][['kids']]
  kids_2 <- node_obj[[j]][['kids']]
  kid_1 <- findCommonName(kids_1)
  kid_2 <- findCommonName(kids_2)
  getTTOL(kid_1, kid_2)
}


if(!file.exists("timetree_cache")) {
  dir.create("timetree_cache")
}

nt_dir <- "x_ncbi_taxonomy"

load(file.path(nt_dir, "res.Rdata"))

txids <- unlist(lapply(node_obj, function(x) x[['txid']]))

for(i in which(!ignore_bool)) {
  print(i)
  # get time since split
  sids <- node_obj[[i]][['sstr']]
  if(length(sids) == 0) {
    next
  }
  txid <- node_obj[[i]][['txid']]
  tmsplt <- getDivergence(txid, sids[1])
  if(length(sids) > 1) {
    for(j in 2:length(sids)) {
      tmp <- getDivergence(txid, sids[j])
      bool <- tmsplt[['mean_ttol']] > tmp[['mean_ttol']]
      if(!is.na(bool) && bool) {
        tmsplt <- tmp
      }
    }
  }
  if(is.na(tmsplt['mean_ttol']) || tmsplt['mean_ttol'] < 50) {
    next
  }
  # get age
  ptids <- node_obj[[i]][['ptid']]
  if(length(ptids) > 1) {
    cmbs <- combn(ptids, 2)
    age <- getDivergence(cmbs[1,1], cmbs[2,1])
    if(ncol(cmbs) > 1) {
      for(j in 2:length(sids)) {
        tmp <- getDivergence(cmbs[1,j], cmbs[2,j])
        bool <- age[['mean_ttol']] > tmp[['mean_ttol']]
        if(!is.na(bool) && bool) {
          age <- tmp
        }
      }
    }
  } else {
    next
  }
  node_obj[[i]][['age']] <- age
  node_obj[[i]][['tmsplt']] <- tmsplt
  node_obj[[i]][['spn']] <- tmsplt - age
}
