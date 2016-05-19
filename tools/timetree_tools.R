
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
  nnms <- unlist(lapply(ids, function(x) length(node_obj[[x]][['nm']])))
  i <- which.max(nnms)
  if(length(i) > 1) {
    i <- sample(i, 1)
  }
  node_obj[[ids[i]]][['nm']][['scientific name']]
}

getDivergence <- function(id1, id2) {
  kids_1 <- node_obj[[id1]][['kids']]
  kids_2 <- node_obj[[id2]][['kids']]
  if(length(kids_1) > 1) {
    nm1 <- findCommonName(kids_1)
  } else {
    nm1 <- node_obj[[id1]][['nm']][['scientific name']]
  }
  if(length(kids_2) > 1) {
    nm2 <- findCommonName(kids_2)
  } else {
    nm2 <- node_obj[[id2]][['nm']][['scientific name']]
  }
  getTTOL(nm1, nm2)
}

getTmsplt <- function(txid) {
  # Return time since split based on divergence from sister
  sstrs <- node_obj[[txid]][['sstr']]
  tmsplt <- getDivergence(txid, sstrs[1])
  if(length(sstrs) > 1) {
    if(length(sstrs) > 10) {
      sstrs <- sample(sstrs, 10)
    }
    for(j in 2:length(sstrs)) {
      tmp <- getDivergence(txid, sstrs[j])
      # always seek greatest time
      bool <- tmsplt[['mean_ttol']] > tmp[['mean_ttol']]
      if(!is.na(bool) && bool) {
        tmsplt <- tmp
      }
    }
  }
  tmsplt
}

getAge <- function(txid) {
  # Return age based on divergence of ptids
  ptids <- node_obj[[txid]][['ptid']]
  cmbs <- combn(ptids, 2)
  age <- getDivergence(cmbs[1,1], cmbs[2,1])
  if(ncol(cmbs) > 1) {
    if(ncol(cmbs) > 10) {
      cmbs <- cmbs[ ,sample(2:ncol(cmbs), 10)]
    }
    for(j in 2:ncol(cmbs)) {
      tmp <- getDivergence(cmbs[1,j], cmbs[2,j])
      # always seek greatest time
      bool <- age[['mean_ttol']] > tmp[['mean_ttol']]
      if(!is.na(bool) && bool) {
        age <- tmp
      }
    }
  }
  age
}