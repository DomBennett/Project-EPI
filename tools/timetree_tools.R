
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
    id1 <- findCommonName(kids_1)
  }
  if(length(kids_2) > 1) {
    id2 <- findCommonName(kids_2)
  }
  getTTOL(id1, id2)
}