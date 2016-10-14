
# DIRS
cache_dir <- "caches"
if(!file.exists(cache_dir)) {
  dir.create(cache_dir)
}
cache_dir <- file.path("caches", "timetree")
if(!file.exists(cache_dir)) {
  dir.create(cache_dir)
}
tmsplt_dir <- file.path("caches", "tmsplt")
if(!file.exists(tmsplt_dir)) {
  dir.create(tmsplt_dir)
}

# FUNCTIONS
getMeanVal <- function(lns) {
  res <- NA
  for(i in 2:length(lns)) {
    if(grepl('pairwise-results', lns[i])) {
      for(j in i:(i+20)) {
        if(grepl('TTOL', lns[j])) {
          break
        }
      }
      for(k in j:i) {
        if(grepl('MYA', lns[k], ignore.case = FALSE)) {
          res <- gsub("<[^>]*>", "", lns[k])
          res <- as.numeric(gsub("[^\\.0-9]", "", res))
          break
        }
      }
    }
  }
  res
}

getTTOL <- function(id1, id2) {
  # search timetree via html
  fl <- file.path(cache_dir, paste0(id1, "_", id2, ".RData"))
  if(file.exists(fl)) {
    load(fl)
  } else {
    res <- searchByID(id1, id2)
    save(res, file=fl)
  }
  res
}

searchByID <- function(id1, id2) {
  url <- "http://www.timetree.org/ajax/pairwise/"
  qry <- paste0(url, id1, '/', id2)
  lns <- searchURL(qry, site="timetree")
  getMeanVal(lns)
}

getTmsplt <- function(txid) {
  # Return time since split based on divergence from sister
  fl <- file.path(tmsplt_dir, paste0(txid, ".RData"))
  if(file.exists(fl)) {
    load(fl)
    if(!try_again | !is.na(tmsplt)) {
      return(tmsplt)
    }
  }
  sstrs <- node_obj[[txid]][['sstr']]
  tmsplt <- getTTOL(txid, sstrs[1])
  if(length(sstrs) > 1) {
    if(length(sstrs) > (max_tt + 1)) {
      sstrs <- sample(sstrs[-1], max_tt)
    }
    for(j in 2:length(sstrs)) {
      tmp <- getTTOL(txid, sstrs[j])
      # always seek greatest time
      bool <- tmsplt > tmp
      if(is.na(tmsplt) & is.numeric(tmp) |
         !is.na(bool) && bool) {
        tmsplt <- tmp
      }
    }
  }
  save(tmsplt, file=fl)
  tmsplt
}