
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
  for(ln in lns) {
    if(grepl('class="time"', ln) &
      grepl('TTOL', ln)) {
      res <- gsub("<[^>]*>", "", ln)
      res <- as.numeric(gsub("[^\\.0-9]", "", res))
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
    url <- "http://www.timetree.org/search/pairwise/"
    qry <- paste0(url, id1, '/', id2)
    mean_ttol <- getMeanVal(searchURL(qry, site))
    res <- mean_ttol
    save(res, file=fl)
  }
  res
}

getTmsplt <- function(txid) {
  # Return time since split based on divergence from sister
  fl <- file.path(tmsplt_dir, paste0(txid, ".RData"))
  if(file.exists(fl) & !try_again) {
    load(fl)
    return(tmsplt)
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