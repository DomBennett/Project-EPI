
# DIRS
cache_dir <- "caches"
if(!file.exists(cache_dir)) {
  dir.create(cache_dir)
}
cache_dir <- file.path("caches", "timetree")
if(!file.exists(cache_dir)) {
  dir.create(cache_dir)
}

# IGNR FUNCTIONS
# functions that will prevent searching for names
# that have failed to be found
.ignr_scr_max <- 3 # maximum number of times a nm fails to be found before it is ignored
.ignr_file <- file.path(cache_dir, "ignr.RData")
if(!file.exists(.ignr_file)) {
  .ignr <- list("scrs"=list(), 'nms'=NULL)
} else {
  load(.ignr_file)
}
.updateIgnr <- function() {
  bool <- unlist(lapply(.ignr[['scrs']],
                        function(x) x > .ignr_scr_max))
  if(!is.null(bool)) {
    nms <- names(bool)[bool]
    .ignr[['nms']] <- c(.ignr[['nms']], nms)
    .ignr[['scrs']] <- .ignr[['scrs']][!bool]
    save(.ignr, file=.ignr_file)
    .ignr <<- .ignr
  }
  NULL
}
.ignrChck <- function(res, nm1, nm2) {
  if(all(is.na(res))) {
    if(!is.null(.ignr[['scrs']][[nm1]])) {
      .ignr[['scrs']][[nm1]] <-
        .ignr[['scrs']][[nm1]] + 1
    } else {
      .ignr[['scrs']][[nm1]] <- 1
    }
    if(!is.null(.ignr[['scrs']][[nm2]])) {
      .ignr[['scrs']][[nm2]] <-
        .ignr[['scrs']][[nm2]] + 1
    } else {
      .ignr[['scrs']][[nm2]] <- 1
    }
  }
}

# FUNCTIONS
getVal <- function(strng) {
  stop <- regexpr(" Mya", strng) - 1
  as.numeric(substr(strng, start=26, stop=stop))
}

getTTOL <- function(sp1, sp2) {
  # search timetree via html
  # TODO: avoid using line numbers
  fl <- file.path(cache_dir, paste0(sp1, "_", sp2, ".RData"))
  if(file.exists(fl)) {
    load(fl)
  } else {
    Sys.sleep(4)  # prevent overloading timetree
    ping <- suppressWarnings(try(readLines("http://www.timetree.org/", n=1),
                                 silent=TRUE))
    if(class(ping) == "try-error") {
      stop("---- Server is down ----")
    }
    url <- "http://www.timetree.org/search/pairwise/"
    sp1 <- gsub("\\s+", "%20", sp1)
    sp2 <- gsub("\\s+", "%20", sp2)
    qry <- paste0(url, sp1, '/', sp2)
    res <- suppressWarnings(try(expr=readLines(qry),
                                silent=TRUE))
    if(class(res) != 'try-error') {
      unlink(qry)
      mean_ttol <- getVal(res[[195]])
      median_ttol <- getVal(res[[199]])
    } else {
      mean_ttol <- NA
      median_ttol <- NA
    }
    res <- data.frame(mean_ttol, median_ttol)
    save(res, file=fl)
  }
  .ignrChck(res, sp1, sp2)
  .updateIgnr()
  res
}

findCommonName <- function(ids) {
  # Return the most likely name to be found in timetree
  # (Based on number of name entries in genbank)
  nms <- unlist(lapply(ids, function(x) node_obj[[x]][['nm']][['scientific name']]))
  ids <- ids[!nms %in% .ignr[['nms']]]
  nnms <- unlist(lapply(ids, function(x) length(node_obj[[x]][['nm']])))
  i <- which.max(nnms)
  if(length(i) > 1) {
    i <- sample(i, 1)
  }
  nm <- node_obj[[ids[i]]][['nm']][['scientific name']]
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