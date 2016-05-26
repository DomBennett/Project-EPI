getGrpTxids <- function(txids, grp) {
  # Return txids for certain group
  .get <- function(txid, bool) {
    prid <- node_obj[[txid]][['prid']]
    if(prid != grpid) {
      bool <- (txids == prid) + bool
      bool <- .get(prid, bool)
    }
    bool
  }
  fl <- file.path('2_node_obj', paste0(grp, "_txids.RData"))
  if(file.exists(fl)) {
    load(fl)
    return(txids)
  }
  ignr <- NULL
  nms <- names(anlyss_grps)[names(anlyss_grps) != grp]
  for(nm in nms) {
    ignr <- c(ignr, node_obj[[anlyss_grps[[nm]]]][['kids']])
  }
  txids <- txids[!txids %in% ignr]
  grpid <- anlyss_grps[[grp]]
  bool <- txids == grpid
  kids <- node_obj[[grpid]][['kids']]
  for(txid in kids) {
    bool <- .get(txid, bool)
  }
  txids <- c(kids, txids[bool > 0])
  save(txids, file=fl)
  txids
}

getSppTxids <- function(txids) {
  # return all txids for species
  txids[unlist(lapply(txids, function(x) node_obj[[x]][['rank']] == "species"))]
}

addKid <- function(prid, kid) {
  # Recursviely add to kids slot
  if(prid %in% txids) {
    node_obj[[prid]][['kids']] <-
      c(node_obj[[prid]][['kids']], kid)
    prid <- node_obj[[prid]][['prid']]
    addKid(prid, kid)
  }
}

# dropNamesLines <- function(names_lines, txids, parallel) {
#   # Use vectorisation to remove unnecessary
#   # lines from names, to make looping faster
#   check <- function(i) {
#     txid <- as.numeric(names_lines[[i]][1])
#     bool <- txid %in% txids
#     data.frame(bool=bool)
#   }
#   res <- plyr::mdply(.data=data.frame(i=1:length(names_lines)),
#                      .fun=check, .parallel=parallel)
#   which.is <- res[res[ ,'bool'],'i']
#   names_lines[which.is]
# }

# fetchName <- function(txid) {
#   # Return scientific name using NCBI taxonomic ID
#   fp <- file.path("ncbi_cache", paste0(txid, ".RData"))
#   if(file.exists(fp)) {
#     load(fp)
#   } else {
#     res <- rentrez::entrez_fetch(db="taxonomy", id=txid, rettype="xml")
#     f1 <- regexpr("<ScientificName>", res)
#     f2 <- regexpr("</ScientificName>", res)
#     nm <- substr(res, f1+16, f2-1)
#     save(nm, file=fp)
#   }
#   nm
# }