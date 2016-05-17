addKid <- function(node_obj, prid, kid) {
  # Recursviely add to kids slot
  i <- which(txids == prid)
  if(length(i) > 0) {
    node_obj[[i]][['kids']] <-
      c(node_obj[[i]][['kids']], kid)
    prid <- node_obj[[i]][['prid']]
    node_obj <- addKid(node_obj, prid, kid)
  }
  node_obj
}

dropNamesLines <- function(names_lines, txids, parallel) {
  # Use vectorisation to remove unnecessary
  # lines from names, to make looping faster
  check <- function(i) {
    txid <- as.numeric(names_lines[[i]][1])
    bool <- txid %in% txids
    data.frame(bool=bool)
  }
  res <- plyr::mdply(.data=data.frame(i=1:length(names_lines)),
                     .fun=check, .parallel=parallel)
  which.is <- res[res[ ,'bool'],'i']
  names_lines[which.is]
}
  

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