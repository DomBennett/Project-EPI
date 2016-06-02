# FOLDERS
cache_dir <- "caches"
if(!file.exists(cache_dir)) {
  dir.create(cache_dir)
}
cache_dir <- file.path("caches", "iucn")
if(!file.exists(cache_dir)) {
  dir.create(cache_dir)
}
nrrtv_dir <- file.path("caches", "iucn", "narratives")
if(!file.exists(nrrtv_dir)) {
  dir.create(nrrtv_dir)
}
hbbts_dir <- file.path("caches", "iucn", "habitats")
if(!file.exists(hbbts_dir)) {
  dir.create(hbbts_dir)
}
cntrs_dir <- file.path("caches", "iucn", "countries")
if(!file.exists(cntrs_dir)) {
  dir.create(cntrs_dir)
}
cat_dir <- file.path("caches", "iucn", "category")
if(!file.exists(cat_dir)) {
  dir.create(cat_dir)
}

# FUNCTIONS
getToken <- function() {
  if(file.exists("iucn_token.R")) {
    source("iucn_token.R")
  } else {
    msg <- "No token found!
Apply for one here: http://apiv3.iucnredlist.org/api/v3/token
Save the token in an R script called `iucn_token.R` in the working dir:
    `token <- [USER TOKEN]`\n"
    stop(msg)
  }
  token
}

cleanNm <- function(nm) {
  # Return name that is html safe and API safe
  nm <- sub("\\/.*", "", nm)
  nm <- gsub("[0-9]", "", nm)
  nm <- sub("sp\\.", "", nm)
  nm <- paste0(toupper(substring(nm, 1,1)),
               tolower(substring(nm, 2)), collapse="")
  nm
}

getIUCNCat <- function(nm, token) {
  # Get category data for species from IUCN API
  # first make sure nm is html safe
  nm <- cleanNm(nm)
  # second check if not already downloaded
  fl <- file.path(cat_dir, paste0(gsub(" ", "_", nm), '.RData'))
  if(file.exists(fl)) {
    load(fl)
  } else {
    url <- paste0("http://apiv3.iucnredlist.org/api/v3/species/", nm,"?token=", token)
    res <- MoreTreeTools:::.safeFromJSON(url)
    save(res, file=fl)
  }
  res
}

getIUCNNrrtv <- function(nm, token) {
  # Get narrative data for species from IUCN API
  # first make sure nm is html safe
  nm <- cleanNm(nm)
  # second check if not already downloaded
  fl <- file.path(nrrtv_dir, paste0(gsub(" ", "_", nm), '.RData'))
  if(file.exists(fl)) {
    load(fl)
  } else {
    url <- paste0("http://apiv3.iucnredlist.org/api/v3/species/narrative/", nm,"?token=", token)
    res <- MoreTreeTools:::.safeFromJSON(url)
    save(res, file=fl)
  }
  res
}

getIUCNHbbts <- function(nm, token) {
  # Get habitat data for species from IUCN API
  # first make sure nm is html safe
  nm <- cleanNm(nm)
  # second check if not already downloaded
  fl <- file.path(hbbts_dir, paste0(gsub(" ", "_", nm), '.RData'))
  if(file.exists(fl)) {
    load(fl)
  } else {
    url <- paste0("http://apiv3.iucnredlist.org/api/v3/habitats/species/name/",
                  nm,"?token=", token)
    res <- MoreTreeTools:::.safeFromJSON(url)
    save(res, file=fl)
  }
  res
}

getIUCNCntrs <- function(nm, token) {
  # Get country data for species from IUCN API
  # first make sure nm is html safe
  nm <- cleanNm(nm)
  # second check if not already downloaded
  fl <- file.path(cntrs_dir, paste0(gsub(" ", "_", nm), '.RData'))
  if(file.exists(fl)) {
    load(fl)
  } else {
    url <- paste0("http://apiv3.iucnredlist.org/api/v3/species/countries/name/",
                  nm,"?token=", token)
    res <- MoreTreeTools:::.safeFromJSON(url)
    save(res, file=fl)
  }
  res
}

getKidNms <- function(txid) {
  # Return the names of the children of a txid
  kds <- node_obj[[txid]][['kids']]
  if(kds[1] == "none") {
    return(node_obj[[txid]][['nm']][['scientific name']])
  }
  nms <- vector(length=length(kds))
  for(i in 1:length(kds)) {
    nms[i] <- node_obj[[kds[i]]][['nm']][['scientific name']]
  }
  nms
}