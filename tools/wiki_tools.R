# DIRS
cache_dir <- "caches"
if(!file.exists(cache_dir)) {
  dir.create(cache_dir)
}
cache_dir <- file.path("caches", "wiki")
if(!file.exists(cache_dir)) {
  dir.create(cache_dir)
}

getWikiUrl <- function(ncbi_tax_page) {
  # Retrieve the Wiki URL from an ncbi_tax_page
  res <- NA
  for(ln in ncbi_tax_page) {
    if(grepl('wikipedia', ln)) {
      strt <- regexec('href=', ln)
      end <- regexec('\\sref=', ln)
      res <- substr(ln, strt[[1]]+6, end[[1]]-2)
    }
  }
  res
}

getWikiLFMention <- function(txid) {
  # T/F -- does wiki article mention living fossil?
  # First search NCBI taxonomy
  fl <- file.path(cache_dir, paste0(txid, ".RData"))
  if(file.exists(fl)) {
    load(fl)
    return(res)
  }
  res <- NA
  url <- "http://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Info&id="
  qry <- paste0(url, txid)
  ncbi_page <- searchURL(qry)
  wiki_url <- getWikiUrl(ncbi_page)
  if(!is.na(wiki_url)) {
    res <- FALSE
    # ensure https
    if(grepl('http:', wiki_url)) {
      wiki_url <- sub('http:', 'https:', wiki_url)
    }
    wiki_page <- searchURL(wiki_url)
    for(ln in wiki_page) {
      if(grepl('living fossil', ln, ignore.case=TRUE)) {
        res <- TRUE
      }
    }
  }
  save(res, file=fl)
  res
}