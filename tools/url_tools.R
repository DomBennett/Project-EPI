searchURL <- function(url, site) {
  # Safely search for a url, re-try if not accessible
  attmpts <- 1
  while(attmpts <= length(wt)) {
    Sys.sleep(wt[attmpts])
    res <- suppressWarnings(try(expr=readLines(url),
                                silent=TRUE))
    if(class(res) != 'try-error') {
      if(attmpts > 1) {
        cat('---- Reconnected ----')
      }
      unlink(url)
      break
    }
    attmpts <- attmpts + 1
    cat('----- Failed to reach, [', site, '] trying again in [',
        wt[attmpts], 's] -----\n', sep='')
  }
  if(attmpts > length(wt)) {
    stop("---- Max attempts made, server must be down ----")
  }
  res
}