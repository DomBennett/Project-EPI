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
    cat('----- Failed to reach, [', url, '] trying again in [',
        wt[attmpts], 's] -----\n', sep='')
    attmpts <- attmpts + 1
  }
  if(attmpts > length(wt)) {
    cat("---- Max attempts made, server must be down ----")
  }
  res
}