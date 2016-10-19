searchURL <- function(qry_url, timeout=10) {
  # Safely search for a url, re-try if not accessible
  attmpts <- 1
  res <- NA
  while(attmpts <= length(wt)) {
    Sys.sleep(wt[attmpts])
    # why does R insist on making this the ugliest and most complicated of code?!?!?
    suppressWarnings(try_err <- try(expr={
      res <- R.utils::withTimeout(expr={
        readLines(con=qry_url)
      }, timeout=timeout, onTimeout='error')
    }, silent=TRUE))
    if(!grepl('cannot open the connection', try_err[[1]])) {
      if(attmpts > 1) {
        cat('---- Reconnected ----')
      }
      unlink(res)
      break
    }
    cat('----- Failed to reach, [', qry_url, '] trying again in [',
        wt[attmpts], 's] -----\n', sep='')
    attmpts <- attmpts + 1
  }
  if(attmpts > length(wt)) {
    cat("---- Max attempts made, server must be down ----")
  }
  res
}
