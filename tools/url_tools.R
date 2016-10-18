searchURL <- function(qry_url, timeout=10) {
  # Safely search for a url, re-try if not accessible
  attmpts <- 1
  res <- NA
  while(attmpts <= length(wt)) {
    Sys.sleep(wt[attmpts])
    # why does R insist on making this the ugliest and most complicated of code?!?!?
    suppressWarnings(try_err <- try(expr={
      cnnctn <- R.utils::withTimeout(expr={
        url(description=qry_url, open='r')
      }, timeout=timeout, onTimeout='error')
    }, silent=TRUE))
    if(grepl('reached elapsed time limit', try_err[[1]])) {
      return(NA)
    }
    if(exists('cnnctn') && !is.null(cnnctn)) {
      is_open <- try(expr={
        isOpen(cnnctn)
      }, silent=TRUE)
    }
    if(class(is_open) != 'try-error' && is_open) {
      if(attmpts > 1) {
        cat('---- Reconnected ----')
      }
      res <- readLines(con=cnnctn)
      close(cnnctn)
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