searchURL <- function(qry_url, site) {
  # Safely search for a url, re-try if not accessible
  attmpts <- 1
  res <- NA
  while(attmpts <= length(wt)) {
    Sys.sleep(wt[attmpts])
    suppressWarnings(try_err <- try(expr={
      cnnctn <- R.utils::withTimeout(expr={
        url(description=qry_url, open='r')
      }, timeout=10, onTimeout='error')
    }, silent=TRUE))
    if(grepl('reached elapsed time limit', try_err[[1]])) {
      return(NA)
    }
    if(class(try_err) != 'try-error') {
      if(file.exists('cnnctn') && !is.null(cnnctn) && isOpen(cnnctn)) {
        if(attmpts > 1) {
          cat('---- Reconnected ----')
        }
        res <- readLines(con=cnnctn)
        close(cnnctn)
      }
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

searchURL(qry_url='http://www.timetree.org/ajax/pairwise/1839748/1031332')
rm(try_err, res, cnnctn)
