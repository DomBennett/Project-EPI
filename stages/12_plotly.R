# PLOT LF RANKINGS RESULTS

# FUNCTION
simpleCap <- function(x) {
  # http://stackoverflow.com/questions/6364783/capitalize-the-first-letter-of-both-words-in-a-two-word-string
  s <- strsplit(x, " ")[[1]]
  paste(toupper(substring(s, 1,1)), substring(s, 2),
        sep="", collapse=" ")
}

# START
cat(paste0('\nStage `plotly` started at [', Sys.time(), ']\n'))

# PARAMETERS
source('parameters.R')

# FUNCTIONS
library(plotly)

# DIRS
input_file <- file.path('8_epi', "res.RData")

# INPUT
load(input_file)

# IDENTIFY PLOTTABLE RESULTS
txids <- epi[['txid']]
rnkngs <-data.frame(ID=txids,
                    Time=rep(NA, length(txids)),
                    Success=rep(NA, length(txids)),
                    Change=rep(NA, length(txids)),
                    Sci=rep(NA, length(txids)),
                    Common=rep("Not found", length(txids)),
                    stringsAsFactors=FALSE)
for(i in 1:length(txids)) {
  txid <- txids[i]
  tm <- node_obj[[txid]][["tmsplt"]]
  if(!is.null(tm)) {
    rnkngs[i, 'Time'] <- log(tm)
  }
  sccss <- node_obj[[txid]][["cntrst_n"]]
  if(!is.null(sccss)) {
    rnkngs[i, 'Success'] <- log(sccss)
  }
  chng <- node_obj[[txid]][["chng"]]
  if(!is.null(chng)) {
    rnkngs[i, 'Change'] <- chng
  }
  nms <- node_obj[[txid]][["nm"]]
  rnkngs[i, 'Sci'] <- nms[['scientific name']]
  bool <- grepl("common", names(nms))
  if(any(bool)) {
    rnkngs[i, 'Common'] <- simpleCap(nms[bool][1])
  }
}
rnkngs <- rnkngs[!is.na(rnkngs[['Time']]), ]
rnkngs <- rnkngs[!is.na(rnkngs[['Change']]), ]
p <- ggplot(data=rnkngs, aes(x=Time, y=Success, colour=Change)) +
  geom_point(aes(Sci=Sci, Common=Common))
p <- ggplotly(p, tooltip=c("Sci", "Common"))
plotly_POST(p, sharing = "secret")