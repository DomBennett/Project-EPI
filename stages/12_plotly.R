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
input_file <- file.path('6_timetree', "res.RData")

# INPUT
load(input_file)

# IDENTIFY PLOTTABLE RESULTS
rnkngs <-data.frame(ID=cnddts,
                    Time=rep(NA, length(cnddts)),
                    Success=rep(NA, length(cnddts)),
                    Change=rnorm(length(cnddts)),
                    Sci=rep(NA, length(cnddts)),
                    Common=rep("Not found", length(cnddts)),
                    stringsAsFactors=FALSE)
for(i in 1:length(cnddts)) {
  txid <- cnddts[i]
  tm <- node_obj[[txid]][["tmsplt"]]
  if(!is.null(tm)) {
    rnkngs[i, 'Time'] <- log(tm)
  }
  sccss <- node_obj[[txid]][["cntrst_n"]]
  if(!is.null(sccss)) {
    rnkngs[i, 'Success'] <- log(sccss)
  }
  chng <- node_obj[[txid]][["cntrst_chng"]]
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
p <- ggplot(data=rnkngs, aes(x=Time, y=Success, colour=Change)) +
  stat_smooth(method="lm", se=FALSE, colour="red") + geom_point(aes(Sci=Sci, Common=Common))
p <- ggplotly(p, tooltip=c("Sci", "Common"))
plotly_POST(p, sharing = "secret")