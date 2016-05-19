# PLOT LF RANKINGS RESULTS

# START
cat(paste0('\nStage `rankings_results` started at [', Sys.time(), ']\n'))

# PARAMETERS
source('parameters.R')

# FUNCTIONS
library(plotly)

# DIRS
input_file <- file.path('4_timetree', "res.RData")

# INPUT
load(input_file)

# IDENTIFY PLOTTABLE RESULTS
rnkngs <-data.frame(ID=cnddts,
                    Time=rep(NA, length(cnddts)),
                    Success=rep(NA, length(cnddts)),
                    Change=rnorm(length(cnddts)),
                    Sci=rep(NA, length(cnddts)),
                    Common=rep("Unknown", length(cnddts)),
                    stringsAsFactors=FALSE)
for(i in 1:length(cnddts)) {
  txid <- cnddts[i]
  tm <- node_obj[[txid]][["cntrst_pe"]]
  if(!is.null(tm)) {
    rnkngs[i, 'Time'] <- log(tm)
  }
  sccss <- node_obj[[txid]][["cntrst_n"]]
  if(!is.null(sccss)) {
    rnkngs[i, 'Success'] <- sccss
  }
  chng <- node_obj[[txid]][["cntrst_chng"]]
  if(!is.null(chng)) {
    rnkngs[i, 'Change'] <- chng
  }
  nms <- node_obj[[txid]][["nm"]]
  rnkngs[i, 'Sci'] <- nms[['scientific name']]
  bool <- grepl("common", names(nms))
  if(any(bool)) {
    rnkngs[i, 'Common'] <- nms[bool][1]
  }
}
rnkngs <- rnkngs[!is.na(rnkngs[['Time']]), ]
p <- ggplot(data=rnkngs, aes(x=Time, y=Success, colour=Change,
                             Sci=Sci, Common=Common)) +
  geom_point()
ggplotly(tooltip=c("Sci", "Common"))
