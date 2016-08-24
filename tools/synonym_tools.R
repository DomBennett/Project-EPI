library(ggplot2)

cleanWrds <- function(txt) {
  # Take text and turn into list of unique words
  wrds <- strsplit(txt, " ")[[1]]
  wrds <- tolower(wrds)
  wrds <- tm::removePunctuation(wrds)
  wrds <- tm::removeNumbers(wrds)
  wrds <- gsub("[[:space:]]", "", wrds)
  wrds <- unique(wrds)
  wrds <- wrds[sapply(wrds, nchar) > min_nchar]
  wrds
}

senseScore <- function(slt_nm, pstvs, ngtvs) {
  # Calculate the sense score
  # A measure of whether descriptions match a given meaning
  # based on synonyms (pstvs) and antonyms (ngtvs)
  pstv_counts <- rep(0, length(pstvs))
  ngtv_counts <- rep(0, length(ngtvs))
  names(pstv_counts) <- pstvs
  names(ngtv_counts) <- ngtvs
  epi[[slt_nm]] <<- NA
  for(i in wwrds) {
    pscr <- 0
    for(pstv in pstvs) {
      wc <- any(grepl(pstv, wrds[[i]]))
      pscr <- pscr + wc
      pstv_counts[[pstv]] <- pstv_counts[[pstv]] + wc
    }
    nscr <- 0
    for(ngtv in ngtvs) {
      wc <- any(grepl(ngtv, wrds[[i]]))
      nscr <- nscr + wc
      ngtv_counts[[ngtv]] <- ngtv_counts[[ngtv]] + wc
    }
    epi[[slt_nm]][i] <<- (pscr > 0) - (nscr > 0)
  }
  pstv_counts <- data.frame(wrd=names(pstv_counts),
                            count=pstv_counts,
                            type="synonym")
  ngtv_counts <- data.frame(wrd=names(ngtv_counts),
                            count=ngtv_counts,
                            type="antonym")
  ordrd_levels <- c(as.character(pstv_counts$wrd[
    order(pstv_counts$count, decreasing=FALSE)]),
    as.character(ngtv_counts$wrd[
      order(ngtv_counts$count, decreasing=TRUE)]))
  wrd_counts <- rbind(pstv_counts, ngtv_counts)
  wrd_counts$wrd <- factor(wrd_counts$wrd, levels=ordrd_levels)
  wrd_counts
}

getExampleTextFrmI <- function(epi, i) {
  nms <- getKidNms(epi[i,'txid'])
  res <- NULL
  for(nm in nms) {
    nrrtv <- getIUCNNrrtv(nm, token)
    if(class(nrrtv) == "list" && length(nrrtv[['result']]) > 0 &&
       !is.null(nrrtv[['result']][[1]][['habitat']])) {
      nrrtv <- nrrtv[['result']][[1]][['habitat']]
      nrrtv <- gsub("<.*?>", "", nrrtv)  # remove html tags
      res <- c(res, nrrtv)
    }
  }
  cat("EPI: [", epi[i,'nm']," @ ", i, "]\n", sep="")
  res
}

getExampleTextFrmWrd <- function(epi, pttrn) {
  # Get example text from a word a pattern
  bool <- sapply(wrds, function(x) any(grepl(pttrn, x)))
  i <- sample(which(bool), 1)
  getExampleTextFrmI(epi, i)
}

ggBoxplot <- function(score, x="pepi") {
  p_data <- data.frame(score=epi[[score]],
                       epi=epi[[x]])
  p_data <- na.omit(p_data)
  p <- ggplot(p_data, aes(factor(score), epi)) +
    geom_boxplot()
  p + theme_bw(fill="cornflowerblue")
}

ggViolin <- function(score, x="pepi") {
  p_data <- data.frame(score=epi[[score]],
                       epi=epi[[x]])
  p_data <- p_data[p_data$score %in% c(-1, 1), ]
  p_data$score[p_data$score == -1] <- "antonym"
  p_data$score[p_data$score == 1] <- "synonym"
  p_data$score <- factor(p_data$score, levels=c('synonym',
                                                "antonym"))
  p <- ggplot(p_data, aes(score, epi, fill=score)) + geom_violin()
  p <- p + theme_bw() + theme(axis.text.x=element_blank(),
                              legend.position="none")
  p
}

ggWrdBars <- function(wrd_counts) {
  ggplot(wrd_counts, aes(x=wrd, y=count, fill=type)) +
    coord_flip() + geom_bar(stat="identity") +
    ylab("Freq.") + xlab("") + theme_bw() +
    theme(legend.title=element_blank())
}