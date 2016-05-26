# TEST IF LIVING FOSSILS SHARE MORE

# START
cat(paste0('\nStage `analysis` started at [', Sys.time(), ']\n'))

# FUNCTIONS
source(file.path('tools', 'analysis_tools.R'))
source(file.path('tools', 'node_obj_tools.R'))

# PARAMETERS
source('parameters.R')
token <- getToken()

# DIRS
output_dir <- '7_analysis'
if (!file.exists(output_dir)) {
  dir.create(output_dir)
}
input_file <- file.path("6_epi", "res.RData")
output_file <- file.path("7_analysis", "res.RData")

# INPUT
load(input_file)

# GET GROUP IDS
txids <- ls(node_obj)
txids <- getGrpTxids(txids, grp="mammals")
spp <- getSppTxids(txids)

# GET LIVING FOSSILS
epi$ntlg_indx <- log(epi$cntrst_n) - log(epi$tmsplt)
epi <- epi[order(epi[['ntlg_indx']], decreasing=FALSE), ]
lf_txids <- epi[['txid']][epi[['txid']] %in% txids][1:40]
lf_data <- vector("list", length=length(lf_txids))
names(lf_data) <- lf_txids

# SEARCH IUCN FOR LIVING FOSSILS
for(i in 1:length(lf_data)) {
  txid <- names(lf_data)[i]
  nms <- getKidNms(txid)
  res <- NULL
  for(j in 1:length(nms)) {
    cat('Searching [', nms[j], '] ....\n', sep="")
    nrrtv <- getIUCNNrrtv(nms[j], token)
    if(class(nrrtv) == "list" && length(nrrtv[['result']]) > 0 &&
       !is.null(nrrtv[['result']][[1]][['habitat']])) {
      res <- c(res, nrrtv[['result']][[1]][['habitat']])
    }
  }
  if(!is.null(res)) {
    lf_data[[txid]][['hbbts']] <- res
  }
}
lf_data <- lf_data[unlist(lapply(lf_data, function(x) length(x[[1]]) > 0))]

# GET STRING DISTANCES USING PERMUTATION
itrntns <- 99
lf_dst <- matrix(nrow=itrntns, ncol=5)
for(i in 1:itrntns) {
  txts <- vector(length=length(lf_data))
  for(j in 1:length(lf_data)) {
    n <- length(lf_data[[j]][[1]])
    txts[j] <- lf_data[[j]][[1]][[sample(1:n, 1)]]
  }
  res <- calcStrDst(txts, mthd='cosine')
  lf_dst[i, ] <- as.numeric(res)
}
colnames(lf_dst) <- colnames(res)

# GET HBBTS FOR RANDOM SPP
itrtns <- 99
nulls_hbbts <- list()
for(itrtn in 1:itrtns) {
  cat('Iteration [', itrtn, '] ....\n', sep="")
  nulls_hbbts[[itrtn]] <- list()
  cc <- 0
  while(cc < length(lf_data)) {
    null_spp <- sample(spp, size=1)
    nm <- node_obj[[null_spp]][['nm']][['scientific name']]
    nrrtv <- getIUCNNrrtv(nm, token)
    if(class(nrrtv) == "list" && length(nrrtv[['result']]) > 0 &&
       !is.null(nrrtv[['result']][[1]][['habitat']])) {
      nulls_hbbts[[itrtn]][[nm]] <-
        nrrtv[['result']][[1]][['habitat']]
      cc <- cc + 1
    }
  }
}

# CALCULATING STR DST OF NULL
null_dsts <- matrix(nrow=itrntns, ncol=5)
for(i in 1:length(nulls_hbbts)) {
  txts <- as.character(unlist(nulls_hbbts[[i]]))
  res <- calcStrDst(txts, mthd="cosine")
  null_dsts[i, ] <- as.numeric(res)
}
colnames(null_dsts) <- colnames(res)

# TEST SIGNIFICANCE
obs_mean <- mean(lf_dst[ ,"median"])
p_val <- sum(null_dsts[ ,"median"] <= obs_mean)/nrow(null_dsts)
hist(null_dsts[ ,"median"], main=paste0('P = ', signif(p_val, 3)))
abline(v=obs_mean, col='red')

# WORD FREQUENCIES
obs <- getWrdFrq(lfs_hes, min_freq=2)  # observed freqs in living fossils
itrtns <- 999
null_dsts <- list()  # null freqs
for(itrtn in 1:itrtns) {
  cat('Iteration [', itrtn, '] ....\n')
  rnds <- sample(1:nrow(metrics), size=nrow(lfs))
  nulls <- metrics[rnds, ]
  for(i in 1:nrow(nulls)) {
    nm <- nulls[i, 'node.label']
    nulls_nrrtvs[[nm]] <- getIUCNNrrtv(nm, token)
  }
  nulls_hes <- getHbttEclgy(nulls_nrrtvs[nulls[['node.label']]])
  res <- getWrdFrq(nulls_hes, min_freq=2)
  null_dsts <- c(null_dsts, list(res))
}
# how many more times does a term appear in observed?
wrd_res <- data.frame(wrd=NA, obs_freq=NA, exp_freq=NA,
                      p_val=NA, z_score=NA)
for(i in 1:length(obs)) {
  wrd <- names(obs)[i]
  nd <- unlist(lapply(null_dsts, function(x) {
    if(wrd %in% names(x)) {
      res <- x[[wrd]]
    } else {
      res <- 0
    }
    res
  }))
  wrd_res[i, 'wrd'] <- wrd
  wrd_res[i, 'obs_freq'] <- obs[[wrd]]
  wrd_res[i, 'exp_freq'] <- mean(nd)
  wrd_res[i, 'z_score'] <- (obs[[wrd]] - mean(nd))/sd(nd)
  wrd_res[i, 'p_val'] <- sum(nd >= obs[[wrd]])/length(nd)
}
temp_res <- wrd_res[wrd_res$p_val < 0.05, ]
temp_res <- temp_res[temp_res$z_score < 5, ]
pdf(file.path('figures', paste0('mammals', '_wordcloud.pdf')), w=14, h=14)
wordcloud(temp_res$wrd, temp_res$z_score,
          colors=brewer.pal(8, 'Dark2'), max.words=45)
dev.off()


wrd <- "solitary"

for(i in 1:length(lfs_hes)) {
  res <- gregexpr(wrd, lfs_hes[i])[[1]]
  if(res[[1]] != -1) {
    start <- res[[1]] - 10
    end <- res[[1]] + attr(res, 'match.length') + 10
    txt <- substr(lfs_hes[[i]], start=start, stop=end)
    cat(names(lfs_hes)[i], ": ", txt, "\n", sep="")
  }
}

# END
cat(paste0('\nStage `analysis` finished at [', Sys.time(), ']\n'))