# Test whether living fossils habitat and ecology is significantly different
# from non-living fossils using text analysis of IUCN red list

# START
cat(paste0('\nStage `text` started at [', Sys.time(), ']\n'))

# PARAMETERS
cutoff <- -0.1  # highest EPI for a living fossil

# LIBS
source(file.path('tools', 'text_tools.R'))

# DIRS
output_dir <- 'text'
cache_dir <- 'iucn_data'
if (!file.exists(output_dir)) {
  dir.create(output_dir)
}
if (!file.exists(cache_dir)) {
  dir.create(cache_dir)
}

# GET IUCN TOKEN
source('token.R')

# LOAD EPI DATA
metrics <- read.csv(file.path('other', 'old', 'livingfossils.csv'),
                 stringsAsFactors=FALSE)
metrics$epi <- ((metrics[['change']] + metrics[['n']])/2) - 
  metrics[['time']]
# limit to species
metrics <- metrics[metrics[['n']] == 1, ]
# remove _
metrics[['clade']] <- gsub("_", " ", metrics[['clade']])

# PARITION
lfs <- metrics[metrics[['epi']] < cutoff, ]
lfs$clade

# SEARCH IUCN
lfs_nrrtvs <- list()
for(i in 1:nrow(lfs)) {
  nm <- lfs[i, 'clade']
  cat('Searching [', nm, '] ....\n')
  lfs_nrrtvs[[nm]] <- getIUCNNrrtv(nm, token)
}

# GET HABITAT AND ECOLOGY DISTS
lfs_hes <- getHbttEclgy(lfs_nrrtvs)
obs <- calcStrDst(lfs_hes, mthd='cosine')

# ITERATE
itrtns <- 99
null_dsts <- obs
nulls_nrrtvs <- list()
for(itrtn in 1:itrtns) {
  cat('Iteration [', itrtn, '] ....\n')
  rnds <- sample(1:nrow(metrics), size=nrow(lfs))
  nulls <- metrics[rnds, ]
  for(i in 1:nrow(nulls)) {
    nm <- nulls[i, 'clade']
    nulls_nrrtvs[[nm]] <- getIUCNNrrtv(nm, token)
  }
  nulls_hes <- getHbttEclgy(nulls_nrrtvs[nulls[['clade']]])
  res <- calcStrDst(nulls_hes, mthd='cosine')
  null_dsts <- rbind(null_dsts, res)
}
null_dsts <- null_dsts[-1, ]

# TEST SIGNIFICANCE
p_val <- sum(null_dsts$mean <= obs$mean)/nrow(null_dsts)
hist(null_dsts$mean, main=paste0('P = ', signif(p_val, 3)))
abline(v=obs$mean, col='red')

# WORD FREQUENCIES
obs <- getWrdFrq(lfs_hes, min_freq=2)  # observed freqs in living fossils
itrtns <- 999
null_dsts <- list()  # null freqs
for(itrtn in 1:itrtns) {
  cat('Iteration [', itrtn, '] ....\n')
  rnds <- sample(1:nrow(metrics), size=nrow(lfs))
  nulls <- metrics[rnds, ]
  for(i in 1:nrow(nulls)) {
    nm <- nulls[i, 'clade']
    nulls_nrrtvs[[nm]] <- getIUCNNrrtv(nm, token)
  }
  nulls_hes <- getHbttEclgy(nulls_nrrtvs[nulls[['clade']]])
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
cat(paste0('\nStage `text` finished at [', Sys.time(), ']\n'))