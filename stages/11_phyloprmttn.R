# PERMUTATION OF PD

# START
cat(paste0('\nStage `phyloprmttn` started at [', Sys.time(), ']\n'))

# PARAMETERS
source('parameters.R')

# FUNCTIONS
library(treeman)
source(file.path('tools', 'phyloprmttn_tools.R'))
source(file.path('tools', 'node_obj_tools.R'))
source(file.path('tools', 'i_tools.R'))

# DIRS
output_dir <- '11_phyloprmttn'
if(!file.exists(output_dir)) {
  dir.create(output_dir)
}
tree_dir <- file.path("0_data", "trees")
input_file <- file.path('7_epi', "res.RData")

# INPUT
load(input_file)
tree_files <- list.files(tree_dir)

# LOOP THROUGH TREE FILES
cat("Looping though group trees ....\n")
for(tree_file in tree_files) {
  # INPUT
  grp <- sub("\\.tre", "", tree_file)
  cat('    Reading in [', grp, '] ....\n', sep="")
  txids <- ls(node_obj)
  txids <- getGrpTxids(txids, grp=grp)
  trees <- readTree(file.path(tree_dir, tree_file))
  cat("    Done.\n")
  if(class(trees) == "TreeMan") {
    ntrees <- 1
    trees <- list(trees)
  } else {
    ntrees <- trees['ntrees']
  }
  
  # GET LIVING FOSSILS
  cat("    Finding [", nlfs, "] top most living fossil clades ....",
      sep="")
  epi <- epi[order(epi[['nt_indx_lg']], decreasing=FALSE), ]
  lf_txids <- epi[['txid']][epi[['txid']] %in% txids][1:nlfs]
  lf_txids <- lf_txids[!is.na(lf_txids)]
  if(length(lf_txids) < 10) {
    cat("    \nToo few living fossils for this group!\n")
    next
  }
  cat("Done.\n")
  
  # OBSERVED PD
  cat("    Calculating PD for living fossils ....")
  lf_nms <- unlist(lapply(lf_txids, function(x) node_obj[[x]][['nm']][['scientific name']]))
  lf_nms <- gsub(" ", "_", lf_nms)
  obs_pds <- vector(length=ntrees)
  for(i in 1:ntrees) {
    lf_ids <- NULL
    for(j in 1:length(lf_txids)) {
      if(!lf_nms[j] %in% trees[[i]]['tips']) {
        # if not in tree, use children to find internal node
        kids <- node_obj[[lf_txids[j]]][['kids']]
        kids_nms <- unlist(lapply(kids,
                                  function(x) node_obj[[x]][['nm']][['scientific name']]))
        kids_nms <- gsub(" ", "_", kids_nms)
        tval <- sum(kids_nms %in% tree['tips'])/length(kids_nms)
        if(length(kids_nms) > 0 && tval > 0.5) {
          kids_nms <- kids_nms[kids_nms %in% tree['tips']]
          lf_ids <- c(lf_ids, getPrnt(trees[[i]], kids_nms))
        }
      } else {
        lf_ids <- c(lf_ids, lf_nms[j])
      }
    }
    dst_mtrx <- calcDstMtrx(trees[[i]], ids=lf_ids)
    obs_pds[i] <- mean(dst_mtrx)
  }
  obs_pd <- mean(obs_pds)
  cat("Done, [", length(lf_ids), '] in tree.\n', sep="")
  
  # ITERATION
  cat("    Iterating for null .... ")
  smpl_sz <- length(lf_ids)
  null_pds <- vector(length=nitrtns)
  for(i in 1:nitrtns) {
    iPrnt(i, nitrtns)
    null_ids <- sample(tree['all'], size=smpl_sz)
    null_pds[i] <- meanPD(null_ids, trees, ntrees)
  }
  cat("Done.\n")
  
  # STATS
  cat("    Outputting stats .... ")
  null_mean <- mean(null_pds)
  null_sd <- sd(null_pds)
  p_val <- sum(null_mean <= obs_pd)/length(null_pds)
  z_score <- (obs_pd - null_mean)/null_sd
  pdf(file.path(output_dir, paste0(grp, ".pdf")))
  hist(null_pds,
       main=paste0('P = ', signif(p_val, 3)),
       xlab=paste0(grp, " -- PD"))
  abline(v=obs_pd, col='red')
  dev.off()
  res <- list(obs_pd=obs_pd, null_mean=null_mean,
              null_sd=null_sd, z_score=z_score,
              p_val=p_val, smpl_sz=smpl_sz)
  save(res, file=file.path(output_dir, paste0(grp, '.RData')))
  cat("Done.\n")
}