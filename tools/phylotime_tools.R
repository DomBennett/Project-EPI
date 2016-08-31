assgnWMean <- function(val, nm) {
  # assign a value to node_obj multiple times
  # if value already present, work out mean
  if(is.null(node_obj[[txid]][[nm]])) {
    node_obj[[txid]][[nm]] <- val
  } else {
    # otherwise get mean
    node_obj[[txid]][[nm]] <-
      (node_obj[[txid]][[nm]] + val)/2
  }
}

getTree <- function(i, tree_file) {
  # looks up ith tree in tree_file
  # trees are saved as RData to save processing
  # if tree is not updated, tree is updated
  load(file=tree_file)
  if('tree' %in% ls()) {
    return(tree)
  }
  tree <- trees[[i]]
  rm(trees)
  tree
}

getNtrees <- function(tree_file) {
  load(file=tree_file)
  if('tree' %in% ls()) {
    return(1)
  }
  trees['ntrees']
}

calcFrPrp2 <- function(tree, tids, progress="none") {
  # treeman function without bigmemory
  .calc <- function(i) {
    id <- tree@all[i]
    spn <- getNdSlt(tree, "spn", id)
    kids <- getNdKids(tree, id)
    if(length(kids) == 0) {
      spn_shres[i, id] <<- spn
    } else {
      spn_shre <- spn/length(kids)
      spn_shres[i, kids] <<- spn_shre
    }
  }
  spn_shres <- matrix(0, ncol=tree@ntips, nrow=tree@nall)
  colnames(spn_shres) <- tree@tips
  plyr::m_ply(.data=data.frame(i=1:tree@nall), .fun = .calc,
              .progress=progress)
  colSums(spn_shres[, tids])
}