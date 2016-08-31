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
  load(tree_file)
  if('tree' %in% ls()) {
    return(tree)
  }
  tree <- trees[[i]]
  rm(trees)
  updateTree(tree)
}

getNtrees <- function(tree_file) {
  load(tree_file)
  if('tree' %in% ls()) {
    return(1)
  }
  trees['ntrees']
}