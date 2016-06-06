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