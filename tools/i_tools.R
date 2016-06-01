.i_cntr <- -1

iPrnt <- function(i, max_i, print_out=10) {
  # Report progress in a for loop
  # Cats to console the print_out number of iterations
  if(i == max_i) {
    .i_cntr <- -1
  }
  dvs <- i %/% (max_i/print_out)
  if(dvs > .i_cntr) {
    cat(i, ".... ")
    .i_cntr <<- .i_cntr + 1
  }
}
