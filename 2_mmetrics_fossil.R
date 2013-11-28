## 18/11/2013
## D.J. Bennett
## Using the O'Leary 2013 morpho data

## Libraries
source(file.path("Functions", "EcoDataTools.R"))
require("reshape2")

## Dirs

input.dir <- "0_data"
output.dir <- "2_mmetrics"

## Data
data <- read.nexus.data(file.path(input.dir, "test_nexus.nex"))

file <- file.path(input.dir, "test_nexus.nex")

X <- scan(file = file, what = character(), sep = "\n", quiet = TRUE, 
          comment.char = "[", strip.white = TRUE)

# line 47, fist element of the matrix
Xj <- X[47]
# this looks for bloacks of white space that are followed and preceded by characters
unlist(strsplit(Xj, "(?<=\\S)(\\s+)(?=\\S)", perl = TRUE))
# at the beginning of the function test to see if the names are enclosed by quotes e.g. ' "
# if this is the case, find all names and convert them
grep("(\"|').*(\"|')", Xj)


gsub("('|\")", "", sub("\\s+", "_", "'species name'"))

# read in all the names into a list
# match each name to matrix, read matrix into list
# convert names into R friendly.

findTaxaline <- function(x) {
  taxaline <- rep(NA, 2)
  taxlabel.section <- FALSE
  for (i in 1:NROW(x)) {
    if (grepl("taxlabels", x[i], ignore.case = TRUE)) {
      taxaline[1] <- as.numeric(i) + 1
      taxlabel.section <- TRUE
    }
    if (grepl(";", x[i], ignore.case = TRUE) & taxlabel.section) {
      taxaline[2] <- as.numeric(i) - 1
      break
    }
  }
  return(taxaline)
}

taxaline <- findTaxaline(X)

taxa <- vector()
for (j in taxaline[1]:taxaline[2]) {
  #Xj <- trim.semicolon(X[j])
  taxa <- c(taxa, X[j])
}