readNexusData <- function(file) {
  # internal functions
  getBlock <- function (block.header, data) {
    start.line <- 0
    for (i in 1:length(data)) {
      if (grepl(block.header, data[i])) {
        start.line <- i + 1
        break
      }
    }
    end.line <- 0
    for (i in start.line:length(data)) {
      if (grepl("ENDBLOCK", data[i])) {
        end.line <- i - 1
        break
      }
    }
    return (data[start.line:end.line])
  }
  getBlockData <- function (label, block) {
    start.line <- 0
    for (i in 1:length(block)) {
      if (grepl(label, block[i])) {
        start.line <- i
        break
      }
    }
    if (start.line == 0) {
      stop (paste0 ("Block data error: no label ", label))
    }
    end.line <- 0
    for (i in start.line:length(block)) {
      if (grepl(';$', block[i])) {
        end.line <- i
        break
      }
    }
    if (start.line == end.line) {
      block.data <- strsplit (block[start.line], '\\s')[[1]][-1]
      res <- matrix (ncol = 2, nrow = length(block.data))
      colnames(res) <- c("name", "value")
      for (i in 1:length(block.data)) {
        element <- strsplit(block.data[i], "=")[[1]]
        # strip excessive characters caused by quotes and semi-colon at end of line
        element[2] <- gsub('\\"|;', "", element[2])
        res[i,1] <- element[1]
        res[i,2] <- element[2]
      }
    } else {
      res <- vector()
      for (i in (start.line+1):(end.line-1)) {
        res <- c(res, block[i])
      }
    }
    return (res)
  }
  getDatum <- function(block, name, numeric = FALSE) {
    res <- block[grep(name, block[,'name'], ignore.case = TRUE), 'value']
    if (numeric) {
      res <- as.numeric(res)
    } else {
      names(res) <- NULL
    }
    return (res)
  }
  readMatrix <- function(matrix, nchars, chars, tax.labels, gap, missing) {
    chars <- c(chars, gap, missing)
    res <- matrix(nrow = length(tax.labels), ncol = nchars)
    for (i in 1:length (tax.labels)) {
      # loop through each tax label given
      tax.label.found <- FALSE
      for (j in 1:length (matrix)) {
        # loop through each element of the matrix, test if name is present
        if (grepl (tax.labels[i], matrix[j])) {
          tax.label.found <- TRUE
          # if taxlabel is in matrix, add its elements to the res
          characters <- sub(tax.labels[i], "", matrix[j])
          characters <- gsub("\\s+", "", characters)
          characters <- strsplit(characters, "")[[1]]
          while ("{" %in% characters) {
            start <- grep("\\{", characters)
            end <- grep("\\}", characters)
            # convert "{1,2,3}" into "123"
            multichar <- characters[start[1]:end[1]]
            multichar <- paste(multichar[multichar %in% chars], collapse = "")
            # replace in character string
            if (start[1] == 1) {
              characters <- c(multichar, characters[(end[1] + 1):length(characters)])
            } else if (end[1] == length(characters)) {
              characters <- c(characters[1:(start[1] - 1)], multichar)
            } else {
              characters <- c(characters[1:(start[1] - 1)],
                              multichar, characters[(end[1] + 1):length(characters)])
            }
          }
          # test if right number of chars
          if (length(characters) != nchars) {
            error.message <- "Number of characters in matrix does not match number
          specified in metadata."
            stop(paste0("Error! [", tax.labels[i],"] ", error.message))
          }
          # make sure all characters are valid
          # (need to first separate out the potential multichars)
          test.characters <- unlist(strsplit(characters, ""))
          if (!all(test.characters %in% chars)) {
            invalid.character <- test.characters[!test.characters %in% chars][1]
            error.message <- "Character not listed in metadata found!"
            stop(paste0("Error! Label: [", tax.labels[i],"], Character: [", invalid.character,"] "
                        , error.message))
          }
          # convert missing and gap symbols
          if (any(characters %in% missing)) {
            characters[characters %in% missing] <- rep(NA, sum(missing %in% characters))
          }
          if (any(characters %in% gap)) {
            characters[characters %in% gap] <- rep("-1", sum(gap %in% characters))
          }
          # add characters to res
          res[i, ] <- characters
        }
      }
      if (!tax.label.found) {
        # if matrix loop finishes without breaking, tax label is not in matrix
        stop (paste0("Error! Tax label not found in matrix: [", tax.labels[i], "]"))
      }
    }
    rownames(res) <- tax.labels
    return (res)
  }
  # read in data
  data <- scan(file = file, what = character(), sep = "\n", quiet = TRUE, 
               strip.white = TRUE, quote = '\'')
  # extract blocks from nexus
  taxa.block <- getBlock(block.header = 'BEGIN TAXA', data)
  char.block <- getBlock(block.header = 'BEGIN CHARACTERS', data)
  # extract data from blocks
  taxa.block.dim <- getBlockData(label = "DIMENSIONS", taxa.block)
  char.block.dim <- getBlockData(label = "DIMENSIONS", char.block)
  char.block.format <- getBlockData(label = "FORMAT", char.block)
  tax.labels <- getBlockData(label = 'TAXLABELS', taxa.block)
  matrix <- getBlockData(label = 'MATRIX', char.block)
  # extract datum from metadata blocks
  nchars <- getDatum(block = char.block.dim, name = 'nchar', numeric = TRUE)
  m.symbol <- getDatum(block = char.block.format, name = 'missing')
  g.symbol <- getDatum(block = char.block.format, name = 'gap')
  chars <- strsplit(getDatum(block = char.block.format, name = 'symbols'), "")[[1]]
  # convert matrix
  matrix <- readMatrix(matrix = matrix, nchars = nchars, chars = chars, tax.labels = tax.labels,
                       gap = g.symbol, missing = m.symbol)
  return (matrix)
}

#matrix  <- readNexusData ("0_data/test_nexus.nex")
#matrix <- readNexusData ("0_data/mbank_X1766_1-8-2014_922_no_notes.nex")