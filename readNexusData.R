readNexusData <- function (file)
{
  "find.ntax" <- function(x) {
    for (i in 1:NROW(x)) {
      if (any(f <- grep("\\bntax", x[i], ignore.case = TRUE))) {
        ntax <- as.numeric(sub("(.+?)(ntax\\s*\\=\\s*)(\\d+)(.+)", 
                               "\\3", x[i], perl = TRUE, ignore.case = TRUE))
        break
      }
    }
    ntax
  }
  "find.nchar" <- function(x) {
    for (i in 1:NROW(x)) {
      if (any(f <- grep("\\bnchar", x[i], ignore.case = TRUE))) {
        nchar <- as.numeric(sub("(.+?)(nchar\\s*\\=\\s*)(\\d+)(.+)", 
                                "\\3", x[i], perl = TRUE, ignore.case = TRUE))
        break
      }
    }
    nchar
  }
  "find.matrix.line" <- function(x) {
    for (i in 1:NROW(x)) {
      if (any(f <- grep("\\bmatrix\\b", x[i], ignore.case = TRUE))) {
        matrix.line <- as.numeric(i)
        break
      }
    }
    matrix.line
  }
  "trim.whitespace" <- function(x) {
    gsub("\\s+", "", x)
  }
  "trim.semicolon" <- function(x) {
    gsub(";", "", x)
  }
  X <- scan(file = file, what = character(), sep = "\n", quiet = TRUE, 
            comment.char = "[", strip.white = TRUE)
  ntax <- find.ntax(X)
  nchar <- find.nchar(X)
  matrix.line <- find.matrix.line(X)
  start.reading <- matrix.line + 1
  Obj <- list()
  length(Obj) <- ntax
  i <- 1
  pos <- 0
  tot.nchar <- 0
  tot.ntax <- 0
  for (j in start.reading:NROW(X)) {
    Xj <- trim.semicolon(X[j])
    if (Xj == "") {
      break
    }
    if (any(jtmp <- grep("\\bend\\b", X[j], perl = TRUE, 
                         ignore.case = TRUE))) {
      break
    }
    ts <- unlist(strsplit(Xj, "(?<=\\S)(\\s+)(?=\\S)", perl = TRUE))
    if (length(ts) > 2) {
      stop("nexus parser does not handle spaces in sequences or taxon names (ts>2)")
    }
    if (length(ts) != 2) {
      stop("nexus parser failed to read the sequences (ts!=2)")
    }
    Seq <- trim.whitespace(ts[2])
    Name <- trim.whitespace(ts[1])
    nAME <- paste(c("\\b", Name, "\\b"), collapse = "")
    if (any(l <- grep(nAME, names(Obj)))) {
      tsp <- strsplit(Seq, NULL)[[1]]
      for (k in 1:length(tsp)) {
        p <- k + pos
        Obj[[l]][p] <- tsp[k]
        chars.done <- k
      }
    }
    else {
      names(Obj)[i] <- Name
      tsp <- strsplit(Seq, NULL)[[1]]
      for (k in 1:length(tsp)) {
        p <- k + pos
        Obj[[i]][p] <- tsp[k]
        chars.done <- k
      }
    }
    tot.ntax <- tot.ntax + 1
    if (tot.ntax == ntax) {
      i <- 1
      tot.ntax <- 0
      tot.nchar <- tot.nchar + chars.done
      if (tot.nchar == nchar * ntax) {
        print("ntot was more than nchar*ntax")
        break
      }
      pos <- tot.nchar
    }
    else {
      i <- i + 1
    }
  }
  if (tot.ntax != 0) {
    cat("ntax:", ntax, "differ from actual number of taxa in file?\n")
    stop("nexus parser did not read names correctly (tot.ntax!=0)")
  }
  for (i in 1:length(Obj)) {
    if (length(Obj[[i]]) != nchar) {
      cat(names(Obj[i]), "has", length(Obj[[i]]), "characters\n")
      stop("nchar differ from sequence length (length(Obj[[i]])!=nchar)")
    }
  }
  Obj <- lapply(Obj, tolower)
  Obj
}

data <- readNexusData("0_data/morpho_matrix_forR.nex")

data[[1]]


file <- "0_data/test_nexus.nex"

data <- scan(file = file, what = character(), sep = "\n", quiet = TRUE, 
             strip.white = TRUE, quote = '\'')

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

matrixReader <- function(matrix, nchars, chars, tax.labels) {
  temp.matrix <- matrix
  for (i in 1:length (tax.labels)) {
    for (j in 1:length (matrix)) {
      if (grepl (tax.labels[i], matrix[j])) {
        characters <- sub(tax.labels[i], "", matrix[j])
        characters <- gsub("\\s+", "", characters)
        break
      }
    }
    characters.frame <- rep(NA, nchars)
    characters <- "{013}{0,1,3}111111"
    characters <- strsplit(characters, "")[[1]]
    if ("{" %in% characters) {
      start <- grep("\\{", characters)
      end <- grep("\\}", characters)
      for (j in 1:length(start)) {
        multichar <- characters[start[j]:end[j]]
        multichar <- paste(multichar[multichar %in% chars], collapse = "")
        # how to remove the {0,1,3} once i've reduced the multistate?
        characters <- c(characters[1:start[j]], multichar, characters[end[j]:length(characters)])
      }
    }
    characters.bool <- characters %in% chars
    grep("\\{", characters)
    grep("\\}", characters)
    match(TRUE, characters.bool)
    for (j in 1:length(nchars)) {
      # check if character is allowed, if not assume it is a multiple state character
      strsplit(characters, "")[[1]]
      if (characters[j] %in% chars) {
        characters.frame[j] <- chars[j]
      } else {
        
      }
    }
  }
}

taxa.block <- getBlock(block.header = 'BEGIN TAXA', data)
char.block <- getBlock(block.header = 'BEGIN CHARACTERS', data)
taxa.block.dim <- getBlockData(label = "DIMENSIONS", taxa.block)
char.block.dim <- getBlockData(label = "DIMENSIONS", char.block)
char.block.format <- getBlockData(label = "FORMAT", char.block)
tax.labels <- getBlockData(label = 'TAXLABELS', taxa.block)
matrix <- getBlockData(label = 'MATRIX', char.block)