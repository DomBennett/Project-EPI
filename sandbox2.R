data <- read.csv("2_mmetrics/LFI_mammal_masses.csv", stringsAsFactors = FALSE)
data <- na.omit(data)
data$N <- as.vector(sapply(data$taxon, function(x) length(gregexpr("\\|", x)[[1]])))

#data <- data[data$mass.change != 0, ] #exclude species with no mass change
hist(log(data$mass.by.edge))
data <- data[data$N > 1, ] # remove single species
head(data)
data$log.mass.by.edge <- log(data$mass.by.edge)
hist(data$log.mass.by.edge) # non-normal distribution

head(data)