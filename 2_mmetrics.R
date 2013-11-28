## 18/11/2013
## D.J. Bennett
## Using the body mass data to determine morphological stasis
## TODO: Expand to include clades, not just species.

## Libraries
source(file.path("Functions", "EcoDataTools.R"))
require("reshape2")

## Dirs
input.dir <- "0_data"
output.dir <- "2_mmetrics"

## Data
data <- read.delim(file.path(input.dir, "qbodydata.txt"), na.strings = -999)
data <- na.omit(data)
#head(data)

## Data manip
data$Binomial <- paste(data$Genus, data$Species)
head(data)
molten <- melt(data, id.vars = c("Continent", "Status", "Order", "Family", "Genus",
                                 "Species", "References", "Binomial"))
data <- dcast(molten, formula = Status + Order + Family + Genus + Binomial ~ variable, mean)
#head(data)


## Calculating character change
# rough and ready and this stage, calc mean for extinct genus, family or order
# subtract this number from extant species in the same clade
extinct <- data[data$Status == "extinct", ]
extant <- data[data$Status == "extant", ]
mass.orders <- tapply(extinct$Log_mass, extinct$Order, mean)
mass.orders <- mass.orders[!is.na(mass.orders)]
mass.family <- tapply(extinct$Log_mass, extinct$Family, mean)
mass.family <- mass.family[!is.na(mass.family)]
mass.genus <- tapply(extinct$Log_mass, extinct$Genus, mean)
mass.genus <- mass.genus[!is.na(mass.genus)]
# calc difference
ancestral.mass <- rep(NA, nrow(extant))
for (i in 1:nrow(extant)) {
  #if (extant$Genus[i] %in% names(mass.genus)) {
  #  ancestral.mass[i] <- mass.genus[names(mass.genus) %in% extant$Genus[i]]
  #} else if (extant$Family[i] %in% names(mass.family)) {
  #  ancestral.mass[i] <- mass.family[names(mass.family) %in% extant$Family[i]]
  #} else 
  if (extant$Order[i] %in% names(mass.orders)) { # trying just orders
    ancestral.mass[i] <- mass.orders[names(mass.orders) %in% extant$Order[i]]
  } else {
    ancestral.mass[i] <- NA
  }
}
extant$Ancestral.mass <- ancestral.mass
extant <- extant[!is.na(ancestral.mass),]
extant$Mass.change <- abs(extant$Log_mass - extant$Ancestral.mass)
extant <- extant[order(extant$Mass.change, decreasing = FALSE),]

## Output
write.csv(x = extant, file = file.path(output.dir, "mmetrics.csv"), row.names = FALSE)

## Explore
hist(extant$Mass.change)
names(extant)
extant$Binomial[order(abs(extant$Mass.change), decreasing = TRUE)]
hist(abs(extant$Mass.change))