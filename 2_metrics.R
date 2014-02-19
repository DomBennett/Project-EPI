## Testing Metrics

## Metric 1 (rescale to 0-1)
# removing big clades
#lfi.data <- lfi.data[lfi.data$n < 40, ]
time <- lfi.data$time / max (lfi.data$time)
change <- lfi.data$s.edge.change + lfi.data$d.edge.change
change <- change / max (change)
performance <- lfi.data$d.edge.length
performance <- performance / max (performance)
lfi <- lfiChecker (time, change, performance)

## Exploring LFI
living.fossils <- lfi.data[lfi > cutoff, ]
living.fossils <- living.fossils[order (living.fossils$lfi, decreasing = TRUE), ]
plot (phylo)
nodelabels (text = round (living.fossils$lfi, digits = 3), node = living.fossils$node)