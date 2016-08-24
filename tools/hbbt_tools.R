library(ggplot2)

ggBinomial <- function(data) {
  #http://docs.ggplot2.org/current/geom_smooth.html
  binomial_smooth <- function(...) {
    geom_smooth(method = "glm", method.args = list(family = "binomial"), ...)
  }
  # To fit a logistic regression, you need to coerce the values to
  # a numeric vector lying between 0 and 1.
  # add quotes to htypes
  data$type <- paste0("'", data$type, "'")
  p <- ggplot(data, aes(x=pepi, y=scr, colour=type, group=type)) +
    binomial_smooth(aes(fill=type)) +
    ylab("Prop. in habitat") + xlab("pEPI") +
    scale_fill_discrete(name="Habitat") +
    scale_colour_discrete(name="Habitat") +
    theme_bw()
  p
}