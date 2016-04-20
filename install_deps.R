# Ensure deps are installed

# START
user.input <- readline('WARNING: this script will automaticaly install packages.
                       Are you sure you want to continue? (y/n) ')
if (user.input != 'y') {
  stop ('Execution halted')
}

# GET INSTALLED PACKAGES
packages <- installed.packages ()[ ,1]

# CHECK AND INSTALL
counter <- 0
if (!'ape' %in% packages) {
  install.packages('ape')
  counter <- counter + 1
}
if (!'plyr' %in% packages) {
  install.packages('plyr')
  counter <- counter + 1
}
if (!'geiger' %in% packages) {
  install.packages('geiger')
  counter <- counter + 1
}
if (!'caper' %in% packages) {
  install.packages('caper')
  counter <- counter + 1
}
if (!'ggplot2' %in% packages) {
  install.packages('ggplot2')
  counter <- counter + 1
}
if (!'devtools' %in% packages) {
  install.packages ('devtools')
  counter <- counter + 1
}
if (!'MoreTreeTools' %in% packages) {
  library (devtools)
  install_github('DomBennett/MoreTreeTools')
  counter <- counter + 1
}

# END
cat ('\nComplete! Installed [', counter, '] packages.', sep = '')