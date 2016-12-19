# prerequisite packages
requiredPackages.cran <- c(
  "devtools", 
  "reshape",  
  "RColorBrewer", 
  "ggplot2", 
  "gplots", 
  "cowplot")

requiredPackages.bioconductor <- c(
  "genefu",
  "ComplexHeatmap")

# install packages if not installed on your system
# CRAN
install.packages(requiredPackages.cran)
# bioconductor 
source("https://bioconductor.org/biocLite.R")
biocLite(requiredPackages.bioconductor)
