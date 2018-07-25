#Load necessary Packages
library("vegan")
## CPL_5 is species cover data for CPL estuarine sites that only includes 
## species that occured on 5% (16) or more of the plots.
cpl5 <- read.csv("CPL_5.csv")
# First two columns are UID and wetland class, must be removed before running NMDS 
cpl5 <- cpl5[,-c(1:2)]
# Run NMDS using metaMDS function in the vegan package.  
# Set trymax to 1000 in order to reach convergence, dimensions can be added to k
cpl5.ord <- metaMDS(cpl5, trymax = 1000)
# Summary of Ordination
cpl5.ord