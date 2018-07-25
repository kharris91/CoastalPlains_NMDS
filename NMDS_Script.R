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
## Plot cpl5.ord
plot(cpl5.ord, type = "n")
points(cpl5.ord, display = "sites", cex = 0.8, pch=21, col="red", bg="yellow")
text(cpl5.ord, display = "spec", cex=0.7, col="blue")
## View as stress plot
stressplot(cpl5.ord)
## By adding a third dimension to the ordination, we will be able to create a 3d plot
## This also lowers the stress values
cpl5.ord.3d <- metaMDS(cpl5, trymax = 1000, k = 3)
# View summary
cpl5.ord.3d
# Plot 2d again
plot(cpl5.ord.3d, type = "n")
points(cpl5.ord.3d, display = "sites", cex = 0.8, pch=21, col="red", bg="yellow")
text(cpl5.ord.3d, display = "spec", cex=0.7, col="blue")
# Plot in 3d
library("vegan3d")
library("rgl")
ordirgl(cpl5.ord.3d, scaling = 3, display = "sites")
orglpoints(cpl5.ord.3d, display = "sites", choices = 1:3, col = "red")
orgltext(cpl5.ord.3d, display = "species", choices = 1:3, adj = 0.5,
         col = "black")