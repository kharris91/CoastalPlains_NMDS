#Load necessary Packages
library("vegan")
library("plyr")

############### Bring in Species Data ###############
## CPL_5 is species cover data for CPL estuarine sites that only includes 
## species that occured on 5% (16) or more of the plots.
cpl5 <- read.csv("cpl_v1_soilplots_soilonly.csv")

# First three columns are UID and wetland class, must be removed before running NMDS 
cpl5 <- cpl5[,-c(2:3)]

############## Calculate Species Measures ##############################
######Species Richness ############
### Sum of non-zero entities per row
s.richness <- apply(cpl5[,-1] > 0, 1, sum)
s.richness
#or
s.richness2 <- ddply(cpl5,~uidplot,function(x) {
  data.frame(RICHNESS = sum(x[-1]>0))
} )
s.richness2 <- s.richness2[order(s.richness2[,2],decreasing = TRUE),]
barplot(s.richness2[,2],names.arg=s.richness2[,1])
###Menhinicks Index (number of species divided by sqrt of numbre of individuals)
menhinick <- function (x) {
  sum(x>0)/sqrt(sum(x))
}
cpl.menhinick <- ddply(cpl5,~uidplot,function(x) {
  data.frame(RICHNESS=menhinick(x[-1]))
  })
cpl.menhinick
############################# NMDS #################################
# Run NMDS using metaMDS function in the vegan package.  
# Set trymax to 1000 in order to reach convergence, dimensions can be added to k
cpl5.ord <- metaMDS(cpl5, trymax = 1000)

# Summary of Ordination
cpl5.ord
summary(cpl5.ord)


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

## Load environmental file
ord <- cpl5.ord
env <- read.csv("cpl_v1_soilplots_soilenvonly.csv")
env <- env[,-c(1,2)]
env <- scale(env)
write.csv(env, "soilplot_env_scaled2.csv")
env <- read.csv("soilplot_env_scaled2.csv")
attach(env)

###
env <- env[,-c(1)]
variable <- STRESS_HARD
plot(ord, disp="sites", type="n")
ordihull(ord, variable, col=1:4, lwd=3)
ordiellipse(ord, variable, col=1:4, kind = "ehull", lwd=3)
ordiellipse(ord, variable, col=1:4, draw="polygon")
ordispider(ord, variable, col=1:4, label = TRUE)
#points(ord, disp="sites", pch=21, col="red", bg="yellow", cex=1.3)

## CCA
#env <- env[,c(3:8)]
#env <- scale(env)
## 
cca <- cca(cpl5 ~ SAND_10 + SILT_10 + CLAY_max_10 + CLAY_mean_10 + pH_H2O_max_10
           + pH_H2O_min_10 + pH_CaCl2_max_10 + pH_CaCl2_min_10 + EC_max_10
           + EC_mean_10 + TOT_CARBON_gkg_10 + TOTC_mmol_10 + InorgC_gkg_10 + OrgC_gkg_mean_10
           + OrgC_gkg_max_10 + OrgC_mmol_10 + TOTN_gkg_10 + TOTN_mmol_10 + TOTP_mgkg_10
           + TOTP_mmol_10 + OLSEN_P_mgkg_10 + OLSEN_P_mmol_10 + MEHLICH_P_mgkg_10 + MEHLICH_P_mmol_10
           + Pox_mgkg_10 + Pox_mmol_10 + S_gkg_10 + TOTS_mmol_10 + CEC_max_10 + CEC_mean_10
           + BaseSat_10 + Ca_cmolkg_10 + K_cmolkg_10 + Mg_cmolkg_10 + Na_cmolkg_10
           + Alox_mmol_10 + Feox_mmol_10 + BD_mean_10 + BD_max_10 + BD_min_10 + TOT_CN_10 
           + Org_CN_10 + TOTC_TOTP_10 + OrgC_TOTP_10 + OrgC_Pox_10 + OrgC_MehlichP_10
           + OrgC_OlsenP_10 + TOTC_TOTS_10 + OrgC_TOTS_10 + TOTN_TOTP_10 + TOTN_Pox_10 
           + TOTN_MehlichP_10 + TOTN_OlsenP_10 + MehlichP_TOTP_10 + OlsenP_TOTP_10
           + Pox_TOTP_10 + PSRox_10 + DPS_10, na.action = na.omit, data=env)
plot(cca, display = c("species", "bp"))
cca2 <- cca(cpl5 ~ SAND_10 + SILT_10 + pH_CaCl2_max_10  
           + EC_mean_10 + InorgC_gkg_10 + OrgC_gkg_mean_10
           + OrgC_mmol_10 + TOTN_gkg_10 + TOTP_mgkg_10
           + OLSEN_P_mgkg_10 + MEHLICH_P_mgkg_10 
           + S_gkg_10 
           + BaseSat_10 + Ca_cmolkg_10 + K_cmolkg_10 + Mg_cmolkg_10 + Na_cmolkg_10
           + Feox_mmol_10 + BD_mean_10 + TOT_CN_10 
           + Org_CN_10 + TOTC_TOTP_10 + OrgC_TOTP_10 +  OrgC_MehlichP_10
           + TOTC_TOTS_10 +  TOTN_TOTP_10 + TOTN_Pox_10 
           + TOTN_MehlichP_10 + TOTN_OlsenP_10 + MehlichP_TOTP_10 + OlsenP_TOTP_10
           + Pox_TOTP_10 + PSRox_10 + DPS_10, na.action = na.omit, data=env)
plot(cca2, display = c("species", "bp"))
cca
summary(cca)

######## Model fitting with ordistep ###############
mod0 <- rda(cpl5 ~ 1, env)  # Model with intercept only
mod0
plot(mod0)
mod1 <- rda(cpl5 ~ ., env, na.action = na.omit)  # Model with all explanatory variables
mod1
plot(mod1)

## With scope present, the default direction is "both"
ordistep(mod0, scope = formula(mod1), perm.max = 200, na.action = na.omit)

## Example without scope. Default direction is "backward"
ordistep(mod1, perm.max = 200) 

###############
ord <- cca(cpl5~  CA_MG_ratio ,env)
fit <- envfit(ord~  CA_MG_ratio,env,perm=999,display="species")
plot(ord, display = c("species","bp"))
ord

##
# VIFs give a measure of the extend of multicollinearity in the predictors of a regression.
# If the VIF of a predictor is high, it indicates that the predictor is highly correlated 
# with other predictors, it contains little or no unique info, and there is 
# redundancy in the set of predictors.
vif.cca(ord)
anova(ord, by="terms")
drop1(ord, test="perm")

barplot(ord$CA$eig/ord$tot.chi, names.arg = 1:ord$CA$rank, cex.names = 0.5, ylab="Proportion of variance explained", xlab="CCA axis")

anova(ord)


