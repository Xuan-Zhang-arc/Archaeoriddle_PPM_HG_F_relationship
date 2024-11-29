
# Title: Using Point Process Modelling to detect the peaceful vs hostile relationship
# Author: Xuan Zhang, PhD student in School of Architecture, Tianjin University, China; visiting PhD student in Department of Archaeology, University of Cambridge, UK (2022-2023)

# Note:
# The file paths have been standardized using the here package, so users do not need to adjust the paths.
# A file named all_sites_for_R.csv is provided, which includes all data obtained after selecting five additional grids and processing the radiocarbon dates via Oxcal Online. Users may choose to skip these steps.
# The maptools pacakge has retired and removed from the CRAN repository, so please visit 'https://cran.r-project.org/src/contrib/Archive/maptools/' to obtain former available versions and install the package. The users can also use other packages as suggested on 'https://cran.r-project.org/web/packages/maptools/index.html' as a substitution. However, switching to other packages might have an impact on the results.

# prepare the packages
library(readr)
library(spatstat)
library(spatstat.model)
library(raster)
library(maptools)
library(sf)
library(here)

# package versions
packageVersion("readr")
#[1] '2.1.5'
packageVersion("spatstat")
#[1] '3.0.2' 
packageVersion("spatstat.model")
#[1] '3.2.1.7'
packageVersion("raster")
#[1] '3.6.3'
packageVersion("maptools")
#[1] '1.1.5'
packageVersion("sf")
#[1] '1.0.9'
packageVersion("here")
#[1] '1.0.1'

#1.select aother 5 squares
#1.1 prepare the spatstat objects
#read the csv of sites
FHG <- read_csv(here("Sites","Biblio_data.csv"))
F <- subset(FHG,FHG$economy=="F")
HG <- subset(FHG,FHG$economy=="HG")
# read the maps
map <- raster(here("Layers","east_narnia4x.tif"))
resource <- raster(here("Layers","resources.tiff"))
target_crs <- '+proj=longlat +datum=WGS84 +no_defs'
# create the window
ext <- extent(map)
polygon <- as(ext, "SpatialPolygons")
owin=as.owin(polygon)
# transfer the sites and maps to spatstat formats
F_ppp<-ppp(x=F$lon,y=F$lat,crs=target_crs,window=owin)
HG_ppp<-ppp(x=HG$lon,y=HG$lat,crs=target_crs,window=owin)
FHG_ppp<-ppp(x=FHG$lon,y=FHG$lat,crs=target_crs,window=owin)
map_im <- as.im(map)
resource_im <- as.im(resource)
#1.2 PPM
#1.2.1 Farmers
ppm_F_map<-ppm(F_ppp~map_im)
plot(predict(ppm_F_map))
AIC(ppm_F_map)
# [1] 50.99587
ppm_F_map_Kres_env <-envelope(ppm_F_map,fun=Kres,nsim=99,correction='best')
plot(ppm_F_map_Kres_env)

ppm_F_resource<-ppm(F_ppp~resource_im)
# Warning message:
# Values of the covariate 'resource_im' were NA or undefined at 29% (306 out of 10
# 53) of the quadrature points. Occurred while executing: ppm.ppp(Q = F_ppp, trend
#  = ~resource_im, data = NULL, interaction = NULL)
# Seen from the warning message, the 29% missing value is too high so the resouce_im was not used for the grid selection.

gam_F_map<-ppm(F_ppp~s(map_im),use.gam=TRUE)
plot(predict(gam_F_map))
AIC(gam_F_map)
# [1] 20.94797
gam_F_map_Kres_env <-envelope(gam_F_map,fun=Kres,nsim=99,correction='best')
plot(gam_F_map_Kres_env)

# the predictive map for grid selection
jpeg(here("Output", "1_select_5_more", "predict(gam_F_map).jpg"),width = 7, height = 7, units = "in", res = 1200)
plot(predict(gam_F_map))
dev.off()

#1.2.2 Hunter-gatherers
ppm_HG_map<-ppm(HG_ppp~map_im)
plot(predict(ppm_HG_map))
AIC(ppm_HG_map)
# [1] 46.55376
ppm_HG_map_Kres_env <-envelope(ppm_HG_map,fun=Kres,nsim=99,correction='best')
plot(ppm_HG_map_Kres_env)

gam_HG_map<-ppm(HG_ppp~s(map_im),use.gam=TRUE)
plot(predict(gam_HG_map))
AIC(gam_HG_map)
# [1] 27.1717
gam_HG_map_Kres_env <-envelope(gam_HG_map,fun=Kres,nsim=99,correction='best')
plot(gam_HG_map_Kres_env)

# the predictive map for grid selection
jpeg(here("Output", "1_select_5_more", "predict(gam_HG_map).jpg"),width = 7, height = 7, units = "in", res = 1200)
plot(predict(gam_HG_map))
dev.off()

#1.2.3 all sites
ppm_FHG_map<-ppm(FHG_ppp~map_im)
plot(predict(ppm_FHG_map))
AIC(ppm_FHG_map)
# [1]42.62654 
ppm_FHG_map_Kres_env <-envelope(ppm_FHG_map,fun=Kres,nsim=99,correction='best')
plot(ppm_FHG_map_Kres_env)

gam_FHG_map<-ppm(FHG_ppp~s(map_im),use.gam=TRUE)
plot(predict(gam_FHG_map))
AIC(gam_FHG_map)
# [1]-18.54415 
gam_FHG_map_Kres_env <-envelope(gam_FHG_map,fun=Kres,nsim=99,correction='best')
plot(gam_FHG_map_Kres_env)

# the predictive map for grid selection
jpeg(here("Output", "1_select_5_more", "predict(gam_FHG_map).jpg"),width = 7, height = 7, units = "in", res = 1200)
plot(predict(gam_FHG_map))
dev.off()


#2.Analysis
#2.1 prepare the spatstat objects
#read the csv of sites
FHG <- read_csv(here("Sites", "all_sites_for_R.csv"))
early.FHG <- subset(FHG, FHG$early_period==TRUE)
early.F <- subset(early.FHG,early.FHG$economy=="F")
early.HG <- subset(early.FHG,early.FHG$economy=="HG")

middle.FHG <- subset(FHG, FHG$middle_period==TRUE)
middle.F <- subset(middle.FHG,middle.FHG$economy=="F")
duplicated_rows <- duplicated(cbind(middle.F$longitude, middle.F$latitude))
middle.F <- middle.F[!duplicated_rows, ]
middle.HG <- subset(middle.FHG,middle.FHG$economy=="HG")

late.FHG <- subset(FHG, FHG$late_period==TRUE)
late.F <- subset(late.FHG,late.FHG$economy=="F")
late.HG <- subset(late.FHG,late.FHG$economy=="HG")

middle.early.FHG <- subset(FHG,FHG$middle_early_period==TRUE)
middle.early.F <- subset(middle.early.FHG,middle.early.FHG$economy=="F")
duplicated_rows <- duplicated(cbind(middle.early.F$longitude, middle.early.F$latitude))
middle.early.F <- middle.early.F[!duplicated_rows, ]
middle.early.HG<- subset(middle.early.FHG,middle.early.FHG$economy=="HG")

middle.late.FHG <- subset(FHG,FHG$middle_late_period==TRUE)
middle.late.F <- subset(middle.late.FHG,middle.late.FHG$economy=="F")
middle.late.HG<- subset(middle.late.FHG,middle.late.FHG$economy=="HG")
# read the maps
map <- raster(here("Layers","east_narnia4x.tif"))
resource <- raster(here("Layers","resources.tiff"))
target_crs <- crs(map)
# create the window
ext <- extent(map)
polygon <- as(ext, "SpatialPolygons")
owin=as.owin(polygon,crs=target_crs)
# transfer the sites and maps to spatstat formats
map.im <- as.im(map)
resource.im <- as.im(resource)
plot(resource.im)
early.F.ppp<-ppp(x=early.F$longitude,y=early.F$latitude,crs=target_crs,window=owin)
early.HG.ppp<-ppp(x=early.HG$longitude,y=early.HG$latitude,crs=target_crs,window=owin)
early.FHG.multitype<-superimpose(Farmer=early.F.ppp,Hunter_gatherer=early.HG.ppp)
plot(early.FHG.multitype)
middle.F.ppp<-ppp(x=middle.F$longitude,y=middle.F$latitude,crs=target_crs,window=owin)
middle.HG.ppp<-ppp(x=middle.HG$longitude,y=middle.HG$latitude,crs=target_crs,window=owin)
middle.FHG.multitype<-superimpose(Farmer=middle.F.ppp,Hunter_gatherer=middle.HG.ppp)
plot(middle.FHG.multitype)

middle.early.F.ppp<-ppp(x=middle.early.F$longitude,y=middle.early.F$latitude,crs=target_crs,window=owin)
middle.early.HG.ppp<-ppp(x=middle.early.HG$longitude,y=middle.early.HG$latitude,crs=target_crs,window=owin)
middle.early.FHG.multitype<-superimpose(Farmer=middle.early.F.ppp,Hunter_gatherer=middle.early.HG.ppp)
plot(middle.early.FHG.multitype)

middle.late.F.ppp<-ppp(x=middle.late.F$longitude,y=middle.late.F$latitude,crs=target_crs,window=owin)
middle.late.HG.ppp<-ppp(x=middle.late.HG$longitude,y=middle.late.HG$latitude,crs=target_crs,window=owin)
middle.late.FHG.multitype<-superimpose(Farmer=middle.late.F.ppp,Hunter_gatherer=middle.late.HG.ppp)
plot(middle.late.FHG.multitype)

late.F.ppp<-ppp(x=late.F$longitude,y=late.F$latitude,crs=target_crs,window=owin)
late.HG.ppp<-ppp(x=late.HG$longitude,y=late.HG$latitude,crs=target_crs,window=owin)
late.FHG.multitype<-superimpose(Farmer=late.F.ppp,Hunter_gatherer=late.HG.ppp)
plot(late.FHG.multitype)

#2.2 early stage
#Cross-type K-function
Kcross.early.FHG<-envelope(early.FHG.multitype,fun=Kcross,nsim=99)
jpeg(here("Output", "2_analysis", "Kcross_early_FHG.jpg"),width = 7, height = 7, units = "in", res = 1200)
plot(Kcross.early.FHG)
dev.off()

# PPM
early.hf <- hyperframe(p=early.FHG.multitype)
early.hf$m <- with(early.hf,map.im)
early.hf$r <- with(early.hf,resource.im)
early.mppm1 <- mppm(p~m+r,data=early.hf)
early.dist.1 <- distfun(early.FHG.multitype)
summary(early.dist.1)
# Distance function for point pattern
# defined in a polygonal window inside the frame [-4, 1] x [-1, 4] units
# 
# Distance function values:
#         range = [0.001513771, 2.628423]
#         mean = 0.8335595
early.dist.2 <- pairdist(early.FHG.multitype)
summary(early.dist.2)
max(early.dist.2)
# [1] 4.462367
mean(early.dist.2)
# [1] 1.814042

#multitype strauss
early.RR <- data.frame(R=seq(0.01, 4.47, by=0.01))
early.MS <- function(R) { MultiStrauss(radii=diag(c(R,R))) }
early.pm <- profilepl(early.RR, early.MS,early.FHG.multitype~marks*polynom(x,y,3),rbord=NULL,correction="Ripley")
early.pm
# profile log pseudolikelihood for model:
# ppm(early.FHG.multitype ~ marks * polynom(x,  y,  3),  correction = "Ripley",   
#     interaction = early.MS,  rbord = NULL)
# interaction: Multitype Strauss process
# irregular parameter: R in [0.01, 4.47]
# optimum value of irregular parameter:  R = 1.38
jpeg(here("Output", "2_analysis", "early_pm.jpg"),width = 7, height = 7, units = "in", res = 1200)
plot(early.pm)
dev.off()

as.ppm(early.pm)
#  Nonstationary Multitype Strauss process
# Fitted to point pattern dataset 'early.FHG.multitype'
# 
# Possible marks: 'Farmer' and 'Hunter_gatherer'
# 
# Log trend:  ~marks * (x + y + I(x^2) + I(x * y) + I(y^2) + I(x^3) + I(x^2 * y)  
# + I(x * y^2) + I(y^3))
# 
# Fitted trend coefficients:
#                     (Intercept)            marksHunter_gatherer
#                      10.4412578                      -6.3761620
#                               x                               y
#                       6.1434358                       8.5209893 
#                          I(x^2)                        I(x * y)
#                       4.1728062                      -3.3570407
#                          I(y^2)                          I(x^3)
#                      -3.2781182                       1.4156157 
#                      I(x^2 * y)                      I(x * y^2)
#                      -3.0128819                      -2.8069921
#                          I(y^3)          marksHunter_gatherer:x
#                      -0.9011120                       6.7694587 
#          marksHunter_gatherer:y     marksHunter_gatherer:I(x^2)
#                      10.3112259                      -8.9345371
#   marksHunter_gatherer:I(x * y)     marksHunter_gatherer:I(y^2)
#                     -11.1450437                      -8.4206572 
#     marksHunter_gatherer:I(x^3) marksHunter_gatherer:I(x^2 * y)
#                      -0.9691066                       3.3180110
# marksHunter_gatherer:I(x * y^2)     marksHunter_gatherer:I(y^3) 
#                       4.8446464                       2.5790241
# 
# 2 types of points
# Possible types:
# [1] Farmer          Hunter_gatherer
# Interaction radii:
#                 Farmer Hunter_gatherer
# Farmer            1.38              NA
# Hunter_gatherer     NA            1.38
# Fitted interaction parameters gamma_ij
#                    Farmer Hunter_gatherer
# Farmer          0.1136299              NA
# Hunter_gatherer        NA       0.5183687
# 
# Relevant coefficients:
#                   markFarmerxFarmer markHunter_gathererxHunter_gatherer
#                          -2.1748084                          -0.6570685 
# 
# For standard errors, type coef(summary(x))
# 
anova(as.ppm(early.pm),test="LR")
# Analysis of Deviance Table
# 
# Terms added sequentially (first to last)
# 
# Model 1: ~marks + x + y + I(x^2) + I(x * y) + I(y^2) + I(x^3) + I(x^2 * y) + I(x
#  * y^2) + I(y^3)         MultiStrauss
# Model 2: ~marks + x + y + I(x^2) + I(x * y) + I(y^2) + I(x^3) + I(x^2 * y) + I(x
#  * y^2) + I(y^3) + marks:x       MultiStrauss
# Model 3: ~marks + x + y + I(x^2) + I(x * y) + I(y^2) + I(x^3) + I(x^2 * y) + I(x
#  * y^2) + I(y^3) + marks:x + marks:y     MultiStrauss
# Model 4: ~marks + x + y + I(x^2) + I(x * y) + I(y^2) + I(x^3) + I(x^2 * y) + I(x
#  * y^2) + I(y^3) + marks:x + marks:y + marks:I(x^2)      MultiStrauss
# Model 5: ~marks + x + y + I(x^2) + I(x * y) + I(y^2) + I(x^3) + I(x^2 * y) + I(x
#  * y^2) + I(y^3) + marks:x + marks:y + marks:I(x^2) + marks:I(x * y)     MultiSt
# rauss
# Model 6: ~marks + x + y + I(x^2) + I(x * y) + I(y^2) + I(x^3) + I(x^2 * y) + I(x
#  * y^2) + I(y^3) + marks:x + marks:y + marks:I(x^2) + marks:I(x * y) + marks:I(y
# ^2)      MultiStrauss
# Model 7: ~marks + x + y + I(x^2) + I(x * y) + I(y^2) + I(x^3) + I(x^2 * y) + I(x
#  * y^2) + I(y^3) + marks:x + marks:y + marks:I(x^2) + marks:I(x * y) + marks:I(y
# ^2) + marks:I(x^3)       MultiStrauss
# Model 8: ~marks + x + y + I(x^2) + I(x * y) + I(y^2) + I(x^3) + I(x^2 * y) + I(x
#  * y^2) + I(y^3) + marks:x + marks:y + marks:I(x^2) + marks:I(x * y) + marks:I(y
# ^2) + marks:I(x^3) + marks:I(x^2 * y)    MultiStrauss
# Model 9: ~marks + x + y + I(x^2) + I(x * y) + I(y^2) + I(x^3) + I(x^2 * y) + I(x
#  * y^2) + I(y^3) + marks:x + marks:y + marks:I(x^2) + marks:I(x * y) + marks:I(y
# ^2) + marks:I(x^3) + marks:I(x^2 * y) + marks:I(x * y^2)         MultiStrauss   
# Model 10: ~marks * (x + y + I(x^2) + I(x * y) + I(y^2) + I(x^3) + I(x^2 * y) + I
# (x * y^2) + I(y^3))      MultiStrauss
#    Npar Df AdjDeviance  Pr(>Chi)
# 1    13
# 2    14  1        8.17  0.004263 ** 
# 3    15  1      362.09 < 2.2e-16 ***
# 4    16  1       41.07 1.466e-10 ***
# 5    17  1       23.50 1.248e-06 ***
# 6    18  1        4.46  0.034696 *
# 7    19  1       16.68 4.414e-05 ***
# 8    20  1        1.26  0.261180
# 9    21  1        1.16  0.282233
# 10   22  1      131.67 < 2.2e-16 ***
# ---
# Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

#2.3 middle stage
#2.3.1 middle.early stage
#Cross-type K-function
Kcross.middle.early.FHG<-envelope(middle.early.FHG.multitype,fun=Kcross,nsim=99)
jpeg(here("Output", "2_analysis", "Kcross_middle_early_FHG.jpg"),width = 7, height = 7, units = "in", res = 1200)
plot(Kcross.middle.early.FHG)
dev.off()

# PPM
middle.early.hf <- hyperframe(p=middle.early.FHG.multitype)
middle.early.hf$m <- with(middle.early.hf,map.im)
middle.early.hf$r <- with(middle.early.hf,resource.im)
middle.early.mppm1 <- mppm(p~m+r,data=middle.early.hf)
plot(predict(middle.early.mppm1))
middle.early.dist1 <- distfun(middle.early.FHG.multitype)
summary(middle.early.dist1)
# Distance function for point pattern
# defined in a polygonal window inside the frame [-4, 1] x [-1, 4] units
# 
# Distance function values:
#         range = [0.001513771, 2.628423]
#         mean = 0.7620134
middle.early.dist2 <- pairdist(middle.early.FHG.multitype)
summary(middle.early.dist2)
max(middle.early.dist2)
# [1] 4.559198
mean(middle.early.dist2)
# [1] 1.840964
ks.test(early.dist.2,middle.early.dist2)
# Asymptotic two-sample Kolmogorov-Smirnov test
# 
# data:  early.dist.2 and middle.early.dist2
# D = 0.021377, p-value = 0.4391
# alternative hypothesis: two-sided
# 
# Warning message:
# In ks.test.default(early.dist.2, middle.early.dist2) :
#   p-value will be approximate in the presence of ties

#multitype-strauss
middle.early.RR <- data.frame(R=seq(0.01, 4.56, by=0.01))
middle.early.MS <- function(R) { MultiStrauss(radii=diag(c(R,R))) }
middle.early.pm <- profilepl(middle.early.RR, middle.early.MS,middle.early.FHG.multitype~marks*polynom(x,y,3), rbord = NULL, correction="Ripley")
middle.early.pm
# profile log pseudolikelihood for model:
# ppm(middle.early.FHG.multitype ~ marks * polynom(x,  y,  3),  correction =      
# "Ripley",       interaction = middle.early.MS,  rbord = NULL)
# interaction: Multitype Strauss process
# irregular parameter: R in [0.01, 4.56]
# optimum value of irregular parameter:  R = 1.24
jpeg(here("Output", "2_analysis", "middle_early_pm.jpg"),width = 7, height = 7, units = "in", res = 1200)
plot(middle.early.pm)
dev.off()

as.ppm(middle.early.pm)
# Nonstationary Multitype Strauss process
# Fitted to point pattern dataset 'middle.early.FHG.multitype'
# 
# Possible marks: 'Farmer' and 'Hunter_gatherer'
# 
# Log trend:  ~marks * (x + y + I(x^2) + I(x * y) + I(y^2) + I(x^3) + I(x^2 * y)  
# + I(x * y^2) + I(y^3))
# 
# Fitted trend coefficients:
#                     (Intercept)            marksHunter_gatherer 
#                       3.9899822                      -1.2176194
#                               x                               y
#                       3.2188311                       4.7110527
#                          I(x^2)                        I(x * y) 
#                       2.8992131                      -1.7036887
#                          I(y^2)                          I(x^3)
#                      -1.3950067                       0.8689269
#                      I(x^2 * y)                      I(x * y^2) 
#                      -1.6996094                      -1.0422251
#                          I(y^3)          marksHunter_gatherer:x
#                      -0.3060087                       9.7969151
#          marksHunter_gatherer:y     marksHunter_gatherer:I(x^2) 
#                      14.6951409                      -7.3724131
#   marksHunter_gatherer:I(x * y)     marksHunter_gatherer:I(y^2)
#                     -12.2598512                     -10.5548084
#     marksHunter_gatherer:I(x^3) marksHunter_gatherer:I(x^2 * y) 
#                      -1.4221105                       0.8388868
# marksHunter_gatherer:I(x * y^2)     marksHunter_gatherer:I(y^3)
#                       2.7411405                       2.0721103
# 
# 2 types of points
# Possible types:
# [1] Farmer          Hunter_gatherer
# Interaction radii:
#                 Farmer Hunter_gatherer
# Farmer            1.24              NA
# Hunter_gatherer     NA            1.24
# Fitted interaction parameters gamma_ij
#                    Farmer Hunter_gatherer
# Farmer          0.6852975              NA
# Hunter_gatherer        NA       0.5347298
# 
# Relevant coefficients:
#                   markFarmerxFarmer markHunter_gathererxHunter_gatherer
#                          -0.3779022                          -0.6259936 
# 
# For standard errors, type coef(summary(x))
anova((as.ppm(middle.early.pm)),test="LR")
# Analysis of Deviance Table
# 
# Terms added sequentially (first to last)
# 
# Model 1: ~marks + x + y + I(x^2) + I(x * y) + I(y^2) + I(x^3) + I(x^2 * y) + I(x
#  * y^2) + I(y^3)         MultiStrauss
# Model 2: ~marks + x + y + I(x^2) + I(x * y) + I(y^2) + I(x^3) + I(x^2 * y) + I(x
#  * y^2) + I(y^3) + marks:x       MultiStrauss
# Model 3: ~marks + x + y + I(x^2) + I(x * y) + I(y^2) + I(x^3) + I(x^2 * y) + I(x
#  * y^2) + I(y^3) + marks:x + marks:y     MultiStrauss
# Model 4: ~marks + x + y + I(x^2) + I(x * y) + I(y^2) + I(x^3) + I(x^2 * y) + I(x
#  * y^2) + I(y^3) + marks:x + marks:y + marks:I(x^2)      MultiStrauss
# Model 5: ~marks + x + y + I(x^2) + I(x * y) + I(y^2) + I(x^3) + I(x^2 * y) + I(x
#  * y^2) + I(y^3) + marks:x + marks:y + marks:I(x^2) + marks:I(x * y)     MultiSt
# rauss
# Model 6: ~marks + x + y + I(x^2) + I(x * y) + I(y^2) + I(x^3) + I(x^2 * y) + I(x
#  * y^2) + I(y^3) + marks:x + marks:y + marks:I(x^2) + marks:I(x * y) + marks:I(y
# ^2)      MultiStrauss
# Model 7: ~marks + x + y + I(x^2) + I(x * y) + I(y^2) + I(x^3) + I(x^2 * y) + I(x
#  * y^2) + I(y^3) + marks:x + marks:y + marks:I(x^2) + marks:I(x * y) + marks:I(y
# ^2) + marks:I(x^3)       MultiStrauss
# Model 8: ~marks + x + y + I(x^2) + I(x * y) + I(y^2) + I(x^3) + I(x^2 * y) + I(x
#  * y^2) + I(y^3) + marks:x + marks:y + marks:I(x^2) + marks:I(x * y) + marks:I(y
# ^2) + marks:I(x^3) + marks:I(x^2 * y)    MultiStrauss
# Model 9: ~marks + x + y + I(x^2) + I(x * y) + I(y^2) + I(x^3) + I(x^2 * y) + I(x
#  * y^2) + I(y^3) + marks:x + marks:y + marks:I(x^2) + marks:I(x * y) + marks:I(y
# ^2) + marks:I(x^3) + marks:I(x^2 * y) + marks:I(x * y^2)         MultiStrauss   
# Model 10: ~marks * (x + y + I(x^2) + I(x * y) + I(y^2) + I(x^3) + I(x^2 * y) + I
# (x * y^2) + I(y^3))      MultiStrauss
#    Npar Df AdjDeviance  Pr(>Chi)
# 1    13
# 2    14  1      12.938 0.0003220 ***
# 3    15  1     245.386 < 2.2e-16 ***
# 4    16  1      17.187 3.388e-05 ***
# 5    17  1      35.628 2.388e-09 ***
# 6    18  1       4.413 0.0356669 *
# 7    19  1      -2.912 0.0879171 .
# 8    20  1       0.661 0.4161415
# 9    21  1      11.528 0.0006854 ***
# 10   22  1    -206.515 < 2.2e-16 ***
# ---
# Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

#2.3.2 middle.late stage
#Cross-type K-function
Kcross.middle.late.FHG<-envelope(middle.late.FHG.multitype,fun=Kcross,nsim=99)
jpeg(here("Output", "2_analysis", "Kcross_middle_late_FHG.jpg"),width = 7, height = 7, units = "in", res = 1200)
plot(Kcross.middle.late.FHG)
dev.off()

# PPM
middle.late.hf <- hyperframe(p=middle.late.FHG.multitype)
middle.late.hf$m <- with(middle.late.hf,map.im)
middle.late.hf$r <- with(middle.late.hf,resource.im)
middle.late.mppm1 <- mppm(p~m+r,data=middle.late.hf)
middle.late.dist1 <- distfun(middle.late.FHG.multitype)
summary(middle.late.dist1)
# Distance function for point pattern
# defined in a polygonal window inside the frame [-4, 1] x [-1, 4] units
# 
# Distance function values:
#         range = [0.01381068, 2.628423]
#         mean = 0.8290538
middle.late.dist2 <- pairdist(middle.late.FHG.multitype)
summary(middle.late.dist2)
max(middle.late.dist2)
# [1] 4.384972
mean(middle.late.dist2)
# [1] 1.817461
ks.test(middle.early.dist2,middle.late.dist2)
# Asymptotic two-sample Kolmogorov-Smirnov test
# 
# data:  middle.early.dist2 and middle.late.dist2
# D = 0.053897, p-value = 0.0001859
# alternative hypothesis: two-sided
# 
# Warning message:
# In ks.test.default(middle.early.dist2, middle.late.dist2) :
#   p-value will be approximate in the presence of ties
ks.test(early.dist.2,middle.late.dist2)
# Asymptotic two-sample Kolmogorov-Smirnov test
# 
# data:  early.dist.2 and middle.late.dist2
# D = 0.049088, p-value = 0.008723
# alternative hypothesis: two-sided
# 
# Warning message:
# In ks.test.default(early.dist.2, middle.late.dist2) :
#   p-value will be approximate in the presence of ties

#multitype strauss
middle.late.RR <- data.frame(R=seq(0.01, 4.39, by=0.01))
middle.late.MS <- function(R) { MultiStrauss(radii=diag(c(R,R))) }
middle.late.pm <- profilepl(middle.late.RR, middle.late.MS,middle.late.FHG.multitype~marks*polynom(x,y,3), rbord = NULL, correction="Ripley")
middle.late.pm
# profile log pseudolikelihood for model:
# ppm(middle.late.FHG.multitype ~ marks * polynom(x,  y,  3),  correction = 
# "Ripley",       interaction = middle.late.MS,  rbord = NULL)
# interaction: Multitype Strauss process
# irregular parameter: R in [0.01, 4.39]
# optimum value of irregular parameter:  R = 1.25
jpeg(here("Output", "2_analysis", "middle_late_pm.jpg"),width = 7, height = 7, units = "in", res = 1200)
plot(middle.late.pm)
dev.off()

as.ppm(middle.late.pm)
# Nonstationary Multitype Strauss process
# Fitted to point pattern dataset 'middle.late.FHG.multitype'
# 
# Possible marks: 'Farmer' and 'Hunter_gatherer'
# 
# Log trend:  ~marks * (x + y + I(x^2) + I(x * y) + I(y^2) + I(x^3) + I(x^2 * y) 
# + I(x * y^2) + I(y^3))
# 
# Fitted trend coefficients:
#                     (Intercept)            marksHunter_gatherer
#                       6.5422430                     -15.0417177 
#                               x                               y
#                       4.9113544                       6.5033586
#                          I(x^2)                        I(x * y)
#                       2.3796293                      -3.4166658 
#                          I(y^2)                          I(x^3)
#                      -2.1317860                       0.6921642
#                      I(x^2 * y)                      I(x * y^2) 
#                      -2.3032571                      -1.1956764
#                          I(y^3)          marksHunter_gatherer:x
#                      -0.3609944                       3.4220299
#          marksHunter_gatherer:y     marksHunter_gatherer:I(x^2) 
#                      32.3235942                      -7.5128501
#   marksHunter_gatherer:I(x * y)     marksHunter_gatherer:I(y^2)
#                      -3.9862308                     -18.7082635 
#     marksHunter_gatherer:I(x^3) marksHunter_gatherer:I(x^2 * y)
#                      -0.3347197                       1.4673152
# marksHunter_gatherer:I(x * y^2)     marksHunter_gatherer:I(y^3)
#                       1.1170386                       3.4090225 
# 
# 2 types of points
# Possible types:
# [1] Farmer          Hunter_gatherer
# Interaction radii:
#                 Farmer Hunter_gatherer
# Farmer            1.25              NA
# Hunter_gatherer     NA            1.25
# Fitted interaction parameters gamma_ij
#                    Farmer Hunter_gatherer
# Farmer          0.4715164              NA
# Hunter_gatherer        NA       0.1740358
# 
# Relevant coefficients:
#                   markFarmerxFarmer markHunter_gathererxHunter_gatherer
#                          -0.7518014                          -1.7484941
# 
# For standard errors, type coef(summary(x))
anova((as.ppm(middle.late.pm)),test="LR")
# Analysis of Deviance Table
# 
# Terms added sequentially (first to last)
# 
# Model 1: ~marks + x + y + I(x^2) + I(x * y) + I(y^2) + I(x^3) + I(x^2 * y) + I(x
#  * y^2) + I(y^3)         MultiStrauss
# Model 2: ~marks + x + y + I(x^2) + I(x * y) + I(y^2) + I(x^3) + I(x^2 * y) + I(x
#  * y^2) + I(y^3) + marks:x       MultiStrauss
# Model 3: ~marks + x + y + I(x^2) + I(x * y) + I(y^2) + I(x^3) + I(x^2 * y) + I(x
#  * y^2) + I(y^3) + marks:x + marks:y     MultiStrauss
# Model 4: ~marks + x + y + I(x^2) + I(x * y) + I(y^2) + I(x^3) + I(x^2 * y) + I(x
#  * y^2) + I(y^3) + marks:x + marks:y + marks:I(x^2)      MultiStrauss
# Model 5: ~marks + x + y + I(x^2) + I(x * y) + I(y^2) + I(x^3) + I(x^2 * y) + I(x
#  * y^2) + I(y^3) + marks:x + marks:y + marks:I(x^2) + marks:I(x * y)     MultiSt
# rauss
# Model 6: ~marks + x + y + I(x^2) + I(x * y) + I(y^2) + I(x^3) + I(x^2 * y) + I(x
#  * y^2) + I(y^3) + marks:x + marks:y + marks:I(x^2) + marks:I(x * y) + marks:I(y
# ^2)      MultiStrauss
# Model 7: ~marks + x + y + I(x^2) + I(x * y) + I(y^2) + I(x^3) + I(x^2 * y) + I(x
#  * y^2) + I(y^3) + marks:x + marks:y + marks:I(x^2) + marks:I(x * y) + marks:I(y
# ^2) + marks:I(x^3)       MultiStrauss
# Model 8: ~marks + x + y + I(x^2) + I(x * y) + I(y^2) + I(x^3) + I(x^2 * y) + I(x
#  * y^2) + I(y^3) + marks:x + marks:y + marks:I(x^2) + marks:I(x * y) + marks:I(y
# ^2) + marks:I(x^3) + marks:I(x^2 * y)    MultiStrauss
# Model 9: ~marks + x + y + I(x^2) + I(x * y) + I(y^2) + I(x^3) + I(x^2 * y) + I(x
#  * y^2) + I(y^3) + marks:x + marks:y + marks:I(x^2) + marks:I(x * y) + marks:I(y
# ^2) + marks:I(x^3) + marks:I(x^2 * y) + marks:I(x * y^2)         MultiStrauss   
# Model 10: ~marks * (x + y + I(x^2) + I(x * y) + I(y^2) + I(x^3) + I(x^2 * y) + I
# (x * y^2) + I(y^3))      MultiStrauss
#    Npar Df AdjDeviance  Pr(>Chi)
# 1    13
# 2    14  1       2.609  0.106256
# 3    15  1      76.083 < 2.2e-16 ***
# 4    16  1      -2.379  0.122986    
# 5    17  1      18.187 2.003e-05 ***
# 6    18  1      35.328 2.786e-09 ***
# 7    19  1       6.796  0.009137 **
# 8    20  1      15.250 9.419e-05 ***
# 9    21  1      36.968 1.201e-09 ***
# 10   22  1     117.126 < 2.2e-16 ***
# ---
# Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

#2.4 late stage
#Cross-type K-function
Kcross.late.FHG<-envelope(late.FHG.multitype,fun=Kcross,nsim=99)
jpeg(here("Output", "2_analysis", "Kcross_late_FHG.jpg"),width = 7, height = 7, units = "in", res = 1200)
plot(Kcross.late.FHG)
dev.off()

# PPM
late.hf <- hyperframe(p=late.FHG.multitype)
late.hf$m <- with(late.hf,map.im)
late.hf$r <- with(late.hf,resource.im)
late.mppm1 <- mppm(p~m+r,data=late.hf)
late.dist1 <- distfun(late.FHG.multitype)
summary(late.dist1)
# Distance function for point pattern
# defined in a polygonal window inside the frame [-4, 1] x [-1, 4] units
# 
# Distance function values:
#         range = [0.01381068, 4.241228]
#         mean = 1.403091
late.dist2 <- pairdist(late.FHG.multitype)
summary(late.dist2)
max(late.dist2)
# [1] 4.107511
mean(late.dist2)
# [1] 1.803547
ks.test(middle.late.dist2,late.dist2)
# Asymptotic two-sample Kolmogorov-Smirnov test
# 
# data:  middle.late.dist2 and late.dist2
# D = 0.13416, p-value = 0.3541
# alternative hypothesis: two-sided
# 
# Warning message:
# In ks.test.default(middle.late.dist2, late.dist2) :
#   p-value will be approximate in the presence of ties
ks.test(early.dist.2,late.dist2)
# Asymptotic two-sample Kolmogorov-Smirnov test
# 
# data:  early.dist.2 and late.dist2
# D = 0.12202, p-value = 0.4726
# alternative hypothesis: two-sided
# 
# Warning message:
# In ks.test.default(early.dist.2, late.dist2) :
#   p-value will be approximate in the presence of ties
ks.test(middle.early.dist2,late.dist2)
# Asymptotic two-sample Kolmogorov-Smirnov test
# 
# data:  middle.early.dist2 and late.dist2
# D = 0.12993, p-value = 0.385
# alternative hypothesis: two-sided
# 
# Warning message:
# In ks.test.default(middle.early.dist2, late.dist2) :
#   p-value will be approximate in the presence of ties

#multitype strauss
late.RR <- data.frame(R=seq(0.01, 4.11, by=0.01))
late.MS <- function(R) { MultiStrauss(radii=diag(c(R,R))) }
late.pm <- profilepl(late.RR, late.MS,late.FHG.multitype~marks*polynom(x,y,2), rbord = NULL,correction="Ripley")
warnings(late.pm)
# Warning messages:
# 1: glm.fit: algorithm did not converge Error in cat("1: glm.fit: algorithm did not c
# onverge", list(param = list( :
#   argument 2 (type 'list') cannot be handled by 'cat'
late.pm
# profile log pseudolikelihood for model:
# ppm(late.FHG.multitype ~ marks * polynom(x,  y,  2),  correction = "Ripley",    
#    interaction = late.MS,  rbord = NULL)
# interaction: Multitype Strauss process
# irregular parameter: R in [0.01, 4.11]
# optimum value of irregular parameter:  R = 1.6
plot(late.pm)
jpeg(here("Output", "2_analysis", "late_pm.jpg"),width = 7, height = 7, units = "in", res = 1200)
plot(late.pm)
dev.off()
as.ppm(late.pm)
# Nonstationary Multitype Strauss process
# Fitted to point pattern dataset 'late.FHG.multitype'
# 
# Possible marks: 'Farmer' and 'Hunter_gatherer'
# 
# Log trend:  ~marks * (x + y + I(x^2) + I(x * y) + I(y^2))
# 
# Fitted trend coefficients:
#                   (Intercept)          marksHunter_gatherer 
#                      9.680035                    -84.608875
#                             x                             y
#                     15.293028                     18.035504
#                        I(x^2)                      I(x * y)
#                     -9.332770                    -15.249145
#                        I(y^2)        marksHunter_gatherer:x 
#                     -9.219586                    -58.106315
#        marksHunter_gatherer:y   marksHunter_gatherer:I(x^2)
#                    -42.842600                    -11.125352
# marksHunter_gatherer:I(x * y)   marksHunter_gatherer:I(y^2)
#                     14.664136                     21.698545
# 
# 2 types of points
# Possible types:
# [1] Farmer          Hunter_gatherer
# Interaction radii:
#                 Farmer Hunter_gatherer
# Farmer             1.6              NA
# Hunter_gatherer     NA             1.6
# Fitted interaction parameters gamma_ij
#                    Farmer Hunter_gatherer
# Farmer          0.0035656              NA
# Hunter_gatherer        NA               0
# 
# Relevant coefficients:
#                   markFarmerxFarmer markHunter_gathererxHunter_gatherer
#                           -5.636435                          -24.057489
# 
# For standard errors, type coef(summary(x))
anova((as.ppm(late.pm)),test="LR")

#area-interaction
late.RR.AI <- data.frame(r=seq(0.01, 4.11, by=0.01))
late.pm.AI <- profilepl(late.RR.AI,AreaInter,late.FHG.multitype~marks+polynom(x,y,2),rbord=NULL, correction="Ripley")
late.pm.AI
# profile log pseudolikelihood for model:
# ppm(late.FHG.multitype ~ marks + polynom(x,  y,  2),  correction = "Ripley",    
#    interaction = AreaInter,  rbord = NULL)
# interaction: Area-interaction process
# irregular parameter: r in [0.01, 4.11]
# optimum value of irregular parameter:  r = 2.52
warnings()
# Warning messages:
# 1: In pot(X, U, EqualPairs, pars, correction, ...) :
#   Correction 'isotropic' is not supported and was ignored
jpeg(here("Output", "2_analysis", "late_pm_AI.jpg"),width = 7, height = 7, units = "in", res = 1200)
plot(late.pm.AI)
dev.off()

as.ppm(late.pm.AI)
# Nonstationary Area-interaction process
# Fitted to point pattern dataset 'late.FHG.multitype'
# 
# Possible marks: 'Farmer' and 'Hunter_gatherer'
# 
# Log trend:  ~marks + (x + y + I(x^2) + I(x * y) + I(y^2))
# 
# Fitted trend coefficients:
#          (Intercept) marksHunter_gatherer                    x
#          155.5660442           -0.8782445           11.3842737
#                    y               I(x^2)             I(x * y)
#            6.8573047           -5.8690519           -7.2092667
#               I(y^2)
#           -2.1871316 
# 
# Disc radius:    2.52
# Fitted interaction parameter eta:       1.38197e-70
# 
# Relevant coefficients:
# Interaction
#   -160.8574
# 
# For standard errors, type coef(summary(x))

#the end




