## Author: Peter J. Galante
## This script performs niche comparisons using Schoener's D-value and Proportion Predicted Present (PPP)
## It also calculates niche similarity and equivalency using the ecospat package
## Inputs: from "canescens.R", or "sumichrasti.R"

###########################################################################################################################
###########################################################################################################################
###########################################    NICHE COMPARISONS    #######################################################
###########################################################################################################################
###########################################################################################################################
###########################################################################################################################
library(spThin);library(ENMeval);library(rgdal);library(rgeos)
## D-value comparisons between homochrous and sumichrasti for each time period
csCurrentFull <- raster("canescens/CanescensModel/CurrentCanFull.tif")
csCurrent <- raster("canescens/CanescensModel/CurrentCan.tif")
csMH <- raster("canescens/CanescensModel/MHCan.tif")
csLGM <- raster("canescens/CanescensModel/LGMCan.tif")
csLIG <- raster("canescens/CanescensModel/LIGCan.tif")
##
scCurrentFull<- raster("sumichrasti/SumichrastiModel/CurrentSumFull.tif")
scCurrent <- raster("sumichrasti/SumichrastiModel/CurrentSum.tif")
scMH <- raster("sumichrasti/SumichrastiModel/MHSum.tif")
scLGM <- raster("sumichrasti/SumichrastiModel/LGMSum.tif")
scLIG <-raster("sumichrasti/SumichrastiModel/LIGSum.tif")
# D values for each time period
Dcurrent <- nicheOverlap(csCurrent, scCurrent, stat = "D")
DMH <- nicheOverlap(csMH, scMH, stat = "D")
DLGM <- nicheOverlap(csLGM, scLGM, stat = "D")
DLIG <- nicheOverlap(csLIG, scLIG, stat = "D")
## bargraph of D through time
e <- extent(c(-82.6, -77.5, 6, 11))
curcan <- crop(csCurrent, e)
MHcan <- crop(csMH, e)
LGMcan <- crop(csLGM, e)
LIGcan <- crop(csLIG, e)
cursum <- crop(scCurrent, e)
MHsum <- crop(scMH, e)
LGMsum <- crop(scLGM, e)
LIGsum <- crop(scLIG, e)

dCur<-nicheOverlap(curcan, cursum, stat = "D")
dMH<-nicheOverlap(MHcan, MHsum, stat = "D")
dLGM<-nicheOverlap(LGMcan, LGMsum, stat = "D")
dLIG<-nicheOverlap(LIGcan, LIGsum, stat = "D")

## PPP: get total potential pixels (loaded from timeSlicePlot.R)
nCurrent <- ncell(curcan[!is.na(curcan)])
PPPcurCan <- ncell(CurCan[!is.na(CurCan)]) / nCurrent *100
PPPcurSum <- ncell(CurSum[!is.na(CurSum)]) / nCurrent *100
nMH <- ncell(MHcan[!is.na(MHcan)])
PPPMHcan <- ncell(MHCan[!is.na(MHCan)]) / nMH *100
PPPMHsum <- ncell(MHSum[!is.na(MHSum)]) / nMH *100
nLGM <- ncell(LGMcan[!is.na(LGMcan)])
PPPLGMcan <- ncell(LGMCan[!is.na(LGMCan)]) / nLGM *100
PPPLGMsum <- ncell(LGMSum[!is.na(LGMSum)]) / nLGM *100
nLIG <- ncell(LIGcan[!is.na(LIGcan)])
PPPLIGcan <- ncell(LIGCan[!is.na(LIGCan)]) / nLIG *100
PPPLIGsum <- ncell(LIGSum[!is.na(LIGSum)]) / nLIG *100

##### Ecospat niche equivalency and similarity tests
library(ecospat)
# env layers
env <- stack(list.files("/layers", pattern = "\\.tif$", full.names = T))
# Occurrence points
occd1 <- read.csv("/sumichrastiThin1/thinned_data_thin4.csv")
occd2 <- read.csv("/canescensThin1/thinned_data_thin1.csv")
# background points
bg1 <- read.csv("sumichrasti/sumichrastiBackg.csv")
colnames(bg1) <- c("Long", "Lat")
bg2 <- read.csv("canescens/canescensBackg.csv")
colnames(bg2) <- c("Long", "Lat")
# environmental data
extract1 <- na.omit(cbind(occd1[,2:3], extract(env, occd1[,2:3]), rep(1, nrow(occd1))))
extract2 <- na.omit(cbind(occd2[,2:3], extract(env, occd2[,2:3]), rep(1, nrow(occd2))))
colnames(extract1)[ncol(extract1)] = "occ"
colnames(extract2)[ncol(extract2)] = "occ"
extbg1 <- na.omit(cbind(bg1, extract(env, bg1), rep(0, nrow(bg1))))
extbg2 <- na.omit(cbind(bg2, extract(env, bg2), rep(0, nrow(bg2))))
colnames(extbg1)[ncol(extbg1)] = "occ"
colnames(extbg2)[ncol(extbg2)] = "occ"
# merge occ and bg data
dat1 <- rbind(extract1, extbg1)
dat2 <- rbind(extract2, extbg2)

pca.env <- dudi.pca(
  rbind(dat1, dat2)[,3:21],
  scannf = F,
  nf = 2)
# Variable contribution
ecospat.plot.contrib(contrib = pca.env$co, eigen = pca.env$eig)
scores.globclim <- pca.env$li
scores.sp1 <- suprow(pca.env, extract1[which(extract1[,22]==1), 3:21])$li
scores.sp2 <- suprow(pca.env, extract2[which(extract2[,22]==1), 3:21])$li
scores.clim1 <- suprow(pca.env, dat1[,3:21])$li
scores.clim2 <- suprow(pca.env, dat2[,3:21])$li
grid.clim1 <- ecospat.grid.clim.dyn(
  glob = scores.globclim,
  glob1 = scores.clim1,
  sp = scores.sp1,
  R = 100,
  th.sp = 0)
grid.clim2 <- ecospat.grid.clim.dyn(
  glob = scores.globclim,
  glob1 = scores.clim2,
  sp = scores.sp2,
  R = 100,
  th.sp = 0)
D.overlap <- ecospat.niche.overlap(grid.clim1, grid.clim2, cor = T)$D
### Niche Equivalency test
# H1: Is the overlap between sp1 and sp2 niche higher than two random niches?
eq.test <- ecospat.niche.equivalency.test(grid.clim1, grid.clim2, rep = 1000, alternative = "greater")
# H1: Is the overlap between sp1 and sp2 higher than when sp2 niche is randomly introduced in the area of sp1
sim.test <- ecospat.niche.similarity.test(grid.clim1, grid.clim2, rep = 1000, alternative = "greater", rand.type = 2)

# Plotting
# Bars represent the frequency of niche overlap D value. Red line represents observed niche overlap D value. p-value is significance test
# The niche equivalency test evaluates, through random permutations of occurrences between ranges, if two niches are significantly different from each other.
pdf("birdEquivalency.pdf")
ecospat.plot.overlap.test(eq.test, "D", "Equivalency") # 0.26; p=0.001
dev.off()
# Bars represent the frequency of niche overlap D value. Red line represents observed niche overlap D value. p-value is significance test
# The niche similarity test assesses, through random shifts of the niches within available conditions in the study area, if two niches are more or less similar than expected by chance
pdf("birdSimilarity.pdf")
ecospat.plot.overlap.test(sim.test, "D", "Similarity") # 0.26; p=0.04695
dev.off()

