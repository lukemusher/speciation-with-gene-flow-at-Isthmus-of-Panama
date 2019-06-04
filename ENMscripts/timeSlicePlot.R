## Author: Peter J. Galante
## This script loads thresholded rasters, creates barplots of Schoener's D-values, calculates proportion concordance,
## creates similarity histograms, calculates neighborhood statistics in both cloglog and raw outputs.
## Inputs from: "canescens.R", or "sumichrasti.R"

library(raster); library(ggplot2); library(ggspatial); library(rasterVis); library(sf);library(dismo)
e <- extent(c(-82.6, -77.5, 6, 11))
# Load canescens ETSS thresholded models
CurCan <- crop(raster('canescens/ETSS/curCanThresh.tif'), e)
LGMCan <- crop(raster('canescens/ETSS/LGMCanThresh.tif'), e)
LIGCan <- crop(raster('canescens/ETSS/LIGCanThresh.tif'), e)
MHCan <- crop(raster('canescens/ETSS/MHCanThresh.tif'), e)
# Load sumichrasti ETSS thresholded models
CurSum <- crop(raster('sumichrasti/ETSS/curSumThresh.tif'), e)
LGMSum <- crop(raster('sumichrasti/ETSS/LGMSumThresh.tif'), e)
LIGSum <- crop(raster('sumichrasti/ETSS/LIGSumThresh.tif'), e)
MHSum <- crop(raster('sumichrasti/ETSS/MHSUMThresh.tif'), e)
## Load polygon outline
poly <- crop(readOGR(dsn = "layers/10mStatesProvince/ne_10m_admin_1_states_provinces.shp"), e)
poly <- st_as_sf(poly)
poly <- st_union(poly)

## Set up plot area
dev.off()
png("bothBinary.png", height=420, res = 100, width = 1000)
pdf("bothBinary.pdf", height = 4, width=10)
plot.new()
par(mfrow=(c(1,4)), oma = c(4,4,2,2), mar = c(1,1,0,0))
plot.window(xlim = extent(CurCan)[1:2], ylim = extent(CurCan)[3:4])

plot(CurCan, col=heat.colors(1, alpha=1), legend=F)
plot(CurSum, col=heat.colors(2, alpha=0.65), add=T, legend=F)
mtext("Current", side=3, line=0.25, cex=1.5, col="black", outer=F)
legend("topright", legend=c(expression(italic("P. homochrous")), expression(italic("P. aglaiae"))), fill = c(heat.colors(2, alpha=1)[[1]],heat.colors(2, alpha=0.65)[[2]]))
plot(poly, add=T)
legend("bottom", legend = "Schoener's D-value = 83.28% \n PPP P. homochrous: 71.12% \n PPP P. aglaiae: 94.97%", inset = 0.03, bty = "n")

plot(MHCan, col=heat.colors(1, alpha=1), legend=F)
plot(MHSum, col=heat.colors(2, alpha=0.65), add=T, legend=F)
mtext("Mid-Holocene", side=3, line=0.25, cex=1.5, col="black", outer=F)
legend("topright", legend=c(expression(italic("P. homochrous")), expression(italic("P. aglaiae"))), fill = c(heat.colors(2, alpha=1),heat.colors(2, alpha=0.65)))
plot(poly, add=T)
legend("bottom", legend = "Schoener's D-value = 77.67% \n PPP P. homochrous: 59.87% \n PPP P. aglaiae: 93.54%", inset = 0.03, bty = "n")

plot(LGMCan, col=heat.colors(1, alpha=1), legend=F)
plot(LGMSum, col=heat.colors(2, alpha=0.65), add=T, legend=F)
mtext("LGM", side=3, line=0.25, cex=1.5, col="black", outer=F)
legend("topright", legend=c(expression(italic("P. homochrous")), expression(italic("P. aglaiae"))), fill = c(heat.colors(2, alpha=1)[[1]],heat.colors(2, alpha=0.65)[[2]]))
plot(poly, add=T)
legend("bottom", legend = "Schoener's D-value = 70.34% \n PPP P. homochrous: 37.17% \n PPP P. aglaiae: 89.06%", inset = 0.03, bty = "n")

plot(LIGCan, col=heat.colors(1, alpha=1), legend=F)
plot(LIGSum, col=heat.colors(2, alpha=0.65), add=T, legend=F)
mtext("LIG", side=3, line=0.25, cex=1.5, col="black", outer=F)
legend("topright", legend=c(expression(italic("P. homochrous")), expression(italic("P. aglaiae"))), fill = c(heat.colors(2, alpha=1)[[1]],heat.colors(2, alpha=0.65)[[2]]))
plot(poly, add=T)
legend("bottom", legend = "Schoener's D-value = 81.71% \n PPP P. homochrous: 62.08% \n PPP P. aglaiae: 87.40%", inset = 0.03, bty = "n")
dev.off()

##########################################
###  Barplot of d-values through time for both species
library(spThin);library(ENMeval);library(rgdal);library(rgeos)
e <- extent(c(-82.6, -77.5, 6, 11))
## D-value comparisons between homochrous and sumichrasti for each time period
csCurrent <- crop(raster("canescens/CanescensModel/CurrentCan.tif"), e)
csMH <- crop(raster("canescens/CanescensModel/MHCan.tif"), e)
csLGM <- crop(raster("canescens/CanescensModel/LGMCan.tif"), e)
csLIG <- crop(raster("canescens/CanescensModel/LIGCan.tif"), e)

##
scCurrent <- crop(raster("sumichrasti/SumichrastiModel/CurrentSum.tif"), e)
scMH <- crop(raster("sumichrasti/SumichrastiModel/MHSum.tif"), e)
scLGM <- crop(raster("sumichrasti/SumichrastiModel/LGMSum.tif"), e)
scLIG <-crop(raster("sumichrasti/SumichrastiModel/LIGSum.tif"), e)
# D values for each time period
Dcurrent <- nicheOverlap(csCurrent, scCurrent, stat = "D")
DMH <- nicheOverlap(csMH, scMH, stat = "D")
DLGM <- nicheOverlap(csLGM, scLGM, stat = "D")
DLIG <- nicheOverlap(csLIG, scLIG, stat = "D")
## bargraph of D through time
ddf<-data.frame(matrix(data = c(Dcurrent, DMH, DLGM, DLIG)))
ddf$TimePeriod <- c("Current", "Mid-Holocene", "LGM", "LIG")
colnames(ddf) <- c("Dvalues", "TimePeriod")

pdf("DvaluesTime.pdf")
png("DvaluesTime.png")
ggplot(data = ddf, aes(x=TimePeriod, y=Dvalues))+
  geom_bar(stat="identity", width=0.75, fill="steelblue")+
  geom_text(aes(label=paste0(round(ddf$Dvalues, 4) * 100, "%")), vjust = 1.5, color = "white", size = 3.5)+
  theme_minimal()+
  ylim(0,1)+ ylab("Schoener's D-value") + xlab("Time Period")
dev.off()


#################################################################
## niche projection plots
library(rasterVis)
canescens <- raster('canescens/CanescensModel/CurrentCanFull.tif')
sumichrasti <- raster('sumichrasti/SumichrastiModel/CurrentSumFull.tif')
plot(canescens)
plot(sumichrasti)

mytheme <- rasterTheme(region = rev(brewer.pal(9, "RdYlBu")))
png("Phomochrous.png")
levelplot(canescens, layer=1, par.settings=mytheme, main = expression(italic("P. homochrous")), margin = F)
dev.off()
png("Paglaiae.png")
levelplot(sumichrasti, layer=1, par.settings=mytheme, main = expression(italic("P. aglaiae")), margin = F)
dev.off()

#########################################################################
####################  proportion concordance  ###########################
e <- extent(c(-82.6, -77.5, 6, 11))
# Load canescens ETSS thresholded models
CurCan <- crop(raster('canescens/ETSS/curCanThresh.tif'), e)
LGMCan <- crop(raster('canescens/ETSS/LGMCanThresh.tif'), e)
LIGCan <- crop(raster('canescens/ETSS/LIGCanThresh.tif'), e)
MHCan <- crop(raster('canescens/ETSS/MHCanThresh.tif'), e)
# Load sumichrasti ETSS thresholded models
CurSum <- crop(raster('sumichrasti/ETSS/curSumThresh.tif'), e)
LGMSum <- crop(raster('sumichrasti/ETSS/LGMSumThresh.tif'), e)
LIGSum <- crop(raster('sumichrasti/ETSS/LIGSumThresh.tif'), e)
MHSum <- crop(raster('sumichrasti/ETSS/MHSUMThresh.tif'), e)


CurCan[is.na(CurCan)] <- 0
CurSum[is.na(CurSum)] <- 0
current <- CurCan + CurSum
current[current == 0] <- NA
propCur <- ncell(current[current == 2]) / ncell(current[!is.na(current)]) * 100
CurOverlap <- tapply(area(current), current[], sum) / 1000

LGMCan[is.na(LGMCan)] <- 0
LGMSum[is.na(LGMSum)] <- 0
LGM <- LGMCan + LGMSum
LGM[LGM == 0] <- NA
propLGM <- ncell(LGM[LGM == 2]) / ncell(LGM[!is.na(LGM)]) * 100
LGMOverlap <- tapply(area(LGM), LGM[], sum) / 1000

LIGCan[is.na(LIGCan)] <- 0
LIGSum[is.na(LIGSum)] <- 0
LIG <- LIGCan + LIGSum
LIG[LIG == 0] <- NA
propLIG <- ncell(LIG[LIG == 2]) / ncell(LIG[!is.na(LIG)]) * 100
LIGOverlap <- tapply(area(LIG), LIG[], sum) / 1000

MHCan[is.na(MHCan)] <- 0
MHSum[is.na(MHSum)] <- 0
MH <- MHCan + MHSum
MH[MH == 0]<- NA
propMH <- ncell(MH[MH == 2]) / ncell(MH[!is.na(MH)]) * 100
MHOverlap <- tapply(area(MH), MH[], sum) / 1000

props <- matrix(nrow = 2, ncol = 4)
colnames(props) <- c("current", "MH", "LGM", "LIG")
rownames(props) <- c("area concordance (km)", "Proportion concordance (%)")
props[1,] <- cbind(CurOverlap[[2]], MHOverlap[[2]], LGMOverlap[[2]], LIGOverlap[[2]])
props[2,] <- cbind(propCur, propMH, propLGM, propLIG)
write.csv(props, "proportionsOverlap.csv")


#########################################################################
####################  similarity histograms  ###########################
e <- extent(c(-82.6, -77.5, 6, 11))
csCurrent <- crop(raster("canescens/CanescensModel/CurrentCan.tif"), e)
csMH <- crop(raster("canescens/CanescensModel/MHCan.tif"), e)
csLGM <- crop(raster("canescens/CanescensModel/LGMCan.tif"), e)
csLIG <- crop(raster("canescens/CanescensModel/LIGCan.tif"), e)

##
scCurrent <- crop(raster("sumichrasti/SumichrastiModel/CurrentSum.tif"), e)
scMH <- crop(raster("sumichrasti/SumichrastiModel/MHSum.tif"), e)
scLGM <- crop(raster("sumichrasti/SumichrastiModel/LGMSum.tif"), e)
scLIG <-crop(raster("sumichrasti/SumichrastiModel/LIGSum.tif"), e)

dCur<-nicheOverlap(csCurrent, scCurrent, stat = "D")
dMH<-nicheOverlap(csMH, scMH, stat = "D")
dLGM<-nicheOverlap(csLGM, scLGM, stat = "D")
dLIG<-nicheOverlap(csLIG, scLIG, stat = "D")

#################  bootstrap niche overlaps  #######################3
## current 50, 250, 500km diameters
cur.pts <- randomPoints((csCurrent + scCurrent), 1000)
cur.pts.buf <- apply(cur.pts, 1, function(x) buffer(SpatialPoints(as.data.frame(rbind(x))), width = 250000))
scmasks <- lapply(cur.pts.buf, function(x) mask(scCurrent, x))
csmasks <- lapply(cur.pts.buf, function(x) mask(csCurrent, x))
laps <- list()
for (i in 1:length(scmasks)){
 laps[[i]] <- nicheOverlap(x = scmasks[[i]], csmasks[[i]], stat = "D", mask = T)
}
pdf('timeSliceBootstrapCurrent500.pdf')
hist(unlist(laps), main = "current")
abline(v= dCur, col = 'red')
dev.off()
write.csv(unlist(laps),'500kmDvals100RepsCurrent.csv')
rm(list = c("scmasks", "csmasks", "laps"))
gc()
##  MH
MH.pts <- randomPoints((csMH + scMH), 1000)
MH.pts.buf <- apply(MH.pts, 1, function(x) buffer(SpatialPoints(as.data.frame(rbind(x))), width = 250000))
MHscmasks <- lapply(MH.pts.buf, function(x) mask(scMH, x))
MHcsmasks <- lapply(MH.pts.buf, function(x) mask(csMH, x))
MHlaps <- list()
for (i in 1:length(MHscmasks)){
  MHlaps[[i]] <- nicheOverlap(x = MHscmasks[[i]], MHcsmasks[[i]], stat = "D", mask = T)
}
pdf('timeSliceBootstrapMH500.pdf')
hist(unlist(MHlaps), main = "Mid-Holocene")
abline(v= dMH, col = 'red')
dev.off()
write.csv(unlist(MHlaps),'500kmDvals100RepsMH.csv')
rm(list = c("MHscmasks", "MHcsmasks", "MHlaps"))
gc()

##  LGM
LGM.pts <- randomPoints((csLGM + scLGM), 1000)
LGM.pts.buf <- apply(LGM.pts, 1, function(x) buffer(SpatialPoints(as.data.frame(rbind(x))), width = 250000))
LGMscmasks <- lapply(LGM.pts.buf, function(x) mask(scLGM, x))
LGMcsmasks <- lapply(LGM.pts.buf, function(x) mask(csLGM, x))
LGMlaps <- list()
for (i in 1:length(LGMscmasks)){
  LGMlaps[[i]] <- nicheOverlap(x = LGMscmasks[[i]], LGMcsmasks[[i]], stat = "D", mask = T)
}
pdf('timeSliceBootstrapLGM500.pdf')
hist(unlist(LGMlaps), main = "LGM")
abline(v= dLGM, col = 'red')
dev.off()
write.csv(unlist(LGMlaps),'500kmDvals100RepsLGM.csv')
rm(list = c("LGMscmasks", "LGMcsmasks", "LGMlaps"))
gc()

##  LIG
LIG.pts <- randomPoints((csLIG + scLIG), 1000)
LIG.pts.buf <- apply(LIG.pts, 1, function(x) buffer(SpatialPoints(as.data.frame(rbind(x))), width = 250000))
LIGscmasks <- lapply(LIG.pts.buf, function(x) mask(scLIG, x))
LIGcsmasks <- lapply(LIG.pts.buf, function(x) mask(csLIG, x))
LIGlaps <- list()
for (i in 1:length(LIGscmasks)){
  LIGlaps[[i]] <- nicheOverlap(x = LIGscmasks[[i]], LIGcsmasks[[i]], stat = "D", mask = T)
}
pdf('timeSliceBootstrapLIG500.pdf')
hist(unlist(LIGlaps), main="LIG")
abline(v= dLIG, col = 'red')
dev.off()
write.csv(unlist(LIGlaps),'500kmDvals100RepsLIG.csv')
rm(list = c("LIGscmasks", "LIGcsmasks", "LIGlaps"))
gc()

#########################################################################################
############################  Neigborhood % concordances  ###############################
e <- extent(c(-82.6, -77.5, 6, 11))
csCurrent <- crop(raster("canescens/ETSS/curCanThresh.tif"), e)
csMH <- crop(raster("canescens/ETSS/MHCanThresh.tif"), e)
csLGM <- crop(raster("canescens/ETSS/LGMCanThresh.tif"), e)
csLIG <- crop(raster("canescens/ETSS/LIGCanThresh.tif"), e)

scCurrent <- crop(raster("sumichrasti/ETSS/curSumThresh.tif"), e)
scMH <- crop(raster("sumichrasti/ETSS/MHSUMThresh.tif"), e)
scLGM <- crop(raster("sumichrasti/ETSS/LGMSumThresh.tif"), e)
scLIG <-crop(raster("sumichrasti/ETSS/LIGSumThresh.tif"), e)
#########################################
## % concordance of neighborhoods of 50, 250, 500km
cur.Both <- sum(scCurrent, csCurrent, na.rm=T)
cur.Both[cur.Both == 0] <- NA
cur.pts <- randomPoints(cur.Both, 1000)
cur.pts.buf <- apply(cur.pts, 1, function(x) buffer(SpatialPoints(as.data.frame(rbind(x))), width = 25000))
scmasks <- lapply(cur.pts.buf, function(x) mask(scCurrent, x))
csmasks <- lapply(cur.pts.buf, function(x) mask(csCurrent, x))
prop.current <- list()
for (i in 1:length(scmasks)){
  scmasks[[i]][is.na(scmasks[[i]])] <- 0
  csmasks[[i]][is.na(csmasks[[i]])] <- 0
  current <- scmasks[[i]] + csmasks[[i]]
  current[current == 0] <- NA
  prop.current[[i]] <- ncell(current[current == 2]) / ncell(current[!is.na(current)]) * 100
}
write.csv(unlist(prop.current), "50kmProportion1000RepsCurrent.csv")
rm(cur.Both, cur.pts, scmasks, csmasks, prop.current)
gc()

MH.Both <- sum(scMH, csMH, na.rm=T)
MH.Both[MH.Both == 0] <- NA
MH.pts <- randomPoints(MH.Both, 1000)
MH.pts.buf <- apply(MH.pts, 1, function(x) buffer(SpatialPoints(as.data.frame(rbind(x))), width = 25000))
MHscmasks <- lapply(MH.pts.buf, function(x) mask(scCurrent, x))
MHcsmasks <- lapply(MH.pts.buf, function(x) mask(csCurrent, x))
prop.MH <- list()
for (i in 1:length(MHscmasks)){
  MHscmasks[[i]][is.na(MHscmasks[[i]])] <- 0
  MHcsmasks[[i]][is.na(MHcsmasks[[i]])] <- 0
  MH <- MHscmasks[[i]] + MHcsmasks[[i]]
  MH[MH == 0] <- NA
  prop.MH[[i]] <- ncell(MH[MH == 2]) / ncell(MH[!is.na(MH)]) * 100
}
write.csv(unlist(prop.MH), "50kmProportion1000RepsMH.csv")
rm(MH.Both, MH.pts, MHscmasks, MHcsmasks, prop.MH)
gc()

LGM.Both <- sum(scLGM, csLGM, na.rm=T)
LGM.Both[LGM.Both == 0] <- NA
LGM.pts <- randomPoints(LGM.Both, 1000)
LGM.pts.buf <- apply(LGM.pts, 1, function(x) buffer(SpatialPoints(as.data.frame(rbind(x))), width = 25000))
LGMscmasks <- lapply(LGM.pts.buf, function(x) mask(scCurrent, x))
LGMcsmasks <- lapply(LGM.pts.buf, function(x) mask(csCurrent, x))
prop.LGM <- list()
for (i in 1:length(LGMscmasks)){
  LGMscmasks[[i]][is.na(LGMscmasks[[i]])] <- 0
  LGMcsmasks[[i]][is.na(LGMcsmasks[[i]])] <- 0
  LGM <- LGMscmasks[[i]] + LGMcsmasks[[i]]
  LGM[LGM == 0] <- NA
  prop.LGM[[i]] <- ncell(LGM[LGM == 2]) / ncell(LGM[!is.na(LGM)]) * 100
}
write.csv(unlist(prop.LGM), "50kmProportion1000RepsLGM.csv")
rm(LGM.Both, LGM.pts, LGMscmasks, LGMcsmasks, prop.LGM)
gc()

LIG.Both <- sum(scLIG, csLIG, na.rm=T)
LIG.Both[LIG.Both == 0] <- NA
LIG.pts <- randomPoints(LIG.Both, 1000)
LIG.pts.buf <- apply(LIG.pts, 1, function(x) buffer(SpatialPoints(as.data.frame(rbind(x))), width = 25000))
LIGscmasks <- lapply(LIG.pts.buf, function(x) mask(scCurrent, x))
LIGcsmasks <- lapply(LIG.pts.buf, function(x) mask(csCurrent, x))
prop.LIG <- list()
for (i in 1:length(LIGscmasks)){
  LIGscmasks[[i]][is.na(LIGscmasks[[i]])] <- 0
  LIGcsmasks[[i]][is.na(LIGcsmasks[[i]])] <- 0
  LIG <- LIGscmasks[[i]] + LIGcsmasks[[i]]
  LIG[LIG == 0] <- NA
  prop.LIG[[i]] <- ncell(LIG[LIG == 2]) / ncell(LIG[!is.na(LIG)]) * 100
}
write.csv(unlist(prop.LIG), "50kmProportion1000RepsLIG.csv")
rm(LIG.Both, LIG.pts, LIGscmasks, LIGcsmasks, prop.LIG)
gc()

##################################################################################################################
##################################################################################################################
######################################  Bootstrapping raw joint probabilities  ################################### 50, 250, 500km
## Load raw preds - created in "sumichrasti.R" or "canescens.R"
# Canescens
CurcanRaw <- raster('canescens/CanescensModel/CurrentCanRaw.tif')
MHcanRaw <- raster('canescens/CanescensModel/MHCanRaw.tif')
LGMcanRaw <- raster('canescens/CanescensModel/LGMCanRaw.tif')
LIGcanRaw <- raster('canescens/CanescensModel/LIGCanRaw.tif')
# aglaiae
CursumRaw <- raster('sumichrasti/SumichrastiModel/CurrentSumRaw.tif')
MHsumRaw <- raster('sumichrasti/SumichrastiModel/MHSumRaw.tif')
LGMsumRaw <- raster('sumichrasti/SumichrastiModel/LGMSumRaw.tif')
LIGsumRaw <- raster('sumichrasti/SumichrastiModel/LIGSumRaw.tif')
# current
cur.pts <- randomPoints((CurcanRaw + CursumRaw), 1000)
cur.pts.buf <- apply(cur.pts, 1, function(x) buffer(SpatialPoints(as.data.frame(rbind(x))), width = 500000))
scmasks <- lapply(cur.pts.buf, function(x) mask(CursumRaw, x))
csmasks <- lapply(cur.pts.buf, function(x) mask(CurcanRaw, x))
CurProbs <- list()
for (i in 1:length(scmasks)){
  CurProbs[[i]] <- mean(values(scmasks[[i]] * csmasks[[i]]), na.rm=T)
}
write.csv(unlist(CurProbs), "500kmProbabilitiesCurrent.csv")
# MH
MH.pts <- randomPoints((MHcanRaw + MHsumRaw), 1000)
MH.pts.buf <- apply(MH.pts, 1, function(x) buffer(SpatialPoints(as.data.frame(rbind(x))), width = 500000))
scmasksMH <- lapply(MH.pts.buf, function(x) mask(MHsumRaw, x))
csmasksMH <- lapply(MH.pts.buf, function(x) mask(MHcanRaw, x))
MHProbs <- list()
for (i in 1:length(scmasksMH)){
  MHProbs[[i]] <- mean(values(scmasksMH[[i]] * csmasksMH[[i]]), na.rm=T)
}
write.csv(unlist(MHProbs), "500kmProbabilitiesMH.csv")
# LGM
LGM.pts <- randomPoints((LGMcanRaw + LGMsumRaw), 1000)
LGM.pts.buf <- apply(LGM.pts, 1, function(x) buffer(SpatialPoints(as.data.frame(rbind(x))), width = 500000))
scmasksLGM <- lapply(LGM.pts.buf, function(x) mask(LGMsumRaw, x))
csmasksLGM <- lapply(LGM.pts.buf, function(x) mask(LGMcanRaw, x))
LGMProbs <- list()
for (i in 1:length(scmasksLGM)){
  LGMProbs[[i]] <- mean(values(scmasksLGM[[i]] * csmasksLGM[[i]]), na.rm=T)
}
write.csv(unlist(LGMProbs), "500kmProbabilitiesLGM.csv")
# LIG
LIG.pts <- randomPoints((LIGcanRaw + LIGsumRaw), 1000)
LIG.pts.buf <- apply(LIG.pts, 1, function(x) buffer(SpatialPoints(as.data.frame(rbind(x))), width = 500000))
scmasksLIG <- lapply(LIG.pts.buf, function(x) mask(LIGsumRaw, x))
csmasksLIG <- lapply(LIG.pts.buf, function(x) mask(LIGcanRaw, x))
LIGProbs <- list()
for (i in 1:length(scmasksLIG)){
  LIGProbs[[i]] <- mean(values(scmasksLIG[[i]] * csmasksLIG[[i]]), na.rm=T)
}
write.csv(unlist(LIGProbs), "500kmProbabilitiesLIG.csv")


