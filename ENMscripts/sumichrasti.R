## Author: Peter J. Galante
## This loads occurrence records, creates a background region, sets up, and performs model tuning and selection.
## Then models are projected to time periods and thresholded using the ETSS. Finally they are projected in raw output for use in "timeSlicePlot.R"

library(spThin);library(ENMeval);library(rgdal);library(rgeos)
bio<-stack(list.files(path = 'layers', pattern = '\\.tif$', full.names = T))
##optimize function
optimize <- function(res) {
  ###Remove any candidate model which has an AUC less than 0.51= models with no discrimination
  opt.auc <- res[res$train.AUC >= 0.5,]
  ###Remove any candidates which have no parameters
  no.param <- opt.auc[opt.auc$parameters > 1,]
  ###Remove any candidates where the AIC score was NA (too many parameters)
  noAICNA<- no.param[!is.na(no.param$parameters),]
  #noAICNA<- no.param[which(!is.na(no.param$AICc)),]
  ###Remove any models which have an OR of zero
  noOR0 <- noAICNA[noAICNA$avg.test.or10pct != 0,]
  ###Order the remaining list by lowest OR then highest AUC, sequentially
  ordered<-noOR0[with(noOR0, order(avg.test.or10pct, -avg.test.AUC)), ]
  ###Grab the settings of that first model (the optimal model)
  ordered[1,]
}
######################################################################
#####################  sumichrasti  ##################################
######################################################################
## sumichrasti - thin points and select dataset
s1<-read.csv('/sumichrastiThin1/thinned_data_thin1.csv')[,2:3]
s2<-read.csv('/sumichrastiThin1/thinned_data_thin2.csv')[,2:3]
s3<-read.csv('/sumichrastiThin1/thinned_data_thin3.csv')[,2:3]
s4<-read.csv('/sumichrastiThin1/thinned_data_thin4.csv')[,2:3]
s5<-read.csv('/sumichrastiThin1/thinned_data_thin5.csv')[,2:3]
length(is.na(extract(bio[[1]], s1))[is.na(extract(bio[[1]], s1))==T]) # 7
length(is.na(extract(bio[[1]], s2))[is.na(extract(bio[[1]], s2))==T]) # 6
length(is.na(extract(bio[[1]], s3))[is.na(extract(bio[[1]], s3))==T]) # 5
length(is.na(extract(bio[[1]], s4))[is.na(extract(bio[[1]], s4))==T]) # 5
length(is.na(extract(bio[[1]], s5))[is.na(extract(bio[[1]], s5))==T]) # 7
## 3 and 4 have the fewest NA points
# load points
randomThin <- sample(3:4,1) # 4
sumichrasti<-read.csv('/sumichrastiThin1/thinned_data_thin4.csv')[,2:3]
# Generate background data
mcp <- function (xy) {
  xy <- as.data.frame(coordinates(xy))
  coords.t <- chull(xy[, 1], xy[, 2])
  xy.bord <- xy[coords.t, ]
  xy.bord <- rbind(xy.bord[nrow(xy.bord), ], xy.bord)
  return(SpatialPolygons(list(Polygons(list(Polygon(as.matrix(xy.bord))), 1))))
}
##  Set a 10km buffer distance
buff.dist <- 10000
## Creating A MCP background for thinned canescens data
sc.MCP<-mcp(sumichrasti)
crs(sc.MCP)<-"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
sc.MCP<-spTransform(sc.MCP, CRSobj = crs('+proj=lcc +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs'))
sc.shp <- gBuffer(sc.MCP, width = buff.dist)
sc.shp<-spTransform(sc.shp, CRSobj = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
sc.env.crp<-crop(bio, sc.shp)
sc.env<-mask(sc.env.crp, sc.shp)
plot(sc.env[[1]])
sc.backg<-randomPoints(sc.env, 10000)
write.csv(sc.backg, 'sumichrasti/sumichrastiBackg.csv', row.names = F)
scback<-read.csv('sumichrasti/sumichrastiBackg.csv')
## Tune in ENMeval
scres<-ENMevaluate(occ = sumichrasti, env = sc.env, bg.coords = scback, RMvalues = seq(5,10,0.5), fc = c("L", "LQ", "H", "LQH", "LQHP"),
                   method = 'block', rasterPreds = F, parallel = T, algorithm = "maxent.jar", numCores = 4)

write.csv(scres@results, "sumichrasti/SumichrastiResults.csv", row.names=F)
Sum<-read.csv('sumichrasti/SumichrastiResults.csv')
SumOpt<-optimize(Sum) # LQ_8.5
dir.create("sumichrasti/SumichrastiModel")
SumMod<-maxent(x = sc.env,
               p = sumichrasti,
               a = scback,
               path = "sumichrasti/SumichrastiModel",
               args=c(
                 'betamultiplier=8.5',
                 'linear=true',
                 'quadratic=true',
                 'product=false',
                 'threshold=false',
                 'hinge=false',
                 'threads=2',
                 'responsecurves=true',
                 'jackknife=true',
                 'askoverwrite=false'
               )
)
saveRDS(SumMod, "sumichrasti/SumichrastiModel/SumichrastiMaxent.RDS")
SumMod<-readRDS("sumichrasti/SumichrastiModel/SumichrastiMaxent.RDS")
# Read env data dn crop to appropriate area
e <- extent(c(-82.6, -77.5, 6, 11))
# Current
cur <- crop(bio, e)
# Mid-holocene
MH <- stack(list.files(path = "layers/mid_holo0-5", pattern = "\\.tif$", full.names = T))
MH <- crop(MH, e)
names(MH) <- names(bio)
# LGM
LGM <- stack(list.files(path = "layers/LGM0-5/", pattern = "\\.tif$", full.names = T))
LGM <- crop(LGM, e)
names(LGM) <- names(bio)
# LIG
LIG <- stack(list.files(path = "layers/LIG0-5/", pattern = "\\.tif$", full.names = T))
LIG <- crop(LIG, e)
names(LIG) <- names(bio)

## Projecting: Current
SumPred<-predict(
  object = SumMod,
  x = cur,
  filename = "sumichrasti/SumichrastiModel/CurrentSum",
  na.rm=T,
  format = "GTiff",
  overwrite=T,
  args = 'cloglog'
)
## Projecting: MH
SumMH<-predict(
  object = SumMod,
  x = MH,
  filename = "sumichrasti/SumichrastiModel/MHSum",
  na.rm=T,
  format = "GTiff",
  overwrite=F,
  args = 'cloglog'
)
## Projecting: LGM
SumLGM<-predict(
  object = SumMod,
  x = LGM,
  filename = "sumichrasti/SumichrastiModel/LGMSum",
  na.rm=T,
  format = "GTiff",
  overwrite=F,
  args = 'cloglog'
)
## Projecting: LIG
SumLIG<-predict(
  object = SumMod,
  x = LIG,
  filename = "sumichrasti/SumichrastiModel/LIGSum",
  na.rm=T,
  format = "GTiff",
  overwrite=F,
  args = 'cloglog'
)

## Predict current to full extent
bio<-stack(list.files(path = 'layers', pattern = '\\.tif$', full.names = T))
colnames(sumichrasti) <- c("x","y")
e1 <- c(extent(sumichrasti)[1]-2, extent(sumichrasti)[2]+2, extent(sumichrasti)[3]-2, extent(sumichrasti)[4]+2)
bio <- crop(bio, e1)

SumCurFull<-predict(
  object = SumMod,
  x = bio,
  filename = "sumichrasti/SumichrastiModel/CurrentSumFull.tif",
  na.rm=T,
  format = "GTiff",
  overwrite=F,
  args = 'cloglog'
)

## Thresholding models
#threshold
threshor10 <- subset(SumMod@results, rownames(SumMod@results) %in% "X10.percentile.training.presence.Cloglog.threshold")[[1]]
threshETSS <- subset(SumMod@results, rownames(SumMod@results) %in% "Equal.training.sensitivity.and.specificity.Cloglog.threshold")[[1]]
curSumFull <- raster('sumichrasti/SumichrastiModel/CurrentSumFull.tif')
curSum <- raster('sumichrasti/SumichrastiModel/CurrentSum.tif')
LGMSum <- raster('sumichrasti/SumichrastiModel/LGMSum.tif')
LIGSum <- raster('sumichrasti/SumichrastiModel/LIGSum.tif')
MHSum <- raster('sumichrasti/SumichrastiModel/MHSum.tif')

curSumFull[curSumFull >= threshETSS] <- 1
curSumFull[curSumFull < threshETSS] <- NA
writeRaster(curSumFull, 'sumichrasti/ETSS/CurSumFullThresh.tif')
curSum[curSum >= threshETSS] <- 1
curSum[curSum < threshETSS] <- NA
writeRaster(curSum, 'sumichrasti/ETSS/curSumThresh.tif')
LGMSum[LGMSum >= threshETSS] <- 1
LGMSum[LGMSum < threshETSS] <- NA
writeRaster(LGMSum, 'sumichrasti/ETSS/LGMSumThresh.tif')
LIGSum[LIGSum >= threshETSS] <- 1
LIGSum[LIGSum < threshETSS] <- NA
writeRaster(LIGSum, 'sumichrasti/ETSS/LIGSumThresh.tif')
MHSum[MHSum >= threshETSS] <- 1
MHSum[MHSum < threshETSS] <- NA
writeRaster(MHSum, 'sumichrasti/ETSS/MHSUMThresh.tif')


############################################################################################################
############################################################################################################
#############################  Projecting models in Raw output format  #####################################
saveRDS(SumMod, "sumichrasti/SumichrastiModel/SumichrastiMaxent.RDS")
SumMod<-readRDS("sumichrasti/SumichrastiModel/SumichrastiMaxent.RDS")
bio<-stack(list.files(path = 'layers', pattern = '\\.tif$', full.names = T))
# Read env data dn crop to appropriate area
e <- extent(c(-82.6, -77.5, 6, 11))
# Current
cur <- crop(bio, e)
# Mid-holocene
MH <- stack(list.files(path = "layers/mid_holo0-5", pattern = "\\.tif$", full.names = T))
MH <- crop(MH, e)
names(MH) <- names(bio)
# LGM
LGM <- stack(list.files(path = "layers/LGM0-5/", pattern = "\\.tif$", full.names = T))
LGM <- crop(LGM, e)
names(LGM) <- names(bio)
# LIG
LIG <- stack(list.files(path = "layers/LIG0-5/", pattern = "\\.tif$", full.names = T))
LIG <- crop(LIG, e)
names(LIG) <- names(bio)

## Projecting: Current
SumPred<-predict(
  object = SumMod,
  x = cur,
  filename = "sumichrasti/SumichrastiModel/CurrentSumRaw",
  na.rm=T,
  format = "GTiff",
  overwrite=T,
  args = 'raw'
)
## Projecting: MH
SumMH<-predict(
  object = SumMod,
  x = MH,
  filename = "sumichrasti/SumichrastiModel/MHSumRaw",
  na.rm=T,
  format = "GTiff",
  overwrite=F,
  args = 'raw'
)
## Projecting: LGM
SumLGM<-predict(
  object = SumMod,
  x = LGM,
  filename = "sumichrasti/SumichrastiModel/LGMSumRaw",
  na.rm=T,
  format = "GTiff",
  overwrite=F,
  args = 'raw'
)
## Projecting: LIG
SumLIG<-predict(
  object = SumMod,
  x = LIG,
  filename = "sumichrasti/SumichrastiModel/LIGSumRaw",
  na.rm=T,
  format = "GTiff",
  overwrite=F,
  args = 'raw'
)
