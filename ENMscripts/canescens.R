## Author: Peter J. Galante
## This loads occurrence records, creates a background region, sets up, and performs model tuning and selection.
## Then models are projected to time periods and thresholded using the ETSS. Finally they are projected in raw output for use in "timeSlicePlot.R"

library(spThin);library(ENMeval);library(rgdal);library(rgeos)
bio<-stack(list.files(path = 'layers/', pattern = '\\.tif$', full.names = T))
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
s1<-read.csv('/canescensThin1/thinned_data_thin1.csv')[,2:3]
s2<-read.csv('/canescensThin1/thinned_data_thin2.csv')[,2:3]
s3<-read.csv('/canescensThin1/thinned_data_thin3.csv')[,2:3]
s4<-read.csv('/canescensThin1/thinned_data_thin4.csv')[,2:3]
s5<-read.csv('/canescensThin1/thinned_data_thin5.csv')[,2:3]
length(is.na(extract(bio[[1]], s1))[is.na(extract(bio[[1]], s1))==T]) # 0
length(is.na(extract(bio[[1]], s2))[is.na(extract(bio[[1]], s2))==T]) # 0
length(is.na(extract(bio[[1]], s3))[is.na(extract(bio[[1]], s3))==T]) # 0
length(is.na(extract(bio[[1]], s4))[is.na(extract(bio[[1]], s4))==T]) # 0
length(is.na(extract(bio[[1]], s5))[is.na(extract(bio[[1]], s5))==T]) # 0
## 3 and 4 have the fewest NA points
# load points
randomThin <- sample(1:5,1) # 4
canescens<-read.csv('/canescensThin1/thinned_data_thin1.csv')[,2:3]
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
ca.MCP<-mcp(canescens)
crs(ca.MCP)<-"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
ca.MCP<-spTransform(ca.MCP, CRSobj = crs('+proj=lcc +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs'))
ca.shp <- gBuffer(ca.MCP, width = buff.dist)
ca.shp<-spTransform(ca.shp, CRSobj = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
ca.env.crp<-crop(bio, ca.shp)
ca.env<-mask(ca.env.crp, ca.shp)
plot(ca.env[[1]])
ca.backg<-randomPoints(ca.env, 10000)
write.csv(ca.backg, 'canescens/canescensBackg.csv', row.names = F)
caback<-read.csv('canescens/canescensBackg.csv')
## Tune in ENMeval
cares<-ENMevaluate(occ = canescens, env = ca.env, bg.coords = caback, RMvalues = seq(1,6,0.5), fc = c("L", "LQ", "H", "LQH", "LQHP"),
                   method = 'block', rasterPreds = F, parallel = T, algorithm = "maxent.jar", numCores = 4)

write.csv(cares@results, "canescens/CanescensResults.csv", row.names=F)
Can<-read.csv('canescens/CanescensResults.csv')
CanOpt<-optimize(Can) # L_3.5
dir.create("sumichrasti/SumichrastiModel")
CanMod<-maxent(x = ca.env,
               p = canescens,
               a = caback,
               path = "canescens/CanescensModel",
               args=c(
                 'betamultiplier=3.5',
                 'linear=true',
                 'quadratic=false',
                 'product=false',
                 'threshold=false',
                 'hinge=false',
                 'threads=2',
                 'responsecurves=true',
                 'jackknife=true',
                 'askoverwrite=false'
               )
)
saveRDS(CanMod, "canescens/CanescensModel/CanescensMaxent.RDS")
CanMod<-readRDS("canescens/CanescensModel/CanescensMaxent.RDS")
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
CanPred<-predict(
  object = CanMod,
  x = cur,
  filename = "canescens/CanescensModel/CurrentCan",
  na.rm=T,
  format = "GTiff",
  overwrite=T,
  args = 'cloglog'
)
## Projecting: MH
CanMH<-predict(
  object = CanMod,
  x = MH,
  filename = "canescens/CanescensModel/MHCan",
  na.rm=T,
  format = "GTiff",
  overwrite=F,
  args = 'cloglog'
)
## Projecting: LGM
CanLGM<-predict(
  object = CanMod,
  x = LGM,
  filename = "canescens/CanescensModel/LGMCan",
  na.rm=T,
  format = "GTiff",
  overwrite=F,
  args = 'cloglog'
)
## Projecting: LIG
CanLIG<-predict(
  object = CanMod,
  x = LIG,
  filename = "canescens/CanescensModel/LIGCan",
  na.rm=T,
  format = "GTiff",
  overwrite=F,
  args = 'cloglog'
)

## Predict current to full extent
bio<-stack(list.files(path = 'layers', pattern = '\\.tif$', full.names = T))
colnames(canescens) <- c("x","y")
e1 <- c(extent(canescens)[1]-2, extent(canescens)[2]+2, extent(canescens)[3]-2, extent(canescens)[4]+2)
bio <- crop(bio, e1)

CanCurFull<-predict(
  object = CanMod,
  x = bio,
  filename = "canescens/CanescensModel/CurrentCanFull.tif",
  na.rm=T,
  format = "GTiff",
  overwrite=F,
  args = 'cloglog'
)

## Thresholding models
#threshold
threshor10 <- subset(CanMod@results, rownames(CanMod@results) %in% "X10.percentile.training.presence.Cloglog.threshold")[[1]]
threshETSS <- subset(CanMod@results, rownames(CanMod@results) %in% "Equal.training.sensitivity.and.specificity.Cloglog.threshold")[[1]]
## load rasters
curCanFull <- raster('canescens/CanescensModel/CurrentCanFull.tif')
curCan <- raster('canescens/CanescensModel/CurrentCan.tif')
LGMCan <- raster('canescens/CanescensModel/LGMCan.tif')
LIGCan <- raster('canescens/CanescensModel/LIGCan.tif')
MHCan <- raster('canescens/CanescensModel/MHCan.tif')

curCanFull[curCanFull >= threshETSS] <- 1
curCanFull[curCanFull < threshETSS] <- NA
writeRaster(curCanFull, 'canescens/ETSS/CurCanFullThresh.tif')
curCan[curCan >= threshETSS] <- 1
curCan[curCan < threshETSS] <- NA
writeRaster(curCan, 'canescens/ETSS/curCanThresh.tif')
LGMCan[LGMCan >= threshETSS] <- 1
LGMCan[LGMCan < threshETSS] <- NA
writeRaster(LGMCan, 'canescens/ETSS/LGMCanThresh.tif')
LIGCan[LIGCan >= threshETSS] <- 1
LIGCan[LIGCan < threshETSS] <- NA
writeRaster(LIGCan, 'canescens/ETSS/LIGCanThresh.tif')
MHCan[MHCan >= threshETSS] <- 1
MHCan[MHCan < threshETSS] <- NA
writeRaster(MHCan, 'canescens/ETSS/MHCanThresh.tif')

############################################################################################################
############################################################################################################
#############################  Projecting models in Raw output format  #####################################
CanMod <- readRDS("canescens/CanescensMaxent.RDS")
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
CanPred<-predict(
  object = CanMod,
  x = cur,
  filename = "canescens/CanescensModel/CurrentCanRaw",
  na.rm=T,
  format = "GTiff",
  overwrite=T,
  args = 'raw'
)
## Projecting: MH
CanMH<-predict(
  object = CanMod,
  x = MH,
  filename = "canescens/CanescensModel/MHCanRaw",
  na.rm=T,
  format = "GTiff",
  overwrite=F,
  args = 'raw'
)
## Projecting: LGM
CanLGM<-predict(
  object = CanMod,
  x = LGM,
  filename = "canescens/CanescensModel/LGMCanRaw",
  na.rm=T,
  format = "GTiff",
  overwrite=F,
  args = 'raw'
)
## Projecting: LIG
CanLIG<-predict(
  object = CanMod,
  x = LIG,
  filename = "canescens/CanescensModel/LIGCanRaw",
  na.rm=T,
  format = "GTiff",
  overwrite=F,
  args = 'raw'
)



