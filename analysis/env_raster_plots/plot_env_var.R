setwd("/Users/matthew/Box\ Sync/Projects/pimpGEA/associations/dataanalysis/envplots")

library(raster)
library(adegenet)
library(maps)
library(scales)
library(spatialEco)
pimp_info <- read.csv("./pimp_info.csv")



ann_pet <- raster("/Volumes/thevoid/Projects/pimpGEA/PET/PET_he_annual/pet_he_yr/w001001.adf")

pal <- colorRampPalette(c("white","black"))

peru <- getData("GADM",country="Peru",level=0)
ecu <- getData("GADM",country="Ecuador",level=0)
bra <- getData("GADM",country="Brazil",level=0)
col <- getData("GADM",country="Colombia",level=0)

#png("pet.png", width = 7, height = 7, units="in", res=400)
plot(ann_pet, xlim=c(-85,-68), ylim=c(-17,2),  col = bluepal(1000))

plot(peru, add=T)
plot(ecu, add=T)
#plot(bra, add=T, col="white")
#plot(col, add=T, col="white")

points(pimp_info$lon, pimp_info$lat, col="black", pch=21, bg=funky(3)[2])
#dev.off()





ann_prec <- raster("/Volumes/thevoid/Projects/pimpGEA/worldclim/wc2.0_30s_bio/wc2.0_bio_30s_12.tif")
ann_prec <- raster.transformation(ann_prec, trans = "nl", smin = 0, smax = 255)


png("prec.png", width = 7, height = 7, units="in", res=400)
plot(ann_prec, xlim=c(-85,-68), ylim=c(-17,2),  col = greenpal(1000))

plot(peru, add=T)
plot(ecu, add=T)
#plot(bra, add=T, col="white")
#plot(col, add=T, col="white")

points(pimp_info$lon, pimp_info$lat, col="black", pch=21, bg=funky(3)[2])
map.scale()
dev.off()





wind_1 <- raster("/Volumes/thevoid/Projects/pimpGEA/worldclim/wc2.0_30s_wind/wc2.0_30s_wind_01.tif")
wind_2 <- raster("/Volumes/thevoid/Projects/pimpGEA/worldclim/wc2.0_30s_wind/wc2.0_30s_wind_02.tif")
wind_3 <- raster("/Volumes/thevoid/Projects/pimpGEA/worldclim/wc2.0_30s_wind/wc2.0_30s_wind_03.tif")
wind_4 <- raster("/Volumes/thevoid/Projects/pimpGEA/worldclim/wc2.0_30s_wind/wc2.0_30s_wind_04.tif")
wind_5 <- raster("/Volumes/thevoid/Projects/pimpGEA/worldclim/wc2.0_30s_wind/wc2.0_30s_wind_05.tif")
wind_6 <- raster("/Volumes/thevoid/Projects/pimpGEA/worldclim/wc2.0_30s_wind/wc2.0_30s_wind_06.tif")
wind_7 <- raster("/Volumes/thevoid/Projects/pimpGEA/worldclim/wc2.0_30s_wind/wc2.0_30s_wind_07.tif")
wind_8 <- raster("/Volumes/thevoid/Projects/pimpGEA/worldclim/wc2.0_30s_wind/wc2.0_30s_wind_08.tif")
wind_9 <- raster("/Volumes/thevoid/Projects/pimpGEA/worldclim/wc2.0_30s_wind/wc2.0_30s_wind_09.tif")
wind_10 <- raster("/Volumes/thevoid/Projects/pimpGEA/worldclim/wc2.0_30s_wind/wc2.0_30s_wind_10.tif")
wind_11 <- raster("/Volumes/thevoid/Projects/pimpGEA/worldclim/wc2.0_30s_wind/wc2.0_30s_wind_11.tif")
wind_12 <- raster("/Volumes/thevoid/Projects/pimpGEA/worldclim/wc2.0_30s_wind/wc2.0_30s_wind_12.tif")
rr <- list(wind_1, wind_2, wind_3, wind_4, wind_5, wind_6, wind_7, wind_8, wind_9, wind_10, wind_11,
           wind_12)
stack  <-  stack(rr)

#ann_max_wnd <- raster.transformation(ann_prec, trans = "nl", smin = 0, smax = 255)
maxs <- max(stack)

png("wind.png", width = 7, height = 7, units="in", res=400)
plot(maxs, xlim=c(-85,-68), ylim=c(-17,2),  col = bluepal(1000))

plot(peru, add=T)
plot(ecu, add=T)
#plot(bra, add=T, col="white")
#plot(col, add=T, col="white")

points(pimp_info$lon, pimp_info$lat, col="black", pch=21, bg=funky(3)[2])
map.scale()
dev.off()





peru <- getData("GADM",country="Peru",level=0)
ecu <- getData("GADM",country="Ecuador",level=0)
bra <- getData("GADM",country="Brazil",level=0)
col <- getData("GADM",country="Colombia",level=0)


vapr_1 <- raster("/Volumes/thevoid/Projects/pimpGEA/worldclim/wc2.0_30s_vapr/wc2.0_30s_vapr_01.tif")
vapr_2 <- raster("/Volumes/thevoid/Projects/pimpGEA/worldclim/wc2.0_30s_vapr/wc2.0_30s_vapr_02.tif")
vapr_3 <- raster("/Volumes/thevoid/Projects/pimpGEA/worldclim/wc2.0_30s_vapr/wc2.0_30s_vapr_03.tif")
vapr_4 <- raster("/Volumes/thevoid/Projects/pimpGEA/worldclim/wc2.0_30s_vapr/wc2.0_30s_vapr_04.tif")
vapr_5 <- raster("/Volumes/thevoid/Projects/pimpGEA/worldclim/wc2.0_30s_vapr/wc2.0_30s_vapr_05.tif")
vapr_6 <- raster("/Volumes/thevoid/Projects/pimpGEA/worldclim/wc2.0_30s_vapr/wc2.0_30s_vapr_06.tif")
vapr_7 <- raster("/Volumes/thevoid/Projects/pimpGEA/worldclim/wc2.0_30s_vapr/wc2.0_30s_vapr_07.tif")
vapr_8 <- raster("/Volumes/thevoid/Projects/pimpGEA/worldclim/wc2.0_30s_vapr/wc2.0_30s_vapr_08.tif")
vapr_9 <- raster("/Volumes/thevoid/Projects/pimpGEA/worldclim/wc2.0_30s_vapr/wc2.0_30s_vapr_09.tif")
vapr_10 <- raster("/Volumes/thevoid/Projects/pimpGEA/worldclim/wc2.0_30s_vapr/wc2.0_30s_vapr_10.tif")
vapr_11 <- raster("/Volumes/thevoid/Projects/pimpGEA/worldclim/wc2.0_30s_vapr/wc2.0_30s_vapr_11.tif")
vapr_12 <- raster("/Volumes/thevoid/Projects/pimpGEA/worldclim/wc2.0_30s_vapr/wc2.0_30s_vapr_12.tif")
rr <- list(vapr_1, vapr_2, vapr_3, vapr_4, vapr_5, vapr_6, vapr_7, vapr_8, vapr_9, vapr_10, vapr_11,
           vapr_12)
stack  <-  stack(rr)
mean_vapr <- mean(stack)

png("vapr.png", width = 7, height = 7, units="in", res=400)
plot(mean_vapr, xlim=c(-85,-68), ylim=c(-17,2),  col = greenpal(1000))

plot(peru, add=T)
plot(ecu, add=T)
#plot(bra, add=T, col="white")
#plot(col, add=T, col="white")

points(pimp_info$lon, pimp_info$lat, col="black", pch=21, bg=funky(3)[2])
map.scale(ratio=F)
dev.off()





seasonality <- raster("/Volumes/thevoid/Projects/pimpGEA/worldclim/wc2.0_30s_bio/wc2.0_bio_30s_15.tif")

png("season.png", width = 7, height = 7, units="in", res=400)
plot(seasonality, xlim=c(-85,-68), ylim=c(-17,2),  col = bluepal(1000))

plot(peru, add=T)
plot(ecu, add=T)
#plot(bra, add=T, col="white")
#plot(col, add=T, col="white")

points(pimp_info$lon, pimp_info$lat, col="black", pch=21, bg=funky(3)[2])
map.scale(ratio=F)
dev.off()