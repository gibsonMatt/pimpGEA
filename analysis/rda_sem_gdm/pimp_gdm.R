setwd("/Users/matthew/Box\ Sync/Projects/pimpGEA/associations/dataanalysis/stacks_analysis/rda")

library(gdm)
library(psych)
library(vegan)
library(adegenet)
library(caret)
library(rgdal)
library(lavaan)
library(simba)
library(ggpubr)
filter_7_data_pimp = read.PLINK("/Users/matthew/Box\ Sync/Projects/pimpGEA/associations/data/plinks/stacks/v2/populations.snps.filter7.pimp.recode.raw")
env_info <- read.csv("/Users/matthew/Box\ Sync/Projects/pimpGEA/associations/data/environmental_data/pimp_all_env_variables.csv")
rem <- c("365911-101", "1601-101", "1416-101",  "1987-101", "2649-101", "2933-101", "390695-101", "503517-101", "2857-101", "3123-101")
#rem <- c("365911-101", "1601-101", "1416-101",  "1987-101", "2649-101", "2933-101", "390695-101", "503517-101")
a <- env_info$acc %in% rem
env_info <- env_info[!a,]
filter_7_data_pimp <- filter_7_data_pimp[indNames(filter_7_data_pimp) != "3123-101"]
filter_7_data_pimp <- filter_7_data_pimp[indNames(filter_7_data_pimp) != "2857-101"]

filter_7_data_pimp_imp <- tab(filter_7_data_pimp, NA.method = 'mean')
filter_7_data_pimp_df <- as.data.frame(filter_7_data_pimp_imp)
rem_col <- c('max_aet', 'min_aet', 'cv_aet')
b <- names(env_info) %in% rem_col

env_info <- env_info[!b]
rnames <- env_info$acc

env_info <- data.frame(sapply(env_info, as.numeric))

rownames(env_info) <- rnames


full_env  <- env_info[,c(6:ncol(env_info))]
full_env$txt <- as.factor(full_env$txt)
one_hot <- onehot::onehot(full_env)

full_env <- data.frame(predict(one_hot, full_env))
rownames(full_env) <- rownames(env_info)

full_env$bio_1 <- log(full_env$bio_1)
full_env$bio_2 <- log(full_env$bio_2)
full_env$bio_3 <- log(full_env$bio_3)
#bio4
#bio5
#bio6
#bio7
#bio8
#bio9
#bio10
#bio11
full_env$bio_12 <- log(full_env$bio_12)
full_env$bio_13 <- log(full_env$bio_13)
full_env$bio_14 <- log(full_env$bio_14 + 1)
#bio15
full_env$bio_16 <- log(full_env$bio_16)
full_env$bio_17 <- log(full_env$bio_17 +0.2)
full_env$bio_18 <- log(full_env$bio_18)
full_env$bio_19 <- log(full_env$bio_19 + 0.1)
full_env$ann_avg_rad <- log(full_env$ann_avg_rad)
full_env$ann_max_rad <- log(full_env$ann_max_rad)
#ann_min_rad
full_env$cv_rad <- log(full_env$cv_rad)
#ann_avg_vapr
#ann_max_vapr
full_env$ann_min_vapr <- log(full_env$ann_min_vapr)
#cv_vapr
#ann_avg_wind
#ann_max_wind
#ann_min_wind
#cv_wind
full_env$AHM <- log(full_env$AHM)
#DD5
full_env$DD_18 <- log(full_env$DD_18)
#DD18
#NFFD
#EMT
#Eref
#CMD
full_env$ann_ai <- log(full_env$ann_ai)
#ann_pet
#ann_pet <- log(full_env$ann_pet + 0.01)
#max_pet
#cv_pet
full_env$ann_aet <- log(full_env$ann_aet)
full_env$ann_swt <- log(full_env$ann_swt)
#max_swt
full_env$min_swt <- log(full_env$min_swt)
#cv_swt
full_env$ptac <- log(full_env$ptac + 0.1)
#bdfe
full_env$cec <- log(full_env$cec)
full_env$swc <- log(full_env$swc)
#ph


unscaled_env <- full_env
full_env[,c(1:54)] <- scale(full_env[,c(1:54)])



for(i in 1:ncol(full_env)){ #replace the few missing values with mean for column
  full_env[is.na(full_env[,i]), i] <- mean(full_env[,i], na.rm = TRUE)
}


for(i in 1:ncol(unscaled_env)){ #replace the few missing values with mean for column
  unscaled_env[is.na(unscaled_env[,i]), i] <- mean(unscaled_env[,i], na.rm = TRUE)
}

countNA <- function(row){
  return(sum(is.na(row))/44064)
}

filter_7 <- as.matrix(filter_7_data_pimp)
nas <- apply(filter_7, 2, countNA)
filter_7 <- as.data.frame(filter_7)
filter_7 <- filter_7[names(nas[nas == 0])]


##Env variables from forward selection
#[1] "DD18"         "bio_7"        "ann_max_rad"  "max_swt"      "cv_wind"      "bio_4"        "ann_avg_wind" "ann_max_wind"
#[9] "max_pet"      "CMD"          "cv_vapr"      "EMT"          "bio_3"        "bio_2"   
#Will need to manually select these when the forward selection variable is no longer saved to the environment
good_env <- full_env[filter_7_env_fwd$variables]
red_env <- good_env
#####Variance partitioning#####
###############################

#Projections##########################
latlon <- SpatialPoints(cbind(env_info$lon, -env_info$lat), proj4string=CRS("+proj=longlat"))

pimp_xy <- as.data.frame(spTransform(latlon, CRS("+init=epsg:31973")))
env_info$x <- pimp_xy$coords.x1
env_info$y <- pimp_xy$coords.x2




#####GDM#####
#Format data#
sitespecies <- data.frame(env_info$lon, env_info$lat, filter_7)
sitespecies$site <- row.names(filter_7)
red_env$site <- row.names(filter_7)
red_env  <- data.frame(red_env, env_info$lon, env_info$lat)

no_env <- red_env[,c(15,16,17)]

gdmTab <- formatsitepair(sitespecies, bioFormat = 1, XColumn = 'env_info.lon', YColumn = 'env_info.lat', predData = red_env, siteColumn = "site")
gdmTab_g <- formatsitepair(sitespecies, bioFormat = 1, XColumn = 'env_info.lon', YColumn = 'env_info.lat', predData = no_env, siteColumn = "site")


#Run GDM#
gdm.e <- gdm(gdmTab, geo=F)
gdm.eg <- gdm(gdmTab, geo=T)
gdm.g <- gdm(gdmTab_g, geo=T)


plot(gdm.eg, plot.layout=c(4,3))




gdm.1.imp <- gdm.varImp(gdmTab, geo=T, nPerm=50, cores=12, parallel = T, fullModelOnly = F)

gdm.1.imp.full <- gdm.varImp(gdmTab, geo=T, nPerm=999, cores=12, parallel = T, fullModelOnly = T)
gdm.1.imp.full_nogeo <- gdm.varImp(gdmTab, geo=F, nPerm=999, cores=12, parallel = T, fullModelOnly = T)
