setwd("/Users/matthew/Box\ Sync/Projects/pimpGEA/associations/dataanalysis/stacks_analysis/rda")
library(psych)
library(vegan)
library(adegenet)
library(caret)
library(rgdal)
library(lavaan)
library(simba)
library(ggpubr)
library(adespatial)


load_funcs <- function(){
  # CA.newr
  source ('https://raw.githubusercontent.com/zdealveindy/anadat-r/master/scripts/NumEcolR2/CA.newr.R')
  
  # PCA.newr
  source ('https://raw.githubusercontent.com/zdealveindy/anadat-r/master/scripts/NumEcolR2/PCA.newr.R')
  
  # Rao
  source ('https://raw.githubusercontent.com/zdealveindy/anadat-r/master/scripts/NumEcolR2/Rao.R')
  
  # bartlett.perm
  source ('https://raw.githubusercontent.com/zdealveindy/anadat-r/master/scripts/NumEcolR2/bartlett.perm.R')
  
  # boxplerk
  source ('https://raw.githubusercontent.com/zdealveindy/anadat-r/master/scripts/NumEcolR2/boxplerk.R')
  
  # boxplert
  source ('https://raw.githubusercontent.com/zdealveindy/anadat-r/master/scripts/NumEcolR2/boxplert.R')
  
  # cleanplot.pca
  source ('https://raw.githubusercontent.com/zdealveindy/anadat-r/master/scripts/NumEcolR2/cleanplot.pca.R')
  
  # coldiss
  source ('https://raw.githubusercontent.com/zdealveindy/anadat-r/master/scripts/NumEcolR2/coldiss.R')
  
  # drawmap
  source ('https://raw.githubusercontent.com/zdealveindy/anadat-r/master/scripts/NumEcolR2/drawmap.R')
  
  # drawmap3
  source ('https://raw.githubusercontent.com/zdealveindy/anadat-r/master/scripts/NumEcolR2/drawmap3.R')
  
  # hcoplot
  source ('https://raw.githubusercontent.com/zdealveindy/anadat-r/master/scripts/NumEcolR2/hcoplot.R')
  
  # panelutils
  source ('https://raw.githubusercontent.com/zdealveindy/anadat-r/master/scripts/NumEcolR2/panelutils.R')
  
  # plot.lda
  source ('https://raw.githubusercontent.com/zdealveindy/anadat-r/master/scripts/NumEcolR2/plot.lda.R')
  
  # plot.links
  source ('https://raw.githubusercontent.com/zdealveindy/anadat-r/master/scripts/NumEcolR2/plot.links.R')
  
  # polyvars
  source ('https://raw.githubusercontent.com/zdealveindy/anadat-r/master/scripts/NumEcolR2/polyvars.R')
  
  # quickMEM
  source ('https://raw.githubusercontent.com/zdealveindy/anadat-r/master/scripts/NumEcolR2/quickMEM.R')
  
  # scalog
  source ('https://raw.githubusercontent.com/zdealveindy/anadat-r/master/scripts/NumEcolR2/scalog.R')
  
  # screestick
  source ('https://raw.githubusercontent.com/zdealveindy/anadat-r/master/scripts/NumEcolR2/screestick.R')
  
  # sr.value
  source ('https://raw.githubusercontent.com/zdealveindy/anadat-r/master/scripts/NumEcolR2/sr.value.R')
  
  # triplot.rda
  source ('https://raw.githubusercontent.com/zdealveindy/anadat-r/master/scripts/NumEcolR2/triplot.rda.R')
  
}
load_funcs()

filter_7_data_pimp = read.PLINK("/Users/matthew/Box\ Sync/Projects/pimpGEA/associations/data/plinks/stacks/v2/populations.snps.filter7.pimp.recode.raw")
#filter_5_data_pimp = read.PLINK("/Users/matthew/Box\ Sync/Projects/pimpGEA/associations/data/plinks/stacks/v2/populations.snps.filter5.pimp.raw")



env_info <- read.csv("/Users/matthew/Box\ Sync/Projects/pimpGEA/associations/data/environmental_data/pimp_all_env_variables.csv")

rem <- c("365911-101", "1601-101", "1416-101",  "1987-101", "2649-101", "2933-101", "390695-101", "503517-101", "2857-101", "3123-101")

a <- env_info$acc %in% rem
env_info <- env_info[!a,]

filter_7_data_pimp <- filter_7_data_pimp[indNames(filter_7_data_pimp) != "3123-101"]
filter_7_data_pimp <- filter_7_data_pimp[indNames(filter_7_data_pimp) != "2857-101"]

#filter_5_data_pimp <- filter_5_data_pimp[indNames(filter_5_data_pimp) != "3123-101"]
#filter_5_data_pimp <- filter_5_data_pimp[indNames(filter_5_data_pimp) != "2857-101"]

filter_7_data_pimp_imp <- tab(filter_7_data_pimp, NA.method = 'mean')
filter_7_data_pimp_df <- as.data.frame(filter_7_data_pimp_imp)


countNA_7 <- function(row){
  return(sum(is.na(row))/44064)
}
countNA_5 <- function(row){
  return(sum(is.na(row))/44243)
}


filter_7 <- as.matrix(filter_7_data_pimp)
nas <- apply(filter_7, 2, countNA_7)
filter_7 <- as.data.frame(filter_7)
filter_7 <- filter_7[names(nas[nas == 0])]


#filter_5 <- as.matrix(filter_5_data_pimp)
#nas <- apply(filter_5, 2, countNA_5)
#filter_5 <- as.data.frame(filter_5)
#filter_5 <- filter_5[names(nas[nas < 0.001])]
#dim(filter_5)

#filter_5_data_pimp_imp <- tab(as.genlight(filter_5), NA.method = 'mean')
#filter_5_data_pimp_df <- as.data.frame(filter_5_data_pimp_imp)



####dbMEMs#####
##############

#Projections##########################
latlon <- SpatialPoints(cbind(env_info$lon, -env_info$lat), proj4string=CRS("+proj=longlat"))

pimp_xy <- as.data.frame(spTransform(latlon, CRS("+init=epsg:31973")))
env_info$x <- pimp_xy$coords.x1
env_info$y <- pimp_xy$coords.x2
#######################################


#Choosing dbMEMs#

#Forward selection of dbMEMs

dbmem <- quickMEM(filter_7, pimp_xy, perm.max=999)

dbmem_all <- dbmem(pimp_xy)

#These were identified with quickMEM()
#dbmem_original11 <- dbmem_all[,c(1,2,3,4,5,6,7,8,9,17)]

dbmem_10 <- dbmem_all[,c(1,2,3,4,5,6,7,8,9,17)]

#png("sig_dbmems.png", width = 10, height = 8, units="in", res=600)
par(mfrow = c(2, 5))
for(i in 1 : ncol(dbmem$dbMEM_red_model))
{
  sr.value(pimp_xy, dbmem$dbMEM_red_model[,i], sub=names(dbmem$dbMEM_red_model)[i])
}
#dev.off()
#1-4 = Broad scale structure
#5, 6, 7, 8, 9, 11, 17 = Fine scale structure



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

full_env <- full_env[,c(-60)]
full_env <- full_env[,c(-59)]

#######Forward selection with env variables########

filter_7_env_rda <- rda(filter_7 ~., full_env)
(filter_7_env_r2a <- RsquareAdj(filter_7_env_rda)$adj.r.squared)

filter_7_env_fwd <- forward.sel(filter_7, full_env,
                                adjR2thresh = filter_7_env_r2a,
                                nperm = 9999, alpha=0.01)
#This selection above identified these variables...
#[1] "DD18"         "bio_7"        "ann_max_rad"  "max_swt"      "cv_wind"      "bio_4"        "ann_avg_wind" "ann_max_wind"
#[9] "max_pet"      "CMD"          "cv_vapr"      "EMT"          "bio_3"        "bio_2"   

variables <- c("DD18", "bio_7", "ann_max_rad", "max_swt", "cv_wind", "bio_4", "ann_avg_wind", "ann_max_wind",
                                "max_pet", "CMD", "cv_vapr" ,"EMT", "bio_3", "bio_2" )
pc_env <- prcomp(full_env, scale. = T)


filter_7_pca_rda <- rda(filter_7 ~., data.frame(pc_env$x))
(filter_7_pca_r2a <- RsquareAdj(filter_7_pca_rda)$adj.r.squared)

filter_7_pca_fwd <- forward.sel(filter_7, pc_env$x,
                                adjR2thresh = filter_7_pca_r2a,
                                nperm = 99)



#####Variance partitioning#####
###############################

good_env <- full_env[variables]

vp <- varpart(filter_7, good_env, dbmem$dbMEM_red_model[,c(1:4)], dbmem$dbMEM_red_model[,c(5:9)])

vp_comb <- varpart(filter_7, good_env, dbmem$dbMEM_red_model)
png("venn_partitions.png", width = 6, height = 6, units="in", res=600)
plot(vp_comb, bg=c("red", "blue"), digits=1, Xnames=c('Env', 'Space'))
dev.off()

pc_env_selected <- pc_env$x[,filter_7_pca_fwd$variables]

vp_comb_pcs <- varpart(filter_7, pc_env_selected, dbmem$dbMEM_red_model)
plot(vp_comb_pcs, bg=c("red", "blue"), digits=1, Xnames=c('Env', 'Space'))

anova(rda(filter_7, full_env[filter_7_env_fwd$variables], dbmem))
      
anova(rda(filter_7, dbmem, full_env[filter_7_env_fwd$variables]))



###SEM
filter_7_gl <- new("genlight", filter_7)
#filter_5_gl <- new("genlight", filter_5)

geo.dist <- dist(pimp_xy)
env.dist.1 <- base::scale(liste(vegdist(full_env[,1], method = "euclidean"), entry="Env1")$Env1)
env.dist.2 <- base::scale(liste(vegdist(full_env[,2], method = "euclidean"), entry="Env2")$Env2)
env.dist.3 <- base::scale(liste(vegdist(full_env[,3], method = "euclidean"), entry="Env3")$Env3)
env.dist.4 <- base::scale(liste(vegdist(full_env[,4], method = "euclidean"), entry="Env4")$Env4)
env.dist.5 <- base::scale(liste(vegdist(full_env[,5], method = "euclidean"), entry="Env5")$Env5)
env.dist.6 <- base::scale(liste(vegdist(full_env[,6], method = "euclidean"), entry="Env6")$Env6)
env.dist.7 <- base::scale(liste(vegdist(full_env[,7], method = "euclidean"), entry="Env7")$Env7)
env.dist.8 <- base::scale(liste(vegdist(full_env[,8], method = "euclidean"), entry="Env8")$Env8)
env.dist.9 <- base::scale(liste(vegdist(full_env[,9], method = "euclidean"), entry="Env9")$Env9)
env.dist.10 <- base::scale(liste(vegdist(full_env[,10], method = "euclidean"), entry="Env10")$Env10)
env.dist.11 <- base::scale(liste(vegdist(full_env[,11], method = "euclidean"), entry="Env11")$Env11)
env.dist.12 <- base::scale(liste(vegdist(full_env[,12], method = "euclidean"), entry="Env12")$Env12)
env.dist.13 <- base::scale(liste(vegdist(full_env[,13], method = "euclidean"), entry="Env13")$Env13)
env.dist.14 <- base::scale(liste(vegdist(full_env[,14], method = "euclidean"), entry="Env14")$Env14)
env.dist.15 <- base::scale(liste(vegdist(full_env[,15], method = "euclidean"), entry="Env15")$Env15)
env.dist.16 <- base::scale(liste(vegdist(full_env[,16], method = "euclidean"), entry="Env16")$Env16)
env.dist.17 <- base::scale(liste(vegdist(full_env[,17], method = "euclidean"), entry="Env17")$Env17)
env.dist.18 <- base::scale(liste(vegdist(full_env[,18], method = "euclidean"), entry="Env18")$Env18)
env.dist.19 <- base::scale(liste(vegdist(full_env[,19], method = "euclidean"), entry="Env19")$Env19)
env.dist.20 <- base::scale(liste(vegdist(full_env[,20], method = "euclidean"), entry="Env20")$Env20)
env.dist.21 <- base::scale(liste(vegdist(full_env[,21], method = "euclidean"), entry="Env21")$Env21)
env.dist.22 <- base::scale(liste(vegdist(full_env[,22], method = "euclidean"), entry="Env22")$Env22)
env.dist.23 <- base::scale(liste(vegdist(full_env[,23], method = "euclidean"), entry="Env23")$Env23)
env.dist.24 <- base::scale(liste(vegdist(full_env[,24], method = "euclidean"), entry="Env24")$Env24)
env.dist.25 <- base::scale(liste(vegdist(full_env[,25], method = "euclidean"), entry="Env25")$Env25)
env.dist.26 <- base::scale(liste(vegdist(full_env[,26], method = "euclidean"), entry="Env26")$Env26)
env.dist.27 <- base::scale(liste(vegdist(full_env[,27], method = "euclidean"), entry="Env27")$Env27)
env.dist.28 <- base::scale(liste(vegdist(full_env[,28], method = "euclidean"), entry="Env28")$Env28)
env.dist.29 <- base::scale(liste(vegdist(full_env[,29], method = "euclidean"), entry="Env29")$Env29)
env.dist.30 <- base::scale(liste(vegdist(full_env[,30], method = "euclidean"), entry="Env30")$Env30)
env.dist.31 <- base::scale(liste(vegdist(full_env[,31], method = "euclidean"), entry="Env31")$Env31)
env.dist.32 <- base::scale(liste(vegdist(full_env[,32], method = "euclidean"), entry="Env32")$Env32)
env.dist.33 <- base::scale(liste(vegdist(full_env[,33], method = "euclidean"), entry="Env33")$Env33)
env.dist.34 <- base::scale(liste(vegdist(full_env[,34], method = "euclidean"), entry="Env34")$Env34)
env.dist.35 <- base::scale(liste(vegdist(full_env[,35], method = "euclidean"), entry="Env35")$Env35)
env.dist.36 <- base::scale(liste(vegdist(full_env[,36], method = "euclidean"), entry="Env36")$Env36)
env.dist.37 <- base::scale(liste(vegdist(full_env[,37], method = "euclidean"), entry="Env37")$Env37)
env.dist.38 <- base::scale(liste(vegdist(full_env[,38], method = "euclidean"), entry="Env38")$Env38)
env.dist.39 <- base::scale(liste(vegdist(full_env[,39], method = "euclidean"), entry="Env39")$Env39)
env.dist.40 <- base::scale(liste(vegdist(full_env[,40], method = "euclidean"), entry="Env40")$Env40)
env.dist.41 <- base::scale(liste(vegdist(full_env[,41], method = "euclidean"), entry="Env41")$Env41)
env.dist.42 <- base::scale(liste(vegdist(full_env[,42], method = "euclidean"), entry="Env42")$Env42)
env.dist.43 <- base::scale(liste(vegdist(full_env[,43], method = "euclidean"), entry="Env43")$Env43)
env.dist.44 <- base::scale(liste(vegdist(full_env[,44], method = "euclidean"), entry="Env44")$Env44)
env.dist.45 <- base::scale(liste(vegdist(full_env[,45], method = "euclidean"), entry="Env45")$Env45)
env.dist.46 <- base::scale(liste(vegdist(full_env[,46], method = "euclidean"), entry="Env46")$Env46)
env.dist.47 <- base::scale(liste(vegdist(full_env[,47], method = "euclidean"), entry="Env47")$Env47)
env.dist.48 <- base::scale(liste(vegdist(full_env[,48], method = "euclidean"), entry="Env48")$Env48)
env.dist.49 <- base::scale(liste(vegdist(full_env[,49], method = "euclidean"), entry="Env49")$Env49)
env.dist.50 <- base::scale(liste(vegdist(full_env[,50], method = "euclidean"), entry="Env50")$Env50)
env.dist.51 <- base::scale(liste(vegdist(full_env[,51], method = "euclidean"), entry="Env51")$Env51)
env.dist.52 <- base::scale(liste(vegdist(full_env[,52], method = "euclidean"), entry="Env52")$Env52)
env.dist.53 <- base::scale(liste(vegdist(full_env[,53], method = "euclidean"), entry="Env53")$Env53)
env.dist.54 <- base::scale(liste(vegdist(full_env[,54], method = "euclidean"), entry="Env54")$Env54)
env.dist.55 <- base::scale(liste(vegdist(full_env[,55], method = "euclidean"), entry="Env55")$Env55)
env.dist.56 <- base::scale(liste(vegdist(full_env[,56], method = "euclidean"), entry="Env56")$Env56)
env.dist.57 <- base::scale(liste(vegdist(full_env[,57], method = "euclidean"), entry="Env57")$Env57)
env.dist.58 <- base::scale(liste(vegdist(full_env[,58], method = "euclidean"), entry="Env58")$Env58)

gen_dist  <- dist(filter_7_gl)
#gen_dist  <- dist(filter_5_gl)


geo.dist <- base::scale(liste(geo.dist, entry="Geo")$Geo)
gen.dist <- base::scale(liste(gen_dist, entry="Gen")$Gen)


full  <- data.frame(gen.dist, env.dist.1, env.dist.2, env.dist.3, env.dist.4, env.dist.5, 
                    env.dist.6, env.dist.7, env.dist.8, env.dist.9, env.dist.10, env.dist.11, 
                    env.dist.12, env.dist.13, env.dist.14, env.dist.15, env.dist.16, env.dist.17,
                    env.dist.18, env.dist.19, env.dist.20, env.dist.21, env.dist.22, env.dist.23,
                    env.dist.24, env.dist.25, env.dist.26, env.dist.27, env.dist.28, env.dist.29,
                    env.dist.30, env.dist.31, env.dist.32, env.dist.33, env.dist.34, env.dist.35,
                    env.dist.36, env.dist.37, env.dist.38, env.dist.39, env.dist.40, env.dist.41, 
                    env.dist.42, env.dist.43, env.dist.44, env.dist.45, env.dist.46, env.dist.47, 
                    env.dist.48, env.dist.49, env.dist.50, env.dist.51, env.dist.52, env.dist.53, 
                    env.dist.54, env.dist.55, env.dist.56, env.dist.57, env.dist.58, geo.dist)

model_all <- '
env =~ env.dist.1+env.dist.2+env.dist.3+env.dist.4+env.dist.5+env.dist.6+env.dist.7+env.dist.8+env.dist.9+env.dist.10+env.dist.11+env.dist.12+env.dist.13+env.dist.14+env.dist.15+env.dist.16+env.dist.17+env.dist.18+env.dist.19+env.dist.20+env.dist.21+env.dist.22+env.dist.23+env.dist.24+env.dist.25+env.dist.26+env.dist.27+env.dist.28+env.dist.29+env.dist.30+env.dist.31+env.dist.32+env.dist.33+env.dist.34+env.dist.35+env.dist.36+env.dist.37+env.dist.38+env.dist.39+env.dist.40+env.dist.41+env.dist.42+env.dist.43+env.dist.44+env.dist.45+env.dist.46+env.dist.47+env.dist.48+env.dist.49+env.dist.50+env.dist.51+env.dist.52+env.dist.53+env.dist.54+env.dist.55+env.dist.56+env.dist.57+env.dist.58
geo =~ geo.dist

gen.dist ~ env + geo
geo ~~ env
env ~~ env
geo ~~ geo
'

fit_all <-  sem(model_all, data=full, estimator = "MLR", se="robust.huber.white")
s <- summary(fit_all, standardized=T, fit.measures=T, rsquare=T)
p <- parameterEstimates(fit_all)

summary(fit_all, rsquare=T)

##Env fit###
env_rda <- rda(filter_7, good_env)
e1 <- envfit(env_rda, good_env, choi=c(1:4), permutations = 9999)

env_rda_cond <- rda(filter_7, good_env, dbmem)
envfit(env_rda_cond, good_env, choi=c(1:4), permutations = 9999)

df <- data.frame(e1$vectors$arrows, e1$vectors$pvals)





####Variable selection, constraining on space####
filter_7_env_rda_cond <- rda(filter_7, full_env, dbmem)
(filter_7_env_rda_cond_r2a <- RsquareAdj(filter_7_env_rda_cond)$r.squared)

null_m <- rda(filter_7 ~ Condition(as.matrix(dbmem)), data = full_env)

full_m <- rda(filter_7 ~ . + Condition(as.matrix(dbmem)), full_env)

ordi <- ordistep(null_m, scope = formula(full_m), direction="forward")



cond_envs <- data.frame(full_env$cv_vapr, full_env$bio_15, full_env$txt.7, full_env$ann_max_rad, full_env$min_pet, full_env$max_pet)



grp_pimp <- find.clusters(filter_8_data_pimp_replace, max.n.clust = 140, n.pca = 141, n.clust = 5)


env_rda_cond <- rda(filter_7_data_pimp_imp, cond_envs, dbmem_10)

env_rda <- rda(filter_7_data_pimp_imp, cond_envs)

palette(alpha(funky(5)[c(2,5,1,3,4)], 0.8))

vec_col = alpha("red",0.6)

pdf("biplot_sites_nocond.pdf", width = 8, height = 7)

par(mar=c(6,6,6,6))

plot(scores(env_rda_cond, display='wa'), type="none", scaling=1, cex.axis=1.4, cex.lab=1.4, xlab="RDA1 site score", ylab = "RDA2 site score",las=1)
points(scores(env_rda_cond, display='wa'), display="sites", pch=19, cex=2.5, col=grp_pimp$grp, scaling=1)
abline(0,0, lty=2, col="gray50")
abline(v=0, lty=2, col="gray50")
labels <- c('CV vapr','Prec Seas.','Soil texture','Ann max rad','Min PET','Max PET')
vectors <- scores(env_rda_cond, display='bp')
t_vectors <- vectors

t_vectors[1,1] <- t_vectors[1,1] +0.1
t_vectors[1,2] <- t_vectors[1,2] +0.05
t_vectors[4,1] <- t_vectors[4,1] +0.2
t_vectors[4,2] <- t_vectors[4,2] -0.1
t_vectors[5,1] <- t_vectors[5,1] +0.1
t_vectors[5,2] <- t_vectors[5,2] -0.1


text(env_rda_cond, scaling=2, display="bp", col="black", cex=1.2, labels=NULL, font=2, lwd=2)
text(t_vectors[,1]*5, t_vectors[,2]*5, labels = labels, font=2)

axis(side=3, at = pretty(range(vectors[,1]*4),n=3), labels = pretty(range(vectors[,1]*4),n=3)/4, cex.axis=1.4, col=vec_col, col.ticks=vec_col,col.axis=vec_col)
mtext('Env Correlation RDA1', side=3, cex=1.4, adj=0.28, padj=-3.7, col=vec_col, font=1)
axis(side=4, at = pretty(range(vectors[,2]*4), n=3), labels = pretty(range(vectors[,2]*4),n=3)/4, cex.axis=1.4, las=1,col=vec_col, col.ticks=vec_col,col.axis=vec_col)
mtext('Env Correlation RDA2', side=4, cex=1.4, adj=0.67, padj=5.5, col=vec_col,font=1)
legend(5,7.5, legend=c("1", "2", "3", "4", "5"), title = "Genetic Cluster", fill=alpha(funky(5)[c(4,3,1,2,5)], 0.8),ncol=2)
dev.off()

#points(env_rda_cond, display="sites", pch=21, cex=1, col="gray70", bg='gray70', scaling=1)

library(maps)
library(maptools)

pdf("map_grp.pdf", width = 5, height = 5)
map("world","Peru", col="gray95", xlim=c(-93,-68.9), ylim=c(-18.2,1.4), fill=TRUE)
map("world","Ecuador", col="gray95",fill=TRUE, add=TRUE)
points(env_info$lon, env_info$lat, col=grp_pimp$grp, pch=19, cex=1.3)
dev.off()



snp_positions <- read.table("/Users/matthew/Box\ Sync/Projects/pimpGEA/associations/data/vcfs/stacks/v2/populations.snps.filter7.pimp.recode.snp.positions", sep="\t", header = T)


plot_rda <- function(model, snp_positions, outliers, ...){
  loadings_RDA1 <- data.frame(scores(model, choices = 1)$species)
  loadings_RDA1$SNP <- rownames(loadings_RDA1)
  loadings_RDA1$SNP <- gsub('.{2}$', '', loadings_RDA1$SNP)
  snp_positions <- snp_positions[order(snp_positions$ID),]
  
  loadings_RDA1 <- loadings_RDA1[order(loadings_RDA1$SNP),]
  
  snp_loadings_positions <- data.frame(snp_positions$CHROM, snp_positions$POS, loadings_RDA1$RDA1)
  snp_loadings_positions <- snp_loadings_positions[order(snp_loadings_positions$snp_positions.CHROM, snp_loadings_positions$snp_positions.POS),]
  snp_loadings_positions$loadings_RDA1.RDA1 <-snp_loadings_positions$loadings_RDA1.RDA1^2
  
  palette(c('gray80', 'gray60','gray80', 'gray60','gray80', 'gray60','gray80',
            'gray60','gray80', 'gray60','gray80', 'gray60'))
  
  tops <- which.maxn(loadings_RDA1$RDA1^2, outliers)
  tops <- loadings_RDA1[tops,]
  m <- min(tops$RDA1^2)
  ma <- max(tops$RDA1^2)
  snp_loadings_positions$snp_positions.CHROM <- ifelse(snp_loadings_positions$loadings_RDA1.RDA1 <= ma & snp_loadings_positions$loadings_RDA1.RDA1 >= m, funky(4)[3], snp_loadings_positions$snp_positions.CHROM)
  print(tops)
  
  #snp_loadings_positions$cex <- scales::rescale(snp_loadings_positions$loadings_RDA1.RDA1, to = c(0.5, 1))
  plot(seq(1,length(snp_loadings_positions$snp_positions.CHROM)), snp_loadings_positions$loadings_RDA1.RDA1, 
       col=snp_loadings_positions$snp_positions.CHROM, pch=19, cex=0.5,xaxt='n', yaxt='n', ...)
  axis(side=1, at=c(2645,7393,11550,15503,18841,22138,25678,28906,32153,35780,39327,42640), mgp=c(3,0.6,0), 
       labels = c("Ch 1", "Ch 2", "Ch 3", "Ch 4", "Ch 5", "Ch 6", "Ch 7", "Ch 8", "Ch 9", "Ch 10", "Ch 11", "Ch 12" ))
  print(pretty(c(0,max(snp_loadings_positions$loadings_RDA1.RDA1)), n = 3))
  axis(side=2, at=pretty(c(0,max(snp_loadings_positions$loadings_RDA1.RDA1)), n = 3), labels = pretty(c(0,max(snp_loadings_positions$loadings_RDA1.RDA1)), n = 3))
  #abline(cutoff, 0, lty=2, cex=2)
}




pdf("RDA_1_by_snp.pdf", width = 10, height = 6)
par(mfrow=c(2,1), mar=c(3,4,2,1))
plot_rda(env_rda, snp_positions, 157, xlab="", ylab="Squared RDA1 loading", main="Full RDA")
abline(h = 0.0373, lty = 2)
plot_rda(env_rda_cond, snp_positions, 252, xlab= "", ylab="Squared RDA1 loading", main = "RDA constrained on space",
         ylim=c(0,0.015))
abline(h = 0.00777, lty = 2)

dev.off()


loadings_cond <- data.frame(scores(env_rda_cond, choices = 1)$species)

loadings_cond$SNP <- rownames(loadings_cond)
loadings_cond$SNP <- gsub('.{2}$', '', loadings_cond$SNP)
snp_positions <- snp_positions[order(snp_positions$ID),]

loadings_cond <- loadings_cond[order(loadings_cond$SNP),]

snp_loadings_positions <- data.frame(snp_positions$CHROM, snp_positions$POS, snp_positions$ID, loadings_cond$RDA1)
names(snp_loadings_positions) <- c('CHROM', 'POS', 'SNP', 'RDA1')



###OUTLIERS
outliers <- function(x,z){
  lims <- mean(x$RDA1) + c(-1, 1) * z * sd(x$RDA1)     # find loadings +/-z sd from mean loading     
  
  x[x$RDA1 < lims[1] | x$RDA1 > lims[2],]     # locus names in these tails
}


#Cond
loadings_cond$SNP <- rownames(loadings_cond)

RDA1_candidates1 <- outliers(snp_loadings_positions,4)
RDA1_candidates1$RDA1 <- RDA1_candidates1$RDA1^2
RDA1_candidates1 <- RDA1_candidates1[order(RDA1_candidates1$RDA1, decreasing=T),]




#No cond
loadings_cond <- data.frame(scores(env_rda, choices = 1)$species)

loadings_cond$SNP <- rownames(loadings_cond)
loadings_cond$SNP <- gsub('.{2}$', '', loadings_cond$SNP)
snp_positions <- snp_positions[order(snp_positions$ID),]

loadings_cond <- loadings_cond[order(loadings_cond$SNP),]

snp_loadings_positions <- data.frame(snp_positions$CHROM, snp_positions$POS, snp_positions$ID, loadings_cond$RDA1)
names(snp_loadings_positions) <- c('CHROM', 'POS', 'SNP', 'RDA1')

loadings_cond$SNP <- rownames(loadings_cond)

RDA1_candidates2 <- outliers(snp_loadings_positions,4)
RDA1_candidates2$RDA1 <- RDA1_candidates2$RDA1^2
RDA1_candidates2 <- RDA1_candidates2[order(RDA1_candidates2$RDA1, decreasing=T),]



#write.table(RDA1_candidates, 'candidate_loci.txt', sep='\t')


foo <- matrix(nrow=(dim(RDA1_candidates)[1]), ncol=6)  # 8 columns for 8 predictors
colnames(foo) <- colnames(cond_envs)

for (i in 1:length(RDA1_candidates[,3])) {
  nam <- RDA1_candidates[i,3]
  snp.gen <- filter_7_data_pimp_imp[,nam]
  foo[i,] <- apply(cond_envs,2,function(x) cor(x,snp.gen))
}

cand <- cbind.data.frame(RDA1_candidates,foo)  
head(cand)

#dups <- cand$snp[duplicated(cand$snp)]
#cand <- cand[!duplicated(cand$snp),]



ns <- c()
cor <- c()

for (i in 1:length(cand$SNP)) {
  bar <- cand[i,]
  bar <- bar[4:ncol(bar)]
  #bar <- bar[,c(-9,-10)]
  ns <-  rbind(ns, names(which.max(abs(bar)))) # gives the variable
  cor <- rbind(cor, max(abs(bar)))           # gives the correlation
}

cand$predictor <- ns
cand$correlation <- cor

