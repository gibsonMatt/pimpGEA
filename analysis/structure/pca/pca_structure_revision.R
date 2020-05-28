setwd("/Users/matthew/Box\ Sync/Projects/pimpGEA/associations/dataanalysis/stacks_analysis/structure/pca")


library(adegenet)
library(ape)
library(maps)
library(maptools)
library(scales)
library(ellipse)

###Keeping only sites with <5% missing data

filter_8_data_pimp = read.PLINK("/Users/matthew/Box\ Sync/Projects/pimpGEA/associations/data/plinks/stacks/v2/populations.snps.filter8.pimp.raw")

#population assignments
pops_2 <- read.csv("./popmap.csv")

#Remove suffix from names
#names <- as.vector(filter_3_data_pimp$ind.names)
#names <- gsub(".filtered.sorted", "", names)
#filter_3_data_pimp$ind.names <- names


countNA <- function(row){
  return(sum(is.na(row))/142)
}

gen <- as.matrix(filter_8_data_pimp)
nas <- apply(gen, 2, countNA)

keep <- names(nas[nas < 0.05])

remove <- c('2933-101', '365911-101', '503517-101')
filter_8_data_pimp <- filter_8_data_pimp[!rownames(filter_8_data_pimp) %in% remove, ]

df <- filter_8_data_pimp[,which(filter_8_data_pimp$loc.names %in% keep)]

#pdf('missing_data_hist.pdf', 5,4)
hist(nas, breaks=25, col = "gray80", border = "gray40", main="", xlab = "Proportion of individuals with missing data")
#dev.off()

df <- new("genlight", df)

df <- tab(df, NA.method = 'mean')
df <- new("genlight", df)
df$pop <- pops_2$Species


##PCA within pimp
pca_pimp <-glPca(df)

myCol <- colorplot(pca_pimp$scores, pca_pimp$scores, transp=T, cex=3)
abline(h=0,v=0, col="grey")
#add.scatter.eig(pca1$eig[1:40],2,1,2, posi="topright", inset=.05, ratio=.3)

plot(pca_pimp$scores[,1], pca_pimp$scores[,2], col=funky(6)[3], pch=19, cex=2, ylab = "PC2", xlab="PC1")

pca_pimp$eig[1]/sum(pca_pimp$eig)
pca_pimp$eig[2]/sum(pca_pimp$eig)

##DAPC within pimp

grp_pimp <- find.clusters(df, max.n.clust = 140)


pimp_info <- read.csv("./pimp_info.csv")

pdf('struc_map2.pdf', 7,10)
map("world","Peru", col="gray90", xlim=c(-84,-69), ylim=c(-18.4,1.4), fill=TRUE)
map("world","Ecuador", col="gray90",fill=TRUE, add=TRUE)
#map("world","Brazil", col="gray90",fill=TRUE, add=TRUE) 
#map("world","Colombia", col="gray90",fill=TRUE, add=TRUE) 
palette(alpha(funky(7)[c(4,6,2,7,1)], 0.99))
#palette(alpha(funky(6), 0.99))
#palette(alpha(c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99'),0.99))
points(pimp_info$lon, pimp_info$lat, pch=21, bg=grp_pimp$grp, cex=3)
#legend("bottomright", levels(grp_pimp$grp), col=palette(), pch=19)
dev.off()

#pdf('pc_2.pdf', 6,6)
plot(pca_pimp$scores[,1], pca_pimp$scores[,2], pch=21, bg=grp_pimp$grp, cex=2.6, xlab="PC1 (18.03%)", ylab="PC2 (6.53%)", xlim = c(-14,27), ylim=c(-12,11))
dataEllipse(pca_pimp$scores[,1], pca_pimp$scores[,2], grp_pimp$grp, levels = 0.90, add=T, group.labels = NULL, 
            center.pch = NULL, plot.points = F, col = palette())
#dev.off()



setwd("/Users/matthew/Box\ Sync/Projects/pimpGEA/associations/dataanalysis/stacks_analysis/structure/faststructure/rev_0.05missing")


k4 <- read.table('./output/revstruc.4.meanQ')
n <- k4$V1
groups <- match(k4$V1, names(grp_pimp$grp))

k4$group <- grp_pimp$grp[groups]
group <- grp_pimp$grp[groups]

k4 <- k4[order(k4$group),]
a <- k4
k4 <- k4[,c(-6, -1)]
k4 <- t(k4)

pdf('struc_plot.pdf', 10,3)
palette(alpha(funky(7)[c(2,1,4,6,7)], 0.99))
barplot(k4, col = palette(), border=NA, space = 0)
dev.off()


setwd("/Users/matthew/Box\ Sync/Projects/pimpGEA/associations/dataanalysis/stacks_analysis/structure/faststructure/rev_0.05missing_withld")

k4a <- read.table('./outputs/revstruc.3.meanQ')
names(k4a) <- c('V2', 'V3', 'V4')
k4a$V1 <- n

groups <- match(k4a$V1, names(grp_pimp$grp))

k4a$group <- grp_pimp$grp[groups]
group <- grp_pimp$grp[groups]

k4a <- k4a[order(k4a$group),]
b <- k4a
k4a <- k4a[,c(-4,-5)]
k4a <- t(k4a)

pdf('struc_plot_comb.pdf', 10,6)
par(mfrow=c(2,1))
palette(alpha(funky(7)[c(2,1,4,6,7)], 0.99))
barplot(k4, col = palette(), border=NA, space = 0)
palette(alpha(bluepal(4), 0.99))
barplot(k4a, col = palette(), border=NA, space = 0)
dev.off()
