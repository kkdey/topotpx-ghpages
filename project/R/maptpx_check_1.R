

## maptpx fit for Alex Data

library(HimalayanBirdsAbundance)
data(HimalayanBirdsAbundance)
counts <- t(exprs(HimalayanBirdsAbundance));
metadata <- pData(HimalayanBirdsAbundance)

counts <- counts[, which(colSums(counts)!=0)];

library(maptpx)
library(slam)

source("../../../maptpx/R/topics.R")
source("../../../maptpx/R/tpx.R")
source("../../../maptpx/R/count.R")

K=3
system.time(suppressWarnings(Topic_clus <- topics(counts, K=K, tol=0.00001, kill=0, bf=T, light=TRUE)))

elevation_metadata=metadata$Elevation;
east_west_dir = metadata$WorE;

docweights <- Topic_clus$omega

par(mfrow=c(1,1))
par(mar=c(5,2,2,2))
barplot(t(Topic_clus$omega[order(elevation_metadata),]),col=2:(K+1),axisnames=F,space=0,border=NA,main="",las=1,ylim=c(0,1),cex.axis=1.5,cex.main=1.4, xlab="elevation")
title(main=paste("Structure Plot topic proportions,k=",K))
combo_patch_dir = paste0(round(elevation_metadata));
combo_patch_dir_ordered = combo_patch_dir[order(elevation_metadata)];

match_labs=match(unique(combo_patch_dir_ordered),combo_patch_dir_ordered);
match_labs_suffix=c(match_labs[2:length(unique(combo_patch_dir_ordered))],35);
match_labs_prefix=match_labs[1:(length(unique(combo_patch_dir_ordered)))];
labs=match_labs_prefix + 0.5*(match_labs_suffix - match_labs_prefix);

axis(1,at=labs,unique(combo_patch_dir_ordered),las=2);
abline(v=match_labs[2:length(match_labs)]);


docweights <- Topic_clus$omega
east_west_elevation = paste0(metadata$WorE, "_", metadata$Elevation);

index1 <- which(metadata$WorE=="E");
index2 <- which(metadata$WorE=="W");
elevation1 <- metadata$Elevation[index1]; elevation2 <- metadata$Elevation[index2];
index_WE <- c(index1[order(elevation1)], index2[order(elevation2)]);

barplot(t(docweights[index_WE,]),col=2:(K+1),axisnames=F,space=0,border=NA,main=paste("No. of clusters=",K),las=1,ylim=c(0,1),cex.axis=1.5,cex.main=1.4)

combo_patch_dir = paste0(east_west_elevation);
combo_patch_dir_ordered = combo_patch_dir[index_WE];

match_labs=match(unique(combo_patch_dir_ordered),combo_patch_dir_ordered);
match_labs_suffix=c(match_labs[2:length(unique(combo_patch_dir_ordered))],35);
match_labs_prefix=match_labs[1:(length(unique(combo_patch_dir_ordered)))];
labs=match_labs_prefix + 0.5*(match_labs_suffix - match_labs_prefix);

axis(1,at=labs,unique(combo_patch_dir_ordered),las=2);
abline(v=match_labs[2:length(match_labs)]);



#####################  Ordering affects the results?   #########################

library(slam)
library(maptpx)

for (i in 1:20){
  rand_counts <- counts[,sample(ncol(counts), replace=FALSE)]
  topics.4 <- topics(rand_counts, K=4, tol=0.000001, kill=0, bf=T, light=0, 
                     init.adapt=TRUE)
  docweights <- topics.4$omega
  east_west_elevation = paste0(metadata$WorE, "_", metadata$Elevation);
  K <- 4
  index1 <- which(metadata$WorE=="E");
  index2 <- which(metadata$WorE=="W");
  elevation1 <- metadata$Elevation[index1]; elevation2 <- metadata$Elevation[index2];
  index_WE <- c(index1[order(elevation1)], index2[order(elevation2)]);
  
  barplot(t(docweights[index_WE,]),col=2:(K+1),axisnames=F,space=0,border=NA,main=paste("No. of clusters=",K),las=1,ylim=c(0,1),cex.axis=1.5,cex.main=1.4)
  
  combo_patch_dir = paste0(east_west_elevation);
  combo_patch_dir_ordered = combo_patch_dir[index_WE];
  
  match_labs=match(unique(combo_patch_dir_ordered),combo_patch_dir_ordered);
  match_labs_suffix=c(match_labs[2:length(unique(combo_patch_dir_ordered))],35);
  match_labs_prefix=match_labs[1:(length(unique(combo_patch_dir_ordered)))];
  labs=match_labs_prefix + 0.5*(match_labs_suffix - match_labs_prefix);
  
  axis(1,at=labs,unique(combo_patch_dir_ordered),las=2);
  abline(v=match_labs[2:length(match_labs)]);
  
}


