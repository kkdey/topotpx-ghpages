rm(list=ls())
library(CountClust)
N=360;
M=560;
lat <- rep(1:N, M);
long <- rep(1:M, each=N)

comp <- cbind(lat, long);

Topic_clus <- readRDS("../rdas/topic_clus_2_maps_no_abundance.rds")
center_gravity <- cg_topics(Topic_clus$theta, comp);

Topic_clus <- readRDS("../rdas/topic_clus_3_maps_no_abundance.rds")
center_gravity <- cg_topics(Topic_clus$theta, comp);

Topic_clus <- readRDS("../rdas/topic_clus_2_maps_abundance.rds")
center_gravity <- cg_topics(Topic_clus$theta, comp);

Topic_clus <- readRDS("../rdas/topic_clus_3_maps_abundance.rds")
center_gravity <- cg_topics(Topic_clus$theta, comp);

