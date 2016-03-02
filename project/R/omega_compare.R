
## test reorder_topics

Topic_clus_1 <- readRDS("../rdas/topic_clus_3_maps_abundance.rds")
omega1 <- Topic_clus_1$omega;

Topic_clus_2 <- readRDS("../rdas/topic_clus_3_maps_no_abundance.rds")
omega2 <- Topic_clus_2$omega;

compare.omega  <- compare_omega(omega1, omega2)

