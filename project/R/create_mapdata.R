
## Create a data matrix across cells along with metadata for easy handling of the dispersion field maps

files <- list.files(path=Sys.glob("../external_data/Dispersion_Fields_with_Abundance/"));
# files <- list.files(path=Sys.glob("../external_data/Dispersion_Fields_No_Abund/"));

map_data <- vector()
for(f in 1:length(files))
{
  data <- read.csv(paste0("../external_data/Dispersion_Fields_with_Abundance/",files[f]));
 # data <- read.csv(paste0("../external_data/Dispersion_Fields_No_Abund/",files[f]));

  data[is.na(data)] <- 0;
  data <- data[,-1];
  data_vec <- matrix(as.matrix(data), ncol=dim(data)[1]*dim(data)[2], byrow=TRUE);
  map_data <- rbind(map_data, data_vec);
}

spots <- unlist(lapply(1:length(files), function(f) strsplit(strsplit(files[f],"_")[[1]][2],"[.]")[[1]][1]))
rownames(map_data) <- spots;
write.table(map_data, "../external_data/map_data_abundance.txt")
# write.table(map_data, "../external_data/map_data_no_abundance.txt")
write.csv(map_data, "../external_data/map_data_abundance.csv")
# write.csv(map_data, "../external_data/map_data_no_abundance.csv")
