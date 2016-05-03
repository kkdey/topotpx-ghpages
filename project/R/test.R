dim <- 1000
pop_cov <- generate_cov(1000);
plot(sort(pop_cov$eigen, decreasing = TRUE), type="l")
plot(sort(pop_cov$eigen, decreasing = TRUE), type="l")

generate_sample <- mvtnorm::rmvnorm(200, rep(0, dim), pop_cov$Sigma);
cov_sample <- cov(generate_sample)
eigen_sample <- eigen(cov_sample, only.values = TRUE)
plot(sort(as.numeric(eigen_sample$values), decreasing = TRUE), type="l")
cor_sample <- cov2cor(cov_sample)
image(nearPD(as.matrix(cor_sample), conv.tol = 1e-06)$mat)

shafer_mat <- shafer_shrink_master(generate_sample);
shafer_eigen <- eigen(shafer_mat, only.values = TRUE)

plot(sort(pop_cov$eigen, decreasing = TRUE), type="l", col="black")
lines(sort(as.numeric(shafer_eigen$values), decreasing = TRUE), type="l", col="blue")
lines(sort(as.numeric(eigen_sample$values), decreasing = TRUE), type="l", col="green")

ash_cor_sample <- ash_cor(cor_sample, 200);
#ash_cor_sample[abs(ash_cor_sample)  < 0.0001] <- 0
image(nearPD(as.matrix(ash_cor_sample), conv.tol = 1e-06)$mat)

# samp_length <- c(10, 50, 200, 300, 500, 5);
# cols <- c("red", "green", "yellow", "blue", "cyan", "brown")

# plot(sort(as.numeric(pop_cov$eigen), decreasing = TRUE), type="l", col="black")

# for(num in 1:length(samp_length)){
#  ash_cor_sample <- ash_cor(cor_sample, samp_length[num]);
#  ash_cov_sample <- diag(sqrt(diag(cov_sample)))%*%ash_cor_sample%*%diag(sqrt(diag(cov_sample)))
#  eigen_ash_sample <- eigen(ash_cov_sample, only.values = TRUE)
#  lines(sort(as.numeric(eigen_ash_sample$values), decreasing = TRUE), type="l", col=cols[num])
#}


ash_cov_sample <- diag(sqrt(diag(cov_sample)))%*%ash_cor_sample%*%diag(sqrt(diag(cov_sample)))
eigen_ash_sample <- eigen(ash_cov_sample, only.values = TRUE)
plot(sort(as.numeric(eigen_ash_sample$values), decreasing = TRUE), type="l", col="red")
lines(sort(as.numeric(eigen_sample$values), decreasing = TRUE), type="l", col="green")
lines(sort(as.numeric(shafer_eigen$values), decreasing = TRUE), type="l", col="blue")
lines(sort(as.numeric(pop_cov$eigen), decreasing = TRUE), type="l", col="black")

ash_cor_master_sample <- ash_cor_master(generate_sample);
shrinkage <- ash_cor_master_sample$shrink_intensity
ash_cov_master <- ash_cor_master_sample$ash_cov_master;
eigen_ash_cov_master <- eigen(ash_cov_master);
plot(sort(as.numeric(eigen_ash_cov_master$values), decreasing = TRUE), type="l", col="red")
lines(sort(as.numeric(pop_cov$eigen), decreasing = TRUE), type="l", col="black")
lines(sort(as.numeric(eigen_sample$values), decreasing = TRUE), type="l", col="green")
lines(sort(as.numeric(eigen_ash_sample$values), decreasing = TRUE), type="l", col="yellow")
