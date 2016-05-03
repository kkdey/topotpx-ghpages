
library(ashr)
library(reshape)
library(reshape2)
library(Matrix)

## The function for correlation matrix ash framework

ash_cor <- function(cormat, nsamples)
{
  cor_table <- melt(cormat);
  cor_table_non_diag <- cor_table[which(cor_table[,1] !=cor_table[,2]),];

  cor_table_non_diag.val <- cor_table_non_diag[,3];
  cor_table_non_diag.val[which(cor_table_non_diag.val==1)]=(1- 1e-7);
  #cor_table_non_diag.val[which(cor_table_non_diag.val==0)]=(1e-7);

  cor_transform_mean_vec=0.5*log((1+cor_table_non_diag.val)/(1-cor_table_non_diag.val))
  cor_transform_sd_vec=rep(sqrt(1/(nsamples-3)), dim(cor_table_non_diag)[1]);
  options(warn=-1)
  fit=ash(cor_transform_mean_vec,cor_transform_sd_vec,mixcompdist="normal");
  ash_cor_vec=(exp(2*fit$PosteriorMean)-1)/(exp(2*fit$PosteriorMean)+1);

  newdata.table <- cor_table_non_diag;
  newdata.table[,3] <- ash_cor_vec;
  new_mat <- dcast(newdata.table, X1~X2, value.var = "value")[,-1];
  new_mat[is.na(new_mat)]=1;
  pd_completion <- nearPD(as.matrix(new_mat), conv.tol=1e-06);
  new_mat <- sweep(pd_completion$mat,diag(pd_completion$mat), MARGIN=1,"/")
  return(new_mat)
}

shafer_shrink_master <- function(data){
  nsamples <- dim(data)[1];
  nfeat <- dim(data)[2];
  colmean <- colMeans(data);
  covmat <- cov(data);
  cormat <- cov2cor(covmat);
  w <- array(0, c(nsamples, nfeat, nfeat));
  for(n in 1:nsamples){
    w[n,,] <- (data[n,] - colmean)%*% t(data[n,] - colmean);
  }
  var_hat_s <- (nsamples/(nsamples-1)^2) * apply(w, c(2,3), var);
  sum_var_hat_s <- sum(var_hat_s[row(var_hat_s)!=col(var_hat_s)])
  square_cor <- covmat^2;
  sum_s_square <- sum(square_cor[row(square_cor)!=col(square_cor)]);
  shrink_intensity <- sum_var_hat_s/sum_s_square;

  if(shrink_intensity > 1){
    shafer_shrink <- diag(diag(covmat));
  }
  if(shrink_intensity < 0){
    shafer_shrink <- covmat;
  }
  else{
    shafer_shrink <- diag(rep(1, nfeat))
    shafer_shrink[row(shafer_shrink)!=col(shafer_shrink)] <- (1- shrink_intensity)*covmat[row(covmat)!=col(covmat)]
  }
  return(shafer_shrink)
}


ash_cor_master <- function(data){
  nsamples <- dim(data)[1];
  nfeat <- dim(data)[2];
  colmean <- colMeans(data);
  covmat <- cov(data);
  cormat <- cov2cor(covmat);
  w <- array(0, c(nsamples, nfeat, nfeat));
  for(n in 1:nsamples){
    w[n,,] <- (data[n,] - colmean)%*% t(data[n,] - colmean);
  }
  var_hat_s <- (nsamples/(nsamples-1)^2) * apply(w, c(2,3), var);
  sum_var_hat_s <- sum(var_hat_s[row(var_hat_s)!=col(var_hat_s)])
  square_cor <- covmat^2;
  sum_s_square <- sum(square_cor[row(square_cor)!=col(square_cor)]);
  shrink_intensity <- sum_var_hat_s/sum_s_square;
  ash_cormat <- ash_cor(cormat, nsamples);

  if(shrink_intensity < 0){
    ash_cor_master <- cormat;
  }
  if(shrink_intensity > 1){
    ash_cor_master <- ash_cormat
  }
  else{
    ash_cormat_master <- shrink_intensity*ash_cormat + (1- shrink_intensity)*cormat;

  }
  ash_covmat_master <- diag(sqrt(diag(covmat)))%*%ash_cormat_master%*%diag(sqrt(diag(covmat)))
  ll <- list("ash_cov_master"=ash_covmat_master,
             "ash_cor_master"=ash_cormat_master,
             "shrink_intensity"=shrink_intensity,
             "ash_cormat"=ash_cormat,
             "sample_cormat"=cormat)
  return(ll)
}


generate_cov <- function(dim){
  pop_cov_class <- clusterGeneration::genPositiveDefMat(dim, covMethod="eigen")
  pop_cov_eigen <- pop_cov_class$egvalues;
  pop_cov_Sigma <- pop_cov_class$Sigma;
  ll <- list("Sigma"=pop_cov_Sigma, "eigen"=pop_cov_eigen);
  return(ll)
}

