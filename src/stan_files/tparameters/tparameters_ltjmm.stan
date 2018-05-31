  // impose consraint on bMat1 such that random intercepts sum to zero for each subject
  // this allow identifiability of latent time shift parameter
  for(i in 1:bN1) {
    for(j in 1:(M-1)){
      bMat10[i,j] = bMat1[i, bK1_idx[j,1]];
    }
    bMat1[i, bK1_idx[M,1]] = -sum(bMat10[i,]);
  }
  
  sigma_lt = make_aux(sigma_lt_unscaled, y_prior_dist_for_sigma_lt,
                      y_prior_mean_for_sigma_lt, y_prior_scale_for_sigma_lt);