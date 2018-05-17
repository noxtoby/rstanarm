  // impose consraint on bMat1 such that random intercepts sum to zero for each subject
  // this allow identifiability of latent time shift parameter
  for(i in 1:bN1) {
    for(j in 1:(M-1)){
      bMat10[i,j] = bMat1[i, bK1_idx[j,1]];
    }
    bMat1[i, bK1_idx[M,1]] = -sum(bMat10[i,]);
  }