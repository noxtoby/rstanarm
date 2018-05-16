real mean_PPD[M];
real yAlpha1[intercept_type[1] > 0];
real yAlpha2[intercept_type[2] > 0];
real yAlpha3[intercept_type[3] > 0];
vector[prior_dist_for_cov == 2 && bK1 > 0 ? size(bCov1_idx) : 0] bCov1;
vector[prior_dist_for_cov == 2 && bK2 > 0 ? size(bCov2_idx) : 0] bCov2;
vector[bN1 * bK1] b1 = to_vector(bMat1'); // ensures same order as stan_glmer (make_b)
vector[bN2 * bK2] b2 = to_vector(bMat2');

matrix[yNeta[1],yK[1]] yXlt1; // fe design matrix with latent time shift
matrix[yNeta[2],yK[2]] yXlt2;
matrix[yNeta[3],yK[3]] yXlt3;

vector[yK[1]] yXltbar1; // predictor (with latent time) means
vector[yK[2]] yXltbar2;
vector[yK[3]] yXltbar3;

// Evaluate mean_PPD
{
  int bMat1_colshift = 0; // column shift in bMat1
  int bMat2_colshift = 0; // column shift in bMat2

  // Linear predictor for submodel 1
  if (M > 0) {
    vector[yNeta[1]] yEta1; // linear predictor
    for(i in 1:yK[1]){
      yXlt1[,i] = yX1[,i];
    }
    for(j in 1:yNeta[1]){
      yXlt1[j,lt_idx[1]] = yXlt1[j,lt_idx[1]] + Delta[y1_Z1_id[j]];
    }
    yEta1 = evaluate_eta(yXlt1, y1_Z1, y1_Z2, y1_Z1_id, y1_Z2_id, yGamma1, yBeta1,
                         bMat1, bMat2, bMat1_colshift, bMat2_colshift, intercept_type[1]);
    yEta1 = evaluate_mu(yEta1, family[1], link[1]);
    mean_PPD[1] = mean_PPD_rng(yEta1, yAux1, family[1]);
  }

  // Linear predictor for submodel 2
  if (M > 1) {
    vector[yNeta[2]] yEta2;
    bMat1_colshift = bMat1_colshift + bK1_len[1];
    bMat2_colshift = bMat2_colshift + bK2_len[1];
    for(i in 1:yK[2]){
      yXlt2[,i] = yX2[,i];
    }
    for(j in 1:yNeta[2]){
      yXlt2[j,lt_idx[2]] = yXlt2[j,lt_idx[2]] + Delta[y2_Z1_id[j]];
    }
    yEta2 = evaluate_eta(yXlt2, y2_Z1, y2_Z2, y2_Z1_id, y2_Z2_id, yGamma2, yBeta2,
                         bMat1, bMat2, bMat1_colshift, bMat2_colshift, intercept_type[2]);
    yEta2 = evaluate_mu(yEta2, family[2], link[2]);
    mean_PPD[2] = mean_PPD_rng(yEta2, yAux2, family[2]);
  }

  // Linear predictor for submodel 3
  if (M > 2) {
    vector[yNeta[3]] yEta3;
    bMat1_colshift = bMat1_colshift + bK1_len[2];
    bMat2_colshift = bMat2_colshift + bK2_len[2];
    for(i in 1:yK[3]){
      yXlt3[,i] = yX3[,i];
    }
    for(j in 1:yNeta[3]){
      yXlt3[j,lt_idx[3]] = yXlt3[j,lt_idx[3]] + Delta[y3_Z1_id[j]];
    }
    yEta3 = evaluate_eta(yXlt3, y3_Z1, y3_Z2, y3_Z1_id, y3_Z2_id, yGamma3, yBeta3,
                         bMat1, bMat2, bMat1_colshift, bMat2_colshift, intercept_type[3]);
    yEta3 = evaluate_mu(yEta3, family[3], link[3]);
    mean_PPD[3] = mean_PPD_rng(yEta3, yAux3, family[3]);
  }
}

// Transform intercept parameters
if (intercept_type[1] > 0){
  for(i in 1:yK[1]){
    yXltbar1[i] = mean(yXlt1[,i]);
  }
  yAlpha1[1] = yGamma1[1] - dot_product(yXltbar1, yBeta1);
}
if (M > 1 && intercept_type[2] > 0){
  for(i in 1:yK[2]){
    yXltbar2[i] = mean(yXlt2[,i]);
  }
  yAlpha2[1] = yGamma2[1] - dot_product(yXltbar2, yBeta2);
}
if (M > 2 && intercept_type[3] > 0){
  for(i in 1:yK[3]){
    yXltbar3[i] = mean(yXlt3[,i]);
  }
  yAlpha3[1] = yGamma3[1] - dot_product(yXltbar3, yBeta3);
}

  // Transform variance-covariance matrices

    // Grouping factor 1
  if (prior_dist_for_cov == 2 && bK1 == 1) {
    bCov1[1] = bSd1[1] * bSd1[1];
  }
  else if (prior_dist_for_cov == 2 && bK1 > 1) {
    bCov1 = to_vector(quad_form_diag(
      multiply_lower_tri_self_transpose(bCholesky1), bSd1))[bCov1_idx];
  }

  // Grouping factor 2
  if (prior_dist_for_cov == 2 && bK2 == 1) {
    bCov2[1] = bSd2[1] * bSd2[1];
  }
  else if (prior_dist_for_cov == 2 && bK2 > 1) {
    bCov2 = to_vector(quad_form_diag(
      multiply_lower_tri_self_transpose(bCholesky2), bSd2))[bCov2_idx];
  }
