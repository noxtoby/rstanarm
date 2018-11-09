  vector[yK[1]] yBeta1; // population level params
  vector[yK[2]] yBeta2;
  vector[yK[3]] yBeta3;
  vector[yK[4]] yBeta4;
  vector[yK[5]] yBeta5;
  vector[yK[6]] yBeta6;
  vector[yK[7]] yBeta7;
  vector[yK[8]] yBeta8;
  vector[yK[9]] yBeta9;
  vector[yK[10]] yBeta10;
  vector[yK[11]] yBeta11;
  vector[yK[12]] yBeta12;
  vector[yK[13]] yBeta13;
  vector[yK[14]] yBeta14;
  vector[yK[15]] yBeta15;
  vector[yK[16]] yBeta16;
  vector[yK[17]] yBeta17;
  vector[yK[18]] yBeta18;
  vector[yK[19]] yBeta19;
  vector[yK[20]] yBeta20;
  real yAux1[has_aux[1]]; // auxiliary params
  real yAux2[has_aux[2]];
  real yAux3[has_aux[3]];
  real yAux4[has_aux[4]];
  real yAux5[has_aux[5]];
  real yAux6[has_aux[6]];
  real yAux7[has_aux[7]];
  real yAux8[has_aux[8]];
  real yAux9[has_aux[9]];
  real yAux10[has_aux[10]];
  real yAux11[has_aux[11]];
  real yAux12[has_aux[12]];
  real yAux13[has_aux[13]];
  real yAux14[has_aux[14]];
  real yAux15[has_aux[15]];
  real yAux16[has_aux[16]];
  real yAux17[has_aux[17]];
  real yAux18[has_aux[18]];
  real yAux19[has_aux[19]];
  real yAux20[has_aux[20]];
  vector[len_theta_L] theta_L; // cov matrix for decov prior
  real yAuxMaximum = 1.0; // used for scaling in theta_L

  // group level params
  matrix[bK1 >  0 ? bN1 : 0, bK1] bMat1; // for grouping factor 1
  matrix[bK2 >  0 ? bN2 : 0, bK2] bMat2; // for grouping factor 2

  // population level params, auxiliary params
  if (has_aux[1] == 1) {
    yAux1[1] = make_aux(yAux1_unscaled[1], y_prior_dist_for_aux[1],
                        y_prior_mean_for_aux[1], y_prior_scale_for_aux[1]);
    if (yAux1[1] > yAuxMaximum)
      yAuxMaximum = yAux1[1];
  }

  if (yK[1] > 0){
    yBeta1 = make_beta(z_yBeta1, y_prior_dist[1], y_prior_mean1,
                       y_prior_scale1, y_prior_df1, y_global_prior_scale[1],
                       yGlobal1, yLocal1, yOol1, yMix1, yAux1, family[1],
                       y_slab_scale[1], y_caux1);
    yBeta1[lt_idx[1]] = abs(yBeta1[lt_idx[1]]);
  }
  if (M > 1) {
    if (has_aux[2] == 1) {
      yAux2[1] = make_aux(yAux2_unscaled[1], y_prior_dist_for_aux[2],
                          y_prior_mean_for_aux[2], y_prior_scale_for_aux[2]);
      if (yAux2[1] > yAuxMaximum)
        yAuxMaximum = yAux2[1];
    }
    if (yK[2] > 0){
      yBeta2 = make_beta(z_yBeta2, y_prior_dist[2], y_prior_mean2,
                         y_prior_scale2, y_prior_df2, y_global_prior_scale[2],
                         yGlobal2, yLocal2, yOol2, yMix2, yAux2, family[2],
                         y_slab_scale[2], y_caux2);
      yBeta2[lt_idx[2]] = abs(yBeta2[lt_idx[2]]);
    }
  }
  if (M > 2) {
    if (has_aux[3] == 1) {
      yAux3[1] = make_aux(yAux3_unscaled[1], y_prior_dist_for_aux[3],
        y_prior_mean_for_aux[3], y_prior_scale_for_aux[3]);
      if (yAux3[1] > yAuxMaximum)
        yAuxMaximum = yAux3[1];
    }
    if (yK[3] > 0){
      yBeta3 = make_beta(z_yBeta3, y_prior_dist[3], y_prior_mean3,
        y_prior_scale3, y_prior_df3, y_global_prior_scale[3],
        yGlobal3, yLocal3, yOol3, yMix3, yAux3, family[3],
        y_slab_scale[3], y_caux3);
      yBeta3[lt_idx[3]] = abs(yBeta3[lt_idx[3]]);
    }
  } 
  if (M > 3) {
    if (has_aux[4] == 1) {
      yAux4[1] = make_aux(yAux4_unscaled[1], y_prior_dist_for_aux[4],
        y_prior_mean_for_aux[4], y_prior_scale_for_aux[4]);
      if (yAux4[1] > yAuxMaximum)
        yAuxMaximum = yAux4[1];
    }
    if (yK[4] > 0){
      yBeta4 = make_beta(z_yBeta4, y_prior_dist[4], y_prior_mean4,
        y_prior_scale4, y_prior_df4, y_global_prior_scale[4],
        yGlobal4, yLocal4, yOol4, yMix4, yAux4, family[4],
        y_slab_scale[4], y_caux4);
      yBeta4[lt_idx[4]] = abs(yBeta4[lt_idx[4]]);
    }
  } 
  if (M > 4) {
    if (has_aux[5] == 1) {
      yAux5[1] = make_aux(yAux5_unscaled[1], y_prior_dist_for_aux[5],
        y_prior_mean_for_aux[5], y_prior_scale_for_aux[5]);
      if (yAux5[1] > yAuxMaximum)
        yAuxMaximum = yAux5[1];
    }
    if (yK[5] > 0){
      yBeta5 = make_beta(z_yBeta5, y_prior_dist[5], y_prior_mean5,
        y_prior_scale5, y_prior_df5, y_global_prior_scale[5],
        yGlobal5, yLocal5, yOol5, yMix5, yAux5, family[5],
        y_slab_scale[5], y_caux5);
      yBeta5[lt_idx[5]] = abs(yBeta5[lt_idx[5]]);
    }
  } 
  if (M > 5) {
    if (has_aux[6] == 1) {
      yAux6[1] = make_aux(yAux6_unscaled[1], y_prior_dist_for_aux[6],
        y_prior_mean_for_aux[6], y_prior_scale_for_aux[6]);
      if (yAux6[1] > yAuxMaximum)
        yAuxMaximum = yAux6[1];
    }
    if (yK[6] > 0){
      yBeta6 = make_beta(z_yBeta6, y_prior_dist[6], y_prior_mean6,
        y_prior_scale6, y_prior_df6, y_global_prior_scale[6],
        yGlobal6, yLocal6, yOol6, yMix6, yAux6, family[6],
        y_slab_scale[6], y_caux6);
      yBeta6[lt_idx[6]] = abs(yBeta6[lt_idx[6]]);
    }
  } 
  if (M > 6) {
    if (has_aux[7] == 1) {
      yAux7[1] = make_aux(yAux7_unscaled[1], y_prior_dist_for_aux[7],
        y_prior_mean_for_aux[7], y_prior_scale_for_aux[7]);
      if (yAux7[1] > yAuxMaximum)
        yAuxMaximum = yAux7[1];
    }
    if (yK[7] > 0){
      yBeta7 = make_beta(z_yBeta7, y_prior_dist[7], y_prior_mean7,
        y_prior_scale7, y_prior_df7, y_global_prior_scale[7],
        yGlobal7, yLocal7, yOol7, yMix7, yAux7, family[7],
        y_slab_scale[7], y_caux7);
      yBeta7[lt_idx[7]] = abs(yBeta7[lt_idx[7]]);
    }
  } 
  if (M > 7) {
    if (has_aux[8] == 1) {
      yAux8[1] = make_aux(yAux8_unscaled[1], y_prior_dist_for_aux[8],
        y_prior_mean_for_aux[8], y_prior_scale_for_aux[8]);
      if (yAux8[1] > yAuxMaximum)
        yAuxMaximum = yAux8[1];
    }
    if (yK[8] > 0){
      yBeta8 = make_beta(z_yBeta8, y_prior_dist[8], y_prior_mean8,
        y_prior_scale8, y_prior_df8, y_global_prior_scale[8],
        yGlobal8, yLocal8, yOol8, yMix8, yAux8, family[8],
        y_slab_scale[8], y_caux8);
      yBeta8[lt_idx[8]] = abs(yBeta8[lt_idx[8]]);
    }
  } 
  if (M > 8) {
    if (has_aux[9] == 1) {
      yAux9[1] = make_aux(yAux9_unscaled[1], y_prior_dist_for_aux[9],
        y_prior_mean_for_aux[9], y_prior_scale_for_aux[9]);
      if (yAux9[1] > yAuxMaximum)
        yAuxMaximum = yAux9[1];
    }
    if (yK[9] > 0){
      yBeta9 = make_beta(z_yBeta9, y_prior_dist[9], y_prior_mean9,
        y_prior_scale9, y_prior_df9, y_global_prior_scale[9],
        yGlobal9, yLocal9, yOol9, yMix9, yAux9, family[9],
        y_slab_scale[9], y_caux9);
      yBeta9[lt_idx[9]] = abs(yBeta9[lt_idx[9]]);
    }
  } 
  if (M > 9) {
    if (has_aux[10] == 1) {
      yAux10[1] = make_aux(yAux10_unscaled[1], y_prior_dist_for_aux[10],
        y_prior_mean_for_aux[10], y_prior_scale_for_aux[10]);
      if (yAux10[1] > yAuxMaximum)
        yAuxMaximum = yAux10[1];
    }
    if (yK[10] > 0){
      yBeta10 = make_beta(z_yBeta10, y_prior_dist[10], y_prior_mean10,
        y_prior_scale10, y_prior_df10, y_global_prior_scale[10],
        yGlobal10, yLocal10, yOol10, yMix10, yAux10, family[10],
        y_slab_scale[10], y_caux10);
      yBeta10[lt_idx[10]] = abs(yBeta10[lt_idx[10]]);
    }
  } 
  if (M > 10) {
    if (has_aux[11] == 1) {
      yAux11[1] = make_aux(yAux11_unscaled[1], y_prior_dist_for_aux[11],
        y_prior_mean_for_aux[11], y_prior_scale_for_aux[11]);
      if (yAux11[1] > yAuxMaximum)
        yAuxMaximum = yAux11[1];
    }
    if (yK[11] > 0){
      yBeta11 = make_beta(z_yBeta11, y_prior_dist[11], y_prior_mean11,
        y_prior_scale11, y_prior_df11, y_global_prior_scale[11],
        yGlobal11, yLocal11, yOol11, yMix11, yAux11, family[11],
        y_slab_scale[11], y_caux11);
      yBeta11[lt_idx[11]] = abs(yBeta11[lt_idx[11]]);
    }
  } 
  if (M > 11) {
    if (has_aux[12] == 1) {
      yAux12[1] = make_aux(yAux12_unscaled[1], y_prior_dist_for_aux[12],
        y_prior_mean_for_aux[12], y_prior_scale_for_aux[12]);
      if (yAux12[1] > yAuxMaximum)
        yAuxMaximum = yAux12[1];
    }
    if (yK[12] > 0){
      yBeta12 = make_beta(z_yBeta12, y_prior_dist[12], y_prior_mean12,
        y_prior_scale12, y_prior_df12, y_global_prior_scale[12],
        yGlobal12, yLocal12, yOol12, yMix12, yAux12, family[12],
        y_slab_scale[12], y_caux12);
      yBeta12[lt_idx[12]] = abs(yBeta12[lt_idx[12]]);
    }
  } 
  if (M > 12) {
    if (has_aux[13] == 1) {
      yAux13[1] = make_aux(yAux13_unscaled[1], y_prior_dist_for_aux[13],
        y_prior_mean_for_aux[13], y_prior_scale_for_aux[13]);
      if (yAux13[1] > yAuxMaximum)
        yAuxMaximum = yAux13[1];
    }
    if (yK[13] > 0){
      yBeta13 = make_beta(z_yBeta13, y_prior_dist[13], y_prior_mean13,
        y_prior_scale13, y_prior_df13, y_global_prior_scale[13],
        yGlobal13, yLocal13, yOol13, yMix13, yAux13, family[13],
        y_slab_scale[13], y_caux13);
      yBeta13[lt_idx[13]] = abs(yBeta13[lt_idx[13]]);
    }
  } 
  if (M > 13) {
    if (has_aux[14] == 1) {
      yAux14[1] = make_aux(yAux14_unscaled[1], y_prior_dist_for_aux[14],
        y_prior_mean_for_aux[14], y_prior_scale_for_aux[14]);
      if (yAux14[1] > yAuxMaximum)
        yAuxMaximum = yAux14[1];
    }
    if (yK[14] > 0){
      yBeta14 = make_beta(z_yBeta14, y_prior_dist[14], y_prior_mean14,
        y_prior_scale14, y_prior_df14, y_global_prior_scale[14],
        yGlobal14, yLocal14, yOol14, yMix14, yAux14, family[14],
        y_slab_scale[14], y_caux14);
      yBeta14[lt_idx[14]] = abs(yBeta14[lt_idx[14]]);
    }
  } 
  if (M > 14) {
    if (has_aux[15] == 1) {
      yAux15[1] = make_aux(yAux15_unscaled[1], y_prior_dist_for_aux[15],
        y_prior_mean_for_aux[15], y_prior_scale_for_aux[15]);
      if (yAux15[1] > yAuxMaximum)
        yAuxMaximum = yAux15[1];
    }
    if (yK[15] > 0){
      yBeta15 = make_beta(z_yBeta15, y_prior_dist[15], y_prior_mean15,
        y_prior_scale15, y_prior_df15, y_global_prior_scale[15],
        yGlobal15, yLocal15, yOol15, yMix15, yAux15, family[15],
        y_slab_scale[15], y_caux15);
      yBeta15[lt_idx[15]] = abs(yBeta15[lt_idx[15]]);
    }
  } 
  if (M > 15) {
    if (has_aux[16] == 1) {
      yAux16[1] = make_aux(yAux16_unscaled[1], y_prior_dist_for_aux[16],
        y_prior_mean_for_aux[16], y_prior_scale_for_aux[16]);
      if (yAux16[1] > yAuxMaximum)
        yAuxMaximum = yAux16[1];
    }
    if (yK[16] > 0){
      yBeta16 = make_beta(z_yBeta16, y_prior_dist[16], y_prior_mean16,
        y_prior_scale16, y_prior_df16, y_global_prior_scale[16],
        yGlobal16, yLocal16, yOol16, yMix16, yAux16, family[16],
        y_slab_scale[16], y_caux16);
      yBeta16[lt_idx[16]] = abs(yBeta16[lt_idx[16]]);
    }
  } 
  if (M > 16) {
    if (has_aux[17] == 1) {
      yAux17[1] = make_aux(yAux17_unscaled[1], y_prior_dist_for_aux[17],
        y_prior_mean_for_aux[17], y_prior_scale_for_aux[17]);
      if (yAux17[1] > yAuxMaximum)
        yAuxMaximum = yAux17[1];
    }
    if (yK[17] > 0){
      yBeta17 = make_beta(z_yBeta17, y_prior_dist[17], y_prior_mean17,
        y_prior_scale17, y_prior_df17, y_global_prior_scale[17],
        yGlobal17, yLocal17, yOol17, yMix17, yAux17, family[17],
        y_slab_scale[17], y_caux17);
      yBeta17[lt_idx[17]] = abs(yBeta17[lt_idx[17]]);
    }
  } 
  if (M > 17) {
    if (has_aux[18] == 1) {
      yAux18[1] = make_aux(yAux18_unscaled[1], y_prior_dist_for_aux[18],
        y_prior_mean_for_aux[18], y_prior_scale_for_aux[18]);
      if (yAux18[1] > yAuxMaximum)
        yAuxMaximum = yAux18[1];
    }
    if (yK[18] > 0){
      yBeta18 = make_beta(z_yBeta18, y_prior_dist[18], y_prior_mean18,
        y_prior_scale18, y_prior_df18, y_global_prior_scale[18],
        yGlobal18, yLocal18, yOol18, yMix18, yAux18, family[18],
        y_slab_scale[18], y_caux18);
      yBeta18[lt_idx[18]] = abs(yBeta18[lt_idx[18]]);
    }
  } 
  if (M > 18) {
    if (has_aux[19] == 1) {
      yAux19[1] = make_aux(yAux19_unscaled[1], y_prior_dist_for_aux[19],
        y_prior_mean_for_aux[19], y_prior_scale_for_aux[19]);
      if (yAux19[1] > yAuxMaximum)
        yAuxMaximum = yAux19[1];
    }
    if (yK[19] > 0){
      yBeta19 = make_beta(z_yBeta19, y_prior_dist[19], y_prior_mean19,
        y_prior_scale19, y_prior_df19, y_global_prior_scale[19],
        yGlobal19, yLocal19, yOol19, yMix19, yAux19, family[19],
        y_slab_scale[19], y_caux19);
      yBeta19[lt_idx[19]] = abs(yBeta19[lt_idx[19]]);
    }
  } 
  if (M > 19) {
    if (has_aux[20] == 1) {
      yAux20[1] = make_aux(yAux20_unscaled[1], y_prior_dist_for_aux[20],
        y_prior_mean_for_aux[20], y_prior_scale_for_aux[20]);
      if (yAux20[1] > yAuxMaximum)
        yAuxMaximum = yAux20[1];
    }
    if (yK[20] > 0){
      yBeta20 = make_beta(z_yBeta20, y_prior_dist[20], y_prior_mean20,
        y_prior_scale20, y_prior_df20, y_global_prior_scale[20],
        yGlobal20, yLocal20, yOol20, yMix20, yAux20, family[20],
        y_slab_scale[20], y_caux20);
      yBeta20[lt_idx[20]] = abs(yBeta20[lt_idx[20]]);
    }
  }

  // group level params, under decov prior
  if (prior_dist_for_cov == 1) {
    int mark = 1;
    // cov matrix
    theta_L = make_theta_L(len_theta_L, p, yAuxMaximum, tau,
                           b_prior_scale, zeta, rho, z_T);
    // group-level params for first grouping factor
    if (bK1 > 0)
      bMat1 = make_b_matrix(z_b, theta_L, p, l, 1);
    // group level params for second grouping factor
    if (bK2 > 0)
      bMat2 = make_b_matrix(z_b, theta_L, p, l, 2);
  }

  // group-level params, under lkj prior
  else if (prior_dist_for_cov == 2) {
    // group-level params for first grouping factor
    if (bK1 == 1)
      bMat1 = (bSd1[1] * z_bMat1)';
    else if (bK1 > 1)
      bMat1 = (diag_pre_multiply(bSd1, bCholesky1) * z_bMat1)';
    // group level params for second grouping factor
    if (bK2 == 1)
      bMat2 = (bSd2[1] * z_bMat2)';
    else if (bK2 > 1)
      bMat2 = (diag_pre_multiply(bSd2, bCholesky2) * z_bMat2)';
  }

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