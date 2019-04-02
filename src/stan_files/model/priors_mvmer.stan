  // Log-priors, auxiliary params
  if (has_aux[1] == 1)
    aux_lp(yAux1_unscaled[1], y_prior_dist_for_aux[1],
           y_prior_scale_for_aux[1], y_prior_df_for_aux[1]);
  if (M > 1 && has_aux[2] == 1)
    aux_lp(yAux2_unscaled[1], y_prior_dist_for_aux[2],
           y_prior_scale_for_aux[2], y_prior_df_for_aux[2]);
  if (M > 2 && has_aux[3] == 1)
    aux_lp(yAux3_unscaled[1], y_prior_dist_for_aux[3],
           y_prior_scale_for_aux[3], y_prior_df_for_aux[3]);
  if (M > 3 && has_aux[4] == 1)
    aux_lp(yAux4_unscaled[1], y_prior_dist_for_aux[4],
           y_prior_scale_for_aux[4], y_prior_df_for_aux[4]);
  if (M > 4 && has_aux[5] == 1)
    aux_lp(yAux5_unscaled[1], y_prior_dist_for_aux[5],
           y_prior_scale_for_aux[5], y_prior_df_for_aux[5]);
  if (M > 5 && has_aux[6] == 1)
    aux_lp(yAux6_unscaled[1], y_prior_dist_for_aux[6],
           y_prior_scale_for_aux[6], y_prior_df_for_aux[6]);
  if (M > 6 && has_aux[7] == 1)
    aux_lp(yAux7_unscaled[1], y_prior_dist_for_aux[7],
           y_prior_scale_for_aux[7], y_prior_df_for_aux[7]);
  if (M > 7 && has_aux[8] == 1)
    aux_lp(yAux8_unscaled[1], y_prior_dist_for_aux[8],
           y_prior_scale_for_aux[8], y_prior_df_for_aux[8]);
  if (M > 8 && has_aux[9] == 1)
    aux_lp(yAux9_unscaled[1], y_prior_dist_for_aux[9],
           y_prior_scale_for_aux[9], y_prior_df_for_aux[9]);
  if (M > 9 && has_aux[10] == 1)
    aux_lp(yAux10_unscaled[1], y_prior_dist_for_aux[10],
           y_prior_scale_for_aux[10], y_prior_df_for_aux[10]);
  if (M > 10 && has_aux[11] == 1)
    aux_lp(yAux11_unscaled[1], y_prior_dist_for_aux[11],
           y_prior_scale_for_aux[11], y_prior_df_for_aux[11]);
  if (M > 11 && has_aux[12] == 1)
    aux_lp(yAux12_unscaled[1], y_prior_dist_for_aux[12],
           y_prior_scale_for_aux[12], y_prior_df_for_aux[12]);
  if (M > 12 && has_aux[13] == 1)
    aux_lp(yAux13_unscaled[1], y_prior_dist_for_aux[13],
           y_prior_scale_for_aux[13], y_prior_df_for_aux[13]);
  if (M > 13 && has_aux[14] == 1)
    aux_lp(yAux14_unscaled[1], y_prior_dist_for_aux[14],
           y_prior_scale_for_aux[14], y_prior_df_for_aux[14]);
  if (M > 14 && has_aux[15] == 1)
    aux_lp(yAux15_unscaled[1], y_prior_dist_for_aux[15],
           y_prior_scale_for_aux[15], y_prior_df_for_aux[15]);
  if (M > 15 && has_aux[16] == 1)
    aux_lp(yAux16_unscaled[1], y_prior_dist_for_aux[16],
           y_prior_scale_for_aux[16], y_prior_df_for_aux[16]);
  if (M > 16 && has_aux[17] == 1)
    aux_lp(yAux17_unscaled[1], y_prior_dist_for_aux[17],
           y_prior_scale_for_aux[17], y_prior_df_for_aux[17]);
  if (M > 17 && has_aux[18] == 1)
    aux_lp(yAux18_unscaled[1], y_prior_dist_for_aux[18],
           y_prior_scale_for_aux[18], y_prior_df_for_aux[18]);
  if (M > 18 && has_aux[19] == 1)
    aux_lp(yAux19_unscaled[1], y_prior_dist_for_aux[19],
           y_prior_scale_for_aux[19], y_prior_df_for_aux[19]);
  if (M > 19 && has_aux[20] == 1)
    aux_lp(yAux20_unscaled[1], y_prior_dist_for_aux[20],
           y_prior_scale_for_aux[20], y_prior_df_for_aux[20]);
  // Log priors, intercepts
  if (intercept_type[1] > 0)
    gamma_lp(yGamma1[1], y_prior_dist_for_intercept[1], y_prior_mean_for_intercept[1],
             y_prior_scale_for_intercept[1], y_prior_df_for_intercept[1]);
  if (M > 1 && intercept_type[2] > 0)
    gamma_lp(yGamma2[1], y_prior_dist_for_intercept[2], y_prior_mean_for_intercept[2],
             y_prior_scale_for_intercept[2], y_prior_df_for_intercept[2]);
  if (M > 2 && intercept_type[3] > 0)
    gamma_lp(yGamma3[1], y_prior_dist_for_intercept[3], y_prior_mean_for_intercept[3],
             y_prior_scale_for_intercept[3], y_prior_df_for_intercept[3]);
  if (M > 3 && intercept_type[4] > 0)
    gamma_lp(yGamma4[1], y_prior_dist_for_intercept[4], y_prior_mean_for_intercept[4],
             y_prior_scale_for_intercept[4], y_prior_df_for_intercept[4]);
  if (M > 4 && intercept_type[5] > 0)
    gamma_lp(yGamma5[1], y_prior_dist_for_intercept[5], y_prior_mean_for_intercept[5],
             y_prior_scale_for_intercept[5], y_prior_df_for_intercept[5]);
  if (M > 5 && intercept_type[6] > 0)
    gamma_lp(yGamma6[1], y_prior_dist_for_intercept[6], y_prior_mean_for_intercept[6],
             y_prior_scale_for_intercept[6], y_prior_df_for_intercept[6]);
  if (M > 6 && intercept_type[7] > 0)
    gamma_lp(yGamma7[1], y_prior_dist_for_intercept[7], y_prior_mean_for_intercept[7],
             y_prior_scale_for_intercept[7], y_prior_df_for_intercept[7]);
  if (M > 7 && intercept_type[8] > 0)
    gamma_lp(yGamma8[1], y_prior_dist_for_intercept[8], y_prior_mean_for_intercept[8],
             y_prior_scale_for_intercept[8], y_prior_df_for_intercept[8]);
  if (M > 8 && intercept_type[9] > 0)
    gamma_lp(yGamma9[1], y_prior_dist_for_intercept[9], y_prior_mean_for_intercept[9],
             y_prior_scale_for_intercept[9], y_prior_df_for_intercept[9]);
  if (M > 9 && intercept_type[10] > 0)
    gamma_lp(yGamma10[1], y_prior_dist_for_intercept[10], y_prior_mean_for_intercept[10],
             y_prior_scale_for_intercept[10], y_prior_df_for_intercept[10]);
  if (M > 10 && intercept_type[11] > 0)
    gamma_lp(yGamma11[1], y_prior_dist_for_intercept[11], y_prior_mean_for_intercept[11],
             y_prior_scale_for_intercept[11], y_prior_df_for_intercept[11]);
  if (M > 11 && intercept_type[12] > 0)
    gamma_lp(yGamma12[1], y_prior_dist_for_intercept[12], y_prior_mean_for_intercept[12],
             y_prior_scale_for_intercept[12], y_prior_df_for_intercept[12]);
  if (M > 12 && intercept_type[13] > 0)
    gamma_lp(yGamma13[1], y_prior_dist_for_intercept[13], y_prior_mean_for_intercept[13],
             y_prior_scale_for_intercept[13], y_prior_df_for_intercept[13]);
  if (M > 13 && intercept_type[14] > 0)
    gamma_lp(yGamma14[1], y_prior_dist_for_intercept[14], y_prior_mean_for_intercept[14],
             y_prior_scale_for_intercept[14], y_prior_df_for_intercept[14]);
  if (M > 14 && intercept_type[15] > 0)
    gamma_lp(yGamma15[1], y_prior_dist_for_intercept[15], y_prior_mean_for_intercept[15],
             y_prior_scale_for_intercept[15], y_prior_df_for_intercept[15]);
  if (M > 15 && intercept_type[16] > 0)
    gamma_lp(yGamma16[1], y_prior_dist_for_intercept[16], y_prior_mean_for_intercept[16],
             y_prior_scale_for_intercept[16], y_prior_df_for_intercept[16]);
  if (M > 16 && intercept_type[17] > 0)
    gamma_lp(yGamma17[1], y_prior_dist_for_intercept[17], y_prior_mean_for_intercept[17],
             y_prior_scale_for_intercept[17], y_prior_df_for_intercept[17]);
  if (M > 17 && intercept_type[18] > 0)
    gamma_lp(yGamma18[1], y_prior_dist_for_intercept[18], y_prior_mean_for_intercept[18],
             y_prior_scale_for_intercept[18], y_prior_df_for_intercept[18]);
  if (M > 18 && intercept_type[19] > 0)
    gamma_lp(yGamma19[1], y_prior_dist_for_intercept[19], y_prior_mean_for_intercept[19],
             y_prior_scale_for_intercept[19], y_prior_df_for_intercept[19]);
  if (M > 19 && intercept_type[20] > 0)
    gamma_lp(yGamma20[1], y_prior_dist_for_intercept[20], y_prior_mean_for_intercept[20],
             y_prior_scale_for_intercept[20], y_prior_df_for_intercept[20]);
  
  // Log priors, population level params
  if (yK[1] > 0)
    beta_lp(z_yBeta1, y_prior_dist[1], y_prior_scale1, y_prior_df1,
            y_global_prior_df[1], yLocal1, yGlobal1, yMix1, yOol1,
            y_slab_df[1], y_caux1);
  if (M > 1 && yK[2] > 0)
    beta_lp(z_yBeta2, y_prior_dist[2], y_prior_scale2, y_prior_df2,
            y_global_prior_df[2], yLocal2, yGlobal2, yMix2, yOol2,
            y_slab_df[2], y_caux2);
  if (M > 2 && yK[3] > 0)
    beta_lp(z_yBeta3, y_prior_dist[3], y_prior_scale3, y_prior_df3,
            y_global_prior_df[3], yLocal3, yGlobal3, yMix3, yOol3,
            y_slab_df[3], y_caux3);
  if (M > 3 && yK[4] > 0)
    beta_lp(z_yBeta4, y_prior_dist[4], y_prior_scale4, y_prior_df4,
            y_global_prior_df[4], yLocal4, yGlobal4, yMix4, yOol4,
            y_slab_df[4], y_caux4);
  if (M > 4 && yK[5] > 0)
    beta_lp(z_yBeta5, y_prior_dist[5], y_prior_scale5, y_prior_df5,
            y_global_prior_df[5], yLocal5, yGlobal5, yMix5, yOol5,
            y_slab_df[5], y_caux5);
  if (M > 5 && yK[6] > 0)
    beta_lp(z_yBeta6, y_prior_dist[6], y_prior_scale6, y_prior_df6,
            y_global_prior_df[6], yLocal6, yGlobal6, yMix6, yOol6,
            y_slab_df[6], y_caux6);
  if (M > 6 && yK[7] > 0)
    beta_lp(z_yBeta7, y_prior_dist[7], y_prior_scale7, y_prior_df7,
            y_global_prior_df[7], yLocal7, yGlobal7, yMix7, yOol7,
            y_slab_df[7], y_caux7);
  if (M > 7 && yK[8] > 0)
    beta_lp(z_yBeta8, y_prior_dist[8], y_prior_scale8, y_prior_df8,
            y_global_prior_df[8], yLocal8, yGlobal8, yMix8, yOol8,
            y_slab_df[8], y_caux8);
  if (M > 8 && yK[9] > 0)
    beta_lp(z_yBeta9, y_prior_dist[9], y_prior_scale9, y_prior_df9,
            y_global_prior_df[9], yLocal9, yGlobal9, yMix9, yOol9,
            y_slab_df[9], y_caux9);
  if (M > 9 && yK[10] > 0)
    beta_lp(z_yBeta10, y_prior_dist[10], y_prior_scale10, y_prior_df10,
            y_global_prior_df[10], yLocal10, yGlobal10, yMix10, yOol10,
            y_slab_df[10], y_caux10);
  if (M > 10 && yK[11] > 0)
    beta_lp(z_yBeta11, y_prior_dist[11], y_prior_scale11, y_prior_df11,
            y_global_prior_df[11], yLocal11, yGlobal11, yMix11, yOol11,
            y_slab_df[11], y_caux11);
  if (M > 11 && yK[12] > 0)
    beta_lp(z_yBeta12, y_prior_dist[12], y_prior_scale12, y_prior_df12,
            y_global_prior_df[12], yLocal12, yGlobal12, yMix12, yOol12,
            y_slab_df[12], y_caux12);
  if (M > 12 && yK[13] > 0)
    beta_lp(z_yBeta13, y_prior_dist[13], y_prior_scale13, y_prior_df13,
            y_global_prior_df[13], yLocal13, yGlobal13, yMix13, yOol13,
            y_slab_df[13], y_caux13);
  if (M > 13 && yK[14] > 0)
    beta_lp(z_yBeta14, y_prior_dist[14], y_prior_scale14, y_prior_df14,
            y_global_prior_df[14], yLocal14, yGlobal14, yMix14, yOol14,
            y_slab_df[14], y_caux14);
  if (M > 14 && yK[15] > 0)
    beta_lp(z_yBeta15, y_prior_dist[15], y_prior_scale15, y_prior_df15,
            y_global_prior_df[15], yLocal15, yGlobal15, yMix15, yOol15,
            y_slab_df[15], y_caux15);
  if (M > 15 && yK[16] > 0)
    beta_lp(z_yBeta16, y_prior_dist[16], y_prior_scale16, y_prior_df16,
            y_global_prior_df[16], yLocal16, yGlobal16, yMix16, yOol16,
            y_slab_df[16], y_caux16);
  if (M > 16 && yK[17] > 0)
    beta_lp(z_yBeta17, y_prior_dist[17], y_prior_scale17, y_prior_df17,
            y_global_prior_df[17], yLocal17, yGlobal17, yMix17, yOol17,
            y_slab_df[17], y_caux17);
  if (M > 17 && yK[18] > 0)
    beta_lp(z_yBeta18, y_prior_dist[18], y_prior_scale18, y_prior_df18,
            y_global_prior_df[18], yLocal18, yGlobal18, yMix18, yOol18,
            y_slab_df[18], y_caux18);
  if (M > 18 && yK[19] > 0)
    beta_lp(z_yBeta19, y_prior_dist[19], y_prior_scale19, y_prior_df19,
            y_global_prior_df[19], yLocal19, yGlobal19, yMix19, yOol19,
            y_slab_df[19], y_caux19);
  if (M > 19 && yK[20] > 0)
    beta_lp(z_yBeta20, y_prior_dist[20], y_prior_scale20, y_prior_df20,
            y_global_prior_df[20], yLocal20, yGlobal20, yMix20, yOol20,
            y_slab_df[20], y_caux20);
            
  // Log priors, group level terms
  if (prior_dist_for_cov == 1) { // decov
    real dummy = decov_lp(z_b, z_T, rho, zeta, tau, b_prior_regularization,
                          delta, b_prior_shape, t, p);
  }
  else if (prior_dist_for_cov == 2) { // lkj
    if (bK1 > 0) {
      // sds for group factor 1
      target += student_t_lpdf(bSd1 | b1_prior_df, 0, b1_prior_scale);
      // primitive coefs for group factor 1
      target += normal_lpdf(to_vector(z_bMat1) | 0, 1);
      // corr matrix for group factor 1
      if (bK1 > 1)
        target += lkj_corr_cholesky_lpdf(bCholesky1 | b1_prior_regularization);
    }
    if (bK2 > 0) {
      // sds for group factor 2
      target += student_t_lpdf(bSd2 | b2_prior_df, 0, b2_prior_scale);
      // primitive coefs for group factor 2
      target += normal_lpdf(to_vector(z_bMat2) | 0, 1);
      // corr matrix for group factor 2
      if (bK2 > 1)
        target += lkj_corr_cholesky_lpdf(bCholesky2 | b2_prior_regularization);
    }
  }
