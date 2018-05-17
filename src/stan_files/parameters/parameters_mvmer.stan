  // intercepts
  real<lower=lb(intercept_type[1]),upper=ub(intercept_type[1])>
    yGamma1[intercept_type[1] > 0];
  real<lower=lb(intercept_type[2]),upper=ub(intercept_type[2])>
    yGamma2[intercept_type[2] > 0];
  real<lower=lb(intercept_type[3]),upper=ub(intercept_type[3])>
    yGamma3[intercept_type[3] > 0];
  real<lower=lb(intercept_type[4]),upper=ub(intercept_type[4])>
    yGamma4[intercept_type[4] > 0];
  real<lower=lb(intercept_type[5]),upper=ub(intercept_type[5])>
    yGamma5[intercept_type[5] > 0];
  real<lower=lb(intercept_type[6]),upper=ub(intercept_type[6])>
    yGamma6[intercept_type[6] > 0];
  real<lower=lb(intercept_type[7]),upper=ub(intercept_type[7])>
    yGamma7[intercept_type[7] > 0];
  real<lower=lb(intercept_type[8]),upper=ub(intercept_type[8])>
    yGamma8[intercept_type[8] > 0];
  real<lower=lb(intercept_type[9]),upper=ub(intercept_type[9])>
    yGamma9[intercept_type[9] > 0];
  real<lower=lb(intercept_type[10]),upper=ub(intercept_type[10])>
    yGamma10[intercept_type[10] > 0];
  real<lower=lb(intercept_type[11]),upper=ub(intercept_type[11])>
    yGamma11[intercept_type[11] > 0];
  real<lower=lb(intercept_type[12]),upper=ub(intercept_type[12])>
    yGamma12[intercept_type[12] > 0];
  real<lower=lb(intercept_type[13]),upper=ub(intercept_type[13])>
    yGamma13[intercept_type[13] > 0];
  real<lower=lb(intercept_type[14]),upper=ub(intercept_type[14])>
    yGamma14[intercept_type[14] > 0];
  real<lower=lb(intercept_type[15]),upper=ub(intercept_type[15])>
    yGamma15[intercept_type[15] > 0];
  real<lower=lb(intercept_type[16]),upper=ub(intercept_type[16])>
    yGamma16[intercept_type[16] > 0];
  real<lower=lb(intercept_type[17]),upper=ub(intercept_type[17])>
    yGamma17[intercept_type[17] > 0];
  real<lower=lb(intercept_type[18]),upper=ub(intercept_type[18])>
    yGamma18[intercept_type[18] > 0];
  real<lower=lb(intercept_type[19]),upper=ub(intercept_type[19])>
    yGamma19[intercept_type[19] > 0];
  real<lower=lb(intercept_type[20]),upper=ub(intercept_type[20])>
    yGamma20[intercept_type[20] > 0];
  // population level primitive params
  vector[yK[1]] z_yBeta1;
  vector[yK[2]] z_yBeta2;
  vector[yK[3]] z_yBeta3;
  vector[yK[4]] z_yBeta4;
  vector[yK[5]] z_yBeta5;
  vector[yK[6]] z_yBeta6;
  vector[yK[7]] z_yBeta7;
  vector[yK[8]] z_yBeta8;
  vector[yK[9]] z_yBeta9;
  vector[yK[10]] z_yBeta10;
  vector[yK[11]] z_yBeta11;
  vector[yK[12]] z_yBeta12;
  vector[yK[13]] z_yBeta13;
  vector[yK[14]] z_yBeta14;
  vector[yK[15]] z_yBeta15;
  vector[yK[16]] z_yBeta16;
  vector[yK[17]] z_yBeta17;
  vector[yK[18]] z_yBeta18;
  vector[yK[19]] z_yBeta19;
  vector[yK[20]] z_yBeta20;

  // group level params, decov prior
  vector[prior_dist_for_cov == 1 ? q : 0] z_b;
  vector[prior_dist_for_cov == 1 ? len_z_T : 0] z_T;
  vector<lower=0,upper=1>[prior_dist_for_cov == 1 ? len_rho : 0] rho;
  vector<lower=0>[prior_dist_for_cov == 1 ? len_concentration : 0] zeta;
  vector<lower=0>[prior_dist_for_cov == 1 ? t : 0] tau;

  // group level params for first grouping factor
    // group-level sds
    vector<lower=0>[prior_dist_for_cov == 2 ? bK1 : 0] bSd1;
    // unscaled group-level params
    matrix[prior_dist_for_cov == 2 && bK1 >  0 ? bK1 : 0, bK1 >  0 ? bN1 : 0] z_bMat1;
    // cholesky factor of corr matrix (if > 1 random effect)
    cholesky_factor_corr[prior_dist_for_cov == 2 && bK1 > 1 ? bK1 : 0] bCholesky1;

  // group level params for second grouping factor
    // group-level sds
    vector<lower=0>[prior_dist_for_cov == 2 ? bK2 : 0] bSd2;
    // unscaled group-level params
    matrix[prior_dist_for_cov == 2 && bK2 >  0 ? bK2 : 0, bK2 >  0 ? bN2 : 0] z_bMat2;
    // cholesky factor of corr matrix (if > 1 random effect)
    cholesky_factor_corr[prior_dist_for_cov == 2 && bK2 > 1 ? bK2 : 0] bCholesky2;

  // auxiliary params, interpretation depends on family
  real<lower=0> yAux1_unscaled[has_aux[1]];
  real<lower=0> yAux2_unscaled[has_aux[2]];
  real<lower=0> yAux3_unscaled[has_aux[3]];
  real<lower=0> yAux4_unscaled[has_aux[4]];
  real<lower=0> yAux5_unscaled[has_aux[5]];
  real<lower=0> yAux6_unscaled[has_aux[6]];
  real<lower=0> yAux7_unscaled[has_aux[7]];
  real<lower=0> yAux8_unscaled[has_aux[8]];
  real<lower=0> yAux9_unscaled[has_aux[9]];
  real<lower=0> yAux10_unscaled[has_aux[10]];
  real<lower=0> yAux11_unscaled[has_aux[11]];
  real<lower=0> yAux12_unscaled[has_aux[12]];
  real<lower=0> yAux13_unscaled[has_aux[13]];
  real<lower=0> yAux14_unscaled[has_aux[14]];
  real<lower=0> yAux15_unscaled[has_aux[15]];
  real<lower=0> yAux16_unscaled[has_aux[16]];
  real<lower=0> yAux17_unscaled[has_aux[17]];
  real<lower=0> yAux18_unscaled[has_aux[18]];
  real<lower=0> yAux19_unscaled[has_aux[19]];
  real<lower=0> yAux20_unscaled[has_aux[20]];
  
  // params for priors
  real<lower=0> yGlobal1[yHs1];
  real<lower=0> yGlobal2[yHs2];
  real<lower=0> yGlobal3[yHs3];
  real<lower=0> yGlobal4[yHs4];
  real<lower=0> yGlobal5[yHs5];
  real<lower=0> yGlobal6[yHs6];
  real<lower=0> yGlobal7[yHs7];
  real<lower=0> yGlobal8[yHs8];
  real<lower=0> yGlobal9[yHs9];
  real<lower=0> yGlobal10[yHs10];
  real<lower=0> yGlobal11[yHs11];
  real<lower=0> yGlobal12[yHs12];
  real<lower=0> yGlobal13[yHs13];
  real<lower=0> yGlobal14[yHs14];
  real<lower=0> yGlobal15[yHs15];
  real<lower=0> yGlobal16[yHs16];
  real<lower=0> yGlobal17[yHs17];
  real<lower=0> yGlobal18[yHs18];
  real<lower=0> yGlobal19[yHs19];
  real<lower=0> yGlobal20[yHs20];
  vector<lower=0>[yK[1]] yLocal1[yHs1];
  vector<lower=0>[yK[2]] yLocal2[yHs2];
  vector<lower=0>[yK[3]] yLocal3[yHs3];
  vector<lower=0>[yK[4]] yLocal4[yHs4];
  vector<lower=0>[yK[5]] yLocal5[yHs5];
  vector<lower=0>[yK[6]] yLocal6[yHs6];
  vector<lower=0>[yK[7]] yLocal7[yHs7];
  vector<lower=0>[yK[8]] yLocal8[yHs8];
  vector<lower=0>[yK[9]] yLocal9[yHs9];
  vector<lower=0>[yK[10]] yLocal10[yHs10];
  vector<lower=0>[yK[11]] yLocal11[yHs11];
  vector<lower=0>[yK[12]] yLocal12[yHs12];
  vector<lower=0>[yK[13]] yLocal13[yHs13];
  vector<lower=0>[yK[14]] yLocal14[yHs14];
  vector<lower=0>[yK[15]] yLocal15[yHs15];
  vector<lower=0>[yK[16]] yLocal16[yHs16];
  vector<lower=0>[yK[17]] yLocal17[yHs17];
  vector<lower=0>[yK[18]] yLocal18[yHs18];
  vector<lower=0>[yK[19]] yLocal19[yHs19];
  vector<lower=0>[yK[20]] yLocal20[yHs20];
  real<lower=0> y_caux1[yHs1 > 0];
  real<lower=0> y_caux2[yHs2 > 0];
  real<lower=0> y_caux3[yHs3 > 0];
  real<lower=0> y_caux4[yHs4 > 0];
  real<lower=0> y_caux5[yHs5 > 0];
  real<lower=0> y_caux6[yHs6 > 0];
  real<lower=0> y_caux7[yHs7 > 0];
  real<lower=0> y_caux8[yHs8 > 0];
  real<lower=0> y_caux9[yHs9 > 0];
  real<lower=0> y_caux10[yHs10 > 0];
  real<lower=0> y_caux11[yHs11 > 0];
  real<lower=0> y_caux12[yHs12 > 0];
  real<lower=0> y_caux13[yHs13 > 0];
  real<lower=0> y_caux14[yHs14 > 0];
  real<lower=0> y_caux15[yHs15 > 0];
  real<lower=0> y_caux16[yHs16 > 0];
  real<lower=0> y_caux17[yHs17 > 0];
  real<lower=0> y_caux18[yHs18 > 0];
  real<lower=0> y_caux19[yHs19 > 0];
  real<lower=0> y_caux20[yHs20 > 0];
  real<lower=0> yOol1[y_prior_dist[1] == 6]; // one_over_lambda
  real<lower=0> yOol2[y_prior_dist[2] == 6];
  real<lower=0> yOol3[y_prior_dist[3] == 6];
  real<lower=0> yOol4[y_prior_dist[4] == 6];
  real<lower=0> yOol5[y_prior_dist[5] == 6];
  real<lower=0> yOol6[y_prior_dist[6] == 6];
  real<lower=0> yOol7[y_prior_dist[7] == 6];
  real<lower=0> yOol8[y_prior_dist[8] == 6];
  real<lower=0> yOol9[y_prior_dist[9] == 6];
  real<lower=0> yOol10[y_prior_dist[10] == 6];
  real<lower=0> yOol11[y_prior_dist[11] == 6];
  real<lower=0> yOol12[y_prior_dist[12] == 6];
  real<lower=0> yOol13[y_prior_dist[13] == 6];
  real<lower=0> yOol14[y_prior_dist[14] == 6];
  real<lower=0> yOol15[y_prior_dist[15] == 6];
  real<lower=0> yOol16[y_prior_dist[16] == 6];
  real<lower=0> yOol17[y_prior_dist[17] == 6];
  real<lower=0> yOol18[y_prior_dist[18] == 6];
  real<lower=0> yOol19[y_prior_dist[19] == 6];
  real<lower=0> yOol20[y_prior_dist[20] == 6];
  vector<lower=0>[yK[1]] yMix1[y_prior_dist[1] == 5 || y_prior_dist[1] == 6];
  vector<lower=0>[yK[2]] yMix2[y_prior_dist[2] == 5 || y_prior_dist[2] == 6];
  vector<lower=0>[yK[3]] yMix3[y_prior_dist[3] == 5 || y_prior_dist[3] == 6];
  vector<lower=0>[yK[4]] yMix4[y_prior_dist[4] == 5 || y_prior_dist[4] == 6];
  vector<lower=0>[yK[5]] yMix5[y_prior_dist[5] == 5 || y_prior_dist[5] == 6];
  vector<lower=0>[yK[6]] yMix6[y_prior_dist[6] == 5 || y_prior_dist[6] == 6];
  vector<lower=0>[yK[7]] yMix7[y_prior_dist[7] == 5 || y_prior_dist[7] == 6];
  vector<lower=0>[yK[8]] yMix8[y_prior_dist[8] == 5 || y_prior_dist[8] == 6];
  vector<lower=0>[yK[9]] yMix9[y_prior_dist[9] == 5 || y_prior_dist[9] == 6];
  vector<lower=0>[yK[10]] yMix10[y_prior_dist[10] == 5 || y_prior_dist[10] == 6];
  vector<lower=0>[yK[11]] yMix11[y_prior_dist[11] == 5 || y_prior_dist[11] == 6];
  vector<lower=0>[yK[12]] yMix12[y_prior_dist[12] == 5 || y_prior_dist[12] == 6];
  vector<lower=0>[yK[13]] yMix13[y_prior_dist[13] == 5 || y_prior_dist[13] == 6];
  vector<lower=0>[yK[14]] yMix14[y_prior_dist[14] == 5 || y_prior_dist[14] == 6];
  vector<lower=0>[yK[15]] yMix15[y_prior_dist[15] == 5 || y_prior_dist[15] == 6];
  vector<lower=0>[yK[16]] yMix16[y_prior_dist[16] == 5 || y_prior_dist[16] == 6];
  vector<lower=0>[yK[17]] yMix17[y_prior_dist[17] == 5 || y_prior_dist[17] == 6];
  vector<lower=0>[yK[18]] yMix18[y_prior_dist[18] == 5 || y_prior_dist[18] == 6];
  vector<lower=0>[yK[19]] yMix19[y_prior_dist[19] == 5 || y_prior_dist[19] == 6];
  vector<lower=0>[yK[20]] yMix20[y_prior_dist[20] == 5 || y_prior_dist[20] == 6];