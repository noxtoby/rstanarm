  // population level dimensions
  int<lower=1,upper=20> M; // num submodels with data (limit of 20)
  int<lower=0,upper=1> has_aux[20]; // has auxiliary param
  int<lower=0,upper=1> has_weights; // has observation weights
  int<lower=0,upper=2> resp_type[20]; // 1=real,2=integer,0=none
  int<lower=0,upper=20> intercept_type[20]; // 1=unbounded,2=lob,20=upb,0=none
  int<lower=0> yNobs[20]; // num observations
  int<lower=0> yNeta[20]; // required length of eta
  int<lower=0> yK[20]; // num predictors

  // group level dimensions, for decov prior
  int<lower=0> t;    // num. terms (maybe 0) with a | in the glmer formula
  int<lower=1> p[t]; // num. variables on the LHS of each |
  int<lower=1> l[t]; // num. levels for the factor(s) on the RHS of each |
  int<lower=0> q;    // conceptually equals \sum_{i=1}^t p_i \times l_i
  int<lower=0> len_theta_L; // length of the theta_L vector

  // group level dimensions, for lkj prior

    // group factor 1
    int<lower=0> bN1; // num groups
    int<lower=0> bK1; // total num params
    int<lower=0> bK1_len[20]; // num params in each submodel
    int<lower=0> bK1_idx[20,2]; // beg/end index for group params

    // group factor 2
    int<lower=0> bN2; // num groups
    int<lower=0> bK2; // total num params
    int<lower=0> bK2_len[20]; // num params in each submodel
    int<lower=0> bK2_idx[20,2]; // beg/end index for group params

  int<lower=0> lt_idx[M]; // index for latent time term of X matrix