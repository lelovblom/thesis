// PARTS OF THIS CODE FROM THE RSTANARM PACKAGE, MVMER FUNCTION, BY SAM BRILLEMAN
#include RSTANARM_STAN_FILES/pre/Columbia_copyright.stan
#include RSTANARM_STAN_FILES/pre/Brilleman_copyright.stan
#include RSTANARM_STAN_FILES/pre/license.stan

// Multivariate GLM with correlated group-specific terms
functions {
#include RSTANARM_STAN_FILES/functions/common_functions.stan
#include RSTANARM_STAN_FILES/functions/bernoulli_likelihoods.stan
#include RSTANARM_STAN_FILES/functions/binomial_likelihoods.stan
#include RSTANARM_STAN_FILES/functions/continuous_likelihoods.stan
#include RSTANARM_STAN_FILES/functions/count_likelihoods.stan
#include RSTANARM_STAN_FILES/functions/mvmer_functions.stan

  // CREATE FUNCTIONS TO EVALUATE THE HAZARDS FOR THE ITERATED INTEGRALS
  real evaluate_h12 (vector bspline_basis,  // A VECTOR OF 6 POINTS (INCLUDING THE INTERCEPT)
                    vector bspline_deriv_basis,
                    real x,                 // the time point itself (not the spline expansion)
                    vector beta_basehaz,    // regression coefficient for the baseline hazards
                    vector beta_eta,        // regression coefficients for the fixed covariates
                    real X1, real X2, real X3, // fixed covariates to be used in the different components
                                               // X1 = A1c, X2 = randomization group, X3 = intervention group
                    vector alpha,        // regression coefficient(s) for the longitudinal outcome - meant for "alpha12"
                    real beta_y1_int, real beta_y2_int, real beta_y3_int,       // longitudinal submodel intercepts (from the MLMM)
                    vector beta_y1,        // longitudinal submodel parameters
                    vector beta_y2,                            
                    vector beta_y3,
                    vector bi,
                    vector y1_xbar, vector y2_xbar, vector y3_xbar 
                    ) {
    real h12;
    h12 = exp(beta_basehaz[1] + dot_product(beta_basehaz[2:6],bspline_basis[2:6]) +  
              beta_eta[1]*X1 + beta_eta[2]*X2 +
              alpha[1]*(beta_y1_int - dot_product(y1_xbar, beta_y1) + // ETDRS CURRENT VALUE
                        dot_product(beta_y1[1:5],bspline_basis[2:6]) +  
                        beta_y1[6] * X3 +
                        X2 * (dot_product(beta_y1[7:11],bspline_basis[2:6])) + 
                        bi[1] + dot_product(bi[2:6],bspline_basis[2:6]) 
                        ) +
              alpha[2]*( dot_product(beta_y1[1:5],bspline_deriv_basis[2:6]) + // ETDRS CURRENT SLOPE
                         X2 * (dot_product(beta_y1[7:11],bspline_deriv_basis[2:6])) +
                         dot_product(bi[2:6],bspline_deriv_basis[2:6])
                       ) +
              alpha[3]*(beta_y2_int - dot_product(y2_xbar, beta_y2) + // AER CURRENT VALUE
                        dot_product(beta_y2[1:5],bspline_basis[2:6]) + 
                        beta_y2[6] * X3 +
                        X2 * (dot_product(beta_y2[7:11],bspline_basis[2:6])) + 
                        bi[7] + dot_product(bi[8:12],bspline_basis[2:6]) 
                        ) +
              alpha[4]*( dot_product(beta_y2[1:5],bspline_deriv_basis[2:6]) + // AER CURRENT SLOPE
                         X2 * (dot_product(beta_y2[7:11],bspline_deriv_basis[2:6])) +
                         dot_product(bi[8:12],bspline_deriv_basis[2:6])
                       ) +
              alpha[5]*(beta_y3[1] + beta_y3[2]*X2 + // EGFR CURRENT SLOPE
                        bi[14] + 2*bi[15]*x
                        )
          );
    return h12;
  }
  real evaluate_h23 (vector bspline_basis,  // A VECTOR OF 6 POINTS (INCLUDING THE INTERCEPT)
                    real x,                 // the time point itself (not the spline expansion)
                    vector beta_basehaz,    // regression coefficient for the baseline hazards
                    vector beta_eta,        // regression coefficients for the fixed covariates
                    real X1, real X2, real X3, // fixed covariates to be used in the different components
                                               // X1 = A1c, X2 = randomization group, X3 = intervention group
                    vector alpha,        // regression coefficient(s) for the longitudinal outcome - meant for "alpha12"
                    real beta_y1_int, real beta_y2_int, real beta_y3_int,       // longitudinal submodel intercepts (from the MLMM)
                    vector beta_y1,        // longitudinal submodel parameters
                    vector beta_y2,                            
                    vector beta_y3,
                    vector bi,
                    vector y1_xbar, vector y2_xbar, vector y3_xbar 
                    ) {
    real h23;
    h23 = exp(beta_basehaz[1] + dot_product(beta_basehaz[2:6],bspline_basis[2:6]) + 
              beta_eta[1]*X1 + //beta_eta[2]*X2 +
              alpha[1]*(beta_y1_int - dot_product(y1_xbar, beta_y1) + // ETDRS CURRENT VALUE
                        dot_product(beta_y1[1:5],bspline_basis[2:6]) + 
                        beta_y1[6] * X3 +
                        X2 * (dot_product(beta_y1[7:11],bspline_basis[2:6])) + 
                        bi[1] + dot_product(bi[2:6],bspline_basis[2:6])  
                        ) +
              alpha[2]*(beta_y2_int - dot_product(y2_xbar, beta_y2) + // AER CURRENT VALUE
                        dot_product(beta_y2[1:5],bspline_basis[2:6]) + 
                        beta_y2[6] * X3 +
                        X2 * (dot_product(beta_y2[7:11],bspline_basis[2:6])) + 
                        bi[7] + dot_product(bi[8:12],bspline_basis[2:6]) 
                        ) +
              alpha[3]*(beta_y3[1] + beta_y3[2]*X2 + // EGFR CURRENT SLOPE
                        bi[14] + 2*bi[15]*x
                        )
          );
    return h23;
  }
}
data {
#include RSTANARM_STAN_FILES/data/dimensions_mvmer.stan
#include RSTANARM_STAN_FILES/data/data_mvmer.stan

  // EVENT DATA
  // start with parameter info
  int<lower=0> e_K[2];           // num. of predictors in event submodel, [1] is 1->2 transition, [2] is 2->3 transition
  int<lower=0> basehaz_df;    // "df" for baseline hazard, so 6 for each 1->2 and 2->3
  int<lower=0> a12_K;           // num. of association parameters
  int<lower=0> a23_K;
  
  //type 1
  int<lower=0> type1_N_GK;
  int<lower=1, upper=bN1> type1_id[type1_N_GK]; //b_N replaced by bN1
  vector[type1_N_GK] type1_qpoints;//quadrature points
  vector[type1_N_GK] type1_wpoints;//quadrature weights
  int<lower=0> nrow_e_Xq_type1;     // num. rows in event submodel predictor matrix
  matrix[nrow_e_Xq_type1,e_K[1]] e_Xq_type1; // predictor matrix (event submodel)
  matrix[nrow_e_Xq_type1,basehaz_df] basehaz_X_type1; // design matrix for baseline hazard
  int<lower=0> nrow_y_Xq_type1; // num. rows in long. predictor matrix at quadpoints
  int<lower=0> y_Zq_id_type1[nrow_y_Xq_type1]; // group indexing for re design matrix
  matrix[nrow_y_Xq_type1,yK[1]] y1_Xq_type1; // fe design matrix at quadpoints      // NB: Y1 AND Y2 SHARE THE SAME MATRICES
  vector[nrow_y_Xq_type1] y1_Zq_type1[bK1_len[1]]; // re design matrix at quadpoints
  matrix[nrow_y_Xq_type1,yK[1]] y1_Xq_slope_type1; // fe design matrix for slope at quadpoints (with zeroes inserted appropriately, e.g. for retbase)
  vector[nrow_y_Xq_type1] y1_Zq_slope_type1[bK1_len[1]]; // re design matrix at quadpoints (with zeroes inserted appropriately, e.g. for intercept)
  matrix[nrow_y_Xq_type1,yK[3]] y3_Xq_slope_type1; // fe design matrix FOR EGFR SLOPE at quadpoints      
  vector[nrow_y_Xq_type1] y3_Zq_slope_type1[bK1_len[3]]; // re design matrix FOR EGFR SLOPE at quadpoints
  
  //type 2
  int<lower=0> type2_N;
  int<lower=1, upper=bN1> type2_id[type2_N];
  vector[type2_N] type2_HBAEL;
  vector[type2_N] type2_GROUP2;
  vector[type2_N] type2_RETBASE2;
  real x_scaled_b_type2[6,15,type2_N];
  real x2_scaled_b_type2[6,15,15,type2_N];
  real x3_scaled_b_type2[6,15,15,type2_N];
  real x_scaled_type2[15,type2_N];
  real x2_scaled_type2[15,15,type2_N];
  real x3_scaled_type2[15,15,type2_N];
  real wt_scaled_type2[15,type2_N];
  real wt2_scaled_type2[15,15,type2_N];
  real wt3_scaled_type2[15,15,type2_N];
  real drv_x_scaled_b_type2[6,15,type2_N];
  real drv_x2_scaled_b_type2[6,15,15,type2_N];
  real drv_x3_scaled_b_type2[6,15,15,type2_N];

  //type 3
  int<lower=0> type3_N_GK;
  int<lower=1, upper=bN1> type3_id[type3_N_GK];
  vector[type3_N_GK] type3_qpoints;// PRE-GENERATED QUADRATURE POINTS AND WEIGHTS
  vector[type3_N_GK] type3_wpoints;
  int<lower=0> nrow_e_Xq_type3;     // num. rows in event submodel predictor matrix
  matrix[nrow_e_Xq_type3,e_K[2]] e_Xq_type3; // predictor matrix (event submodel)  // CAREFUL, e_K[2] NOT e_K[1]
  matrix[nrow_e_Xq_type3,basehaz_df] basehaz_X_type3; // design matrix for baseline hazard
  int<lower=0> nrow_y_Xq_type3; // num. rows in long. predictor matrix at quadpoints
  int<lower=0> y_Zq_id_type3[nrow_y_Xq_type3]; // group indexing for re design matrix
  matrix[nrow_y_Xq_type3,yK[1]] y1_Xq_type3; // fe design matrix at quadpoints      // NB: Y1 AND Y2 SHARE THE SAME MATRICES
  vector[nrow_y_Xq_type3] y1_Zq_type3[bK1_len[1]]; // re design matrix at quadpoints
  matrix[nrow_y_Xq_type3,yK[3]] y3_Xq_slope_type3; // fe design matrix FOR EGFR SLOPE at quadpoints      
  vector[nrow_y_Xq_type3] y3_Zq_slope_type3[bK1_len[3]]; // re design matrix FOR EGFR SLOPE at quadpoints
  
  //type 4
  int<lower=0> type4_N;
  int<lower=1, upper=bN1> type4_id[type4_N];
  vector[type4_N] type4_HBAEL;
  vector[type4_N] type4_GROUP2;
  vector[type4_N] type4_RETBASE2;
  real x_scaled_b_type4[6,15,type4_N];
  real x2_scaled_b_type4[6,15,15,type4_N];
  real x3_scaled_b_type4[6,15,15,type4_N];
  real x_scaled_type4[15,type4_N];
  real x2_scaled_type4[15,15,type4_N];
  real x3_scaled_type4[15,15,type4_N];
  real wt_scaled_type4[15,type4_N];
  real wt2_scaled_type4[15,15,type4_N];
  real wt3_scaled_type4[15,15,type4_N];
  real drv_x_scaled_b_type4[6,15,type4_N];
  real drv_x2_scaled_b_type4[6,15,15,type4_N];
  real drv_x3_scaled_b_type4[6,15,15,type4_N];
  
  //type 5
  int<lower=0> type5_N_GK;
  int<lower=1, upper=bN1> type5_id[type5_N_GK];
  vector[type5_N_GK] type5_qpoints;// PRE-GENERATED QUADRATURE POINTS AND WEIGHTS
  vector[type5_N_GK] type5_wpoints;
  int<lower=0> nrow_e_Xq_type5;     // num. rows in event submodel predictor matrix
  matrix[nrow_e_Xq_type5,e_K[2]] e_Xq_type5; // predictor matrix (event submodel)
  matrix[nrow_e_Xq_type5,basehaz_df] basehaz_X_type5; // design matrix for baseline hazard
  int<lower=0> nrow_y_Xq_type5; // num. rows in long. predictor matrix at quadpoints
  int<lower=0> y_Zq_id_type5[nrow_y_Xq_type5]; // group indexing for re design matrix
  matrix[nrow_y_Xq_type5,yK[1]] y1_Xq_type5; // fe design matrix at quadpoints      // NB: Y1 AND Y2 SHARE THE SAME MATRICES
  vector[nrow_y_Xq_type5] y1_Zq_type5[bK1_len[1]]; // re design matrix at quadpoints
  matrix[nrow_y_Xq_type5,yK[3]] y3_Xq_slope_type5; // fe design matrix FOR EGFR SLOPE at quadpoints      
  vector[nrow_y_Xq_type5] y3_Zq_slope_type5[bK1_len[3]]; // re design matrix FOR EGFR SLOPE at quadpoints
  
  //type 6
  int<lower=0> type6_N;
  int<lower=1, upper=bN1> type6_id[type6_N];
  vector[type6_N] type6_HBAEL;
  vector[type6_N] type6_GROUP2;
  vector[type6_N] type6_RETBASE2;
  real x_scaled_b_type6[6,15,type6_N];
  real x2_scaled_b_type6[6,15,15,type6_N];
  real x3_scaled_b_type6[6,15,15,type6_N];
  real x_scaled_type6[15,type6_N];
  real x2_scaled_type6[15,15,type6_N];
  real x3_scaled_type6[15,15,type6_N];
  real wt_scaled_type6[15,type6_N];
  real wt2_scaled_type6[15,15,type6_N];
  real wt3_scaled_type6[15,15,type6_N];
  real drv_x_scaled_b_type6[6,15,type6_N];
  real drv_x2_scaled_b_type6[6,15,15,type6_N];
  real drv_x3_scaled_b_type6[6,15,15,type6_N];
  
  //type 7
  int<lower=0> type7_N_GK;
  int<lower=1, upper=bN1> type7_id[type7_N_GK];
  vector[type7_N_GK] type7_qpoints;// PRE-GENERATED QUADRATURE POINTS AND WEIGHTS
  vector[type7_N_GK] type7_wpoints;
  int<lower=0> nrow_e_Xq_type7;     // num. rows in event submodel predictor matrix
  matrix[nrow_e_Xq_type7,e_K[2]] e_Xq_type7; // predictor matrix (event submodel)
  matrix[nrow_e_Xq_type7,basehaz_df] basehaz_X_type7; // design matrix for baseline hazard
  int<lower=0> nrow_y_Xq_type7; // num. rows in long. predictor matrix at quadpoints
  int<lower=0> y_Zq_id_type7[nrow_y_Xq_type7]; // group indexing for re design matrix
  matrix[nrow_y_Xq_type7,yK[1]] y1_Xq_type7; // fe design matrix at quadpoints      // NB: Y1 AND Y2 SHARE THE SAME MATRICES
  vector[nrow_y_Xq_type7] y1_Zq_type7[bK1_len[1]]; // re design matrix at quadpoints
  matrix[nrow_y_Xq_type7,yK[3]] y3_Xq_slope_type7; // fe design matrix FOR EGFR SLOPE at quadpoints      
  vector[nrow_y_Xq_type7] y3_Zq_slope_type7[bK1_len[3]]; // re design matrix FOR EGFR SLOPE at quadpoints

  //events
  int<lower=0> events_N;
  vector[events_N] events_year;
  int<lower=1, upper=bN1> events_id[events_N];
  int<lower=0> nrow_e_Xq_events;     // num. rows in event submodel predictor matrix
  matrix[nrow_e_Xq_events,e_K[2]] e_Xq_events; // predictor matrix (event submodel)
  matrix[nrow_e_Xq_events,basehaz_df] basehaz_X_events; // design matrix for baseline hazard
  int<lower=0> nrow_y_Xq_events; // num. rows in long. predictor matrix at quadpoints
  int<lower=0> y_Zq_id_events[nrow_y_Xq_events]; // group indexing for re design matrix
  matrix[nrow_y_Xq_events,yK[1]] y1_Xq_events; // fe design matrix at quadpoints      // NB: Y1 AND Y2 SHARE THE SAME MATRICES
  vector[nrow_y_Xq_events] y1_Zq_events[bK1_len[1]]; // re design matrix at quadpoints
  matrix[nrow_y_Xq_events,yK[3]] y3_Xq_slope_events; // fe design matrix FOR EGFR SLOPE at quadpoints      
  vector[nrow_y_Xq_events] y3_Zq_slope_events[bK1_len[3]]; // re design matrix FOR EGFR SLOPE at quadpoints

  // NEW DATA TYPES FOR THE "EVALUATE_ETA" FUNCTION 
  vector[0] y1_offset_eta;
  vector[0] y2_offset_eta;
  vector[0] y3_offset_eta;
  vector[0] y1_z2q_eta[0];
  vector[0] y2_z2q_eta[0];
  vector[0] y3_z2q_eta[0];
  int<lower=0> y1_z2q_id_eta[0];
  int<lower=0> y2_z2q_id_eta[0];
  int<lower=0> y3_z2q_id_eta[0];

  // NOW THE HYPERPARAMETERS
#include RSTANARM_STAN_FILES/data/hyperparameters_mvmer.stan

  vector[e_K[1]]          e12_prior_mean;
  vector[e_K[2]]          e23_prior_mean;
  vector[a12_K]          a12_prior_mean;
  vector[a23_K]          a23_prior_mean;
  vector[basehaz_df]   e12_prior_mean_for_aux;
  vector[basehaz_df]   e23_prior_mean_for_aux;
  vector<lower=0>[e_K[1]]    e12_prior_scale;
  vector<lower=0>[e_K[2]]    e23_prior_scale;
  vector<lower=0>[a12_K]    a12_prior_scale;
  vector<lower=0>[a23_K]    a23_prior_scale;
  vector<lower=0>[basehaz_df] e12_prior_scale_for_aux;
  vector<lower=0>[basehaz_df] e23_prior_scale_for_aux;
}
transformed data {
#include RSTANARM_STAN_FILES/tdata/tdata_mvmer.stan
}
parameters {
#include RSTANARM_STAN_FILES/parameters/parameters_mvmer.stan

  // NOW FOR THE EVENT MODEL
  vector[e_K[1]] e12_z_beta; // primitive coefs in event submodel (log hazard ratios)
  vector[e_K[2]] e23_z_beta;
  vector[a12_K] a12_z_beta; // primitive assoc params (log hazard ratios)
  vector[a23_K] a23_z_beta;
  vector[basehaz_df] e12_aux_unscaled; // unscaled coefs for baseline hazard   
  vector[basehaz_df] e23_aux_unscaled;
}
transformed parameters { 
  // DECLARE FOR THE EVENT MODEL
  vector[e_K[1]] e12_beta;       // coefs in event submodel (log hazard ratios)
  vector[e_K[2]] e23_beta;
  vector[a12_K] a12_beta;       // assoc params in event submodel (log hazard ratios) 
  vector[a23_K] a23_beta;
  vector[basehaz_df] e12_aux; // b-spline coefs for baseline hazard
  vector[basehaz_df] e23_aux;

#include RSTANARM_STAN_FILES/tparameters/tparameters_mvmer.stan

  // coefs for event submodel (incl. association parameters)
  e12_beta = e12_z_beta .* e12_prior_scale + e12_prior_mean;
  e23_beta = e23_z_beta .* e23_prior_scale + e23_prior_mean;
  a12_beta = a12_z_beta .* a12_prior_scale + a12_prior_mean;
  a23_beta = a23_z_beta .* a23_prior_scale + a23_prior_mean;

  // b-spline coefs for baseline hazard
  e12_aux = e12_aux_unscaled .* e12_prior_scale_for_aux + e12_prior_mean_for_aux;
  e23_aux = e23_aux_unscaled .* e23_prior_scale_for_aux + e23_prior_mean_for_aux;
}
model {
  // Log likelihoods
  // increments target with mvmer log liks
#include RSTANARM_STAN_FILES/model/mvmer_lp.stan 

// EVENT SUBMODELS
  { // Type 1 transitions - PRE-DEFINED QUADRATURE POINTS
    vector[nrow_e_Xq_type1] eta;
    vector[nrow_y_Xq_type1] y1_eta_q;
    vector[nrow_y_Xq_type1] y2_eta_q;
    vector[nrow_y_Xq_type1] y1_etaslope_q;
    vector[nrow_y_Xq_type1] y2_etaslope_q;
    vector[nrow_y_Xq_type1] y3_etaslope_q;

    y1_eta_q = evaluate_eta(y1_Xq_type1, y1_Zq_type1, y1_z2q_eta, //y1_xq_type1, y1_z1q_type1, 
                                    y_Zq_id_type1, y1_z2q_id_eta, //y1_z1q_id_type1 these are all the same for my application (across the 3 outcomes)
                                    yGamma1, yBeta1, bMat1, bMat2,
                                    0, 0, 1, //bMat1_colshift, bMat2_colshift, 0,    //NB: the last of these is intercept type, which we have as 1 (unbounded)
                                              y1_offset_eta);
    y1_etaslope_q = evaluate_eta(y1_Xq_slope_type1, y1_Zq_slope_type1, y1_z2q_eta,  
                                y_Zq_id_type1, y1_z2q_id_eta, 
                                yGamma1, yBeta1, bMat1, bMat2,
                                0, 0, 0, //bMat1_colshift, bMat2_colshift, 0,    //NB: SWITCH INTERCEPT TO ZERO SINCE WE DON'T NEED IT FOR SLOPE
                                          y1_offset_eta);
    y2_eta_q = evaluate_eta(y1_Xq_type1, y1_Zq_type1, y2_z2q_eta, // Y1 AND Y2 HAVE SAME FE AND RE DEESIGN MATRICES
                                    y_Zq_id_type1, y2_z2q_id_eta,
                                    yGamma2, yBeta2, bMat1, bMat2,
                                    6, 0, 1, //bMat1_colshift, bMat2_colshift, 0,
                                              y2_offset_eta);
    y2_etaslope_q = evaluate_eta(y1_Xq_slope_type1, y1_Zq_slope_type1, y1_z2q_eta,  
                                y_Zq_id_type1, y1_z2q_id_eta, 
                                yGamma2, yBeta2, bMat1, bMat2,
                                6, 0, 0, //bMat1_colshift, bMat2_colshift, 0,    //NB: SWITCH INTERCEPT TO ZERO SINCE WE DON'T NEED IT FOR SLOPE
                                          y2_offset_eta);
    y3_etaslope_q = evaluate_eta(y3_Xq_slope_type1, y3_Zq_slope_type1, y3_z2q_eta,
                                    y_Zq_id_type1, y3_z2q_id_eta,
                                    yGamma3, yBeta3, bMat1, bMat2,
                                    12, 0, 0, //bMat1_colshift, bMat2_colshift, 0, //NB: SWITCH INTERCEPT OFF
                                              y3_offset_eta);
    
    eta = basehaz_X_type1 * e12_aux + e_Xq_type1 * e12_beta 
          + a12_beta[1] * y1_eta_q 
          + a12_beta[2] * y1_etaslope_q
          + a12_beta[3] * y2_eta_q 
          + a12_beta[4] * y2_etaslope_q 
          + a12_beta[5] * y3_etaslope_q 
          - a12_beta[1] * dot_product(yXbar1, yBeta1) 
          - a12_beta[3] * dot_product(yXbar2, yBeta2) 
          ;
    target += - dot_product(type1_wpoints,exp(eta));
  }
  
  { // Type 3 transitions - PRE-DEFINED QUADRATURE POINTS
    vector[nrow_e_Xq_type3] eta;
    vector[nrow_y_Xq_type3] y1_eta_q;
    vector[nrow_y_Xq_type3] y2_eta_q;
    vector[nrow_y_Xq_type3] y3_etaslope_q;
    
    y1_eta_q = evaluate_eta(y1_Xq_type3, y1_Zq_type3, y1_z2q_eta, //y1_xq_type3, y1_z1q_type3, 
                                    y_Zq_id_type3, y1_z2q_id_eta, //y1_z1q_id_type3 these are all the same for my application (across the 3 outcomes)
                                    yGamma1, yBeta1, bMat1, bMat2,
                                    0, 0, 1, //bMat1_colshift, bMat2_colshift, 0, //intercept_tpye should be 1!
                                              y1_offset_eta);
    y2_eta_q = evaluate_eta(y1_Xq_type3, y1_Zq_type3, y2_z2q_eta,  // Y1 AND Y2 HAVE SAME FE AND RE DEESIGN MATRICES
                                    y_Zq_id_type3, y2_z2q_id_eta,
                                    yGamma2, yBeta2, bMat1, bMat2,
                                    6, 0, 1, //bMat1_colshift, bMat2_colshift, 0,
                                              y2_offset_eta);
    y3_etaslope_q = evaluate_eta(y3_Xq_slope_type3, y3_Zq_slope_type3, y3_z2q_eta,
                                    y_Zq_id_type3, y3_z2q_id_eta,
                                    yGamma3, yBeta3, bMat1, bMat2,
                                    12, 0, 0, //bMat1_colshift, bMat2_colshift, 0, //NB: SWITCH INTERCEPT OFF
                                              y3_offset_eta);
    
    eta = basehaz_X_type3 * e23_aux + e_Xq_type3 * e23_beta 
          + a23_beta[1] * y1_eta_q 
          + a23_beta[2] * y2_eta_q /*+ a23_beta[3] * y3_eta_q*/
          - a23_beta[1] * dot_product(yXbar1, yBeta1) 
          - a23_beta[2] * dot_product(yXbar2, yBeta2) /*- a23_beta[3] * dot_product(yXbar3, yBeta3)*/
          + a23_beta[3] * y3_etaslope_q
          ;
    target += - dot_product(type3_wpoints,exp(eta));
  }
  
  { // Type 5 transitions - PRE-DEFINED QUADRATURE POINTS
    vector[nrow_e_Xq_type5] eta;
    vector[nrow_y_Xq_type5] y1_eta_q;
    vector[nrow_y_Xq_type5] y2_eta_q;
    //vector[nrow_y_Xq_type5] y3_eta_q;
    vector[nrow_y_Xq_type5] y3_etaslope_q;

    y1_eta_q = evaluate_eta(y1_Xq_type5, y1_Zq_type5, y1_z2q_eta, //y1_xq_type5, y1_z1q_type5, 
                                    y_Zq_id_type5, y1_z2q_id_eta, //y1_z1q_id_type5 these are all the same for my application (across the 3 outcomes)
                                    yGamma1, yBeta1, bMat1, bMat2,
                                    0, 0, 1, //bMat1_colshift, bMat2_colshift, 0, //NB: the last of these is intercept type, which we have as 1 (unbounded)
                                              y1_offset_eta);
    y2_eta_q = evaluate_eta(y1_Xq_type5, y1_Zq_type5, y2_z2q_eta,
                                    y_Zq_id_type5, y2_z2q_id_eta,
                                    yGamma2, yBeta2, bMat1, bMat2,
                                    6, 0, 1, //bMat1_colshift, bMat2_colshift, 0,
                                              y2_offset_eta);
    y3_etaslope_q = evaluate_eta(y3_Xq_slope_type5, y3_Zq_slope_type5, y3_z2q_eta,
                                    y_Zq_id_type5, y3_z2q_id_eta,
                                    yGamma3, yBeta3, bMat1, bMat2,
                                    12, 0, 0, //bMat1_colshift, bMat2_colshift, 0, //NB: SWITCH INTERCEPT OFF
                                              y3_offset_eta);
    
    eta = basehaz_X_type5 * e23_aux + e_Xq_type5 * e23_beta 
          + a23_beta[1] * y1_eta_q 
          + a23_beta[2] * y2_eta_q 
          - a23_beta[1] * dot_product(yXbar1, yBeta1) 
          - a23_beta[2] * dot_product(yXbar2, yBeta2) 
          + a23_beta[3] * y3_etaslope_q
          ;
    target += - dot_product(type5_wpoints,exp(eta));
  }

  { // Type 7 transitions - PRE-DEFINED QUADRATURE POINTS, AND ONLY THE TRANSITION PROBABILITY NOT THE EVENT TIME
    vector[nrow_e_Xq_type7] eta;
    vector[nrow_y_Xq_type7] y1_eta_q;
    vector[nrow_y_Xq_type7] y2_eta_q;
    vector[nrow_y_Xq_type7] y3_etaslope_q;
    
    y1_eta_q = evaluate_eta(y1_Xq_type7, y1_Zq_type7, y1_z2q_eta, //y1_xq_type7, y1_z1q_type7, 
                                    y_Zq_id_type7, y1_z2q_id_eta, //y1_z1q_id_type7 these are all the same for my application (across the 3 outcomes)
                                    yGamma1, yBeta1, bMat1, bMat2,
                                    0, 0, 1, //bMat1_colshift, bMat2_colshift, 0,
                                              y1_offset_eta);
    y2_eta_q = evaluate_eta(y1_Xq_type7, y1_Zq_type7, y2_z2q_eta,    // Y1 AND Y2 HAVE SAME FE AND RE DEESIGN MATRICES
                                    y_Zq_id_type7, y2_z2q_id_eta,
                                    yGamma2, yBeta2, bMat1, bMat2,
                                    6, 0, 1, //bMat1_colshift, bMat2_colshift, 0,
                                              y2_offset_eta);
    y3_etaslope_q = evaluate_eta(y3_Xq_slope_type7, y3_Zq_slope_type7, y3_z2q_eta,
                                    y_Zq_id_type7, y3_z2q_id_eta,
                                    yGamma3, yBeta3, bMat1, bMat2,
                                    12, 0, 0, //bMat1_colshift, bMat2_colshift, 0, //NB: SWITCH INTERCEPT OFF
                                              y3_offset_eta);
    
    eta = basehaz_X_type7 * e23_aux + e_Xq_type7 * e23_beta 
          + a23_beta[1] * y1_eta_q 
          + a23_beta[2] * y2_eta_q
          - a23_beta[1] * dot_product(yXbar1, yBeta1) 
          - a23_beta[2] * dot_product(yXbar2, yBeta2)
          + a23_beta[3] * y3_etaslope_q
          ;
    target += - dot_product(type7_wpoints,exp(eta));
  }
  
  { // OBSERVED STATE 3 EVENTS
    vector[nrow_e_Xq_events] log_hazard;
    vector[nrow_y_Xq_events] y1_eta_q;
    vector[nrow_y_Xq_events] y2_eta_q;
    vector[nrow_y_Xq_events] y3_etaslope_q;

    y1_eta_q = evaluate_eta(y1_Xq_events, y1_Zq_events, y1_z2q_eta, //y1_xq_events, y1_z1q_events, 
                                    y_Zq_id_events, y1_z2q_id_eta, //y1_z1q_id_events these are all the same for my application (across the 3 outcomes)
                                    yGamma1, yBeta1, bMat1, bMat2,
                                    0, 0, 1, //bMat1_colshift, bMat2_colshift, 0,
                                              y1_offset_eta);
    y2_eta_q = evaluate_eta(y1_Xq_events, y1_Zq_events, y2_z2q_eta,    // Y1 AND Y2 HAVE SAME FE AND RE DEESIGN MATRICES
                                    y_Zq_id_events, y2_z2q_id_eta,
                                    yGamma2, yBeta2, bMat1, bMat2,
                                    6, 0, 1, //bMat1_colshift, bMat2_colshift, 0,
                                              y2_offset_eta);
    y3_etaslope_q = evaluate_eta(y3_Xq_slope_events, y3_Zq_slope_events, y3_z2q_eta,
                                    y_Zq_id_events, y3_z2q_id_eta,
                                    yGamma3, yBeta3, bMat1, bMat2,
                                    12, 0, 0, //bMat1_colshift, bMat2_colshift, 0, //NB: SWITCH INTERCEPT OFF
                                              y3_offset_eta);
    
    log_hazard = basehaz_X_events * e23_aux + e_Xq_events * e23_beta 
                  + a23_beta[1] * y1_eta_q 
                  + a23_beta[2] * y2_eta_q /*+ a23_beta[3] * y3_eta_q*/
                  - a23_beta[1] * dot_product(yXbar1, yBeta1) 
                  - a23_beta[2] * dot_product(yXbar2, yBeta2) /*- a23_beta[3] * dot_product(yXbar3, yBeta3)*/
                  + a23_beta[3] * y3_etaslope_q
                  ;
    target += sum(log_hazard);
  }
  
  { // TYPE 2 TRANSITIONS, USING THE PRE-BUILT SPLINES MATRICES, AND REQUIRING LOOPS
    vector[6] bspline_basis;
    vector[6] bspline_basis2;
    vector[6] bspline_basis3;
    vector[6] bspline_drv_basis;
    vector[6] bspline_drv_basis2;
    vector[6] bspline_drv_basis3;
    real x;
    real x2;
    real x3;
    real wt;
    real wt2;
    real wt3;
    real HBAEL;
    real GROUP2;
    real RETBASE2;
    vector[15] bi;
    
    real gauss_total;
    real q12;
    real gaussH1;
    real H1;
    real gaussH2;
    real H2;
    real piece;
    real weighted_piece;//NEW
    real p12;
    
    for (n in 1:type2_N) {
      HBAEL = type2_HBAEL[n];
      GROUP2 = type2_GROUP2[n];
      RETBASE2 = type2_RETBASE2[n];
      bi[] = to_vector(bMat1[type2_id[n], ]);// use bMat1 not b_mat
      gauss_total = 0;
      for (i in 1:15) {
        bspline_basis[] = to_vector(x_scaled_b_type2[,i,n]);
        bspline_drv_basis[] = to_vector(drv_x_scaled_b_type2[,i,n]);
        x = x_scaled_type2[i,n];
        gaussH1 = 0;
        gaussH2 = 0;
        wt = wt_scaled_type2[i,n];
        q12 = evaluate_h12(bspline_basis, bspline_drv_basis,
                          x, e12_aux, e12_beta,
                          HBAEL, GROUP2, RETBASE2,
                          a12_beta,
                          yGamma1[1], yGamma2[1], yGamma3[1],
                          yBeta1, yBeta2, yBeta3,
                          bi,
                          yXbar1, yXbar2, yXbar3
                          );
        for (j in 1:15) {
          bspline_basis2[] = to_vector(x2_scaled_b_type2[,j,i,n]);
          bspline_basis3[] = to_vector(x3_scaled_b_type2[,j,i,n]);
          bspline_drv_basis2[] = to_vector(drv_x2_scaled_b_type2[,j,i,n]);
          bspline_drv_basis3[] = to_vector(drv_x3_scaled_b_type2[,j,i,n]);
          x2 = x2_scaled_type2[j,i,n];
          x3 = x3_scaled_type2[j,i,n];
          wt2 = wt2_scaled_type2[j,i,n];
          wt3 = wt3_scaled_type2[j,i,n];
          H1 = wt2 * evaluate_h12(bspline_basis2, bspline_drv_basis2,
                                  x2, e12_aux, e12_beta,
                                  HBAEL, GROUP2, RETBASE2,
                                  a12_beta,
                                  yGamma1[1], yGamma2[1], yGamma3[1],
                                  yBeta1, yBeta2, yBeta3,
                                  bi,
                                  yXbar1, yXbar2, yXbar3
                                  );
          H2 = wt3 * evaluate_h23(bspline_basis3, 
                                  x3, e23_aux, e23_beta,
                                  HBAEL, GROUP2, RETBASE2,
                                  a23_beta,
                                  yGamma1[1], yGamma2[1], yGamma3[1],
                                  yBeta1, yBeta2, yBeta3,
                                  bi,
                                  yXbar1, yXbar2, yXbar3
                                  );
          gaussH1 = gaussH1 + H1;
          gaussH2 = gaussH2 + H2;
          //piece = exp(-gaussH1)*q12*exp(-gaussH2);
        }
        piece = exp(-gaussH1)*q12*exp(-gaussH2);
        weighted_piece = wt * piece;
        //gauss_total = gauss_total + wt * piece;
        gauss_total = gauss_total + weighted_piece;
      }
      p12 = gauss_total;
      target += log(p12);
    }
    
  }

  { // TYPE 4 TRANSITIONS, USING THE PRE-BUILT SPLINES MATRICES, AND REQUIRING LOOPS
    vector[6] bspline_basis;
    vector[6] bspline_basis2;
    vector[6] bspline_basis3;
    vector[6] bspline_drv_basis;
    vector[6] bspline_drv_basis2;
    vector[6] bspline_drv_basis3;
    real x;
    real x2;
    real x3;
    real wt;
    real wt2;
    real wt3;
    real HBAEL;
    real GROUP2;
    real RETBASE2;
    vector[15] bi;
    
    real gauss_total_p11;
    real H1_1;
    real p11;
    
    real gauss_total_p12;
    real q12;
    real gaussH1;
    real H1_2;
    real gaussH2;
    real H2_2;
    real piece;
    real weighted_piece;//NEW
    real p12;
    
    for (n in 1:type4_N) {
      HBAEL = type4_HBAEL[n];
      GROUP2 = type4_GROUP2[n];
      RETBASE2 = type4_RETBASE2[n];
      bi[] = to_vector(bMat1[type4_id[n], ]);// use bMat1 not b_mat
      
      gauss_total_p11 = 0;
      for (i in 1:15) {
        bspline_basis[] = to_vector(x_scaled_b_type4[,i,n]);
        bspline_drv_basis[] = to_vector(drv_x_scaled_b_type4[,i,n]);
        x = x_scaled_type4[i,n];
        wt = wt_scaled_type4[i,n];
        H1_1 = wt * evaluate_h12(bspline_basis, bspline_drv_basis,
                                x, e12_aux, e12_beta,
                                HBAEL, GROUP2, RETBASE2,
                                a12_beta,
                                yGamma1[1], yGamma2[1], yGamma3[1],
                                yBeta1, yBeta2, yBeta3,
                                bi,
                                yXbar1, yXbar2, yXbar3
                                );
        gauss_total_p11 = gauss_total_p11 + H1_1;
      }
      p11 = exp(-gauss_total_p11);
      
      gauss_total_p12 = 0;
      for (i in 1:15) {
        bspline_basis[] = to_vector(x_scaled_b_type4[,i,n]);
        bspline_drv_basis[] = to_vector(drv_x_scaled_b_type4[,i,n]);
        x = x_scaled_type4[i,n];
        gaussH1 = 0;
        gaussH2 = 0;
        wt = wt_scaled_type4[i,n];
        q12 = evaluate_h12(bspline_basis, bspline_drv_basis,
                          x, e12_aux, e12_beta,
                          HBAEL, GROUP2, RETBASE2,
                          a12_beta,
                          yGamma1[1], yGamma2[1], yGamma3[1],
                          yBeta1, yBeta2, yBeta3,
                          bi,
                          yXbar1, yXbar2, yXbar3
                          );
        for (j in 1:15) {
          bspline_basis2[] = to_vector(x2_scaled_b_type4[,j,i,n]);
          bspline_basis3[] = to_vector(x3_scaled_b_type4[,j,i,n]);
          bspline_drv_basis2[] = to_vector(drv_x2_scaled_b_type4[,j,i,n]);
          bspline_drv_basis3[] = to_vector(drv_x3_scaled_b_type4[,j,i,n]);
          x2 = x2_scaled_type4[j,i,n];
          x3 = x3_scaled_type4[j,i,n];
          wt2 = wt2_scaled_type4[j,i,n];
          wt3 = wt3_scaled_type4[j,i,n];
          H1_2 = wt2 * evaluate_h12(bspline_basis2, bspline_drv_basis2,
                                  x2, e12_aux, e12_beta,
                                  HBAEL, GROUP2, RETBASE2,
                                  a12_beta,
                                  yGamma1[1], yGamma2[1], yGamma3[1],
                                  yBeta1, yBeta2, yBeta3,
                                  bi,
                                  yXbar1, yXbar2, yXbar3
                                  );
          H2_2 = wt3 * evaluate_h23(bspline_basis3, 
                                  x3, e23_aux, e23_beta,
                                  HBAEL, GROUP2, RETBASE2,
                                  a23_beta,
                                  yGamma1[1], yGamma2[1], yGamma3[1],
                                  yBeta1, yBeta2, yBeta3,
                                  bi,
                                  yXbar1, yXbar2, yXbar3
                                  );
          gaussH1 = gaussH1 + H1_2;
          gaussH2 = gaussH2 + H2_2;
          //piece = exp(-gaussH1)*q12*exp(-gaussH2);
        }
        piece = exp(-gaussH1)*q12*exp(-gaussH2);
        weighted_piece = wt * piece;
        //gauss_total_p12 = gauss_total_p12 + wt * piece;
        gauss_total_p12 = gauss_total_p12 + weighted_piece;
      }
      p12 = gauss_total_p12;
      target += log(p11 + p12);
    }
    
  }

  { // TYPE 6 TRANSITIONS, USING THE PRE-BUILT SPLINES MATRICES, AND REQUIRING LOOPS
    vector[6] bspline_basis;
    vector[6] bspline_basis2;
    vector[6] bspline_basis3;
    vector[6] bspline_drv_basis;
    vector[6] bspline_drv_basis2;
    vector[6] bspline_drv_basis3;
    real x;
    real x2;
    real x3;
    real wt;
    real wt2;
    real wt3;
    real HBAEL;
    real GROUP2;
    real RETBASE2;
    vector[15] bi;
    
    real gauss_total;
    real q12;
    real gaussH1;
    real H1;
    real gaussH2;
    real H2;
    real piece;
    real weighted_piece;//NEW
    real p12;
    
    for (n in 1:type6_N) {
      HBAEL = type6_HBAEL[n];
      GROUP2 = type6_GROUP2[n];
      RETBASE2 = type6_RETBASE2[n];
      bi[] = to_vector(bMat1[type6_id[n], ]);// use bMat1 not b_mat
      gauss_total = 0;
      for (i in 1:15) {
        bspline_basis[] = to_vector(x_scaled_b_type6[,i,n]);
        bspline_drv_basis[] = to_vector(drv_x_scaled_b_type6[,i,n]);
        x = x_scaled_type6[i,n];
        gaussH1 = 0;
        gaussH2 = 0;
        wt = wt_scaled_type6[i,n];
        q12 = evaluate_h12(bspline_basis, bspline_drv_basis,
                          x, e12_aux, e12_beta,
                          HBAEL, GROUP2, RETBASE2,
                          a12_beta,
                          yGamma1[1], yGamma2[1], yGamma3[1],
                          yBeta1, yBeta2, yBeta3,
                          bi,
                          yXbar1, yXbar2, yXbar3
                          );
        for (j in 1:15) {
          bspline_basis2[] = to_vector(x2_scaled_b_type6[,j,i,n]);
          bspline_basis3[] = to_vector(x3_scaled_b_type6[,j,i,n]);
          bspline_drv_basis2[] = to_vector(drv_x2_scaled_b_type6[,j,i,n]);
          bspline_drv_basis3[] = to_vector(drv_x3_scaled_b_type6[,j,i,n]);
          x2 = x2_scaled_type6[j,i,n];
          x3 = x3_scaled_type6[j,i,n];
          wt2 = wt2_scaled_type6[j,i,n];
          wt3 = wt3_scaled_type6[j,i,n];
          H1 = wt2 * evaluate_h12(bspline_basis2, bspline_drv_basis2,
                                  x2, e12_aux, e12_beta,
                                  HBAEL, GROUP2, RETBASE2,
                                  a12_beta,
                                  yGamma1[1], yGamma2[1], yGamma3[1],
                                  yBeta1, yBeta2, yBeta3,
                                  bi,
                                  yXbar1, yXbar2, yXbar3
                                  );
          H2 = wt3 * evaluate_h23(bspline_basis3, 
                                  x3, e23_aux, e23_beta,
                                  HBAEL, GROUP2, RETBASE2,
                                  a23_beta,
                                  yGamma1[1], yGamma2[1], yGamma3[1],
                                  yBeta1, yBeta2, yBeta3,
                                  bi,
                                  yXbar1, yXbar2, yXbar3
                                  );
          gaussH1 = gaussH1 + H1;
          gaussH2 = gaussH2 + H2;
          //piece = exp(-gaussH1)*q12*exp(-gaussH2);
        }
        piece = exp(-gaussH1)*q12*exp(-gaussH2);
        weighted_piece = wt * piece;
        //gauss_total = gauss_total + wt * piece;
        gauss_total = gauss_total + weighted_piece;
      }
      p12 = gauss_total;
      target += log(p12);
    }
    
  }

  // Log priors
  // increments target with mvmer priors
#include RSTANARM_STAN_FILES/model/priors_mvmer.stan

  // coefficients for event submodel
  target += normal_lpdf(e12_z_beta | 0, 1);
  target += normal_lpdf(e23_z_beta | 0, 1);
  target += normal_lpdf(a12_z_beta | 0, 1);
  target += normal_lpdf(a23_z_beta | 0, 1);
  
  // b-spline coefs for baseline hazard
  target += normal_lpdf(e12_aux_unscaled | 0, 1);
  target += normal_lpdf(e23_aux_unscaled | 0, 1);
}
generated quantities {
#include RSTANARM_STAN_FILES/gqs/gen_quantities_mvmer.stan
}
