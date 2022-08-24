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
  // THESE NEED TO INCORPORATE THE LONGITUDINAL OUTCOME
  vector evaluate_h12(vector time,                             // THIS ASSUMES TIME IS ALWAYS A VECTOR OF 15
                      vector beta1,                            // regression coefficient for the baseline hazards
                      vector X,                                // fixed covariates 
                      vector beta2,                            // regression coefficients for the fixed covariates
                      vector beta3,                            // regression coefficient(s) for the longitudinal outcome - meant for "alpha12"
                      
                      vector betay1,
                      vector betay2,                            // longitudinal submodel parameters
                      vector betay3,
                      real betay1_int, real betay2_int, real betay3_int,
                      vector bi,                               // random effects values
                      real pwc_knot1, real pwc_knot2,
                      vector y_xbar){ 
    vector[15] h12;
    
    for (i in 1:15) {
      h12[i] = exp( beta1[1]*(time[i]<=pwc_knot1) + beta1[2]*(pwc_knot1<time[i] && time[i]<=pwc_knot2) + beta1[3]*(pwc_knot2<time[i]) +
                    beta2[1]*X[1] + beta2[2]*X[2] +
                    beta3[1]*( betay1_int - dot_product(y_xbar, betay1) + betay1[1]*time[i] + betay1[2]*X[1] + betay1[3]*X[3] + 
                              bi[1] + bi[2]*time[i] ) +
                    beta3[2]*( betay2_int - dot_product(y_xbar, betay2) + betay2[1]*time[i] + betay2[2]*X[1] + betay2[3]*X[3] + 
                              bi[3] + bi[4]*time[i] ) +
                    beta3[3]*( betay3_int - dot_product(y_xbar, betay3) + betay3[1]*time[i] + betay3[2]*X[1] + betay3[3]*X[3] + 
                              bi[5] + bi[6]*time[i] )
                    );
    }
    return h12;
  }
  vector evaluate_h23(vector time, 
                      vector beta1, 
                      vector X,                                // fixed covariates 
                      vector beta2,                            // regression coefficients for the fixed covariates
                      vector beta3,                            // regression coefficient(s) for the longitudinal outcome
                      
                      vector betay1,                          // longitudinal submodel parameters
                      vector betay2,                            
                      vector betay3,
                      real betay1_int, real betay2_int, real betay3_int,
                      vector bi,                               // random effects values   
                      real pwc_knot1, real pwc_knot2,
                      vector y_xbar){
    vector[15] h23;
    
    for (i in 1:15) {
      h23[i] = exp( beta1[1]*(time[i]<=pwc_knot1) + beta1[2]*(pwc_knot1<time[i] && time[i]<=pwc_knot2) + beta1[3]*(pwc_knot2<time[i]) +
                    beta2[1]*X[1] +
                    beta3[1]*( betay1_int - dot_product(y_xbar, betay1) + betay1[1]*time[i] + betay1[2]*X[1] + betay1[3]*X[2] + // NOTE THE ORDER OF X - ONLY 2 FOR 2->3 TRANSITION
                              bi[1] + bi[2]*time[i] ) +
                    beta3[2]*( betay2_int - dot_product(y_xbar, betay2) + betay2[1]*time[i] + betay2[2]*X[1] + betay2[3]*X[2] + 
                              bi[3] + bi[4]*time[i] ) +
                    beta3[3]*( betay3_int - dot_product(y_xbar, betay3) + betay3[1]*time[i] + betay3[2]*X[1] + betay3[3]*X[2] + 
                              bi[5] + bi[6]*time[i] )
                    );
    }
    return h23;
  }
  
}
data {
#include RSTANARM_STAN_FILES/data/dimensions_mvmer.stan
#include RSTANARM_STAN_FILES/data/data_mvmer.stan

  // EVENT DATA
  // start with parameter info
  int<lower=0> e_K[2];           // num. of predictors in event submodel, [1] is 1->2 transition, [2] is 2->3 transition
  int<lower=0> basehaz_df;    // "df" for baseline hazard, so it's 3 since there are 2 knots in each 1->2 and 2->3
  int<lower=0> a12_K;           // num. of association parameters
  int<lower=0> a23_K;
  
  //type 1
  int<lower=0> type1_N_GK;
  int<lower=1, upper=bN1> type1_id[type1_N_GK]; 
  vector[type1_N_GK] type1_qpoints;// PRE-GENERATED QUADRATURE POINTS AND WEIGHTS
  vector[type1_N_GK] type1_wpoints;
  int<lower=0> nrow_e_Xq_type1;     // num. rows in event submodel predictor matrix
  matrix[nrow_e_Xq_type1,e_K[1]] e_Xq_type1; // predictor matrix (event submodel)
  matrix[nrow_e_Xq_type1,basehaz_df] basehaz_X_type1; // design matrix for baseline hazard

  int<lower=0> nrow_y_Xq_type1; // num. rows in long. predictor matrix at quadpoints

  matrix[nrow_y_Xq_type1,yK[1]] y_Xq_type1; // fe design matrix at quadpoints      // NOTE THE MAJOR SIMPLIFICATION HERE: y_K[1]=y_K[2]=y_K[3] 
  vector[nrow_y_Xq_type1] y_Zq_type1[bK1_len[1]]; // re design matrix at quadpoints
  int<lower=0> y_Zq_id_type1[nrow_y_Xq_type1]; // group indexing for re design matrix
  
  //type 2
  int<lower=0> type2_N;
  vector[type2_N] type2_year;
  vector[type2_N] type2_lag_year;
  int<lower=1, upper=bN1> type2_id[type2_N];
  vector[type2_N] type2_X1;
  vector[type2_N] type2_X3;
  vector[type2_N] type2_X2;
  
  //type 3
  int<lower=0> type3_N_GK;
  int<lower=1, upper=bN1> type3_id[type3_N_GK];
  vector[type3_N_GK] type3_qpoints;// PRE-GENERATED QUADRATURE POINTS AND WEIGHTS
  vector[type3_N_GK] type3_wpoints;
  int<lower=0> nrow_e_Xq_type3;     // num. rows in event submodel predictor matrix
  matrix[nrow_e_Xq_type3,e_K[2]] e_Xq_type3; // predictor matrix (event submodel)  // CAREFUL, e_K[2] NOT e_K[1]
  matrix[nrow_e_Xq_type3,basehaz_df] basehaz_X_type3; // design matrix for baseline hazard
  
  int<lower=0> nrow_y_Xq_type3; // num. rows in long. predictor matrix at quadpoints
  
  matrix[nrow_y_Xq_type3,yK[1]] y_Xq_type3; // fe design matrix at quadpoints
  vector[nrow_y_Xq_type3] y_Zq_type3[bK1_len[1]]; // re design matrix at quadpoints
  int<lower=0> y_Zq_id_type3[nrow_y_Xq_type3]; // group indexing for re design matrix
  
  //type 4
  int<lower=0> type4_N;
  vector[type4_N] type4_year;
  vector[type4_N] type4_lag_year;
  int<lower=1, upper=bN1> type4_id[type4_N];
  vector[type4_N] type4_X1;
  vector[type4_N] type4_X3;
  vector[type4_N] type4_X2;
  
  //type 5
  int<lower=0> type5_N_GK;
  int<lower=1, upper=bN1> type5_id[type5_N_GK];
  vector[type5_N_GK] type5_qpoints;// PRE-GENERATED QUADRATURE POINTS AND WEIGHTS
  vector[type5_N_GK] type5_wpoints;
  int<lower=0> nrow_e_Xq_type5;     // num. rows in event submodel predictor matrix
  matrix[nrow_e_Xq_type5,e_K[2]] e_Xq_type5; // predictor matrix (event submodel)
  matrix[nrow_e_Xq_type5,basehaz_df] basehaz_X_type5; // design matrix for baseline hazard
  
  int<lower=0> nrow_y_Xq_type5; // num. rows in long. predictor matrix at quadpoints
  
  matrix[nrow_y_Xq_type5,yK[1]] y_Xq_type5; // fe design matrix at quadpoints
  vector[nrow_y_Xq_type5] y_Zq_type5[bK1_len[1]]; // re design matrix at quadpoints
  int<lower=0> y_Zq_id_type5[nrow_y_Xq_type5]; // group indexing for re design matrix
  
  //type 6
  int<lower=0> type6_N;
  vector[type6_N] type6_year;
  vector[type6_N] type6_lag_year;
  int<lower=1, upper=bN1> type6_id[type6_N];
  vector[type6_N] type6_X1;
  vector[type6_N] type6_X3;
  vector[type6_N] type6_X2;
  
  //type 7
  int<lower=0> type7_N_GK;
  int<lower=1, upper=bN1> type7_id[type7_N_GK];
  vector[type7_N_GK] type7_qpoints;// PRE-GENERATED QUADRATURE POINTS AND WEIGHTS
  vector[type7_N_GK] type7_wpoints;
  int<lower=0> nrow_e_Xq_type7;     // num. rows in event submodel predictor matrix
  matrix[nrow_e_Xq_type7,e_K[2]] e_Xq_type7; // predictor matrix (event submodel)
  matrix[nrow_e_Xq_type7,basehaz_df] basehaz_X_type7; // design matrix for baseline hazard
  
  int<lower=0> nrow_y_Xq_type7; // num. rows in long. predictor matrix at quadpoints
  
  matrix[nrow_y_Xq_type7,yK[1]] y_Xq_type7; // fe design matrix at quadpoints
  vector[nrow_y_Xq_type7] y_Zq_type7[bK1_len[1]]; // re design matrix at quadpoints
  int<lower=0> y_Zq_id_type7[nrow_y_Xq_type7]; // group indexing for re design matrix

  //EVENTS FOR STATE 3
  int<lower=0> events_N;
  vector[events_N] events_year;
  int<lower=1, upper=bN1> events_id[events_N];
  int<lower=0> nrow_e_Xq_events;     // num. rows in event submodel predictor matrix
  matrix[nrow_e_Xq_events,e_K[2]] e_Xq_events; // predictor matrix (event submodel)
  matrix[nrow_e_Xq_events,basehaz_df] basehaz_X_events; // design matrix for baseline hazard
  
  int<lower=0> nrow_y_Xq_events; // num. rows in long. predictor matrix at quadpoints
  
  matrix[nrow_y_Xq_events,yK[1]] y_Xq_events; // fe design matrix at quadpoints
  vector[nrow_y_Xq_events] y_Zq_events[bK1_len[1]]; // re design matrix at quadpoints
  int<lower=0> y_Zq_id_events[nrow_y_Xq_events]; // group indexing for re design matrix

  // EXTRA DATA FOR MY SPECIFIC APPLICATION (THE GK NODES AND WEIGHTS AND THE PWC KNOTS)
  vector[15] xk;
  vector[15] wt;
  
  real pwc_knot1;
  real pwc_knot2;
  
  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // NEW DATA TYPES FOR THE "EVALUATE_ETA" FUNCTION FROM RSTANARM
  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  vector[0] y1_offset_eta;
  vector[0] y2_offset_eta;
  vector[0] y3_offset_eta;
  //vector[has_offset[2] && assoc_uses[1,2] == 1 ? nrow_y_Xq[2] : 0] y2_offset_eta;
  //vector[has_offset[3] && assoc_uses[1,3] == 1 ? nrow_y_Xq[3] : 0] y3_offset_eta;
      // re design matrix at quadpoints, group factor 2
  vector[0] y1_z2q_eta[0];
  vector[0] y2_z2q_eta[0];
  vector[0] y3_z2q_eta[0];
  //vector[assoc_uses[1,1] == 1 && bK2_len[1] > 0 ? nrow_y_Xq[1] : 0] y1_z2q_eta[bK2_len[1]];//ETC. AS ABOVE
  //vector[assoc_uses[1,2] == 1 && bK2_len[2] > 0 ? nrow_y_Xq[2] : 0] y2_z2q_eta[bK2_len[2]];
  //vector[assoc_uses[1,3] == 1 && bK2_len[3] > 0 ? nrow_y_Xq[3] : 0] y3_z2q_eta[bK2_len[3]];
  int<lower=0> y1_z2q_id_eta[0];
  int<lower=0> y2_z2q_id_eta[0];
  int<lower=0> y3_z2q_id_eta[0];
  //int<lower=0> y1_z2q_id_eta[assoc_uses[1,1] == 1 && bK2_len[1] > 0 ? nrow_y_Xq[1] : 0];//ETC> AS ABOVE
  //int<lower=0> y2_z2q_id_eta[assoc_uses[1,2] == 1 && bK2_len[2] > 0 ? nrow_y_Xq[2] : 0];
  //int<lower=0> y3_z2q_id_eta[assoc_uses[1,3] == 1 && bK2_len[3] > 0 ? nrow_y_Xq[3] : 0];
  
  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // THE HYPERPARAMETERS
  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#include RSTANARM_STAN_FILES/data/hyperparameters_mvmer.stan

// NEXT THE HYPERPARAMETERS FOR THE EVENTS:
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
    vector[nrow_y_Xq_type1] y3_eta_q;
    
    y1_eta_q = evaluate_eta(y_Xq_type1, y_Zq_type1, y1_z2q_eta, //y1_xq_type1, y1_z1q_type1, 
                                    y_Zq_id_type1, y1_z2q_id_eta, //y1_z1q_id_type1 these are all the same for this application (across the 3 outcomes)
                                    yGamma1, yBeta1, bMat1, bMat2,
                                    0, 0, 1, //bMat1_colshift, bMat2_colshift, 0,
                                              y1_offset_eta);
    y2_eta_q = evaluate_eta(y_Xq_type1, y_Zq_type1, y2_z2q_eta,
                                    y_Zq_id_type1, y2_z2q_id_eta,
                                    yGamma2, yBeta2, bMat1, bMat2,
                                    2, 0, 1, //bMat1_colshift, bMat2_colshift, 0,
                                              y2_offset_eta);
    y3_eta_q = evaluate_eta(y_Xq_type1, y_Zq_type1, y3_z2q_eta,
                                    y_Zq_id_type1, y3_z2q_id_eta,
                                    yGamma3, yBeta3, bMat1, bMat2,
                                    4, 0, 1, //bMat1_colshift, bMat2_colshift, 0,
                                              y3_offset_eta);
    
    eta = basehaz_X_type1 * e12_aux + e_Xq_type1 * e12_beta + a12_beta[1] * y1_eta_q + a12_beta[2] * y2_eta_q + a12_beta[3] * y3_eta_q
          - a12_beta[1] * dot_product(yXbar1, yBeta1) - a12_beta[2] * dot_product(yXbar2, yBeta2) - a12_beta[3] * dot_product(yXbar3, yBeta3);
    target += - dot_product(type1_wpoints,exp(eta));
  }
  
  { // Type 3 transitions - PRE-DEFINED QUADRATURE POINTS
    vector[nrow_e_Xq_type3] eta;
    vector[nrow_y_Xq_type3] y1_eta_q;
    vector[nrow_y_Xq_type3] y2_eta_q;
    vector[nrow_y_Xq_type3] y3_eta_q;
    
    y1_eta_q = evaluate_eta(y_Xq_type3, y_Zq_type3, y1_z2q_eta, //y1_xq_type3, y1_z1q_type3, 
                                    y_Zq_id_type3, y1_z2q_id_eta, //y1_z1q_id_type3 these are all the same for this application (across the 3 outcomes)
                                    yGamma1, yBeta1, bMat1, bMat2,
                                    0, 0, 1, //bMat1_colshift, bMat2_colshift, 0, //intercept_type should be 1!
                                              y1_offset_eta);
    y2_eta_q = evaluate_eta(y_Xq_type3, y_Zq_type3, y2_z2q_eta,
                                    y_Zq_id_type3, y2_z2q_id_eta,
                                    yGamma2, yBeta2, bMat1, bMat2,
                                    2, 0, 1, //bMat1_colshift, bMat2_colshift, 0,
                                              y2_offset_eta);
    y3_eta_q = evaluate_eta(y_Xq_type3, y_Zq_type3, y3_z2q_eta,
                                    y_Zq_id_type3, y3_z2q_id_eta,
                                    yGamma3, yBeta3, bMat1, bMat2,
                                    4, 0, 1, //bMat1_colshift, bMat2_colshift, 0,
                                              y3_offset_eta);
    
    eta = basehaz_X_type3 * e23_aux + e_Xq_type3 * e23_beta + a23_beta[1] * y1_eta_q + a23_beta[2] * y2_eta_q + a23_beta[3] * y3_eta_q
          - a23_beta[1] * dot_product(yXbar1, yBeta1) - a23_beta[2] * dot_product(yXbar2, yBeta2) - a23_beta[3] * dot_product(yXbar3, yBeta3);
    target += - dot_product(type3_wpoints,exp(eta));
  }
  
  { // Type 5 transitions - PRE-DEFINED QUADRATURE POINTS
    vector[nrow_e_Xq_type5] eta;
    vector[nrow_y_Xq_type5] y1_eta_q;
    vector[nrow_y_Xq_type5] y2_eta_q;
    vector[nrow_y_Xq_type5] y3_eta_q;
    
    y1_eta_q = evaluate_eta(y_Xq_type5, y_Zq_type5, y1_z2q_eta, //y1_xq_type5, y1_z1q_type5, 
                                    y_Zq_id_type5, y1_z2q_id_eta, //y1_z1q_id_type5 these are all the same for this application (across the 3 outcomes)
                                    yGamma1, yBeta1, bMat1, bMat2,
                                    0, 0, 1, //bMat1_colshift, bMat2_colshift, 0,
                                              y1_offset_eta);
    y2_eta_q = evaluate_eta(y_Xq_type5, y_Zq_type5, y2_z2q_eta,
                                    y_Zq_id_type5, y2_z2q_id_eta,
                                    yGamma2, yBeta2, bMat1, bMat2,
                                    2, 0, 1, //bMat1_colshift, bMat2_colshift, 0,
                                              y2_offset_eta);
    y3_eta_q = evaluate_eta(y_Xq_type5, y_Zq_type5, y3_z2q_eta,
                                    y_Zq_id_type5, y3_z2q_id_eta,
                                    yGamma3, yBeta3, bMat1, bMat2,
                                    4, 0, 1, //bMat1_colshift, bMat2_colshift, 0,
                                              y3_offset_eta);
    
    eta = basehaz_X_type5 * e23_aux + e_Xq_type5 * e23_beta + a23_beta[1] * y1_eta_q + a23_beta[2] * y2_eta_q + a23_beta[3] * y3_eta_q
          - a23_beta[1] * dot_product(yXbar1, yBeta1) - a23_beta[2] * dot_product(yXbar2, yBeta2) - a23_beta[3] * dot_product(yXbar3, yBeta3);
    target += - dot_product(type5_wpoints,exp(eta));
  }

  { // Type 7 transitions - PRE-DEFINED QUADRATURE POINTS, AND ONLY THE TRANSITION PROBABILITY NOT THE EVENT TIME
    vector[nrow_e_Xq_type7] eta;
    vector[nrow_y_Xq_type7] y1_eta_q;
    vector[nrow_y_Xq_type7] y2_eta_q;
    vector[nrow_y_Xq_type7] y3_eta_q;
    
    y1_eta_q = evaluate_eta(y_Xq_type7, y_Zq_type7, y1_z2q_eta, //y1_xq_type7, y1_z1q_type7, 
                                    y_Zq_id_type7, y1_z2q_id_eta, //y1_z1q_id_type7 these are all the same for my application (across the 3 outcomes)
                                    yGamma1, yBeta1, bMat1, bMat2,
                                    0, 0, 1, //bMat1_colshift, bMat2_colshift, 0,
                                              y1_offset_eta);
    y2_eta_q = evaluate_eta(y_Xq_type7, y_Zq_type7, y2_z2q_eta,
                                    y_Zq_id_type7, y2_z2q_id_eta,
                                    yGamma2, yBeta2, bMat1, bMat2,
                                    2, 0, 1, //bMat1_colshift, bMat2_colshift, 0,
                                              y2_offset_eta);
    y3_eta_q = evaluate_eta(y_Xq_type7, y_Zq_type7, y3_z2q_eta,
                                    y_Zq_id_type7, y3_z2q_id_eta,
                                    yGamma3, yBeta3, bMat1, bMat2,
                                    4, 0, 1, //bMat1_colshift, bMat2_colshift, 0,
                                              y3_offset_eta);
    
    eta = basehaz_X_type7 * e23_aux + e_Xq_type7 * e23_beta + a23_beta[1] * y1_eta_q + a23_beta[2] * y2_eta_q + a23_beta[3] * y3_eta_q
          - a23_beta[1] * dot_product(yXbar1, yBeta1) - a23_beta[2] * dot_product(yXbar2, yBeta2) - a23_beta[3] * dot_product(yXbar3, yBeta3);
    target += - dot_product(type7_wpoints,exp(eta));
  }
  
  { // OBSERVED STATE 3 EVENTS
    vector[nrow_e_Xq_events] log_hazard;
    vector[nrow_y_Xq_events] y1_eta_q;
    vector[nrow_y_Xq_events] y2_eta_q;
    vector[nrow_y_Xq_events] y3_eta_q;
    
    y1_eta_q = evaluate_eta(y_Xq_events, y_Zq_events, y1_z2q_eta, //y1_xq_events, y1_z1q_events, 
                                    y_Zq_id_events, y1_z2q_id_eta, //y1_z1q_id_events these are all the same for my application (across the 3 outcomes)
                                    yGamma1, yBeta1, bMat1, bMat2,
                                    0, 0, 1, //bMat1_colshift, bMat2_colshift, 0,
                                              y1_offset_eta);
    y2_eta_q = evaluate_eta(y_Xq_events, y_Zq_events, y2_z2q_eta,
                                    y_Zq_id_events, y2_z2q_id_eta,
                                    yGamma2, yBeta2, bMat1, bMat2,
                                    2, 0, 1, //bMat1_colshift, bMat2_colshift, 0,
                                              y2_offset_eta);
    y3_eta_q = evaluate_eta(y_Xq_events, y_Zq_events, y3_z2q_eta,
                                    y_Zq_id_events, y3_z2q_id_eta,
                                    yGamma3, yBeta3, bMat1, bMat2,
                                    4, 0, 1, //bMat1_colshift, bMat2_colshift, 0,
                                              y3_offset_eta);
    
    log_hazard = basehaz_X_events * e23_aux + e_Xq_events * e23_beta + a23_beta[1] * y1_eta_q + a23_beta[2] * y2_eta_q + a23_beta[3] * y3_eta_q
                  - a23_beta[1] * dot_product(yXbar1, yBeta1) - a23_beta[2] * dot_product(yXbar2, yBeta2) - a23_beta[3] * dot_product(yXbar3, yBeta3);
    target += sum(log_hazard);
  }
  
  { // Type 2 transitions - vectorized version, but still 1 loop required 
    vector[3] X12;
    vector[3] X23;
    vector[6] bi;
    
    vector[15] x_scaled;
    vector[15] wt_scaled;
    vector[15] x2_scaled;
    vector[15] wt2_scaled;
    vector[15] x3_scaled;
    vector[15] wt3_scaled;
    vector[15] H1_iterated;
    vector[15] H2_iterated;
    
    real p12;
    
    for (n in 1:type2_N) {
      X12[1]= type2_X1[n];
      X12[2] = type2_X3[n]; 
      X12[3] = type2_X2[n]; //  NOTE THE ORDER OF COVARIATES
      X23[1]= type2_X1[n];
      X23[2]= type2_X2[n];
      bi[] = to_vector(bMat1[type2_id[n], ]);
      
      x_scaled = ((xk+1)/2)*(type2_year[n]-type2_lag_year[n])+type2_lag_year[n];
      wt_scaled = ((type2_year[n]-type2_lag_year[n])/2)*wt;
      
      for (k in 1:15){
        x2_scaled = ((xk+1)/2)*(x_scaled[k]-type2_lag_year[n])+type2_lag_year[n];
        wt2_scaled = ((x_scaled[k]-type2_lag_year[n])/2)*wt;
        
        x3_scaled = ((xk+1)/2)*(type2_year[n]-x_scaled[k])+x_scaled[k];
        wt3_scaled = ((type2_year[n]-x_scaled[k])/2)*wt;
        
        H1_iterated[k] = sum(wt2_scaled .* evaluate_h12(x2_scaled, e12_aux, X12, e12_beta,
                                                        a12_beta, 
                                                        yBeta1, yBeta2, yBeta3, 
                                                        yGamma1[1], yGamma2[1], yGamma3[1],
                                                        bi,
                                                        pwc_knot1, pwc_knot2,
                                                        yXbar1)); // note the vector multiplication symbol .*
        H2_iterated[k] = sum(wt3_scaled .* evaluate_h23(x3_scaled, e23_aux, X23, e23_beta,
                                                        a23_beta, 
                                                        yBeta1, yBeta2, yBeta3, 
                                                        yGamma1[1], yGamma2[1], yGamma3[1],
                                                        bi,
                                                        pwc_knot1, pwc_knot2,
                                                        yXbar1)); // note the vector multiplication symbol .*
      }
      p12 = sum(wt_scaled .* ( exp(-H1_iterated) .* evaluate_h12(x_scaled, e12_aux, X12, e12_beta,
                                                                a12_beta, 
                                                                yBeta1, yBeta2, yBeta3, 
                                                                yGamma1[1], yGamma2[1], yGamma3[1],
                                                                bi,
                                                                pwc_knot1, pwc_knot2,
                                                                yXbar1) .* 
                exp(-H2_iterated) ) );
      target += log(p12);
    }
  }
  
  { // Type 4 transitions (i.e. censoring from state 1) - vectorized version
    vector[3] X12;
    vector[2] X23;
    vector[6] bi;
    
    vector[15] x_scaled;
    vector[15] wt_scaled;
    vector[15] x2_scaled;
    vector[15] wt2_scaled;
    vector[15] x3_scaled;
    vector[15] wt3_scaled;
    vector[15] H1_iterated;
    vector[15] H2_iterated;
    
    real H1;
    real p11;
    real p12;
    
    for (n in 1:type4_N) {
      X12[1]= type4_X1[n];
      X12[2] = type4_X3[n];
      X12[3] = type4_X2[n]; //
      X23[1]= type4_X1[n];
      X23[2]= type4_X2[n];
      bi[] = to_vector(bMat1[type4_id[n], ]);
      
      x_scaled = ((xk+1)/2)*(type4_year[n]-type4_lag_year[n])+type4_lag_year[n];
      wt_scaled = ((type4_year[n]-type4_lag_year[n])/2)*wt;
      H1 = sum(wt_scaled .* evaluate_h12(x_scaled, e12_aux, X12, e12_beta,
                                                        a12_beta, 
                                                        yBeta1, yBeta2, yBeta3, 
                                                        yGamma1[1], yGamma2[1], yGamma3[1],
                                                        bi,
                                                        pwc_knot1, pwc_knot2,
                                                        yXbar1));
      p11 = exp(-H1);
      for (k in 1:15){
        x2_scaled = ((xk+1)/2)*(x_scaled[k]-type4_lag_year[n])+type4_lag_year[n];
        wt2_scaled = ((x_scaled[k]-type4_lag_year[n])/2)*wt;
        
        x3_scaled = ((xk+1)/2)*(type4_year[n]-x_scaled[k])+x_scaled[k];
        wt3_scaled = ((type4_year[n]-x_scaled[k])/2)*wt;
        
        H1_iterated[k] = sum(wt2_scaled .* evaluate_h12(x2_scaled, e12_aux, X12, e12_beta,
                                                        a12_beta, 
                                                        yBeta1, yBeta2, yBeta3, 
                                                        yGamma1[1], yGamma2[1], yGamma3[1],
                                                        bi,
                                                        pwc_knot1, pwc_knot2,
                                                        yXbar1)); // note the vector multiplication symbol .*
        H2_iterated[k] = sum(wt3_scaled .* evaluate_h23(x3_scaled, e23_aux, X23, e23_beta,
                                                        a23_beta, 
                                                        yBeta1, yBeta2, yBeta3, 
                                                        yGamma1[1], yGamma2[1], yGamma3[1],
                                                        bi,
                                                        pwc_knot1, pwc_knot2,
                                                        yXbar1)); // note the vector multiplication symbol .*
      }
      p12 = sum(wt_scaled .* (exp(-H1_iterated) .* evaluate_h12(x_scaled, e12_aux, X12, e12_beta,
                                                        a12_beta, 
                                                        yBeta1, yBeta2, yBeta3, 
                                                        yGamma1[1], yGamma2[1], yGamma3[1],
                                                        bi,
                                                        pwc_knot1, pwc_knot2,
                                                        yXbar1) .* exp(-H2_iterated)) );
      target += log(p11+p12);
    }
  }
  
  { // type 6 transitions (i.e. entry into state 3, from state 1) - WE OMIT THE EVENTS - vectorized version
    vector[3] X12;
    vector[2] X23;
    vector[6] bi;
    
    vector[15] x_scaled;
    vector[15] wt_scaled;
    vector[15] x2_scaled;
    vector[15] wt2_scaled;
    vector[15] x3_scaled;
    vector[15] wt3_scaled;
    vector[15] H1_iterated;
    vector[15] H2_iterated;
    
    real p12;
    
    for (n in 1:type6_N) {
      X12[1]= type6_X1[n];
      X12[2] = type6_X3[n];
      X12[3] = type6_X2[n];
      X23[1]= type6_X1[n];
      X23[2]= type6_X2[n];
      bi[] = to_vector(bMat1[type6_id[n], ]);
      
      x_scaled = ((xk+1)/2)*(type6_year[n]-type6_lag_year[n])+type6_lag_year[n];
      wt_scaled = ((type6_year[n]-type6_lag_year[n])/2)*wt;
      
      for (k in 1:15){
        x2_scaled = ((xk+1)/2)*(x_scaled[k]-type6_lag_year[n])+type6_lag_year[n];
        wt2_scaled = ((x_scaled[k]-type6_lag_year[n])/2)*wt;
        
        x3_scaled = ((xk+1)/2)*(type6_year[n]-x_scaled[k])+x_scaled[k];
        wt3_scaled = ((type6_year[n]-x_scaled[k])/2)*wt;
        
        H1_iterated[k] = sum(wt2_scaled .* evaluate_h12(x2_scaled, e12_aux, X12, e12_beta,
                                                        a12_beta, 
                                                        yBeta1, yBeta2, yBeta3, 
                                                        yGamma1[1], yGamma2[1], yGamma3[1],
                                                        bi,
                                                        pwc_knot1, pwc_knot2,
                                                        yXbar1)); // note the vector multiplication symbol .*
        H2_iterated[k] = sum(wt3_scaled .* evaluate_h23(x3_scaled, e23_aux, X23, e23_beta,
                                                        a23_beta, 
                                                        yBeta1, yBeta2, yBeta3, 
                                                        yGamma1[1], yGamma2[1], yGamma3[1],
                                                        bi,
                                                        pwc_knot1, pwc_knot2,
                                                        yXbar1)); // note the vector multiplication symbol .*
      }
      p12 = sum(wt_scaled .* (exp(-H1_iterated) .* evaluate_h12(x_scaled, e12_aux, X12, e12_beta,
                                                                a12_beta, 
                                                                yBeta1, yBeta2, yBeta3, 
                                                                yGamma1[1], yGamma2[1], yGamma3[1],
                                                                bi,
                                                                pwc_knot1, pwc_knot2,
                                                                yXbar1) .* exp(-H2_iterated)) );
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
