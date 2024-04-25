//-----------------------------
// Fit independent factor analysis model (LVM)
// Species-specific coefficients are treated as random slopes and drawn from a common (multivariate normal) distribution.
//-----------------------------
#include <math.h>
#include <TMB.hpp> 

template<class Type>


Type objective_function<Type>::operator() () { 
    using namespace density;
    using namespace Eigen;


    // Declare data
    DATA_MATRIX(y); 
    DATA_MATRIX(X); 
    DATA_INTEGER(num_lv);
    int num_units = y.rows(); 
    int num_spp = y.cols();
    int num_X = X.cols();

    // Declare parameters
    PARAMETER_MATRIX(betas);
    PARAMETER_VECTOR(alpha);
    PARAMETER_VECTOR(log_sds);
    PARAMETER_VECTOR(loadings);
    PARAMETER_MATRIX(lvs);

    // Negative log-likelihood and other parameter set up
    Type nll = Type(0.0);
    vector<Type> sds = exp(log_sds);

    // Set up latent variables
    for(int i=0; i<num_units; i++) { for(int k=0; k<num_lv; k++) {
        nll -= dnorm(lvs(i,k), Type(0.0), Type(1.0), true);
        } }

    // Set up shared random effects distribution for the betas 
     for(int i=0; i<num_spp; i++) { for(int k=0; k<num_X; k++) {
       nll -= dnorm(betas(i,k), alpha(k), sds(k), true);
       } }

    //To create loadings as upper triangular matrix 
    matrix<Type> loadings_mat(num_spp, num_lv); 
    int loadings_ind=0;
    for (int j=0; j<num_lv; j++) { for (int i=0; i<num_spp; i++) {
        if (j > i) { loadings_mat(i,j) = 0; } 
        if (j <= i) { loadings_mat(i,j) = loadings(loadings_ind); loadings_ind += 1; } 
        } }

    // Data likelihood
    matrix<Type> eta_X = X*betas.transpose();
    matrix<Type> eta_LV = lvs*loadings_mat.transpose();
    for (int i=0; i<num_units; i++) { for (int j=0; j<num_spp; j++) {
        nll -= dbinom(y(i,j), Type(1.0), invlogit(eta_X(i,j)+eta_LV(i,j)),true);
        } }

    ADREPORT(betas); 
    ADREPORT(alpha); 
    ADREPORT(sds); 
    ADREPORT(loadings_mat);
    ADREPORT(lvs);

    return nll; 
    }
