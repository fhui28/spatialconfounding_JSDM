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
    PARAMETER_VECTOR(log_dispparam);
    PARAMETER_VECTOR(minusonelogit_powerparam);
    PARAMETER_VECTOR(loadings);
    PARAMETER_MATRIX(lvs);

    // Negative log-likelihood and other parameter set up
    Type nll = Type(0.0);
    vector<Type> dispparam = exp(log_dispparam);
    vector<Type> powerparam(num_spp);
    for(int i=0; i<num_spp; i++) {
        powerparam(i) = 1 + exp(minusonelogit_powerparam(i))/(1+exp(minusonelogit_powerparam(i)));
        }

    // Set up latent variables
    for(int i=0; i<num_units; i++) { for(int k=0; k<num_lv; k++) {
        nll -= dnorm(lvs(i,k), Type(0.0), Type(1.0), true);
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
        nll -= dtweedie(y(i,j), exp(eta_X(i,j)+eta_LV(i,j)), dispparam(j), powerparam(j), true);
        } }

    ADREPORT(betas); 
    ADREPORT(loadings_mat);
    ADREPORT(lvs);

    return nll; 
    }
