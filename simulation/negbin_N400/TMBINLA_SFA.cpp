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
    DATA_MATRIX(OLSmatrix_transpose);        
    DATA_MATRIX(residualprojection);        
    DATA_INTEGER(num_lv);
    DATA_FACTOR(meshidxloc);	// Pointers into random effects vector x
    DATA_INTEGER(n_i); // Number of vertices in mesh
    int num_units = y.rows(); 
    int num_spp = y.cols();
    int num_X = X.cols();
    
    // Declare parameters
    PARAMETER_MATRIX(betas);
    PARAMETER_VECTOR(log_dispparam);
    PARAMETER_VECTOR(loadings);
    PARAMETER_MATRIX(lvs);
    PARAMETER_VECTOR(log_kappa);
    
    // SPDE objects
    DATA_SPARSE_MATRIX(G0);
    DATA_SPARSE_MATRIX(G1);
    DATA_SPARSE_MATRIX(G2);

    // Negative log-likelihood and other parameter set up
    Type nll = Type(0.0); 
    vector<Type> kappa(num_lv);
    vector<Type> kappa2(num_lv);
    vector<Type> kappa4(num_lv);
    vector<Type> Range_raw(num_lv);
    vector<Type> dispparam = exp(log_dispparam);

    // Set up latent variables
    for(int k=0; k<num_lv; k++) {
        kappa(k) = exp(log_kappa(k));
        kappa2(k) = exp(2.0*log_kappa(k));
        kappa4(k) = kappa2(k)*kappa2(k);
        Range_raw(k) = sqrt(8) / exp(log_kappa(k));
        }

    vector<Type> log_nu(num_lv);
    matrix<Type> lvs_scale(n_i, num_lv);
    Eigen::SparseMatrix<Type> Q;
    for(int k=0; k<num_lv; k++) {
        Q = kappa4(k)*G0 + Type(2.0)*kappa2(k)*G1 + G2;
        nll += GMRF(Q)(lvs.col(k));
        log_nu(k) = log(1 / (exp(log_kappa(k)) * sqrt(4*3.141592)) );  // Ensures that MargSD = 1
        for(int k2=0; k2<n_i; k2++) 
            lvs_scale(k2,k) = lvs(k2,k) / exp(log_nu(k));
        }

    //To create loadings as upper triangular matrix
    matrix<Type> loadings_mat(num_spp, num_lv); 
    int loadings_ind=0;

    for (int j=0; j<num_lv; j++) { for (int i=0; i<num_spp; i++) {
        if (j > i) { loadings_mat(i,j) = 0; } 
        if (j <= i) { loadings_mat(i,j) = loadings(loadings_ind); loadings_ind += 1; } 
        } }

    // Data likelihood
    matrix<Type> eta_X = X*betas.transpose();
    matrix<Type> lvs_units(num_units, num_lv);
    for (int i=0; i<num_units; i++) { for(int k=0; k<num_lv; k++) {
            lvs_units(i,k) = lvs_scale(meshidxloc(i),k);
            } }        
    matrix<Type> eta_LV = lvs_units*loadings_mat.transpose();
    for (int i=0; i<num_units; i++) { for (int j=0; j<num_spp; j++) {
        nll -= dnbinom2(y(i,j), exp(eta_X(i,j)+eta_LV(i,j)), exp(eta_X(i,j)+eta_LV(i,j)) + dispparam(j)*pow(exp(eta_X(i,j)+eta_LV(i,j)),2), true);
        } }


    matrix<Type> betas_rsr = betas + eta_LV.transpose() * OLSmatrix_transpose;    
    matrix<Type> lvs_rsr = residualprojection * lvs_units;

    ADREPORT(betas); 
    ADREPORT(loadings_mat);
    ADREPORT(lvs_units);
    ADREPORT(betas_rsr); 
    ADREPORT(lvs_rsr);

    return nll; 
    }
