#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::mat calculate_n_kh(const ListOf<IntegerVector>& X, const ListOf<IntegerVector>& Z, int H, int K) 
{
    
    int N = X.size();
    
    arma::mat n_kh(K, H, arma::fill::zeros);
    
    for (int h = 0; h < H; ++h) 
    {
        for (int k = 0; k < K; ++k) 
        {
            for (int i = 0; i < N; ++i) 
            {
                IntegerVector x_i = X[i];
                IntegerVector z_i = Z[i];
                
                for (int j = 0; j < x_i.size(); ++j) 
                {
                    n_kh(k, h) += (x_i[j] == h + 1) && (z_i[j] == k + 1);
                }
            }
        }
    }
    
    return n_kh;
}

// [[Rcpp::export]]
arma::mat calculate_m_ik(const ListOf<IntegerVector>& Z, int K) 
{
    
    int N = Z.size();
    
    arma::mat m_ik(N, K, arma::fill::zeros);
    
    for (int k = 0; k < K; ++k) 
    {
        for (int i = 0; i < N; ++i) 
        {
            IntegerVector z_i = Z[i];
            
            for (int j = 0; j < z_i.size(); ++j) 
            {
                m_ik(i, k) += (z_i[j] == k + 1);
            }
        }
    }
    
    return m_ik;
}

// [[Rcpp::export]]
List update_Z(const ListOf<IntegerVector>& Z, const ListOf<IntegerVector>& X, const arma::mat& theta, const arma::mat& n_kh, const arma::mat& m_ik, const arma::vec& gamma) 
{
    
    int N = Z.size();
    int K = n_kh.n_rows;
    int H = n_kh.n_cols;
    
    List Z_new = clone(Z);
    arma::mat n_kh_new(n_kh.begin(), K, H);
    arma::mat m_ik_new(m_ik.begin(), N, K);
    
    for (int i = 0; i < N; ++i) 
    {
        IntegerVector x_i = X[i];
        IntegerVector z_i = Z_new[i];
        int Mi = x_i.size();
        for (int j = 0; j < Mi; ++j) 
        {
            arma::mat n0_kh = n_kh_new;
            arma::mat m0_ik = m_ik_new;
            n0_kh(z_i[j] - 1, x_i[j] - 1) = n0_kh(z_i[j] - 1, x_i[j] - 1) - 1;
            m0_ik(i, z_i[j] - 1) = m0_ik(i, z_i[j] - 1) - 1;
            
            arma::vec prob_ij = (n0_kh.col(x_i[j] - 1) + gamma[x_i[j] - 1]) / (arma::sum(n0_kh, 1) + arma::sum(gamma)) %
                (m0_ik.row(i).t() + theta.col(i)) / arma::sum(m0_ik.row(i).t() + theta.col(i));
            
            double    prob_sum = arma::sum(prob_ij);
            arma::vec prob_normalized = prob_ij / prob_sum;
            
            //Rcpp::Rcout << prob_normalized << std::endl; // Uncomment for debugging
            
            z_i[j] = Rcpp::RcppArmadillo::sample(arma::linspace(1, K, K), 1, false, prob_normalized)[0];
            
            n_kh_new = n0_kh;
            m_ik_new = m0_ik;
            n_kh_new(z_i[j] - 1, x_i[j] - 1) = n_kh_new(z_i[j] - 1, x_i[j] - 1) + 1;
            m_ik_new(i, z_i[j] - 1) = m_ik_new(i, z_i[j] - 1) + 1;
        }
        
        Z_new[i] = z_i;
    }
    
    return List::create(Named("Z") = Z_new, Named("n_kh") = n_kh_new, Named("m_ik") = m_ik_new);
}

// [[Rcpp::export]]
arma::mat update_beta(const arma::mat& n_kh, const arma::vec& gamma) 
{
    
    int K = n_kh.n_rows;
    int H = n_kh.n_cols;
    
    arma::mat beta(H, K, arma::fill::zeros);
    
    for (int k = 0; k < K; ++k) 
    {
        for (int h = 0; h < H; ++h)
        {
            beta(h, k) = Rcpp::rgamma(1, n_kh(k, h) + gamma[h], 1.0)[0];
        }
        
        beta.col(k) /= arma::accu(beta.col(k));
        
        // Check and replace zeros
        if (arma::any(beta.col(k) < std::numeric_limits<double>::epsilon())) 
        {
            arma::uvec k_vec(1, arma::fill::value(k));
            beta.submat(arma::find(beta.col(k) < std::numeric_limits<double>::epsilon()), k_vec).fill(std::numeric_limits<double>::epsilon());
        }
    }
    
    return beta;
}

// [[Rcpp::export]]
arma::mat update_alpha(const arma::mat& m_ik, const arma::mat& theta) 
{
    
    int N = m_ik.n_rows;
    int K = m_ik.n_cols;
    
    arma::mat alpha(K, N, arma::fill::zeros);
    
    for (int i = 0; i < N; ++i) 
    {
        for (int k = 0; k < K; ++k)
        {
            alpha(k, i) = Rcpp::rgamma(1, m_ik(i, k) + theta(k, i), 1.0)[0];
        }
        
        alpha.col(i) /= arma::accu(alpha.col(i));
        
        // Check and replace zeros
        if (arma::any(alpha.col(i) < std::numeric_limits<double>::epsilon())) 
        {
            arma::uvec i_vec(1, arma::fill::value(i));
            alpha.submat(arma::find(alpha.col(i) < std::numeric_limits<double>::epsilon()), i_vec).fill(std::numeric_limits<double>::epsilon());
        }
    }
    
    return alpha;
}

// [[Rcpp::export]]
arma::mat dU_tilde(const arma::mat& B, const arma::mat& alpha, const arma::mat& theta, const arma::mat& Bases, const arma::vec& lambda, int n_subsample) 
{
    
    int N = theta.n_cols;
    int L = B.n_rows;
    int K = B.n_cols;
    
    arma::uvec ind_sub   = arma::randi<arma::uvec>(n_subsample, arma::distr_param(0, N - 1));
    arma::mat  theta_sub = theta.cols(ind_sub);
    arma::mat  alpha_sub = alpha.cols(ind_sub);
    arma::mat  Bases_sub_tr = Bases.cols(ind_sub).t();
    arma::rowvec  colsums_theta_sub = arma::sum(theta_sub, 0);
    
    arma::mat out_dU_sub(L, K, arma::fill::zeros);
    
    // Accessing the digamma function from the R environment
    Rcpp::Function digammaR("digamma");
    
    for (int k = 0; k < K; ++k) 
    {
        arma::rowvec theta_sub_k = theta_sub.row(k);
        
        // Using R's digamma function
        arma::rowvec digamma_diff = as<arma::rowvec>(digammaR(theta_sub_k)) - as<arma::rowvec>(digammaR(colsums_theta_sub)) - arma::log(alpha_sub.row(k)) + std::numeric_limits<double>::epsilon();
        
        arma::rowvec colsums_term = arma::sum(Bases_sub_tr.each_col() % (digamma_diff % theta_sub_k).t(), 0);
        out_dU_sub.col(k) = (N / n_subsample) * colsums_term.t() + B.col(k) / lambda(k);
    }
    
    return out_dU_sub;
}

// [[Rcpp::export]]
List update_B(const arma::mat& B, const arma::vec& nu, const arma::mat& alpha, const arma::mat& Bases, const arma::vec& lambda, int n_leapfrog, int n_subsample, double rho, double eta) 
{
    
    int N = Bases.n_cols;
    int L = B.n_rows;
    int K = B.n_cols;
    
    arma::mat B_new(B.begin(), L, K);
    arma::vec nu_new(nu.begin(), nu.size());
    arma::mat theta_new(K, N, arma::fill::zeros);
    
    for (int t = 0; t < n_leapfrog; ++t) 
    {
        B_new    += arma::reshape(nu_new, L, K); 
        theta_new = arma::exp(B_new.t() * Bases);
        arma::vec w  = arma::randn(L * K) * sqrt(2 * rho * eta);
        arma::mat dU = dU_tilde(B_new, alpha, theta_new, Bases, lambda, n_subsample);
        nu_new = (1 - rho) * nu_new - eta * vectorise(dU) + w;
    }
    
    return List::create(Named("B") = B_new, Named("nu") = nu_new, Named("theta") = theta_new);
}

// [[Rcpp::export]]
double calculate_deviance(const ListOf<IntegerVector>& X, const arma::mat& beta, const arma::mat& alpha) 
{
    
    double deviance = 0;
    
    for (int i = 0; i < X.size(); ++i) 
    {
        IntegerVector x_i = X[i];
        arma::vec     beta_times_alpha_i = sum(beta.each_row() % alpha.col(i).t(), 1);
        
        for (int j = 0; j < x_i.size(); ++j) 
        {
            deviance -= 2 * log(beta_times_alpha_i[x_i[j] - 1]);
        }
    }
    
    return deviance;
}
