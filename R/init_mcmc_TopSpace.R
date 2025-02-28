#' Initialize the MCMC algorithm for TopSpace
#'
#' @param X 
#' @param nbd.sizes 
#' @param H 
#' @param K 
#' @param coord 
#' @param kernel.params 
#' @param poly.degree 
#'
#' @return
#' @export
#'
#' @examples
init.mcmc.TopSpace = function(X, nbd.sizes, H, K, coord, kernel.params = c(0.01, 1), poly.degree = 20L)
{
   # initialize Z, beta, alpha using the vanilla LDA
   out_lda = ldamcmc::lda_acgs(K = K, V = H, wid = as.matrix(unlist(X)), doc.N = as.matrix(nbd.sizes), alpha.v = rep(0.5, K), eta = 0.5, 
                               max.iter = 20000, burn.in = 10000, save.z = 1, save.beta = 1, save.theta = 1)
   
   Z_ppm_lda = matrix(NA, sum(nbd.sizes), K)
   for (k in 1 : K)
   {
      Z_ppm_lda[ , k] = rowMeans(out_lda$Z == k)
   }
   
   Z_lda = apply(Z_ppm_lda, 1, which.max)
   Z_lda_list = list()
   Z_lda_list[[1]] = Z_lda[1 : nbd.sizes[1]]
   for (i in 2 : length(X))
   {
      Z_lda_list[[i]] = Z_lda[sum(nbd.sizes[1 : (i - 1)]) + (1 : nbd.sizes[i])]
   }
   
   beta_lda = apply(out_lda$beta, c(1, 2), mean)
   
   alpha_lda = t(apply(out_lda$theta, c(1, 2), mean))
   
   # perform kernel decomposition for the GP with the modified exponential squared correlation kernel
   coord0 = BayesGPfit::GP.std.grids(coord, max_range = 1)
   Psi    = BayesGPfit::GP.eigen.funcs.fast(coord0, poly_degree = poly.degree, a = kernel.params[1], b = kernel.params[2])
   lambda = BayesGPfit::GP.eigen.value(poly_degree = poly.degree, a = kernel.params[1], b = kernel.params[2], d = ncol(coord0))
   Bases  = t(Psi) * sqrt(lambda)
   L      = ncol(Psi)
   
   # return initial values for the MCMC algorithm for GPSTM
   init = list()
   init$Z     = Z_lda_list
   init$beta  = beta_lda
   init$alpha = alpha_lda
   init$B       = matrix(0, L, K)
   init$n.eigen = L
   init$lambda  = lambda
   init$Bases   = Bases
   init$theta   = t(exp(crossprod(init$B, Bases)))
   return(init)
}