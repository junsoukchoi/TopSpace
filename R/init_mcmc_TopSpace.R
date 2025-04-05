#' Initialize the MCMC algorithm for TopSpace
#'
#' @param X A list of vectors, each representing a local neighborhood and specifying the cell types observed within that neighborhood.
#' @param nbd.sizes A vector indicating the number of cells present in each local neighborhood.
#' @param H An integer indicating the total number of cell types observed in the data.
#' @param K An integer indicating the number of topics.
#' @param coord A matrix whose rows correspond to the spatial coordinates of each local neighborhood.
#' @param kernel.params A two-element vector specifying the hyperparameters of the modified exponential squared covariance kernel for the GP over spatially varying Dirichlet hyperparameters. The first element is the concentration parameter, where higher values yield greater concentration around the center. The second element is the smoothness parameter, where lower values result in a smoother process. The default values are (0.01, 1).
#' @param poly.degree An integer specifying the highest degree of Hermite polynomials that determines the dimension (length) of the truncated basis expansion for the GP over spatially varying Dirichlet hyperparameters. The default value is 10L.
#'
#' @return A list of initial values for the MCMC algorithm for TopSpace. Valid tags include 'Z', 'beta', 'alpha', 'B', 'n.eigen', 'lambda', 'Bases', and 'theta'. The value corresponding to each tag contains either the initial MCMC values or the eigendecomposition specification for the GP in TopSpace.
#' 
#' @export
#'
#' @examples # set random seed
#' set.seed(2024)
#' 
#' # load the pre-processed non-small cell lung cancer (NSCLC) dataset
#' data("NSCLC")
#' 
#' # select the image of the second patient
#' img = 2
#' X = NSCLC[[img]]$X
#' nbd_sizes = NSCLC[[img]]$nbd.sizes
#' coord = NSCLC[[img]]$coord
#' 
#' # specify the total number of cell types observed in the data: 
#' # CD19+ B-cells, CD14+ cells, CD8+ T-cells, CD4+ T-cells, CK+ cancer cells, and Others
#' H = 6
#' 
#' # specify the number of topics to be considered
#' K = 4
#'    
#' # obtain initial values for the MCMC algorithm for TopSpace
#' init = init.mcmc.TopSpace(X, nbd_sizes, H, K, coord, kernel.params =  c(0.01, 1))
#' print(init)
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
   
   # perform the eigendecomposition of the modified exponential squared covariance kernel
   coord0 = BayesGPfit::GP.std.grids(coord, max_range = 1)
   Psi    = BayesGPfit::GP.eigen.funcs.fast(coord0, poly_degree = poly.degree, a = kernel.params[1], b = kernel.params[2])
   lambda = BayesGPfit::GP.eigen.value(poly_degree = poly.degree, a = kernel.params[1], b = kernel.params[2], d = ncol(coord0))
   Bases  = t(Psi) * sqrt(lambda)
   L      = ncol(Psi)
   
   # return initial values for the MCMC algorithm for TopSpace
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