#' An efficient implementation of Markov chain Monte Carlo for TopSpace
#'
#' @param X A list of vectors, each representing a local neighborhood and specifying the cell types observed within that neighborhood.
#' @param nbd.sizes A vector indicating the number of cells present in each local neighborhood.
#' @param H An integer indicating the total number of cell types observed in the data.
#' @param K An integer indicating the number of topics.
#' @param init A list containing initial values for the MCMC algorithm for TopSpace. Valid tags include 'Z', 'beta', 'alpha', 'B', 'n.eigen', 'lambda', 'Bases', and 'theta'. The value corresponding to each tag should be either the initial MCMC values or the eigendecomposition specification for the GP in TopSpace.
#' @param priors A list of hyperparameter values. The valid tag is 'gamma', which specifies the Dirichlet prior hyperparameters for the per-topic cell type distributions \eqn{\beta_k}.
#' @param n.mcmc An integer specifying the total number of MCMC iterations.
#' @param n.burnin An integer specifying the number of burn-in samples.
#' @param step.show An integer indicating the interval to report MCMC progress.
#' @param learning.rate A double specifying the learning rate for the SGHMC update applied to the GP over spatially varying Dirichlet hyperparameters.
#' @param friction A double representing the friction term used for the SGHMC update for the GP over spatially varying Dirichlet hyperparameters.
#' @param n.leapfrog An integer specifying the number of leapfrog steps used in the SGHMC update for the GP over spatially varying Dirichlet hyperparameters.
#' @param prop.subsample A double between 0 and 1 indicating the proportion of local neighborhoods to subsample during the SGHMC update for the GP over spatially varying Dirichlet hyperparameters.
#' @param thres.merge A double specifying the threshold for merging of similar topics---specifically, topics \eqn{k} and \eqn{k'} are merged if \eqn{|| \beta_k - \beta_{k'} ||_2 \le} \code{thresh.merge}.
#'
#' @return A list of posterior samples from the MCMC algorithm for TopSpace.  Valid tags include 'Z', 'beta', 'alpha', 'K.eff', and 'DIC', each containing MCMC samples for the corresponding parameter: latent topic assignments (Z), per-topic cell type distributions (beta), and per-neighborhood topic distributions (alpha). The list also includes the effective number of topics (K.eff) and the Deviance Information Criterion (DIC).
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
#' 
#' ### fit our proposed TopSpace to the selected image ###
#' # will choose the optimal number of topics K among {2, 3, 4, 5}
#' DIC   = rep(NA, 4)
#' K_eff = rep(NA, 4)
#' TIME  = rep(NA, 4)
#' out   = list()
#' for (K in 2 : 5)
#' {
#'    # number of cell types: CD19+ B-cells, CD14+ cells, CD8+ T-cells, CD4+ T-cells, CK+ cancer cells, and Others
#'    H = 6
#'    
#'    # initialize MCMC for TopSpace using the standard LDA
#'    init   = init.mcmc.TopSpace(X, nbd_sizes, H, K, coord, kernel.params =  c(0.01, 1))
#'    priors = list(gamma = rep(0.5, H))
#'    
#'    # run MCMC for the proposed TopSpace, with the specified K
#'    p_time = proc.time()
#'    out_GPSTM  = mcmc.TopSpace(X, nbd_sizes, H, K, init, priors, n.mcmc = 20000, n.burnin = 10000)
#'    time_GPSTM = proc.time() - p_time
#'    
#'    # store the deviance information criterion (DIC), effective number of topics K, and runtime
#'    DIC[K - 1]   = out_GPSTM$DIC
#'    K_eff[K - 1] = out_GPSTM$K.eff
#'    TIME[K - 1]  = time_GPSTM[3]
#'    
#'    # summarize posterior samples: compute posterior means for Z, beta, and alpha
#'    Z_ppm = matrix(NA, sum(nbd_sizes), K)
#'    for (k in 1 : K)
#'    {
#'       Z_ppm[ , k] = rowMeans(out_GPSTM$Z == k)
#'    }
#'    
#'    Z_est = apply(Z_ppm, 1, which.max)
#'    Z_est_list = list()
#'    Z_est_list[[1]] = Z_est[1 : nbd_sizes[1]]
#'    for (i in 2 : length(X))
#'    {
#'       Z_est_list[[i]] = Z_est[sum(nbd_sizes[1 : (i - 1)]) + (1 : nbd_sizes[i])]
#'    }
#'    
#'    beta_est = apply(out_GPSTM$beta, c(1, 2), mean)
#'    
#'    alpha_est = apply(out_GPSTM$alpha, c(1, 2), mean)
#'    
#'    # find the dominant topic for each neighborhood by selecting the highest topic probability
#'    nbd_topic_est = apply(alpha_est, 1, which.max)
#'    
#'    # store the posterior inference results
#'    out[[K - 1]] = list()
#'    out[[K - 1]]$Z = Z_est_list
#'    out[[K - 1]]$beta = beta_est
#'    out[[K - 1]]$alpha = alpha_est
#'    out[[K - 1]]$nbd_topic = nbd_topic_est
#' }
#' 
#' # choose the optimal K based on the DIC, and obtain the corresponding results
#' result = list()
#' result$dat = X
#' result$coord = coord
#' result$nbd.sizes = nbd_sizes
#' opt = which.min(DIC)
#' result$K    = K_eff[opt]
#' result$DIC  = DIC[opt]
#' result$TIME = TIME
#' result$est  = out[[opt]]
#' 
#' 
#' ### visualize the TopSpace results ###
#' library(ggplot2)
#' library(ggforce)
#' library(MBA)
#' library(fields)
#' 
#' # illustrate topics (beta), i.e., per-topic cell type distribution
#' x = paste("Topic", 1 : result$K)
#' y = c("CD19+ B-cell", "CD14+ cell", "CD8+ T-cell", "CD4+ T-cell", "CK+ cancer cell","Other")
#' df_beta = expand.grid(Topics = x, Phenotypes = y)
#' df_beta$beta = c(result$est$beta)
#' ggplot(data=df_beta, aes(x=Topics, y=beta, fill=Phenotypes)) + geom_bar(stat="identity", position = position_fill(reverse = TRUE)) + 
#'    xlab("") + ylab("Probability") +
#'    theme(legend.title = element_text(face = "bold"), legend.text = element_text(face = "bold"), 
#'          axis.title.x = element_text(face = "bold"), axis.title.y = element_text(face = "bold"), 
#'          axis.text.x = element_text(face = "bold"), axis.text.y = element_text(face = "bold"))
#' 
#' # illustrate spatial clustering defined by the dominant topics
#' df_clustering = data.frame(x = result$coord[ , 1], y = result$coord[ , 2], r = rep(19, nrow(result$coord)), 
#'                            Topic = factor(paste("Topic", apply(result$est$alpha, 1, which.max))))
#' ggplot(df_clustering, aes(x0 = x, y0 = y, r = r, color = Topic, fill = Topic)) + geom_circle() + scale_fill_brewer(palette = "Set2") +
#'    scale_color_brewer(palette = "Set2") + xlab("") + ylab("") + 
#'    theme(legend.title = element_text(face = "bold"), legend.position="top", legend.text = element_text(size = 12, face = "bold"),
#'          axis.title.x = element_text(face = "bold"), axis.title.y = element_text(face = "bold"), 
#'          axis.text.x = element_text(face = "bold"), axis.text.y = element_text(face = "bold"))
#' 
#' # generatea spatial probability map for tertiary lymphoid structures (TLS)
#' result$cl_Bcells = which(result$est$beta[ , 1] > 0.5) 
#' result$cl_Tcells = which(result$est$beta[ , 3] + result$est$beta[ , 4] > 0.5) 
#' surf = mba.surf(cbind(result$coord, data = rowSums(result$est$alpha[ , c(result$cl_Bcells, result$cl_Tcells)])), 
#'                 no.X = 200, no.Y = 200)$xyz.est
#' image.plot(surf, xlab = "", ylab ="", zlim = c(-0.1, 1.1), legend.mar = 6, main = "Spatial Probability Map of TLS")
mcmc.TopSpace = function(X, nbd.sizes, H, K, init, priors, n.mcmc, n.burnin, step.show = 100, learning.rate = 0.01, friction = 0.01, n.leapfrog = 10, prop.subsample = 0.5, thres.merge = 0.01)
{
   
   M = length(X)
   
   # initialize parameters and set the basis functions for the GP
   Z      = init$Z
   beta   = t(init$beta)
   alpha  = t(init$alpha)
   B      = init$B
   L      = init$n.eigen
   lambda = init$lambda
   Bases  = init$Bases
   theta  = t(init$theta)
   
   # initialize MCMC samples
   Z_mcmc        = matrix(NA, sum(nbd.sizes), n.mcmc)
   beta_mcmc     = array(NA, dim = c(K, H, n.mcmc))
   alpha_mcmc    = array(NA, dim = c(M, K, n.mcmc))
   #B_mcmc        = array(NA, dim = c(L, K, n.mcmc))
   #theta_mcmc    = array(NA, dim = c(M, K, n.mcmc))
   deviance_mcmc = rep(NA, n.mcmc - n.burnin)
   
   # calculate n_kh's and m_ik's
   n_kh = calculate_n_kh(X, Z, H, K)
   m_ik = calculate_m_ik(Z, K)
   
   # update the Markov chain
   iter = 1
   while(iter <= n.mcmc)
   {
      # update Z
      out_update_Z = update_Z(Z, X, theta, n_kh, m_ik, priors$gamma)
      Z    = out_update_Z$Z
      n_kh = out_update_Z$n_kh
      m_ik = out_update_Z$m_ik
      
      # update beta
      beta = update_beta(n_kh, priors$gamma)
      
      # update alpha
      alpha = update_alpha(m_ik, theta)
      
      if (iter %% 50 == 1)
         nu = rnorm(L * K) * sqrt(learning.rate / M)
      
      # update B
      out_update_B = update_B(B, nu, alpha, Bases, lambda, n.leapfrog, floor(prop.subsample * M), friction, learning.rate / M)
      B     = out_update_B$B
      nu    = out_update_B$nu
      theta = out_update_B$theta
      
      # save MCMC samples
      Z_mcmc[ , iter]       = unlist(Z)
      beta_mcmc[ , , iter]  = t(beta)
      alpha_mcmc[ , , iter] = t(alpha)
      #B_mcmc[ , , iter]     = B
      #theta_mcmc[ , , iter] = t(theta)
      
      # after burn-in period, calculate the deviance for the current parameters 
      if (iter > n.burnin)
      {
         deviance_mcmc[iter - n.burnin] = calculate_deviance(X, beta, alpha)
      }
      
      # display the current iteration count (at each step.show iteration)
      if (iter %% step.show == 0) cat("iter =", iter, "\n")
      
      # merge topics if they are identical
      if (iter %% 1000 == 0)
      {
         beta_temp = t(apply(beta_mcmc[ , , (iter - 999) : iter], c(1, 2), mean))
         do_merge  = TRUE
         while (do_merge)
         {
            diff_beta = as.matrix(dist(t(beta_temp), method = "euclidean"))
            diff_beta[lower.tri(diff_beta)] = Inf
            diag(diff_beta) = Inf
            ids_merge = which(diff_beta < thres.merge, arr.ind = TRUE)
            
            if (nrow(ids_merge) == 0)
            {
               do_merge = FALSE
            } 
            else
            {
               id_keep = ids_merge[1, 1]
               id_drop = ids_merge[1, 2]
               
               beta[ , id_keep] = (beta[ , id_keep] + beta[ , id_drop]) / 2
               beta = beta[ , -id_drop, drop = FALSE]
               
               alpha[id_keep, ] = alpha[id_keep, ] + alpha[id_drop, ]
               alpha = alpha[-id_drop, , drop = FALSE]
               
               n_kh[id_keep, ] = n_kh[id_keep, ] + n_kh[id_drop, ]
               n_kh = n_kh[-id_drop, , drop = FALSE]
               
               m_ik[ , id_keep] = m_ik[ , id_keep] + m_ik[ , id_drop]
               m_ik = m_ik[ , -id_drop, drop = FALSE]
               
               for (i in 1 : M)
               {
                  Z[[i]][Z[[i]] == id_drop] = id_keep
                  Z[[i]][Z[[i]] > id_drop] = Z[[i]][Z[[i]] > id_drop] - 1
               }
               
               beta_temp[ , id_keep] = (beta_temp[ , id_keep] + beta_temp[ , id_drop]) / 2
               beta_temp = beta_temp[ , -id_drop, drop = FALSE]
               
               K = ncol(beta)
               
               B = matrix(0, L, K)
               theta = exp(crossprod(B, Bases))
               
               Z_mcmc     = matrix(NA, sum(nbd.sizes), n.mcmc)
               beta_mcmc  = array(NA, dim = c(K, H, n.mcmc))
               alpha_mcmc = array(NA, dim = c(M, K, n.mcmc))
               #B_mcmc     = array(NA, dim = c(L, K, n.mcmc))
               #theta_mcmc = array(NA, dim = c(M, K, n.mcmc))
               
               iter = 0
            }
         }
      }
      
      iter = iter + 1
   }
   
   # return MCMC samples of Z, beta, alpha, B, and theta, along with deviances
   out = list()
   out$Z        = Z_mcmc[ , (n.burnin + 1) : n.mcmc]
   out$beta     = beta_mcmc[ , , (n.burnin + 1) : n.mcmc]
   out$alpha    = alpha_mcmc[ , , (n.burnin + 1) : n.mcmc]
   #out$B        = B_mcmc[ , , (n.burnin + 1) : n.mcmc]
   #out$theta    = theta_mcmc[ , , (n.burnin + 1) : n.mcmc]
   out$K.eff    = K
   out$DIC      = mean(deviance_mcmc) + 0.5 * var(deviance_mcmc)
   
   return(out)
}
