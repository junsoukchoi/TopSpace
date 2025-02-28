#' An efficient implementation of Markov chain Monte Carlo for TopSpace
#'
#' @param X 
#' @param nbd.sizes 
#' @param H 
#' @param K 
#' @param init 
#' @param priors 
#' @param n.mcmc 
#' @param n.burnin 
#' @param step.show 
#' @param learning.rate 
#' @param decay 
#' @param n.leapfrog 
#' @param prop.subsample 
#' @param thres.merge 
#'
#' @return
#' @export
#'
#' @examples
mcmc.TopSpace = function(X, nbd.sizes, H, K, init, priors, n.mcmc, n.burnin, step.show = 100, learning.rate = 0.01, decay = 0.01, n.leapfrog = 10, prop.subsample = 0.5, thres.merge = 0.01)
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
      out_update_B = update_B(B, nu, alpha, Bases, lambda, n.leapfrog, floor(prop.subsample * M), decay, learning.rate / M)
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
