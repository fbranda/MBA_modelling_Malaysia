crp_gibbs <- function(data,alpha=0.01,mu0,sigma0,sigma_y,c_init,maxIters=1000)
 {
data_dim <- ncol(data)
N <- nrow(data)
tau0 <- solve(sigma0)
tau_y <- solve(sigma_y)      
z <- c_init

n_k <- as.vector(table(z)) 
Nclust <- length(n_k) 

res <- matrix(NA, nrow = N, ncol = maxIters) 
pb <- txtProgressBar(min = 0, max = maxIters, style = 3)
for(iter in 1:maxIters) { # maxIters also prevents endless loops
 for(n in 1:N) { # one data point (customer) at a time
   
    c_i <- z[n] # what is the nth persons table assignment?
    n_k[c_i] <- n_k[c_i] - 1 # remove the nth person from table
    
    if( n_k[c_i] == 0 )
       {
         n_k[c_i] <- n_k[Nclust] # last cluster to replace this empty cluster
         loc_z <- ( z == Nclust ) # who are in the last cluster?
         z[loc_z] <- c_i # move them up to fill just emptied cluster
         n_k <- n_k[ -Nclust ] # take out the last cluster, now empty
         Nclust <- Nclust - 1 # decrease total number of clusters by 1
       }
  z[n] <-  -1 
  logp <- rep( NA, Nclust + 1 )
   
  for( c_i in 1:Nclust ) {
   tau_p <- tau0 + n_k[c_i] * tau_y # cluster precision as per Eq (4)
   sig_p <- solve(tau_p) 
   loc_z <- which(z == c_i)
     if(length(loc_z) > 1) {
       sum_data <-  colSums(data[z == c_i, , drop = FALSE]) }
      else {
       sum_data <- as.numeric(data[z == c_i, , drop = FALSE])
      }
    sum_data <- matrix(sum_data, ncol = 1)
    mean_p <- sig_p %*% (tau_y %*% sum_data + tau0 %*% t(mu0))
    logp[c_i] <- log(n_k[c_i]) +
    dmvnorm(data[n,], mean = mean_p, sigma = sig_p + sigma_y, log = TRUE) }
  
  logp[ Nclust+1 ] <- log(alpha) +
 dmvnorm(data[n,], mean = mu0, sigma = sigma0 + sigma_y, log = TRUE)
 max_logp <- max(logp)
  logp <- logp - max_logp
loc_probs <- exp(logp)
 loc_probs <- loc_probs / sum(loc_probs)
newz <- sample(1:(Nclust+1), 1, replace = TRUE, prob = loc_probs)
 if(newz == Nclust + 1) {
  n_k <- c(n_k, 0)
 Nclust <- Nclust + 1
  }
z[n] <- newz
 n_k[newz] <- n_k[newz] + 1 # update the cluster n_k
 
}
 setTxtProgressBar(pb, iter) # update text progress bar after each iter
 
 res[, iter] <- z # cluster membership of N observations
}

close(pb) # close text progress bar
invisible(res) # return results, N by maxIters matrix

}
      