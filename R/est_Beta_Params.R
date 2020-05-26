

## Estimates starting parameters for the Beta distribution calculation
##
## @param mu The mean of the data.
## @param var The variance of the data.
##
## @return A list containing the estimated beta parameters alpha and beta.
##
est_Beta_Params <- function(mu, var){
  alpha <- ((1-mu)/var - 1/mu) * mu^2
  beta <- alpha * (1/mu - 1)
  return(params = list(alpha = alpha, beta = beta))
}
