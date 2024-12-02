# Generate data following the Negative Binomial model
## INPUTS: n = number of observations
##         m = such that p = Bm
##         p = number of parameters
##         Gamma = sqrt of covariate matrix Sigma^{1/2} = O D^{1/2} O^T
##         beta = coefficients of the model
##         data_gen = data generation method (either normal or uniform)
data.nb = function(n, m, p, Gamma, beta, data_gen = 'normal') {
  # Set up matrix for data generation
  Z = NULL
  # Get data generation method string in lower case
  data_gen = tolower(data_gen)
  # Generate data based on the method specified
  ## - normal: Z ~ N(0,1)
  ## - uniform: Z ~ Uniform[-√3, √3]
  ## - otherwise: print error message and break
  if (data_gen == 'normal') {
    Z = matrix(rnorm(n*m, 0, 1), n, m)
  } else if (data_gen = 'uniform') {
    Z = matrix(runif(n*m, -sqrt(3), sqrt(3)), n, m)
  }
  else {
    print("ERROR: Data Generation Method not specified properly")
    break
  }
  # Generate X as Sigma^0.5 * Z
  X = Z %*% Gamma
}