# Compute covariance matrix and coefficient vectors for data generation
# PARAMS:
# n: Number of observations we want
# m: Number of coefficients. p=Bm, and this model only runs when p>n.
# B: Number of "blocks" that will form the block-diagonal Sigma matrix
# n_sparse: Number of non-zero components of delta_1
# b2: b^2
# RETURN:
# List of 3 items: Gamma: sqrt of covariance matrix
#                  beta_1, beta_2: coef vectors for the two settings
setting <- function(n, m, B, n_sparse, b2, seed){
  set.seed(seed)
  p<-B*m 
  # Compute D: diagonal matrix in the decomposition of Sigma
  s<-ceiling(n^0.8)
  ## define weights as described in paper
  w = numeric(p)
  w[(s+1):p] = (1:(p-s))^-4
  W<-sum(w)
  ## Start as a vector, easier to use logic with
  ## then compute D from the weights and s
  ## finally take the sqrt as we are working with
  ## Gamma = O D^{1/2} O^T 
  diag.vec<-rep(1, p) 
  diag.vec[(s+1):p]<-((n-s)*w/W)[(s+1):p]
  D.sqrt <-diag(sqrt(diag.vec), nrow=p)
  
  # Compute O matrix
  ## If B = 1, then we do not have a blockwise matrix
  ## so generate p x p matrix from standard normal
  ## then make it orthogonal
  if(B==1) { 
    temp = matrix(rnorm(p^2,0,1),p,p)
    temp.sym = temp + t(temp)
    O = as.matrix(eigen(temp.sym, TRUE)$vectors)
    ## otherwise, generate a block-wise diagonal matrix
    ## with each block being mxm from standard normal
    ## then make it orthogonal
  } else {
    temp_full = matrix(0,p,p)
    for(i in 1:B){
      temp_sub = matrix(rnorm(m^2,0,1),m,m)
      temp_full[((i-1)*m+1):(i*m), ((i-1)*m+1):(i*m)] = temp_sub
    }
    temp.sym = temp_full + t(temp_full)
    O = as.matrix(eigen(temp.sym, TRUE)$vectors)
  }
  
  # Generate Gamma which is Gamma = O D^{1/2} O^T = Sigma^{1/2}
  Gamma = O%*%D.sqrt%*%t(O)
  Sigma = Gamma %*% t(Gamma)
  
  # Generate beta for both the sparse case (delta_1) and
  # the random case (delta_2). Denote beta_i where i = delta_i.
  ## For delta_1: Generate delta where n_sparse of them = 1,
  ## rest = 0.
  delta_1 = rep(0, p)
  delta_1[sample(1:p, n_sparse, replace=FALSE)] = rep(1, n_sparse)
  ## For delta_2: take a linear combination of the first 100 
  ## columns of O and utilize those as delta.
  delta_2 = as.vector(O[,1:(min(p, 100))] %*% rnorm(min(p, 100)))
  ## Generate beta = b*delta/sqrt{delta^T Sigma delta}
  b = b2^0.5
  ### beta_1
  denom_1 = as.numeric(sqrt(t(delta_1) %*% Sigma %*% delta_1))
  beta_1 = b*delta_1 / denom_1
  ### beta_2
  denom_2 = as.numeric(sqrt(t(delta_2) %*% Sigma %*% delta_2))
  beta_2 = b*delta_2 / denom_2
  
  # Return Gamma and betas
  return(list(Gamma = Gamma, beta_1 = beta_1, beta_2 = beta_2))
}

# Generate data following the Poisson model
## INPUTS: n = number of observations
##         m = such that p = Bm
##         p = number of parameters
##         Gamma = sqrt of covariate matrix Sigma^{1/2} = O D^{1/2} O^T
##         beta = coefficients of the model
##         data_gen = data generation method (either normal or uniform)
data.poisson = function(n, m, p, Gamma, beta, data_gen = 'normal') {
  # Set up matrix for data generation
  Z = NULL
  # Get data generation method string in lower case
  data_gen = tolower(data_gen)
  # Generate data based on the method specified
  ## - normal: Z ~ N(0,1)
  ## - uniform: Z ~ Uniform[-√3, √3]
  ## - otherwise: print error message and break
  if (data_gen == 'normal') {
    Z = matrix(rnorm(n*p, 0, 1), n, p)
  } else if (data_gen == 'uniform') {
    Z = matrix(runif(n*p, -sqrt(3), sqrt(3)), n, p)
  }
  else {
    print("ERROR: Data Generation Method not specified properly")
    break
  }
  # Generate X as Sigma^0.5 * Z
  X = Z %*% Gamma
  # Generate Y from Poisson(mu)
  mu = as.vector(exp(X %*% beta))
  Y = rpois(n, mu)
  
  # Return X and Y
  return(list(X=X, Y=Y))
}
