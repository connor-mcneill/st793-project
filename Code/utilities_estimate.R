score_fn = function(Y, X, beta) {
  score = matrix(numeric(ncol(X)), nrow=1, ncol=ncol(X))
  for (i in 1:length(Y)) {
    score = score + as.vector((Y[i] - exp(X[i,]%*%beta)))*X[i,]
  }
  return(score)
}

## RP test (code modified from supplemental material, comments are our own)
## Additional tests added for LRT, Score, and Wald based on RP data
##  X:  n x p matrix of covariates.
##  Y:  n dimensional vector.
##  rho: projection ratio.
##  D: random projection times.
RP.tests<-function(X,Y,rho,D=10){
  # Get number of rows n, columns p, and our projection dim k
  n = nrow(X)
  p = ncol(X)
  k = ceiling(rho*n)
  # Obtain I - P_1
  I.P1 = diag(n) - (1/n)*matrix(1,n,n)
  # Obtain P_k for each random projection, where P_k is a
  # random projection matrix with random entries, drawn
  # independently of the data. Perform D times and take 
  # mean of each time.
  Pk.array = array(rnorm((p*k*D),0,1),dim=c(p,k,D))
  Pk=(1/sqrt(p))*as.matrix(apply(Pk.array,c(1,2),mean))
  # Obtain the U matrix, (I - P1)*X*Pk
  Uk = I.P1 %*% X %*% Pk
  # Obtain the H matrix, the Hat matrix
  Hk = Uk %*% solve(t(Uk) %*% Uk) %*% t(Uk)
  # Compute the test statistic
  T.multi = ((t(Y)%*%Hk%*%Y)/k)/((t(Y)%*%( I.P1-Hk)%*%Y)/(n-k-1))
  # Return the standardized test statistic as given in 
  # Section 4.1
  deno = sqrt(2/(n*rho*(1-rho)))
  multi = as.numeric((T.multi-1)/deno)
  
  # Now utilize U_k as X to do random-project classical tests
  ## estimate the betas via MLE
  fit = glm(Y ~ Uk + 0, family = poisson(link = 'log'))
  beta_hat = as.vector(fit$coefficients)
  ## Wald = (beta_hat)^T {n I(beta_hat)} beta_hat
  W = diag(as.vector(exp(Uk %*% beta_hat)))
  obs_FI = (1/n)* t(Uk) %*% W %*% Uk
  Wald.stat = as.numeric(t(beta_hat) %*% (n * obs_FI) %*% beta_hat)
  ## Score = Score(beta_0)^T {n I(beta_0)}^{-1} Score(beta_0)
  score = score_fn(Y, Uk, rep(0,k))
  null_FI = (1/n)* t(Uk) %*% Uk
  Score.stat = as.numeric(score %*% solve(n * null_FI) %*% t(score))
  ## LRT is the deviance test
  LRT.stat = fit$null.deviance - fit$deviance
  
  # Return the four test statistics
  return(c(RP = multi,
           LRT = LRT.stat,
           Wald = Wald.stat,
           Score = Score.stat))
}

simulation = function(n, p, m, beta, Gamma, data.gen, rho, index,
                      #g.assume = NA, g.assume.deri = NA, V.assume = NA,
                      #RP.Flag = TRUE, Classical.Flag = TRUE, 
                      #CSX.Flag = FALSE,
                      sign = 0.95, L = 1000) {
  # obtain projection dimension k
  k = ceiling(rho * n)
  # initialize the storage of estimators
  estimator_storage = matrix(nrow = L, ncol = 4)
  # set up progress bar
  pb = txtProgressBar(min = 0, max = L, style=3)
  # Perform simulation
  for (i in 1:L) {
    # Set the seed
    set.seed(i)
    # Generate the data set from the given distribution
    data = data.gen(n, m, p, Gamma, beta)
    # Get the X and Y from the data
    X = data$X; Y = data$Y
    # Compute the test statistics and record them
    test_stats = RP.tests(X, Y, rho)
    estimator_storage[i,] = test_stats
    ## Update progress bar
    setTxtProgressBar(pb, i)
  }
  # Close progress bar
  close(pb)
  # Compute rejection rate of each estimator
  clas_rej_rate = apply(estimator_storage[,2:4], 2, 
                        function(x) {mean(ifelse(x > qchisq(sign,k), T, F))})
  rp_rej_rate = mean(ifelse(estimator_storage[,1] > qnorm(sign), T, F))
  
  # Return rejection rates
  return(c(rp_rej_rate, clas_rej_rate))
}