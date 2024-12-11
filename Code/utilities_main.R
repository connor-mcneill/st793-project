require(MASS)

# Initialize n, m, and B
n = 400
m = 10
B = 100
p = B*m
# SETTING 1: null hypothesis is true
set_null = setting(n, m, B, n_sparse = 10, b2=0, seed=793)
results_null = simulation(n, p, m, beta = set_null$beta_1, set_null$Gamma,
                          data.gen = data.poisson, rho = 0.4,
                          index = 0)
# SETTING 2: b^2 = 0.1
set_0.1 = setting(n, m, B, n_sparse = 10, b2 = 0.1, seed=793)
results_d1_0.1 = simulation(n, p, m, beta = set_0.1$beta_1,
                            set_0.1$Gamma, data.gen = data.poisson,
                            rho = 0.4, index=0)
results_d2_0.1 = simulation(n, p, m, beta = set_0.1$beta_2,
                            set_0.1$Gamma, data.gen = data.poisson,
                            rho = 0.4, index=0)
# SETTING 3: b^2 = 0.2
set_0.2 = setting(n, m, B, n_sparse = 10, b2 = 0.2, seed=793)
results_d1_0.2 = simulation(n, p, m, beta = set_0.2$beta_1,
                            set_0.2$Gamma, data.gen = data.poisson,
                            rho = 0.4, index=0)
results_d2_0.2 = simulation(n, p, m, beta = set_0.2$beta_2,
                            set_0.2$Gamma, data.gen = data.poisson,
                            rho = 0.4, index=0)
# Display results
results.table = rbind(results_null, results_d1_0.1, results_d2_0.1,
                      results_d1_0.2, results_d2_0.2)
print(results.table, digits=3)
# Save image of results in R
save.image(file = 'simresults.RData')