if(!require(VarianceGamma)){install.packages("VarianceGamma")}; library(VarianceGamma)
if(!require(e1071)){install.packages("e1071")}; library(e1071)

centralized_moments = function(data)
{
  m1 = mean(data)
  m2 = var(data)
  d1 = skewness(data)
  d2 = kurtosis(data)
  
  return(c(m1, m2, d1, d2))
}

vg_mom_estimates = function(data)
{
  moments = centralized_moments(data)
  m1 = moments[1]; m2 = moments[2]; d1 = moments[3]; d2 = moments[4]
  
  sigma = sqrt(m2)
  nu = (d2-3)/(3)
  theta = d1*m2^(1/2)/(d2-3)
  mu = m1-theta
  
  return(c(mu, sigma, theta, nu))
}

run_one = function(n, true_param) 
{
  start = Sys.time()
  x = rvg(n, param = true_param)
  
  # Fit VG with MOM
  est = tryCatch({
    vg_mom_estimates(x)
  }, error = function(e) NULL)
  
  runtime = as.numeric(Sys.time() - start, units = "secs")
  
  if (is.null(est)) {
    return(list(est = rep(NA, 4), runtime = runtime))
  }
  
  return(list(est = est, runtime = runtime))
}

sim_study = function(R = 200, ns = c(100, 500, 1000),
                      true_param = c(0, 1, 0.5, 0.8)) {
  results = list()
  
  for (n in ns) {
    est_mat = matrix(NA, nrow = R, ncol = 4)
    runtimes = numeric(R)
    
    for (r in 1:R) {
      sim = run_one(n, true_param)
      est_mat[r, ] = sim$est
      runtimes[r] = sim$runtime
    }
    
    colnames(est_mat) = c("mu", "sigma", "nu", "theta")
    
    # Performance metrics
    biases = colMeans(est_mat, na.rm = TRUE) - true_param
    rmse   = sqrt(colMeans((t(t(est_mat) - true_param))^2, na.rm = TRUE))
    avg_runtime = mean(runtimes, na.rm = TRUE)
    
    results[[as.character(n)]] = list(
      estimates = est_mat,
      bias = biases,
      rmse = rmse,
      runtime = avg_runtime
    )
  }
  
  return(results)
}

# --- Run the study ---
true_theta = .5; true_sigma = 1; true_nu = .8; true_mu = 0
true_param = c(true_mu, true_sigma, true_theta, true_nu)
ns = c(100, 500, 1000)

set.seed(123)
res = sim_study(R = 200, ns = ns, true_param = true_param)

print(res[["500"]]$bias)
print(res[["500"]]$rmse)
print(res[["500"]]$runtime)


