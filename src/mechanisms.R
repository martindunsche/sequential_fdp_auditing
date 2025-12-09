make_sum_gauss <- function(sigma) {
  function(x) sum(x) + rnorm(1, 0, sigma)
}

make_sum_laplace <- function(sigma) {
  function(x) {
    sum(x) + rmutil::rlaplace(1, m = 0, s = sigma)
  }
}

make_noisy_sgd <- function(theta_0, eta_learn, sigma, T_, m) {
  function(x){
    theta = theta_0
    for (i in 1:T_){
      x_subsample = sample(x, size = m)
      theta = theta - eta_learn * (sum(theta - x_subsample) / m +
                                   rnorm(1, 0, sigma))
    }
    theta
  }
}
