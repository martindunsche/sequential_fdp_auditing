
library(KernSmooth)


inner_int <- function(x, eta, grid_width, p, log_lr) {
  thresh <- eta + x
  if (thresh <= 0) return(sum(p * grid_width))
  log_thresh <- log(thresh)
  p_vec <- p
  p_vec[log_lr <= log_thresh] <- 0
  sum(p_vec * grid_width)
}

alpha_value <- function(eta, x_vector, h, grid_width, hat_p, log_lr) {
  dx <- x_vector[2] - x_vector[1]
  val <- sum(
    sapply(x_vector, inner_int,
           eta = eta,
           grid_width = grid_width,
           p = hat_p,
           log_lr = log_lr) * dx / h
  )
  1 - val
}

beta_value <- function(eta, x_vector, h, grid_width, hat_q, log_lr) {
  dx <- x_vector[2] - x_vector[1]
  sum(
    sapply(x_vector, inner_int,
           eta = eta,
           grid_width = grid_width,
           p = hat_q,
           log_lr = log_lr) * dx / h
  )
}


KDE_Estimator <- function(eta_max, Mechanism, x1, x2, N, h){
  
  p_sample <- sapply(rep(list(x1),N), Mechanism)
  q_sample <- sapply(rep(list(x2),N), Mechanism)
  
  t_min = min(c(p_sample,q_sample))
  t_max = max(c(p_sample,q_sample))
  
  bw_p = dpik(p_sample, kernel = 'normal', gridsize = min(1000L,N/10), range.x = c(t_min,t_max))
  bw_q = dpik(q_sample, kernel = 'normal', gridsize = min(1000L,N/10), range.x = c(t_min,t_max))
  
  KDE_p <- bkde(p_sample, kernel = 'normal', bandwidth = bw_p, gridsize = min(1000L,N/10), range.x = c(t_min,t_max))
  KDE_q <- bkde(q_sample, kernel = 'normal', bandwidth = bw_q, gridsize = min(1000L,N/10), range.x = c(t_min,t_max))
  
  hat_p <- KDE_p$y
  hat_q <- KDE_q$y
  grid_width <- KDE_p$x[2] - KDE_p$x[1]
  
  hat_p_q_ratio <- hat_p / hat_q
  hat_p_q_ratio[is.infinite(hat_p_q_ratio)] <- 0
  hat_p_q_ratio[is.na(hat_p_q_ratio)] <- 0
  
  
  x_vector = seq(from = - h/2, to = h/2, length.out = 1000)
  eta_vector = seq(from = 0, to = eta_max, length.out = 1000)
  alpha = sapply(eta_vector, alpha_value, x_vector,h,grid_width,hat_p,hat_p_q_ratio)
  beta = sapply(eta_vector, beta_value, x_vector,h,grid_width,hat_q,hat_p_q_ratio)
  xd<- data.frame(alpha,beta,eta_vector)
  return(list(p=p_sample,q=q_sample,xd=xd, hat_p=KDE_p, hat_q=KDE_q, hat_p_q_ratio=hat_p_q_ratio))
}













