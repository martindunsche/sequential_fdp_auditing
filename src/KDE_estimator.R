
library(KernSmooth)



inner_int <- function(x,eta,grid_width,p,p_q_ratio){
  p_vector <- p
  p_vector[which(p_q_ratio <= eta + x)] <- 0
  int_value <- sum(p_vector * grid_width)
  return(int_value)
}


alpha_value <- function(eta,x_vector,h,grid_width,hat_p,hat_p_q_ratio){
  alpha = 1 - sum(sapply(x_vector,inner_int, eta = eta, grid_width = grid_width,
                         p = hat_p, p_q_ratio = hat_p_q_ratio) / h * (x_vector[2] - x_vector[1]))
  return(alpha)
}

beta_value <- function(eta,x_vector,h,grid_width,hat_q,hat_p_q_ratio){
  beta = sum(sapply(x_vector,inner_int, eta = eta, grid_width = grid_width,
                    p = hat_q, p_q_ratio = hat_p_q_ratio) / h * (x_vector[2] - x_vector[1]))
  return(beta)
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













