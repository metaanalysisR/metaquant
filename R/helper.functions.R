################################################################################

Sfkml_l3l4 <- function(p,
                       l3,
                       l4 
){ 
  (p^l3 - 1)/l3 - ((1 - p)^l4 - 1)/l4
}

################################################################################

Fmin<- function(opt_pars, 
                rho3.hat,
                rho4.hat, 
                u){
  l3<-opt_pars[1]
  l4<-opt_pars[2]
  
  if (l3==0)
    return (10^10)
  if (l4==0)
    return (10^10)
  
  rho2 <- Sfkml_l3l4(1 - u, l3, l4) - 
    Sfkml_l3l4(u, l3, l4)  
  rho3 <- (Sfkml_l3l4(1/2, l3, l4) - 
             Sfkml_l3l4(u, l3, l4))/
    (Sfkml_l3l4(1 - u, l3, l4) - 
       Sfkml_l3l4(1/2, l3, l4))
  rho4 <- (Sfkml_l3l4(3/4, l3, l4) - Sfkml_l3l4(1/4, l3, l4))/rho2 
  
  if ((Sfkml_l3l4(1 - u, l3, l4) - Sfkml_l3l4(u, l3, l4))<0) 
    return(10^10)
  else
    return((rho3.hat - rho3)^2 + (rho4.hat - rho4)^2)
}


################################################################################

objective_function_gl <- function(par,n,empirical_quantiles) {
  lambda1 <- par[1]
  lambda2 <- par[2]
  lambda3 <- par[3]
  lambda4 <- par[4]
  u<-0.5/n
  
  if (lambda2<=0)
    return(1e10)
  
  theoretical_quantiles <- qgl(p = c(u, 0.25, 0.5, 0.75, 1-u),
                               lambda1 = c(lambda1, lambda2, lambda3, lambda4),
                               param = "fkml")
  
  return(sum((theoretical_quantiles - empirical_quantiles)^2))
  
  
}


################################################################################

objective_function_sl <- function(par,empirical_quantiles) { #n,
  
  alpha<- par[1] #lambda 
  beta <- par[2] #eta
  delta <- par[3] #delta
  
  if (delta<0 | delta >1 | beta <= 0)
    return(1e10)
  
  theoretical_quantiles <- qsl(p = c(0.25, 0.5, 0.75),
                               parameters = c(alpha,beta, delta)
  )
  
  return(sum((theoretical_quantiles - empirical_quantiles)^2))
  
  
}


################################################################################

objective_function_sl_minmax <- function(par,empirical_quantiles,n) { #n,
  
  alpha<- par[1] #lambda 
  beta <- par[2] #eta
  delta <- par[3] #delta
  u<-0.5/n
  
  if (delta<0 | delta >1 | beta <= 0)
    return(1e10)
  
  
  theoretical_quantiles <- qsl(p = c(u, 0.5,1-u),
                               parameters = c(alpha,beta, delta)
  )
  
  return(sum((theoretical_quantiles - empirical_quantiles)^2))
  
  
}


################################################################################

true_summary<- function(dist, param){ 
  qf <- paste0("q", dist)
  qfp <- Vectorize(function(p, dist, param) 
    do.call(qf, c(p = p, param)), "p")
  qfp.sq <- Vectorize(function(p, dist, param) 
    do.call(qf, c(p = p, param))^2, "p")
  
  mu<-NA
  if(qf=="qgl"){
    
    gld_moments<-gld.moments(unlist(param))
    mu<-gld_moments[1]
    sigma2<-gld_moments[2]
    
  }

  if (is.na(mu))
    mu <- tryCatch(integrate(qfp, 0, 1, 
                             dist = dist, 
                             param = param)$value, error = function(e) NA)
  if (is.na(mu))
    mu<- tryCatch(
      integrate(qfp, 0 , 0.5, dist = dist, param = param)$value + 
        integrate(qfp, 0.5, 0.99999999, dist = dist, param = param)$value, 
      error = function(e) NA)
  
  if (is.na(mu))
    mu<- tryCatch(
      integrate(qfp, 0 , 0.5, dist = dist, param = param)$value + 
        integrate(qfp, 0.5, 0.9999999, dist = dist, param = param)$value, 
      error = function(e) NA)
  
  if (is.na(mu))
    mu<- tryCatch(
      integrate(qfp, 0 , 0.5, dist = dist, param = param)$value + 
        integrate(qfp, 0.5, 0.999999, dist = dist, param = param)$value, 
      error = function(e) NA)
  
  if (is.na(mu)){
    n <- 1000000
    p_values <- runif(n,min=0,max=1)
    if(qf=="qgl")
      samples <-qgl(p_values,unlist(param)) 
    else if(qf=="qsl")
      samples <-qsl(p_values,unlist(param)) 
    mu<-mean(samples)
  }
  
  sigma2 <- tryCatch(integrate(qfp.sq, 0, 1,
                               dist = dist, 
                               param = param)$value, 
                     error = function(e) NA) - mu^2
  if (is.na(sigma2)){
    sigma2 <- tryCatch(
      integrate(qfp.sq, 0, 0.5 , dist = dist, param = param)$value+
        integrate(qfp.sq, 0.5, 0.999999 , dist = dist, param = param)$value, 
      error = function(e) NA) - mu^2
  }
  
  if (sigma2<=0 | is.na(sigma2)){
    n <- 1000000
    p_values <- runif(n)
    if(qf=="qgl")
      samples <-  qgl(p_values,unlist(param)) 
    if(qf=="qsl")
      samples <-  qsl(p_values,unlist(param))   
    sigma2<-var(samples)
  }

  res <- list(mean = mu, variance = sigma2)
  class(res) <- c("table")
  res
}


################################################################################

luo_using5<-function(min, q1, med, q3, max, n_size){
  a<-min
  q1<-q1
  m<-med
  q3<-q3
  b<-max
  
  weight_est1<-2.2/(2.2+n_size^0.75)
  weight_est2<-0.7-0.72/n_size^0.55
  
  mean_est<-weight_est1*(a+b)/2+weight_est2*(q1+q3)/2+(1-weight_est1-weight_est2)*m

  return(as.numeric(mean_est))
}

luo_using_minq2max<-function(min, med, max, n_size){ 
  
  weight=4/(4+n_size^0.75)
  
  mean_est=weight*(min+max)/2+(1-weight)*med
  return(as.numeric(mean_est))
}

luo_using_q1q2q3<-function(q1, med, q3, n_size){ 
  
  weight=0.7+0.39/n_size
  
  mean_est<-weight*(q1+q3)/2+(1-weight)*med
  
  return(as.numeric(mean_est))
}



################################################################################

shi_using5<-function(min, q1, med, q3, max, n_size){

  weight=1/(1+0.07*n_size^0.6)
  sd_est=weight*(max-min)/(2*qnorm((n_size-0.375)/(n_size+0.25),0,1))+
    (1-weight)*(q3-q1)/(2*qnorm((0.75*n_size-0.125)/(n_size+0.25),0,1))
  
  return(as.numeric(sd_est))
}

################################################################################

wan_using5<-function(min, q1, med, q3, max, n_size){ 
  
  mean_est<-(min+2*q1+2*med+2*q3+max)/8 
  sd_est<-((max-min)/(4*qnorm((n_size-0.375)/(n_size+0.25))))+
    ((q3-q1)/(4*qnorm((0.75*n_size-0.125)/(n_size+0.25))))
  
  return(list(mean_est=mean_est, sd_est=sd_est))
}

wan_using_minq2max<-function(min, med, max, n_size){ 
  
  weight=4/(4+n_size^0.75)
  
  mean_est=(min+2*med+max)/4
  sd_est=(max-min)/(2*qnorm((n_size-0.375)/(n_size+0.25)))
  
  return(list(mean_est=mean_est,sd_est=sd_est))
}

wan_using_q1q2q3<-function(q1, med, q3, n_size){ 
  
  mean_est=(q1+med+q3)/3
  sd_est=(q3-q1)/(2*qnorm((0.75*n_size-0.125)/(n_size+0.25)))
  
  return(list(mean_est=mean_est, sd_est=sd_est))
}