########################### GLD and SLD functions ##############################

Sfkml_l3l4 <- function(p,
                       l3,
                       l4 
){ 
  (p^l3 - 1)/l3 - ((1 - p)^l4 - 1)/l4
}


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


######################## Other methods (mean estimation) ########################

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



####################### Other methods (SD estimation) ########################

shi_using5<-function(min, q1, med, q3, max, n_size){

  weight=1/(1+0.07*n_size^0.6)
  sd_est=weight*(max-min)/(2*qnorm((n_size-0.375)/(n_size+0.25),0,1))+
    (1-weight)*(q3-q1)/(2*qnorm((0.75*n_size-0.125)/(n_size+0.25),0,1))
  
  return(as.numeric(sd_est))
}

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


########################## metaquant() helper functions ########################

get.fp.gld <- function(x, p) {
  
  params <- x$parameters
  p_est <- gld::qgl(p, params)
  fp <- gld::dgl(gld::qgl(p, params), params)
  
  return(list("est_perc" = p_est, "est_f"= fp))
}


get.fp.qe <- function(x, family, p) {
  if (!(family %in% c("normal", "log-normal", "gamma", "weibull"))) {
    stop("family must be either normal, log-normal, gamma, or Weibull.")
  }
  if (family == "normal") {
    par <- x$norm.par
    p_est <- stats::qnorm(p, par[1], par[2])
    fp <- stats::dnorm(stats::qnorm(p, par[1], par[2]), par[1], par[2])
  }
  else if (family == "log-normal") {
    par <- x$lnorm.par
    p_est <- stats::qlnorm(p, par[1], par[2])
    fp <- stats::dlnorm(stats::qlnorm(p, par[1], par[2]), par[1], par[2])
  }
  else if (family == "gamma") {
    par <- x$gamma.par
    p_est <- stats::qgamma(p, par[1], par[2])
    fp <- stats::dgamma(stats::qgamma(p, par[1], par[2]), par[1], par[2])
  }
  else if (family == "weibull") {
    par <- x$weibull.par
    p_est <- stats::qweibull(p, par[1], par[2])
    fp <- stats::dweibull(stats::qweibull(p, par[1], par[2]), par[1], par[2])
  }
  
  return(list("est_perc" = p_est, "est_f"= fp))
}


check.meta.df <- function(df){
  if (!is.data.frame(df)){
    stop("Argument 'data' must be a data frame.")
  } else {
    if ("tbl_df" %in% class(df)) {
      df <- as.data.frame(df)
    }
  }
  
  possible.cols <- c('min.g1', 'q1.g1', 'med.g1', 'q3.g1', 'max.g1', 'n.g1',
                     'min.g2', 'q1.g2', 'med.g2', 'q3.g2', 'max.g2', 'n.g2')
  
  missing.cols <- possible.cols[which(!possible.cols %in% colnames(df))]
  df[, missing.cols] <- NA
  
  for (col in possible.cols){
    if (!class(df[, col]) %in% c('numeric', 'integer')){
      if (class(df[, col]) %in% 'logical' | all(is.na(df[, col]))){
        next()
      }
      stop(sprintf(
        "Summary columns must be numeric or integer; column '%s' is not one of these types.",
        col
      ), call. = FALSE)
    }
  }
  
  temp <- df[, possible.cols, drop = FALSE]
  all.na.row <- rowSums(!is.na(temp)) == 0
  
  if (any(all.na.row)) {
    if (all(all.na.row)) {
      stop(
        "All rows in the input dataset have all summary data set to NA. 
        This often indicates incorrect or missing column names.",
        call. = FALSE
      )
    } else {
      stop(
        sprintf(
          "The following rows have all summary data set to NA: %s",
          paste(which(all.na.row), collapse = " ")
        ),
        call. = FALSE
      )
    }
  }
  
  return(df)
}



############################# plotdist() heper functions #######################

check.df <- function(df) {
  #Assign 'study.index' if missing
  if (!"study.index" %in% colnames(df)) {
    warning("The 'study.index' column is missing. It has been added automatically as Study1, Study2, etc..")
    df$study.index <- paste0("study", seq_len(nrow(df)))
  }
  #Check for Group 1 columns and Group 2 columns
  s1_g1 <- c("study.index", "min.g1", "med.g1", "max.g1", "n.g1")
  s2_g1 <- c("study.index", "q1.g1", "med.g1", "q3.g1")
  s3_g1 <- c("study.index", "min.g1", "q1.g1", "med.g1", "q3.g1", "max.g1", "n.g1")
  
  valid_group1 <- all(s1_g1 %in% colnames(df)) ||
    all(s2_g1 %in% colnames(df)) ||
    all(s3_g1 %in% colnames(df))
  
  if (!valid_group1) {
    stop("The group 1 does not have the required quantile columns for estimation.")
  }
  g2_present <- any(grepl("\\.g2$", colnames(df)))
  if (g2_present) {
    s1_g2 <- c("study.index", "min.g2", "med.g2", "max.g2", "n.g2")
    s2_g2 <- c("study.index", "q1.g2", "med.g2", "q3.g2")
    s3_g2 <- c("study.index", "min.g2", "q1.g2", "med.g2", "q3.g2", "max.g2", "n.g2")
    
    valid_group2 <- all(s1_g2 %in% colnames(df)) ||
      all(s2_g2 %in% colnames(df)) ||
      all(s3_g2 %in% colnames(df))
    
    if (!valid_group2) {
      stop("The group 2 does not have the required quantile columns for estimation.")
    }
  }
  
  #Check for non-numeric values in quantile columns
  invalid_rows_g1 <-NULL
  invalid_rows_g2 <-NULL
  quantile_cols_g1 <-grep("\\.g1$", colnames(df), value = TRUE)
  quantile_cols_g2 <-grep("\\.g2$", colnames(df), value = TRUE)
  
  if (length(quantile_cols_g1)>0) {
    invalid_rows_g1 <- which(rowSums(sapply(df[quantile_cols_g1], function(col) {
      !is.finite(as.numeric(col))
    })) > 0)
  }
  if (length(quantile_cols_g2)>0) {
    invalid_rows_g2 <- which(rowSums(sapply(df[quantile_cols_g2], function(col) {
      !is.finite(as.numeric(col))
    })) > 0)
  }
  if (length(invalid_rows_g1)>0) {
    warning(paste("Rows", paste(invalid_rows_g1, collapse = ", "), 
                  "contain invalid data for Group 1 and will be skipped from density visualisation."))
  }
  if (length(invalid_rows_g2)>0) {
    warning(paste("Rows", paste(invalid_rows_g2, collapse = ", "), 
                  "contain invalid data for Group 2 and will be skipped from density visualisation."))
  }
  
  return(list(df = df, invalid_rows_g1 = invalid_rows_g1, invalid_rows_g2 = invalid_rows_g2))
}


get.scenario <- function(data, gi) {
  if (all(c(paste0("min.", gi), paste0("q1.", gi), paste0("med.", gi), 
            paste0("q3.", gi), paste0("max.", gi)) %in% colnames(data))) {
    return("5-number")
  } else if (all(c(paste0("min.", gi), paste0("med.", gi), paste0("max.", gi)) %in% colnames(data))) {
    return("min-med-max")
  } else if (all(c(paste0("q1.", gi), paste0("med.", gi), paste0("q3.", gi)) %in% colnames(data))) {
    return("q1-med-q3")
  } else {
    stop(paste("Invalid or incomplete data for group", gi))
  }
}


get.density <- function(data, gi, summary_type, group_label, xmin, xmax, length.out, opt) {
  group_data <- data.frame()
  
  if (is.null(xmin) || is.null(xmax)) {
    if (summary_type == "5-number" || summary_type == "min-med-max") {
      if (is.null(xmin)) xmin <-min(data[[paste0("min.", gi)]], na.rm=TRUE)
      if (is.null(xmax)) xmax <-max(data[[paste0("max.", gi)]], na.rm=TRUE)
    } else if (summary_type == "q1-med-q3") {
      if (is.null(xmin) || is.null(xmax)) {
        stop("Please set argument 'xmin' and 'xmax' to compute density data (as default values cannot be computed for the given quantile scenario)")
      }
    }
  }
  
  for (i in seq_len(nrow(data))) {
    n <- data[[paste0("n.", gi)]][i]
    
    if (summary_type == "5-number") {
      q <- as.numeric(data[i, c(paste0("min.", gi), paste0("q1.", gi), paste0("med.", gi), 
                                paste0("q3.", gi), paste0("max.", gi))])
      dens_out <- est.gld.five(min = q[1], q1 = q[2], med = q[3], q3 = q[4], max = q[5], n = n, opt = opt)
      x_vals <- seq(xmin, xmax, length.out = length.out)
      y_vals <- dgl(x_vals, lambda1 = dens_out$parameters[1], lambda2 = dens_out$parameters[2], 
                    lambda3 = dens_out$parameters[3], lambda4 = dens_out$parameters[4])
    } else if (summary_type == "min-med-max") {
      q <- as.numeric(data[i, c(paste0("min.", gi), paste0("med.", gi), paste0("max.", gi))])
      dens_out <- est.sld.minq2max(min = q[1], med = q[2], max = q[3], n = n, opt = opt)
      x_vals <- seq(xmin, xmax, length.out = length.out)
      y_vals <- dsl(x_vals, dens_out$parameters)
    } else if (summary_type == "q1-med-q3") {
      q <- as.numeric(data[i, c(paste0("q1.", gi), paste0("med.", gi), paste0("q3.", gi))])
      dens_out <- est.sld.q1q2q3(q1 = q[1], med = q[2], q3 = q[3], n = n, opt = opt)
      x_vals <- seq(xmin, xmax, length.out = length.out)
      y_vals <- dsl(x_vals, dens_out$parameters)
    }
    
    group_density <- data.frame(x=x_vals, y=y_vals, Study=data$study.index[i], Group=group_label)
    group_data <- rbind(group_data, group_density)
  }
  
  return(group_data)
}


get.pooled.density <- function(data, gi, summary_type, group_label, xmin, xmax, length.out, opt) {
  
  if (is.null(xmin) || is.null(xmax)) {
    if (summary_type == "5-number" || summary_type == "min-med-max") {
      if (is.null(xmin)) xmin <-min(data[[paste0("min.", gi)]], na.rm=TRUE)
      if (is.null(xmax)) xmax <-max(data[[paste0("max.", gi)]], na.rm=TRUE)
    } else if (summary_type == "q1-med-q3") {
      if (is.null(xmin) || is.null(xmax)) {
        stop("Please set argument 'xmin' and 'xmax' to compute density data (as default values cannot be computed for the given quantile scenario)")
      }
    }
  }
  
  if (summary_type == "5-number") {
    data$l1 <-NA
    data$l2 <-NA
    data$l3 <-NA
    data$l4 <-NA
  } else {
    data$lambda <-NA
    data$eta <-NA
    data$delta <-NA
  }
  
  #Store parameters for each row
  for (i in seq_len(nrow(data))) {
    n <- data[[paste0("n.", gi)]][i]
    if (summary_type == "5-number") {
      q <- as.numeric(data[i, c(paste0("min.", gi), paste0("q1.", gi), paste0("med.", gi), 
                                paste0("q3.", gi), paste0("max.", gi))])
      dens_out <- est.gld.five(min = q[1], q1 = q[2], med = q[3], q3 = q[4], max = q[5], n = n, opt = opt)
      data$l1[i] <- dens_out$parameters[1]
      data$l2[i] <- dens_out$parameters[2]
      data$l3[i] <- dens_out$parameters[3]
      data$l4[i] <- dens_out$parameters[4]
    } else if(summary_type == "min-med-max"){
      q <- as.numeric(data[i, c(paste0("min.", gi), paste0("med.", gi), paste0("max.", gi))])
      dens_out <- est.sld.minq2max(min = q[1], med = q[2], max = q[3], n = n, opt = opt)
      data$lambda[i] <- dens_out$parameters[1]
      data$eta[i] <- dens_out$parameters[2]
      data$delta[i] <- dens_out$parameters[3]
    } else if(summary_type == "q1-med-q3"){
      q <- as.numeric(data[i, c(paste0("q1.", gi), paste0("med.", gi), paste0("q3.", gi))])
      dens_out <- est.sld.q1q2q3(q1 = q[1], med = q[2], q3 = q[3], n = n, opt = opt)
      data$lambda[i] <- dens_out$parameters[1]
      data$eta[i] <- dens_out$parameters[2]
      data$delta[i] <- dens_out$parameters[3]
    }
  }
  
  # weights
  data$weight <- data[[paste0("n.", gi)]]/sum(data[[paste0("n.", gi)]])
  
  # pooled parameters
  if (summary_type == "5-number"){
    pooled.l1 <-sum(data$weight*data$l1,na.rm =TRUE)
    pooled.l2 <- sum(data$weight*data$l2,na.rm =TRUE)
    pooled.l3 <-sum(data$weight*data$l3,na.rm =TRUE)
    pooled.l4 <- sum(data$weight*data$l4,na.rm =TRUE)
    x_vals <-seq(xmin, xmax, length.out = length.out)
    y_vals <- dgl(x_vals, lambda1 =pooled.l1, lambda2 =pooled.l2, lambda3 =pooled.l3, lambda4 =pooled.l4)
  } else{
    pooled.lambda <-sum(data$weight*data$lambda, na.rm =TRUE)
    pooled.eta <-sum(data$weight*data$eta, na.rm =TRUE)
    pooled.delta <-sum(data$weight*data$delta, na.rm =TRUE)
    x_vals <-seq(xmin, xmax, length.out =length.out)
    y_vals <-dsl(x_vals, c(pooled.lambda, pooled.eta, pooled.delta))
  }
  
  return(data.frame(x=x_vals, y=y_vals, Group = paste0(group_label, " (Pooled)"), Study="Pooled"))
}
