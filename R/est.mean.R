#' Estimating Sample Mean using Quantiles
#'
#' @description
#' This function estimates the sample mean from a study presenting quantile summary measures with the sample size (\eqn{n}). The quantile summaries can fall into one of the following categories:
#' \itemize{
#'   \item \eqn{S_1}: \{ minimum, median, maximum \}
#'   \item \eqn{S_2}: \{ first quartile, median, third quartile \}
#'   \item \eqn{S_3}: \{ minimum, first quartile, median, third quartile, maximum \}
#' }
#'
#' The \code{est.mean} function implements newly proposed flexible quantile-based distribution methods for estimating sample mean (De Livera et al., 2024). 
#' It also incorporates existing methods for estimating sample means as described by Luo et al. (2018) and McGrath et al. (2020).
#' 
#' 
#' @usage est.mean(
#'    min = NULL, 
#'    q1 = NULL, 
#'    med = NULL, 
#'    q3 = NULL, 
#'    max = NULL, 
#'    n = NULL, 
#'    method = "gld/sld", 
#'    opt = TRUE
#'    )
#' 
#' @param min numeric value representing the sample minimum.
#' @param q1 numeric value representing the first quartile of the sample.
#' @param med numeric value representing the median of the sample.
#' @param q3 numeric value representing the third quartile of the sample.
#' @param max numeric value representing the sample maximum.
#' @param n numeric value specifying the sample size.
#' @param method character string specifying the approach used to estimate the sample means. The options are the following:
#' \describe{
#'   \item{\code{'gld/sld'}}{The default option. The method proposed by De Livera et al. (2024). Estimation using the Generalized Lambda Distribution (GLD) for 5-number summaries (\eqn{S_3}), and the Skew Logistic Distribution (SLD)  for 3-number summaries (\eqn{S_1} and \eqn{S_2}).}
#'   \item{\code{'luo'}}{Method of Luo et al. (2018).}
#'   \item{\code{'hozo/wan/bland'}}{The method proposed by Wan et al. (2014). i.e., the method of Hozo et al. (2005) for \eqn{S_1}, method of Wan et al. (2014) for \eqn{S_2}, and method of Bland (2015) for \eqn{S_3}.} 
#'   \item{\code{'bc'}}{Box-Cox method proposed by McGrath et al. (2020).}
#'   \item{\code{'qe'}}{Quantile Matching Estimation method proposed by McGrath et al. (2020).}
#' }
#' @param opt logical value indicating whether to apply the optimization step of \code{'gld/sld'} method, in estimating their parameters using theoretical quantiles. 
#'   The default value is \code{TRUE}.
#' 
#' @details
#' The \code{'gld/sld'} method (i.e., the method of De Livera et al., (2024)) of \code{est.mean} uses the following quantile based distributions:
#' \itemize{
#'   \item Generalized Lambda Distribution (GLD) for estimating the sample mean using 5-number summaries (\eqn{S_3}).
#'   \item Skew Logistic Distribution (SLD) for estimating the sample mean using 3-number summaries (\eqn{S_1} and \eqn{S_2}).
#' } 
#' The generalised lambda distribution (GLD) is a four parameter family of distributions defined by its quantile function under the FKML parameterisation (Freimer et al., 1988).
#' De Livera et al. propose that the GLD quantlie function can be used to approximate a sample's distribution using 5-point summaries. 
#' The four parameters of GLD quantile function include: a location parameter (\eqn{\lambda_1}), an inverse scale parameter (\eqn{\lambda_2}>0), and two shape parameters (\eqn{\lambda_3} and \eqn{\lambda_4}).
#' 
#' The quantile-based skew logistic distribution (SLD), introduced by Gilchrist (2000) and further modified by van Staden and King (2015) 
#' is used to approximate the sample's distribution using 3-point summaries.
#' The SLD quantile function is defined using three parameters: a location parameter (\eqn{\lambda}), a scale parameter (\eqn{\eta}), and a skewing parameter (\eqn{\delta}).
#' 
#' For \code{'gld/sld'} method, the parameters of the generalized lambda distribution (GLD) and skew logistic distribution (SLD) are estimated 
#' by formulating and solving a set of simultaneous equations. These equations relate the estimated sample quantiles to their theoretical counterparts
#' of the respective distribution (GLD or SLD). Finally, the mean for each scenario is calculated by integrating functions of the estimated quantile function.
#'
#' @return \code{mean}: numeric value representing the estimated mean of the sample.
#'   
#' @examples
#' #Generate 5-point summary data
#' set.seed(123)
#' n <- 1000
#' x <- stats::rlnorm(n, 4, 0.3)
#' quants <- c(min(x), stats::quantile(x, probs = c(0.25, 0.5, 0.75)), max(x))
#' obs_mean <- mean(x)
#' 
#' #Estimate sample mean using s3 (5 number summary)
#' est_mean_s3 <- est.mean(min = quants[1], q1 = quants[2], med = quants[3], q3 = quants[4], 
#'                         max = quants[5], n=n, method = "gld/sld")
#' est_mean_s3
#' 
#' #Estimate sample mean using s1 (min, median, max)
#' est_mean_s1 <- est.mean(min = quants[1], med = quants[3], max = quants[5],
#'                         n=n, method = "gld/sld")
#' est_mean_s1
#'
#' #Estimate sample mean using s2 (q1, median, q3)
#' est_mean_s2 <- est.mean(q1 = quants[2], med = quants[3], q3 = quants[4],
#'                         n=n, method = "gld/sld")
#' est_mean_s2
#'
#' @references Alysha De Livera, Luke Prendergast, and Udara Kumaranathunga. A novel density-based approach for estimating unknown means, distribution visualisations, and meta-analyses of quantiles. \emph{Submitted for Review}, 2024, pre-print available here: <https://arxiv.org/abs/2411.10971>
#' @references Dehui Luo, Xiang Wan, Jiming Liu, and Tiejun Tong. Optimally estimating the sample mean from the sample size, median, mid-range, and/or mid-quartile range. \emph{Statistical methods in medical research}, 27(6):1785–1805,2018.
#' @references Xiang Wan, Wenqian Wang, Jiming Liu, and Tiejun Tong. Estimating the sample mean and standard deviation from the sample size, median, range and/or interquartile range. \emph{BMC medical research methodology}, 14:1–13, 2014.
#' @references Sean McGrath, XiaoFei Zhao, Russell Steele, Brett D Thombs, Andrea Benedetti, and DEPRESsion Screening Data (DEPRESSD) Collaboration. Estimating the sample mean and standard deviation from commonly reported quantiles in meta-analysis. \emph{Statistical methods in medical research}, 29(9):2520–2537, 2020b.
#' @references Marshall Freimer, Georgia Kollia, Govind S Mudholkar, and C Thomas Lin. A study of the generalized tukey lambda family. \emph{Communications in Statistics-Theory and Methods}, 17(10):3547–3567, 1988.
#' @references Warren Gilchrist. \emph{Statistical modelling with quantile functions}. Chapman and Hall/CRC, 2000.
#' @references P. J. van Staden and R. A. R. King. The quantile-based skew logistic distribution.  \emph{Statistics & Probability Letters}, 96:109–116, 2015.
#' @export 
#' 
#' @importFrom gld qgl gld.moments
#' @importFrom sld qsl
#' @importFrom stats integrate optim integrate runif var qnorm
#' @importFrom estmeansd bc.mean.sd qe.mean.sd


est.mean <- function(min = NULL, 
                     q1 = NULL, 
                     med = NULL, 
                     q3 = NULL, 
                     max = NULL, 
                     n = NULL, 
                     method = "gld/sld", 
                     opt = TRUE) {
  
  #5-number summary
  if (!is.null(min) && !is.null(q1) && !is.null(med) && !is.null(q3) && !is.null(max)) {
    
    if (method == "gld/sld") {
      glsl_est <- est.density.five(min, q1, med, q3, max, n, opt)
      mean_est <- glsl_est$mean
      return(list("mean" = mean_est))
      
    } else if (method == "hozo/wan/bland") {
      wan_est <- wan_using5(min, q1, med, q3, max, n)
      mean_est <- as.numeric(wan_est$mean_est)
      return(list("mean" = mean_est))
      
    }else if (method == "luo") {
      mean_est <- luo_using5(min, q1, med, q3, max, n)
      return(list("mean" = mean_est))
      
    } else if (method == "bc") {
      if (any(c(min, q1, med, q3, max) <= 0)) {
        add <- abs(min(c(min, q1, med, q3, max))) + 0.5
      } else {
        add <- 0
      }
      q_bc <- c(min, q1, med, q3, max) + add
      mean_est_bc <- bc.mean.sd(min.val = q_bc[1], q1.val = q_bc[2], med.val = q_bc[3], q3.val = q_bc[4], max.val = q_bc[5], n = n)
      mean_est <- mean_est_bc$est.mean - add
      return(list("mean" = mean_est))
      
    } else if (method == "qe") {
      if (any(c(min, q1, med, q3, max) <= 0)) {
        add <- abs(min(c(min, q1, med, q3, max))) + 0.5
      } else {
        add <- 0
      }
      q_qe <- c(min, q1, med, q3, max) + add
      mean_est_qe <- qe.mean.sd(min.val = q_qe[1], q1.val = q_qe[2], med.val = q_qe[3], q3.val = q_qe[4], max.val = q_qe[5], n = n)
      mean_est <- mean_est_qe$est.mean - add
      return(list("mean" = mean_est))
      
    } else {
      stop("Unsupported method.")
    }
  
  #3-number summary: min,med,max  
  } else if (!is.null(min) && !is.null(med) && !is.null(max)) {
    
    if (method == "gld/sld"){
      glsl_est <- est.density.three1(min=min, med=med, max=max, n=n, opt=opt)
      mean_est <- glsl_est$mean
      return(list("mean" = mean_est))
      
    } else if (method == "hozo/wan/bland"){
      wan_est <- wan_using_minq2max(min, med, max, n)
      mean_est <- as.numeric(wan_est$mean_est)
      return(list("mean" = mean_est))
      
    }else if (method == "luo"){
      mean_est <- luo_using_minq2max(min, med, max, n)
      return(list("mean" = mean_est))
      
    } else if (method == "bc"){
      if (any(c(min, med, max) <= 0)) {
        add <- abs(min(c(min, med, max))) + 0.5
      } else {
        add <- 0
      }
      q_bc <- c(min,med,max) + add
      mean_est_bc <- bc.mean.sd(min.val = q_bc[1], med.val = q_bc[2], max.val = q_bc[3], n = n)
      mean_est <- mean_est_bc$est.mean - add
      return(list("mean" = mean_est))
      
    } else if (method == "qe"){
      if (any(c(min, med, max) <= 0)) {
        add <- abs(min(c(min, med, max))) + 0.5
      } else {
        add <- 0
      }
      q_qe <- c(min,med,max) + add
      mean_est_qe <- qe.mean.sd(min.val = q_qe[1], med.val = q_qe[2], max.val = q_qe[3], n = n)
      mean_est <- mean_est_qe$est.mean - add
      return(list("mean" = mean_est))
      
    } else {
      stop("Unsupported method.")
    }
  
  #3-number summary: q1,med,q3  
  } else if (!is.null(q1) && !is.null(med) && !is.null(q3)) {
    
    if (method == "gld/sld") {
      glsl_est <- est.density.three2(q1=q1, med=med, q3=q3, opt=opt)
      mean_est <- glsl_est$mean
      return(list("mean" = mean_est))
      
    } else if (method == "hozo/wan/bland") {
      wan_est <- wan_using_q1q2q3(q1, med, q3, n)
      mean_est <- as.numeric(wan_est$mean_est)
      return(list("mean" = mean_est))
      
    } else if (method == "luo") {
      mean_est <- luo_using_q1q2q3(q1, med, q3, n)
      return(list("mean" = mean_est))
      
    } else if (method == "bc") {
      if (any(c(q1, med, q3) <= 0)) {
        add <- abs(min(c(q1, med, q3))) + 0.5
      } else {
        add <- 0
      }
      q_bc <- c(q1, med, q3) + add
      mean_est_bc <- bc.mean.sd(q1.val = q_bc[1], med.val = q_bc[2], q3.val = q_bc[3], n = n)
      mean_est <- mean_est_bc$est.mean - add
      return(list("mean" = mean_est))
      
    } else if (method == "qe") {
      if (any(c(q1, med, q3) <= 0)) {
        add <- abs(min(c(q1, med, q3))) + 0.5
      } else {
        add <- 0
      }
      q_qe <- c(q1, med, q3) + add
      mean_est_qe <- qe.mean.sd(q1.val = q_qe[1], med.val = q_qe[2], q3.val = q_qe[3], n = n)
      mean_est <- mean_est_qe$est.mean - add
      return(list("mean" = mean_est))
      
    } else {
      stop("Unsupported method.")
    }
    
  } else {
    stop("Invalid input: Please provide either {min, q1, med, q3, max} or {min, med, max} or {q1, med, q3}.\n")
  }
}




est.density.five <- function(min = NULL, 
                             q1 = NULL, 
                             med = NULL, 
                             q3 = NULL, 
                             max = NULL, 
                             n = NULL, 
                             opt = TRUE) {
    
  u <- 0.5 / n
  rho1.hat <- med
  rho2.hat <- max - min
    
  rho3.hat <- (med - min) / (max - med)
  rho4.hat <- (q3 - q1) / rho2.hat
    
  res <- optim(par = c(0.5, 0.5),
               Fmin,
               rho3.hat = rho3.hat,
               rho4.hat = rho4.hat,
               u = u)
    
   l3 <- res$par[1]
   l4 <- res$par[2]
   l2 <- (Sfkml_l3l4(1 - u, l3, l4) - Sfkml_l3l4(u, l3, l4)) / rho2.hat
   l1 <- rho1.hat - Sfkml_l3l4(1 / 2, l3, l4) / l2
    
   params_temp <- c(l1, l2, l3, l4)
    
   if (opt) {
      params <- optim(par = params_temp, n = n, empirical_quantiles = c(min, q1, med, q3, max), fn = objective_function_gl)$par
   } else {
      params <- params_temp
   }
    
   names(params) <- c("location", "inverse scale", "shape 1", "shape 2")
   res <- true_summary("gl", list(params))
   mean_est <- as.numeric(res$mean)
   sd_est <- as.numeric(sqrt(res$variance))
    
   return(list("parameters" = params, "mean"= mean_est, "sd"= sd_est))
}


est.density.three1 <- function(min = NULL, 
                               q1 = NULL, 
                               med = NULL, 
                               q3 = NULL, 
                               max = NULL, 
                               n = NULL, 
                               opt = TRUE) {
  a <- min
  m <- med
  b <- max
    
  v <- (b - m) / (m - a)
  d <- log((2 * n - 1) / n)
  delta <- max(((d - (log(n) * v)) / ((d - log(n)) * (v + 1))), 0)
  delta <- min(delta, 0.99)
  eta <- (b - a) / log(2 * n - 1)
  lambda <- m + eta * log(2) * (1 - 2 * delta)
    
  params_temp <- c(lambda, eta, delta)
    
  if (opt) {
    params <- optim(par = params_temp, n = n, empirical_quantiles = c(min, med, max), fn = objective_function_sl_minmax)$par
  } else {
    params <- params_temp
  }
    
  names(params) <- c("location", "scale", "skewing")
  res <- true_summary("sl", list(params))
  mean_est <- as.numeric(res$mean)
  sd_est <- as.numeric(sqrt(res$variance))
    
  return(list("parameters" = params, "mean"= mean_est, "sd"= sd_est))
    
}



est.density.three2 <- function(min = NULL, 
                               q1 = NULL, 
                               med = NULL, 
                               q3 = NULL, 
                               max = NULL, 
                               n = NULL, 
                               opt = TRUE) {
  
  rho <- (q3 - med) / (med - q1)
  delta <- max((log(3 / 2) - log(2) * rho) / log(3 / 4) / (rho + 1), 0)
  delta <- min(delta, 0.99)
  eta <- (q3 - q1) / log(3)
  lambda <- med + eta * log(2) * (1 - 2 * delta)
    
  params_temp <- c(lambda, eta, delta)
    
  if (opt) {
    params <- optim(par = params_temp, empirical_quantiles = c(q1, med, q3), fn = objective_function_sl)$par
  } else {
      params <- params_temp
  }
    
  names(params) <- c("location", "scale", "skewing")
  res <- true_summary("sl", list(params))
  mean_est <- as.numeric(res$mean)
  sd_est <- as.numeric(sqrt(res$variance))
    
  return(list("parameters" = params, "mean"= mean_est, "sd"= sd_est))
    
}

