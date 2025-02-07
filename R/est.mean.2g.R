#' Estimating Sample Means of Two Groups using Quantiles 
#'
#' @description
#' This function estimates the sample means from a two group study presenting quantile summary measures with the sample size (\eqn{n}). 
#' The quantile summaries of each group can fall into one of the following categories:
#' \itemize{
#'   \item \eqn{S_1}: \{ minimum, median, maximum \} 
#'   \item \eqn{S_2}: \{ first quartile, median, third quartile \}
#'   \item \eqn{S_3}: \{ minimum, first quartile, median, third quartile, maximum \}
#' }
#'
#' The \code{est.mean.2g} function uses a novel quantile-based distribution methods for estimating sample mean for two groups such as 'Treatment' and 'Control' (De Livera et al., 2024). 
#' The method is based on the following quantile-based distributions for 
#' estimating sample means:
#' \itemize{
#'   \item Generalized Lambda Distribution (GLD) for estimating sample means using 
#'         5-number summaries (\eqn{S_3}).
#'   \item Skew Logistic Distribution (SLD) for estimating sample means using 
#'         3-number summaries (\eqn{S_1} and \eqn{S_2}).
#' }
#' 
#' @usage
#' est.mean.2g(
#'    min.g1 = NULL, 
#'    q1.g1 = NULL, 
#'    med.g1 = NULL, 
#'    q3.g1 = NULL, 
#'    max.g1 = NULL,
#'    min.g2 = NULL, 
#'    q1.g2 = NULL, 
#'    med.g2 = NULL, 
#'    q3.g2 = NULL, 
#'    max.g2 = NULL,
#'    n.g1, 
#'    n.g2, 
#'    opt = TRUE
#' )
#'
#' @param min.g1 numeric value representing the sample minimum of group 1.
#' @param q1.g1 numeric value representing the first quartile of group 1.
#' @param med.g1 numeric value representing the median of group 1.
#' @param q3.g1 numeric value representing the third quartile of group 1.
#' @param max.g1 numeric value representing the sample maximum of group 1.
#' @param min.g2 numeric value representing the sample minimum of group 2.
#' @param q1.g2 numeric value representing the first quartile of group 2.
#' @param med.g2 numeric value representing the median of group 2.
#' @param q3.g2 numeric value representing the third quartile of group 2.
#' @param max.g2 numeric value representing the sample maximum of group 2.
#' @param n.g1 numeric value specifying the sample size of group 1.
#' @param n.g2 numeric value specifying the sample size of group 2.
#' @param opt logical value indicating whether to apply the optimization step in estimating 
#'      the parameters of GLD or SLD. Default is \code{TRUE}.
#'
#' @details
#' The \code{est.mean.2g} function implement the methods proposed by De Livera et al. (2024) for the two group case by incorporating shared information across the two groups 
#' to improve the accuracy of the  estimates.
#' 
#' The generalised lambda distribution (GLD) is a four parameter family of distributions defined by its quantile function under the FKML parameterisation (Freimer et al., 1988).
#' De Livera et al. propose that the GLD quantlie function can be used to approximate a sample's distribution using 5-point summaries (\eqn{S_3}). 
#' The four parameters of GLD quantile function include: a location parameter (\eqn{\lambda_1}), an inverse scale parameter (\eqn{\lambda_2}>0), and two shape parameters (\eqn{\lambda_3} and \eqn{\lambda_4}).
#' The \code{est.mean.sld.2g} function considers the case where the underlying distribution in each group has 
#' the same shape  (i.e., common \eqn{\lambda_3} and \eqn{\lambda_4}), and differ only in location and scale. Weights are used in the optimisation step in estimating 
#' \eqn{\lambda_3} and \eqn{\lambda_4} to put more emphasis on the group with the larger sample size.
#' 
#' The quantile-based skew logistic distribution (SLD), introduced by Gilchrist (2000) and further modified by van Staden and King (2015) 
#' is used to approximate the sample's distribution using 3-point summaries (\eqn{S_1} and \eqn{S_2}).
#' The SLD quantile function is defined using three parameters: a location parameter (\eqn{\lambda}), a scale parameter (\eqn{\eta}), and a skewing parameter (\eqn{\delta}).
#' In \code{est.mean.2g}, an assumption of a common skewing parameter (\eqn{\delta}) is used for the two groups, so a pooled estimate of \eqn{\delta} is computed using weights based on the sample sizes. 
#' 
#' Under each scenario, the parameters of the respective distributions are estimated by formulating and solving a series of simultaneous equations which relate the estimated quantiles 
#' with the population counterparts. The estimated mean is then obtained via integration of functions of the estimated quantile function.
#'
#' @return
#' A list containing the estimated sample means for the two groups:
#' \itemize{
#'   \item \code{mean.g1}: numeric value representing the estimated mean of group 1.
#'   \item \code{mean.g2}: numeric value representing the estimated mean of group 2.
#' }
#'
#' @examples
#' #Generate 5-point summary data for two groups
#' set.seed(123)
#' n_t <- 1000
#' n_c <- 1500
#' x_t <- stats::rlnorm(n_t, 4, 0.3)
#' x_c <- 1.1*(stats::rlnorm(n_c, 4, 0.3))
#' q_t <- c(min(x_t), stats::quantile(x_t, probs = c(0.25, 0.5, 0.75)), max(x_t))
#' q_c <- c(min(x_c), stats::quantile(x_c, probs = c(0.25, 0.5, 0.75)), max(x_c))
#' obs_mean_t <- mean(x_t)
#' obs_mean_c <- mean(x_c)
#' 
#' #Estimate sample mean using s3 (5 number summary)
#' est_means_s3 <- est.mean.2g(q_t[1],q_t[2],q_t[3],q_t[4],q_t[5],
#'                             q_c[1],q_c[2],q_c[3],q_c[4],q_c[5],
#'                             n.g1 = n_t,
#'                             n.g2 = n_c)
#' est_means_s3
#' 
#' #Estimate sample mean using s1 (min, med, max)
#' est_means_s1 <- est.mean.2g(min.g1=q_t[1], med.g1=q_t[3], max.g1=q_t[5],
#'                             min.g2=q_c[1], med.g2=q_c[3], max.g2=q_c[5],
#'                             n.g1 = n_t,
#'                             n.g2 = n_c)
#' est_means_s1
#' 
#' #Estimate sample mean using s2 (q1, med, q3)
#' est_means_s2 <- est.mean.2g(q1.g1=q_t[2], med.g1=q_t[3], q3.g1=q_t[4],
#'                             q1.g2=q_c[2], med.g2=q_c[3], q3.g2=q_c[4],
#'                             n.g1 = n_t,
#'                             n.g2 = n_c)
#' est_means_s2
#'
#' @references Alysha De Livera, Luke Prendergast, and Udara Kumaranathunga. A novel density-based approach for estimating unknown means, distribution visualisations, and meta-analyses of quantiles. \emph{Submitted for Review}, 2024, pre-print available here: <https://arxiv.org/abs/2411.10971>
#' @references Marshall Freimer, Georgia Kollia, Govind S Mudholkar, and C Thomas Lin. A study of the generalized tukey lambda family. \emph{Communications in Statistics-Theory and Methods}, 17(10):3547–3567, 1988.
#' @references Warren Gilchrist. \emph{Statistical modelling with quantile functions}. Chapman and Hall/CRC, 2000.
#' @references P. J. van Staden and R. A. R. King. The quantile-based skew logistic distribution.  \emph{Statistics & Probability Letters}, 96:109–116, 2015.
#'
#' @seealso
#' \code{\link{est.mean}} for estimating means from one-group quantile data.
#'
#' @export
#' 
#' @importFrom gld qgl gld.moments
#' @importFrom sld qsl
#' @importFrom stats integrate optim integrate runif var
#' @importFrom estmeansd bc.mean.sd qe.mean.sd


est.mean.2g <- function(
    min.g1 = NULL, q1.g1 = NULL, med.g1=NULL, q3.g1 = NULL, max.g1 = NULL,
    min.g2 = NULL, q1.g2 = NULL, med.g2=NULL, q3.g2 = NULL, max.g2 = NULL,
    n.g1, 
    n.g2, 
    opt = TRUE) {
  
  #5-number summary
  if (!is.null(min.g1) && !is.null(q1.g1) && !is.null(med.g1) && !is.null(q3.g1) && !is.null(max.g1) &&
      !is.null(min.g2) && !is.null(q1.g2) && !is.null(med.g2) && !is.null(q3.g2) && !is.null(max.g2)) {
    
      glsl_est <- est.density.five.2g(min_t = min.g1, q1_t = q1.g1, med_t=med.g1, q3_t = q3.g1, max_t = max.g1,
                                      min_c = min.g2, q1_c = q1.g2, med_c=med.g2, q3_c = q3.g2, max_c = max.g2,
                                      n_t=n.g1, 
                                      n_c=n.g2, 
                                      opt = opt)
      mean_est_g1 <- glsl_est$mean_t
      mean_est_g2 <- glsl_est$mean_c
      return(list("mean.g1" = mean_est_g1,
                  "mean.g2" = mean_est_g2))
    
  #3-number summary: min,med,max  
  } else if (!is.null(min.g1) && !is.null(med.g1) && !is.null(max.g1) &&
             !is.null(min.g2) && !is.null(med.g2) && !is.null(max.g2)) {
    
      glsl_est <- est.density.three1.2g(min_t = min.g1, med_t=med.g1, max_t = max.g1,
                                        min_c = min.g2, med_c=med.g2, max_c = max.g2,
                                        n_t=n.g1, 
                                        n_c=n.g2, 
                                        opt = opt)
      mean_est_g1 <- glsl_est$mean_t
      mean_est_g2 <- glsl_est$mean_c
      return(list("mean.g1" = mean_est_g1,
                  "mean.g2" = mean_est_g2))

  #3-number summary: q1,med,q3  
  } else if (!is.null(q1.g1) && !is.null(med.g1) && !is.null(q3.g1) &&
             !is.null(q1.g2) && !is.null(med.g2) && !is.null(q3.g2)) {
    
      glsl_est <- est.density.three2.2g(q1_t = q1.g1, med_t=med.g1, q3_t = q3.g1,
                                        q1_c = q1.g2, med_c=med.g2, q3_c = q3.g2,
                                        n_t=n.g1, 
                                        n_c=n.g2, 
                                        opt = opt)
      mean_est_g1 <- glsl_est$mean_t
      mean_est_g2 <- glsl_est$mean_c
      return(list("mean.g1" = mean_est_g1,
                  "mean.g2" = mean_est_g2))

  } else {
    stop("Invalid input: Please provide either {min, q1, med, q3, max} or {min, med, max} or {q1, med, q3}.\n")
  }
}



est.density.five.2g <- function(
    min_t = NULL, q1_t = NULL, med_t=NULL, q3_t = NULL, max_t = NULL,
    min_c = NULL, q1_c = NULL, med_c=NULL, q3_c = NULL, max_c = NULL,
    n_t, n_c, 
    opt = TRUE) {

    u_t <- 0.5 / n_t
    w_t<-n_t/(n_t+n_c)
    
    u_c <- 0.5 / n_c
    w_c<-n_c/(n_t+n_c)
    
    rho1.hat_t <- med_t
    rho2.hat_t <- max_t - min_t
    rho3.hat_t <- (med_t - min_t) / (max_t - med_t)
    rho4.hat_t <- (q3_t - q1_t) / rho2.hat_t
    
    rho1.hat_c <- med_c
    rho2.hat_c <- max_c - min_c
    rho3.hat_c <- (med_c - min_c) / (max_c - med_c)
    rho4.hat_c <- (q3_c - q1_c) / rho2.hat_c
    
    res <- optim(par = c(0.5, 0.5),
                 Fmin_2g,
                 rho3.hat_c = rho3.hat_c,
                 rho4.hat_c = rho4.hat_c,
                 u_c = u_c,
                 w_c = w_c,
                 rho3.hat_t = rho3.hat_t,
                 rho4.hat_t = rho4.hat_t,
                 u_t = u_t,
                 w_t = w_t)
    
    l3 <- res$par[1]
    l4 <- res$par[2]
    
    l2_t <- (Sfkml_l3l4(1 - u_t, l3, l4) - Sfkml_l3l4(u_t, l3, l4)) / rho2.hat_t
    l1_t <- rho1.hat_t - Sfkml_l3l4(1 / 2, l3, l4) / l2_t
    
    l2_c <- (Sfkml_l3l4(1 - u_c, l3, l4) - Sfkml_l3l4(u_c, l3, l4)) / rho2.hat_c
    l1_c <- rho1.hat_c - Sfkml_l3l4(1 / 2, l3, l4) / l2_c
    
    params_temp <- c(l1_c, l2_c, l1_t, l2_t, l3, l4)
    
    if (opt) {
      params <- optim(par = params_temp, 
                      n_c = n_c,
                      n_t = n_t,
                      emp_quant_c = c(min_c, q1_c, med_c, q3_c, max_c),
                      emp_quant_t = c(min_t, q1_t, med_t, q3_t, max_t),
                      fn = objective_function_gl_2g)$par #c(l1_c, l2_c, l1_t, l2_t, l3, l4)
    } else {
      params <- params_temp
    }
    
    names(params) <- c("location_c", "inverse scale_c", "location_t", "inverse scale_t", "shape 1", "shape 2")
    
    res_c <- true_summary("gl", list(params[c(1, 2, 5, 6)]))
    res_t <- true_summary("gl", list(params[c(3, 4, 5, 6)]))
    
    mean_est_c <- as.numeric(res_c$mean)
    mean_est_t <- as.numeric(res_t$mean)
    
    sd_est_c <- as.numeric(sqrt(res_c$variance))
    sd_est_t <- as.numeric(sqrt(res_t$variance))
    
    return(list("parameters" = params,
                "mean_t" = mean_est_t,
                "mean_c" = mean_est_c,
                "sd_t" = sd_est_t,
                "sd_c" = sd_est_c))
    
}



est.density.three1.2g <- function(
    min_t = NULL, q1_t = NULL, med_t=NULL, q3_t = NULL, max_t = NULL,
    min_c = NULL, q1_c = NULL, med_c=NULL, q3_c = NULL, max_c = NULL,
    n_t, n_c, 
    opt = TRUE) {

    a_t <- min_t
    m_t <- med_t
    b_t <- max_t
    
    v_t <- (b_t - m_t) / (m_t - a_t)
    d_t <- log((2 * n_t - 1) / n_t)
    delta_t <- max(((d_t - (log(n_t) * v_t)) / ((d_t - log(n_t)) * (v_t + 1))), 0)
    delta_t <- min(delta_t, 0.99)
    eta_t <- (b_t - a_t) / log(2 * n_t - 1)
    
    a_c <- min_c
    m_c <- med_c
    b_c <- max_c
    
    v_c <- (b_c - m_c) / (m_c - a_c)
    d_c <- log((2 * n_c - 1) / n_c)
    delta_c <- max(((d_c - (log(n_c) * v_c)) / ((d_c - log(n_c)) * (v_c + 1))), 0)
    delta_c <- min(delta_c, 0.99)
    eta_c <- (b_c - a_c) / log(2 * n_c - 1)
    
    #weighted delta
    w_c<-n_c/(n_t+n_c)
    w_t<-n_t/(n_t+n_c)
    delta_w<-w_c*delta_c+w_t*delta_t
    
    #lambda using weighted delta
    lambda_t <- m_t + eta_t * log(2) * (1 - 2*delta_w)
    lambda_c <- m_c + eta_c * log(2) * (1 - 2*delta_w)
    
    params_temp <- c(lambda_c, eta_c, lambda_t, eta_t, delta_w)
    
    if (opt) {
      params <- optim(par = params_temp, 
                      n_c = n_c,
                      n_t = n_t,
                      empirical_quantiles_c = c(min_c, med_c, max_c),
                      empirical_quantiles_t = c(min_t, med_t, max_t),
                      fn = objective_function_sl_minmax_2g)$par
    } else {
      params <- params_temp
    }
    
    names(params) <- c("location_c", "scale_c", "location_t", "scale_t", "skewing")
    
    res_c <- true_summary("sl", list(params[c(1, 2, 5)]))
    res_t <- true_summary("sl", list(params[c(3, 4, 5)]))
    
    mean_est_c <- as.numeric(res_c$mean)
    mean_est_t <- as.numeric(res_t$mean)
    
    sd_est_c <- as.numeric(sqrt(res_c$variance))
    sd_est_t <- as.numeric(sqrt(res_t$variance))
    
    return(list("parameters" = params,
                "mean_t" = mean_est_t,
                "mean_c" = mean_est_c,
                "sd_t" = sd_est_t,
                "sd_c" = sd_est_c))
    
}


est.density.three2.2g <- function(
    min_t = NULL, q1_t = NULL, med_t=NULL, q3_t = NULL, max_t = NULL,
    min_c = NULL, q1_c = NULL, med_c=NULL, q3_c = NULL, max_c = NULL,
    n_t, n_c, 
    opt = TRUE) {
  
    rho_t <- (q3_t - med_t)/(med_t-q1_t)
    delta_t <- max((log(3/2) - log(2)*rho_t)/log(3/4)/(rho_t + 1), 0)
    delta_t <- min(delta_t, 0.99)
    eta_t <- (q3_t - q1_t)/(log(3))
    
    rho_c <- (q3_c - med_c) / (med_c-q1_c)
    delta_c <- max((log(3/2) - log(2)*rho_c) / log(3/4) / (rho_c + 1), 0)
    delta_c <- min(delta_c, 0.99)
    eta_c <- (q3_c - q1_c)/(log(3))
    
    #weighted delta
    w_c<-n_c/(n_t+n_c)
    w_t<-n_t/(n_t+n_c)
    delta_w<-w_c*delta_c+w_t*delta_t
    
    #lambda using weighted delta
    lambda_t <- med_t + eta_t*log(2)*(1-2*delta_w)
    lambda_c <- med_c + eta_c*log(2)*(1-2*delta_w)
    
    params_temp <- c(lambda_c, eta_c, lambda_t, eta_t, delta_w)
    
    if (opt) {
      params <- optim(par = params_temp, 
                      n_c = n_c,
                      n_t = n_t,
                      empirical_quantiles_c = c(q1_c, med_c, q3_c),
                      empirical_quantiles_t = c(q1_t, med_t, q3_t),
                      fn = objective_function_sl_2g)$par
    } else {
      params <- params_temp
    }
    
    names(params) <- c("location_c", "scale_c", "location_t", "scale_t", "skewing")
    
    res_c <- true_summary("sl", list(params[c(1, 2, 5)]))
    res_t <- true_summary("sl", list(params[c(3, 4, 5)]))
    
    mean_est_c <- as.numeric(res_c$mean)
    mean_est_t <- as.numeric(res_t$mean)
    
    sd_est_c <- as.numeric(sqrt(res_c$variance))
    sd_est_t <- as.numeric(sqrt(res_t$variance))
    
    return(list("parameters" = params,
                "mean_t" = mean_est_t,
                "mean_c" = mean_est_c,
                "sd_t" = sd_est_t,
                "sd_c" = sd_est_c))
}


