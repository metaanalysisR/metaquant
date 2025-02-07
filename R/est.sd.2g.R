#' Estimating Sample Standard Deviations of Two Groups using Quantiles 
#'
#' @description
#' This function estimates the sample standard deviations (SD) from a two group study presenting quantile summary measures with the sample size (\eqn{n}). 
#' The quantile summaries of each group can fall into one of the following categories:
#' \itemize{
#'   \item \eqn{S_1}: \{ minimum, median, maximum \} 
#'   \item \eqn{S_2}: \{ first quartile, median, third quartile \}
#'   \item \eqn{S_3}: \{ minimum, first quartile, median, third quartile, maximum \}
#' }
#'
#' The \code{est.sd.2g} function uses a novel quantile-based distribution methods for estimating sample SD for two groups such as 'Treatment' and 'Control' (De Livera et al., 2024). 
#' The method is based on the following quantile-based distributions:
#' \itemize{
#'   \item Generalized Lambda Distribution (GLD) for estimating sample SDs using 
#'         5-number summaries (\eqn{S_3}).
#'   \item Skew Logistic Distribution (SLD) for estimating sample SDs using 
#'         3-number summaries (\eqn{S_1} and \eqn{S_2}).
#' }
#' 
#' @usage
#' est.sd.2g(
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
#' For details explaining the method of estimating using GLD or SLD, check \code{\link{est.mean.2g}}. 
#'
#' @return
#' A list containing the estimated sample SDs for the two groups:
#' \itemize{
#'   \item \code{sd.g1}: numeric value representing the estimated standard deviation of group 1.
#'   \item \code{sd.g2}: numeric value representing the estimated standard deviation of group 2.
#' }
#'
#' @examples
#' #Generate 5-point summary data for two groups
#' set.seed(123)
#' n_t <- 1000
#' n_c <- 1500
#' x_t <- stats::rlnorm(n_t, 5, 0.5)
#' x_c <- 1.1*(stats::rlnorm(n_c, 5, 0.5))
#' q_t <- c(min(x_t), stats::quantile(x_t, probs = c(0.25, 0.5, 0.75)), max(x_t))
#' q_c <- c(min(x_c), stats::quantile(x_c, probs = c(0.25, 0.5, 0.75)), max(x_c))
#' obs_sd_t <- sd(x_t)
#' obs_sd_c <- sd(x_c)
#' 
#' #Estimate sample SD using s3 (5 number summary)
#' est_sds_s3 <- est.sd.2g(q_t[1],q_t[2],q_t[3],q_t[4],q_t[5],
#'                         q_c[1],q_c[2],q_c[3],q_c[4],q_c[5],
#'                         n.g1 = n_t,
#'                         n.g2 = n_c)
#' est_sds_s3
#' 
#' #Estimate sample SD using s1 (min, med, max)
#' est_sds_s1 <- est.sd.2g(min.g1=q_t[1], med.g1=q_t[3], max.g1=q_t[5],
#'                         min.g2=q_c[1], med.g2=q_c[3], max.g2=q_c[5],
#'                         n.g1 = n_t,
#'                         n.g2 = n_c)
#' est_sds_s1
#' 
#' #Estimate sample SD using s2 (q1, med, q3)
#' est_sds_s2 <- est.sd.2g(q1.g1=q_t[2], med.g1=q_t[3], q3.g1=q_t[4],
#'                         q1.g2=q_c[2], med.g2=q_c[3], q3.g2=q_c[4],
#'                         n.g1 = n_t,
#'                         n.g2 = n_c)
#' est_sds_s2
#'
#' @references Alysha De Livera, Luke Prendergast, and Udara Kumaranathunga. A novel density-based approach for estimating unknown means, distribution visualisations, and meta-analyses of quantiles. \emph{Submitted for Review}, 2024, pre-print available here: <https://arxiv.org/abs/2411.10971>
#'
#' @seealso
#' \code{\link{est.sd}} for estimating standard deviation from one-group quantile data.
#'
#' @export


est.sd.2g <- function(
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
    sd_est_g1 <- glsl_est$sd_t
    sd_est_g2 <- glsl_est$sd_c
    return(list("sd.g1" = sd_est_g1,
                "sd.g2" = sd_est_g2))
    
    #3-number summary:min, med,max  
  } else if (!is.null(min.g1) && !is.null(med.g1) && !is.null(max.g1) &&
             !is.null(min.g2) && !is.null(med.g2) && !is.null(max.g2)) {
    
    glsl_est <- est.density.three1.2g(min_t = min.g1, med_t=med.g1, max_t = max.g1,
                            min_c = min.g2, med_c=med.g2, max_c = max.g2,
                            n_t=n.g1, 
                            n_c=n.g2, 
                            opt = opt)
    sd_est_g1 <- glsl_est$sd_t
    sd_est_g2 <- glsl_est$sd_c
    return(list("sd.g1" = sd_est_g1,
                "sd.g2" = sd_est_g2))
  
    #3-number summary:q1, med,q3  
  } else if (!is.null(q1.g1) && !is.null(med.g1) && !is.null(q3.g1) &&
             !is.null(q1.g2) && !is.null(med.g2) && !is.null(q3.g2)) {
    
    glsl_est <- est.density.three2.2g(q1_t = q1.g1, med_t=med.g1, q3_t = q3.g1,
                            q1_c = q1.g2, med_c=med.g2, q3_c = q3.g2,
                            n_t=n.g1, 
                            n_c=n.g2, 
                            opt = opt)
    sd_est_g1 <- glsl_est$sd_t
    sd_est_g2 <- glsl_est$sd_c
    return(list("sd.g1" = sd_est_g1,
                "sd.g2" = sd_est_g2))
    
  } else {
    stop("Invalid input: Please provide either {min, q1, med, q3, max} or {min, med, max} or {q1, med, q3}.\n")
  }
}
