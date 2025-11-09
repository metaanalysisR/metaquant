#' Estimating Sample Standard Deviation using Quantiles
#'
#' @description
#' This function estimates the sample standard deviation from a study presenting quantile summary measures with the sample size (\eqn{n}). The quantile summaries can fall into one of the following categories:
#' \itemize{
#'   \item \eqn{S_1}: \{ minimum, median, maximum \}
#'   \item \eqn{S_2}: \{ first quartile, median, third quartile \}
#'   \item \eqn{S_3}: \{ minimum, first quartile, median, third quartile, maximum \}
#' }
#'
#' The \code{est.sd} function implements newly proposed flexible quantile-based distribution methods for estimating sample standard deviation by De Livera et al. (2024) 
#' as well as other existing methods for estimating sample standard deviations by Shi et al. (2020) and McGrath et al. (2020).
#' 
#' 
#' @usage est.sd(
#'    min = NULL, 
#'    q1 = NULL, 
#'    med = NULL, 
#'    q3 = NULL, 
#'    max = NULL, 
#'    n = NULL, 
#'    method = "shi/wan", 
#'    opt = TRUE
#'    )
#' 
#' @param min numeric value representing the sample minimum.
#' @param q1 numeric value representing the first quartile of the sample.
#' @param med numeric value representing the median of the sample.
#' @param q3 numeric value representing the third quartile of the sample.
#' @param max numeric value representing the sample maximum.
#' @param n numeric value specifying the sample size.
#' @param method character string specifying the approach used to estimate the sample standard deviations. The options are the following:
#' \describe{
#'   \item{\code{'shi/wan'}}{The default option. Method of Shi et al. (2020).}
#'   \item{\code{'gld/sld'}}{The method proposed by De Livera et al. (2024). Estimation using the generalised lambda distribution (GLD) for 5-number summaries (\eqn{S_3}), and the skew logistic distribution (SLD) for 3-number summaries (\eqn{S_1} and \eqn{S_2}).}
#'   \item{\code{'wan'}}{The method proposed by Wan et al. (2014).} 
#'   \item{\code{'bc'}}{Box-Cox method proposed by McGrath et al. (2020).}
#'   \item{\code{'qe'}}{Quantile Matching Estimation method proposed by McGrath et al. (2020).}
#' }
#' @param opt logical value indicating whether to apply the optimisation step of \code{'gld/sld'} method, in estimating their parameters using theoretical quantiles. 
#'   The default value is \code{TRUE}.
#' 
#' @details
#' For details explaining the new method \code{'gld/sld'}, check \code{\link{est.mean}}. 
#' 
#' @return \code{sd}: numeric value representing the estimated standard deviation of the sample.
#' 
#' @examples
#' #Generate 5-point summary data
#' set.seed(123)
#' n <- 1000
#' x <- stats::rlnorm(n, 5, 0.5)
#' quants <- c(min(x), stats::quantile(x, probs = c(0.25, 0.5, 0.75)), max(x))
#' obs_sd <- sd(x)
#' 
#' #Estimate sample SD using s3 (5 number summary)
#' est_sd_s3 <- est.sd(min = quants[1], q1 = quants[2], med = quants[3], q3 = quants[4], 
#'                     max = quants[5], n=n, method = "gld/sld")
#' est_sd_s3
#' 
#' #Estimate sample SD using s1 (min, median, max)
#' est_sd_s1 <- est.sd(min = quants[1], med = quants[3], max = quants[5],
#'                     n=n, method = "gld/sld")
#' est_sd_s1
#'
#' #Estimate sample SD using s2 (q1, median, q3)
#' est_sd_s2 <- est.sd(q1 = quants[2], med = quants[3], q3 = quants[4],
#'                     n=n, method = "gld/sld")
#' est_sd_s2
#'
#' 
#' @references De Livera, A. M., Prendergast, L., & Kumaranathunga, U. (2024). A novel density-based approach for estimating unknown means, distribution visualisations and meta-analyses of quantiles. \emph{arXiv preprint arXiv:2411.10971}. <https://arxiv.org/abs/2411.10971>.
#' @references Shi, J., Luo, D., Weng, H., Zeng, X.-T., Lin, L., Chu, H., & Tong, T. (2020). Optimally estimating the sample standard deviation from the five-number summary. \emph{Research Synthesis Methods, 11}(5), 641–654.
#' @references Wan, X., Wang, W., Liu, J., & Tong, T. (2014). Estimating the sample mean and standard deviation from the sample size, median, range and/or interquartile range. \emph{BMC Medical Research Methodology, 14}, 1–13.
#' @references McGrath, S., Zhao, X., Steele, R., Thombs, B. D., Benedetti, A., & the DEPRESSD Collaboration. (2020b). Estimating the sample mean and standard deviation from commonly reported quantiles in meta-analysis. \emph{Statistical Methods in Medical Research, 29}(9), 2520–2537.
#' @export 


est.sd <- function(min = NULL, 
                   q1 = NULL, 
                   med = NULL, 
                   q3 = NULL, 
                   max = NULL, 
                   n = NULL, 
                   method = "shi/wan", 
                   opt = TRUE) {
  
  #5-number summary
  if (!is.null(min) && !is.null(q1) && !is.null(med) && !is.null(q3) && !is.null(max)) {
    
    if (method == "gld/sld") {
      glsl_est <- est.gld.five(min, q1, med, q3, max, n, opt)
      sd_est <- glsl_est$sd
      return(list("sd" = sd_est))
      
    } else if (method == "wan") {
      wan_est <- wan_using5(min, q1, med, q3, max, n)
      sd_est <- as.numeric(wan_est$sd_est)
      return(list("sd" = sd_est))
      
    } else if (method == "shi/wan") {
      sd_est <- shi_using5(min, q1, med, q3, max, n)
      return(list("sd" = sd_est))
      
    } else if (method == "bc") {
      if (any(c(min, q1, med, q3, max) <= 0)) {
        add <- abs(min(c(min, q1, med, q3, max))) + 0.5
      } else {
        add <- 0
      }
      q_bc <- c(min, q1, med, q3, max) + add
      sd_est_bc <- estmeansd::bc.mean.sd(min.val = q_bc[1], q1.val = q_bc[2], med.val = q_bc[3], q3.val = q_bc[4], max.val = q_bc[5], n = n)
      sd_est <- sd_est_bc$est.sd
      return(list("sd" = sd_est))
      
    } else if (method == "qe") {
      if (any(c(min, q1, med, q3, max) <= 0)) {
        add <- abs(min(c(min, q1, med, q3, max))) + 0.5
      } else {
        add <- 0
      }
      q_qe <- c(min, q1, med, q3, max) + add
      sd_est_qe <- estmeansd::qe.mean.sd(min.val = q_qe[1], q1.val = q_qe[2], med.val = q_qe[3], q3.val = q_qe[4], max.val = q_qe[5], n = n)
      sd_est <- sd_est_qe$est.sd
      return(list("sd" = sd_est))
      
    } else {
      stop("Unsupported method.")
    }
    
    #3-number summary:min,med,max  
  } else if (!is.null(min) && !is.null(med) && !is.null(max)) {
    
    if (method == "gld/sld") {
      glsl_est <- est.sld.minq2max(min=min, med=med, max=max, n=n, opt=opt)
      sd_est <- glsl_est$sd
      return(list("sd" = sd_est))
      
    } else if (method == "wan") {
      wan_est <- wan_using_minq2max(min, med, max, n)
      sd_est <- as.numeric(wan_est$sd_est)
      return(list("sd" = sd_est))
      
    } else if (method == "shi/wan") {
      wan_est <- wan_using_minq2max(min, med, max, n)
      sd_est <- as.numeric(wan_est$sd_est)
      return(list("sd" = sd_est))
      
    } else if (method == "bc") {
      if (any(c(min, med, max) <= 0)) {
        add <- abs(min(c(min, med, max))) + 0.5
      } else {
        add <- 0
      }
      q_bc <- c(min, med, max) + add
      sd_est_bc <- estmeansd::bc.mean.sd(min.val = q_bc[1], med.val = q_bc[2], max.val = q_bc[3], n = n)
      sd_est <- sd_est_bc$est.sd
      return(list("sd" = sd_est))
      
    } else if (method == "qe") {
      if (any(c(min, med, max) <= 0)) {
        add <- abs(min(c(min, med, max))) + 0.5
      } else {
        add <- 0
      }
      q_qe <- c(min, med, max) + add
      sd_est_qe <- estmeansd::qe.mean.sd(min.val = q_qe[1], med.val = q_qe[2], max.val = q_qe[3], n = n)
      sd_est <- sd_est_qe$est.sd
      return(list("sd" = sd_est))
      
    } else {
      stop("Unsupported method.")
    }
    
    #3-number summary:q1,med,q3  
  } else if (!is.null(q1) && !is.null(med) && !is.null(q3)) {
    
    if (method == "gld/sld") {
      glsl_est <- est.sld.q1q2q3(q1=q1, med=med, q3=q3, opt=opt)
      sd_est <- glsl_est$sd
      return(list("sd" = sd_est))
      
    } else if (method == "wan") {
      wan_est <- wan_using_q1q2q3(q1, med, q3, n)
      sd_est <- as.numeric(wan_est$sd_est)
      return(list("sd" = sd_est))
      
    } else if (method == "shi/wan") {
      wan_est <- wan_using_q1q2q3(q1, med, q3, n)
      sd_est <- as.numeric(wan_est$sd_est)
      return(list("sd" = sd_est))
      
    } else if (method == "bc") {
      if (any(c(q1, med, q3) <= 0)) {
        add <- abs(min(c(q1, med, q3))) + 0.5
      } else {
        add <- 0
      }
      q_bc <- c(q1, med, q3) + add
      sd_est_bc <- estmeansd::bc.mean.sd(q1.val = q_bc[1], med.val = q_bc[2], q3.val = q_bc[3], n = n)
      sd_est <- sd_est_bc$est.sd
      return(list("sd" = sd_est))
      
    } else if (method == "qe") {
      if (any(c(q1, med, q3) <= 0)) {
        add <- abs(min(c(q1, med, q3))) + 0.5
      } else {
        add <- 0
      }
      q_qe <- c(q1, med, q3) + add
      sd_est_qe <- estmeansd::qe.mean.sd(q1.val = q_qe[1], med.val = q_qe[2], q3.val = q_qe[3], n = n)
      sd_est <- sd_est_qe$est.sd
      return(list("sd" = sd_est))
      
    } else {
      stop("Unsupported method.")
    }
    
  } else {
    stop("Invalid input: Please provide either {min, q1, med, q3, max} or {min, med, max} or {q1, med, q3}.\n")
  }
}
