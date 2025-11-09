#' Meta-Analysis of Quantiles and Functions of Quantiles
#'
#' @description
#'
#' This function implements statistical methods for meta-analysis of quantiles 
#' and functions of quantiles for single-group and two-group studies. The 
#' function uses inverse-variance weighting to synthesise information from studies 
#' that report five-number summaries (minimum, first quartile, median, third 
#' quartile, maximum) and sample sizes—particularly useful for skewed outcomes.
#'
#' The `metaquant` function currently supports two density-based frameworks:
#' (1) a Generalized Lambda Distribution (GLD) fitted via percentile matching,
#' following De Livera et al. (2024), to estimate parameters for meta-analysis of
#' medians and other quantiles; and (2) an extension of the Quantile Estimation
#' (QE) method of McGrath et al. (2020) to additional quantiles and functions of
#' quantiles, with derived standard errors for inverse-variance pooling.
#'
#' The function facilitates meta-analyses of the following effect sizes:
#' \itemize{
#'   \item Single-group quantiles: median (\eqn{m}), first quartile (\eqn{q_1}),
#'         third quartile (\eqn{q_3}).
#'   \item Two-group differences in quantiles: difference in medians
#'         (\eqn{m_{g1}-m_{g2}}), difference in first quartiles
#'         (\eqn{q_{1g1}-q_{1g2}}), difference in third quartiles
#'         (\eqn{q_{3g1}-q_{3g2}}).
#'   \item Ratio of squared interquartile ranges (IQRs) between two groups:
#'         \eqn{r = (q_{3g1}-q_{1g1})^2 / (q_{3g2}-q_{1g2})^2}.
#' }
#'
#' Portions of this implementation are adapted from
#' \code{\link[metamedian:metamedian]{metamedian}} for the QE
#' method, and have been extended to support other quantiles and functions of
#' quantiles beyond the median.
#' 
#' @usage metaquant(
#'    data,
#'    method = "gld",
#'    effect.size.type = "median",
#'    opt = TRUE,
#'    single.family = FALSE,
#'    pool.studies = TRUE,
#'    ...)
#' 
#' @param data a data frame with one row per study containing five-number 
#' summaries and sample sizes. 
#'  For one-group studies, the input should contain the following columns:
#' \describe{
#'   \item{\code{'min.g1'}}{minimum value}
#'   \item{\code{'q1.g1'}}{first quartile}
#'   \item{\code{'med.g1'}}{median}
#'   \item{\code{'q3.g1'}}{third quartile}
#'   \item{\code{'max.g1'}}{maximum value}
#'   \item{\code{'n.g1'}}{sample size}
#' }
#' For two-group studies, also include the corresponding columns for the second
#' group: \code{min.g2}, \code{q1.g2}, \code{med.g2}, \code{q3.g2}, 
#' \code{max.g2}, and \code{n.g2}.
#' 
#' @param method character string specifying the density-based approach used to 
#' perform the meta analysis of quantiles or their functions. Options:
#' \describe{
#'   \item{\code{'gld'}}{The default option. Estimation method proposed by 
#'   De Livera et al. (2024) using the generalised lambda distribution (GLD).}
#'   \item{\code{'qe'}}{Quantile Matching Estimation method proposed by McGrath 
#'   et al. (2020).}
#' }
#' @param effect.size.type character string specifying the quantile-based effect
#'   size for the meta-analysis. Options:
#' \describe{
#'   \item{\code{'median'}}{The default option. Median for single-group studies; 
#'       difference in medians for two-group studies.}
#'   \item{\code{'q1'}}{First quartile for single-group studies; difference in
#'       first quartiles for two-group studies.}
#'   \item{\code{'q3'}}{Third quartile for single-group studies; difference in
#'       third quartiles for two-group studies.}
#'   \item{\code{'logr2'}}{Log ratio of squared IQRs between two groups. Only
#'       applicable when the input data frame provides columns for both groups.}
#' }
#' @param opt logical; whether to apply the optimisation step of the
#'   \code{"gld"} method when estimating its parameters. Default is \code{TRUE}.
#'
#' @param single.family logical; for two-group studies using the \code{"qe"}
#'   method, whether to assume the same parametric family of distributions for
#'   both groups. Default is \code{FALSE}. 
#'   See \code{\link[metamedian:qe.study.level]{qe.study.level}}.
#'
#' @param pool.studies logical; whether to pool study-specific effect sizes via
#'   inverse-variance–weighted meta-analysis. Default is \code{TRUE}. If
#'   \code{FALSE}, the function returns a list of effect sizes and their
#'   within-study variance estimates. See \code{\link[metafor:rma.uni]{rma.uni}}.
#'
#' @param ... additional arguments passed to
#'   \code{\link[metafor:rma.uni]{rma.uni}} for pooling.
#' 
#' 
#' @return An object of class "rma.uni" or a list of effect sizes and their 
#' estimated variances.
#' 
#' @seealso [est.mean()], [est.sd()] 
#' 
#' @examples
#' # Example dataset of 5-number summaries (min, q1, med, q3, max) for 2 groups
#' data_2g <- data.frame(
#'   study.index = c("Study1", "Study2", "Study3"),
#'   min.g1 = c(15, 15, 13),
#'   q1.g1  = c(57, 59, 55),
#'   med.g1 = c(66, 68, 60),
#'   q3.g1  = c(74, 72, 69),
#'   max.g1 = c(108, 101, 100),
#'   n.g1   = c(226, 230, 200),
#'   min.g2 = c(18, 19, 15),
#'   q1.g2  = c(66, 71, 69),
#'   med.g2 = c(73, 82, 81),
#'   q3.g2  = c(80, 93, 89),
#'   max.g2 = c(110, 115, 100),
#'   n.g2   = c(226, 230, 200)
#' )
#' print(data_2g)
#'
#' # Meta-analysis of difference in first quartiles
#' metaquant(data = data_2g, method = "gld", effect.size.type = "q1")
#' metaquant(data = data_2g, method = "qe",  effect.size.type = "q1")
#'
#' # Meta-analysis of log ratio of squared IQRs
#' ma_lr <- metaquant(data = data_2g, method = "gld", effect.size.type = "logr2")
#' # Back-transform to original scale (ratio of squared IQRs)
#' est_r <- exp(ma_lr$b)                                       # pooled estimate
#' ci_r  <- exp(c(ma_lr$ci.lb, ma_lr$ci.ub))                   # confidence interval
#' pi_r  <- exp(c(predict(ma_lr)$pi.lb, predict(ma_lr)$pi.ub)) # prediction interval
#' est_r; ci_r; pi_r
#'
#' @references De Livera, A. M., Prendergast, L., & Kumaranathunga, U. (2024). A novel 
#' density-based approach for estimating unknown means, distribution visualisations and meta-analyses 
#' of quantiles. \emph{arXiv preprint arXiv:2411.10971}. <https://arxiv.org/abs/2411.10971>.
#' @references King, R., Dean, B., Klinke, S., & van Staden, P. (2025). gld: Estimation 
#' and use of the Generalised (Tukey) Lambda Distribution (R package Version 2.6.7). Comprehensive R 
#' Archive Network (CRAN). https://doi.org/10.32614/CRAN.package.gld. <https://CRAN.R-project.org/package=gld>.
#' @references McGrath, S., Sohn, H., Steele, R., & Benedetti, A. (2020). Meta‐analysis of 
#' the difference of medians. \emph{Biometrical Journal}, 62(1), 69-98.
#' @references McGrath, S., Zhao, X., Ozturk, O., Katzenschlager, S., Steele, R., & Benedetti, A. (2024). 
#' Metamedian: an R package for meta‐analyzing studies reporting medians. \emph{Research Synthesis Methods}, 15(2), 332-346.
#' 
#' @export 
#' 
#' @importFrom metafor rma.uni
#'


metaquant <- function(data,
                      method = "gld",
                      effect.size.type = "median",
                      opt = TRUE,
                      single.family = FALSE,
                      pool.studies = TRUE,
                      ...) {

  if (!method %in% c("gld", "qe")) {
    stop("The 'method' argument must be 'gld' or 'qe'.")
  }
  if (!effect.size.type %in% c("median", "q1", "q3", "logr2")) {
    stop("The 'effect.size.type' argument must be 'median', 'q1', 'q3', or 'logr2'.")
  }
  
  df <- check.meta.df(data)
  
  if (effect.size.type %in% c("median", "q1", "q3")) {
    res <- mapply(est.q.study.level,
                  df$min.g1, df$q1.g1, df$med.g1, df$q3.g1, df$max.g1, df$n.g1,
                  df$min.g2, df$q1.g2, df$med.g2, df$q3.g2, df$max.g2, df$n.g2,
                  method, effect.size.type, opt, single.family)
    
    group.counts <- as.integer(unlist(res["number.of.groups", ]))
    if (length(unique(group.counts)) != 1L) {
      stop("All studies must share the same design (either all studies must have 
            one group or all must have two groups) for effect.size.type 'median', 
           'q1', or 'q3'.")
    }
    
    if (pool.studies) {
      return( metafor::rma.uni( yi = unlist(res["effect.size", ]),
                                vi = unlist(res["estvar", ]), ... ))
    } else {
      return(list( effect_sizes  = unlist(res["effect.size", ]),
                   estimated_variances = unlist(res["estvar", ])))
    }
    
  } else { # effect.size.type == "logr2"
    res <- mapply(est.r.study.level,
                  df$min.g1, df$q1.g1, df$med.g1, df$q3.g1, df$max.g1, df$n.g1,
                  df$min.g2, df$q1.g2, df$med.g2, df$q3.g2, df$max.g2, df$n.g2,
                  method, opt, single.family)
    
    group.counts <- as.integer(unlist(res["number.of.groups", ]))
    if (length(unique(group.counts)) != 1L) {
      stop("All studies must have two groups for effect.size.type 'logr2'.")
    }
    
    if (pool.studies) {
      return(metafor::rma.uni(yi = unlist(res["effect.size.log", ]),
                              vi = unlist(res["estvar.log", ]), ... ))
    } else {
      return(list( effect_sizes  = unlist(res["effect.size.log", ]),
                   estimated_variances = unlist(res["estvar.log", ])))
    }
  }
}


#' Estimating Quantile-Based Effect Sizes and Variances 
#'
#' @description
#'
#' This function estimates the variances of quantiles and the differences of 
#' quantiles for single-group and two-group studies, respectively, from studies 
#' that report five-number summaries (minimum, first quartile, median, third 
#' quartile, maximum) and sample sizes, using density-based approaches.
#'
#' The `est.q.study.level` function currently supports two density-based frameworks:
#' (1) a Generalized Lambda Distribution (GLD) fitted via percentile matching,
#' following De Livera et al. (2024); and (2) an extension of the Quantile Estimation
#' (QE) method of McGrath et al. (2020) to additional quantiles and functions of
#' quantiles.
#'
#' The function estimates the asymptotic variances of the following effect sizes:
#' \itemize{
#'   \item Single-group quantiles: median (\eqn{m}), first quartile (\eqn{q_1}),
#'         third quartile (\eqn{q_3}).
#'   \item Two-group differences in quantiles: difference in medians
#'         (\eqn{m_{g1}-m_{g2}}), difference in first quartiles
#'         (\eqn{q_{1g1}-q_{1g2}}), difference in third quartiles
#'         (\eqn{q_{3g1}-q_{3g2}}).
#' }
#'
#' Portions of this implementation are adapted from
#' \code{\link[metamedian:qe.study.level]{qe.study.level}}
#' for the QE method, and have been extended to support other quantiles and 
#' functions of quantiles beyond the median.
#' 
#' @usage est.q.study.level(
#'    min.g1, 
#'    q1.g1, 
#'    med.g1, 
#'    q3.g1, 
#'    max.g1, 
#'    n.g1, 
#'    min.g2, 
#'    q1.g2, 
#'    med.g2, 
#'    q3.g2, 
#'    max.g2, 
#'    n.g2,
#'    method, 
#'    effect.size.type, 
#'    opt = TRUE, 
#'    single.family = FALSE, 
#'    qe.fit.control.g1 = list(), 
#'    qe.fit.control.g2 = list()
#'  )
#' 
#' @param min.g1 numeric value representing the sample minimum (of group one for two-group studies).
#' @param q1.g1 numeric value representing the first quartile of the sample (of group one for two-group studies).
#' @param med.g1 numeric value representing the median of the sample (of group one for two-group studies).
#' @param q3.g1 numeric value representing the third quartile of the sample (of group one for two-group studies).
#' @param max.g1 numeric value representing the sample maximum (of group one for two-group studies).
#' @param n.g1 numeric value specifying the sample size (of group one for two-group studies).
#' @param min.g2 numeric value representing the sample minimum of group two for two-group studies.
#' @param q1.g2 numeric value representing the first quartile of the sample of group two for two-group studies.
#' @param med.g2 numeric value representing the median of the sample of group two for two-group studies.
#' @param q3.g2 numeric value representing the third quartile of the sample of group two for two-group studies.
#' @param max.g2 numeric value representing the sample maximum of group two for two-group studies.
#' @param n.g2 numeric value specifying the sample size of group two for two-group studies.
#' 
#' @param method character string specifying the density-based approach used to 
#' estimate variances of quantiles or their functions. Options:
#' \describe{
#'   \item{\code{'gld'}}{The default option. Estimation method proposed by 
#'   De Livera et al. (2024) using the generalised lambda distribution (GLD).}
#'   \item{\code{'qe'}}{Quantile Matching Estimation method proposed by McGrath 
#'   et al. (2020).}
#' }
#' @param effect.size.type character string specifying the quantile-based effect
#'   size for the meta-analysis. Options:
#' \describe{
#'   \item{\code{'median'}}{The default option. Median for single-group studies; 
#'       difference in medians for two-group studies.}
#'   \item{\code{'q1'}}{First quartile for single-group studies; difference in
#'       first quartiles for two-group studies.}
#'   \item{\code{'q3'}}{Third quartile for single-group studies; difference in
#'       third quartiles for two-group studies.}
#' }
#' 
#' @param opt logical; whether to apply the optimisation step of the
#'   \code{"gld"} method when estimating its parameters. Default is \code{TRUE}.
#' @param single.family logical; for two-group studies using the \code{"qe"}
#'   method, whether to assume the same parametric family of distributions for
#'   both groups. Default is \code{FALSE}. 
#'   See \code{\link[metamedian:qe.study.level]{qe.study.level}}
#'   
#' @param qe.fit.control.g1	 optional list of control parameters for 
#'   \code{\link[estmeansd:qe.fit]{qe.fit}} (of group one for two-group studies).
#' @param qe.fit.control.g2	 optional list of control parameters for 
#'   \code{\link[estmeansd:qe.fit]{qe.fit}} of group two for two-group studies.
#'
#' @return A list containing following components:
#' - \code{effect.size}: numeric value of quantile-based effect size of the study 
#' based on the input of \code{effect.size.type} argument.
#' - \code{estvar}: numeric value of the estimated variance of the effect size.
#' - \code{number.of.groups}: integer indicating the number of groups 
#' in the input study data.
#' - \code{effect.size.name}: character string specifying a label for the effect 
#' size depending on \code{number.of.groups} and \code{effect.size.type}.
#' 
#' @seealso [est.r.study.level()]
#' 
#' @examples
#' #Generate 5-number summary data (group one)
#' set.seed(123)
#' n1 <- 100
#' x1 <- stats::rlnorm(n1, 4, 0.3)
#' quants1 <- c(min(x1), stats::quantile(x1, probs = c(0.25, 0.5, 0.75)), max(x1))
#' 
#' #Estimate variance of the first quartile
#' est.q.study.level(min.g1 = quants1[1], q1.g1 = quants1[2], med.g1 = quants1[3], 
#'                   q3.g1 = quants1[4],max.g1 = quants1[5], n.g1=n1, 
#'                   method = "gld", effect.size.type = "q1")
#'
#' #Generate 5-number summary data (group two)
#' set.seed(123) 
#' n2 <- 120
#' x2 <- stats::rlnorm(n2, 3, 0.5)
#' quants2 <- c(min(x2), stats::quantile(x2, probs = c(0.25, 0.5, 0.75)), max(x2))
#' 
#' #Estimate variance of the difference in first quartiles (for two groups)
#' est.q.study.level(min.g1 = quants1[1], q1.g1 = quants1[2], med.g1 = quants1[3], 
#'                   q3.g1 = quants1[4], max.g1 = quants1[5], n.g1=n1, 
#'                   min.g2 = quants2[1], q1.g2 = quants2[2], med.g2 = quants2[3], 
#'                   q3.g2 = quants2[4],  max.g2 = quants2[5], n.g2=n2, 
#'                   method = "gld", effect.size.type = "q1")
#'                  
#'                  
#' @references De Livera, A. M., Prendergast, L., & Kumaranathunga, U. (2024). A novel 
#' density-based approach for estimating unknown means, distribution visualisations and meta-analyses 
#' of quantiles. \emph{arXiv preprint arXiv:2411.10971}. <https://arxiv.org/abs/2411.10971>.
#' @references King, R., Dean, B., Klinke, S., & van Staden, P. (2025). gld: Estimation 
#' and use of the Generalised (Tukey) Lambda Distribution (R package Version 2.6.7). Comprehensive R 
#' Archive Network (CRAN). https://doi.org/10.32614/CRAN.package.gld. <https://CRAN.R-project.org/package=gld>.
#' @references McGrath, S., Sohn, H., Steele, R., & Benedetti, A. (2020). Meta‐analysis of 
#' the difference of medians. \emph{Biometrical Journal}, 62(1), 69-98.
#' @references McGrath, S., Zhao, X., Ozturk, O., Katzenschlager, S., Steele, R., & Benedetti, A. (2024). 
#' Metamedian: an R package for meta‐analyzing studies reporting medians. \emph{Research Synthesis Methods}, 15(2), 332-346.
#' 
#' @export
#' 
#' @importFrom estmeansd qe.fit
#' 

################################################################################
# This function contains code adapted from the qe.study.level() function of 
# 'metamedian' R package (https://CRAN.R-project.org/package=metamedian)
################################################################################

est.q.study.level <- function(min.g1, q1.g1, med.g1, q3.g1, max.g1, n.g1, 
                              min.g2, q1.g2, med.g2, q3.g2, max.g2, n.g2,
                              method, effect.size.type, opt = TRUE, single.family = FALSE, 
                              qe.fit.control.g1 = list(), qe.fit.control.g2 = list()) {
  
  absent <- function(x) missing(x) || isTRUE(is.na(x))
  one.group <- absent(min.g2)  && absent(q1.g2) && absent(med.g2) && 
               absent(q3.g2) && absent(max.g2)
  
  if (effect.size.type == "median") {
    p=0.5
    effect.size.g1 <- med.g1
  } else if (effect.size.type == "q1") {
    p=0.25
    effect.size.g1 <- q1.g1
  } else if (effect.size.type == "q3"){
    p=0.75
    effect.size.g1 <- q3.g1
  }
  
  if (method == "gld") {
    gld.fit.g1 <- est.gld.five(min = min.g1, q1 = q1.g1, med = med.g1, 
                                   q3 = q3.g1,  max = max.g1, n = n.g1, opt = opt) 
    fp.g1<- get.fp.gld(gld.fit.g1, p=p)$est_f

  
  } else {
    qe.fit.g1 <- estmeansd::qe.fit(min.val = min.g1, q1.val = q1.g1, med.val = med.g1, 
                                   q3.val = q3.g1, max.val = max.g1, n = n.g1,
                                   qe.fit.control = qe.fit.control.g1,
                                   two.sample.default = TRUE)
    qe.dist.g1 <- names(which.min(qe.fit.g1$values))
    fp.g1<- get.fp.qe(qe.fit.g1, qe.dist.g1, p=p)$est_f
  }
  
  estvar.g1 <- p * (1-p) / n.g1 / fp.g1^2
  
  if (one.group) {
    effect.size <- effect.size.g1
    estvar <- estvar.g1 
  } else {
      if (effect.size.type == "median") {
        p=0.5
        effect.size.g2 <- med.g2
      } else if (effect.size.type == "q1") {
        p=0.25
        effect.size.g2 <- q1.g2
      } else if (effect.size.type == "q3"){
        p=0.75
        effect.size.g2 <- q3.g2
      }
      
      if (method == "gld") {
        gld.fit.g2 <- est.gld.five(min = min.g2, q1 = q1.g2, med = med.g2, 
                                   q3 = q3.g2,  max = max.g2, n = n.g2, opt = opt) 
        fp.g2<- get.fp.gld(gld.fit.g2, p=p)$est_f
        
      } else {
        qe.fit.g2 <- estmeansd::qe.fit(min.val = min.g2, q1.val = q1.g2,
                                    med.val = med.g2, q3.val = q3.g2,
                                    max.val = max.g2, n = n.g2,
                                    qe.fit.control = qe.fit.control.g2,
                                    two.sample.default = TRUE)
        qe.dist.g2 <- names(which.min(qe.fit.g2$values))
        fp.g2<- get.fp.qe(qe.fit.g2, qe.dist.g2, p=p)$est_f
      }
      estvar.g2 <- p * (1-p) / n.g2 / fp.g2^2
      
    effect.size <- effect.size.g1 - effect.size.g2
    estvar <- estvar.g1 + estvar.g2
    
    if ( method == "qe" & single.family ){

        qe.dist <- names(which.min(qe.fit.g1$values + qe.fit.g2$values))
        fp.qe.g1 <- get.fp.qe(qe.fit.g1, qe.dist, p=p)$est_f
        fp.qe.g2 <- get.fp.qe(qe.fit.g2, qe.dist, p=p)$est_f
        
        estvar <- p * (1-p) * (1 / (fp.qe.g1^2 * n.g1) +
                                 1 / (fp.qe.g2^2 * n.g2))
    }
    
  } 
  
  number.of.groups <- if (one.group) 1L else 2L
  
  effect.size.name <- if (number.of.groups == 1L) {
    switch(effect.size.type,
           median = "Median",
           q1     = "First quartile",
           q3     = "Third quartile")
  } else {
    switch(effect.size.type,
           median = "Difference in medians",
           q1     = "Difference in first quartiles",
           q3     = "Difference in third quartiles")
  }
  
  return(list(effect.size = effect.size,
              estvar = estvar,
              number.of.groups = number.of.groups,
              effect.size.name = effect.size.name))
}


#' Estimating Variances of Squared IQR Ratio and its Natural Logarithm 
#'
#' @description
#'
#' This function estimates the variances of squared IQR ratio and its logarithm
#' for two-group studies, from studies that report five-number summaries 
#' (minimum, first quartile, median, third quartile, maximum) and sample sizes, 
#' using density-based approaches.
#'
#' The `est.r.study.level` function currently supports two density-based frameworks:
#' (1) a Generalized Lambda Distribution (GLD) fitted via percentile matching,
#' following De Livera et al. (2024); and (2) an extension of the Quantile Estimation
#' (QE) method of McGrath et al. (2020) to additional quantiles and functions of
#' quantiles.
#'
#' The function estimates the asymptotic variances of the following effect sizes:
#' \itemize{
#'   \item Ratio of squared interquartile ranges (IQRs) between two groups:
#'         \eqn{r = (q_{3g1}-q_{1g1})^2 / (q_{3g2}-q_{1g2})^2}.
#'   \item Log ratio of squared IQRs between two groups:
#'         \eqn{log(r)}.
#'         
#' }
#'
#' Portions of this implementation are adapted from
#' \code{\link[metamedian:qe.study.level]{qe.study.level}} 
#' for the QE method, and have been extended to support functions of quantiles 
#' beyond the median.
#' 
#' @usage est.r.study.level(
#'    min.g1, 
#'    q1.g1, 
#'    med.g1, 
#'    q3.g1, 
#'    max.g1, 
#'    n.g1, 
#'    min.g2, 
#'    q1.g2, 
#'    med.g2, 
#'    q3.g2, 
#'    max.g2, 
#'    n.g2,
#'    method, 
#'    opt = TRUE, 
#'    single.family = FALSE, 
#'    qe.fit.control.g1 = list(), 
#'    qe.fit.control.g2 = list()
#'  )
#' 
#' @param min.g1 numeric value representing the sample minimum (of group one for two-group studies).
#' @param q1.g1 numeric value representing the first quartile of the sample (of group one for two-group studies).
#' @param med.g1 numeric value representing the median of the sample (of group one for two-group studies).
#' @param q3.g1 numeric value representing the third quartile of the sample (of group one for two-group studies).
#' @param max.g1 numeric value representing the sample maximum (of group one for two-group studies).
#' @param n.g1 numeric value specifying the sample size (of group one for two-group studies).
#' @param min.g2 numeric value representing the sample minimum of group two for two-group studies.
#' @param q1.g2 numeric value representing the first quartile of the sample of group two for two-group studies.
#' @param med.g2 numeric value representing the median of the sample of group two for two-group studies.
#' @param q3.g2 numeric value representing the third quartile of the sample of group two for two-group studies.
#' @param max.g2 numeric value representing the sample maximum of group two for two-group studies.
#' @param n.g2 numeric value specifying the sample size of group two for two-group studies.
#' 
#' @param method character string specifying the density-based approach used to 
#' estimate variances of squared IQR ratio and its natural logarithm. Options:
#' \describe{
#'   \item{\code{'gld'}}{The default option. Estimation method proposed by 
#'   De Livera et al. (2024) using the generalised lambda distribution (GLD).}
#'   \item{\code{'qe'}}{Quantile Matching Estimation method proposed by McGrath 
#'   et al. (2020).}
#' }
#' 
#' @param opt logical; whether to apply the optimisation step of the
#'   \code{"gld"} method when estimating its parameters. Default is \code{TRUE}.
#' @param single.family logical; for two-group studies using the \code{"qe"}
#'   method, whether to assume the same parametric family of distributions for
#'   both groups. Default is \code{FALSE}. 
#'   See \code{\link[metamedian:qe.study.level]{qe.study.level}}
#'   
#' @param qe.fit.control.g1	 optional list of control parameters for 
#'   \code{\link[estmeansd:qe.fit]{qe.fit}} (of group one for two-group studies).
#' @param qe.fit.control.g2	 optional list of control parameters for 
#'   \code{\link[estmeansd:qe.fit]{qe.fit}} of group two for two-group studies.
#'
#' @return A list containing following components:
#' - \code{effect.size}: numeric value of the effect size of the study (ratio of squared IQRs).
#' - \code{estvar}: estimated variance of the effect size (ratio of squared IQRs).
#' - \code{effect.size.log}: numeric value of log ratio of squared IQRs.
#' - \code{estvar.log}: estimated variance of log ratio of squared IQRs.
#' - \code{number.of.groups}: integer indicating the number of groups in the
#'  input study data.
#' 
#' @seealso [est.q.study.level()]
#' 
#' @examples
#' #Generate 5-number summary data (group one)
#' set.seed(123)
#' n1 <- 100
#' x1 <- stats::rlnorm(n1, 4, 0.3)
#' quants1 <- c(min(x1), stats::quantile(x1, probs = c(0.25, 0.5, 0.75)), max(x1))
#' 
#' #Generate 5-number summary data (group two)
#' set.seed(123) 
#' n2 <- 120
#' x2 <- stats::rlnorm(n2, 3, 0.5)
#' quants2 <- c(min(x2), stats::quantile(x2, probs = c(0.25, 0.5, 0.75)), max(x2))
#' 
#' #Estimate variance of the squared IQR ratio and its natural logarithm (for two groups)
#' est.r.study.level(min.g1 = quants1[1], q1.g1 = quants1[2], med.g1 = quants1[3], 
#'                   q3.g1 = quants1[4], max.g1 = quants1[5], n.g1=n1, 
#'                   min.g2 = quants2[1], q1.g2 = quants2[2], med.g2 = quants2[3], 
#'                   q3.g2 = quants2[4],  max.g2 = quants2[5], n.g2=n2, 
#'                   method = "gld")
#'                  
#'                  
#' @references De Livera, A. M., Prendergast, L., & Kumaranathunga, U. (2024). A novel 
#' density-based approach for estimating unknown means, distribution visualisations and meta-analyses 
#' of quantiles. \emph{arXiv preprint arXiv:2411.10971}. <https://arxiv.org/abs/2411.10971>.
#' @references King, R., Dean, B., Klinke, S., & van Staden, P. (2025). gld: Estimation 
#' and use of the Generalised (Tukey) Lambda Distribution (R package Version 2.6.7). Comprehensive R 
#' Archive Network (CRAN). https://doi.org/10.32614/CRAN.package.gld. <https://CRAN.R-project.org/package=gld>.
#' @references McGrath, S., Sohn, H., Steele, R., & Benedetti, A. (2020). Meta‐analysis of 
#' the difference of medians. \emph{Biometrical Journal}, 62(1), 69-98.
#' @references McGrath, S., Zhao, X., Ozturk, O., Katzenschlager, S., Steele, R., & Benedetti, A. (2024). 
#' Metamedian: an R package for meta‐analyzing studies reporting medians. \emph{Research Synthesis Methods}, 15(2), 332-346.
#' 
#' @export
#' 
#' @importFrom estmeansd qe.fit
#' 

est.r.study.level <- function(min.g1, q1.g1, med.g1, q3.g1, max.g1, n.g1,
                            min.g2, q1.g2, med.g2, q3.g2, max.g2, n.g2,
                            method, opt = TRUE, single.family = FALSE, 
                            qe.fit.control.g1 = list(), qe.fit.control.g2 = list()) {
  
  absent <- function(x) missing(x) || isTRUE(is.na(x))
  one.group <- absent(min.g2)  && absent(q1.g2) && absent(med.g2) && 
               absent(q3.g2) && absent(max.g2)

  if (one.group) {
    stop("Studies must have two groups for effect.size.type of 'logr2'")
    
  } else {

    r2 <- (q3.g1-q1.g1)^2 / (q3.g2-q1.g2)^2
    effect.size <- r2
      
    if (method == "gld") {

      gld.fit.g1 <- est.gld.five(min = min.g1, q1 = q1.g1, med = med.g1, 
                                              q3 = q3.g1,  max = max.g1, n = n.g1, opt = opt) 
      gld.fit.g2 <- est.gld.five(min = min.g2, q1 = q1.g2, med = med.g2, 
                                              q3 = q3.g2,  max = max.g2, n = n.g2, opt = opt) 
      
      fq1.g1<- get.fp.gld(gld.fit.g1, p=0.25)$est_f
      fq3.g1<- get.fp.gld(gld.fit.g1, p=0.75)$est_f
      
      fq1.g2<- get.fp.gld(gld.fit.g2, p=0.25)$est_f
      fq3.g2<- get.fp.gld(gld.fit.g2, p=0.75)$est_f

      
    } else {

      qe.fit.g1 <- estmeansd::qe.fit(min.val = min.g1, q1.val = q1.g1,
                                     med.val = med.g1, q3.val = q3.g1,
                                     max.val = max.g1, n = n.g1,
                                     qe.fit.control = qe.fit.control.g1,
                                     two.sample.default = TRUE)
      qe.fit.g2 <- estmeansd::qe.fit(min.val = min.g2, q1.val = q1.g2,
                                     med.val = med.g2, q3.val = q3.g2,
                                     max.val = max.g2, n = n.g2,
                                     qe.fit.control = qe.fit.control.g2,
                                     two.sample.default = TRUE)
      
      qe.dist.g1 <- names(which.min(qe.fit.g1$values))
      qe.dist.g2 <- names(which.min(qe.fit.g2$values))
      
      fq1.g1<- get.fp.qe(qe.fit.g1, qe.dist.g1, p=0.25)$est_f
      fq3.g1<- get.fp.qe(qe.fit.g1, qe.dist.g1, p=0.75)$est_f
      
      fq1.g2<- get.fp.qe(qe.fit.g2, qe.dist.g2, p=0.25)$est_f
      fq3.g2<- get.fp.qe(qe.fit.g2, qe.dist.g2, p=0.75)$est_f
      
    }
    
    gq1.g1 <- 1/fq1.g1 ;  gq3.g1 <- 1/fq3.g1
    gq1.g2 <- 1/fq1.g2 ;  gq3.g2 <- 1/fq3.g2
    
    w1 <- n.g1/(n.g1+n.g2)
    w2 <- n.g2/(n.g1+n.g2)
    
    estvar <- (4*0.25*r2^2 / (n.g1 + n.g2)) * (
      ((gq1.g1^2 + gq3.g1^2 - 0.25*(gq1.g1 + gq3.g1)^2) / ( w1 * (q3.g1 - q1.g1)^2 )) +
        ((gq1.g2^2 + gq3.g2^2 - 0.25*(gq1.g2 + gq3.g2)^2) / ( w2 * (q3.g2 - q1.g2)^2 )))

    if (method == "qe" & single.family ){

        qe.dist <- names(which.min(qe.fit.g1$values + qe.fit.g2$values))
        
        fq1.qe.g1<- get.fp.qe(qe.fit.g1, qe.dist, p=0.25)$est_f
        fq3.qe.g1<- get.fp.qe(qe.fit.g1, qe.dist, p=0.75)$est_f
        fq1.qe.g2<- get.fp.qe(qe.fit.g2, qe.dist, p=0.25)$est_f
        fq3.qe.g2<- get.fp.qe(qe.fit.g2, qe.dist, p=0.75)$est_f
        
        gq1.qe.g1 <- 1/fq1.qe.g1 ;  gq3.qe.g1 <- 1/fq3.qe.g1
        gq1.qe.g2 <- 1/fq1.qe.g2 ;  gq3.qe.g2 <- 1/fq3.qe.g2
        
        w1 <- n.g1/(n.g1+n.g2)
        w2 <- n.g2/(n.g1+n.g2)
        
        estvar <- (4*0.25*r2^2 / (n.g1 + n.g2)) * (
          ((gq1.qe.g1^2 + gq3.qe.g1^2 - 0.25*(gq1.qe.g1 + gq3.qe.g1)^2) / ( w1 * (q3.g1 - q1.g1)^2 )) +
            ((gq1.qe.g2^2 + gq3.qe.g2^2 - 0.25*(gq1.qe.g2 + gq3.qe.g2)^2) / ( w2 * (q3.g2 - q1.g2)^2 )))
    }
    
    effect.size.log <- log(r2)
    estvar.log <- estvar/r2^2
    
    number.of.groups <- if (one.group) 1L else 2L
    
  } 
  
  return(list(effect.size = effect.size, estvar = estvar, 
              effect.size.log = effect.size.log, estvar.log = estvar.log, 
              number.of.groups= number.of.groups))
}



