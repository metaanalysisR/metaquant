#' Visualising Densities using Quantiles
#'
#' @description
#' The function estimates and visualizes the density curves of one-group or two-group studies presenting quantile summary measures with the sample size (\eqn{n}). The quantile summaries can fall into one of the following categories:
#' \itemize{
#'   \item \eqn{S_1}: \{ minimum, median, maximum \}
#'   \item \eqn{S_2}: \{ first quartile, median, third quartile \}
#'   \item \eqn{S_3}: \{ minimum, first quartile, median, third quartile, maximum \}
#' }
#'
#' The \code{plotdist} function uses the following quantile-based distribution methods for visualising densities using qantiles (De Livera et al., 2024). 
#' \itemize{
#'   \item Generalized Lambda Distribution (GLD) when 5-number summaries present (\eqn{S_3}).
#'   \item Skew Logistic Distribution (SLD) when 3-number summaries present (\eqn{S_1} and \eqn{S_2}).
#' } 
#' 
#' 
#' @usage plotdist(
#'    data, 
#'    xmin, 
#'    xmax, 
#'    length.out = 1000, 
#'    title = "", 
#'    xlab = "x", 
#'    ylab = "Density",
#'    line.size = 0.5,
#'    title.size = 12,
#'    lab.size = 10,
#'    color.g1 = "red", 
#'    color.g2 = "blue",
#'    label.g1 = NULL, 
#'    label.g2 = NULL,
#'    display.index = FALSE,
#'    display.legend = FALSE,
#'    opt = TRUE
#' )
#'
#' @param data data frame containing the quantile summary data. For one-group studies, the input dataset can contain the following columns:
#' \describe{
#'   \item{\code{'stduy.index'}}{stduy index or name}
#'   \item{\code{'min.g1'}}{minimum value}
#'   \item{\code{'q1.g1'}}{first quartile}
#'   \item{\code{'med.g1'}}{median}
#'   \item{\code{'q3.g1'}}{third quartile}
#'   \item{\code{'max.g1'}}{maximum value}
#' } For two-group studies, the data frame can also contain the following columns for the second group: \code{min.g2, q1.g2, med.g2, q3.g2, max.g2} and \code{ n.g2}.
#' @param xmin numeric value for the lower limit of the x-axis for density calculation. It is recommended to set this to a value smaller than the smallest value across the quantile summaries to ensure the density curve is fully captured.
#' @param xmax numeric value for the upper limit of the x-axis for density calculation. It is recommended to set this to a value larger than the largest value across the quantile summaries to ensure the density curve is fully captured.
#' @param length.out integer specifying the number of points along the x-axis for density calculation. Default is \code{1000}.
#' @param title character string for the plot title. Default is an empty string.
#' @param xlab character string for the x-axis label. Default is \code{"x"}. 
#' @param ylab character string for the y-axis label. Default is \code{"Density"}.
#' @param line.size numeric. Thickness of the density curve lines. Default is \code{0.5}.
#' @param title.size numeric. Font size for the plot title. Default is \code{12}.
#' @param lab.size numeric. Font size for axis labels. Default is \code{10}.
#' @param color.g1 character string specifying the color for density curves of group 1. Default is \code{"red"}.
#' @param color.g2 character string specifying the color for density curves of group 2. Default is \code{"blue"}.
#' @param label.g1 character string indicating label or name for group 1 (eg., 'Treatment')
#' @param label.g2 character string indicating label or name for group 2 (eg., 'Control'). 
#'   
#'    If \code{'label.g1'} and \code{'label.g2'} are not provided, the function will assign labels as 'Group 1' and 'Group 2'.
#' @param display.index logical. If \code{TRUE}, the \code{'study.index'} of each quantile set (row) will be displayed alongside the corresponding density curve. The default is \code{FALSE}, meaning no labels will be shown. The label text size is controlled by the \code{lab.size} parameter.
#' @param display.legend logical. If \code{TRUE}, legends (\code{'label.g1'} and/or \code{'label.g2'}) will be displayed on the right side of the plot. The default is \code{FALSE}. The legend text size is controlled by the \code{lab.size} parameter.
#' @param opt logical value indicating whether to apply the optimization step when estimating GLD or SLD parameters. The default value is \code{TRUE}.
#'   
#' @details
#' The generalised lambda distribution (GLD) is a four parameter family of distributions defined by its quantile function under the FKML parameterisation (Freimer et al., 1988).
#' De Livera et al. propose that the GLD quantlie function can be used to approximate a sample's distribution using 5-point summaries. 
#' The four parameters of GLD quantile function include: a location parameter (\eqn{\lambda_1}), an inverse scale parameter (\eqn{\lambda_2}>0), and two shape parameters (\eqn{\lambda_3} and \eqn{\lambda_4}).
#' 
#' The quantile-based skew logistic distribution (SLD), introduced by Gilchrist (2000) and further modified by van Staden and King (2015) 
#' is used to approximate the sample's distribution using 3-point summaries.
#' The SLD quantile function is defined using three parameters: a location parameter (\eqn{\lambda}), a scale parameter (\eqn{\eta}), and a skewing parameter (\eqn{\delta}).
#' 
#' These parameters of GLD and SLD are estimated by formulating and solving a series of simultaneous equations which relate the estimated quantiles
#' with the population counterparts of respective distribution (GLD or SLD). The \code{plotdist} uses these estimated parameters, to compute the density data 
#' using \code{\link{dgl}} function in the "gld" package and \code{\link{dsl}} function in the "sld" package, based on the 5-number summaries and 3-number summaries, respectively. 
#' 
#'
#' @return A ggplot object visualizing estimated densitiy curve/s for one or two groups.
#'
#' @examples 
#' #Example dataset of 3-point summaries (min, med, max) for 2 groups
#' data_3num_2g <- data.frame(
#'   study.index = c("Study1", "Study2", "Study3"),
#'   min.g1 = c(15, 15, 13),
#'   med.g1 = c(66, 68, 63),
#'   max.g1 = c(108, 101, 100),
#'   n.g1 = c(226, 230, 200),
#'   min.g2 = c(18, 19, 15),
#'   med.g2 = c(73, 82, 81),
#'   max.g2 = c(110, 115, 100),
#'   n.g2 = c(226, 230, 200)
#'  )
#' print(data_3num_2g)
#' 
#' plot_3num_2g <- plotdist(
#'   data_3num_2g,
#'   xmin = 10,
#'   xmax = 125,
#'   title = "Example Density Plots of Two Groups",
#'   xlab = "x data",
#'   color.g1 = "blue",
#'   color.g2 = "red",
#'   label.g1 = "Treatment", 
#'   label.g2 = "Control",
#'   display.index = TRUE,
#'   display.legend = TRUE
#' )
#' print(plot_3num_2g)
#'
#' @references Alysha De Livera, Luke Prendergast, and Udara Kumaranathunga. A novel density-based approach for estimating unknown means, distribution visualisations, and meta-analyses of quantiles. \emph{Submitted for Review}, 2024. (Article available on request to authors.)
#' @references Marshall Freimer, Georgia Kollia, Govind S Mudholkar, and C Thomas Lin. A study of the generalized tukey lambda family. \emph{Communications in Statistics-Theory and Methods}, 17(10):3547–3567, 1988.
#' @references Warren Gilchrist. \emph{Statistical modelling with quantile functions}. Chapman and Hall/CRC, 2000.
#' @references P. J. van Staden and R. A. R. King. The quantile-based skew logistic distribution.  \emph{Statistics & Probability Letters}, 96:109–116, 2015.
#' @export 
#' 
#' @importFrom gld dgl
#' @importFrom sld dsl
#' @importFrom ggplot2 ggplot aes geom_line scale_color_manual labs theme_minimal theme element_text geom_text
#' @importFrom stats setNames
#' @importFrom magrittr %>%
#' @importFrom dplyr group_by slice_max ungroup 


plotdist <- function(data, 
                     xmin, 
                     xmax, 
                     length.out = 1000, 
                     title = "", 
                     xlab = "x", 
                     ylab = "Density",
                     line.size = 0.5,
                     title.size = 12,
                     lab.size = 10,
                     color.g1 = "red", 
                     color.g2 = "blue",
                     label.g1 = NULL, 
                     label.g2 = NULL,
                     display.index = FALSE,
                     display.legend = FALSE,
                     opt = TRUE) {
  
  check_result <- check.df(data)
  data <- check_result$df
  invalid_rows_g1 <- check_result$invalid_rows_g1
  invalid_rows_g2 <- check_result$invalid_rows_g2
  
  g1_type <- get.scenario(data, "g1")
  g2_present <- any(grepl("\\.g2$", colnames(data)))
  g2_type <- if (g2_present) get.scenario(data, "g2") else NULL
  
  if (is.null(label.g1)) label.g1 <- "Group 1"
  if (g2_present && is.null(label.g2)) label.g2 <- "Group 2"
  
  density_data <- data.frame()
  
  valid_data_g1 <- if (length(invalid_rows_g1) > 0) {
    data[-invalid_rows_g1, c("study.index", grep("\\.g1$", colnames(data), value = TRUE))]
  } else {
    data[, c("study.index", grep("\\.g1$", colnames(data), value = TRUE))]
  }

  density_data <- rbind(density_data, 
                        get.density(valid_data_g1, "g1", g1_type, label.g1, xmin, xmax, length.out, opt))
  
  if (g2_present) {
    valid_data_g2 <- if (length(invalid_rows_g2) > 0) {
      data[-invalid_rows_g2, c("study.index", grep("\\.g2$", colnames(data), value = TRUE))]
    } else {
      data[, c("study.index", grep("\\.g2$", colnames(data), value = TRUE))]
    }
    density_data <- rbind(density_data, 
                          get.density(valid_data_g2, "g2", g2_type, label.g2, xmin, xmax, length.out, opt))
  }

  group_colors <- if (g2_present) {
    setNames(c(color.g1, color.g2), c(label.g1, label.g2))
  } else {
    setNames(c(color.g1), c(label.g1))
  }
  
  x <- density_data$x
  y <- density_data$y
  Group <- density_data$Group
  Study <- density_data$Study

  plot <- ggplot(density_data, aes(x = x, y = y, color = Group, group = interaction(Group, Study))) +
    geom_line(size=line.size) +
    scale_color_manual(values = group_colors) +
    labs(
      title = title,
      x = xlab,
      y = ylab,
      color = "Group"
    ) +
    theme_minimal() +
    theme(
      legend.position = ifelse(display.legend, "right", "none"),
      legend.title = element_text(size = lab.size),
      legend.text = element_text(size = lab.size*0.9),
      plot.title = element_text(size = title.size, hjust = 0.5),  
      axis.title.x = element_text(size = lab.size), 
      axis.title.y = element_text(size = lab.size)
    )
  
  if (display.index) {
    label_data <- density_data %>%
      group_by(Study, Group) %>%
      slice_max(order_by = y, n = 1, with_ties = FALSE) %>%
      ungroup()
    
    plot <- plot + 
      geom_text(
        data = label_data,
        aes(x = x, y = y, label = Study),
        hjust = -0.2, vjust = -0.5,
        size = lab.size/3,
        inherit.aes = FALSE
      )
  }
  return(plot)
}



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
  invalid_rows_g1 <- NULL
  invalid_rows_g2 <- NULL
  quantile_cols_g1 <- grep("\\.g1$", colnames(df), value = TRUE)
  quantile_cols_g2 <- grep("\\.g2$", colnames(df), value = TRUE)
  
  if (length(quantile_cols_g1) > 0) {
    invalid_rows_g1 <- which(rowSums(sapply(df[quantile_cols_g1], function(col) {
      !is.finite(as.numeric(col))
    })) > 0)
  }
  if (length(quantile_cols_g2) > 0) {
    invalid_rows_g2 <- which(rowSums(sapply(df[quantile_cols_g2], function(col) {
      !is.finite(as.numeric(col))
    })) > 0)
  }
  
  if (length(invalid_rows_g1) > 0) {
    warning(paste("Rows", paste(invalid_rows_g1, collapse = ", "), 
                  "contain invalid data for Group 1 and will be skipped from density visualisation."))
  }
  if (length(invalid_rows_g2) > 0) {
    warning(paste("Rows", paste(invalid_rows_g2, collapse = ", "), 
                  "contain invalid data for Group 2 and will be skipped from density visualisation."))
  }
  
  return(list(df = df, invalid_rows_g1 = invalid_rows_g1, invalid_rows_g2 = invalid_rows_g2))
}




get.scenario <- function(data, prefix) {
  if (all(c(paste0("min.", prefix), paste0("q1.", prefix), paste0("med.", prefix), 
            paste0("q3.", prefix), paste0("max.", prefix)) %in% colnames(data))) {
    return("5-number")
  } else if (all(c(paste0("min.", prefix), paste0("med.", prefix), paste0("max.", prefix)) %in% colnames(data))) {
    return("min-med-max")
  } else if (all(c(paste0("q1.", prefix), paste0("med.", prefix), paste0("q3.", prefix)) %in% colnames(data))) {
    return("q1-med-q3")
  } else {
    stop(paste("Invalid or incomplete data for group", prefix))
  }
}

 
get.density <- function(data, prefix, summary_type, group_label, xmin, xmax, length.out, opt) {
  group_data <- data.frame()
  
  for (i in seq_len(nrow(data))) {
    n <- data[[paste0("n.", prefix)]][i]
    
    if (summary_type == "5-number") {
      q <- as.numeric(data[i, c(paste0("min.", prefix), paste0("q1.", prefix), paste0("med.", prefix), 
                                paste0("q3.", prefix), paste0("max.", prefix))])
      dens_out <- est.density.five(min = q[1], q1 = q[2], med = q[3], q3 = q[4], max = q[5], n = n, opt = opt)
      x_vals <- seq(xmin, xmax, length.out = length.out)
      y_vals <- dgl(x_vals, lambda1 = dens_out$parameters[1], lambda2 = dens_out$parameters[2], 
                    lambda3 = dens_out$parameters[3], lambda4 = dens_out$parameters[4])
      
    } else if (summary_type == "min-med-max") {
      q <- as.numeric(data[i, c(paste0("min.", prefix), paste0("med.", prefix), paste0("max.", prefix))])
      dens_out <- est.density.three1(min = q[1], med = q[2], max = q[3], n = n, opt = opt)
      x_vals <- seq(xmin, xmax, length.out = length.out)
      y_vals <- dsl(x_vals, dens_out$parameters)
      
    } else if (summary_type == "q1-med-q3") {
      q <- as.numeric(data[i, c(paste0("q1.", prefix), paste0("med.", prefix), paste0("q3.", prefix))])
      dens_out <- est.density.three2(q1 = q[1], med = q[2], q3 = q[3], n = n, opt = opt)
      x_vals <- seq(xmin, xmax, length.out = length.out)
      y_vals <- dsl(x_vals, dens_out$parameters)
    }
    
    group_density <- data.frame(x = x_vals, y = y_vals, Study = data$study.index[i], Group = group_label)
    group_data <- rbind(group_data, group_density)
  }
  
  return(group_data)
}

