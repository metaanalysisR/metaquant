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
#'   \item Generalised Lambda Distribution (GLD) when 5-number summaries present (\eqn{S_3}).
#'   \item Skew Logistic Distribution (SLD) when 3-number summaries present (\eqn{S_1} and \eqn{S_2}).
#' } 
#' 
#' 
#' @usage plotdist(
#'    data, 
#'    xmin = NULL, 
#'    xmax = NULL, 
#'    ymax = NULL,
#'    length.out = 1000, 
#'    title = "", 
#'    xlab = "x", 
#'    ylab = "Density",
#'    line.size = 0.5,
#'    title.size = 12,
#'    lab.size = 10,
#'    color.g1 = "pink",
#'    color.g2 = "skyblue",
#'    color.g1.pooled = "red",
#'    color.g2.pooled = "blue",
#'    label.g1 = NULL, 
#'    label.g2 = NULL,
#'    display.index = FALSE,
#'    display.legend = FALSE,
#'    pooled.dist = FALSE, 
#'    pooled.only = FALSE, 
#'    opt = TRUE
#' )
#'
#' @param data data frame containing the quantile summary data. For one-group studies, the input may contain the following columns depending on the quantile scenario:
#' \describe{
#'   \item{\code{'stduy.index'}}{stduy index or name}
#'   \item{\code{'min.g1'}}{minimum value}
#'   \item{\code{'q1.g1'}}{first quartile}
#'   \item{\code{'med.g1'}}{median}
#'   \item{\code{'q3.g1'}}{third quartile}
#'   \item{\code{'max.g1'}}{maximum value}
#'   \item{\code{'n.g1'}}{sample size}
#' } For two-group studies, the data frame may also contain the following columns for the second group: \code{min.g2, q1.g2, med.g2, q3.g2, max.g2} and \code{ n.g2}. 
#'    Note that, for three-point summaries (\eqn{S_1} and \eqn{S_2}), only the relevant columns should be included.
#' @param xmin numeric value for the lower limit of the x-axis for density calculation. It is recommended to set this to a value smaller than the smallest value across the quantile summaries to ensure the density curve is fully captured. 
#'    If \code{xmin} is not provided, the minimum value of the 'min.' columns will be used for scenario \eqn{S_1} or \eqn{S_3}. Note that for scenario \eqn{S_2}, no default calculation is performed for \code{xmin}. 
#' @param xmax numeric value for the upper limit of the x-axis for density calculation. It is recommended to set this to a value larger than the largest value across the quantile summaries to ensure the density curve is fully captured.
#'    If \code{xmax} is not provided, the maximum value of the 'max.' columns will be used for scenario \eqn{S_1} or \eqn{S_3}. Similarly, for scenario \eqn{S_2}, no default calculation is performed for \code{xmax}. 
#' @param ymax numeric value for the upper limit of the y-axis. If NULL, the highest density value will be used.
#' @param length.out integer specifying the number of points along the x-axis for density calculation. Default is \code{1000}.
#' @param title character string for the plot title. Default is an empty string.
#' @param xlab character string for the x-axis label. Default is \code{"x"}. 
#' @param ylab character string for the y-axis label. Default is \code{"Density"}.
#' @param line.size numeric. Thickness of the density curve lines. Default is \code{0.5}.
#' @param title.size numeric. Font size for the plot title. Default is \code{12}.
#' @param lab.size numeric. Font size for axis labels. Default is \code{10}.
#' @param color.g1 character string specifying the color for individual density curves of group 1 for each study (row). Default is \code{"pink"}.
#' @param color.g2 character string specifying the color for individual density curves of group 2 for each study (row). Default is \code{"skyblue"}.
#' @param color.g1.pooled character string specifying the color for pooled density curve of group 1. Default is \code{"red"}.
#' @param color.g2.pooled character string specifying the color for pooled density curve of group 2. Default is \code{"blue"}.
#' @param label.g1 character string indicating label or name for group 1 (eg., 'Treatment')
#' @param label.g2 character string indicating label or name for group 2 (eg., 'Control'). 
#'   
#'    If \code{'label.g1'} and \code{'label.g2'} are not provided, the function will assign labels as 'Group 1' and 'Group 2'.
#' @param display.index logical. If \code{TRUE}, the \code{'study.index'} of each quantile set (row) will be displayed alongside the corresponding density curve. The default is \code{FALSE}, meaning no labels will be shown. The label text size is controlled by the \code{lab.size} parameter.
#' @param display.legend logical. If \code{TRUE}, legends (\code{'label.g1'} and/or \code{'label.g2'}) will be displayed on the right side of the plot. The default is \code{FALSE}. The legend text size is controlled by the \code{lab.size} parameter.
#' @param pooled.dist logical. If \code{TRUE}, pooled density curves for group 1 and/or group 2 will be plotted along with the individual density curves. The default is \code{FALSE}. 
#' @param pooled.only logical. If \code{TRUE}, only the pooled density curves of group 1 and/or group 2 will be plotted, excluding the individual density curves. The default is \code{FALSE}. 
#' @param opt logical value indicating whether to apply the optimization step when estimating GLD or SLD parameters. The default value is \code{TRUE}.
#'   
#' @details
#' The generalised lambda distribution (GLD) is a four parameter family of distributions defined by its quantile function under the FKML parameterisation (Freimer et al., 1988).
#' De Livera et al. propose that the GLD quantile function can be used to approximate a sample's distribution using 5-point summaries. 
#' The four parameters of GLD quantile function include: a location parameter (\eqn{\lambda_1}), an inverse scale parameter (\eqn{\lambda_2}>0), and two shape parameters (\eqn{\lambda_3} and \eqn{\lambda_4}).
#' 
#' The quantile-based skew logistic distribution (SLD), introduced by Gilchrist (2000) and further modified by van Staden and King (2015) 
#' is used to approximate the sample's distribution using 3-point summaries.
#' The SLD quantile function is defined using three parameters: a location parameter (\eqn{\lambda}), a scale parameter (\eqn{\eta}), and a skewing parameter (\eqn{\delta}).
#' 
#' These parameters of GLD and SLD are estimated by formulating and solving a series of simultaneous equations which relate the estimated quantiles
#' with the population counterparts of respective distribution (GLD or SLD). The \code{plotdist} uses these estimated parameters, to compute the density data 
#' using \code{\link[gld]{dgl}} function from the \pkg{gld} package and \code{\link[sld]{dsl}} function from the \pkg{sld} package.
#' 
#' If one needs to generate pooled density plots, they can use the \code{pooled.dist} or \code{pooled.only} arguments as described in the *Arguments* section. 
#' The pooled density curves represent a weighted average of individual study densities, with weights determined by sample sizes. The method is similar to obtaining pooled 
#' estimates of effects in a standard meta-analysis and it serves as a way to visualize combined estimated distributional information across studies.
#'
#' @return An interactive plotly object visualizing the estimated density curve(s) for one or two groups.
#'
#' @examples 
#' #Example dataset of 3-point summaries (min, med, max) for 2 groups
#' data_3num_2g <- data.frame(
#'   study.index = c("Study 1", "Study 2", "Study 3"),
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
#' #Density plots of two groups along with the pooled plots
#' plot_2g <- plotdist(
#'   data_3num_2g,
#'   xmin = 10,
#'   xmax = 125,
#'   title = "Example Density Plots of Two Groups",
#'   xlab = "x data",
#'   color.g1 = "skyblue",
#'   color.g2 = "pink",
#'   color.g1.pooled = "blue",
#'   color.g2.pooled = "red",
#'   label.g1 = "Treatment", 
#'   label.g2 = "Control",
#'   display.legend = TRUE,
#'   pooled.dist = TRUE
#' )
#' print(plot_2g)
#'
#' @references De Livera, A. M., Prendergast, L., & Kumaranathunga, U. (2024). A novel density-based approach for estimating unknown means, distribution visualisations and meta-analyses of quantiles. \emph{arXiv preprint arXiv:2411.10971}. <https://arxiv.org/abs/2411.10971>.
#' @references Freimer, M., Kollia, G., Mudholkar, G. S., & Lin, C. T. (1988). A study of the generalized Tukey lambda family. \emph{Communications in Statistics—Theory and Methods, 17}(10), 3547–3567.
#' @references Gilchrist, W. (2000). \emph{Statistical modelling with quantile functions}. Chapman & Hall/CRC.
#' @references van Staden, P. J., & King, R. A. R. (2015). The quantile-based skew logistic distribution. \emph{Statistics & Probability Letters, 96}, 109–116.
#' @export 
#' 
#' @importFrom gld dgl
#' @importFrom sld dsl
#' @importFrom ggplot2 ggplot aes geom_line scale_color_manual labs theme_minimal theme element_text geom_text xlim ylim
#' @importFrom plotly ggplotly
#' @importFrom stats setNames
#' @importFrom magrittr %>%
#' @importFrom dplyr group_by slice_max ungroup mutate


plotdist <- function(data, 
                     xmin = NULL, 
                     xmax = NULL, 
                     ymax = NULL,
                     length.out = 1000, 
                     title = "", 
                     xlab = "x", 
                     ylab = "Density",
                     line.size = 0.5,
                     title.size = 12,
                     lab.size = 10,
                     color.g1 = "pink",
                     color.g2 = "skyblue",
                     color.g1.pooled = "red",
                     color.g2.pooled = "blue",
                     label.g1 = NULL, 
                     label.g2 = NULL,
                     display.index = FALSE,
                     display.legend = FALSE,
                     pooled.dist = FALSE,
                     pooled.only = FALSE,
                     opt = TRUE) {
  
  check_result <-check.df(data)
  data <-check_result$df
  invalid_rows_g1 <-check_result$invalid_rows_g1
  invalid_rows_g2 <-check_result$invalid_rows_g2
  
  g1_type <-get.scenario(data, "g1")
  g2_present <-any(grepl("\\.g2$", colnames(data)))
  g2_type <- if(g2_present) get.scenario(data, "g2") else NULL
  
  if (is.null(label.g1)) label.g1 <-"Group 1"
  if (g2_present && is.null(label.g2)) label.g2 <-"Group 2"
  
  density_data <- data.frame()
  pooled_data <- data.frame()
  
  valid_data_g1 <- if (length(invalid_rows_g1) > 0) {
    data[-invalid_rows_g1, c("study.index", grep("\\.g1$", colnames(data), value = TRUE))]
  } else {
    data[, c("study.index", grep("\\.g1$", colnames(data), value = TRUE))]
  }
  
  if (pooled.only || pooled.dist) {
    pooled_data <- rbind(
      pooled_data,
      get.pooled.density(valid_data_g1, "g1", g1_type, label.g1, xmin, xmax, length.out, opt)
    )
  }
  
  if (!pooled.only) {
    density_data <- rbind(density_data, 
                          get.density(valid_data_g1, "g1", g1_type, label.g1, xmin, xmax, length.out, opt))
  }
  
  if (g2_present) {
    valid_data_g2 <- if (length(invalid_rows_g2) > 0) {
      data[-invalid_rows_g2, c("study.index", grep("\\.g2$", colnames(data), value = TRUE))]
    } else {
      data[, c("study.index", grep("\\.g2$", colnames(data), value = TRUE))]
    }
    
    if (pooled.only || pooled.dist) {
      pooled_data <- rbind(
        pooled_data,
        get.pooled.density(valid_data_g2, "g2", g2_type, label.g2, xmin, xmax, length.out, opt)
      )
    }
    
    if (!pooled.only) {
      density_data <- rbind(density_data, 
                            get.density(valid_data_g2, "g2", g2_type, label.g2, xmin, xmax, length.out, opt))
    }
  }
  
  group_colors <- setNames(c(color.g1, color.g2), c(label.g1, label.g2))
  if (pooled.dist || pooled.only) {
    group_colors <- c(group_colors, setNames(c(color.g1.pooled, color.g2.pooled), 
                                             c(paste0(label.g1, " (Pooled)"), paste0(label.g2, " (Pooled)"))))
  }
  
  plot_data <- if (pooled.only) pooled_data else rbind(density_data, pooled_data)
  plot_data$Group <- factor(plot_data$Group, levels = names(group_colors))
  plot_data <- plot_data %>% mutate(Interaction = paste(Group, Study, sep = ", "))
  
  x <-plot_data$x
  y <-plot_data$y
  Group <-plot_data$Group
  Study <-plot_data$Study
  Interaction <-plot_data$Interaction
  
  plot <- ggplot(plot_data, aes(x = x, y = y, color = Group, group = Interaction)) +
    geom_line(size = line.size) +
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
      legend.title = element_text(size=lab.size),
      legend.text = element_text(size=lab.size*0.9),
      plot.title = element_text(size=title.size, hjust=0.5),  
      axis.title.x = element_text(size=lab.size), 
      axis.title.y = element_text(size=lab.size)
    )+
    xlim(if (is.null(xmin)) min(x) else xmin, if (is.null(xmax)) max(plot_data$x) else xmax) +
    ylim(0, if (is.null(ymax)) max(plot_data$y) else ymax)
  
  if (display.index && !pooled.only) {
    label_data <- plot_data %>% 
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
  
  interact_plot <- ggplotly(plot, tooltip = c("x", "y", "Interaction"))
  return(interact_plot)
}


