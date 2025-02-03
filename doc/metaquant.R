## -----------------------------------------------------------------------------
#Using CRAN
#install.packages("metaquant")
library(metaquant)

## -----------------------------------------------------------------------------
#install.packages("devtools")
#library(devtools)
#devtools::install_github("metaanalysisR/metaquant")

## -----------------------------------------------------------------------------
# Load the libraries
library(metaquant)
library(stats)

## -----------------------------------------------------------------------------
#Generate quantile summary data
set.seed(123)
n <- 100
x <- rlnorm(n, 4, 0.3)
quants <- c(min(x), quantile(x, probs = c(0.25, 0.5, 0.75)), max(x))
quants

## -----------------------------------------------------------------------------
#Estimate sample mean of S3 using 'gld/sld'
estmean_gl <- est.mean(min = quants[1], 
                       q1 = quants[2], 
                       med = quants[3],
                       q3 = quants[4],
                       max = quants[5],
                       n=n)
estmean_gl

## -----------------------------------------------------------------------------
#Estimate sample mean of S3 using the method 'luo'
estmean_luo <- est.mean(min = quants[1], 
                        q1 = quants[2], 
                        med = quants[3],
                        q3 = quants[4],
                        max = quants[5],
                        n=n,
                        method = "luo")
estmean_luo

## -----------------------------------------------------------------------------
# 3-point summary data for S1
quants1 <- c(min(x), quantile(x, probs = 0.5), max(x))
quants1

## -----------------------------------------------------------------------------
#Estimate sample mean for S1
estmean_sl_1 <- est.mean(min = quants1[1], 
                        med = quants1[2],
                        max = quants1[3],
                        n=n,
                        method = "gld/sld")
estmean_sl_1

## -----------------------------------------------------------------------------
# 3-point summary data for S2
quants2 <- quantile(x, probs = c(0.25, 0.5, 0.75))
quants2

## -----------------------------------------------------------------------------
#Estimate sample mean for S2
estmean_sl_2 <- est.mean(q1 = quants2[1], 
                        med = quants2[2],
                        q3 = quants2[3],
                        method = "gld/sld")
estmean_sl_2

## -----------------------------------------------------------------------------
#Estimate sample SD of S3 using 'shi/wan' method
estsd_shi <- est.sd(min = quants[1], 
                    q1 = quants[2], 
                    med = quants[3],
                    q3 = quants[4],
                    max = quants[5],
                    n=n)
estsd_shi

## -----------------------------------------------------------------------------
#Generate 5-point summary data for two groups
set.seed(123)
n_t <- 100
n_c <- 120
x_t <- rexp(n_t, 5)
x_c <- rexp(n_c, 10)
q_t <- c(min(x_t), quantile(x_t, probs = c(0.25, 0.5, 0.75)), max(x_t))
q_c <- c(min(x_c), quantile(x_c, probs = c(0.25, 0.5, 0.75)), max(x_c))


## -----------------------------------------------------------------------------
#Estimate sample mean of S3 
estmean_2g <- est.mean.2g(q_t[1],q_t[2],q_t[3],q_t[4],q_t[5],
                          q_c[1],q_c[2],q_c[3],q_c[4],q_c[5],
                          n.g1 = n_t,
                          n.g2 = n_c)
estmean_2g

## -----------------------------------------------------------------------------
#Estimate sample SD of S3 
estsd_2g <- est.sd.2g(q_t[1],q_t[2],q_t[3],q_t[4],q_t[5],
                      q_c[1],q_c[2],q_c[3],q_c[4],q_c[5],
                      n.g1 = n_t,
                      n.g2 = n_c)
estsd_2g

## -----------------------------------------------------------------------------
# Dataset of 5-point summaries for 1 group
data_s3 <- data.frame(
  study.index = c("Study1", "Study2", "Study3"),
  min.g1 = c(18, 19, 15),
  q1.g1 = c(66, 71, 69),
  med.g1 = c(73, 82, 81),
  q3.g1 = c(80, 93, 89),
  max.g1 = c(110, 115, 100),
  n.g1 = c(226, 230, 200)
)
data_s3

## -----------------------------------------------------------------------------
# Plot densities 
plot_s3 <- plotdist(
  data_s3,
  xmin = 10,
  xmax = 125,
  title = "Example Density Plot of S3",
  xlab = "x data",
  title.size = 11,
  lab.size = 10,
  color.g1 = "blue",
  display.index = FALSE,
  display.legend = FALSE
)
plot_s3

## -----------------------------------------------------------------------------
# Dataset of 3-point summaries for 1 group
data_s1 <- data.frame(
  study.index = c("Study1", "Study2", "Study3"),
  min.g1 = c(18, 19, 15),
  med.g1 = c(73, 82, 81),
  max.g1 = c(110, 115, 100),
  n.g1 = c(226, 230, 200)
)
data_s1

## -----------------------------------------------------------------------------
# Plot densities 
plot_s1 <- plotdist(
  data_s1,
  xmin = 10,
  xmax = 125,
  title = "Example Density Plot of S1",
  xlab = "x data",
  color.g1 = "purple",
  display.index = FALSE,
  display.legend = FALSE
)
plot_s1

## -----------------------------------------------------------------------------
# Dataset of 5-point summaries for 2 groups
data_2g <- data.frame(
  study.index = c("Study1", "Study2", "Study3"),
  min.g1 = c(18, 19, 15),
  q1.g1 = c(66, 71, 69),
  med.g1 = c(73, 82, 81),
  q3.g1 = c(80, 93, 89),
  max.g1 = c(110, 115, 100),
  n.g1 = c(226, 230, 200),
  min.g2 = c(15, 15, 13),
  q1.g2 = c(57, 59, 55),
  med.g2 = c(66, 68, 60),
  q3.g2 = c(74, 72, 69),
  max.g2 = c(108, 101, 100),
  n.g2 = c(201, 223, 198)
)
data_2g

## -----------------------------------------------------------------------------
# Plot densities 
plot_2g <- plotdist(
  data_2g,
  xmin = 10,
  xmax = 125,
  title = "Example Density Plot of Two Groups",
  xlab = "x data",
  title.size = 11,
  label.g1 = "Treatment", 
  label.g2 = "Control",
  display.index = FALSE,
  display.legend = TRUE
)
plot_2g

## -----------------------------------------------------------------------------
# Plot densities with pooled curves
plot_2g <- plotdist(
  data_2g,
  xmin = 10,
  xmax = 125,
  title = "Example Density Plot with Pooled Densities",
  xlab = "x data",
  title.size = 11,
  label.g1 = "Treatment", 
  label.g2 = "Control",
  display.index = FALSE,
  display.legend = FALSE,
  pooled.dist = TRUE
)
plot_2g

## -----------------------------------------------------------------------------
sessionInfo()

