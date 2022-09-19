
library(ggplot2)
library(reshape2)

set.seed(100)
data <- data.frame(x = rnorm(100, 2, 1), y = rnorm(100, 1, 1))
data$id <- 1:nrow(data)
data2 <- melt(data, id.vars = "id")
lapply(data, function(x) get_summary_stats(data.frame(x)))
# $x
# # A tibble: 1 x 13
# variable     n    min   max median    q1    q3   iqr   mad  mean    sd    se    ci
# <chr>    <dbl>  <dbl> <dbl>  <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>
#   1 x          100 -0.272  4.58   1.94  1.39  2.66  1.26 0.974  2.00  1.02 0.102 0.203
# 
# $y
# # A tibble: 1 x 13
# variable     n   min   max median    q1    q3   iqr   mad  mean    sd    se    ci
# <chr>    <dbl> <dbl> <dbl>  <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>
#   1 x          100 -1.14  3.17  0.927 0.568  1.45 0.878 0.648  1.01 0.796  0.08 0.158

data$diff <- data$x - data$y
shapiro.test(data$diff)
# Shapiro-Wilk normality test
# 
# data:  data$diff
# W = 0.98792, p-value = 0.5024

t.test(data$x, data$y, paired = T)
# Paired t-test
# 
# data:  data$x and data$y
# t = 7.2039, df = 99, p-value = 1.165e-10
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   0.7186033 1.2649401
# sample estimates:
#   mean of the differences 
# 0.9917717 

wilcox.test(data$x, data$y, paired = T)
# Wilcoxon signed rank test with continuity correction
# 
# data:  data$x and data$y
# V = 4287, p-value = 1.39e-09
# alternative hypothesis: true location shift is not equal to 0

ggplot(data = data2, aes(x = variable, y  = value)) +
  geom_line(aes(group = id))+
  geom_point(aes(color = variable))