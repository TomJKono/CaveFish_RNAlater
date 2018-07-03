# This script will do the linear model fitting of the RNAlater data. It will
# read in the merged predictor data and plot some summaries of the variables,
# fit the linear model, and test for goodness of fit.

# Read in the covariate data
covariates <- read.csv("/Users/tomkono/Dropbox/GitHub/RIS/SEM_RNALater/Data/RNALater_Predictors.csv", header=TRUE)

# And then the l2fc data
l2fc <- read.csv("/Users/tomkono/Dropbox/GitHub/RIS/SEM_RNALater/Results/RNAlater_v_LN_Pvalues.csv", header=TRUE)

# Merge the two tables on gene ID
merged <- merge(covariates, l2fc, by.x="gene_id", by.y="ENSID")

# Trim down the merged file to just the predictors we care about
merged <- merged[,c("baseMean", "gene_id", "gene_length", "GC_prop", "exon_number", "SSR", "HPR", "log2FoldChange")]

# Convert SSR and HPR into a factor
merged$SSR <- factor(merged$SSR, levels=c("0", "1"))
merged$HPR <- factor(merged$HPR, levels=c("0", "1"))

# Then, fit a model that has everything:
#   Mean expression
#   Length
#   GC content
#   N. exons
#   SSR pres/abs
#   Homopolymer repeat pres/abs
#   GC:SSR
#   GC:HPR
full_model <- lm(
    log2FoldChange ~ baseMean + GC_prop + gene_length + exon_number + SSR + HPR + GC_prop:SSR + GC_prop:HPR,
    data=merged)
summary(full_model)

# Call:
# lm(formula = log2FoldChange ~ baseMean + GC_prop + gene_length +
#     exon_number + SSR + HPR + GC_prop:SSR + GC_prop:HPR, data = merged)

# Residuals:
#      Min       1Q   Median       3Q      Max
# -25.6429  -0.3399   0.1743   0.5788  26.4914

# Coefficients:
#                Estimate Std. Error t value Pr(>|t|)
# (Intercept)  -2.467e-01  4.976e-01  -0.496   0.6201
# baseMean      3.857e-06  1.643e-06   2.348   0.0189 *
# GC_prop      -4.528e-01  9.505e-01  -0.476   0.6338
# gene_length  -1.204e-05  1.284e-05  -0.938   0.3484
# exon_number   2.097e-02  2.142e-03   9.789   <2e-16 ***
# SSR1          1.273e-01  5.130e-01   0.248   0.8040
# HPR1          4.902e-01  2.983e-01   1.643   0.1003
# GC_prop:SSR1 -1.806e-01  9.797e-01  -0.184   0.8537
# GC_prop:HPR1 -8.449e-01  5.691e-01  -1.485   0.1376
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Residual standard error: 1.556 on 14602 degrees of freedom
#   (533 observations deleted due to missingness)
# Multiple R-squared:  0.017, Adjusted R-squared:  0.01646
# F-statistic: 31.57 on 8 and 14602 DF,  p-value: < 2.2e-16

# It looks like the interaction terms aren't significant. Remove them and then
# test the additional variance explained
red1_mod <- lm(
    log2FoldChange ~ baseMean + gene_length + GC_prop + exon_number + SSR + HPR,
    data=merged)
summary(red1_mod)

# Call:
# lm(formula = log2FoldChange ~ baseMean + gene_length + GC_prop +
#     exon_number + SSR + HPR, data = merged)

# Residuals:
#      Min       1Q   Median       3Q      Max
# -25.6465  -0.3426   0.1744   0.5785  26.4995

# Coefficients:
#               Estimate Std. Error t value Pr(>|t|)
# (Intercept)  7.689e-02  1.562e-01   0.492 0.622471
# baseMean     3.855e-06  1.643e-06   2.347 0.018944 *
# gene_length -1.182e-05  1.284e-05  -0.921 0.357201
# GC_prop     -1.075e+00  2.844e-01  -3.781 0.000157 ***
# exon_number  2.094e-02  2.142e-03   9.777  < 2e-16 ***
# SSR1         3.400e-02  4.466e-02   0.761 0.446487
# HPR1         4.912e-02  2.611e-02   1.881 0.059984 .
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Residual standard error: 1.556 on 14604 degrees of freedom
#   (533 observations deleted due to missingness)
# Multiple R-squared:  0.01685,   Adjusted R-squared:  0.01644
# F-statistic:  41.7 on 6 and 14604 DF,  p-value: < 2.2e-16

# How do the models compare?
anova(red1_mod, full_model, test="LRT")
# Analysis of Variance Table

# Model 1: log2FoldChange ~ baseMean + gene_length + GC_prop + exon_number +
#     SSR + HPR
# Model 2: log2FoldChange ~ baseMean + GC_prop + gene_length + exon_number +
#     SSR + HPR + GC_prop:SSR + GC_prop:HPR
#   Res.Df   RSS Df Sum of Sq Pr(>Chi)
# 1  14604 35377
# 2  14602 35371  2    5.6537   0.3113

# The SSR term doesn't help us, I think. Let's drop it.
red2_mod <- lm(
    log2FoldChange ~ baseMean + gene_length + GC_prop + exon_number + HPR,
    data=merged)
summary(red2_mod)
# Call:
# lm(formula = log2FoldChange ~ baseMean + gene_length + GC_prop +
#     exon_number + HPR, data = merged)

# Residuals:
#      Min       1Q   Median       3Q      Max
# -25.6447  -0.3413   0.1747   0.5782  26.5043

# Coefficients:
#               Estimate Std. Error t value Pr(>|t|)
# (Intercept)  1.058e-01  1.515e-01   0.698 0.485106
# baseMean     3.851e-06  1.642e-06   2.345 0.019056 *
# gene_length -1.182e-05  1.284e-05  -0.921 0.357118
# GC_prop     -1.074e+00  2.844e-01  -3.778 0.000158 ***
# exon_number  2.095e-02  2.142e-03   9.778  < 2e-16 ***
# HPR1         5.209e-02  2.582e-02   2.017 0.043673 *
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Residual standard error: 1.556 on 14605 degrees of freedom
#   (533 observations deleted due to missingness)
# Multiple R-squared:  0.01681,   Adjusted R-squared:  0.01647
# F-statistic: 49.93 on 5 and 14605 DF,  p-value: < 2.2e-16

# It does not improve model fit
anova(red1_mod, red2_mod, test="LRT")
# Analysis of Variance Table

# Model 1: log2FoldChange ~ baseMean + gene_length + GC_prop + exon_number +
#     SSR + HPR
# Model 2: log2FoldChange ~ baseMean + gene_length + GC_prop + exon_number +
#     HPR
#   Res.Df   RSS Df Sum of Sq Pr(>Chi)
# 1  14604 35377
# 2  14605 35378 -1    -1.404   0.4465

# What about length?
red3_mod <- lm(
    log2FoldChange ~ baseMean + GC_prop + exon_number + HPR,
    data=merged)
summary(red3_mod)
# Call:
# lm(formula = log2FoldChange ~ baseMean + GC_prop + exon_number +
#     HPR, data = merged)

# Residuals:
#      Min       1Q   Median       3Q      Max
# -25.6417  -0.3403   0.1757   0.5789  26.5112

# Coefficients:
#               Estimate Std. Error t value Pr(>|t|)
# (Intercept)  1.114e-01  1.514e-01   0.736 0.461736
# baseMean     3.893e-06  1.642e-06   2.371 0.017761 *
# GC_prop     -1.092e+00  2.837e-01  -3.850 0.000119 ***
# exon_number  1.941e-02  1.340e-03  14.481  < 2e-16 ***
# HPR1         5.196e-02  2.582e-02   2.012 0.044202 *

anova(red3_mod)

anova(red2_mod, red3_mod, test="LRT")
# Analysis of Deviance Table
# Model 1: log2FoldChange ~ gene_length + GC_prop + exon_number + HPR
# Model 2: log2FoldChange ~ GC_prop + exon_number + HPR
#   Resid. Df Resid. Dev Df Deviance Pr(>Chi)
# 1     14606      35391
# 2     14607      35394 -1  -2.3526   0.3245

anova(red3_mod, full_model, test="LRT")
