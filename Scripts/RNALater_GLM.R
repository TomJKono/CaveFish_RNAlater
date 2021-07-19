# This script will do the linear model fitting of the RNAlater data. It will
# read in the merged predictor data and plot some summaries of the variables,
# fit the linear model, and test for goodness of fit.

# Read in the covariate data
covariates <- read.csv("/Users/tomkono/Dropbox/GitHub/RIS/SEM_RNALater/Data/RNALater_Predictors.csv", header=TRUE)

# And then the l2fc data
l2fc <- read.csv("/Users/tomkono/Dropbox/GitHub/RIS/SEM_RNALater/Results/RNAlater_v_LN_Pvalues.csv", header=TRUE)

# Merge the two tables on gene ID
merged <- merge(covariates, l2fc, by.x="GeneID", by.y="ENSID")

# Trim down the merged file to just the predictors we care about
merged <- merged[,c("baseMean", "GeneID", "Longest.TxLen", "Longest.GCProp", "Longest.NExon", "SSR", "HPR", "log2FoldChange")]

merged$baseMean <- log2(merged$baseMean)

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
    log2FoldChange ~ baseMean + Longest.GCProp + Longest.TxLen + Longest.NExon + SSR + HPR + Longest.GCProp:SSR + Longest.GCProp:HPR,
    data=merged)
summary(full_model)

# Call:
# lm(formula = log2FoldChange ~ baseMean + Longest.GCProp + Longest.TxLen +
#     Longest.NExon + SSR + HPR + Longest.GCProp:SSR + Longest.GCProp:HPR,
#     data = merged)

# Residuals:
#      Min       1Q   Median       3Q      Max
# -25.4940  -0.4429   0.1086   0.5925  26.3636

# Coefficients:
#                       Estimate Std. Error t value Pr(>|t|)
# (Intercept)          9.269e-01  7.046e-01   1.316 0.188371
# baseMean             1.552e-01  7.329e-03  21.175  < 2e-16 ***
# Longest.GCProp      -5.068e+00  1.380e+00  -3.674 0.000241 ***
# Longest.TxLen       -4.264e-06  1.074e-05  -0.397 0.691268
# Longest.NExon        2.768e-02  2.307e-03  11.995  < 2e-16 ***
# SSR1                -1.706e+00  7.148e-01  -2.387 0.017006 *
# HPR1                 2.988e-01  3.128e-01   0.955 0.339504
# Longest.GCProp:SSR1  3.440e+00  1.406e+00   2.447 0.014423 *
# Longest.GCProp:HPR1 -6.069e-01  6.233e-01  -0.974 0.330248
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Residual standard error: 1.481 on 8772 degrees of freedom
# Multiple R-squared:  0.08651,   Adjusted R-squared:  0.08568
# F-statistic: 103.8 on 8 and 8772 DF,  p-value: < 2.2e-16


# It looks like the GC:HPR interaction term isn't significant. Remove it then
# test the additional variance explained
red1_mod <- lm(
    log2FoldChange ~ baseMean + Longest.TxLen + Longest.GCProp + Longest.NExon + SSR + HPR + Longest.GCProp:SSR,
    data=merged)
summary(red1_mod)

# Call:
# lm(formula = log2FoldChange ~ baseMean + Longest.TxLen + Longest.GCProp +
#     Longest.NExon + SSR + HPR + Longest.GCProp:SSR, data = merged)

# Residuals:
#      Min       1Q   Median       3Q      Max
# -25.4836  -0.4444   0.1067   0.5902  26.3761

# Coefficients:
#                       Estimate Std. Error t value Pr(>|t|)
# (Intercept)          1.040e+00  6.950e-01   1.496   0.1348
# baseMean             1.554e-01  7.327e-03  21.203  < 2e-16 ***
# Longest.TxLen       -3.671e-06  1.072e-05  -0.342   0.7320
# Longest.GCProp      -5.295e+00  1.360e+00  -3.893 9.97e-05 ***
# Longest.NExon        2.739e-02  2.289e-03  11.968  < 2e-16 ***
# SSR1                -1.600e+00  7.063e-01  -2.265   0.0236 *
# HPR1                -3.540e-03  3.783e-02  -0.094   0.9254
# Longest.GCProp:SSR1  3.235e+00  1.390e+00   2.327   0.0200 *
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Residual standard error: 1.481 on 8773 degrees of freedom
# Multiple R-squared:  0.08641,   Adjusted R-squared:  0.08568
# F-statistic: 118.5 on 7 and 8773 DF,  p-value: < 2.2e-16

# How do the models compare?
anova(red1_mod, full_model, test="LRT")

# Analysis of Variance Table

# Model 1: log2FoldChange ~ baseMean + Longest.TxLen + Longest.GCProp +
#     Longest.NExon + SSR + HPR + Longest.GCProp:SSR
# Model 2: log2FoldChange ~ baseMean + Longest.GCProp + Longest.TxLen +
#     Longest.NExon + SSR + HPR + Longest.GCProp:SSR + Longest.GCProp:HPR
#   Res.Df   RSS Df Sum of Sq Pr(>Chi)
# 1   8773 19251
# 2   8772 19249  1    2.0803   0.3302

# Next term to drop: HPR
red2_mod <- lm(
    log2FoldChange ~ baseMean + Longest.TxLen + Longest.GCProp + Longest.NExon + SSR + Longest.GCProp:SSR,
    data=merged)
summary(red2_mod)

# Call:
# lm(formula = log2FoldChange ~ baseMean + Longest.TxLen + Longest.GCProp +
#     Longest.NExon + SSR + Longest.GCProp:SSR, data = merged)

# Residuals:
#      Min       1Q   Median       3Q      Max
# -25.4844  -0.4448   0.1067   0.5903  26.3780

# Coefficients:
#                       Estimate Std. Error t value Pr(>|t|)
# (Intercept)          1.039e+00  6.949e-01   1.495   0.1350
# baseMean             1.554e-01  7.325e-03  21.212  < 2e-16 ***
# Longest.TxLen       -3.798e-06  1.063e-05  -0.357   0.7210
# Longest.GCProp      -5.295e+00  1.360e+00  -3.894 9.94e-05 ***
# Longest.NExon        2.738e-02  2.286e-03  11.978  < 2e-16 ***
# SSR1                -1.603e+00  7.053e-01  -2.273   0.0231 *
# Longest.GCProp:SSR1  3.240e+00  1.389e+00   2.333   0.0197 *
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Residual standard error: 1.481 on 8774 degrees of freedom
# Multiple R-squared:  0.08641,   Adjusted R-squared:  0.0857

anova(red2_mod, red1_mod, test="LRT")

# Analysis of Variance Table

# Model 1: log2FoldChange ~ baseMean + Longest.TxLen + Longest.GCProp +
#     Longest.NExon + SSR + Longest.GCProp:SSR
# Model 2: log2FoldChange ~ baseMean + Longest.TxLen + Longest.GCProp +
#     Longest.NExon + SSR + HPR + Longest.GCProp:SSR
#   Res.Df   RSS Df Sum of Sq Pr(>Chi)
# 1   8774 19251
# 2   8773 19251  1  0.019215   0.9254


# Let's get rid of length
red3_mod <- lm(
    log2FoldChange ~ baseMean + Longest.GCProp + Longest.NExon + SSR + Longest.GCProp:SSR,
    data=merged)
summary(red3_mod)
# Call:
# lm(formula = log2FoldChange ~ baseMean + Longest.GCProp + Longest.NExon +
#     SSR + Longest.GCProp:SSR, data = merged)

# Residuals:
#      Min       1Q   Median       3Q      Max
# -25.4766  -0.4432   0.1066   0.5908  26.3817

# Coefficients:
#                      Estimate Std. Error t value Pr(>|t|)
# (Intercept)          1.027018   0.694111   1.480 0.139012
# baseMean             0.155547   0.007308  21.283  < 2e-16 ***
# Longest.GCProp      -5.277778   1.358944  -3.884 0.000104 ***
# Longest.NExon        0.026825   0.001670  16.059  < 2e-16 ***
# SSR1                -1.620474   0.703584  -2.303 0.021292 *
# Longest.GCProp:SSR1  3.269607   1.386380   2.358 0.018377 *
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Residual standard error: 1.481 on 8775 degrees of freedom
# Multiple R-squared:  0.0864,    Adjusted R-squared:  0.08588
# F-statistic:   166 on 5 and 8775 DF,  p-value: < 2.2e-16

anova(red3_mod, red2_mod, test="LRT")
# Analysis of Variance Table

# Model 1: log2FoldChange ~ baseMean + Longest.GCProp + Longest.NExon +
#     SSR + Longest.GCProp:SSR
# Model 2: log2FoldChange ~ baseMean + Longest.TxLen + Longest.GCProp +
#     Longest.NExon + SSR + Longest.GCProp:SSR
#   Res.Df   RSS Df Sum of Sq Pr(>Chi)
# 1   8775 19251
# 2   8774 19251  1   0.27985    0.721

anova(red3_mod)


red4_mod <- lm(
    log2FoldChange ~ baseMean + Longest.GCProp + Longest.NExon + SSR,
    data=merged)
summary(red4_mod)
# Call:
# lm(formula = log2FoldChange ~ baseMean + Longest.GCProp + Longest.NExon +
#     SSR, data = merged)

# Residuals:
#      Min       1Q   Median       3Q      Max
# -25.4849  -0.4440   0.1062   0.5920  26.3895

# Coefficients:
#                Estimate Std. Error t value Pr(>|t|)
# (Intercept)    -0.55693    0.17529  -3.177  0.00149 **
# baseMean        0.15538    0.00731  21.256  < 2e-16 ***
# Longest.GCProp -2.13822    0.27308  -7.830 5.45e-15 ***
# Longest.NExon   0.02693    0.00167  16.120  < 2e-16 ***
# SSR1            0.02637    0.08611   0.306  0.75941
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Residual standard error: 1.482 on 8776 degrees of freedom
# Multiple R-squared:  0.08582,   Adjusted R-squared:  0.0854
# F-statistic:   206 on 4 and 8776 DF,  p-value: < 2.2e-16
anova(red4_mod)
# Analysis of Variance Table

# Response: log2FoldChange
#                  Df  Sum Sq Mean Sq  F value    Pr(>F)
# baseMean          1  1088.8 1088.76 496.0140 < 2.2e-16 ***
# Longest.GCProp    1   134.5  134.50  61.2750  5.54e-15 ***
# Longest.NExon     1   584.9  584.94 266.4833 < 2.2e-16 ***
# SSR               1     0.2    0.21   0.0938    0.7594
# Residuals      8776 19263.5    2.20
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

red5_mod <- lm(
    log2FoldChange ~ baseMean + Longest.GCProp + Longest.NExon,
    data=merged)
summary(red5_mod)
# Call:
# lm(formula = log2FoldChange ~ baseMean + Longest.GCProp + Longest.NExon,
#     data = merged)

# Residuals:
#      Min       1Q   Median       3Q      Max
# -25.4840  -0.4436   0.1061   0.5927  26.3914

# Coefficients:
#                 Estimate Std. Error t value Pr(>|t|)
# (Intercept)    -0.528941   0.149587  -3.536 0.000408 ***
# baseMean        0.155294   0.007304  21.262  < 2e-16 ***
# Longest.GCProp -2.144052   0.272405  -7.871 3.94e-15 ***
# Longest.NExon   0.026997   0.001654  16.325  < 2e-16 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Residual standard error: 1.481 on 8777 degrees of freedom
# Multiple R-squared:  0.08581,   Adjusted R-squared:  0.0855
# F-statistic: 274.6 on 3 and 8777 DF,  p-value: < 2.2e-16

anova(red5_mod)
# Analysis of Variance Table

# Response: log2FoldChange
#                  Df  Sum Sq Mean Sq F value    Pr(>F)
# baseMean          1  1088.8 1088.76 496.065 < 2.2e-16 ***
# Longest.GCProp    1   134.5  134.50  61.281 5.522e-15 ***
# Longest.NExon     1   584.9  584.94 266.511 < 2.2e-16 ***
# Residuals      8777 19263.7    2.19
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
