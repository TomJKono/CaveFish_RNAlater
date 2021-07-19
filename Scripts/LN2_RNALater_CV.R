# Generate plots of the coefficient of variation for RNALater vs. LN2 samples.

norm_counts <- read.csv("/Users/tomkono/Dropbox/GitHub/RIS/SEM_RNALater/Results/RNALater_Normalized_Counts.csv", header=TRUE, row.names=1)

# Calculate the coefficient of variation for the LN2 samples
ln_cv <- apply(norm_counts, 1, function(x) {
    m <- mean(as.numeric(x[1:6]))
    s <- sd(as.numeric(x[1:6]))
    return(s/m)
    })

# and the same for the RNAlater
rl_cv <- apply(norm_counts, 1, function(x) {
    m <- mean(as.numeric(x[7:11]))
    s <- sd(as.numeric(x[7:11]))
    return(s/m)
    })

# Do a correlation test, with a kendall correlation to test rank order
cor.test(ln_cv, rl_cv, method="kendall")

# And make a plot
pdf(file="RNALater_LN2_CV_Plot.pdf", height=6, width=6)
plot(
    ln_cv ~ rl_cv,
    pch=19,
    cex=0.25,
    xlab="RNAlater CV",
    ylab="LN2 CV",
    main="Coefficient of Variation Across Storage Conditions")
abline(c(0, 0), c(1, 1), lwd=2, col="red")
dev.off()

comb_cor <- as.data.frame(cbind(ln_cv, rl_cv))
sum(comb_cor$ln_cv < comb_cor$rl_cv)
sum(comb_cor$ln_cv > comb_cor$rl_cv)
