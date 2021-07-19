# Make a plot of the RNAlater predictor variables. We will make one that is
# potentially "publication ready" that just has the significant predictors
# and one that has all the predictors that were tested.

library(GGally)

# Read the data
dat <- read.csv("/Users/tomkono/Dropbox/GitHub/RIS/SEM_RNALater/Data/RNALater_Predictors.csv", header=TRUE)
l2fc <- read.csv("/Users/tomkono/Dropbox/GitHub/RIS/SEM_RNALater/Results/RNAlater_v_LN_Pvalues.csv", header=TRUE)

dat <- merge(dat, l2fc, by.x="GeneID", by.y="ENSID")

# Slice down the data frame to just the variables for the main text plot
full_dat <- dat[,c("log2FoldChange", "baseMean", "Longest.TxLen", "Longest.GCProp", "Longest.NExon", "SSR", "HPR")]
full_dat$baseMean <- log2(full_dat$baseMean)

# Make the SSR and HPR values factors
full_dat$SSR <- as.factor(full_dat$SSR)
full_dat$HPR <- as.factor(full_dat$HPR)

full_dat <- full_dat[!(is.na(full_dat$HPR) | is.na(full_dat$SSR)),]

small_cont <- function(data, mapping, ...) {
    ggplot(data=data, mapping=mapping) +
    geom_point(size=0.25)
}

# Now, make a plot
pdf(file="Predictor_Scatterplot.pdf", height=6, width=6)
par(mar=c(0, 0, 0, 0), mgp=c(1.5, 0.5, 0))
ggpairs(
    full_dat,
    columns=c("log2FoldChange", "baseMean", "Longest.GCProp", "Longest.NExon", "SSR"),
    columnLabels=c("L2FC", "log2(M)", "G", "E", "S"),
    axisLabels="show") +
theme_bw() +
theme(
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank()
    )
dev.off()

pdf(file="Predictor_Scatterplot_Supp.pdf", height=8, width=8)
par(mar=c(0, 0, 0, 0), mgp=c(1.5, 0.5, 0))
ggpairs(
    full_dat,
    columns=c("log2FoldChange", "baseMean", "Longest.TxLen", "Longest.GCProp", "Longest.NExon", "SSR", "HPR"),
    columnLabels=c("L2FC", "log2(M)", "L", "G", "E", "S", "H"),
    axisLabels="show") +
theme_bw() +
theme(
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank()
    )
dev.off()
