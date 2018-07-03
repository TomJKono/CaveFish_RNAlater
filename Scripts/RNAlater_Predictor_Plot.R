# Make a plot of the RNAlater predictor variables. We will make one that is
# potentially "publication ready" that just has the significant predictors
# and one that has all the predictors that were tested.

# Read the data
dat <- read.csv("/Users/tomkono/Dropbox/GitHub/RIS/SEM_RNALater/Data/RNALater_Predictors.csv", header=TRUE)
l2fc <- read.csv("/Users/tomkono/Dropbox/GitHub/RIS/SEM_RNALater/Results/RNAlater_v_LN_Pvalues.csv", header=TRUE)

dat <- merge(dat, l2fc, by.x="gene_id", by.y="ENSID")

# Slice down the data frame to just the variables for the main text plot
full_dat <- dat[,c("baseMean", "gene_length", "GC_prop", "exon_number", "SSR", "HPR")]

# Now, make a plot
pdf(file="Predictor_Scatterplot.pdf", height=3, width=3)
par(mar=c(0, 0, 0, 0), mgp=c(1.5, 0.5, 0))
pairs(
    ~baseMean + GC_prop + exon_number + jitter(HPR),
    data=full_dat,
    pch=19,
    cex=0.25,
    col=rgb(0, 0, 0, 0.25),
    labels=c("M", "G", "E", "H"),
    cex.labels=1.5,
    lower.panel=NULL,
    oma=c(1, 1, 1, 1))
dev.off()

pdf(file="Predictor_Scatterplot_Supp.pdf", height=6, width=6)
par(mar=c(0, 0, 0, 0), mgp=c(1.5, 0.5, 0))
pairs(
    ~baseMean + gene_length + GC_prop + exon_number + jitter(SSR) + jitter(HPR),
    data=full_dat,
    pch=19,
    cex=0.25,
    col=rgb(0, 0, 0, 0.25),
    labels=c("M", "L", "G", "E", "S", "H"),
    cex.labels=1.5,
    lower.panel=NULL,
    oma=c(1, 1, 1, 1))
dev.off()
