# A slight re-write of the CNP RNALater.r script for current/modern packages.
# This script will summarize raw counts and then drop a .CSV file with the
# normalized log2(fold change) values for future testing.
library("DESeq2")
library("pheatmap")
library("ggplot2")
library("dplyr")

# Set colors for plots
rnalater_col <- "#e41a1c"
ln_col <- "#377eb8"

# Read the data
#Expression was calculated based on Stringtie/prepDE.py and GCcontent, gene length and exon number were extracted from Astyanax GTF.
counts <- read.csv("/Users/tomkono/Dropbox/GitHub/RIS/SEM_RNALater/Data/RNAlater_exp.csv",header=TRUE, row.names=1)
exp_dat <- read.csv("/Users/tomkono/Dropbox/GitHub/RIS/SEM_RNALater/Data/coldata.csv", header=TRUE)
preds <- read.csv("/Users/tomkono/Dropbox/GitHub/RIS/SEM_RNALater/Data/RNALater_Predictors.csv", header=TRUE)

# Subset the expression matrix to just the samples we are comparing
sub_exp <- counts[, as.character(exp_dat$samplename)]

# Make a DESeqDataSet
dds <- DESeqDataSetFromMatrix(
    countData=sub_exp,
    colData=exp_dat,
    design=~Condition)

# Make summary plots of the read counts
pdf(file="RNAlater_RawCounts.pdf", height=8, width=11)
par(mar=c(8, 4, 0.1, 0.1), mgp=c(2, 1, 0))
toplot <- data.frame(
    Val=c(
        log(1+counts(dds)[, "CHOY.16.01"]), log(1+counts(dds)[, "CHOY.16.04"]),
        log(1+counts(dds)[, "CHOY.16.05"]), log(1+counts(dds)[, "CHOY.16.08"]),
        log(1+counts(dds)[, "CHOY.16.11"]), log(1+counts(dds)[, "CHOY.16.12"]),
        log(1+counts(dds)[, "CHOY.16.R.01"]), log(1+counts(dds)[, "CHOY.16.R.03"]),
        log(1+counts(dds)[, "CHOY.16.R.04"]), log(1+counts(dds)[, "CHOY.16.R.05"]),
        log(1+counts(dds)[, "CHOY.16.R.2"])),
    Lab=c(
        rep("CHOY.16.01", nrow(dds)), rep("CHOY.16.04", nrow(dds)),
        rep("CHOY.16.05", nrow(dds)), rep("CHOY.16.08", nrow(dds)),
        rep("CHOY.16.11", nrow(dds)), rep("CHOY.16.12", nrow(dds)),
        rep("CHOY.16.R.01", nrow(dds)), rep("CHOY.16.R.03", nrow(dds)),
        rep("CHOY.16.R.04", nrow(dds)), rep("CHOY.16.R.05", nrow(dds)),
        rep("CHOY.16.R.2", nrow(dds)))
    )
# Reorder the labels
toplot$Lab <- factor(toplot$Lab, 
    levels=c("CHOY.16.01", "CHOY.16.04", "CHOY.16.05", "CHOY.16.08", 
             "CHOY.16.11", "CHOY.16.12", "CHOY.16.R.01", "CHOY.16.R.03",
             "CHOY.16.R.04", "CHOY.16.R.05", "CHOY.16.R.2"))
# Plot it
boxplot(
    toplot$Val ~ toplot$Lab,
    at=c(
        seq(1, 6),
        seq(8, 12)
        ),
    las=2,
    border=c(rep(ln_col, 6), rep(rnalater_col, 5)),
    main="\nRaw Counts",
    ylab="log(1+Counts)",
    xlab="",
    ylim=c(0, 16),
    lwd=3)
legend("topright", c("LN2", "RNAlater"), col=c(ln_col, rnalater_col), lwd=3)
dev.off()

# Trim down the counts matrix to just genes with at least 100 counts across
# all samples
keep <- rowSums(counts(dds)) >= 100
dim(dds)
dds <- dds[keep,]
dim(dds)

# Then re-plot the data for the filtered set
pdf(file="RNAlater_FilteredCounts.pdf", height=8, width=10)
par(mar=c(7.5, 4, 0, 0), mgp=c(2, 1, 0))
toplot <- data.frame(
    Val=c(
        log(1+counts(dds)[, "CHOY.16.01"]), log(1+counts(dds)[, "CHOY.16.04"]),
        log(1+counts(dds)[, "CHOY.16.05"]), log(1+counts(dds)[, "CHOY.16.08"]),
        log(1+counts(dds)[, "CHOY.16.11"]), log(1+counts(dds)[, "CHOY.16.12"]),
        log(1+counts(dds)[, "CHOY.16.R.01"]), log(1+counts(dds)[, "CHOY.16.R.03"]),
        log(1+counts(dds)[, "CHOY.16.R.04"]), log(1+counts(dds)[, "CHOY.16.R.05"]),
        log(1+counts(dds)[, "CHOY.16.R.2"])),
    Lab=c(
        rep("CHOY.16.01", nrow(dds)), rep("CHOY.16.04", nrow(dds)),
        rep("CHOY.16.05", nrow(dds)), rep("CHOY.16.08", nrow(dds)),
        rep("CHOY.16.11", nrow(dds)), rep("CHOY.16.12", nrow(dds)),
        rep("CHOY.16.R.01", nrow(dds)), rep("CHOY.16.R.03", nrow(dds)),
        rep("CHOY.16.R.04", nrow(dds)), rep("CHOY.16.R.05", nrow(dds)),
        rep("CHOY.16.R.2", nrow(dds)))
    )
# Reorder the labels
toplot$Lab <- factor(toplot$Lab, 
    levels=c("CHOY.16.01", "CHOY.16.04", "CHOY.16.05", "CHOY.16.08", 
             "CHOY.16.11", "CHOY.16.12", "CHOY.16.R.01", "CHOY.16.R.03",
             "CHOY.16.R.04", "CHOY.16.R.05", "CHOY.16.R.2"))
# Plot it
boxplot(
    toplot$Val ~ toplot$Lab,
    at=c(
        seq(1, 6),
        seq(8, 12)
        ),
    las=2,
    border=c(rep(ln_col, 6), rep(rnalater_col, 5)),
    main="\nFiltered Counts",
    ylab="log(1+Counts)",
    xlab="",
    ylim=c(0, 16),
    lwd=3)
legend("topright", c("LN2", "RNAlater"), col=c(ln_col, rnalater_col), lwd=3)
dev.off()

# Make a PC
pdf(file="RNAlater_FilteredCounts_PCA_Pub.pdf", height=3, width=3)
par(mar=c(2.5, 2.5, 0, 0), mgp=c(1.5, 0.5, 0))
cts <- t(assay(vst(dds, blind=FALSE), normalized=TRUE))
pc <- prcomp(cts)
liquidn <- exp_dat$samplename[exp_dat$Condition == "LN2"]
rnalater <- exp_dat$samplename[exp_dat$Condition == "RNAlater"]
plot(
    pc$x[liquidn, "PC1"], pc$x[liquidn, "PC2"],
    col=ln_col,
    pch=19,
    xlim=c(-50, 50),
    ylim=c(-50, 50),
    xlab="PC1 (27.2% Var Explained)",
    ylab="PC2 (14.9% Var Explained)")
points(
    pc$x[rnalater, "PC1"], pc$x[rnalater, "PC2"],
    col=rnalater_col,
    pch=19)
text(pc$x[liquidn, "PC1"], pc$x[liquidn, "PC2"], cex=0.5, srt=-30, pos=2, offset=0.25, labels=liquidn, col=ln_col)
text(pc$x[rnalater, "PC1"], pc$x[rnalater, "PC2"], cex=0.5, srt=-30, pos=2, offset=0.25, labels=rnalater, col=rnalater_col)
legend("topright", c("LN2", "RNAlater"), col=c(ln_col, rnalater_col), pch=19, cex=0.5)

dev.off()
pdf(file="RNAlater_FilteredCounts_PCA.pdf", 6, 6)
cts <- t(assay(vst(dds, blind=FALSE), normalized=TRUE))
pc <- prcomp(cts)
summary(pc)
liquidn <- exp_dat$samplename[exp_dat$Condition == "LN2"]
rnalater <- exp_dat$samplename[exp_dat$Condition == "RNAlater"]
plot(
    pc$x[liquidn, "PC1"], pc$x[liquidn, "PC2"],
    col=ln_col,
    pch=19,
    xlim=c(-50, 50),
    ylim=c(-50, 50),
    xlab="PC1",
    ylab="PC2")
points(
    pc$x[rnalater, "PC1"], pc$x[rnalater, "PC2"],
    col=rnalater_col,
    pch=19)
text(pc$x[liquidn, "PC1"], pc$x[liquidn, "PC2"], cex=0.75, adj=c(1, 1), labels=liquidn, col=ln_col)
text(pc$x[rnalater, "PC1"], pc$x[rnalater, "PC2"], cex=0.75, adj=c(1, 1), labels=rnalater, col=rnalater_col)
legend("topright", c("LN2", "RNAlater"), col=c(ln_col, rnalater_col), pch=19)

plot(
    pc$x[liquidn, "PC2"], pc$x[liquidn, "PC3"],
    col=ln_col,
    pch=19,
    xlim=c(-50, 50),
    ylim=c(-50, 50),
    xlab="PC2",
    ylab="PC3")
points(
    pc$x[rnalater, "PC2"], pc$x[rnalater, "PC3"],
    col=rnalater_col,
    pch=19)
text(pc$x[liquidn, "PC2"], pc$x[liquidn, "PC3"], cex=0.75, adj=c(1, 1), labels=liquidn, col=ln_col)
text(pc$x[rnalater, "PC2"], pc$x[rnalater, "PC3"], cex=0.75, adj=c(1, 1), labels=rnalater, col=rnalater_col)
legend("topright", c("LN2", "RNAlater"), col=c(ln_col, rnalater_col), pch=19)

plot(
    pc$x[liquidn, "PC1"], pc$x[liquidn, "PC3"],
    col=ln_col,
    pch=19,
    xlim=c(-50, 50),
    ylim=c(-50, 50),
    xlab="PC1",
    ylab="PC3")
points(
    pc$x[rnalater, "PC1"], pc$x[rnalater, "PC3"],
    col=rnalater_col,
    pch=19)
text(pc$x[liquidn, "PC1"], pc$x[liquidn, "PC3"], cex=0.75, adj=c(1, 1), labels=liquidn, col=ln_col)
text(pc$x[rnalater, "PC1"], pc$x[rnalater, "PC3"], cex=0.75, adj=c(1, 1), labels=rnalater, col=rnalater_col)
legend("topright", c("LN2", "RNAlater"), col=c(ln_col, rnalater_col), pch=19)

# Make a barplot of the proportion of variance explained
eig <- pc$sdev^2
propvar <- eig/sum(eig)
at <- barplot(propvar, xlab="PC", ylab="Variance Explained", axes=FALSE, ylim=c(0, 0.6))
text(at, 0.5, labels=round(propvar, 3), srt=90, adj=0)
axis(side=1, at=at, labels=1:length(at))
axis(side=2, at=seq(0, 0.5, by=0.05), labels=seq(0, 0.5, by=0.05))
dev.off()

# Then, do the differential expression analysis and get l2fc values
dds$Condition <- relevel(dds$Condition, ref="LN2")
dds <- estimateSizeFactors(dds)
dds <- DESeq(dds)
res <- results(dds, name="Condition_RNAlater_vs_LN2")
summary(res)

# Make a heatmap because why not!
degs <- dds[res$padj < 0.05 & !is.na(res$padj),]
dim(degs)
if(nrow(degs) > 0) {
    pdf(file="RNAlater_DEG_Heatmap.pdf", height=6, width=6, onefile=FALSE)
    par(mar=c(2.5, 2.5, 0, 0), mgp=c(1.5, 0.5, 0))
    df <- as.data.frame(colData(degs)[,c("Condition")])
    annotation_colors <- list(
        Condition=c(LN2=ln_col, RNAlater=rnalater_col)
        )
    rownames(df) <- colnames(degs)
    colnames(df) <- c("Condition")
    vs_degs <- rlog(degs, blind=FALSE)
    mat <- assay(vs_degs) - rowMeans(assay(vs_degs))
    pheatmap(mat, cluster_rows=TRUE, show_rownames=FALSE, cluster_cols=TRUE, annotation_col=df, annotation_colors=annotation_colors, breaks=seq(-2, 2, length.out=101))
    dev.off()
}

# Save the list, just in case we want to have actual p-values for these
deg <- as.data.frame(res)
deg$ENSID <- rownames(deg)
rownames(deg) <- NULL
deg_ordered <- deg[order(deg$padj),]
deg_ordered$RNAlater.v.LN <- ifelse(deg_ordered$log2FoldChange>0, "Up", "Down")
deg_ordered <- deg_ordered[,c( "ENSID", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj", "RNAlater.v.LN")]
write.csv(deg_ordered, file="RNAlater_v_LN_Pvalues.csv", quote=FALSE, row.names=FALSE)
