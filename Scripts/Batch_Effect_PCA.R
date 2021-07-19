# Make a PC plot of the entire dataset, color-coded by batch and symbol-coded
# by population.

# We have already saved the normalized and transformed counts in this data file
load("/Users/tomkono/Dropbox/GitHub/RIS/SEM_RNALater/Results/VarExplained/RNAlater_VarExplained.rda")

# Decompose the counts with PCA
pc <- prcomp(t(vst_dds))

# Set up symbols for the populations:
#   - Choy N=42
#   - Molino N=42
#   - Pachon N=42
#   - Tinaja N=42
pchs <- c(rep(19, 42), rep(17, 42), rep(15, 42), rep(18, 42))
# Set up the colors for the batches.
cols <- rainbow(20)[as.numeric(metadata$ExtractBatch)]

# Make the plot
pdf(file="RNAseq_Batch_Effect_PCA.pdf", height=10, width=10)
plot(
    pc$x[,"PC1"], pc$x[,"PC2"],
    pch=pchs,
    col=cols,
    xlab="PC1",
    ylab="PC2")
legend("topleft",
    c("Choy", "Molino", "Pachon", "Tinaja"),
    pch=c(19, 17, 15, 18),
    ncol=1,
    title="Population")
legend("topright",
    as.character(levels(metadata$ExtractBatch)),
    pch=19,
    col=rainbow(20),
    ncol=4,
    title="Extraction Batch")
dev.off()
