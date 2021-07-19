# Estimate the variance explained by various experimental variables, like
# extraction batch and sequencing lane. We will be using the
# variancePartition package from Bioconductor

library("variancePartition")
library("DESeq2")
library("doParallel")

# Run 4 threads
cl <- makeCluster(3)
registerDoParallel(cl)

metadata <- read.csv("/Users/tomkono/Dropbox/GitHub/RIS/SEM_RNALater/Data/MasterSheet.test.for.Batch.effects_Rev.csv", header=TRUE)
gene_counts <- read.csv("/Users/tomkono/Dropbox/GitHub/RIS/SEM_RNALater/Data/Circadian_RawGeneCounts.csv", header=TRUE)

# Fix the samplenames in the metadata sheet so that they match the names in
# the gene counts matrix
metadata$Sample <- gsub("-", ".", metadata$Sample)
metadata$Sample <- gsub("[.]N", "", metadata$Sample)

# Extract the time point and population from the sample name
metadata$Pop <- unlist(lapply(strsplit(as.character(metadata$Sample), "[.]"), function(x) {x[1]}))
metadata$Time <- unlist(lapply(strsplit(as.character(metadata$Sample), "[.]"), function(x) {x[2]}))

# And make sure all our experimental variables are factors
metadata$SequencingLane <- as.factor(metadata$SequencingLane)
metadata$ExtractBatch <- as.factor(metadata$ExtractBatch)
metadata$Pop <- as.factor(metadata$Pop)
metadata$Time <- as.factor(metadata$Time)

# Trim down the counts matrix to only samples that are present in the metadata
# sheet.
keep_exp <- names(gene_counts) %in% metadata$Sample
flt_gene_counts <- gene_counts[,keep_exp]
rownames(flt_gene_counts) <- as.character(gene_counts$X)

# Then, we normalize the counts using only library size
dds <- DESeqDataSetFromMatrix(
    countData=as.matrix(flt_gene_counts),
    colData=metadata,
    design=~1)
dds <- estimateSizeFactors(dds)

# Remove genes with less than 100 counts across all samples
keep_gene <- rowSums(counts(dds)) >= 100
dds <- dds[keep_gene,]

# And extract variance-stabilized normalized counts
vst_dds <- assay(vst(dds, blind=TRUE))

# Define the forumla with which to test variance. We will include interactions
# in this first, full, model.
fm_i <- ~ (1|ExtractBatch) + (1|SequencingLane) + (1|Pop) + (1|Time) + 
        (1|ExtractBatch:SequencingLane) + (1|ExtractBatch:Pop) + (1|ExtractBatch:Time) +
        (1|SequencingLane:Pop) + (1|SequencingLane:Time) +
        (1|Pop:Time)
fm_me <- ~ (1|ExtractBatch) + (1|SequencingLane) + (1|Pop) + (1|Time)
fm_fixed <- ~ ExtractBatch + SequencingLane + Pop + Time
# Fit the model
vp_me <- fitExtractVarPartModel(vst_dds, fm_me, metadata)
vp_i <- fitExtractVarPartModel(vst_dds, fm_i, metadata)
vp_fixed <- fitExtractVarPartModel(vst_dds, fm_fixed, metadata)

# Plot it
pdf(file="RNAseq_Variance_Explained_MainEffects.pdf", height=6, width=6)
plotVarPart(vp_me)
dev.off()
pdf(file="RNAseq_Variance_Explained_FullModel.pdf", height=6, width=6)
plotVarPart(vp_i)
dev.off()
pdf(file="RNAseq_Variance_Explained_FixedEffects.pdf", height=6, width=6)
plotVarPart(vp_fixed)
dev.off()

# Save the model fit summaries
vp_me_fit <- fitVarPartModel(vst_dds, fm_me, metadata, fxn=summary)
vp_i_fit <- fitVarPartModel(vst_dds, fm_i, metadata, fxn=summary)
vp_fixed_fit <- fitVarPartModel(vst_dds, fm_fixed, metadata, fxn=summary)

print(vp_me_fit)
print(vp_i_fit)
print(vp_fixed_fit)

save(
    vst_dds,
    fm_i,
    fm_me,
    fm_fixed,
    vp_me,
    vp_i,
    vp_fixed,
    metadata,
    gene_counts,
    vp_me_fit,
    vp_i_fit,
    vp_i_fixed,
    file="RNAlater_VarExplained.rda")
