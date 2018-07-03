#Load libraries 
library(reshape2)
library(reshape)
library(stats)
library(heplots)
library(mvpart)
library(car)
library(nlme)
library(lme4)
library(ggplot2)
library(pscl)
library(MASS)
library(MuMIn)
library(labdsv)
library(effects)
library(limma)
library(DESeq2)
library(edgeR)
library(splines)
library(RColorBrewer)
library(BiocParallel)
library(glmm)
library(MASS)

## Set working directory
setwd("/Users/courtney/Desktop/cavefish/RNAlater")

################################################################################################
#
#           RNASEQ LIQUID N AND RNALATER EXPRESSION
#           1. Bring in files and trim expession
#
################################################################################################

#Expression was calculated based on Stringtie/prepDE.py and GCcontent, gene length and exon number were extracted from Astyanax GTF.
expression <- read.csv("RNAlater_exp.csv",header=T, row.names=1)
GCcontent <- read.csv("GCcontentAsty.csv",header=T)
#Note, gene length is based on the transcript and gene rather than CDS, currently using the CDS which is in GCcontent
geneLength <- read.csv("gene_length.csv",header=T)
exonNumber <- read.csv("merged_exon.csv",header=T)
logfc <- read.csv("Logfc_expression.csv", header=T)
stats <- read.csv("condition_treated_results.csv", header=T)

# Count number of transcripts and number of individuals
dim(expression)

# Trim data: Each row (i.e. transcript) must have greater than 2 counts per million (cpm) and at least 3 of the individuals must have counts data
keep <- rowSums(cpm(expression)>2) >= 3     

# Keep data that meet the filtering criteria
RNAexpression <- expression[keep,]

# Count number of transcripts left after trimming
dim(RNAexpression)

#############NORMALIZE USING EDGER############################

# Subset quality samples from file
d <- subset(RNAexpression, select = c(CHOY.16.R.01, CHOY.16.R.03, CHOY.16.R.04, CHOY.16.R.05, CHOY.16.R.2,
                                 CHOY.16.01, CHOY.16.04, CHOY.16.05, CHOY.16.08, CHOY.16.11, CHOY.16.12))

# Group repliates
group <-c(rep("RNAlater",5),rep("liquidN",6))

# Trim out any transcripts with no counts
d <- d[rowSums(d) > 0, ]

# Count number of transcripts and individuals
dim(d)

#set up DGEList
dgeg <- DGEList(counts=d, group=group)

#extract counts per million that are normalized
counts.per.m <- cpm(dgeg, normalized.lib.sizes=TRUE)

#extract data in a matrix
A <- as.matrix(counts.per.m)

########################NORMALIZE USING DESEQ2############################
#Bring in files
expression <- read.csv("RNAlater_exp.csv",header=T, row.names=1)
#expression <- subset(expression, select=c("CHOY.16.01","CHOY.16.04","CHOY.16.05","CHOY.16.11","CHOY.16.12",
#                                   "CHOY.16.R.01","CHOY.16.R.03","CHOY.16.R.04","CHOY.16.R.05","CHOY.16.R.2"))
coldata <- read.csv("coldata.csv",header=T, row.names=1)
#check row names are in the same order
all(rownames(coldata) == colnames(expression))
#create DEGlist
dds <- DESeqDataSetFromMatrix(countData = expression, colData = coldata, design = ~ Condition)
#trim out reads that have less than 10countspemillion
keep <- rowSums(counts(dds) >= 10) >= 3
#keep the ones above
dds <- dds[keep,]
#check number of counts
dim(dds)
#normalize reads
normalized <- fpm(dds)

#runDESeq2
dds <- DESeq(dds)

#check results
res <- results(dds)

#Extract description infomation
mcols(res)$description
resOrdered <- res[order(res$pvalue),]
summary(res)
res05 <- results(dds, alpha=0.05)
summary(res05)
par(mar=c(8,5,2,2))
boxplot(log10(assays(dds)[["cooks"]]), range=0, las=2)
write.csv(as.data.frame(resOrdered), 
          file="condition_treated_results.csv")

vsd <- vst(dds, blind=FALSE)
plotPCA(vsd, intgroup=c("Condition"))

#USE LRT to calculate LogFC based on "Condition"
dds <- DESeq(dds, test="LRT", reduced=~1)
res <- results(dds)


#brought file back in as "stats"

#add empty column for treatment
stats$expr <- 0

#Change logFC padjusted values that are < 0.05 to 0 and and > 0.05 to 1
stats$expr <- ifelse(stats$padj <0.05, 0, 1)

#change headers
names(stats) <- c("gene_id","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj","expr")

#remove rows with NAs
final_expression <- stats[complete.cases(stats[ , 7]),]

#combine GC content 
storexp <- merge(final_expression[, c("gene_id","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj","expr")], GCcontent, by = c("gene_id"))

#combine exon number
Expression_stor <- merge(storexp[, c("gene_id","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj","expr", "gene_length","GC_length","GC_prop")], exonNumber, by = c("gene_id"))

#sort data based on gene ID and then variable (sample ID)
Expstorage <- Expression_stor[order(Expression_stor$gene_id),]

model1 <- lm(expr~GC_prop+exon_number+gene_length, data=Expstorage)     

model1 <- glm(expr~GC_prop+exon_number+gene_length, family=binomial, data=Expstorage)
model2 <- glm(expr~exon_number+gene_length, family=binomial, data=Expstorage)
model3 <- glm(expr~gene_length, family=binomial, data=Expstorage)
anova(model1, model2, model3)











d <- subset(Expstorage, select=c(logfc, gene_length))
ggscatter(Expstorage, x = "log2FoldChange", y = "gene_length", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "logfc", ylab = "gene length")
ggscatter(Expstorage, x = "log2FoldChange", y = "GC_prop", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "logfc", ylab = "GC_prop")
ggscatter(Expstorage, x = "log2FoldChange", y = "exon_number", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "logfc", ylab = "exon_number")

#2612
plotMA(res, ylim=c(-2,2))
plotCounts(dds, gene=which.min(res$padj), intgroup="Condition")
d <- plotCounts(dds, gene=which.min(res$padj), intgroup="Condition", 
                returnData=TRUE)
ggplot(d, aes(x=Condition, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10(breaks=c(25,100,400))

library("pheatmap")
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:20]
df <- as.data.frame(colData(dds)["Condition"],
pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df))

################################################################################################
#
#           RNASEQ LIQUID N AND RNALATER EXPRESSION
#           2. Combine files into one dataframe 
#
################################################################################################

#convert each row to a unique ID based on "gene_id"
exp <- melt(normalized) 

#change headers
names(exp) <- c("gene_id","sample_id","expression")

#combine GC content 
storexp <- merge(exp[, c("gene_id","sample_id","expression")], GCcontent, by = c("gene_id"))

#combine exon number
Expression_stor <- merge(storexp[, c("gene_id","sample_id","expression", "gene_length","GC_length","GC_prop")], exonNumber, by = c("gene_id"))

#add empty column for treatment
Expression_stor$treatment <- 0

#Add storage conditions (LiquidN and RNAlater)
Expression_stor$treatment <- ifelse(Expression_stor$sample_id == "CHOY.16.01", "LIQUIDN",
                                    ifelse(Expression_stor$sample_id == "CHOY.16.05", "LIQUIDN",
                                           ifelse(Expression_stor$sample_id == "CHOY.16.11", "LIQUIDN",
                                                  ifelse(Expression_stor$sample_id == "CHOY.16.12", "LIQUIDN",
                                                         ifelse(Expression_stor$sample_id == "CHOY.16.04", "LIQUIDN",
                                                                ifelse(Expression_stor$sample_id == "CHOY.16.08", "LIQUIDN",
                                                                       ifelse(Expression_stor$sample_id == "CHOY.16.R.01", "RNALATER",
                                                                              ifelse(Expression_stor$sample_id == "CHOY.16.R.03", "RNALATER",
                                                                                     ifelse(Expression_stor$sample_id == "CHOY.16.R.04", "RNALATER",
                                                                                            ifelse(Expression_stor$sample_id == "CHOY.16.R.05", "RNALATER",
                                                                                                   ifelse(Expression_stor$sample_id == "CHOY.16.R.2", "RNALATER", NA)))))))))))

#Expression_stor$treatment <- treat$treatment[match(Expression_stor$sample_id, treat$sample_id)]

#sort data based on gene ID and then variable (sample ID)
Expstorage <- Expression_stor[order(Expression_stor$gene_id, Expression_stor$sample_id),]

################################################################################################
#
#           RNASEQ LIQUID N AND RNALATER EXPRESSION
#           4. MDS plot and Heatmap
#
################################################################################################

summary(m1 <- glm.nb(expression ~ treatment + gene_length + GC_length + exon_number, data=Expstorage))
model1 <- glm(expression ~ treatment + gene_length + GC_length + exon_number, family = negative.binomial (2), data=Expstorage)
model2 <- glm(expression ~ gene_length + GC_length + exon_number, data=Expstorage)
anova(model1, model2)


################################################################################################
#
#           RNASEQ LIQUID N AND RNALATER EXPRESSION
#           3. Multiple regression
#
################################################################################################
#the mixed effect model
model1 <- lme(expression ~ treatment + gene_length + GC_length + exon_number, random=~1|gene_id, data=Expstorage)
summary(model1)
anova(model1)
model2 <- lme(expression ~ treatment + gene_length + GC_length + exon_number + treatment*gene_length + treatment*GC_length + treatment*exon_number + gene_length*GC_length + gene_length*exon_number + GC_length*exon_number, random=~1|gene_id, data=Expstorage)
summary(model2)
anova(model2)
model3 <- lme(expression ~ treatment + gene_length + GC_length + exon_number + treatment*gene_length + treatment*GC_length + treatment*exon_number + gene_length*GC_length + gene_length*exon_number + GC_length*exon_number + treatment*gene_length*GC_length + treatment*gene_length*exon_number + treatment*GC_length*exon_number + gene_length*GC_length*exon_number, random=~1|gene_id, data=Expstorage)
summary(model3)
anova(model3)

#glmm main effects + random effects
glmm_model <- MCMCglmm(expression~treatment+gene_length+GC_prop+exon_number,random=~gene_id,data=Expstorage,nitt=1300, burnin=300, thin=1, family="gaussian")
summary(glmm_model)

#Linear model main effects
fit <- lm(expression ~ treatment + gene_length + GC_prop + exon_number,  data=Expstorage)
summary(fit)
#Linear model main effects+ two-way interactions

#Compare model of treatment + full model
fit2 <- lmer(expression ~ treatment + (1 | sample_id), data=Expstorage, REML=FALSE)
fit3 <- lmer(expression ~ gene_length + GC_prop + exon_number + treatment + (1 | sample_id), data=Expstorage, REML=FALSE)
anova(fit2, fit3)

















################################################################################################
#
#           RNASEQ LIQUID N AND RNALATER EXPRESSION
#           2. LOGFC
#
################################################################################################
#change headers
names(logfc) <- c("gene_id","logfc")

#combine GC content 
storexp <- merge(logfc[, c("gene_id","logfc")], GCcontent, by = c("gene_id"))

#combine exon number
Expression_stor <- merge(storexp[, c("gene_id","logfc", "gene_length","GC_length","GC_prop")], exonNumber, by = c("gene_id"))

#add empty column for treatment
Expression_stor$treatment <- 0

#Note, can't add random effects
fit <- lm(logfc ~ gene_length + GC_prop + exon_number,  data=Expression_stor)
summary(fit)

#glmm main effects + random effects
glmm_model <- MCMCglmm(logfc ~ gene_length + GC_prop + exon_number,random=~gene_id,data=Expression_stor,nitt=1300, burnin=300, thin=1, family="gaussian")
summary(glmm_model)

negative_d <- subset(Expression_stor, logfc < 0, select =c(logfc, GC_prop, exon_number, gene_length, GC_length))






library(ggpubr)
#correlation between logfc and gene_length 
d <- subset(Expression_stor, select=c(logfc, gene_length))
ggscatter(d, x = "logfc", y = "gene_length", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "logfc", ylab = "gene length")

#correlation between logfc and exon_number
d <- subset(Expression_stor, select=c(logfc, exon_number))
ggscatter(d, x = "logfc", y = "exon_number", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "logfc", ylab = "exon number")

#correlation between logFC and GC proportion
d <- subset(Expression_stor, select=c(logfc, GC_prop))
ggscatter(d, x = "logfc", y = "GC_prop", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "logfc", ylab = "GC proportion")

#correlation between logFC and GC_length
d <- subset(Expression_stor, select=c(logfc, GC_length))
ggscatter(d, x = "logfc", y = "GC_length", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "logfc", ylab = "GC Length")

##Subset positive and negative logfc
positive_d <- subset(Expression_stor, logfc >0, select =c(logfc, GC_prop, exon_number, gene_length, GC_length))
ggscatter(positive_d, x = "logfc", y = "GC_length", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "logfc", ylab = "GC Length")
ggscatter(positive_d, x = "logfc", y = "GC_prop", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "logfc", ylab = "GC proportion")
ggscatter(positive_d, x = "logfc", y = "exon_number", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "logfc", ylab = "exon number")
ggscatter(positive_d, x = "logfc", y = "gene_length", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "logfc", ylab = "gene length")
negative_d <- subset(Expression_stor, logfc < 0, select =c(logfc, GC_prop, exon_number, gene_length, GC_length))
ggscatter(negative_d, x = "logfc", y = "GC_length", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "logfc", ylab = "GC Length")
ggscatter(negative_d, x = "logfc", y = "GC_prop", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "logfc", ylab = "GC proportion")
ggscatter(negative_d, x = "logfc", y = "exon_number", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "logfc", ylab = "exon number")
ggscatter(negative_d, x = "logfc", y = "gene_length", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "logfc", ylab = "gene length")

#Significant LogFC padjusted Log FC
#subset rows where padjusted > 0.05
significant <- subset(stats, padj < 0.05, select =c(gene_id, log2FoldChange, padj))

#combine GC content 
storexp <- merge(significant[, c("gene_id","log2FoldChange","padj")], GCcontent, by = c("gene_id"))

#combine exon number
Expression_stor <- merge(storexp[, c("gene_id","log2FoldChange","padj","gene_length","GC_length","GC_prop")], exonNumber, by = c("gene_id"))

#sort data based on gene ID and then variable (sample ID)
Expstorage <- Expression_stor[order(Expression_stor$gene_id),]

##Subset positive and negative logfc
positive_d <- subset(Expstorage, log2FoldChange >0, select =c(log2FoldChange, GC_prop, exon_number, gene_length, GC_length))
ggscatter(positive_d, x = "log2FoldChange", y = "GC_length", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "logfc", ylab = "GC Length")
ggscatter(positive_d, x = "log2FoldChange", y = "GC_prop", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "logfc", ylab = "GC proportion")
ggscatter(positive_d, x = "log2FoldChange", y = "exon_number", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "logfc", ylab = "exon number")
ggscatter(positive_d, x = "log2FoldChange", y = "gene_length", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "logfc", ylab = "gene length")
negative_d <- subset(Expstorage, log2FoldChange <0, select =c(log2FoldChange, GC_prop, exon_number, gene_length, GC_length))
ggscatter(negative_d, x = "log2FoldChange", y = "GC_length", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "logfc", ylab = "GC Length")
ggscatter(negative_d, x = "log2FoldChange", y = "GC_prop", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "logfc", ylab = "GC proportion")
ggscatter(negative_d, x = "log2FoldChange", y = "exon_number", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "logfc", ylab = "exon number")
ggscatter(negative_d, x = "log2FoldChange", y = "gene_length", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "logfc", ylab = "gene length")

#glmm main effects + random effects
glmm_model <- MCMCglmm(expr~gene_length+GC_prop+exon_number,random=~gene_id,data=Expstorage,nitt=1300, burnin=300, thin=1, family="gaussian")
summary(glmm_model)

#glmm main effects + random effects
glmm_model2 <- MCMCglmm(expr~gene_length+GC_prop+exon_number,data=Expstorage,nitt=1300, burnin=300, thin=1, family="gaussian")
summary(glmm_model2)

Anova(lm(y~x1*x2),type=2)

aa <- Anova(lm(log2FoldChange~GC_prop*gene_length+GC_prop*exon_number+gene_length*exon_number,data=Expstorage),type=3)












#####################################################################################################
#
#   2. 
#           
#
#####################################################################################################
#### Lab ####
# Group lab by exposure treatment and population
group = c(rep("TacoFreshControl",5), rep("TacoFreshHigh",5), rep("TacoFreshLow",6), rep("TacoFreshMedium",6),
          rep("PuyaFreshControl",6), rep("PuyaFreshHigh",3), rep("PuyaFreshLow",6), rep("PuyaFreshMedium",5),
          rep("PuyaSulfidicControl",6), rep("PuyaSulfidicHigh",3), rep("PuyaSulfidicLow",6), rep("PuyaSulfidicMedium",6),
          rep("TacoSulfidicControl",5), rep("TacoSulfidicHigh",5), rep("TacoSulfidicLow",6), rep("TacoSulfidicMedium",6))

# Generate dataframe with lab sample names and group 
samples = data.frame(cbind(colnames(lab)), as.character(group))
colnames(samples) = c("samples", "group")

# Create a DGEList object to hold the dataset 
m = DGEList(counts = lab, group = samples$group)

# Calculate normalized factors based on raw library sizes
m = calcNormFactors(m)

# Create a design marix 
design <- model.matrix(~0+group)

# Add column names based on sample names in the lab group
colnames(design) <- levels(factor(samples$group))

# Estimate common dispersion and tagwise dispersion 
m <- estimateDisp(m, design)

# Given tagwise dispersion and a design matrix, glmFIT fits the negative binomial GLM for each tag. Will then produce a DGEGLM object with new components. 
fit <- glmFit(m,design)

# Construct contrast matrix of comparisons (compare each exposure treatment to the control treatment for that drainage)
my.contrasts <- makeContrasts(
  Puya_low    = PuyaSulfidicLow - PuyaFreshControl,
  Puya_medium = PuyaSulfidicMedium - PuyaFreshControl,
  Puya_high   = PuyaSulfidicHigh - PuyaFreshControl,
  Taco_low    = TacoSulfidicLow - TacoFreshControl,
  Taco_medium = TacoSulfidicMedium - TacoFreshControl,
  Taco_high   = TacoSulfidicHigh - TacoFreshControl,
  levels = design
)

# Use the DGEGLM object to perform likelihood ratio test based on the contrast matrix of comparisons
PL_lrt = glmLRT(fit, contrast = my.contrasts[,"Puya_low"])
PM_lrt = glmLRT(fit, contrast = my.contrasts[,"Puya_medium"])
PH_lrt = glmLRT(fit, contrast = my.contrasts[,"Puya_high"])
TL_lrt = glmLRT(fit, contrast = my.contrasts[,"Taco_low"])
TM_lrt = glmLRT(fit, contrast = my.contrasts[,"Taco_medium"])
TH_lrt = glmLRT(fit, contrast = my.contrasts[,"Taco_high"])

# Summary of differential expression that was up, down or not significant in each comparison 
summary(decideTestsDGE(PL_lrt, adjust.method="BH", p.value = 0.05))
summary(decideTestsDGE(PM_lrt, adjust.method="BH", p.value = 0.05))
summary(decideTestsDGE(PH_lrt, adjust.method="BH", p.value = 0.05))

summary(decideTestsDGE(TL_lrt, adjust.method="BH", p.value = 0.05))
summary(decideTestsDGE(TM_lrt, adjust.method="BH", p.value = 0.05))
summary(decideTestsDGE(TH_lrt, adjust.method="BH", p.value = 0.05))

# Extract top expressed genes (based on significant upregulated and downregulated genes)
PL = topTags(PL_lrt, n = 473)
PM = topTags(PM_lrt, n = 284)
PH = topTags(PH_lrt, n = 1298)
TL = topTags(TL_lrt, n = 138)
TM = topTags(TM_lrt, n = 1058)
TH = topTags(TH_lrt, n = 71)

# Write logFC, logCPM, LR, Pvalue and FDR statistics for each significant gene to a .csv file 
write.csv(PL$table, file = "PL.csv")
write.csv(PM$table, file = "PM.csv")
write.csv(PH$table, file = "PH.csv")
write.csv(TL$table, file = "TL.csv")
write.csv(TM$table, file = "TM.csv")
write.csv(TH$table, file = "TH.csv")


# Field
# Group field by population
field.group = c(rep("PuyaFieldSulfidic", 5), rep("PuyaFieldFresh",6), rep("TacoFieldFresh",6), rep("TacoFieldSulfidic", 5))

# Generate dataframe with lab sample names and group 
field.samples = data.frame(cbind(colnames(field)), as.character(field.group))
colnames(field.samples) = c("samples", "group")

# Create a DGEList object to hold the dataset 
field.m = DGEList(counts = field, group = field.samples$group)

# Calculate normalized factors based on raw library sizes
field.m = calcNormFactors(field.m)

# Create a design marix 
field.design <- model.matrix(~0+field.group)

# Add column names based on sample names in the lab group
colnames(field.design) <- levels(factor(field.samples$group))

# Estimate common dispersion and tagwise dispersion 
field.m <- estimateDisp(field.m, field.design)

# Given tagwise dispersion and a design matrix, glmFIT fits the negative binomial GLM for each tag. Will then produce a DGEGLM object with new components. 
field.fit <- glmFit(field.m,field.design)

# Construct contrast matrix of comparisons (compare sulfide adapted to nonsulfide adapted populations in each drainage)
my.contrasts <- makeContrasts(
  Field_Puya = PuyaFieldSulfidic - PuyaFieldFresh, 
  Field_Taco = TacoFieldSulfidic - TacoFieldFresh,
  levels = field.design
)

# Use the DGEGLM object to perform likelihood ratio test based on the contrast matrix of comparisons
FieldP_lrt = glmLRT(field.fit, contrast = my.contrasts[,"Field_Puya"])
FieldT_lrt = glmLRT(field.fit, contrast = my.contrasts[,"Field_Taco"])

# Summary of differential expression that was up, down or not significant in each comparison 
summary(decideTestsDGE(FieldP_lrt, adjust.method="BH", p.value = 0.05))
summary(decideTestsDGE(FieldT_lrt, adjust.method="BH", p.value = 0.05))

# Extract top expressed genes (based on significant upregulated and downregulated genes)
FieldP = topTags(FieldP_lrt, n = 1793)
FieldT = topTags(FieldT_lrt, n = 2462)

# Write logFC, logCPM, LR, Pvalue and FDR statistics for each significant gene to a .csv file 
write.csv(FieldP$table, file = "FieldP.csv")
write.csv(FieldT$table, file = "FieldT.csv")