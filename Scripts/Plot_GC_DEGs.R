# Make boxplots of the GC content for significantly up-regulated genes and
# down regulated genes

# Read in the lists of up-regulated and down regulated genes
up <- as.character(read.table("/Users/tomkono/Dropbox/GitHub/RIS/SEM_RNALater/Results/DEGs/RNAlater_v_LN_UP_0.05.txt", header=FALSE)$V1)
down <- as.character(read.table("/Users/tomkono/Dropbox/GitHub/RIS/SEM_RNALater/Results/DEGs/RNAlater_v_LN_DOWN_0.05.txt", header=FALSE)$V1)
# Read in the predictors, which has GC content
gc <- read.csv("/Users/tomkono/Dropbox/GitHub/RIS/SEM_RNALater/Data/RNALater_Predictors.csv", header=TRUE)

# We have to change the T to G and remove the dot from the TX names
up.gene <- unlist(
    lapply(
        strsplit(up, "[.]"),
        function(x) {
            gsub("T", "G", x[1])
        }
        )
    )
down.gene <- unlist(
    lapply(
        strsplit(down, "[.]"),
        function(x) {
            gsub("T", "G", x[1])
        }
        )
    )

toplot <- data.frame(
    Value=c(
        gc[gc$gene_id %in% up.gene, "GC_prop"],
        gc[gc$gene_id %in% down.gene, "GC_prop"],
        gc$GC_prop
        ),
    Label=c(
        rep("R", length(gc[gc$gene_id %in% up.gene, "GC_prop"])),
        rep("L", length(gc[gc$gene_id %in% down.gene, "GC_prop"])),
        rep("A", nrow(gc))
        )
    )
toplot$Label <- factor(toplot$Label, levels=c("R", "L", "A"))

pdf(file="Up_Down_GC_Content.pdf", height=6, width=6)
boxplot(
    toplot$Value ~ toplot$Label,
    xlab="Category",
    ylab="GC Content",
    axes=F)
axis(side=2)
axis(side=1, at=c(1, 2, 3), labels=c("Higher exp.\nin RNAlater", "Higher exp.\nin LN2", "All Genes"), tck=-0.01)
dev.off()
