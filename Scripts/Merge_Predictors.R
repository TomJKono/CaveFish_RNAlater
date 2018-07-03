# This is a dumb script to merge the predictor variables for the RNALater
# project. This is because the values are split across multiple files, and it
# is easier to just have one big file with all the info in it.

GCcontent <- read.csv("/Users/tomkono/Dropbox/GitHub/RIS/SEM_RNALater/Data/GCcontentAsty.csv",header=TRUE)
geneLength <- read.csv("/Users/tomkono/Dropbox/GitHub/RIS/SEM_RNALater/Data/gene_length.csv",header=TRUE)
exonNumber <- read.csv("/Users/tomkono/Dropbox/GitHub/RIS/SEM_RNALater/Data/merged_exon.csv",header=TRUE)
repeats <- read.table("/Users/tomkono/Dropbox/GitHub/RIS/SEM_RNALater/Data/Amex102_Repeats.txt", header=TRUE)

merged_all <- Reduce(
    function(...) merge(..., by="gene_id", all.x=TRUE),
    list(GCcontent, geneLength, exonNumber, repeats)
    )

write.csv(
    merged_all,
    file="RNALater_Predictors.csv",
    quote=FALSE,
    row.names=FALSE)
