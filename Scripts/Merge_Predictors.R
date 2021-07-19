# This is a dumb script to merge the predictor variables for the RNALater
# project.

l_n_g <- read.csv("/Users/tomkono/Dropbox/GitHub/RIS/SEM_RNALater/Data/Amex_Ensembl_Len_Nexon_GC.csv.gz", header=TRUE)
repeats <- read.table("/Users/tomkono/Dropbox/GitHub/RIS/SEM_RNALater/Data/Amex102_Repeats.txt", header=TRUE)

merged <- merge(l_n_g, repeats, by.x="GeneID", by.y="gene_id", all.x=TRUE, all.y=TRUE)

write.csv(
    merged,
    file="RNALater_Predictors.csv",
    quote=FALSE,
    row.names=FALSE)
