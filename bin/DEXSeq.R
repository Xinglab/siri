args <- commandArgs(trailingOnly = TRUE)
### using intron and exon data as input
library(DEXSeq)
library(BiocParallel)
BPPARAM = MulticoreParam(as.integer(args[length(args)]))
##agrs[1] = intron_data
##args[2] = sampleData
##args[3] = alternativeCountData
##args[4] = threads

intron = read.table(args[1],row.names = 1)
rownames(intron) = paste(intron[,1], rownames(intron),sep=":")
intron = intron[,-1]
countData = intron
featureID = gsub(".*:","",rownames(countData))
groupID = gsub(":.*","",rownames(countData))
groupID = gsub(",.*","",groupID)

sampleData = read.table(args[2])
sampleData = as.data.frame(sampleData)
alternativeCountData = read.table(args[3])
colnames(alternativeCountData) = colnames(countData)
alternativeCountData  = as.matrix(alternativeCountData )
dxd = DEXSeqDataSet(countData, sampleData, design= ~ sample + exon + condition:exon , featureID, groupID, featureRanges=NULL, transcripts=NULL, alternativeCountData)

dxd = estimateSizeFactors( dxd )
dxd = estimateDispersions( dxd, BPPARAM = BPPARAM)

dxd = testForDEU( dxd , BPPARAM = BPPARAM)
dxd = estimateExonFoldChanges( dxd, fitExpToVar="condition", BPPARAM = BPPARAM)
dxr1 = DEXSeqResults( dxd )
write.table(dxr1,args[length(args) -1],quote = FALSE, sep="\t")
