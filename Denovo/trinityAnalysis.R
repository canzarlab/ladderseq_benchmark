
require(dplyr)

args = commandArgs(trailingOnly=TRUE)

args1 <- args[1]
args2 <- args[2]
args3 <- args[3]
args4 <- args[4]
args5 <- args[5]
args6 <- args[6]
args7 <- args[7]


# 
# ## hardcoded arguments
# args1 <- "/algbio1/shounak/raw/simulated/human/Resimulation/100MillWithDenovo/simulations/sim_100000000/sim_1/betas_realistic/blat/ladder_seq.psl"
# args2 <- "/algbio1/shounak/raw/simulated/human/Resimulation/100MillWithDenovo/simulations/sim_100000000/sim_1/betas_realistic/blat/full.psl"
# args3 <- "/algbio1/shounak/raw/simulated/human/GroundTruth/CoverageGroundTruth/100Mill/groupedGenes_fullyCov_0.1Tpm.fa"
# args4 <- "/algbio1/shounak/raw/simulated/human/Resimulation/100MillWithDenovo/simulations/sim_100000000/sim_1/betas_realistic/ladderMerge/ladderMergeFinal.fasta"
# args5 <- "/algbio1/shounak/raw/simulated/human/Resimulation/100MillWithDenovo/simulations/sim_100000000/sim_1/betas_realistic/original_trinity/Trinity.fasta"
# args6 <- "/algbio1/shounak/raw/simulated/human/Resimulation/100MillWithDenovo/simulations/sim_100000000/sim_1/betas_realistic/test"
# args7 <- "realistic_betas"



print("Starting BLAT Analysis")

colnames = c("match", "mis.match", "rep.match", "Ns.count", "Q.gap.bases", "Q.gap.count", "T.gap.bases",
	     "T.gap.count", "strand", "Q.name", "Q.size", "Q.start", "Q.end", "T.name", "T.size", "T.start",
	     "T.end", "block.count", "blockSizes", "qStarts", "tStarts")


merged <- read.table(args1, sep="\t", header=T, skip=5, col.names=colnames)
unmerged <- read.table(args2, sep="\t", header=T, skip=5, col.names=colnames)

## size filter - tolerence by percentage
sizeFilterMerged_80 <- merged[which(((merged[["Q.start"]]==0) & (merged[["Q.end"]] == merged[["Q.size"]])) & (((merged[["T.end"]] - merged[["T.start"]]) / merged[["T.size"]])*100 >= 80)),]
sizeFilterUnmerged_80 <- unmerged[which(((unmerged[["Q.start"]]==0) & (unmerged[["Q.end"]] == unmerged[["Q.size"]])) & (((unmerged[["T.end"]] - unmerged[["T.start"]]) / unmerged[["T.size"]])*100 >=  80)),]

sizeFilterMerged_85 <- merged[which(((merged[["Q.start"]]==0) & (merged[["Q.end"]] == merged[["Q.size"]])) & (((merged[["T.end"]] - merged[["T.start"]]) / merged[["T.size"]])*100 >= 85)),]
sizeFilterUnmerged_85 <- unmerged[which(((unmerged[["Q.start"]]==0) & (unmerged[["Q.end"]] == unmerged[["Q.size"]])) & (((unmerged[["T.end"]] - unmerged[["T.start"]]) / unmerged[["T.size"]])*100 >=  85)),]

sizeFilterMerged_90 <- merged[which(((merged[["Q.start"]]==0) & (merged[["Q.end"]] == merged[["Q.size"]])) & (((merged[["T.end"]] - merged[["T.start"]]) / merged[["T.size"]])*100 >= 90)),]
sizeFilterUnmerged_90 <- unmerged[which(((unmerged[["Q.start"]]==0) & (unmerged[["Q.end"]] == unmerged[["Q.size"]])) & (((unmerged[["T.end"]] - unmerged[["T.start"]]) / unmerged[["T.size"]])*100 >=  90)),]

sizeFilterMerged_95 <- merged[which(((merged[["Q.start"]]==0) & (merged[["Q.end"]] == merged[["Q.size"]])) & (((merged[["T.end"]] - merged[["T.start"]]) / merged[["T.size"]])*100 >= 95)),]
sizeFilterUnmerged_95 <- unmerged[which(((unmerged[["Q.start"]]==0) & (unmerged[["Q.end"]] == unmerged[["Q.size"]])) & (((unmerged[["T.end"]] - unmerged[["T.start"]]) / unmerged[["T.size"]])*100 >=  95)),]


## indel filters
indelFilterMerged_80 <- sizeFilterMerged_80[which((sizeFilterMerged_80[["Q.gap.count"]] + sizeFilterMerged_80[["T.gap.count"]])/sizeFilterMerged_80[["Q.size"]] < 0.01),]
indelFilterUnmerged_80 <- sizeFilterUnmerged_80[which((sizeFilterUnmerged_80[["Q.gap.count"]] + sizeFilterUnmerged_80[["T.gap.count"]])/sizeFilterUnmerged_80[["Q.size"]] < 0.01),]

indelFilterMerged_85 <- sizeFilterMerged_85[which((sizeFilterMerged_85[["Q.gap.count"]] + sizeFilterMerged_85[["T.gap.count"]])/sizeFilterMerged_85[["Q.size"]] < 0.01),]
indelFilterUnmerged_85 <- sizeFilterUnmerged_85[which((sizeFilterUnmerged_85[["Q.gap.count"]] + sizeFilterUnmerged_85[["T.gap.count"]])/sizeFilterUnmerged_85[["Q.size"]] < 0.01),]

indelFilterMerged_90 <- sizeFilterMerged_90[which((sizeFilterMerged_90[["Q.gap.count"]] + sizeFilterMerged_90[["T.gap.count"]])/sizeFilterMerged_90[["Q.size"]] < 0.01),]
indelFilterUnmerged_90 <- sizeFilterUnmerged_90[which((sizeFilterUnmerged_90[["Q.gap.count"]] + sizeFilterUnmerged_90[["T.gap.count"]])/sizeFilterUnmerged_90[["Q.size"]] < 0.01),]

indelFilterMerged_95 <- sizeFilterMerged_95[which((sizeFilterMerged_95[["Q.gap.count"]] + sizeFilterMerged_95[["T.gap.count"]])/sizeFilterMerged_95[["Q.size"]] < 0.01),]
indelFilterUnmerged_95 <- sizeFilterUnmerged_95[which((sizeFilterUnmerged_95[["Q.gap.count"]] + sizeFilterUnmerged_95[["T.gap.count"]])/sizeFilterUnmerged_95[["Q.size"]] < 0.01),]



## calculating numbers
numFilteredTranscriptsMerged_80 <- length(unique(indelFilterMerged_80[["T.name"]]))
numFilteredTranscriptsUnmerged_80 <- length(unique(indelFilterUnmerged_80[["T.name"]]))

numFilteredTranscriptsMerged_85 <- length(unique(indelFilterMerged_85[["T.name"]]))
numFilteredTranscriptsUnmerged_85 <- length(unique(indelFilterUnmerged_85[["T.name"]]))

numFilteredTranscriptsMerged_90 <- length(unique(indelFilterMerged_90[["T.name"]]))
numFilteredTranscriptsUnmerged_90 <- length(unique(indelFilterUnmerged_90[["T.name"]]))

numFilteredTranscriptsMerged_95 <- length(unique(indelFilterMerged_95[["T.name"]]))
numFilteredTranscriptsUnmerged_95 <- length(unique(indelFilterUnmerged_95[["T.name"]]))



transcriptomeSize <- as.numeric(args5)
mergedAlign <- as.numeric(args3)
unmergedAlign <- as.numeric(args4)

# transcriptomeSize <- 50533
# mergedAlign <- 71317
# unmergedAlign <- 61530


  

recall_num_merged_80 <- numFilteredTranscriptsMerged_80
recall_num_raw_80 <- numFilteredTranscriptsUnmerged_80
recall_merged_80 <- numFilteredTranscriptsMerged_80/transcriptomeSize
recall_raw_80 <- numFilteredTranscriptsUnmerged_80/transcriptomeSize
prec_merged_80 <- numFilteredTranscriptsMerged_80/mergedAlign
prec_raw_80 <- numFilteredTranscriptsUnmerged_80/unmergedAlign

recall_num_merged_85 <- numFilteredTranscriptsMerged_85
recall_num_raw_85 <- numFilteredTranscriptsUnmerged_85
recall_merged_85 <- numFilteredTranscriptsMerged_85/transcriptomeSize
recall_raw_85 <- numFilteredTranscriptsUnmerged_85/transcriptomeSize
prec_merged_85 <- numFilteredTranscriptsMerged_85/mergedAlign
prec_raw_85 <- numFilteredTranscriptsUnmerged_85/unmergedAlign


recall_num_merged_90 <- numFilteredTranscriptsMerged_90
recall_num_raw_90 <- numFilteredTranscriptsUnmerged_90
recall_merged_90 <- numFilteredTranscriptsMerged_90/transcriptomeSize
recall_raw_90 <- numFilteredTranscriptsUnmerged_90/transcriptomeSize
prec_merged_90 <- numFilteredTranscriptsMerged_90/mergedAlign
prec_raw_90 <- numFilteredTranscriptsUnmerged_90/unmergedAlign


recall_num_merged_95 <- numFilteredTranscriptsMerged_95
recall_num_raw_95 <- numFilteredTranscriptsUnmerged_95
recall_merged_95 <- numFilteredTranscriptsMerged_95/transcriptomeSize
recall_raw_95 <- numFilteredTranscriptsUnmerged_95/transcriptomeSize
prec_merged_95 <- numFilteredTranscriptsMerged_95/mergedAlign
prec_raw_95 <- numFilteredTranscriptsUnmerged_95/unmergedAlign





res <- data.frame(matrix(ncol = 5, nrow = 8))
colnames(res) <- c("Group", "Precision", "Recall", "PercentageLength", "bandProb")
res[["Group"]] <- c("Ladder-seq","Conventional","Ladder-seq","Conventional","Ladder-seq","Conventional","Ladder-seq","Conventional")
res[["Precision"]] <- c(prec_merged_80, prec_raw_80,prec_merged_85, prec_raw_85,prec_merged_90, prec_raw_90,prec_merged_95, prec_raw_95)
res[["Recall"]] <- c(recall_merged_80, recall_raw_80,recall_merged_85, recall_raw_85,recall_merged_90, recall_raw_90,recall_merged_95, recall_raw_95)
res[["Recall_num"]] <- c(recall_num_merged_80, recall_num_raw_80,recall_num_merged_85, recall_num_raw_85,recall_num_merged_90, recall_num_raw_90,recall_num_merged_95, recall_num_raw_95)
res[["PercentageLength"]] <- c(80,80,85,85,90,90,95,95)
res[["bandProb"]] <- args7




write.table(res, file=args6, row.names=FALSE, quote = FALSE, sep="\t")
