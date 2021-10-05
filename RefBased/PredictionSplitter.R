

require(dplyr)

# Hardcoded argument
# grouped_locus_transidFile <- "/algbio1/shounak/raw/simulated/human/ImperfectSep100HighRes/Assemblies/AssemblyByCoverage/tpm0.1/substituted/singleExon_NoAlign/test/locus_transidList.txt"
# groundTransIdFile <- "/algbio1/shounak/raw/simulated/human/ImperfectSep100HighRes/Assemblies/AssemblyByCoverage/tpm0.1/substituted/singleExon_NoAlign/test/groundTruth_transidList.txt"
# outputPath = ""



#### Command line arguments
args = commandArgs(trailingOnly=TRUE)

grouped_locus_transidFile <- args[1]

groundTransIdFile <- args[2]
## Output path
outputPath <- args[3]




filterLoci <- function(grouped_locus_transidFile,groundTransIdFile,outputPath){
  pred_ground_grouped <- read.delim(grouped_locus_transidFile, header = F, comment.char="#",quote = "", sep = "\t")
  names(pred_ground_grouped) <- c("LocusId","TranscriptId")
  #cc <- pred_ground_grouped %>%  subset(LocusId == "Locus_11746")

  groundtruthTransIds <- read.delim(groundTransIdFile, header = F, comment.char="#",quote = "")
  names(groundtruthTransIds) <- c("TranscriptId")
  groundtruthTransIds <- as.data.frame(unique(groundtruthTransIds$TranscriptId))
  names(groundtruthTransIds) <- c("TranscriptId")


  locusToBeKept <- unique(pred_ground_grouped %>% subset(TranscriptId %in% groundtruthTransIds$TranscriptId) %>%select(LocusId))
  length(locusToBeKept$LocusId)
  #dd <- locusToBeKept %>%  subset(LocusId == "Locus_11746")

  outFileName <- paste(outputPath,"locusTobeKept.txt",sep="")

  write.table(locusToBeKept,file = outFileName, sep = "\t",row.names = F, col.names = F, quote = F)
}


filterLoci(grouped_locus_transidFile,groundTransIdFile,outputPath)
