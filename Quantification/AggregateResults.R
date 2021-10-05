


### Command line arguments ###
args <- commandArgs(trailingOnly = TRUE)

numberOfFiles <- strtoi(args[1])
# can be geneEffectiveComp or geneComp or difference
inputFileType <-  args[2]
outputFilePath <- args[3]
outputFileName <- args[4]
bandProb <- args[5]



# numberOfFiles <- 2
# inputFileType <-  "geneEffectiveComp"
# outputFilePath <-  "/algbio1/shounak/raw/simulated/human/JustQuantification/test/"
# outputFileName <- "corr_Mod_File_Aggregated.txt"
# file1 <- "/algbio1/shounak/raw/simulated/human/Resimulation/30Million/simulations/sim_30000000/sim_5/betas_realistic/geneEffectiveComp/Mrd_Orig_File.txt"
# file2 <- "/algbio1/shounak/raw/simulated/human/Resimulation/30Million/simulations/sim_30000000/sim_2/betas_realistic/geneEffectiveComp/Mrd_Orig_File.txt"


aggregateFiles <- function(numberOfFiles,inputFileType,outputFilePath,outputFileName,bandProb){

  if(inputFileType == "geneEffectiveComp"){
    for (i in 1:numberOfFiles){
      currentFile <- read.table(file = args[i + 5], header = T)
      currentFile[3] <- currentFile[3]/numberOfFiles
      if(i == 1){
        finalFile <- currentFile
      }else{
        mergedFile <- merge(currentFile,finalFile, by = c("Complexity","Effective_complexity"))
        mergedFile$addedValue <- mergedFile[3] + mergedFile[6]
        n <- names(finalFile)
        finalFile <- as.data.frame(mergedFile[c(1,2,3,4,5)])
        finalFile[3] <- mergedFile$addedValue
        colnames(finalFile) <- n
      }
    }
    outputFile <- paste(outputFilePath,outputFileName,sep="")
    write.table(finalFile, file = outputFile, col.names = T, row.names = F, quote = F, sep = "\t")
  }else if(inputFileType == "geneComp"){
    for (i in 1:numberOfFiles){
      currentFile <- read.table(file = args[i + 5], header = T)
      currentFile[2] <- currentFile[2]/numberOfFiles
      if(i == 1){
        finalFile <- currentFile
      }else{
        mergedFile <- merge(currentFile,finalFile, by = c("Complexity"))
        mergedFile$addedValue <- mergedFile[2] + mergedFile[5]
        n <- names(finalFile)
        finalFile <- as.data.frame(mergedFile[c(1,2,3,4)])
        finalFile[2] <- mergedFile$addedValue
        colnames(finalFile) <- n
      }
    }
    finalFile["bandProb"] <- bandProb
    outputFile <- paste(outputFilePath,outputFileName,sep="")
    write.table(finalFile, file = outputFile, col.names = T, row.names = F, quote = F, sep = "\t")
  }else if(inputFileType == "difference"){
    for (i in 1:numberOfFiles){
      currentFile <- read.table(file = args[i + 5], header = T)
      currentFile[2] <- currentFile[2]/numberOfFiles
      if(i == 1){
        finalFile <- currentFile
      }else{
        mergedFile <- merge(currentFile,finalFile, by = c("diffGene_EffectiveGene"))
        mergedFile$addedValue <- mergedFile[2] + mergedFile[5]
        n <- names(finalFile)
        finalFile <- as.data.frame(mergedFile[c(1,2,3,4)])
        finalFile[2] <- mergedFile$addedValue
        colnames(finalFile) <- n
      }
    }
    outputFile <- paste(outputFilePath,outputFileName,sep="")
    write.table(finalFile, file = outputFile, col.names = T, row.names = F, quote = F, sep = "\t")
  }

}


aggregateFiles(numberOfFiles,inputFileType,outputFilePath,outputFileName,bandProb)
