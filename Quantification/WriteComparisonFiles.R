# Writes the files containing the correlation and mards and median relative difference (with modifications on original kallisto) based on effective gene complexity between ground truth and original/modified kallisto
# Shounak


# Libraries
require(ggplot2)
require(dplyr)





# ### Command line arguments ###
args <- commandArgs(trailingOnly = TRUE)

modKal <- args[1]
origKal <- args[2]
groundTruth <- args[3]
dataset <- args[4]
geneId_transId_file <- args[5]
outputPath <- args[6]
groundTruthType <- args[7]
# scale <- as.double(args[8])


### Harcoded arguments ###
## Harcoded arguments ###
#modKal <- "/algbio1/shounak/raw/simulated/human/Resimulation/30Million/simulations/sim_30000000/sim_1/betas_realistic/outputMod/abundance.tsv"
#origKal <- "/algbio1/shounak/raw/simulated/human/Resimulation/30Million/simulations/sim_30000000/sim_1/betas_realistic/original/abundance.tsv"
#groundTruth <- "/algbio1/shounak/raw/simulated/human/configData/realData/rsem/out.isoforms.results"
#dataset <- "check"
#geneId_transId_file <- "/algbio1/shounak/raw/simulated/human/configData/geneId_transId.txt"
#outputPath <- "/algbio1/shounak/raw/simulated/human/JustQuantification/test/"
#groundTruthType <- "RSEM"
# scale <- 0.4345646




mards <- function(vector1,vector2){
  numerator = vector1 - vector2
  denom = (vector1 + vector2)
  ards <- abs(numerator)/denom
  ards[which(!is.finite(ards))] <- 0
  return (mean(ards))
}

## Function taken from mamabear
relative_difference <- function(x, y, na_zeroes = FALSE, normalize_counts = TRUE) {

    stopifnot(length(x) == length(y))
  result <- rep(NA_real_, length(x))

  non_zero <- which( x > 0 | y > 0 )
  both_zero <- setdiff(seq_along(x), non_zero)


  if (!na_zeroes) {
    result[both_zero] <- 0.0
  }

  if (normalize_counts) {
    x <- x / sum(x, na.rm = TRUE)
    y <- y / sum(y, na.rm = TRUE)
  }

  result[non_zero] <- 2 * ( ( x[non_zero] - y[non_zero] ) /
                              abs(x[non_zero] + y[non_zero]) )

  median(abs(result))
}


## Function taken from mamabear - modified MRD1
relative_difference_mod1 <- function(x, y, na_zeroes = FALSE, normalize_counts = TRUE) {
  stopifnot(length(x) == length(y))
  result <- rep(NA_real_, length(x))

  non_zero <- which( x > 0 | y > 0 )
  both_zero <- setdiff(seq_along(x), non_zero)


  if (!na_zeroes) {
    result[both_zero] <- 0.0
  }

  if (normalize_counts) {
    x <- x / sum(x, na.rm = TRUE)
    y <- y / sum(y, na.rm = TRUE)
  }
  result[non_zero] <- 2 * ( abs( x[non_zero] - y[non_zero] ) /
                              (x[non_zero] + y[non_zero]) )

  median(abs(result[non_zero]))
}



## Function taken from mamabear - modified MRD2
relative_difference_mod2 <- function(x, y, na_zeroes = FALSE, normalize_counts = TRUE) {
  stopifnot(length(x) == length(y))
  result <- rep(NA_real_, length(x))

  non_zero <- which( x > 0 | y > 0 )
  both_zero <- setdiff(seq_along(x), non_zero)
  above_threshold <- which( x > 50 | y > 50 )


  if (!na_zeroes) {
    result[both_zero] <- 0.0
  }

  if (normalize_counts) {
    x <- x / sum(x, na.rm = TRUE)
    y <- y / sum(y, na.rm = TRUE)
  }
  result[non_zero] <- ( abs( x[non_zero] - y[non_zero] ) /
                          y[non_zero] )

  median(abs(result))
}




correlationPlotter <- function(groundTruth,modKallisto,origKallisto,dataset,geneId_transId_file,outputPath,groundTruthType,scale=1){


  if(groundTruthType=="RSEM"){
    ## Abundances for the ground truth
    groundTruth <- read.table(groundTruth,header = TRUE)
    groundTruth$transcript_id <- gsub("_.*", "", groundTruth$transcript_id)
    groundTruth <- groundTruth[c(1,5,6)]
    names(groundTruth) <- c("target_id","est_counts_ground","tpm_ground")
  }else if(groundTruthType=="KALLISTO"){
    ## Abundances for the ground truth
    groundTruth <- read.csv(groundTruth,header = TRUE)
    groundTruth <- groundTruth[c(2,4,5)]
    names(groundTruth) <- c("target_id","est_counts_ground","tpm_ground")
  }

  ## reading in the original kallisto abundance file here to calculate the scale
  abundanceSample1 <- read.table(origKal, header= T)
  # also reading the modified kallisto file and filtering lowly expressed transcripts
  

  scale <- sum(abundanceSample1$est_counts) / sum(groundTruth$est_counts_ground)
  groundTruth$est_counts_ground <- groundTruth$est_counts_ground*scale
  #groundTruth <- groundTruth %>% subset(groundTruth$est_counts_ground > 0)
  groundTruth$est_counts_ground_log <- (groundTruth$est_counts_ground) + 5
  groundTruth$est_counts_ground_log <- log2(groundTruth$est_counts_ground_log)
  groundTruth$tpm_ground <- groundTruth$tpm_ground + 0.1
  groundTruth$tpm_ground <- log2(groundTruth$tpm_ground)

  #groundTruthGeneId <- merge(groundTruth,geneId_transId,by = "target_id") 
  #geneComplexity1 <- groundTruthGeneId %>% group_by(ensembl_gene_id) %>%
  #  summarise(Complexity=n())
  
  # Reading the modified kallisto file in order to create good transcript list
  abundanceSample2 <- read.table(modKal, header= T)
  mergedTrans <- merge(abundanceSample1,abundanceSample2, by = "target_id")
  # keeping transcripts which are not estimated as 0 by both methods
  mergedTrans <- mergedTrans %>% subset(!(est_counts.x == 0 & est_counts.y == 0))
  mergedTrans <- mergedTrans %>% subset(mergedTrans$target_id %in% groundTruth$target_id) 
  goodTransList <- as.data.frame(mergedTrans$target_id)
  print(length(goodTransList))
  names(goodTransList) <- c("target_id")
  
  
  ## Reading the geneid transcriptid mappings
  geneId_transId <- read.table(file = geneId_transId_file, header = T)
  names(geneId_transId)[2] <- "target_id"
  
  ## Making the genecomplexities here
  groundTruthGeneId <- merge(groundTruth,geneId_transId,by = "target_id") 
  geneComplexity <- groundTruthGeneId %>% group_by(ensembl_gene_id) %>%
    summarise(Complexity=n())


  ## Abundances for the original kallisto sample
  abundanceSample1 <- merge(abundanceSample1,geneId_transId,by = "target_id")
  abundanceSample1["length"] <- abundanceSample1["length"] + 200
  #abundanceSample1 <- abundanceSample1 %>% subset(abundanceSample1$tpm > 0)
  abundanceSample1$tpm <- abundanceSample1$tpm + 0.1
  abundanceSample1$tpm <- log2(abundanceSample1$tpm)
  abundanceSample1$est_counts_log <- abundanceSample1$est_counts + 5
  abundanceSample1$est_counts_log <- log2(abundanceSample1$est_counts_log)



  #### Defining the complexity and effective complexity for original kallisto files

  # Assigning actual band where the transcript should have been ###
  abundanceSample1["ActualBand"] <- abundanceSample1$length %>%
    cut(c(0,1000,1500,2000,3000,4000,6000,10000),labels=FALSE)


  mergedSample1 <- merge(abundanceSample1,groundTruth,by = "target_id")
  mergedSample1 <- mergedSample1[c(1,4,5,6,7,8,9,10,11)]

  #geneComplexity1 <- mergedSample1 %>% group_by(ensembl_gene_id) %>%
  #  summarise(Complexity=n())

  mergedSample1 <- merge(mergedSample1,geneComplexity,by="ensembl_gene_id")

  mergedSample1$Complexity <- mergedSample1$Complexity  %>% replace(., mergedSample1$Complexity>=11, 11)

  effectiveComplexityList1 <- mergedSample1 %>%
    group_by(ensembl_gene_id,Complexity,ActualBand) %>%
    summarise(Effective_complexity=n())

  mergedSample1 <- merge(mergedSample1,effectiveComplexityList1,by=c("ensembl_gene_id","Complexity","ActualBand"))
  mergedSample1 <- mergedSample1 %>% subset(mergedSample1$target_id %in% goodTransList$target_id)
  
  

  
  ###### grouping by both gene and effective complexity ######
  ## Est_counts correlation
  coRRFileCOUNTSample1 <- mergedSample1 %>% group_by(Complexity,Effective_complexity) %>%
    summarise(COR=cor(est_counts_log,est_counts_ground_log, method = "spearman"),numOfTranscripts=n())
  coRRFileCOUNTSample1["Group"] <- "Conventional"
  write.table(coRRFileCOUNTSample1, file = paste(outputPath,"geneEffectiveComp/corr_Orig_File.txt", sep = ""), quote = F, row.names = F, col.names = T)
  ## TPM correlation
  coRRFileTPMSample1 <- mergedSample1 %>% group_by(Complexity,Effective_complexity) %>%
    summarise(COR=cor(tpm,tpm_ground),numOfTranscripts=n())
  coRRFileTPMSample1["Group"] <- "Conventional"
  write.table(coRRFileTPMSample1, file = paste(outputPath,"geneEffectiveComp/tpm_Orig_File.txt", sep = ""), quote = F, row.names = F, col.names = T)
  ## Est_counts MARDS
  mardsFileCountSample1 <-  mergedSample1 %>% group_by(Complexity,Effective_complexity) %>%
    summarise(MARDS=mards(est_counts,est_counts_ground),numOfTranscripts=n())
  mardsFileCountSample1["Group"] <- "Conventional"
  write.table(mardsFileCountSample1, file = paste(outputPath,"geneEffectiveComp/MARDS_Orig_File.txt", sep = ""), quote = F, row.names = F, col.names = T)
  ## Est_counts medianRelativeDiff Conventional kallisto
  medianRelDiffCountsFileSample1 <-  mergedSample1 %>% group_by(Complexity,Effective_complexity) %>%
    summarise(mrd=relative_difference(est_counts_log,est_counts_ground_log),numOfTranscripts=n())
  medianRelDiffCountsFileSample1["Group"] <- "Conventional"
  write.table(medianRelDiffCountsFileSample1, file = paste(outputPath,"geneEffectiveComp/Mrd_Orig_File.txt", sep = ""), quote = F, row.names = F, col.names = T)
  ## Est_counts medianRelativeDiff modified 1
  medianRelDiffCountsFileMod1Sample1 <-  mergedSample1 %>% group_by(Complexity,Effective_complexity) %>%
    summarise(mrd1=relative_difference_mod1(est_counts,est_counts_ground),numOfTranscripts=n())
  medianRelDiffCountsFileMod1Sample1["Group"] <- "Conventional"
  write.table(medianRelDiffCountsFileMod1Sample1, file = paste(outputPath,"geneEffectiveComp/Mrd_mod1_Orig_File.txt", sep = ""), quote = F, row.names = F, col.names = T)
  ## Est_counts medianRelativeDiff modified 2
  medianRelDiffCountsFileMod2Sample1 <-  mergedSample1 %>% group_by(Complexity,Effective_complexity) %>%
    summarise(mrd2=relative_difference_mod2(est_counts,est_counts_ground),numOfTranscripts=n())
  medianRelDiffCountsFileMod2Sample1["Group"] <- "Conventional"
  write.table(medianRelDiffCountsFileMod2Sample1, file = paste(outputPath,"geneEffectiveComp/Mrd_mod2_Orig_File.txt", sep = ""), quote = F, row.names = F, col.names = T)


  ###### grouping by gene complexity without Complexity==EffectiveComplexity######
  ## Est_counts correlation
  geneCompCOUNTSample1 <- mergedSample1 %>%  subset(Complexity!=Effective_complexity) %>% group_by(Complexity) %>%
    summarise(COR=cor(est_counts_log,est_counts_ground_log),numOfTranscripts=n())
  geneCompCOUNTSample1["Group"] <- "Conventional"
  write.table(geneCompCOUNTSample1, file = paste(outputPath,"geneComp/corr_Orig_File.txt", sep = ""), quote = F, row.names = F, col.names = T)
  ## TPM correlation
  geneCompTPMSample1 <- mergedSample1 %>%  subset(Complexity!=Effective_complexity) %>% group_by(Complexity) %>%
    summarise(COR=cor(tpm,tpm_ground),numOfTranscripts=n())
  geneCompTPMSample1["Group"] <- "Conventional"
  write.table(geneCompTPMSample1, file = paste(outputPath,"geneComp/tpm_Orig_File.txt", sep = ""), quote = F, row.names = F, col.names = T)
  ## Est_counts MARDS
  geneCompMardsSample1 <-  mergedSample1 %>%  subset(Complexity!=Effective_complexity) %>% group_by(Complexity) %>%
    summarise(MARDS=mards(est_counts,est_counts_ground),numOfTranscripts=n())
  geneCompMardsSample1["Group"] <- "Conventional"
  write.table(geneCompMardsSample1, file = paste(outputPath,"geneComp/MARDS_Orig_File.txt", sep = ""), quote = F, row.names = F, col.names = T)
  ## Est_counts medianRelativeDiff Conventional kallisto
  geneCompMedianRelDiffSample1 <-  mergedSample1 %>%  subset(Complexity!=Effective_complexity) %>% group_by(Complexity) %>%
    summarise(mrd=relative_difference(est_counts,est_counts_ground),numOfTranscripts=n())
  geneCompMedianRelDiffSample1["Group"] <- "Conventional"
  write.table(geneCompMedianRelDiffSample1, file = paste(outputPath,"geneComp/Mrd_Orig_File.txt", sep = ""), quote = F, row.names = F, col.names = T)
  ## Est_counts medianRelativeDiff modified 1
  geneCompMedianRelDiffMod1Sample1 <-  mergedSample1 %>%  subset(Complexity!=Effective_complexity) %>% group_by(Complexity) %>%
    summarise(mrd1=relative_difference_mod1(est_counts,est_counts_ground),numOfTranscripts=n())
  geneCompMedianRelDiffMod1Sample1["Group"] <- "Conventional"
  write.table(geneCompMedianRelDiffMod1Sample1, file = paste(outputPath,"geneComp/Mrd_mod1_Orig_File.txt", sep = ""), quote = F, row.names = F, col.names = T)
  ## Est_counts medianRelativeDiff modified 2
  geneCompMedianRelDiffMod2Sample1 <-  mergedSample1 %>%  subset(Complexity!=Effective_complexity) %>% group_by(Complexity) %>%
    summarise(mrd2=relative_difference_mod2(est_counts,est_counts_ground),numOfTranscripts=n())
  geneCompMedianRelDiffMod2Sample1["Group"] <- "Conventional"
  write.table(geneCompMedianRelDiffMod2Sample1, file = paste(outputPath,"geneComp/Mrd_mod2_Orig_File.txt", sep = ""), quote = F, row.names = F, col.names = T)


  ###### grouping by gene complexity with Complexity==EffectiveComplexity######
  ## Est_counts correlation
  geneCompCOUNTSample1 <- mergedSample1 %>% group_by(Complexity) %>%
    summarise(COR=cor(est_counts_log,est_counts_ground_log),numOfTranscripts=n())
  geneCompCOUNTSample1["Group"] <- "Conventional"
  write.table(geneCompCOUNTSample1, file = paste(outputPath,"geneComp/allIncluded_corr_Orig_File.txt", sep = ""), quote = F, row.names = F, col.names = T)
  ## TPM correlation
  geneCompTPMSample1 <- mergedSample1 %>% group_by(Complexity) %>%
    summarise(COR=cor(tpm,tpm_ground),numOfTranscripts=n())
  geneCompTPMSample1["Group"] <- "Conventional"
  write.table(geneCompTPMSample1, file = paste(outputPath,"geneComp/allIncluded_tpm_Orig_File.txt", sep = ""), quote = F, row.names = F, col.names = T)
  ## Est_counts MARDS
  geneCompMardsSample1 <-  mergedSample1 %>% group_by(Complexity) %>%
    summarise(MARDS=mards(est_counts,est_counts_ground),numOfTranscripts=n())
  geneCompMardsSample1["Group"] <- "Conventional"
  write.table(geneCompMardsSample1, file = paste(outputPath,"geneComp/allIncluded_MARDS_Orig_File.txt", sep = ""), quote = F, row.names = F, col.names = T)
  ## Est_counts medianRelativeDiff Conventional kallisto
  geneCompMedianRelDiffSample1 <-  mergedSample1 %>% group_by(Complexity) %>%
    summarise(mrd=relative_difference(est_counts,est_counts_ground),numOfTranscripts=n())
  geneCompMedianRelDiffSample1["Group"] <- "Conventional"
  write.table(geneCompMedianRelDiffSample1, file = paste(outputPath,"geneComp/allIncluded_Mrd_Orig_File.txt", sep = ""), quote = F, row.names = F, col.names = T)
  ## Est_counts medianRelativeDiff modified 1
  geneCompMedianRelDiffMod1Sample1 <-  mergedSample1 %>% group_by(Complexity) %>%
    summarise(mrd1=relative_difference_mod1(est_counts,est_counts_ground),numOfTranscripts=n())
  geneCompMedianRelDiffMod1Sample1["Group"] <- "Conventional"
  write.table(geneCompMedianRelDiffMod1Sample1, file = paste(outputPath,"geneComp/allIncluded_Mrd_mod1_Orig_File.txt", sep = ""), quote = F, row.names = F, col.names = T)
  ## Est_counts medianRelativeDiff modified 2
  geneCompMedianRelDiffMod2Sample1 <-  mergedSample1 %>% group_by(Complexity) %>%
    summarise(mrd2=relative_difference_mod2(est_counts,est_counts_ground),numOfTranscripts=n())
  geneCompMedianRelDiffMod2Sample1["Group"] <- "Conventional"
  write.table(geneCompMedianRelDiffMod2Sample1, file = paste(outputPath,"geneComp/allIncluded_Mrd_mod2_Orig_File.txt", sep = ""), quote = F, row.names = F, col.names = T)





  ###### grouping by difference in gene and effective_gene complexity ######
  mergedSample1["diffGene_EffectiveGene"] <- mergedSample1$Complexity - mergedSample1$Effective_complexity

  ## Est_counts correlation
  differenceCOUNTSample1 <- mergedSample1 %>% group_by(diffGene_EffectiveGene) %>%
    summarise(COR=cor(est_counts_log,est_counts_ground_log),numOfTranscripts=n())
  differenceCOUNTSample1["Group"] <- "Conventional"
  write.table(differenceCOUNTSample1, file = paste(outputPath,"difference/corr_Orig_File.txt", sep = ""), quote = F, row.names = F, col.names = T)
  ## TPM correlation
  differenceTPMSample1 <- mergedSample1 %>% group_by(diffGene_EffectiveGene) %>%
    summarise(COR=cor(tpm,tpm_ground),numOfTranscripts=n())
  differenceTPMSample1["Group"] <- "Conventional"
  write.table(differenceTPMSample1, file = paste(outputPath,"difference/tpm_Orig_File.txt", sep = ""), quote = F, row.names = F, col.names = T)
  ## Est_counts MARDS
  differenceMardsSample1 <-  mergedSample1 %>% group_by(diffGene_EffectiveGene) %>%
    summarise(MARDS=mards(est_counts,est_counts_ground),numOfTranscripts=n())
  differenceMardsSample1["Group"] <- "Conventional"
  write.table(differenceMardsSample1, file = paste(outputPath,"difference/MARDS_Orig_File.txt", sep = ""), quote = F, row.names = F, col.names = T)
  ## Est_counts medianRelativeDiff Conventional kallisto
  differenceMedianRelDiffSample1 <-  mergedSample1 %>% group_by(diffGene_EffectiveGene) %>%
    summarise(mrd=relative_difference(est_counts,est_counts_ground),numOfTranscripts=n())
  differenceMedianRelDiffSample1["Group"] <- "Conventional"
  write.table(differenceMedianRelDiffSample1, file = paste(outputPath,"difference/Mrd_Orig_File.txt", sep = ""), quote = F, row.names = F, col.names = T)
  ## Est_counts medianRelativeDiff modified 1
  differenceMedianRelDiffMod1Sample1 <-  mergedSample1 %>% group_by(diffGene_EffectiveGene) %>%
    summarise(mrd1=relative_difference_mod1(est_counts,est_counts_ground),numOfTranscripts=n())
  differenceMedianRelDiffMod1Sample1["Group"] <- "Conventional"
  write.table(differenceMedianRelDiffMod1Sample1, file = paste(outputPath,"difference/Mrd_mod1_Orig_File.txt", sep = ""), quote = F, row.names = F, col.names = T)
  ## Est_counts medianRelativeDiff modified 2
  differenceMedianRelDiffMod2Sample1 <-  mergedSample1 %>% group_by(diffGene_EffectiveGene) %>%
    summarise(mrd2=relative_difference_mod2(est_counts,est_counts_ground),numOfTranscripts=n())
  differenceMedianRelDiffMod2Sample1["Group"] <- "Conventional"
  write.table(differenceMedianRelDiffMod2Sample1, file = paste(outputPath,"difference/Mrd_mod2_Orig_File.txt", sep = ""), quote = F, row.names = F, col.names = T)













  ## Abundances for the modified kallisto sample
  #abundanceSample2 <- read.table(modKal, header= T)
  abundanceSample2 <- merge(abundanceSample2,geneId_transId,by = "target_id")
  abundanceSample2["length"] <- abundanceSample2["length"] + 200
  #abundanceSample2 <- abundanceSample2 %>% subset(abundanceSample2$tpm > 0)
  abundanceSample2$tpm <- abundanceSample2$tpm + 0.1
  abundanceSample2$tpm <- log2(abundanceSample2$tpm)
  abundanceSample2$est_counts_log <- abundanceSample2$est_counts + 5
  abundanceSample2$est_counts_log <- log2(abundanceSample2$est_counts_log)



  #### Defining the complexity and effective complexity for Conventional kallisto files

  # Assigning actual band where the transcript should have been ###
  abundanceSample2["ActualBand"] <- abundanceSample2$length %>%
    cut(c(0,1000,1500,2000,3000,4000,6000,10000),labels=FALSE)


  mergedSample2 <- merge(abundanceSample2,groundTruth,by = "target_id")
  mergedSample2 <- mergedSample2[c(1,4,5,6,7,8,9,10,11)]

  #geneComplexity2 <- mergedSample2 %>% group_by(ensembl_gene_id) %>%
  #  summarise(Complexity=n())

  mergedSample2 <- merge(mergedSample2,geneComplexity,by="ensembl_gene_id") 

  mergedSample2$Complexity <- mergedSample2$Complexity  %>% replace(., mergedSample2$Complexity>=11, 11)

  effectiveComplexityList2 <- mergedSample2 %>%
    group_by(ensembl_gene_id,Complexity,ActualBand) %>%
    summarise(Effective_complexity=n())


  mergedSample2 <- merge(mergedSample2,effectiveComplexityList2,by=c("ensembl_gene_id","Complexity","ActualBand"))
  mergedSample2 <- mergedSample2 %>% subset(mergedSample2$target_id %in% goodTransList$target_id)
  
  

  ###### grouping by both gene and effective complexity ######
  ## Est_counts correlation
  coRRFileCOUNTSample2 <- mergedSample2 %>% group_by(Complexity,Effective_complexity) %>%
    summarise(COR=cor(est_counts_log,est_counts_ground_log),numOfTranscripts=n())
  coRRFileCOUNTSample2["Group"] <- "Ladder-Seq"
  write.table(coRRFileCOUNTSample2, file = paste(outputPath,"geneEffectiveComp/corr_Mod_File.txt", sep = ""), quote = F, row.names = F, col.names = T)
  ## TPM correlation
  coRRFileTPMSample2 <- mergedSample2 %>% group_by(Complexity,Effective_complexity) %>%
    summarise(COR=cor(tpm,tpm_ground),numOfTranscripts=n())
  coRRFileTPMSample2["Group"] <- "Ladder-Seq"
  write.table(coRRFileTPMSample2, file = paste(outputPath,"geneEffectiveComp/tpm_Mod_File.txt", sep = ""), quote = F, row.names = F, col.names = T)
  ## Est_counts MARDS
  mardsFileCountSample2 <-  mergedSample2 %>% group_by(Complexity,Effective_complexity) %>%
    summarise(MARDS=mards(est_counts,est_counts_ground),numOfTranscripts=n())
  mardsFileCountSample2["Group"] <- "Ladder-Seq"
  write.table(mardsFileCountSample2, file = paste(outputPath,"geneEffectiveComp/MARDS_Mod_File.txt", sep = ""), quote = F, row.names = F, col.names = T)
  ## Est_counts medianRelDiff Ladder-Seq kallisto
  medianRelDiffCountsFileSample2 <-  mergedSample2 %>% group_by(Complexity,Effective_complexity) %>%
    summarise(mrd=relative_difference(est_counts,est_counts_ground),numOfTranscripts=n())
  medianRelDiffCountsFileSample2["Group"] <- "Ladder-Seq"
  write.table(medianRelDiffCountsFileSample2, file = paste(outputPath,"geneEffectiveComp/Mrd_Mod_File.txt", sep = ""), quote = F, row.names = F, col.names = T)
  ## Est_counts medianRelDiff Ladder-Seq 1
  medianRelDiffCountsFileMod1Sample2 <-  mergedSample2 %>% group_by(Complexity,Effective_complexity) %>%
    summarise(mrd1=relative_difference_mod2(est_counts,est_counts_ground),numOfTranscripts=n())
  medianRelDiffCountsFileMod1Sample2["Group"] <- "Ladder-Seq"
  write.table(medianRelDiffCountsFileMod1Sample2, file = paste(outputPath,"geneEffectiveComp/Mrd_mod1_Mod_File.txt", sep = ""), quote = F, row.names = F, col.names = T)
  ## Est_counts medianRelDiff Ladder-Seq 2
  medianRelDiffCountsFileMod2Sample2 <-  mergedSample2 %>% group_by(Complexity,Effective_complexity) %>%
    summarise(mrd2=relative_difference_mod2(est_counts,est_counts_ground),numOfTranscripts=n())
  medianRelDiffCountsFileMod2Sample2["Group"] <- "Ladder-Seq"
  write.table(medianRelDiffCountsFileMod2Sample2, file = paste(outputPath,"geneEffectiveComp/Mrd_mod2_Mod_File.txt", sep = ""), quote = F, row.names = F, col.names = T)







  ###### grouping by gene complexity without Complexity==Effective_complexity######
  ## Est_counts correlation
  geneCompCOUNTSample2 <- mergedSample2 %>%  subset(Complexity!=Effective_complexity) %>% group_by(Complexity) %>%
    summarise(COR=cor(est_counts_log,est_counts_ground_log),numOfTranscripts=n())
  geneCompCOUNTSample2["Group"] <- "Ladder-Seq"
  write.table(geneCompCOUNTSample2, file = paste(outputPath,"geneComp/corr_Mod_File.txt", sep = ""), quote = F, row.names = F, col.names = T)
  ## TPM correlation
  geneCompTPMSample2 <- mergedSample2 %>%  subset(Complexity!=Effective_complexity) %>% group_by(Complexity) %>%
    summarise(COR=cor(tpm,tpm_ground),numOfTranscripts=n())
  geneCompTPMSample2["Group"] <- "Ladder-Seq"
  write.table(geneCompTPMSample2, file = paste(outputPath,"geneComp/tpm_Mod_File.txt", sep = ""), quote = F, row.names = F, col.names = T)
  ## Est_counts MARDS
  geneCompMardsSample2 <-  mergedSample2 %>%  subset(Complexity!=Effective_complexity) %>% group_by(Complexity) %>%
    summarise(MARDS=mards(est_counts,est_counts_ground),numOfTranscripts=n())
  geneCompMardsSample2["Group"] <- "Ladder-Seq"
  write.table(geneCompMardsSample2, file = paste(outputPath,"geneComp/MARDS_Mod_File.txt", sep = ""), quote = F, row.names = F, col.names = T)
  ## Est_counts medianRelativeDiff original kallisto
  geneCompMedianRelDiffSample2 <-  mergedSample2 %>%  subset(Complexity!=Effective_complexity) %>% group_by(Complexity) %>%
    summarise(mrd=relative_difference(est_counts,est_counts_ground),numOfTranscripts=n())
  geneCompMedianRelDiffSample2["Group"] <- "Ladder-Seq"
  write.table(geneCompMedianRelDiffSample2, file = paste(outputPath,"geneComp/Mrd_Mod_File.txt", sep = ""), quote = F, row.names = F, col.names = T)
  ## Est_counts medianRelativeDiff Ladder-Seq 1
  geneCompMedianRelDiffMod1Sample2 <-  mergedSample2 %>%  subset(Complexity!=Effective_complexity) %>% group_by(Complexity) %>%
    summarise(mrd1=relative_difference_mod1(est_counts,est_counts_ground),numOfTranscripts=n())
  geneCompMedianRelDiffMod1Sample2["Group"] <- "Ladder-Seq"
  write.table(geneCompMedianRelDiffMod1Sample2, file = paste(outputPath,"geneComp/Mrd_mod1_Mod_File.txt", sep = ""), quote = F, row.names = F, col.names = T)
  ## Est_counts medianRelativeDiff Ladder-Seq 2
  geneCompMedianRelDiffMod2Sample2 <-  mergedSample2 %>%  subset(Complexity!=Effective_complexity) %>% group_by(Complexity) %>%
    summarise(mrd2=relative_difference_mod2(est_counts,est_counts_ground),numOfTranscripts=n())
  geneCompMedianRelDiffMod2Sample2["Group"] <- "Ladder-Seq"
  write.table(geneCompMedianRelDiffMod2Sample2, file = paste(outputPath,"geneComp/Mrd_mod2_Mod_File.txt", sep = ""), quote = F, row.names = F, col.names = T)

  ###### grouping by gene complexity with Complexity==Effective_complexity######
  ## Est_counts correlation
  geneCompCOUNTSample2 <- mergedSample2 %>% group_by(Complexity) %>%
    summarise(COR=cor(est_counts_log,est_counts_ground_log),numOfTranscripts=n())
  geneCompCOUNTSample2["Group"] <- "Ladder-Seq"
  write.table(geneCompCOUNTSample2, file = paste(outputPath,"geneComp/allIncluded_corr_Mod_File.txt", sep = ""), quote = F, row.names = F, col.names = T)
  ## TPM correlation
  geneCompTPMSample2 <- mergedSample2 %>% group_by(Complexity) %>%
    summarise(COR=cor(tpm,tpm_ground),numOfTranscripts=n())
  geneCompTPMSample2["Group"] <- "Ladder-Seq"
  write.table(geneCompTPMSample2, file = paste(outputPath,"geneComp/allIncluded_tpm_Mod_File.txt", sep = ""), quote = F, row.names = F, col.names = T)
  ## Est_counts MARDS
  geneCompMardsSample2 <-  mergedSample2 %>% group_by(Complexity) %>%
    summarise(MARDS=mards(est_counts,est_counts_ground),numOfTranscripts=n())
  geneCompMardsSample2["Group"] <- "Ladder-Seq"
  write.table(geneCompMardsSample2, file = paste(outputPath,"geneComp/allIncluded_MARDS_Mod_File.txt", sep = ""), quote = F, row.names = F, col.names = T)
  ## Est_counts medianRelativeDiff original kallisto
  geneCompMedianRelDiffSample2 <-  mergedSample2 %>% group_by(Complexity) %>%
    summarise(mrd=relative_difference(est_counts,est_counts_ground),numOfTranscripts=n())
  geneCompMedianRelDiffSample2["Group"] <- "Ladder-Seq"
  write.table(geneCompMedianRelDiffSample2, file = paste(outputPath,"geneComp/allIncluded_Mrd_Mod_File.txt", sep = ""), quote = F, row.names = F, col.names = T)
  ## Est_counts medianRelativeDiff Ladder-Seq 1
  geneCompMedianRelDiffMod1Sample2 <-  mergedSample2 %>% group_by(Complexity) %>%
    summarise(mrd1=relative_difference_mod1(est_counts,est_counts_ground),numOfTranscripts=n())
  geneCompMedianRelDiffMod1Sample2["Group"] <- "Ladder-Seq"
  write.table(geneCompMedianRelDiffMod1Sample2, file = paste(outputPath,"geneComp/allIncluded_Mrd_mod1_Mod_File.txt", sep = ""), quote = F, row.names = F, col.names = T)
  ## Est_counts medianRelativeDiff Ladder-Seq 2
  geneCompMedianRelDiffMod2Sample2 <-  mergedSample2 %>% group_by(Complexity) %>%
    summarise(mrd2=relative_difference_mod2(est_counts,est_counts_ground),numOfTranscripts=n())
  geneCompMedianRelDiffMod2Sample2["Group"] <- "Ladder-Seq"
  write.table(geneCompMedianRelDiffMod2Sample2, file = paste(outputPath,"geneComp/allIncluded_Mrd_mod2_Mod_File.txt", sep = ""), quote = F, row.names = F, col.names = T)







  ###### grouping by difference in gene and effective_gene complexity ######
  mergedSample2["diffGene_EffectiveGene"] <- mergedSample2$Complexity - mergedSample2$Effective_complexity

  ## Est_counts correlation
  differenceCOUNTSample2 <- mergedSample2 %>% group_by(diffGene_EffectiveGene) %>%
    summarise(COR=cor(est_counts_log,est_counts_ground_log),numOfTranscripts=n())
  differenceCOUNTSample2["Group"] <- "Ladder-Seq"
  write.table(differenceCOUNTSample2, file = paste(outputPath,"difference/corr_Mod_File.txt", sep = ""), quote = F, row.names = F, col.names = T)
  ## TPM correlation
  differenceTPMSample2 <- mergedSample2 %>% group_by(diffGene_EffectiveGene) %>%
    summarise(COR=cor(tpm,tpm_ground),numOfTranscripts=n())
  differenceTPMSample2["Group"] <- "Ladder-Seq"
  write.table(differenceTPMSample2, file = paste(outputPath,"difference/tpm_Mod_File.txt", sep = ""), quote = F, row.names = F, col.names = T)
  ## Est_counts MARDS
  differenceMardsSample2 <-  mergedSample2 %>% group_by(diffGene_EffectiveGene) %>%
    summarise(MARDS=mards(est_counts,est_counts_ground),numOfTranscripts=n())
  differenceMardsSample2["Group"] <- "Ladder-Seq"
  write.table(differenceMardsSample2, file = paste(outputPath,"difference/MARDS_Mod_File.txt", sep = ""), quote = F, row.names = F, col.names = T)
  ## Est_counts medianRelativeDiff original kallisto
  differenceMedianRelDiffSample2 <-  mergedSample2 %>% group_by(diffGene_EffectiveGene) %>%
    summarise(mrd=relative_difference(est_counts,est_counts_ground),numOfTranscripts=n())
  differenceMedianRelDiffSample2["Group"] <- "Ladder-Seq"
  write.table(differenceMedianRelDiffSample2, file = paste(outputPath,"difference/Mrd_Mod_File.txt", sep = ""), quote = F, row.names = F, col.names = T)
  ## Est_counts medianRelativeDiff Ladder-Seq 1
  differenceMedianRelDiffMod1Sample2 <-  mergedSample2 %>% group_by(diffGene_EffectiveGene) %>%
    summarise(mrd1=relative_difference_mod1(est_counts,est_counts_ground),numOfTranscripts=n())
  differenceMedianRelDiffMod1Sample2["Group"] <- "Ladder-Seq"
  write.table(differenceMedianRelDiffMod1Sample2, file = paste(outputPath,"difference/Mrd_mod1_Mod_File.txt", sep = ""), quote = F, row.names = F, col.names = T)
  ## Est_counts medianRelativeDiff Ladder-Seq 2
  differenceMedianRelDiffMod2Sample2 <-  mergedSample2 %>% group_by(diffGene_EffectiveGene) %>%
    summarise(mrd2=relative_difference_mod2(est_counts,est_counts_ground),numOfTranscripts=n())
  differenceMedianRelDiffMod2Sample2["Group"] <- "Ladder-Seq"
  write.table(differenceMedianRelDiffMod2Sample2, file = paste(outputPath,"difference/Mrd_mod2_Mod_File.txt", sep = ""), quote = F, row.names = F, col.names = T)




}




correlationPlotter(groundTruth,modKal,origKal,dataset,geneId_transId_file,outputPath,groundTruthType,scale)
