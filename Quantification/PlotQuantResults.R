


# Libraries
require(ggplot2)
require(dplyr)


# ### Command line arguments ###
args <- commandArgs(trailingOnly = TRUE)



## hardcoded arguments
MARDSFile1 <- "MARDS_Orig_File_acrossBetas.txt"
MARDSFile2 <- "MARDS_Mod_File_acrossBetas.txt"
tpmCORRFile1 <- "tpm_Orig_File_acrossBetas.txt"
tpmCORRFile2 <- "tpm_Mod_File_acrossBetas.txt"



aspect_ratio <- 1.5
height <- 7
filePath <- args[1]
dataset <- args[2]


## defining colours
ladderColour1 <- rgb(222, 80, 0, 193, names = NULL, maxColorValue = 255)
ladderColour2 <- rgb(235, 40, 0, 90, names = NULL, maxColorValue = 255)
conventionalColour <- rgb(44, 118, 212, 255, names = NULL, maxColorValue = 255)






# this is without the cases where geneComplexity = effective complexity
plotByJustGeneComplexity <- function(countCORRFile1,
                                  countCORRFile2,
                                  tpmCORRFile1,
                                  tpmCORRFile2,
                                  MARDSFile1,
                                  MARDSFile2,
                                  medianRelDiffFile1,
                                  medianRelDiffFile2,
                                  medianRelDiffMod1File1,
                                  medianRelDiffMod1File2,
                                  medianRelDiffMod2File1,
                                  medianRelDiffMod2File2,
                                  filePath,
                                  dataset){

  specificfilePath <- paste(filePath,"geneComp/allIncluded_",sep = "")




  tpmCORRSample1_Unsplitted <- read.table(file = paste(specificfilePath,tpmCORRFile1,sep = ""), header = T)
  tpmCORRSample2_Unsplitted <- read.table(file = paste(specificfilePath,tpmCORRFile2,sep = ""), header = T)


  MARDSSample1_Unsplitted <- read.table(file = paste(specificfilePath,MARDSFile1,sep = ""), header = T)
  MARDSSample2_Unsplitted <- read.table(file = paste(specificfilePath,MARDSFile2,sep = ""), header = T)



  #### plotting tpms ####

  tpmCORRSample2_Unsplitted <- tpmCORRSample2_Unsplitted %>% subset(Complexity!="Complexity")
  tpmCORRSample2_Unsplitted$COR <- as.numeric(as.character(tpmCORRSample2_Unsplitted$COR))
  coRRFileTPMGeneComplex <- rbind(tpmCORRSample2_Unsplitted,tpmCORRSample1_Unsplitted)
  coRRFileTPMGeneComplex$Complexity <- as.numeric(as.character(coRRFileTPMGeneComplex$Complexity))
  coRRFileTPMGeneComplex <- coRRFileTPMGeneComplex %>% subset(Complexity>1 & Complexity<=10)
  coRRFileTPMGeneComplex$Group <- paste(coRRFileTPMGeneComplex$Group,coRRFileTPMGeneComplex$bandProb,sep="_")

  g <- ggplot(coRRFileTPMGeneComplex,aes(x=Complexity,y=COR,colour=Group,group=Group)) +
    theme_classic() +
    theme(axis.title=element_text(size=16),
          legend.title=element_blank(),
          axis.text.x = element_text(size=14),
          axis.text.y = element_text(size=14)) +
    geom_line() +
    geom_point(size=4) +
    scale_x_continuous(breaks = c(2,3,4,5,6,7,8,9,10)) +
    scale_colour_manual("legend", values = c(conventionalColour,ladderColour2, ladderColour1))

  ggsave(device = "svg" ,g, height = 7 , width = 7 * aspect_ratio,filename = paste(specificfilePath,"geneCompTpm.svg",sep = ""))




  ### plotting MARDS

  MARDSSample2_Unsplitted <- MARDSSample2_Unsplitted %>% subset(Complexity!="Complexity")
  MARDSSample2_Unsplitted$MARDS <- as.numeric(as.character(MARDSSample2_Unsplitted$MARDS))
  coRRFileMardsGeneComplex <- rbind(MARDSSample2_Unsplitted,MARDSSample1_Unsplitted)
  coRRFileMardsGeneComplex$Complexity <- as.numeric(as.character(coRRFileMardsGeneComplex$Complexity))
  coRRFileMardsGeneComplex <- coRRFileMardsGeneComplex %>% subset(Complexity>1 & Complexity<=10)
  coRRFileMardsGeneComplex$Group <- paste(coRRFileMardsGeneComplex$Group,coRRFileMardsGeneComplex$bandProb,sep="_")

  g <- ggplot(coRRFileMardsGeneComplex,aes(x=Complexity,y=MARDS,colour=Group,group=Group)) +
    theme_classic() +
    theme(axis.title=element_text(size=16),
          legend.title=element_blank(),
          axis.text.x = element_text(size=14),
          axis.text.y = element_text(size=14)) +
    geom_line() +
    geom_point(size=4) +
    scale_x_continuous(breaks=c(2,3,4,5,6,7,8,9,10)) +
    scale_colour_manual("legend", values = c(conventionalColour,ladderColour2, ladderColour1))

  ggsave(device = "svg" ,g, height = 7 , width = 7 * aspect_ratio,filename = paste(specificfilePath,"geneCompMARDS.svg",sep = ""))


}



plotByJustGeneComplexity(countCORRFile1,
                        countCORRFile2,
                        tpmCORRFile1,
                        tpmCORRFile2,
                        MARDSFile1,
                        MARDSFile2,
                        medianRelDiffFile1,
                        medianRelDiffFile2,
                        medianRelDiffMod1File1,
                        medianRelDiffMod1File2,
                        medianRelDiffMod2File1,
                        medianRelDiffMod2File2,
                        filePath,
                        dataset)
