##### Create precision and recall plots from comparison results for stringtie assembly
##### Input
## 1. comparisonResultModFile : comparison results from modified (Ladder seq: Concat or stringtie merge) assemly runs.
## 2. comparisonResultOrigFile : comparison results from original RNA-seq assemly runs.
## 3. modificationType : Concat or Stringtie merge
## 4. dataset : type of reads used
## 5. comparisonTool : The tool used to compare the assembled transcripts to the ground truths. "gffCompare" or "gtfCompare"
## 6. pipeline: The pipeline used to assemble and analyse the data



## libraries
require(dplyr)
require(ggplot2)
require(reshape2)
require(reshape)
require(stringr)


#### Command line arguments
args = commandArgs(trailingOnly=TRUE)
comparisonResultOrigFile <- args[1]
comparisonResultModFile <- args[2]
outFilePath <- args[3]
dataset <- args[4]
comparisonTool <- args[5]
pipeline <- args[6]
aspect_ratio <- 1.5
height <- 7


# # #### Hardcoded arguments
# comparisonResultModFile <- "/algbio1/shounak/raw/simulated/human/Resimulation/30Million/simulations/sim_30000000/sim_1/gffSingleExonsRemovedAssemblyResults_acrossBetas.txt"
# comparisonResultOrigFile <- "/algbio1/shounak/raw/simulated/human/Resimulation/30Million/simulations/sim_30000000/sim_1/gffOriginalResults_acrossBetas.txt"
# outFilePath <- "/home/schakraborty/Documents/RProjects/RNADecon/stringtie/"
# dataset <- "RSEM Simulations,Imperfect30"
# comparisonTool <- "gffCompare"
# pipeline <- "Shounak'sPipeline"




plotFunction <- function(comparisonResultOrigFile, comparisonResultModFile, outFilePath, dataset, comparisonTool, pipeline){

  #defining colours
  ladderColour1 <- rgb(222, 80, 0, 193, names = NULL, maxColorValue = 255)
  ladderColour2 <- rgb(230, 60, 0, 130, names = NULL, maxColorValue = 255)
  ladderColour3 <- rgb(235, 40, 0, 90, names = NULL, maxColorValue = 255)
  ladderColour4 <- rgb(240, 20, 0, 50, names = NULL, maxColorValue = 255)
  ladderColour5 <- rgb(255, 0, 0, 30, names = NULL, maxColorValue = 255)
  conventionalColour <- rgb(44, 118, 212, 255, names = NULL, maxColorValue = 255)


  ## Read the data
  comparisonResultOriginal <- read.table(comparisonResultOrigFile, header = T)
  comparisonResultOriginal["Group"] <- "Conventional"
  comparisonResultOriginal$Group <- paste(comparisonResultOriginal$Group,comparisonResultOriginal$bandProb,sep="_")

  comparisonResultMod <- read.table(comparisonResultModFile, header = T)
  comparisonResultMod["Group"] <- "LadderSeq"
  comparisonResultMod <- comparisonResultMod %>% subset(Complexity!="Complexity")
  comparisonResultMod$Group <- paste(comparisonResultMod$Group,comparisonResultMod$bandProb,sep="_")



## Calculations for stringtie
  combinedStringtieCompFile <- rbind(comparisonResultOriginal, comparisonResultMod)

  combinedStringtieCompFile$Recall <- as.numeric(combinedStringtieCompFile$Recall)
  combinedStringtieCompFile$Precision <- as.numeric(combinedStringtieCompFile$Precision)
  combinedStringtieCompFile$Complexity <- as.numeric(combinedStringtieCompFile$Complexity)
  combinedStringtieCompFile["F1Score"] <- 2*(combinedStringtieCompFile$Precision*combinedStringtieCompFile$Recall)/(combinedStringtieCompFile$Precision+combinedStringtieCompFile$Recall)

  # Converting all complexities above 10 to 10
  combinedStringtieCompFile$Complexity <- combinedStringtieCompFile$Complexity  %>% replace(., combinedStringtieCompFile$Complexity>10, 10)



  ## Plotting stringtie results
  g <- ggplot(combinedStringtieCompFile, aes(x=Complexity,y=F1Score,colour=Group)) +
    theme_classic() +
    theme(axis.title=element_text(size=16),
          legend.title=element_blank(),
          axis.text.x = element_text(size=14),
          axis.text.y = element_text(size=14)) +
    geom_line() + geom_point(size=4) +
    scale_x_discrete(limits=c(1:10)) +
    scale_colour_manual("legend", values = c(conventionalColour,ladderColour2,ladderColour3,
                                             ladderColour4,ladderColour5,ladderColour1))
  ggsave(device = "svg" ,g, height = 7 , width = 7 * aspect_ratio,filename = paste(outFilePath,"Stringtie_F1_",comparisonTool,"_acrossBetas.svg",sep = ""))


  g <- ggplot(combinedStringtieCompFile, aes(x=Complexity,y=Precision,colour=Group)) +
    theme_classic() +
    theme(axis.title=element_text(size=16),
          legend.title=element_blank(),
          axis.text.x = element_text(size=14),
          axis.text.y = element_text(size=14)) +
    geom_line() + geom_point(size=4) +
    scale_x_discrete(limits=c(1:10)) +
    scale_colour_manual("legend", values = c(conventionalColour,ladderColour2,ladderColour3,
                                             ladderColour4,ladderColour5,ladderColour1))
  ggsave(device = "svg" ,g, height = 7 , width = 7 * aspect_ratio,filename = paste(outFilePath,"Stringtie_Pre_",comparisonTool,"_acrossBetas.svg",sep = ""))


  g <- ggplot(combinedStringtieCompFile, aes(x=Complexity,y=Recall,colour=Group)) +
    theme_classic() +
    theme(axis.title=element_text(size=16),
          legend.title=element_blank(),
          axis.text.x = element_text(size=14),
          axis.text.y = element_text(size=14)) +
    geom_line() + geom_point(size=4) +
    scale_x_discrete(limits=c(1:10)) +
    scale_colour_manual("legend", values = c(conventionalColour,ladderColour2,ladderColour3,
                                             ladderColour4,ladderColour5,ladderColour1))
  ggsave(device = "svg" ,g, height = 7 , width = 7 * aspect_ratio,filename = paste(outFilePath,"Stringtie_Re_",comparisonTool,"_acrossBetas.svg",sep = ""))





  #### Working with the full gtfs which are not divided by complexity
  ## Read the data
  ## Here i am using dummy values for recall and precision in order to plot them properly in the x axis
  ## recall = 85 and precision = 90
  # comparisonResultOriginal_full <- read.table(comparisonResultOrigFile_full, header = T, fill = T)
  # comparisonResultOriginal_full <- comparisonResultOriginal_full[c(1,2,3)]
  # names(comparisonResultOriginal_full) <- c("Recall","Precision","bandProb")
  # comparisonResultOriginal_full["Group"] <- "Conventional"
  # meltedOriginalResult <-  melt(comparisonResultOriginal_full,id.vars = c("Group","bandProb"))
  # meltedOriginalResult$value <- as.numeric(as.character(meltedOriginalResult$value))
  # meltedOriginalResult$variable <- ifelse(str_detect(meltedOriginalResult$variable,"Recall"), 85, 90)
  # 
  # comparisonResultMod_full <- read.table(comparisonResultModFile_full, header = T, fill = T)
  # comparisonResultMod_full <- comparisonResultMod_full[c(1,2,3)]
  # names(comparisonResultMod_full) <- c("Recall","Precision","bandProb")
  # comparisonResultMod_full["Group"] <- paste("LadderSeq",comparisonResultMod_full$bandProb,sep="_")
  # comparisonResultMod_full <- comparisonResultMod_full %>% subset(Recall!="Complexity")
  # meltedModResult <-  melt(comparisonResultMod_full,id.vars = c("Group","bandProb"))
  # meltedModResult$value <- as.numeric(as.character(meltedModResult$value))
  # meltedModResult <- meltedModResult  %>% mutate(Diff = -lag(meltedModResult$value,1)+ value)
  # meltedModResult$Diff <- ifelse(str_detect(meltedModResult$bandProb,"realistic"), meltedModResult$value, meltedModResult$Diff)
  # meltedModResult$variable <- ifelse(str_detect(meltedModResult$variable,"Recall"), 85, 90)
  # 
  # y <- factor(meltedModResult$Group)
  # y <-  factor(y,levels(y)[c(4,3,2,1,5)])
  # barwidth = 2
  # g <- ggplot() +
  #   theme_classic() +
  #   theme(axis.title=element_text(family="Sans",size=20),
  #         legend.title=element_blank(),
  #         legend.text.align = 0 ,
  #         legend.position="bottom" ,
  #         legend.text = element_text(family="Sans",size=20),
  #         axis.text.x = element_text(family="Sans",size=20),
  #         axis.title.x = element_blank(),
  #         axis.text.y = element_text(family="Sans",size=20)) +
  #   geom_bar(data = meltedOriginalResult,
  #            mapping = aes(x = variable - 1.1 , y = value, fill=Group),
  #            stat="identity",
  #            position='dodge',
  #            width = barwidth) +
  #   geom_bar(data = meltedModResult,
  #            mapping = aes(x = variable +1.1 , y = Diff,  fill=y),
  #            stat="identity",
  #            position='stack',
  #            width = barwidth) +
  #   scale_fill_manual("legend",
  #                     labels=c(expression("Conventional Trinity"),expression('Ladder_seq'^'1'),expression("Ladder_seq"^"2"),
  #                              expression("Ladder_seq"^"3"),expression("Ladder_seq"^"4"),"Ladder_seq"),
  #                     values = c(conventionalColour,ladderColour2,ladderColour3,
  #                                ladderColour4,ladderColour5,ladderColour1)) +
  #   xlab("") +
  #   ylab("% of Transcripts") +
  #   scale_x_continuous(breaks=c(85,90),labels=c("Recall","Precision"))+
  #   scale_y_continuous(expand = c(0, 0))
  # 
  # ggsave(device = "svg" ,g, height = 7 , width = 7 * aspect_ratio,filename = paste(outFilePath,"Combined_",comparisonTool,"_acrossBetas.svg",sep = ""))
}


plotFunction(comparisonResultOrigFile, comparisonResultModFile, outFilePath, dataset, comparisonTool, pipeline)
