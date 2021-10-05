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
require(scales)


#### Command line arguments
args = commandArgs(trailingOnly=TRUE)
resultOriginalFile <- args[1]
resultLadderMergeFile <- args[2]
combinedResultsFile <- args[3]
outFilePath <- args[4]
dataset <- args[5]
aspect_ratio <- 1.5
height <- 7
cols <- hue_pal()(2)



# ## Hardcoded arguments
# resultOriginalFile <-  "/algbio1/shounak/raw/simulated/human/Resimulation/30MillWithDenovo/simulations/sim_30000000/sim_1/CompRe_original_acrossBetas.txt"
# resultLadderMergeFile <- "/algbio1/shounak/raw/simulated/human/Resimulation/30MillWithDenovo/simulations/sim_30000000/sim_1/CompRe_ladder_acrossBetas.txt"
# combinedResultsFile <- "/algbio1/shounak/raw/simulated/human/Resimulation/30MillWithDenovo/simulations/sim_30000000/sim_1/denovoResults_acrossBetas.txt"
# outFilePath <- "/algbio1/shounak/raw/simulated/human/Resimulation/30MillWithDenovo/simulations/sim_30000000/sim_1/"
# dataset <- "30Million"



plotFunction <- function(resultOriginalFile, resultLadderMergeFile, combinedResultsFile, outFilePath, dataset){

  ## defining colours
  ladderColour1 <- rgb(222, 80, 0, 193, names = NULL, maxColorValue = 255)
  ladderColour2 <- rgb(230, 60, 0, 130, names = NULL, maxColorValue = 255)
  ladderColour3 <- rgb(235, 40, 0, 90, names = NULL, maxColorValue = 255)
  ladderColour4 <- rgb(240, 20, 0, 50, names = NULL, maxColorValue = 255)
  ladderColour5 <- rgb(255, 0, 0, 30, names = NULL, maxColorValue = 255)
  conventionalColour <- rgb(44, 118, 212, 255, names = NULL, maxColorValue = 255)

  ## Reading data for complexity-wise analysis
  resultOriginal <- read.table(resultOriginalFile, header = T)
  resultOriginal["Group"] <- "Conventional"
  resultOrigList <- split(resultOriginal, f = resultOriginal$PercentageLength)

  resultLadderMerge <- read.table(resultLadderMergeFile, header = T)
  resultLadderMerge["Group"] <- "LadderSeq"
  resultLadderList <- split(resultLadderMerge, f = resultLadderMerge$PercentageLength)

  for (i in c("95","90","85","80")){
    ## Plotting the est counts correlations
    resultOriginal_grouped <- as.data.frame(resultOrigList[i])
    colnames(resultOriginal_grouped) <- c("Complexity","Recall","Recall_num","PercentageLength","bandProb","Group")
    resultOriginal_grouped$Group <- paste(resultOriginal_grouped$Group,resultOriginal_grouped$bandProb,sep="_")

    resultLadder_grouped <- as.data.frame(resultLadderList[i])
    colnames(resultLadder_grouped) <- c("Complexity","Recall","Recall_num","PercentageLength","bandProb","Group")
    resultLadder_grouped$Group <- paste(resultLadder_grouped$Group,resultLadder_grouped$bandProb,sep="_")
    resultLadder_grouped <- resultLadder_grouped %>% subset(Group!="Group")

    percentageLength <- resultOriginal_grouped$PercentageLength[1]

    # Plotting complexity-wise recall
    combinedRecallFile <- rbind(resultOriginal_grouped, resultLadder_grouped)
    combinedRecallFile$Recall <- as.numeric(as.character(combinedRecallFile$Recall))
    combinedRecallFile$Complexity <- as.numeric(as.character(combinedRecallFile$Complexity))


    g <- ggplot(combinedRecallFile, aes(x=Complexity,y=Recall,colour=Group)) +
      theme_classic() +
      theme(axis.title=element_text(family="Sans",size=20),
            legend.title=element_blank(),
            legend.text.align = 0 ,
            legend.text = element_text(family="Sans",size=20),
            axis.text.x = element_text(family="Sans",size=20),
            axis.text.y = element_text(family="Sans",size=20)) +
            geom_line() + geom_point(size=4) +
            scale_x_discrete(limits=c(1:10)) +
      ylab("Sensitivity") +
      scale_colour_manual("legend",
                          labels=c(expression("Conventional Trinity"),expression('Ladder_seq'^'1'),expression("Ladder_seq"^"2"),
                                   expression("Ladder_seq"^"3"),expression("Ladder_seq"^"4"),"Ladder_seq"),
                          values = c(conventionalColour,ladderColour2,ladderColour3,
                                               ladderColour4,ladderColour5,ladderColour1))
    ggsave(device = "svg" ,g, height = 7 , width = 7 * aspect_ratio,filename = paste(outFilePath,"Complexity_Re_",percentageLength,"_acrossBetas.svg",sep = ""))
  }

  # Reading data for combined analysis
  combinedResults <- read.table(combinedResultsFile, header = T)
  names(combinedResults) <- c("Group", "Precision", "Recall","PercentageLength","bandProb","Recall_num")
  combinedResults <- combinedResults %>% subset(Group != "Group")
  meltedCombinedResults <- melt(combinedResults,id.vars = c("Group","PercentageLength","bandProb"))
  meltedCombinedResults$Group <- paste(meltedCombinedResults$Group,meltedCombinedResults$bandProb,sep="_")

  recallPlots <- meltedCombinedResults %>%  subset(variable=="Recall_num")
  recallPlots <- recallPlots %>% subset(!Group %in% c("Conventional_betas_increase_1","Conventional_betas_increase_4","Conventional_betas_increase_7","Conventional_betas_perfect"))
  recallPlots$value <- as.numeric(as.character(recallPlots$value))
  recallPlots$PercentageLength <- as.numeric(as.character(recallPlots$PercentageLength))

  recallPlotsConventional <- recallPlots %>% subset(Group=="Conventional_betas_realistic")
  recallPlotsLadder <- recallPlots %>% subset(Group!="Conventional_betas_realistic")
  recallPlotsLadder <- recallPlotsLadder %>%
                        mutate(Diff = lag(-recallPlotsLadder$value,4)+ value)

  recallPlotsLadder$Diff <- ifelse(is.na(recallPlotsLadder$Diff), recallPlotsLadder$value, recallPlotsLadder$Diff)
  y <- factor(recallPlotsLadder$Group)
  y <-  factor(y,levels(y)[c(4,3,2,1,5)])
  barwidth = 2
  g <- ggplot() +
    theme_classic() +
    theme(axis.title=element_text(family="Sans",size=20),
          legend.title=element_blank(),
          legend.text.align = 0 ,
          legend.position="bottom" ,
          legend.text = element_text(family="Sans",size=20),
          axis.text.x = element_text(family="Sans",size=20),
          axis.text.y = element_text(family="Sans",size=20,angle = 60)) +
    geom_bar(data = recallPlotsConventional,
             mapping = aes(x = PercentageLength - 1.1 , y = value, fill=Group),
             stat="identity",
             position='dodge',
             width = barwidth,
             orientation = 'x') +
    geom_bar(data = recallPlotsLadder,
             mapping = aes(x = PercentageLength + 1.1 , y = Diff,  fill=y),
             stat="identity",
             position='stack',
             width = barwidth,
             orientation = 'x') +
    scale_fill_manual("legend",
                      labels=c(expression("Conventional Trinity"),expression('Ladder_seq'^'1'),expression("Ladder_seq"^"2"),
                               expression("Ladder_seq"^"3"),expression("Ladder_seq"^"4"),"Ladder_seq"),
                      values = c(conventionalColour,ladderColour2,ladderColour3,
                                           ladderColour4,ladderColour5,ladderColour1)) +
    xlab("% Length reconstructed") + ylab("# of Transcripts") +
    scale_y_continuous(expand = c(0, 0))


  ggsave(device = "svg" ,g, height = 7 , width = 7 * aspect_ratio,filename = paste(outFilePath,"RecallComparison_acrossBetas.svg",sep = ""))






  precisionPlots <- meltedCombinedResults %>%  subset(variable=="Precision")
  precisionPlots <- precisionPlots %>% subset(!Group %in% c("Conventional_betas_increase_1","Conventional_betas_increase_4","Conventional_betas_increase_7","Conventional_betas_perfect"))
  precisionPlots$value <- as.numeric(as.character(precisionPlots$value))
  precisionPlots$PercentageLength <- as.numeric(as.character(precisionPlots$PercentageLength))

  precisionPlotsConventional <- precisionPlots %>% subset(Group=="Conventional_betas_realistic")
  precisionPlotsLadder <- precisionPlots %>% subset(Group!="Conventional_betas_realistic")
  precisionPlotsLadder <- precisionPlotsLadder %>%
    mutate(Diff = lag(-precisionPlotsLadder$value,4)+ value)

  precisionPlotsLadder$Diff <- ifelse(is.na(precisionPlotsLadder$Diff), precisionPlotsLadder$value, precisionPlotsLadder$Diff)
  y <- factor(precisionPlotsLadder$Group)
  y <-  factor(y,levels(y)[c(4,3,2,1,5)])
  barwidth = 2
  g <- ggplot() +
    theme_classic() +
    theme(axis.title=element_text(family="Sans",size=20),
          legend.title=element_blank(),
          legend.text.align = 0 ,
          legend.position="bottom" ,
          legend.text = element_text(family="Sans",size=20),
          axis.text.x = element_text(family="Sans",size=20),
          axis.text.y = element_text(family="Sans",size=20)) +
    geom_bar(data = precisionPlotsConventional,
             mapping = aes(x = PercentageLength - 1.1 , y = value, fill=Group),
             stat="identity",
             position='dodge',
             width = barwidth,
             orientation = 'x') +
    geom_bar(data = precisionPlotsLadder,
             mapping = aes(x = PercentageLength + 1.1 , y = Diff,  fill=y),
             stat="identity",
             position='stack',
             width = barwidth,
             orientation = 'x') +
    scale_fill_manual("legend",
                      labels=c(expression("Conventional Trinity"),expression('Ladder_seq'^'1'),expression("Ladder_seq"^"2"),
                               expression("Ladder_seq"^"3"),expression("Ladder_seq"^"4"),"Ladder_seq"),
                      values = c(conventionalColour,ladderColour2,ladderColour3,
                                           ladderColour4,ladderColour5,ladderColour1)) +
    xlab("% Length reconstructed") + ylab("Precision") +
    scale_y_continuous(expand = c(0, 0))

  ggsave(device = "svg" ,g, height = 7 , width = 7 * aspect_ratio,filename = paste(outFilePath,"PrecisionComparison_acrossBetas.svg",sep = ""))

}


plotFunction(resultOriginalFile, resultLadderMergeFile, combinedResultsFile, outFilePath, dataset)
