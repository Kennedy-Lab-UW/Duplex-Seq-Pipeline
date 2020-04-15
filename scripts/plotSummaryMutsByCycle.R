library(tidyr)
library(readr)
library(ggplot2)
library(dplyr)
library(grid)
library(gridExtra)
args = commandArgs(trailingOnly=TRUE)
# Read Config
myConfig <- read_csv(args[1])

# mutsByReadPos
myPlots = list()
colorTrans = c("C>A"='#5abdeb', 
               "C>G"='#050708', 
               "C>T"='#d43c32', 
               "T>A"='#cbcacb', 
               "T>C"='#aacb72', 
               "T>G"='#e7c9c6', 
               "C>N"="#4B0082", 
               "T>N"="#4B0082")
for ( sampIter in seq(length(row.names(myConfig))) ) {
  # load data
  
  inSampName = myConfig$sample[sampIter]
  
  myMutsByCyc <- read_delim(
    paste(
      myConfig$baseDir[sampIter],
      "/Stats/data/",
      inSampName, 
      ".dcs_MutsPerCycle.dat.csv", 
      sep=""
    ), 
    ",", 
    col_names = c(
      "Cycle",
      "C>T", "C>A", "C>G",
      "T>A", "T>C", "T>G",
      "C>N", "T>N", 
      "Count","Base","Count_percent"
    ), 
    skip = 1
  )
  myMutsByCyc$totMuts = myMutsByCyc$`C>T` + 
    myMutsByCyc$`C>G` +
    myMutsByCyc$`C>A` + 
    myMutsByCyc$`T>A` + 
    myMutsByCyc$`T>C` + 
    myMutsByCyc$`T>G`
  myMutsByCyc <- pivot_longer(
    myMutsByCyc, 
    c(`C>T`,`C>A`,`C>G`,`T>A`,`T>C`,`T>G`,`C>N`,`T>N`), 
    names_to = "mutType", values_to = "Number"
  )
  myMutsByCyc <- separate(myMutsByCyc, mutType, c("mutFrom", "mutTo"), sep = ">", remove = FALSE)
  
  maxReads = max(myMutsByCyc$Count)
  readCtr = c()
  insertSizeCtr = 0
  maxCt = max(myMutsByCyc[myMutsByCyc$mutTo != "N",]$totMuts, 1)
  
  myPlot = ggplot(data=myMutsByCyc[myMutsByCyc$mutTo != "N",], mapping = aes(x=Base)) + 
    geom_bar(mapping=aes(y=Number, fill=mutType), stat = "identity", position = "stack") + 
    geom_line(mapping = aes(y=Count_percent * maxCt), linetype = "dashed") + 
    scale_y_continuous(
      sec.axis = sec_axis(~./maxCt, 
                          name="Fraction of Total Reads"), name = "Count")  + 
    scale_fill_manual(values = colorTrans) + 
    ggtitle(inSampName) + 
    theme_bw()
  myPlots[[sampIter]] = ggplotGrob(myPlot)
}

finPlot = marrangeGrob(grobs = myPlots, ncol = 1, nrow = 4, top=NULL)

ggsave(paste(args[1],"summaryMutsByCycle.pdf", sep="."), plot=finPlot, width = 200, 
       height=200, units = "mm", limitsize = FALSE)