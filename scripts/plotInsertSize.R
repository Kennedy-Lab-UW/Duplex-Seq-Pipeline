# Load libraries
library(tidyr)
library(readr)
library(ggplot2)
library(grid)
library(gridExtra)
args = commandArgs(trailingOnly=TRUE)


myPlots = list()
# load data

inSampName = args[1]

myInsertSizes <- read_delim(paste("Stats/data/", inSampName, ".iSize_Metrics.txt", sep=""), 
                            "\t", escape_double = FALSE, comment = "#", 
                            trim_ws = TRUE, skip=10)

numReads = sum(myInsertSizes$All_Reads.fr_count)
totReads = numReads
readCtr = c()
insertSizeCtr = 0
for ( xIter in seq(1,length(myInsertSizes$All_Reads.fr_count)) ) {
  readCtr = c(readCtr, numReads)
  numReads = numReads - myInsertSizes$All_Reads.fr_count[xIter]
  insertSizeCtr = insertSizeCtr + (myInsertSizes$insert_size[xIter] * myInsertSizes$All_Reads.fr_count[xIter])
}
readCtr = readCtr / totReads
myInsertSizes$CumProp = readCtr

myPlot = ggplot(data=myInsertSizes, mapping = aes(x=insert_size)) + 
  geom_bar(mapping=aes(y=All_Reads.fr_count), stat = "identity") + 
  geom_line(mapping = aes(y=CumProp * max(myInsertSizes$All_Reads.fr_count)), color="red") + 
  labs(title = inSampName, x = "Insert Size", y="Count") + 
  scale_y_continuous(sec.axis = sec_axis(~./max(myInsertSizes$All_Reads.fr_count), name="Cumulative Proportion")) + 
  theme_bw()
ggsave(
  filename=paste("Stats/plots/", inSampName, ".iSize_Histogram.png", sep=""), 
  plot=myPlot, height=6, width=6, units="in", 
  device = "tiff", compression="lzw"
  )
