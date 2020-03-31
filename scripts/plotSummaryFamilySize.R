library(tidyr)
library(readr)
library(ggplot2)
library(dplyr)
library(grid)
library(gridExtra)
args = commandArgs(trailingOnly=TRUE)
# Read Config
myConfig <- read_csv(args[1])

# Family size
myPlots = list()
for ( sampIter in seq(length(row.names(myConfig))) ) {
  
  # load data
  
  inSampName = myConfig$sample[sampIter]
  
  myTagstats <- read_delim(
    paste(
      myConfig$baseDir[sampIter],
      "/Stats/data/",
      inSampName, 
      ".tagstats.txt", 
      sep=""
    ), 
    "\t", 
    col_names = c(
      "familySize", 
      "numFamilies", 
      "propReads"
    ), 
  )
  
  myPlot = ggplot(data=myTagstats, mapping = aes(x=familySize, y=propReads)) + 
    geom_bar(stat = "identity", position = "stack") + 
    scale_y_continuous(name = "Proportion of Total Reads") + xlab("Family Size") + 
    ggtitle(inSampName) + 
    theme_bw()
  myPlots[[sampIter]] = ggplotGrob(myPlot)
}

finPlot = marrangeGrob(grobs = myPlots, ncol = 2, nrow = 2, top=NULL)

ggsave(paste(args[1],"summaryFamilySize.pdf", sep="."), plot=finPlot, width = 200, 
       height=200, units = "mm", limitsize = FALSE)