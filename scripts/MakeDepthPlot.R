args = commandArgs(trailingOnly=TRUE)

inBed = args[1]
inSampName = args[2]
inSampType = args[3]

library(ggplot2)
library(readr)
# read in bed file
# input:
#        inBed = get_target_bed,
#        inDepth = "{runPath}/Stats/data/{sample}.{sampType}.region.mutpos.vcf_depth.txt",

#    output:
#        "{runPath}/Stats/plots/{sample}.{sampType}.targetCoverage.png"

myBed <- read_table2(inBed, col_names = FALSE)
myFName = paste("Stats/data/",inSampName, ".region.mutpos.vcf_depth.txt", sep="")
depth <- read_delim(myFName,
                    "\t", escape_double = FALSE, trim_ws = TRUE, 
                    col_types="cicii") 
# Set column names
numCols = length(colnames(myBed))
if (numCols == 3) {
  colnames(myBed) <- c("Chrom","Start","End")
  myBed$Name = row.names(myBed)
} else if (numCols == 4) {
  colnames(myBed) <- c("Chrom","Start","End", "Name")
} else if (numCols == 5) {
  colnames(myBed) <- c("Chrom","Start","End", "Name","Score")
} else if (numCols == 6) {
  colnames(myBed) <- c("Chrom","Start","End", "Name","Score","Strand")
} else if (numCols == 7) {
  colnames(myBed) <- c("Chrom","Start","End", "Name","Score","Strand","thickStart")
} else if (numCols == 8) {
  colnames(myBed) <- c("Chrom","Start","End", "Name","Score","Strand","thickStart","thickEnd")
} else if (numCols == 9) {
  colnames(myBed) <- c("Chrom","Start","End", "Name","Score","Strand","thickStart","thickEnd","itemRgb")
} else if (numCols == 10) {
  colnames(myBed) <- c("Chrom","Start","End", "Name","Score","Strand","thickStart","thickEnd","itemRgb","blockCounts")
}  else if (numCols == 11) {
  colnames(myBed) <- c("Chrom","Start","End", "Name","Score","Strand","thickStart","thickEnd","itemRgb","blockCounts","blockSizes")
} else if (numCols == 12) {
  colnames(myBed) <- c("Chrom","Start","End", "Name","Score","Strand","thickStart","thickEnd","itemRgb","blockCounts","blockSizes","blockStarts")
} 
myBed$Start <- myBed$Start + 1
myBed$End <- myBed$End + 1
myBed$Target = factor(myBed$Name, levels = c(myBed$Name,"Off_Target"))
namesVect = c()
if (length(row.names(depth)) > 0) {
    for (dataIter in seq(1,length(row.names(depth)))) {
      myName = "Off_Target"
      for (bedIter in seq(1, length(row.names(myBed)))) {
        if ( depth$`#Chrom`[dataIter] == myBed$Chrom[bedIter] && 
             depth$Pos[dataIter] >= myBed$Start[bedIter] && 
             depth$Pos[dataIter] < myBed$End[bedIter]) {
          myName = myBed$Name[bedIter]
          break
        }
      }
      namesVect = c(namesVect, myName)
    }
    maxDP = max(depth$DP)
    depth$Nfract = depth$Ns / depth$DP * 100
    depth$Nfract[depth$DP == 0] = 0
    depth$Target = factor(namesVect, levels = c(myBed$Name,"Off_Target"))
    maxNs = ceiling(max(depth$Ns[depth$Target != "Off_Target"])/3)*3
    
} else {
    maxDP = 1
    maxNs = 1
    depth$Target = factor(namesVect, levels = c(myBed$Name,"Off_Target"))
}
if (maxNs == 0) {
  maxNs = 1
}
if (maxDP == 0) {
  maxDP = 1
}
multiplier = -maxDP / maxNs


myFName = paste("Final/",inSampType,"/",inSampName, ".vcf", sep="")
myVCF <- read_delim(
  myFName, "\t", 
  escape_double = FALSE, 
  comment = "##", 
  trim_ws = TRUE
  )
namesVect = c()
if (length(row.names(myVCF)) > 0) {
for (dataIter in seq(1,length(row.names(myVCF)))) {
  myName = "Off_Target"
  for (bedIter in seq(1, length(row.names(myBed)))) {
    if ( myVCF$`#CHROM`[dataIter] == myBed$Chrom[bedIter] && 
         myVCF$POS[dataIter] >= myBed$Start[bedIter] && 
         myVCF$POS[dataIter] < myBed$End[bedIter]) {
      myName = myBed$Name[bedIter]
      break
    }
  }
  namesVect = c(namesVect, myName)
}
myVCF$Target = factor(namesVect, levels = c(myBed$Name,"Off_Target"))
} else {myVCF$Target = factor()}
myPlot=ggplot(data = depth[depth$Target != "Off_Target",]) + 
  geom_area(mapping = aes(x = Pos, y = DP), stat="identity") + 
  geom_area(mapping=aes(x=Pos, y=Ns * multiplier), color='red', alpha=0.7) + 
  geom_segment(data = myBed, mapping = aes(x=Start,xend=End, y=0, yend=0))
if (length(row.names(myVCF)) > 0) {
  myPlot = myPlot + 
  geom_point(data=myVCF[myVCF$Target != "Off_Target",], mapping=aes(x=POS, y=maxDP), shape=1)
}
myPlot = myPlot + 
  scale_y_continuous(
    breaks = c(0, maxDP), 
    sec.axis = sec_axis(~./multiplier, name="Ns (count)", breaks = c(0,maxNs/3, 2*maxNs/3,maxNs), labels = c(0,maxNs/3, 2*maxNs/3,maxNs))) + 
  labs(title = inSampName) + 
  theme_bw() + 
  theme(
    axis.title.y.right = element_text(color='red'), 
    axis.text.x = element_text(angle=90)
  ) + 
  facet_wrap(. ~ Target, scales = "free_x", ncol = min(length(myBed$Name), 4))
ggsave(filename = paste("Stats/plots/", inSampName, ".targetCoverage.png", sep=""), 
       plot=myPlot, 
       width = 200, 
       height=50*ceiling(length(levels(depth$Target)) / 4), 
       units="mm", 
       device = "png")
