import logging
import sys

def Bed_File(fName, mode='r'):
    if mode == 'r':
        return Bed_Reader(fName)
    elif mode == 'w':
        return Bed_Writer(fName)
    else:
        raise ValueError(f"Invalid mode: {mode}")

class Bed_Reader:
    def __init__(self, inFile):
        self.file = open(inFile, 'r')
        
    def __iter__(self):
        return(self)
    
    def __next__(self):
        myLine = self.file.readline()
        if myLine != "":
            return(Bed_Line(*myLine.strip().split()))
        else:
            raise StopIteration

class Bed_Writer:
    def __init__(self, fName):
        self.file = open(fName, 'w')
    def writeline(self, line):
        if type(line) == Bed_Line:
            self.file.write(f"{str(line)}")
        else:
            raise TypeError("Anything written to a bed file must be a Bed_Line")

class Bed_Line:
    def __init__(self, 
                 chrom, 
                 startPos, 
                 endPos, 
                 name="",
                 score="", 
                 strand="+",
                 thickStart="", 
                 thickEnd="",
                 itemRGB="", 
                 blockCount="",
                 blockSizes="", 
                 blockStarts=""):
        self.chrom = chrom
        self.startPos = int(startPos)
        self.endPos = int(endPos)
        if name == "":
            self.name = f"{self.chrom}:{self.startPos}-{self.endPos}"
        else:
            self.name = name
        if score == '':
            self.score = 0
        else:
            self.score = int(score)
        if strand in (".","+","-"):
            self.strand = strand
        else:
            self.strand = "."
        if thickStart == "":
            self.thickStart = self.startPos
        else:
            self.thickStart = int(thickStart)
        if thickEnd == "":
            self.thickEnd = self.endPos
        else:
            self.thickEnd = int(thickEnd)
        if itemRGB == "":
            self.itemRGB = None
        else:
            self.itemRGB = itemRGB
        if blockCount == "":
            self.blockCount = 1
        else:
            self.blockCount = int(blockCount)
        if blockSizes == "":
            self.blockSizes = [abs(self.endPos - self.startPos)]
        else:
            self.blockSizes = list(filter(None,[int(x) for x in blockSizes.split(',')]))
        if blockStarts == "":
            self.blockStarts = [0]
        else:
            self.blockStarts = list(filter([int(x) for x in blockStarts.split(',')]))
        assert self.blockCount == len(self.blockSizes)
        assert self.blockCount == len(self.blockStarts)

    def samtoolsStr(self):
        return f"{self.chrom}:{self.startPos}:{self.endPos}"
    
    def get_subregions(self):
        if self.strand == '+':
            return [Bed_Line(
                self.chrom, 
                self.startPos + self.blockStarts[x], 
                self.startPos + self.blockStarts[x] + self.blockSizes[x], 
                f"{self.name}_block{x + 1}", 
                self.score, 
                self.strand
                ) for x in range(len(self.blockStarts))]
        elif self.strand == '-':
            return [Bed_Line(
                self.chrom, 
                self.startPos + self.blockStarts[-x - 1], 
                self.startPos + self.blockStarts[-x - 1] + self.blockSizes[-x - 1], 
                f"{self.name}_block{x + 1}", 
                self.score, 
                self.strand
                ) for x in range(len(self.blockStarts))]
        
    def contains(self, inChr, inPos):
        if (
                inChr == self.chrom
                and inPos >= self.startPos
                and inPos < self.endPos
                ):
            return True
        else:
            return False
    
    def __len__(self):
        return self.endPos - self.startPos - 1

    def __repr__(self):
        return f"Bed_Line:{str(self)}"

    def __str__(self):
        bSizes = ",".join(self.blockSizes)
        bStarts = ",".join(self.bStarts)
        outString = (
            f"{self.chrom}\t{self.startPos}\t{self.endPos}\t"
            f"{self.name}\t{self.score}\t{self.strans}\t{self.thickStart}\t"
            f"{self.thickEnd}\t{self.itemRGB}\t{self.blockCount}\t{bSizes}\t"
            f"{bStarts}\n")
        return outString