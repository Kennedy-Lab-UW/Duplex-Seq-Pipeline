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
            self.blockSizes = [int(x) for x in blockSizes.split(',') if x != ""]
        if blockStarts == "":
            self.blockStarts = [0]
        else:
            self.blockStarts = [int(x) for x in blockStarts.split(',') if x != ""]
        assert self.blockCount == len(self.blockSizes)
        assert self.blockCount == len(self.blockStarts)
        self.subregions = None

    def samtoolsStr(self):
        return f"{self.chrom}:{self.startPos}:{self.endPos}"
    
    def get_subregions(self):
        if self.subregions is not None:
            return self.subregions
        if self.strand == '+':
            self.subregions = [Bed_Line(
                self.chrom, 
                self.startPos + self.blockStarts[x], 
                self.startPos + self.blockStarts[x] + self.blockSizes[x], 
                f"{self.name}_block{x + 1}", 
                self.score, 
                self.strand
                ) for x in range(len(self.blockStarts))]
        elif self.strand == '-':
            self.subregions =  [Bed_Line(
                self.chrom, 
                self.startPos + self.blockStarts[-x - 1], 
                self.startPos + self.blockStarts[-x - 1] + self.blockSizes[-x - 1], 
                f"{self.name}_block{x + 1}", 
                self.score, 
                self.strand
                ) for x in range(len(self.blockStarts))]
        return self.subregions
        
    def contains(self, inChr, inPos, blocks_only=False):
        if blocks_only:
            if self.subregions is None:
                self.get_subregions()
            for subregion in self.subregions:
                if subregion.contains(inChr, inPos):
                    return True
            return False
        if (
                inChr == self.chrom
                and inPos >= self.startPos
                and inPos < self.endPos
                ):
            return True
        else:
            return False
    
    def __len__(self):
        return self.endPos - self.startPos

    def __repr__(self):
        return f"Bed_Line:{str(self)}"

    def __str__(self):
        bSizes = ",".join([str(x) for x in self.blockSizes])
        bStarts = ",".join([str(x) for x in self.blockStarts])
        outString = (
            f"{self.chrom}\t{self.startPos}\t{self.endPos}\t"
            f"{self.name}\t{self.score}\t{self.strand}\t{self.thickStart}\t"
            f"{self.thickEnd}\t{self.itemRGB}\t{self.blockCount}\t{bSizes}\t"
            f"{bStarts}\n")
        return outString
