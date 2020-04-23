import pysam
from argparse import ArgumentParser
import logging
import os

class iteratorWrapper:
    def __init__(self, inIterator, finalValue):
        self.it = inIterator
        self.finalValue = finalValue
        self.endIter = False
    def __iter__(self):
        return(self)
    def __next__(self):
        try:
            temp = next(self.it)
        except StopIteration:
            if self.endIter == False:
                temp = self.finalValue
                self.endIter = True
            else:
                raise(StopIteration)
        return(temp)
    next = __next__

def main():
    parser = ArgumentParser()
    parser.add_argument("inBam")
    parser.add_argument("outBam")
    parser.add_argument("ambigBam")
    parser.add_argument("goodChr")
    o = parser.parse_args()
    #class o:
    #    pass
    #o.inBam = "/home/data/mitoTest/Mouse202/MouseTest_dcs.ambig.sort.bam"
    #o.outBam = "/home/data/mitoTest/Mouse202/MouseTest_dcs.ambig.recover.bam"
    #o.goodChr = "chrM"
    numeric_level = getattr(logging, "INFO", None)
    if not isinstance(numeric_level, int):
        raise ValueError('Invalid log level: {o.logLvl}')
    logging.basicConfig(
        format='%(levelname)s: %(message)s', 
        level=numeric_level, 
        )
    
    inBam = pysam.AlignmentFile(o.inBam, 'rb')
    resort = False
    
    if inBam.header.as_dict()['HD']['SO'] != "queryname":
        inBam.close()
        logging.info("Input not querry name sorted.  Sorting into query-name order...")
        pysam.sort("-n", "-o", f"{o.inBam}.temp.sort.bam", o.inBam)
        resort = True
        inBam = pysam.AlignmentFile(f"{o.inBam}.temp.sort.bam", 'rb')
    if resort:
        outBam = pysam.AlignmentFile(f"{o.outBam}.temp.bam", 'wb', template=inBam)
        outAmbig = pysam.AlignmentFile(f"{o.ambigBam}.temp.bam", 'wb', template=inBam)
    else:
        outBam = pysam.AlignmentFile(o.outBam, 'wb', template=inBam)
        outAmbig = pysam.AlignmentFile(o.ambigBam, 'wb', template=inBam)
    try:
        firstLine = next(inBam)
        contRecover = True
    except StopIteration:
        logging.info("No reads to recover!")
        contRecover = False
    if contRecover:
        lineStorage = [firstLine]
        FinalValue = pysam.AlignedSegment()
        FinalValue.query_name = "FinalValue#ab:1"
        
        for line in iteratorWrapper(inBam.fetch(until_eof=True), FinalValue):
            if line.query_name == firstLine.query_name:
                lineStorage.append(line)
            else:
                amVals = [x.get_tag("am") for x in lineStorage]
                amSet=set(amVals)
                amTSet = {0, 2, 5}
                amTest = amSet - amTSet
                if len(amTest) == 0:
                    goodPos = False
                    for read in lineStorage:
                        if read.has_tag('am'):
                            if read.get_tag('am') == 2:
                                # test for chrom match
                                chrTest = 1
                                while read.has_tag(f"c{chrTest}"):
                                    if read.get_tag(f"c{chrTest}") == o.goodChr:
                                        goodPos = True
                                    chrTest += 1
                    if goodPos:
                        for read in lineStorage:
                            outBam.write(read)
                    else:
                        for read in lineStorage:
                            outAmbig.write(read)
                firstLine = line
                lineStorage = [firstLine]
    outBam.close()
    if resort:
        logging.info("Resorting into coordinate order...")
        pysam.sort("-o", o.outBam, f"{o.outBam}.temp.bam")
        pysam.sort("-o", o.outAmbig, f"{o.outAmbig}.temp.bam")
        logging.info("Cleaning up temporary files")
        os.remove(f"{o.outBam}.temp.bam")
        os.remove(f"{o.inBam}.temp.sort.bam")
        os.remove(f"{o.outAmbig}.temp.bam")

    logging.info("Done")

if __name__ == "__main__":
    main()