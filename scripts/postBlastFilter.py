import pysam
from argparse import ArgumentParser
import sys
from collections import defaultdict

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

def checkPosition(line, matchingIds):
    bwaChr = line.reference_name
    bwaPos = line.pos
    lineLen = line.infer_read_length()
    matching = []
    for xIter in matchingIds:
        blastChr = line.get_tag(f"c{xIter}")
        blastPos = line.get_tag(f"p{xIter}")
        if (
                blastChr == bwaChr 
                and (blastPos >= bwaPos - 2 * lineLen 
                     or blastPos <= bwaPos + 2 * lineLen
                     )
                ):
            # Maps correctly
            matching.append(xIter)
    if len(matching) == 0:
        return(False)
    else:
        return(True)
    

def main():
    parser = ArgumentParser()
    parser.add_argument('outPrefix')
    parser.add_argument('taxID', type=int)
    o = parser.parse_args()
    
    inBam = pysam.AlignmentFile("-", 'rb')
    # ~ print("Opened inbam")
    outGoodBam = pysam.AlignmentFile("-", 'wb', template=inBam)
    outBadBam = pysam.AlignmentFile(f"{o.outPrefix}.wrongSpecies.bam", 'wb', template=inBam)
    outAmbigBam = pysam.AlignmentFile(f"{o.outPrefix}.ambig.bam", 'wb', template=inBam)
    # ~ print("Opened outBams")
    outTextID = open(f"{o.outPrefix}.speciesComp.txt" ,'w')
    firstLine = next(inBam)
    lineStorage = [firstLine]
    FinalValue = pysam.AlignedSegment()
    FinalValue.query_name = "FinalValue"
    taxID_dict = defaultdict(int)
    for line in iteratorWrapper(inBam.fetch(until_eof=True), FinalValue):
        if line.query_name == firstLine.query_name:
            lineStorage.append(line)
        else:
            # check if either line has a t0 tag
            addedTag = False
            # ~ tagVals1 = [x.get_tag("YT") if x.has_tag("YT") else o.taxID for x in lineStorage]
            # ~ family = [x.query_name for x in lineStorage]
            for xIter in lineStorage:
                if not xIter.has_tag("t0"):
                    xIter.set_tag("t0",o.taxID, 'i')
                    xIter.set_tag("am", 5, 'i')
                    xIter.set_tag("YB", "False", 'Z')
                    addedTag = True
                    
            tagVals = [x.get_tag("t0") for x in lineStorage]
            tagSet=set(tagVals)
            testSet = {o.taxID, -1}
            tagTest = tagSet - testSet
            # ~ sys.stderr.write(f"{family}\n{tagVals}\n{tagSet},{addedTag}\n")
            if len(tagTest) == 0:
                if o.taxID in tagSet:
                    # good pair
                    altTaxID = o.taxID
                else:
                    # blast indeterminate
                    altTaxID = -1
            elif len(tagTest) == 1:
                if o.taxID in tagSet and addedTag == False:
                    # indeterminate
                    altTaxID = -1
                else:
                    # single bad determination
                    altTaxID = next(iter(tagTest))
                    # ~ sys.stderr.write(f"{tagSet}\n")
            else:
                # multiple bad determinations
                altTaxID = -1
            
            # increment counter
            
            for lIter in lineStorage:
                lIter.set_tag("t0",altTaxID, 'i')
            
            # check altTaxID
            if altTaxID == o.taxID:
                # correct species read pair
                # process pseudogene checking
                amVals = [x.get_tag("am") for x in lineStorage]
                amSet=set(amVals)
                amTSet = {0, 5}
                amTest = amSet - amTSet
                if len(amTest) == 1:
                    if 3 in amTest:
                        testInd = [x for x in range(len(amVals)) if amVals[x] == 3]
                        anyMatch = False
                        matchingInds = {}
                        for lIter in testInd:
                            tagIter = 1
                            matchingInds[lIter] = []
                            while lineStorage[lIter].has_tag(f"t{tagIter}"):
                                if lineStorage[lIter].get_tag(f"t{tagIter}") == o.taxID:
                                    matchingInds[lIter].append(tagIter)
                                tagIter += 1
                            if matchingInds[lIter] != []:
                                anyMatch = True
                        if not anyMatch:
                            for lIter in lineStorage:
                                lIter.set_tag("t0",-1, 'i')
                                outBadBam.write(lIter)
                            altTaxID = -1
                            
                        elif True in [len(matchingInds[x]) > 1 for x in matchingInds]:
                            # multiple correct species hit; mark ambiguous
                            for lIter in lineStorage:
                                outAmbigBam.write(lIter)
                        else:
                            ############################################
                            # Revisit in the future to test for if, of 
                            # two reads, one matches and the other doesn't, 
                            # we need to do something else.  
                            ############################################
                            matchingPos = True
                            for lIter in matchingInds:
                                if len(matchingInds[lIter]) == 1:
                                    # Check position
                                    if not checkPosition(lineStorage[lIter], matchingInds[lIter]):
                                        matchingPos = False
                            if matchingPos:
                                for lIter in lineStorage:
                                    outGoodBam.write(lIter)
                            else:
                                for lIter in lineStorage:
                                    if lIter.get_tag('am') != 5:
                                        lIter.set_tag("am", 1,'i')
                                    outAmbigBam.write(lIter)
                    elif 2 in amTest:
                        for lIter in lineStorage:
                            outAmbigBam.write(lIter)
                    else:
                        for lIter in lineStorage:
                            lIter.set_tag("t0",-1, 'i')
                            outBadBam.write(lIter)
                        outTaxID = -2
                elif len(amTest) == 0:
                    testInd = [x for x in range(len(amVals)) if amVals[x] == 0]
                    if len(testInd) == 0:
                        for lIter in lineStorage:
                            outGoodBam.write(lIter)
                    else:
                        matchingInds = {}
                        for lIter in testInd:
                            tagIter = 1
                            matchingInds[lIter] = []
                            while lineStorage[lIter].has_tag(f"c{tagIter}"):
                                matchingInds[lIter].append(tagIter)
                                tagIter += 1
                        if True in [len(matchingInds[x]) > 1 for x in matchingInds]:
                            # multiple correct species hit; mark ambiguous
                            for lIter in lineStorage:
                                outAmbigBam.write(lIter)
                        else:
                            ############################################
                            # Revisit in the future to test for if, of 
                            # two reads, one matches and the other doesn't, 
                            # we need to do something else.  
                            ############################################
                            matchingPos = True
                            for lIter in matchingInds:
                                if len(matchingInds[lIter]) == 1:
                                    # Check position
                                    if not checkPosition(lineStorage[lIter], matchingInds[lIter]):
                                        matchingPos = False
                            if matchingPos:
                                for lIter in lineStorage:
                                    outGoodBam.write(lIter)
                            else:
                                for lIter in lineStorage:
                                    if lIter.get_tag('am') != 5:
                                        lIter.set_tag("am", 1,'i')
                                    outAmbigBam.write(lIter)
                else:
                    for lIter in lineStorage:
                        outAmbigBam.write(lIter)
            else:
                for lIter in lineStorage:
                    outBadBam.write(lIter)
            
            taxID_dict[altTaxID] += len(lineStorage)
            
            firstLine = line
            lineStorage = [firstLine]
    
    for tID in taxID_dict:
        outTextID.write(f"{tID}\t{taxID_dict[tID]}\n")
    inBam.close()
    outGoodBam.close()
    outBadBam.close()
    outTextID.close()

if __name__ == "__main__":
    main()