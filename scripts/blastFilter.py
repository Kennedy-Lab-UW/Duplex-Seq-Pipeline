from argparse import ArgumentParser
from collections import namedtuple, Counter
import pysam
from Bio.Blast import NCBIXML


class multiparser:
    def __init__(self, inBam, inXml):
        self.bamFile=pysam.AlignmentFile(inBam, 'rb')
        self.xml = NCBIXML.parse(open(inXml))
        #self.log = open("testLog.txt", 'w')
        self.bamLineNum = 0
        self.xmlLineNum = 0
        self.bamLine = None
        self.xmlLine = None
        
    def __iter__(self):
        return(self)
    def __next__(self):
        try:
            if self.bamLine is None and self.xmlLine is None:
                self.bamLine = next(self.bamFile)
                self.xmlLine = next(self.xml)
                self.bamLineNum += 1
                self.xmlLineNum += 1
            elif self.bamLine.query_name != self.xmlLine.query.split('/')[0]:
                self.bamLine = next(self.bamFile)
                self.bamLineNum += 1
            else:
                self.bamLineNum += 1
                self.xmlLineNum += 1
                self.bamLine = next(self.bamFile)
                self.xmlLine = next(self.xml)
            
#            while  bamLine.query_name != xmlLine.query.split('/')[0]:
#                bamLine = next(self.bamFile)
#                self.log.write(f"Bam Name:\t{bamLine.query_name}\n")
            
            if self.bamLine.query_name == self.xmlLine.query.split('/')[0]:
                return(self.bamLine, self.xmlLine)
            else:
                return(self.bamLine, None)
        except StopIteration:
            raise StopIteration
        except Exception:
            print(self.bamLine.query_name)
            print(self.xmlLine.query)
            raise

blastAlignment = namedtuple(
    "blastAlignment", 
    [
        "alnNum",
        "ID",
        "start", 
        "end",
        "score"
        ]
    )

def extractBlastData(inRecord):
    alignments = []
    myAlnNum = 0
    # ~ print(inRecord.query)
    for aln in inRecord.alignments:
        for hsp in aln.hsps:
            # ~ print(aln.title)
            # ~ print(hsp.expect)
            alignments.append(
                blastAlignment(
                    myAlnNum, 
                    aln.title, 
                    hsp.sbjct_start, 
                    hsp.sbjct_end, 
                    hsp.expect
                    )
                )
            myAlnNum += 1
    # ~ print()
    alignments.sort(key=lambda aln: aln.score)
    # ~ sys.stderr.write(f"alignments = {alignments}")
    try:
        det = {"queryName":inRecord.query}
        if len(alignments) == 0: 
            # Set other alignment fields <---------------------------------------- ***WORK HERE***
            det["t0"] = -3
            det["am"] = 4
            det["MaxAln"] = 0
        else:
            # Modify to deal with having more than two maximum alignments <------- ***WORK HERE***
            # ~ sys.stderr.write("Entering while...\n")
            
            tieCounter = Counter([x.score for x in alignments])
            tiedAlignments = [alignments[x] for x in range(len(alignments)) if alignments[x].score == min(tieCounter)]
            if len(tiedAlignments) > 1:
                # Figure out which alignments are tied for max value
                # check for type of tie
                taxIDs = {x.ID.split()[1].split('|')[0] for x in tiedAlignments}
                if len(taxIDs) == 1:
                    # within species tie; possible pseudogene
                    det["t0"] = int(tiedAlignments[0].ID.split()[1].split('|')[0])
                    for aln in range(len(tiedAlignments)):
                        det[f"c{aln + 1}"] = tiedAlignments[aln].ID.split()[1].split('|')[1]
                        det[f"p{aln + 1}"] = min(tiedAlignments[aln].start, tiedAlignments[aln].end)
                        det[f"l{aln + 1}"] = abs(tiedAlignments[aln].start - tiedAlignments[aln].end)
                    det["am"] = 2
                    det["MaxAln"] = len(tiedAlignments)
                else:
                    # between species tie; ambiguous.  
                    det["t0"] = -1
                    for aln in range(len(tiedAlignments)):
                        det[f"t{aln + 1}"] = int(tiedAlignments[aln].ID.split()[1].split('|')[0])
                        det[f"c{aln + 1}"] = tiedAlignments[aln].ID.split()[1].split('|')[1]
                        det[f"p{aln + 1}"] = min(tiedAlignments[aln].start, tiedAlignments[aln].end)
                        det[f"l{aln + 1}"] = abs(tiedAlignments[aln].start - tiedAlignments[aln].end)
                    det["am"] = 3
                    det["MaxAln"] = len(tiedAlignments)
            else:
                # Only one best alignment
                det["t0"] = int(tiedAlignments[0].ID.split()[1].split('|')[0])
                det["c1"] = tiedAlignments[0].ID.split()[1].split('|')[1]
                det["p1"] = min(tiedAlignments[0].start, tiedAlignments[0].end)
                det["l1"] = abs(tiedAlignments[0].start - tiedAlignments[0].end)
                det["am"] = 0
                det["MaxAln"] = 1
    except Exception:
        print(alignments)
        print(len(alignments))
        print(tiedAlignments)
        print(len(tiedAlignments))
        raise
    return(det)

def main():
    parser = ArgumentParser()
    parser.add_argument('inBam', action='store')
    parser.add_argument('inXML', action='store')
    parser.add_argument('outPrefix', action='store')
    o = parser.parse_args()

    myIterator = multiparser(o.inBam,o.inXML)
    outBam = pysam.AlignmentFile(f"{o.outPrefix}.speciesLabeled.bam", 'wb', template=myIterator.bamFile)

    records = 0
    for bamLine, blastLine in myIterator:
        if bamLine.query_name == "TGGATACGGTTGATGGCGAAGTGGTTCTCATG":
            print(bamLine.tags)
        # ~ if records %  == 0:
        # ~ sys.stderr.write(f"{records} lines processed...\n")
        if blastLine is not None:
            # ~ sys.stderr.write(f"Determining blast result...\n")
            blastDet = extractBlastData(blastLine)
            # ~ sys.stderr.write(f"Determined blast result...\n")
            # ~ sys.stderr.write(f"{blastDet}")
            bamLine.set_tag("t0",blastDet["t0"], 'i')
            if bamLine.query_name == "TGGATACGGTTGATGGCGAAGTGGTTCTCATG":
                # ~ print(bamLine.tags)
                print(blastDet)
            try:
                
                for aln in range(blastDet["MaxAln"]):
                    if f"t{aln + 1}" in blastDet:
                         bamLine.set_tag(f"t{aln + 1}", blastDet[f"t{aln + 1}"], 'i')
                    bamLine.set_tag(f"c{aln + 1}", blastDet[f"c{aln + 1}"], 'Z')
                    bamLine.set_tag(f"p{aln + 1}", blastDet[f"p{aln + 1}"], 'i')
                    bamLine.set_tag(f"l{aln + 1}", blastDet[f"l{aln + 1}"], 'i')
                    bamLine.set_tag("YB", "True", 'Z')
            except Exception:
                print(blastDet)
                raise
            if "am" in blastDet:
                bamLine.set_tag("am", blastDet["am"], 'i')
        else:
            bamLine.set_tag("t0",-4, 'i')
            bamLine.set_tag("am",4,'i')
        if bamLine.query_name == "TGGATACGGTTGATGGCGAAGTGGTTCTCATG":
            print(bamLine.tags)
        outBam.write(bamLine)
        records += 1
    outBam.close()
    print(records)

if __name__ == "__main__":
    main()