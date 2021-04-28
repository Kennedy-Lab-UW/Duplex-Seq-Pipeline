import datetime
import logging
import sys
from argparse import ArgumentParser
from collections import Counter
import pysam
from BedParser import *

def getParams():
    parser = ArgumentParser(
        description=(
            f"Make a mutpos vcf file from a post-procesing duplex "
            f"sequencing BAM file.  Note that this program will fail if"
            f" the maximum depth in your bam file is > 8000.  "
            )
        )
    parser.add_argument(
        '-i', '--inBam', 
        action='store', 
        dest='inBam', 
        help='An imput bam file. If None, defaults to stdin. [%(default)s]', 
        default=None
        )
    parser.add_argument(
        '-b', '--inBed', 
        action='store', 
        dest='inBed', 
        help='An input bed file. If None, processes all positions. [%(default)s]', 
        default=None
        )
    parser.add_argument(
        '-f', '--inFasta', 
        action='store', 
        dest='in_fasta', 
        help='The reference genome fasta file.  ', 
        required=True
        )
    parser.add_argument(
        '-o', '--outfile', 
        action='store', 
        dest='out_file', 
        help='A filename for the output file.  If None, outputs to stdout.  [%(default)s]', 
        default=None
        )
    parser.add_argument(
        "--round", 
        action="store", 
        type=int, 
        dest="round", 
        help="How many digits to round frequencies to.", 
        default=4
        )
    parser.add_argument(
        '--logLevel', 
        action="store", 
        dest="logLvl", 
        default="Info",
        help=(f"Identification for how much information gets output. "
              f"Acceptable levels are: 'DEBUG', 'INFO', 'WARNING', "
              f"'ERROR', and 'CRITICAL'.  "
              )
        )
    parser.add_argument(
        '--samp_name', 
        action='store', 
        dest='sampName',
        default=None,
        help='A name for the sample in the output VCF file.  '
        )
    return(parser.parse_args())

class MutposEngine: 
    def __init__(
            self, 
            in_bam, 
            out_file, 
            in_fasta, 
            sampName=None, 
            in_bed=None, 
            cmd=" ".join(sys.argv), 
            NsOut=False
            ):
        self.inBam = pysam.AlignmentFile(in_bam, "rb")
        self.outputNs = NsOut
        
        
        # open output file
        self.outDepth = open(out_file, 'w')
        self.outDepth.write("#Chrom\tPos\tRef\tDP\tNs\n")
        # open fasta file
        self.inFasta = pysam.FastaFile(in_fasta)
        # set up constants
        self.rmNumTable = {ord(c): None for c in '1234567890'}
        if in_bed is None:
            self.myBed = None
        else:
            self.myBed = Bed_File(in_bed)
    
    def processLines(self):
        linesProcessed = 0
        if self.myBed is None:
            for pileup_column in self.inBam.pileup(
                    fastafile=self.inFasta, 
                    truncate=True, 
                    stepper="nofilter", 
                    max_depth=1000000, 
                    min_base_quality=0
                    ):
                myLines = self.makeVcfLine(pileup_column)
                logging.debug(f"{str(myLines)}")
                self.outDepth.write(
                    f"{myLines['Chrom']}\t"
                    f"{myLines['Pos']}\t"
                    f"{myLines['Ref']}\t"
                    f"{myLines['DP']}\t"
                    f"{myLines['Ns']}\n"
                    )
                linesProcessed += 1
                if linesProcessed % 1000 == 0:
                    logging.info(f"{linesProcessed} lines processed...")
        else:
            for myRegion in self.myBed:
                for pileup_column in self.inBam.pileup(
                        reference=myRegion.chrom, 
                        start=myRegion.startPos, 
                        end=myRegion.endPos,
                        fastafile=self.inFasta, 
                        truncate=True, 
                        stepper="nofilter", 
                        max_depth=1000000, 
                        min_base_quality=0
                        ):
                    
                    myLines = self.makeVcfLine(pileup_column)
                    logging.debug(f"MyLines are {str(myLines)}")
                    self.outDepth.write(
                        f"{myLines['Chrom']}\t"
                        f"{myLines['Pos']}\t"
                        f"{myLines['Ref']}\t"
                        f"{myLines['DP']}\t"
                        f"{myLines['Ns']}\n"
                        )
                    linesProcessed += 1
                    if linesProcessed % 1000 == 0:
                        logging.info(f"{linesProcessed} lines processed...")
    
    
    def makeVcfLine(self, pileup_column):
        myReads = Counter(
            [x.translate(self.rmNumTable).upper() 
             for x in pileup_column.get_query_sequences(
                add_indels=True
                )
             ]
            )
        myChrom = pileup_column.reference_name
        myPos = pileup_column.reference_pos + 1
        # I'll need to pull the reference base from the fasta file
        
        myRefBase = self.inFasta.fetch(reference=myChrom, start=myPos-1, end=myPos).upper()
        
        logging.debug(f"Refernce is: {myRefBase}")
        myNCount = sum([myReads[x] for x in myReads if x[0] == "N"])
        myTotCount = sum([myReads[x] for x in myReads])
        DepthLine={
            "Chrom": myChrom, 
            "Pos": myPos, 
            "Ref": myRefBase, 
            "DP": myTotCount - myNCount, 
            "Ns": myNCount
            }
        return(DepthLine)

def main():
    o = getParams()
    numeric_level = getattr(logging, o.logLvl.upper(), None)
    if not isinstance(numeric_level, int):
        raise ValueError('Invalid log level: {o.logLvl}')
    logging.basicConfig(
        format='%(levelname)s: %(message)s', 
        level=numeric_level, 
        )
    # Start engine
    logging.info("Starting Engine")
    myEngine = MutposEngine(o.inBam, 
                            o.out_file, 
                            o.in_fasta, 
                            o.sampName, 
                            o.inBed
                            )
    logging.info("Processing Lines")
    myEngine.processLines()

if __name__ == "__main__":
    main()
