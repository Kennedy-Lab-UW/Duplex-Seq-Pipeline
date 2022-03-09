import datetime
import logging
import sys
from argparse import ArgumentParser
from collections import Counter
import pysam
from BedParser import *
from VCF_Parser import SamHeaderToVcfHeader


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
    return parser.parse_args()


class Error(Exception):
    """Base class for exceptions in this module."""
    pass


class VCF_Error(Error):
    def __init__(self, message):
        self.message = f"VCF_ERROR: {message}"
        sys.stderr.write(self.message)


class VCF_Line:
    def __init__(self,
                 chrom=".",
                 pos=".",
                 ref=".",
                 alt=".",
                 refCnt=0,
                 altCnt=0,
                 nCnt=0,
                 depth=0
                 ):
        self.chrom = chrom
        self.pos = pos
        self.id = '.'

        self.qual = '.'
        self.filter = '.'
        self.info = '.'
        self.format = ["AD", "DP", "AF", "NC"]
        # process indels properly
        myRef = ref
        if '+' in alt:
            myAlt = alt.replace('+', '')
        elif '-' in alt:
            myAlt = alt[0]
        else:
            myAlt = alt
        self.alleles = [myRef.upper(), myAlt.upper()]
        self.data = {
            "AD": {myRef.upper(): refCnt, myAlt.upper(): altCnt},
            "DP": depth,
            "AF": {
                myRef.upper(): f"{refCnt / depth:.2E}",
                myAlt.upper(): f"{altCnt / depth:.2E}"
            },
            "NC": nCnt
        }
        self.skippedReads = 0

    def formatStr(self):
        return ":".join(self.format)

    def dataStr(self):
        retStr = (
            f"{','.join([str(self.data['AD'][x]) for x in self.alleles])}"
            f":{self.data['DP']}"
            f":{','.join([str(self.data['AF'][x]) for x in self.alleles])}"
            f":{self.data['NC']}"
        )
        return retStr

    def __str__(self):
        if len(self.alleles) == 1:
            return ""
        else:
            return (f"\n{self.chrom}"
                    f"\t{self.pos}"
                    f"\t{self.id}"
                    f"\t{self.alleles[0]}"
                    f"\t{','.join(self.alleles[1:])}"
                    f"\t{self.qual}"
                    f"\t{self.filter}"
                    f"\t{self.info}"
                    f"\t{self.formatStr()}"
                    f"\t{self.dataStr()}"
                    )


class VCF_File:
    def __init__(self,
                 inFileName,
                 inFileMode,
                 sampName=None,
                 headStart=None
                 ):
        self.mode = inFileMode
        if inFileMode == 'r':
            raise VCF_Error("Reading VCF isn't supported right now.")
        elif inFileMode == 'w':
            if inFileName is None:
                self.file = sys.stdout
            else:
                self.file = open(inFileName, 'w')
            d = datetime.datetime.today()
            if headStart is None:
                self.file.write(
                    f'##fileformat=VCFv4.3\n'
                    f'##fileDate={d.year}{str(d.month).zfill(2)}{str(d.day).zfill(2)}\n'
                    f'##VCF_MutposVersion=v0.2.0\n'
                    f'##VCF_MutposCommand=VCF_Mutpos test command'
                )
            else:
                self.file.write("\n".join(str(headStart).split('\n')[:-2]))
            self.file.write(
                f'\n##ALT=<ID=*,Description="Represents allele(s) other than observed.">\n'
                f'##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">\n'
                f'##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allele Depth">\n'
                f'##FORMAT=<ID=AF,Number=R,Type=Float,Description="Allele frequency">\n'
                f'##FORMAT=<ID=NC,Number=1,Type=Integer,Description="N Count">\n'
                f'#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t'
            )
            if sampName is None:
                self.file.write("SAMPLE")
            else:
                self.file.write(f"{sampName.upper()}")
        else:
            raise VCF_Error(f"invalid mode: {inFileMode}")

    def readline(self):
        if self.mode == 'r':
            return self.file.readline()
        else:
            logging.error(f"This file is in write mode")
            raise ValueError(f"This file is in write mode")

    def write(self, lineToWrite):
        if self.mode == 'w':
            self.file.write(str(lineToWrite))
        else:
            logging.error(f"This file is in read mode")
            raise ValueError(f"This file is in read mode")

    def makeline(self, chrom=".", pos=".", ref="."):
        return ()


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
        # Create vcf header
        vcfHeadStart = SamHeaderToVcfHeader(
            self.inBam.header,
            "Sample",
            "VCF_Mutpos",
            "0.3.0",
            cmd
        )
        self.outFile = VCF_File(out_file, 'w', sampName, vcfHeadStart)

        self.outDepth = open(f"{out_file}_depth.txt", 'w')
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
                for myLine in myLines[0]:
                    self.outFile.write(myLine)
                self.outDepth.write(
                    f"{myLines[1]['Chrom']}\t"
                    f"{myLines[1]['Pos']}\t"
                    f"{myLines[1]['Ref']}\t"
                    f"{myLines[1]['DP']}\t"
                    f"{myLines[1]['Ns']}\n"
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
                    for myLine in myLines[0]:
                        self.outFile.write(myLine)
                    self.outDepth.write(
                        f"{myLines[1]['Chrom']}\t"
                        f"{myLines[1]['Pos']}\t"
                        f"{myLines[1]['Ref']}\t"
                        f"{myLines[1]['DP']}\t"
                        f"{myLines[1]['Ns']}\n"
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

        myRefBase = self.inFasta.fetch(reference=myChrom, start=myPos - 1, end=myPos).upper()

        logging.debug(f"Refernce is: {myRefBase}")
        myVcfLines = []
        myRefCount = myReads[myRefBase]
        myNCount = sum([myReads[x] for x in myReads if x[0] == "N"])
        myTotCount = sum([myReads[x] for x in myReads])
        for readTypeKey in myReads:
            if readTypeKey.upper() != myRefBase and myRefBase != "N":
                # make a VCF line for this sample
                # If alt is "N", depth is treated differently
                if readTypeKey[0] == 'N':
                    if self.outputNs:
                        myVcfLines.append(
                            VCF_Line(
                                chrom=myChrom,
                                pos=myPos,
                                ref=myRefBase,
                                alt=readTypeKey,
                                refCnt=myRefCount,
                                depth=myTotCount,
                                altCnt=myReads[readTypeKey],
                                nCnt=myNCount)
                        )
                elif readTypeKey != "*":
                    if '-' in readTypeKey:
                        myRefBPs = self.inFasta.fetch(
                            reference=myChrom,
                            start=myPos - 1,
                            end=myPos + len(readTypeKey[2:])
                        ).upper()
                    else:
                        myRefBPs = myRefBase
                    logging.debug(
                        f"{readTypeKey}, "
                        f"{myReads[readTypeKey]}, "
                        f"{myRefBPs}, "
                        f"{myReads[myRefBPs]}"
                    )
                    myVcfLines.append(
                        VCF_Line(
                            chrom=myChrom,
                            pos=myPos,
                            ref=myRefBPs,
                            alt=readTypeKey,
                            refCnt=myRefCount,
                            depth=myTotCount - myNCount,
                            altCnt=myReads[readTypeKey],
                            nCnt=myNCount
                        )
                    )
        DepthLine = {
            "Chrom": myChrom,
            "Pos": myPos,
            "Ref": myRefBase,
            "DP": myTotCount - myNCount,
            "Ns": myNCount
        }
        return [myVcfLines, DepthLine]


def delRepFunc(inStr):
    if len(inStr.group(0)) == 2:
        return inStr.group(0)[0]
    else:
        raise IOError(f"Bad input string {inStr}")


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
