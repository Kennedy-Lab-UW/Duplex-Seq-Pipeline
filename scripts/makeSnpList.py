#makeSnpList.py
from argparse import ArgumentParser
# ~ import scipy
# ~ from scipy import stats
from VCF_Parser import *

parser = ArgumentParser()
parser.add_argument(
    '-i', 'inSnps', 
    action='store', 
    dest='in', 
    help='A input SNPs VCF file'
    )
parser.add_argument(
    '-o','--outList',
    action='store', 
    dest='out', 
    help='A name for an output file to append SNPs to.  '
    )
o = parser.parse_args()

snpFile = VariantFile(o.in)
outFile = open(o.out, 'a')
for varLine in variantFile:
    outFile.write(
        f"{varLine.chrom}:{varLine.pos}|{varLine.ref}>{varLine.alts[0]}\n"
        )
snpFile.close()
outFile.close()

