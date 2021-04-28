# makeSnpList.py
from argparse import ArgumentParser
from VCF_Parser import *

parser = ArgumentParser()
parser.add_argument(
    '-i', 'inSnps',
    action='store',
    dest='inSnps',
    help='A input SNPs VCF file'
)
parser.add_argument(
    '-o', '--outList',
    action='store',
    dest='out',
    help='A name for an output file to append SNPs to.  '
)
o = parser.parse_args()

snpFile = VariantFile(o.inSnps)
outFile = open(o.out, 'a')
for varLine in variantFile:
    outFile.write(
        f"{varLine.chrom}:{varLine.pos}|{varLine.ref}>{varLine.alts[0]}\n"
    )
snpFile.close()
outFile.close()
