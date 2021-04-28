import datetime
import logging
import sys
from argparse import ArgumentParser
from collections import namedtuple
import statistics
from BedParser import *

Depth_Line = namedtuple(
    "Depth_Line",
    [
        "Chrom",
        "Pos",
        "Ref",
        "DP",
        "Ns"
    ]
)

def str_or_nonetype(inStr):
    if inStr.upper() == "NONE":
        return(None)
    else:
        return(inStr)

def main():
    # Parse in arguments
    parser = ArgumentParser()
    parser.add_argument(
        '-i','--inFile', 
        action="store", 
        dest="in_file",
        default=None, 
        help=("The input depth file.  "
              "If 'None', defaults to stdin.  [%(default)s]"))
    parser.add_argument(
        '-o','--outFile',
        action="store",
        dest="out_file",
        default=None,
        help=("The output depth summary file.  "
              "If 'None', defaults to stdout.  [%(default)s]"))
    parser.add_argument(
        '-b','--bed_file',
        action="store",
        dest="bed_file",
        required=True,
        help="The bed file to compute stats on.")
    parser.add_argument(
        '-m', '--mask_bed',
        action='store',
        dest='mask_bed',
        help='A bed file with the regions to be masked.',
        type=str_or_nonetype,
        default=None)
    parser.add_argument(
        '--blocks_only',
        action="store_true",
        dest="blocks", 
        help="Use only those sites in blocks for calculating whole-line stats.")
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
    o = parser.parse_args()
    
    cmd=" ".join(sys.argv)
    d = datetime.datetime.today()

    # Set up logging
    numeric_level = getattr(logging, o.logLvl.upper(), None)
    if not isinstance(numeric_level, int):
        raise ValueError('Invalid log level: {o.logLvl}')
    logging.basicConfig(
        format='%(levelname)s: %(message)s', 
        level=numeric_level, 
        )
    logging.info(f"Running DepthSummaryCsv.py on {d} using command:")
    logging.info(cmd)

    # Open bed file
    logging.info(f"Opening bed file {o.bed_file}...")
    in_bed = Bed_File(o.bed_file)
    logging.info("Parsing bed file...")
    # Set up data structures
    bed_dict = []
    for line in in_bed:
        # Add the main line
        bed_dict.append(
            {"region": line,
             "depths": [], 
             "class": "Bed_Line",
             "min": 0,
             "mean": 0,
             "median": 0,
             "max": 0})
        # Add subregions, if they differ from the main region
        subregs = line.get_subregions()
        if subregs[0].samtoolsStr() != line.samtoolsStr():
            for block in subregs:
                bed_dict.append(
                    {"region": block,
                     "depths": [], 
                     "class": "Bed_Block",
                     "min": 0,
                     "mean": 0,
                     "median": 0,
                     "max": 0})

    # Open input file
    if o.in_file is None:
        logging.info("Using input from StdIn")
        f_in = sys.stdin
    else:
        logging.info(f"Using input from {o.in_file}")
        f_in = open(o.in_file,'r')
    # Open masking bed:
    if o.mask_bed is not None:
        mask_bed = Bed_File(o.mask_bed)
        mask_regions = [x for x in mask_bed]
    else:
        mask_bed = False
        mask_regions = []
    # Iterate through the input file
    logging.info("Processing input file...")
    for i, lIn in enumerate(f_in):
        if i % 10000 == 0:
            logging.info(f"Processed {i} lines...")
        if lIn[0] != "#":
            line = Depth_Line(*lIn.strip().split())
            # Iterate through the regions and blocks
            for regIter in bed_dict:
                # check if the line is masked:
                mask_line = False
                for mask_iter in mask_regions:
                    if mask_region.contains(line.Chrom, int(line.Pos) - 1):
                        mask_line = True
                # Count the line in any region that contains it
                if (not mask_line and 
                        regIter["region"].contains(line.Chrom, 
                                                   int(line.Pos) - 1,
                                                   o.blocks)):
                    regIter["depths"].append(int(line.DP))

    # Close input file
    if o.in_file is not None:
        f_in.close()

    # Open output file
    if o.out_file is None:
        logging.info("Writing output to StdOut")
        f_out = sys.stdout
    else:
        logging.info(f"Writing output to {o.out_file}")
        f_out = open(o.out_file, 'w')

    # Write header line
    logging.info("Writing output file...")
    f_out.write(
        f"##Input file: {o.in_file}\n"
        f"##Target bed file: {o.bed_file}\n"
        f"##Masking bed file: {o.mask_bed}\n")
    if o.blocks:
        f_out.write("##Blocks only\n")
    f_out.write(
        "#NAME,"
        "CHROM,"
        "START_POS,"
        "END_POS,"
        "TYPE,"
        "MIN,"
        "MEAN,"
        "MEDIAN,"
        "MAX\n")

    # Calculate average depths and write output
    for regIter in bed_dict:
        if len(regIter["region"]) > 0 and len(regIter["depths"]) > 0:
            #regIter["depths"].extend([0 for x in range(
            #    len(regIter["region"])-len(regIter["depths"]))])
            regIter["min"] = min(regIter["depths"])
            regIter["max"] = max(regIter["depths"])
            regIter["median"] = statistics.median(regIter["depths"])
            regIter["mean"] = statistics.mean(regIter["depths"])

        # Write line to output file
        f_out.write(
            f"{regIter['region'].name},"
            f"{regIter['region'].chrom},"
            f"{regIter['region'].startPos + 1},"
            f"{regIter['region'].endPos},"
            f"{regIter['class']},"
            f"{regIter['min']},"
            f"{regIter['mean']},"
            f"{regIter['median']},"
            f"{regIter['max']}\n")

    # Close output file
    if o.out_file is not None:
        f_out.close()

    # DONE
    logging.info("DONE")

if __name__ == "__main__":
    main()