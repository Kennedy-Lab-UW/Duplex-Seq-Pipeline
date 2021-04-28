import pysam
from argparse import ArgumentParser


def main():
    parser = ArgumentParser()
    parser.add_argument('prefix')
    parser.add_argument('inNonAmbigFile')
    parser.add_argument('inAmbigFile')
    o = parser.parse_args()

    # check non-ambiguous reads
    nonAmbigFile = pysam.AlignmentFile(o.inNonAmbigFile, 'rb')
    nonAmbigCounter = {0: 0, 1: 0, 2: 0, 3: 0, 4: 0, 5: 0, -1: 0}
    for read in nonAmbigFile:
        if read.has_tag("am"):
            nonAmbigCounter[read.get_tag("am")] += 1
        else:
            nonAmbigCounter[-1] += 1
    nonAmbigFile.close()

    # check ambiguous reads
    ambigFile = pysam.AlignmentFile(o.inAmbigFile, 'rb')
    ambigCounter = {0: 0, 1: 0, 2: 0, 3: 0, 4: 0, 5: 0, -1: 0}
    for read in ambigFile:
        if read.has_tag("am"):
            ambigCounter[read.get_tag("am")] += 1
        else:
            ambigCounter[-1] += 1
    ambigFile.close()

    # write outputs
    outFile = open(f"{o.prefix}_ambiguity_counts.txt", 'w')
    # Define ambiguity classes
    outFile.write(
        f"# AMBIGUITY CLASSES: \n"
        f"#   -1: Missing 'am' tag\n"
        f"#   0: BLAST has 1 best match, matches bwa position\n"
        f"#   1: BLAST has 1 best match, does not match bwa position\n"
        f"#   2: BLAST has 2+ best matches, correct species\n"
        f"#   3: BLAST has 2+ best matches, at least one incorrect species\n"
        f"#   4: BLAST did not attempt alignment, or no blast matches\n"
        f"#   5: Not BLASTed\n"
    )

    # nonAmbiguous read counts
    outFile.write("\n# Non-ambiguous file ambiguity class counts: \n")
    outFile.write("Class\tCounts\n")
    for xIter in (0, 1, 2, 3, 4, 5, -1):
        outFile.write(f"{xIter}\t{nonAmbigCounter[xIter]}\n")

    # ambiguous read counts
    outFile.write("\n# Ambiguous file ambiguity class counts: \n")
    outFile.write("Class\tCounts\n")
    for xIter in (0, 1, 2, 3, 4, 5, -1):
        outFile.write(f"{xIter}\t{ambigCounter[xIter]}\n")
    outFile.close()


if __name__ == "__main__":
    main()
