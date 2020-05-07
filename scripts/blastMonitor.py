"""blastMonitor.py
This program will monitor progress of a BLAST run and report on number
of hits returned as the run continues.

Takes input from standard in, writes outputs to standard out and
standard error.
"""
import sys


def main():
    """main function"""
    iterations = 0
    hits = 0
    reported_iteration = True
    sys.stderr.write("Monitoring Blast Output:\n")
    for line in sys.stdin:
        if "<Iteration>" in line:
            iterations += 1
            reported_iteration = False
        if "<Hit>" in line:
            hits += 1
        if iterations % 1000 == 0 and not reported_iteration:
            sys.stderr.write(f"{iterations} reads processed\n\t{hits} hits\n")
            reported_iteration = True
        sys.stdout.write(line)


if __name__ == "__main__":
    main()
