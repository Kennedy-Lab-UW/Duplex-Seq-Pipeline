import sys
iterations = 0
hits = 0
reportedIteration = True
sys.stderr.write("Monitoring Blast Output:\n")
for line in sys.stdin:
    if "<Iteration>" in line:
        iterations += 1
        reportedIteration = False
    if "<Hit>" in line:
        hits += 1
    if iterations % 1000 == 0 and not reportedIteration:
        sys.stderr.write(f"{iterations} reads processed\n\t{hits} hits\n")
        reportedIteration = True
    sys.stdout.write(line)
