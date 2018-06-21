import sys

with open(sys.argv[1]) as a:
    for line in a:
        print(line.rstrip())

