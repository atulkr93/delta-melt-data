import math
import sys
import pandas as pd

filename = str(sys.argv[1])
f = open(filename, 'r')
fout = open(filename[:filename.index('.')] + "_clean.csv", 'w')
line = f.readline()
fout.write(line)
for line in f.readlines():
    line_split = line.split(',')
    if math.fabs(float(line_split[1]) / float(line_split[2])) > 4:
        pass
    else:
        fout.write(line)
fout.close()
f.close()
