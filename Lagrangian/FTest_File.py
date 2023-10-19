import numpy as np
import sys

baseDir = "/mnt/d/FTest/"
name = "Forces"
ext = ".csv"
start_ind = 0
end_ind = 2000
comb_file = baseDir + "FTest_All.csv"
f = open(comb_file, "w")
f.close()
file2 = open(comb_file, "a")
for i in range(start_ind, end_ind+1):
    file = baseDir + name + str(i+1) + ext
    read = open(file, "r")
    lines = read.readlines()


    for index, line in enumerate(lines):
        if index == 0:
            pass
        else:
            file2.write(line)

    # sys.exit()
