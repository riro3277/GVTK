import numpy as np
import sys
#run this to delete repeat time steps in forcemonitor file for each particle. Will create new forcemonitor to be ranin flowVCmonitorfile.py
path = "/mnt/d/lagrangiantest/PoiseuilleFlow/binVC/forceMonitor.mon"
path2 = "/mnt/d/CoW_Clipped/binVC/forceMonitoredit.mon"
file = open(path, "r")
# file2 = open(path2, "w")
lines = file.readlines()
totalp = 10
endt = 0.2
step = 0.0005
savestep = 0.0005




pinterval = int((savestep/step)*5)
totalsteps = int((endt/step) + 1)
row = 0
p = 0
tcalc = 0
# file2.close()
# file2 = open(path2, "a")
indexS =0
badlines = []
for pi in range(totalsteps):
    parr = []
    for index, line in enumerate(lines):
        time = line.split(" ")[1]
        forcename = line.split(" ")[2].replace(":", "")
        force = line.split(" ")[3]


        if tcalc == float(time):
            indexS = index+pinterval
            indexScheck = indexS + 5
            lineS = lines[indexS]
            linesScheck = lines[indexScheck]
            timeScheck = linesScheck.split(" ")[1]
            timeS = lineS.split(" ")[1]
            if float(time) +step  == float(timeS):
                if float(timeScheck) < float(timeS):
                    print(indexS)
                    badlines.append(indexS)


            row+=1
            printline = "p,{},t,{},{},{}".format(p, time, forcename, force)

    tcalc = step*pi
    tcalc = round(tcalc,5)
    # print(tcalc)
# with open(path2, 'w') as fp:
#     # iterate each line
#     for index, line in enumerate(lines):
#         # delete line 5 and 8. or pass any Nth line you want to remove
#         # note list index starts from 0
#         if index not in badlines:
#             fp.write(line)
#         else:
#             pass
