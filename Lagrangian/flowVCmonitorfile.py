import numpy as np
import sys
#will create csv file that arranges particle forces per particle per time. run this using the file created in flowVCmonitorfileeditor.py
path = "/mnt/d/lagrangiantest/PoiseuilleFlow/binVC100/forceMonitor.mon"
path2 = "/mnt/d/lagrangiantest/PoiseuilleFlow/binVC100/forcesVC.csv"
file = open(path, "r")
file2 = open(path2, "w")
lines = file.readlines()
totalp = 9
endt = 0.2
step = 0.00005
savestep = 0.0005


pinterval = int((savestep/step)*5)
totalsteps = int((endt/step) + 1)
row = 0
p = 0
tcalc = 0
file2.close()
file2 = open(path2, "a")
for pi in range(totalsteps):
    parr = []
    for index, line in enumerate(lines):
        time = line.split(" ")[1]
        forcename = line.split(" ")[2].replace(":", "")
        force = line.split(" ")[3]
        if index % pinterval < 1e-10:
            if index < pinterval:
                p = 0
            else:
                p += 1
        if p == totalp:
            p = 0


        if tcalc == float(time):
            row+=1
            printline = "p,{},t,{},{},{}\n".format(p, time, forcename, force)
            file2.write(printline)




    tcalc = step*pi
    tcalc = round(tcalc,5)
    print(tcalc)
