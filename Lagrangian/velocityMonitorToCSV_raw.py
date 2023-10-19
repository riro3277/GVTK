monitor_path = "/mnt/d/lagrangiantest/PoiseuilleFlow/binVC100/velocityMonitor.mon"
csv_path = "/mnt/d/lagrangiantest/PoiseuilleFlow/binVC100/velocityVC_100U_100P.csv"


file1 = open(monitor_path, "r")
file2 = open(csv_path, "w")
file2.close()
file2 = open(csv_path, "a")
lines = file1.readlines()
tstep = 0.00005
timearr = []



p = 97
p = 0
line_count = 1
for index, line in enumerate(lines):

    if index>0:
        if index % 2 < 1e-10:
            p+=1
    if p == 97:
        p=0
    time = line.split(" ")[1]
    forcename = line.split(" ")[2].replace(":", "")
    force = line.split(" ")[3]

    printline = "p,{},t,{},{},{}\n".format(p, time, forcename, force)
    file2.write(printline)


    print(index)

#Below is used for finding duplicate times in monitor file, make sure to not overwrite/clean files as first lines in this code do that
# for index, line in enumerate(lines):
#
#     if index < 485:
#         currentime=0
#     elif index % 485 < 1e-10:
#         currentime = time
#         timearr.append(time)
#
#     time = line.split(",")[3]
#     # forcename = line.split(",")[2].replace(":", "")
#     # force = line.split(" ")[3]
# uniqueList = []
# duplicateList = []
#
# for i in timearr:
#     if i not in uniqueList:
#         uniqueList.append(i)
#     elif i not in duplicateList:
#         duplicateList.append(i)
#
# print(duplicateList)
#
