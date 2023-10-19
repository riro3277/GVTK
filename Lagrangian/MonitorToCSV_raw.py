monitor_path = "/mnt/d/lagrangiantest/PoiseuilleFlow/binVC100/setForceMonitor.mon"
csv_path = "/mnt/d/lagrangiantest/PoiseuilleFlow/binVC100/setForceVC_100U_1P.csv"


file1 = open(monitor_path, "r")
file2 = open(csv_path, "w")
file2.close()
file2 = open(csv_path, "a")
lines = file1.readlines()
tstep = 0.00005
timearr = []
for line in lines:
    p=0
    all = line.split(" ")
    time = line.split(" ")[1]
    forcename = line.split(" ")[2].replace(":", "")
    file2.write("p,0,")
    # if forcename == "Slip":
    #     pass
    # else:
    #     forcex = line.split(" ")[3]
    #     forcey = line.split(" ")[4]
    #     forcez = line.split(" ")[5]
    #     force = line.split(" ")[3]
    #     printline = "p,{},t,{},{},{},{},{}\n".format(p, time, forcename, forcex, forcey, forcez)
    #     # printline = "p,{},t,{},{},{}\n".format(p, time, forcename, force)
    #     file2.write(printline)
    for i in all:
        file2.write(i)
        file2.write(",")
    # file2.write("\n")




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




# p = 97
# lines = file1.readlines()
# p = 0
# line_count = 1
# for index, line in enumerate(lines):
#
#     if index>0:
#         if index % 5 < 1e-10:
#             p+=1
#     if p == 97:
#         p=0
#     time = line.split(" ")[1]
#     forcename = line.split(" ")[2].replace(":", "")
#     force = line.split(" ")[3]
#
#     printline = "p,{},t,{},{},{}\n".format(p, time, forcename, force)
#     file2.write(printline)
#
#
#     print(index)
