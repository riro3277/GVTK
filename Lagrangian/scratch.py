import os
from datetime import date
from testModule import *
# file = open('test.txt', 'a')
# vals = [1,2,3,4,5,6]
# file.write('The vals: ')
# for i in vals:
#     file.write(str(i)+' ')
# file.close()
# p  = 6
# file = open('/mnt/d/lagrangiantest/cylindrical_flow/StandardFileMaster.txt', 'r')
# for line in file:
#     # print(line)
#     vals = line.split(" ")
#     time = vals[0].split(":")[1]
#     force = vals[1].split(":")
#     ForceType = force[0]
#     Forcevals = force[1]
#     # ForceArr = [Forcevals[0], Forcevals[1], Forcevals[2]]
#
#
#     print('abs dirname: ', os.path.dirname(os.path.abspath(__file__)))
#
#     final = vals[1].split(' ')
#     # print(int(final[2])- int(final[1]))
# file.close()

# print(os.path.isfile('test.txt'))
# if os.path.isfile('test.txt') == True:
#     with open('test.txt', 'w'):
#         pass
case = Testing('/mnt/d/lagrangiantest/cylindrical_flow/')
case.Compare()
case.Documentation()
