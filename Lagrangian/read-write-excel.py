# -*- coding: utf-8 -*-
"""
Created on Wed Nov 16 11:24:44 2022

@author: Sreeparna
"""

import matplotlib.pyplot as plt
import pandas as pd
import xlwings as xw
import seaborn as sns
import numpy as np

# dataset-1 Positive LNH
# Specifying a sheet
ws = xw.Book("forceMonitor_Drag.xlsx").sheets['Sheet2']

# Selecting data from
t = ws.range("B1:B6").value
drag = ws.range("D1:D6").value
Tdrag = 0
stoptime = 0

for i in range(len(t)):
    t1 = t[i]
    while (t1 == stoptime):
        Tdrag = Tdrag + drag[i]
    Davg = Tdrag/6.0
    stoptime = stoptime + 0.00005
    
print(Davg)

# for i in range(len(y1)):
#     for j in range(len(y1)):
#         if i!=j:
#             if abs(y1[i]-y1[j]) < eps:
#                 overlap.append(i)
# x1 = [1 for i in range(len(x1))]