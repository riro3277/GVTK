import sys, os, vtk
import numpy as np
import time
import random

from numpy import linalg as la

try:
    from classGVTKGrid import *
    from classGVTKLagrangianData import *
    from classGVTKGenericProblem import *
    from testModule import *
except:
    sys.exit('Could not import user defined modules')


tracer_rootDir = "/mnt/d/lagrangiantest/PoiseuilleFlow/Output_nc/"
tracer_prefix = "tracer_nc_
type = ".vtk"
start = 0
end = 60000
delta = 200

for i in range(start, end+delta, delta):
    file = tracer_rootDir + tracer_prefix + str(i) + type
    ParticleData = GVTKLagrangianData(a_File = file)
    for p in range(ParticleData.numParticles):
        posP         = ParticleData.getX(p)
