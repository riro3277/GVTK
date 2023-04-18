#---------------------------------------------
# Begin system and user-library module imports
#---------------------------------------------
import sys, os
from mpi4py import MPI

try:
    import vtk
except ImportError:
    sys.exit("Module vtk with python bindings is required")

try:
    import classGVTKTracerIntegration as L_INT
    import classGVTKDriftDiffuse as L_DIF
    import classGVTKInertialIntegration as L_INR
    import classGVTKVLPIntegration as L_VLP
    import classGVTKCollision as L_COL
    import moduleTracerInput as L_INP
except ImportError:
    sys.exit("One or more of library modules not accessible")

#-------------------------------------------
# End system and user-library module imports
#-------------------------------------------

#----------------------
# Begin simulation main
#----------------------

#-----------------------------
# parse command-line arguments
#-----------------------------
# if len(sys.argv) != 3:
#     sys.exit("Need input file name and execute mode as argument")
# else:
#     inputFile   = sys.argv[1].strip()
#     executeMode = sys.argv[2].strip()

comm = MPI.COMM_WORLD
rank = comm.rank

# print('rank:', rank)
# print('type:', type(rank))
# sys.exit()

inputFile = 'dats/case-input-'+str(rank)+'.dat'
executeMode = 'diffuse'

#-------------------------------------------------------------
# instantiate simulation input object and configure input data
#-------------------------------------------------------------
simInput    = L_INP.SimInputs(inputFile)
simInput.readInputFile()

#--------------------------------------------------------------
# instantiate Lagrangian physics object and execute computation
#--------------------------------------------------------------
if executeMode == 'tracer':
    ensemble    = L_INT.GVTKTracerIntegration(simInput)
elif executeMode == 'diffuse':
    ensemble    = L_DIF.GVTKDriftDiffuse(simInput)
elif executeMode == 'inertial':
    ensemble    = L_INR.GVTKInertial(simInput)
elif executeMode == 'viral':
    ensemble    = L_VLP.GVTKViral(simInput)
elif executeMode == 'collision':
    ensemble    = L_COL.GVTKCollision(simInput)

ensemble.runCompute()

#-------------
# exit message
#-------------
print("Computation completed. Exiting!")

#--------------------
# End simulation main
#--------------------
