#---------------------------------------------
# Begin system and user-library module imports
#---------------------------------------------
import sys, os

try:
    import vtk
except ImportError:
    sys.exit("Module vtk with python bindings is required")

try:
    import classTracers as L_TRC
    import tracerInput as L_INP
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
if len(sys.argv) != 3:
    printUsage()
    sys.exit("Need input file name and execute mode as argument")
else:
    inputFile   = sys.argv[1].strip()
    executeMode = sys.argv[2].strip()

#-------------------------------------------------------------
# instantiate simulation input object and configure input data
#-------------------------------------------------------------
simInput    = L_INP.SimInputs(inputFile)
simInput.readInputFile()

#--------------------------------------------------------------
# instantiate Lagrangian physics object and execute computation
#--------------------------------------------------------------
ensemble    = L_TRC.SimpleTracer(simInput)
ensemble.runCompute()

#-------------
# exit message
#-------------
print("Computation completed. Exiting!")

#--------------------
# End simulation main
#--------------------
