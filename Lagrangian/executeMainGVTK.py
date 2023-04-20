#---------------------------------------------
# Begin system and user-library module imports
#---------------------------------------------
import sys, os

try:
    import vtk
except ImportError:
    sys.exit("Module vtk with python bindings is required")

try:
    import classGVTKTracerIntegration as L_INT
    import classGVTKDriftDiffuse as L_DIF
    import classGVTKInertialIntegration as L_INR
    import classGVTKVLPIntegration as L_VLP
    import classGVTKCollisionWithMaxeyRiley as L_COL
    import moduleTracerInput as L_INP
    import classGVTKGenericProblem as Prob
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
    print("Modes of Operation:")
    print("\n", "Collision", "\n", "Tracer", "\n", "Diffuse", "\n", "Inertial", "\n", "Viral", "\n", "Test", "\n")
    print("For test mode, give location of -input-tracer-analytical-standard.dat- as the first argument and -test- as the second argument")
    print("-input-tracer-analytical-standard.dat- is located in the GVTK code base, must change paths to ./Cases/cylindrical_flow/ is located on your local machine")
    sys.exit("If not testing, give location of input file as first argument and mode as second arguemnt ")
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
elif executeMode == 'test':
    ensemble = L_COL.GVTKCollision(simInput)
ensemble.runCompute()

#-------------
# exit message
#-------------
print("Computation completed. Exiting!")

#--------------------
# End simulation main
#--------------------
