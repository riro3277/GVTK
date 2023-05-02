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
#Running in test mode automatically passes the input-tracer-analytical-standard.dat as the input file to run the cylindrical_flow simulation
if len(sys.argv) == 2:
    if sys.argv[1].strip() == "test":
        executeMode = sys.argv[1].strip()
        inputFile = os.path.dirname(os.path.abspath(__file__)) +"/input-tracer-analytical-standard.dat"
#General message to show possible modes of operation and how to run the regression test
elif len(sys.argv) != 3:
    print("Modes of Operation:")
    print("\n", "Collision", "\n", "Tracer", "\n", "Diffuse", "\n", "Inertial", "\n", "Viral", "\n", "Test", "\n")
    print("For test mode, pass test as the first argument of the execute")
    print("Test mode completes regression test against standard results of the cylinder_flow field")
    sys.exit("If not testing, give execution mode as first argument and location of input file as second ")
else:
    inputFile   = sys.argv[2].strip()
    executeMode = sys.argv[1].strip()
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
