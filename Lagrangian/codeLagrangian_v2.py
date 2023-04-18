#-------------------------------------------------------------------------------------------------
# This is a single file implementation of a script to compute Lagrangian fields from analytical or 
# advective flow field data, which is designed based on a simple code design philosophy that 
# heavily makes use of the VTK library, numpy, and standalone Fortran modules for all major 
# geometry and compute operations.
#
# The objective of this module is to:
# - perform an FTLE, elongational stretch, and/or longitudinal strain computation for analytical
#   and advective flow data/fields
# - provide a single script for future development and testing of FTLE, elongational stretch,
#   and/or longitudinal strain related capabilities.
#
# Additional detailed features to be included in future versions.
#
# Author:       Zachariah Irwin
# Institution:  University of Colorado, Boulder
# Last Edit:    March 2020
#-------------------------------------------------------------------------------------------------
import sys, os, time

try:
  import numpy as np
except ImportError:
  print("Module Error! Numpy not installed.")
  sys.exit()

try:
  import vtk
except ImportError:
  print("Module Warning! VTK module with Python bindings not installed.")
  sys.exit()

try:
  import tracerInput as INPUT
except ImportError:
  print("Module Warning! Could not import module 'tracerInput', check configuration.")
  sys.exit()

try:
  import preProcess as PRE
except ImportError:
  print("Module Warning! Could not import module 'preProcess', check configuration.")
  sys.exit()

try:
  import postProcess as POST
except ImportError:
  print("Module Warning! Could not import module 'postProcess', check configuration.")
  sys.exit()

try:
  import mesh
except ImportError:
  print("Module Warning! Compiled Fortran 'mesh' module not available.")
  sys.exit()

try:
  import integration
except ImportError:
  print("Module Warning! Compiled Fortran 'integration' module not available.")
  sys.exit()

try:
  import main as MAIN
except ImportError:
  print("Module Warning! Compiled Fortran 'main' module not available.")
  sys.exit()


#--------------------------------
# Begin compute script definition
#--------------------------------
if __name__=="__main__":

  start = time.time()

  EPS   = 1.0e-12

  #----------------------------
  # Parse command line argument
  #----------------------------
  if len(sys.argv) != 2:
    sys.exit("Need Input Filename As An Argument")

  inputFile = sys.argv[1].strip()

  #--------------------------
  # Set-up problem input data
  #--------------------------
  inputData    = INPUT.SimInputs(inputFile)
  inputData.readInputFile()

  #------------------------------
  # Get computational information
  #------------------------------ 
  flowModel = inputData.getFlowModel()
  params    = inputData.getLagrangianParams()

  lagFlags     = np.zeros((3), dtype=np.int32, order='F')
  lagFlags[0]  = params['FTLE']
  lagFlags[1]  = params['Stretch']
  lagFlags[2]  = params['Strain']

  unitVec_1  = inputData.m_UnitVec1
  unitVec_2  = inputData.m_UnitVec2
  unitVec_3  = inputData.m_UnitVec3

  unitVectors       = np.zeros((3,3), dtype=np.float64, order='F')
  unitVectors[0,:]  = unitVec_1
  unitVectors[1,:]  = unitVec_2
  unitVectors[2,:]  = unitVec_3

  #---------------------------
  # Create integration windows
  #---------------------------
  outputTStart     = inputData.m_OutTimeStart
  outputTInterval  = inputData.m_OutTimeInterval
  outputIDStart    = inputData.m_OutIndexStart
  outputIDStop     = inputData.m_OutIndexStop
  outputIDInterval = inputData.m_OutIndexDelta

  if outputIDInterval == 0:
    outputCount = 0
  else:
    outputCount = (outputIDStop - outputIDStart)/outputIDInterval

  integrationWindows = np.linspace(outputTStart, outputTStart + outputCount*outputTInterval, \
                                   int(outputCount+1), endpoint=True)
  
  #-------------------------
  # Get advection parameters
  #-------------------------
  startTime    = inputData.getSimulationStartTime()
  endTime      = inputData.getSimulationStopTime()
  timeDelta    = inputData.getIntegrationTimeStep()
  dataDelta    = inputData.getDataTimeDelta()
  intScheme    = inputData.getIntegrationScheme()

  if intScheme == 'rk4':
    intScheme = 1
  elif intScheme == 'euler':
    intScheme = 2

  if (flowModel == 0):

    #-----------------------------------
    # Parse velocity file ID information
    #-----------------------------------
    fileIn_IDTiming     = np.zeros((3), dtype=np.int32, order='F')
    fileIn_IDTiming[0]  = inputData.m_DataIndexStart
    fileIn_IDTiming[1]  = inputData.m_DataIndexStop
    fileIn_IDTiming[2]  = inputData.m_DataIndexDelta

    #---------------------------------------
    # Parse velocity file timing information
    #---------------------------------------
    fileIn_TTiming      = np.zeros((3), dtype=np.float64, order='F')
    fileIn_TTiming[0]   = inputData.m_DataTimeStart
    fileIn_TTiming[1]   = inputData.m_DataTimeStop
    fileIn_TTiming[2]   = inputData.m_DataTimeDelta

    #-----------------------------------------
    # Convert velocity files from .vtu to .bin
    #-----------------------------------------
    inPrefix = inputData.m_FlowDirectory + inputData.m_FlowFileTag
    print("Generating velocity files...")

    fileIndex = fileIn_IDTiming[0]
    timeIndex = 0
    timeStamp = fileIn_TTiming[0]

    numDataFiles = int((fileIn_IDTiming[1]-fileIn_IDTiming[0])/fileIn_IDTiming[2])+1
    
    for f in range(numDataFiles):
      fileName    = inputData.getFlowDataFileName(a_ID=fileIndex)
      
      PRE.convertVTKDataToBin(fileName, inPrefix, timeIndex, timeStamp, \
              a_DataFieldName = inputData.getVelDataName(), \
              a_LegacyDataType = 'vtu')

      timeIndex   = timeIndex + 1
      fileIndex   = fileIndex + 1
      timeStamp   = timeStamp + fileIn_TTiming[2]
    
    print("Finished generating velocity files.\n")

  #----------------------------------------------
  # Get geometry information (assumes fixed mesh)
  #----------------------------------------------
  flowGeo  = inputData.getStandardFlowDomainGeometryDefinition()
  lagGeo   = inputData.getStandardLagrangianDomainGeometryDefinition()
  numDim   = inputData.getProblemDimension()

  #------------------------------------------------------------------
  # Check to see that Lagrangian domain is constrained by flow domain
  #------------------------------------------------------------------
  if ((lagGeo['X0'] < flowGeo['X0']) or (lagGeo['X1'] > flowGeo['X1']) or \
      (lagGeo['Y0'] < flowGeo['Y0']) or (lagGeo['Y1'] > flowGeo['Y1']) or \
      (lagGeo['Z0'] < flowGeo['Z0']) or (lagGeo['Z1'] > flowGeo['Z1'])):
    print("ERROR: Lagrangian domain must be constrained by flow domain. Check input parameters.")
    sys.exit()

  numX_Lag    = int((lagGeo['X1'] - lagGeo['X0'])/lagGeo['DX']) + 1
  numY_Lag    = int((lagGeo['Y1'] - lagGeo['Y0'])/lagGeo['DY']) + 1
  if inputData.getProblemDimension() < 3:
    numZ_Lag = 1
  else:
    numZ_Lag  = int((lagGeo['Z1'] - lagGeo['Z0'])/lagGeo['DZ']) + 1

  numAlongAxis = np.asarray([numX_Lag, numY_Lag, numZ_Lag], dtype=np.int32, order='F')
  numParticles = numX_Lag*numY_Lag*numZ_Lag

  flowBBox   = np.asarray([flowGeo['X0'], flowGeo['X1'], flowGeo['Y0'], flowGeo['Y1'], flowGeo['Z0'], flowGeo['Z1']], dtype=np.float64, order='F')
  lagBBox    = np.asarray([lagGeo['X0'], lagGeo['X1'], lagGeo['Y0'], lagGeo['Y1'], lagGeo['Z0'], lagGeo['Z1']], dtype=np.float64, order='F')
  bcType     = flowGeo['B0']

  intFlags     = np.zeros((3), dtype=np.int32, order='F')
  intFlags[0]  = numDim
  intFlags[1]  = intScheme
  intFlags[2]  = bcType

  outfileIndex = inputData.m_OutIndexStart
  for integrationWindow in integrationWindows:

    #------------------------------------
    # Initialize integration window times
    #------------------------------------
    T0 = startTime + integrationWindow
    T1 = endTime + integrationWindow

    times     = np.zeros((3), dtype=np.float64, order='F')
    times[0]  = T0
    times[1]  = T1
    times[2]  = timeDelta

    #-----------------------------------------------------------
    # Build NumPy grids (if it is a fixed mesh, built only once)
    #-----------------------------------------------------------
    if inputData.isFixedMesh():

      print("\nBuilding Lagrangian NumPy Grids initiated at %.4fs" %(startTime+integrationWindow))

      Cauchy_Green = np.zeros((numParticles, 3, 3), dtype=np.float64, order='F')
      FTLE_T       = np.zeros((numParticles), dtype=np.float64, order='F')
      FTLE_No_T    = np.zeros((numParticles), dtype=np.float64, order='F')
      Stretch_1    = np.zeros((numParticles), dtype=np.float64, order='F')
      Stretch_2    = np.zeros((numParticles), dtype=np.float64, order='F')
      Stretch_3    = np.zeros((numParticles), dtype=np.float64, order='F')
      Strain_1     = np.zeros((numParticles), dtype=np.float64, order='F')
      Strain_2     = np.zeros((numParticles), dtype=np.float64, order='F')
      Strain_3     = np.zeros((numParticles), dtype=np.float64, order='F')

      ref_dims  = np.asarray([lagGeo['DX'], lagGeo['DY'], lagGeo['DZ']], dtype=np.float64, order='F')

      seedCoordinates = mesh.mesh.seed_points_cartesian(numParticles, ref_dims, lagBBox)

      # Save original configuration of points for mapping/VTK writing    
      referenceCoordinates = np.copy(seedCoordinates, order='F')

      print("Finished Building Lagrangian NumPy grids initiated at %.4fs \n" %(startTime+integrationWindow))

    #---------------------------------
    # Simulated/experimental data sets
    #---------------------------------
    if flowModel == 0: 

      if inputData.isFixedMesh() and integrationWindow == integrationWindows[0]:
        print("Building Cell Locator Maps...")
        locatorObj  = PRE.createCellLocator(inputData.getFlowDataFileName(a_ID=inputData.m_DataIndexStart))
        print("Finished Building Cell Locator Maps \n")

      #-------------------------------------------------------
      # Get mesh data (if it is a fixed mesh, built only once)
      #-------------------------------------------------------
      if inputData.isFixedMesh() and integrationWindow == integrationWindows[0]:
        coordinates, connectivity, adjacency = PRE.getMeshData(inputData.getFlowDataFileName(a_ID=inputData.m_DataIndexStart))
        numNodes      = coordinates.shape[0]
        numCells      = connectivity.shape[0]
      if integrationWindow == integrationWindows[0]:
        print("Generating Cell-IDs for Lagrangian Tracers...")
        cellID_Init   = PRE.getCellIDs_Cartesian(seedCoordinates, locatorObj)
        print("Finished Generating Cell-IDs for Lagrangian Tracers.\n")

      #----------------
      # Run main script
      #----------------
      filePrefix = inPrefix + '.vel'
      cellIDs    = cellID_Init

      print("Now integrating %i tracers from %.4fs to %.4fs" %(numParticles, T0, T1))
      integration.integration.print_integration_parameters(flowModel, bcType, intScheme, T1-T0, timeDelta, flowBBox)

      MAIN.main.lag_fields_mesh(intFlags, lagFlags, numParticles, ref_dims, numAlongAxis, flowBBox, \
                                lagBBox, times, fileIn_IDTiming, fileIn_TTiming, filePrefix, \
                                numNodes, numCells, connectivity, adjacency, coordinates, \
                                unitVectors, FTLE_T, FTLE_No_T, Stretch_1, Stretch_2, Stretch_3, \
                                Strain_1, Strain_2, Strain_3, Cauchy_Green, cellIDs)
      
      #------------------------------------------------------------------
      # At the end of the integration window, write fields into vtk files
      #------------------------------------------------------------------
      print("\nFinished integrating tracers and computing Lagrangian fields.")
      print("\nWriting Lagrangian VTK data...")
      
      #--------------------------------------------------------
      # Flags for which stretch/strain arrays to write (if any)
      #--------------------------------------------------------
      computeDirections = np.ones(3, dtype=np.bool)
      if np.all(unitVectors[0] == 0.0):
        computeDirections[0] = False
      if np.all(unitVectors[1] == 0.0):
        computeDirections[1] = False
      if np.all(unitVectors[2] == 0.0):
        computeDirections[2] = False

      POST.writeLagrangianFieldToVTK(referenceCoordinates, numAlongAxis, FTLE_T, FTLE_No_T,\
                                     Stretch_1, Stretch_2, Stretch_3, \
                                     Strain_1, Strain_2, Strain_3, lagFlags, computeDirections, \
                                     inputData.getLagrangianOutputFile(a_ID=outfileIndex))
      print("Finished writing Lagrangian VTK data into %s" %inputData.getLagrangianOutputFile(a_ID=outfileIndex))
      outfileIndex = outfileIndex + inputData.m_OutIndexDelta
      
    #-----------------------
    # Analytical flow fields
    #-----------------------
    else:
      
      print("Now integrating %i tracers from %.4fs to %.4fs" %(numParticles, T0, T1))
      integration.integration.print_integration_parameters(flowModel, bcType, intScheme, T1-T0, timeDelta, flowBBox)

      MAIN.main.lag_fields_anlyt(flowModel, intFlags, lagFlags, numParticles, ref_dims, numAlongAxis, \
                                 flowBBox, lagBBox, times, unitVectors, FTLE_T, FTLE_No_T, \
                                 Stretch_1, Stretch_2, Stretch_3, Strain_1, Strain_2, Strain_3, Cauchy_Green)

      #------------------------------------------------------------------
      # At the end of the integration window, write fields into vtk files
      #------------------------------------------------------------------
      print("\nFinished integrating tracers and computing Lagrangian fields.")
      print("\nWriting Lagrangian VTK data...")

      #--------------------------------------------------------
      # Flags for which stretch/strain arrays to write (if any)
      #--------------------------------------------------------
      computeDirections = np.ones(3, dtype=np.bool)
      if np.all(unitVectors[0] == 0.0):
        computeDirections[0] = False
      if np.all(unitVectors[1] == 0.0):
        computeDirections[1] = False
      if np.all(unitVectors[2] == 0.0):
        computeDirections[2] = False

      POST.writeLagrangianFieldToVTK(referenceCoordinates, numAlongAxis, FTLE_T, FTLE_No_T,\
                                     Stretch_1, Stretch_2, Stretch_3, \
                                     Strain_1, Strain_2, Strain_3, lagFlags, computeDirections, \
                                     inputData.getLagrangianOutputFile(a_ID=outfileIndex))
      
      print("Finished Writing Lagrangian VTK data into %s" %inputData.getLagrangianOutputFile(a_ID=outfileIndex))
      outfileIndex = outfileIndex + inputData.m_OutIndexDelta
  
  end = time.time()
  elapsed = (end-start)/60.0
  print("Total runtime: %.2f minutes" %elapsed)
