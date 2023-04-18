#-------------------------------------------------------------------------------------------------
# This is a single file implementation of a script to compute tracerrangian fields from analytical or 
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
  print("Module Warning! Could not import module tracerInput, check configuration.")
  sys.exit()

try:
  import preProcess as PRE
except ImportError:
  print("Module Warning! Could not import module preProcess, check configuration.")
  sys.exit()

try:
  import postProcess as POST
except ImportError:
  print("Module Warning! Could not import module postProcess, check configuration.")
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
  print("Module Warning! Compiled Fortran 'main_fortran' module not available.")
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

  #-------------------------
  # Get advection parameters
  #-------------------------
  startTime      = inputData.getSimulationStartTime()
  endTime        = inputData.getSimulationStopTime()
  timeDelta      = inputData.getIntegrationTimeStep()
  dataDelta      = inputData.getDataTimeDelta()
  intScheme      = inputData.getIntegrationScheme()
  injectionStart = inputData.getInjectionStartTime()
  injectionDelta = inputData.getInjectionInterval()
  numInjections  = inputData.getNumberOfInjections()

  bcType      = inputData.m_FlowGeometryPrimitives['B0']
  flowModel   = inputData.m_FlowModel

  if intScheme == 'rk4':
    intScheme = 1
  elif intScheme == 'euler':
    intScheme = 2

  #----------------------------------------------
  # Get geometry information (assumes fixed mesh)
  #----------------------------------------------
  flowGeo     = inputData.getStandardFlowDomainGeometryDefinition()
  tracerGeo   = inputData.getStandardTracerDomainGeometryDefinition()
  numDim      = inputData.getProblemDimension()

  if inputData.getTracerInput() == inputData.m_RootPath + 'none':
    seedCart = True
  else:
    seedCart = False

  #---------------------------------
  # Seed a cartesian grid of tracers
  #---------------------------------
  if seedCart:
    #-------------------------------------------------------------------------
    # Check to see that tracer injection domain is constrained by flow domain.
    #-------------------------------------------------------------------------
    if ((tracerGeo['X0'] < flowGeo['X0']) or (tracerGeo['X1'] > flowGeo['X1']) or \
        (tracerGeo['Y0'] < flowGeo['Y0']) or (tracerGeo['Y1'] > flowGeo['Y1']) or \
        (tracerGeo['Z0'] < flowGeo['Z0']) or (tracerGeo['Z1'] > flowGeo['Z1'])):
      print("ERROR: Tracer domain must be constrained by flow domain. Check input parameters.")
      sys.exit()

    numX_tracer    = int((tracerGeo['X1'] - tracerGeo['X0'])/tracerGeo['DX']) + 1
    numY_tracer    = int((tracerGeo['Y1'] - tracerGeo['Y0'])/tracerGeo['DY']) + 1
    if inputData.getProblemDimension() < 3:
      numZ_tracer = 1
    else:
      numZ_tracer  = int((tracerGeo['Z1'] - tracerGeo['Z0'])/tracerGeo['DZ']) + 1

    tracerBBox   = np.zeros((6), dtype=np.float64, order='F')

    tracerBBox[0]   = tracerGeo['X0']
    tracerBBox[1]   = tracerGeo['X1']
    tracerBBox[2]   = tracerGeo['Y0']
    tracerBBox[3]   = tracerGeo['Y1']
    tracerBBox[4]   = tracerGeo['Z0']
    tracerBBox[5]   = tracerGeo['Z1']

    refDimns     = np.zeros((3), dtype=np.float64, order='F')
    refDimns[0]  = tracerGeo['DX']
    refDimns[1]  = tracerGeo['DY']
    refDimns[2]  = tracerGeo['DZ']

    numParticles = int(numX_tracer*numY_tracer*numZ_tracer)

  #---------------------------
  # Seed a VTK polydata object
  #---------------------------
  else:

    dataObj            = PRE.readVTK(inputData.getTracerInput(), a_Legacy='vtp')
    seedCoordinates    = PRE.getCoordinates(dataObj)
    numParticles       = seedCoordinates.shape[0]

  flowBBox     = np.zeros((6), dtype=np.float64, order='F')

  flowBBox[0]  = flowGeo['X0']
  flowBBox[1]  = flowGeo['X1']
  flowBBox[2]  = flowGeo['Y0']
  flowBBox[3]  = flowGeo['Y1']
  flowBBox[4]  = flowGeo['Z0']
  flowBBox[5]  = flowGeo['Z1']

  intFlags     = np.zeros((3), dtype=np.int32, order='F')
  intFlags[0]  = numDim
  intFlags[1]  = intScheme
  intFlags[2]  = bcType

  times        = np.zeros((5), dtype=np.float64, order='F')
  times[0]     = startTime
  times[1]     = endTime
  times[2]     = timeDelta
  times[3]     = injectionStart
  times[4]     = injectionDelta

  fileOut_IDTiming    = np.zeros((3), dtype=np.int32, order='F')
  fileOut_IDTiming[0] = inputData.m_OutIndexStart
  fileOut_IDTiming[1] = inputData.m_OutIndexStop
  fileOut_IDTiming[2] = inputData.m_OutIndexDelta

  fileOut_TTiming     = np.zeros((2), dtype=np.float64, order='F')
  fileOut_TTiming[0]  = inputData.m_OutTimeStart
  fileOut_TTiming[1]  = inputData.m_OutTimeInterval

  fileOut_Prefix      = inputData.getTracerOutputFile().split('.')[0]

  #---------------------------------
  # Simulated/experimental data sets
  #---------------------------------
  if (flowModel == 0):

    inPrefix = inputData.m_FlowDirectory + inputData.m_FlowFileTag

    fileIn_Prefix       = inPrefix + '.vel'

    fileIn_IDTiming     = np.zeros((3), dtype=np.int32, order='F')
    fileIn_IDTiming[0]  = inputData.m_DataIndexStart
    fileIn_IDTiming[1]  = inputData.m_DataIndexStop
    fileIn_IDTiming[2]  = inputData.m_DataIndexDelta

    fileIn_TTiming      = np.zeros((3), dtype=np.float64, order='F')
    fileIn_TTiming[0]   = inputData.m_DataTimeStart
    fileIn_TTiming[1]   = inputData.m_DataTimeStop
    fileIn_TTiming[2]   = inputData.m_DataTimeDelta

    #-----------------------------------------
    # Convert velocity files from .vtu to .bin
    #-----------------------------------------
    print("Generating velocity files...")

    fileIndex     = fileIn_IDTiming[0]
    timeStamp     = fileIn_TTiming[0]
    numDataFiles  = int((fileIn_IDTiming[1]-fileIn_IDTiming[0])/fileIn_IDTiming[2])+1
    timeIndex     = 0
    
    for f in range(numDataFiles):
      fileName    = inputData.getFlowDataFileName(a_ID=fileIndex)
      
      PRE.convertVTKDataToBin(fileName, inPrefix, timeIndex, timeStamp, \
               a_DataFieldName = inputData.getVelDataName(), \
               a_LegacyDataType = 'vtu')

      timeIndex   = timeIndex + 1
      fileIndex   = fileIndex + 1
      timeStamp   = timeStamp + fileIn_TTiming[2]
    
    print("Finished generating velocity files.\n")

    if inputData.isFixedMesh():
      print("Building Cell Locator Maps...")

      if seedCart:

        locatorObj  = PRE.createCellLocator(inputData.getFlowDataFileName(a_ID=inputData.m_DataIndexStart))
        
        print("Finished Building Cell Locator Maps \n")
        print("\nBuilding Tracer grids for cell-locator...")

        ref_dims  = np.asarray([tracerGeo['DX'], tracerGeo['DY'], tracerGeo['DZ']], dtype=np.float64, order='F')

        seedCoordinates = mesh.mesh.seed_points_cartesian(numParticles, ref_dims, tracerBBox)

        print("Finished Building Tracer grids for cell-locator...")
        print("Generating Cell-IDs for Tracers...")
        
        cellID_Init   = PRE.getCellIDs_Cartesian(seedCoordinates, locatorObj)

      else:

        locatorObj         = PRE.createCellLocator(inputData.getTracerInput(), a_Legacy='vtp')
        print("Finished Building Cell Locator Maps \n")
        print("Generating Cell-IDs for Tracers...")
        cellID_Init        = PRE.getCellIDs_VTK(dataObj, locatorObj)
      
      print("Finished Generating Cell-IDs for Tracers.\n")

    #-------------------------------------------------------
    # Get mesh data (if it is a fixed mesh, built only once)
    #-------------------------------------------------------
    if inputData.isFixedMesh():
      coordinates, connectivity, adjacency = PRE.getMeshData(inputData.getFlowDataFileName(a_ID=inputData.m_DataIndexStart))
      numNodes  = coordinates.shape[0]
      numCells  = connectivity.shape[0]

    #-----------------
    # Run main script.
    #-----------------
    print("Now integrating tracers from %.4fs to %.4fs" %(startTime, endTime))
    integration.integration.print_integration_parameters(flowModel, bcType, intScheme, endTime - startTime, timeDelta, flowBBox)
    
    if seedCart:
      MAIN.main.advect_seed_cart(intFlags, flowBBox, tracerBBox, numParticles, numInjections, \
                                 refDimns, times, fileIn_IDTiming, fileIn_TTiming, fileIn_Prefix, \
                                 fileOut_IDTiming, fileOut_TTiming, fileOut_Prefix, numNodes, \
                                 numCells, connectivity, adjacency, coordinates, cellID_Init)
    else:
      MAIN.main.advect_seed_vtk(intFlags, flowBBox, seedCoordinates, numParticles, numInjections, \
                                times, fileIn_IDTiming, fileIn_TTiming, fileIn_Prefix, \
                                fileOut_IDTiming, fileOut_TTiming, fileOut_Prefix, numNodes, \
                                numCells, connectivity, adjacency, coordinates, cellID_Init)
    print("\nFinished integrating tracers and writing data to files.")
      
  #-----------------------
  # Analytical flow fields
  #-----------------------
  else:

    if seedCart:
      MAIN.main.advect_seed_anlyt(flowModel, intFlags, flowBBox, tracerBBox, numParticles, \
                                  numInjections, refDimns, times, fileOut_IDTiming, fileOut_TTiming, \
                                  fileOut_Prefix)

    if not seedCart:
      print("ERROR: Analytical flow fields with vtk seeding not supported for advection yet.")
      sys.exit()

  end = time.time()
  elapsed = (end-start)/60.0
  print("Total runtime: %.2f minutes" %elapsed)
