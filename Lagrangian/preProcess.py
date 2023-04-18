#-------------------------------------------------------------------------------------------------
# This is a collection of functions to perform pre-processing computations for particle transport
# simulations based on an input file.
#
# Additional detailed features to be included in future versions.
#
# Authors:      Zachariah Irwin, Debanjan Mukherjee
# Institution:  University of Colorado, Boulder
# Last Edit:    March 2020
#-------------------------------------------------------------------------------------------------
import sys, os

try:
  import numpy as np
except ImportError:
  print("Module Error! Numpy not installed.")
  sys.exit()

try:
  import vtk
except ImportError:
  print("Module Error! VTK module with Python bindings not installed.")
  sys.exit()

try:
  import postProcess as POST
except ImportError:
  print("Module Error! Could not import module postProcess, check configuration.")
  sys.exit()

#------------------------------------------------------------------------------------
## @brief A utility function that takes a VTK file as input, and extracts the data 
# from the file based on file type (e.g. polydata, unstructured etc.)
#
# @param[in] a_FileName  Name of the VTK mesh/data file
# @param[in] a_Legacy    For legacy VTK file types, the data type (vtu/vtp) is needed
#
#------------------------------------------------------------------------------------
def readVTK(a_FileName, a_Legacy='none'):

    if a_FileName.endswith('.vtu'):
        reader = vtk.vtkXMLUnstructuredGridReader()
    elif a_FileName.endswith('.vtp'):
        reader = vtk.vtkXMLPolyDataReader()
    elif a_FileName.endswith('.vtk'):
        if a_Legacy == 'none':
            print("Need To Specify Data Type For Legacy Files")
            sys.exit()
        elif a_Legacy == 'vtu':
            reader = vtk.vtkUnstructuredGridReader()
        elif a_Legacy == 'vtp':
            reader = vtk.vtkPolyDataReader()
    else:
        print("Unsupported File Extension")
        sys.exit()
        
    #print("Reading From File", a_FileName)

    reader.SetFileName(a_FileName)
    reader.Update()

    data = reader.GetOutput()

    return data

# --------------------------------------------------------------------------
## @brief A utility unction to create a VTK cell locator object from either
# an unstructured grid or polydata file
#
# @param[in]  a_FileName            Unstructured grid or polydata file name
# @param[in]  a_LocatorType (opt)   Type of vtk locator being employed
# @param[out] a_Legacy              For legacy VTK file types, the data type 
#                                   (vtu/vtp) is needed
#
# --------------------------------------------------------------------------
def createCellLocator(a_FileName, a_LocatorType=None, a_Legacy='vtp'):
    if a_FileName.endswith('vtu'):
        reader = vtk.vtkXMLUnstructuredGridReader()
    elif a_FileName.endswith('vtp'):
        reader = vtk.vtkXMLPolyDataReader()
    elif a_FileName.endswith('.vtk'):
        if a_Legacy == 'none':
            print("Need To Specify Data Type For Legacy Files")
            sys.exit()
        elif a_Legacy == 'vtu':
            reader = vtk.vtkUnstructuredGridReader()
        elif a_Legacy == 'vtp':
            reader = vtk.vtkPolyDataReader()
    else:
        print("Unsupported File Extension")
        sys.exit()

    reader.SetFileName(a_FileName)
    reader.Update()

    if a_LocatorType is None:
        locator = vtk.vtkCellTreeLocator()
    else:
        if a_LocatorType == 'oct':
            locator = vtk.vtkCellLocator()
        elif a_LocatorType == 'tre':
            locator = vtk.vtkCellTreeLocator()
        elif a_LocatorType == 'bsp':
            locator = vtk.vtkModifiedBSPTree()

    locator.SetDataSet(reader.GetOutput())
    locator.BuildLocator()

    return locator

#-------------------------------------------------------------------------------------------------
## @brief Function to take in a vtk mesh-data file and generate the mesh-data arrays in 
# Python/Fortran readable formats for mesh topology
#
# @param[in]   a_FileName               Name of the vtk file to be converted to binary
# @param[in]   a_InputLegacyDataType    If legacy vtk files are used, the data type vtu/vtp has to
#                                       be specified
# @param[in]   a_CreateAdjacency        Boolean flag for creating the adjacency/toplogy array          
# 
#--------------------------------------------------------------------------------------------------
def getMeshData(a_FileName, a_InputLegacyDataType = 'none', a_CreateAdjacency = True):

  #-----------------------------------------
  # Read the input file and extract the data
  #-----------------------------------------
  data    = readVTK(a_FileName)

  numNodes    = data.GetNumberOfPoints()
  numCells    = data.GetNumberOfCells()

  #------------------------------------------------------------------------
  # Create a dummy numpy array to just store the node number as an integer,
  # and an array of floats to store the coordinates of the nodes
  #------------------------------------------------------------------------
  numNodesArray    = numNodes*np.ones(1, dtype=np.int32)
  coordinatesArray = np.zeros((numNodes, 3), dtype=np.float64, order='F')

  #-------------------------------------------------
  # Populate the coordinates array from the vtk file
  #-------------------------------------------------
  print("Converting Coordinate Data")
  
  for p in range(numNodes):
    coordinatesArray[p,:] = data.GetPoint(p)[0], data.GetPoint(p)[1], data.GetPoint(p)[2]

  print("Finished Converting Coordinate Data\n")

  #----------------------------------------------------------------------------
  # Create a dummy numpy array to just store the number of cells as an integer, 
  # and an array of integers to store the cell-nodal connectivity
  # NOTE: ALL ELEMENTS ARE ASSUMED TO BE TETRAHEDRALIZED (4 NODES PER CELL)
  #----------------------------------------------------------------------------
  numCellsArray      = numCells*np.ones(1, dtype=np.int32)
  connectivityArray  = np.ones((numCells, 4), dtype=np.int32, order='F')

  #---------------------------------------------------------------------------
  # Populate the connectivity array from the vtk file, by iterating over cells
  #---------------------------------------------------------------------------
  print("Converting Connectivity")

  for c in range(numCells):
  
    cellNodeIDList = vtk.vtkIdList()

    cell = data.GetCell(c)

    if cell.GetNumberOfPoints() == 4:

      cellNodeList = cell.GetPointIds()
      
      for i in range(4):
        connectivityArray[c, i] = cellNodeList.GetId(i)# + a_Offset

    elif cell.GetNumberOfPoints() == 3:
  
      cellNodeList = cell.GetPointIds()
      for i in range(3):
        connectivityArray[c, i] = cellNodeList.GetId(i)# + a_Offset

      connectivityArray[c,3] = -1
    
    else:

      print("Non triangular and/or non tetrahedral elements not allowed")
      sys.exit()

  print("Finished Converting Connectivity\n")

  #----------------------------------------------------
  # Now we proceed to process the adjacency information
  #----------------------------------------------------
  adjacencyArray          = -1*np.ones((numCells, 4), dtype = np.int32, order='F')

  nodesInCell             = vtk.vtkIdList()
  nodesInFaceOpposite1    = vtk.vtkIdList()
  nodesInFaceOpposite2    = vtk.vtkIdList()
  nodesInFaceOpposite3    = vtk.vtkIdList()
  nodesInFaceOpposite4    = vtk.vtkIdList()

  nodesInEdgeOpposite1    = vtk.vtkIdList()
  nodesInEdgeOpposite2    = vtk.vtkIdList()
  nodesInEdgeOpposite3    = vtk.vtkIdList()
  nodesInEdgeOpposite4    = vtk.vtkIdList()

  #---------------------------------------
  # Each cell assumed to have 4 nodes each
  #---------------------------------------
  nodesInCell.SetNumberOfIds(4)

  #---------------------------------------------------------
  # Each face of a cell thereof assumed to have 3 nodes each 
  #---------------------------------------------------------
  nodesInFaceOpposite1.SetNumberOfIds(3)
  nodesInFaceOpposite2.SetNumberOfIds(3)
  nodesInFaceOpposite3.SetNumberOfIds(3)
  nodesInFaceOpposite4.SetNumberOfIds(3)

  #------------------------------------------
  # Each edge is assumed to have 2 nodes each
  #------------------------------------------
  nodesInEdgeOpposite1.SetNumberOfIds(2)
  nodesInEdgeOpposite2.SetNumberOfIds(2)
  nodesInEdgeOpposite3.SetNumberOfIds(2)
  nodesInEdgeOpposite4.SetNumberOfIds(2)

  print("Converting Adjacency")

  for cell in range(numCells):
      
    data.GetCellPoints(cell, nodesInCell)
  
    if data.GetCellType(cell) == vtk.VTK_TETRA:
      cellAdjacency = np.zeros(4, dtype=np.int32)
    elif data.GetCellType(cell) == vtk.VTK_TRIANGLE:
      cellAdjacency = np.zeros(4, dtype=np.int32)
      
    cellID1 = vtk.vtkIdList()
    cellID2 = vtk.vtkIdList()
    cellID3 = vtk.vtkIdList()
    cellID4 = vtk.vtkIdList()

    #--------------------------------------------------------------------------
    # NOTE: IN THE FOLLOWING LINES OF CODE, THE ORDERING OF NODES, AND FACES IS 
    # EXTREMELY CRUCIAL. THE FIRST THREE LINES ARE BASED ON OLDER CODE THAT
    # ASSUMES THAT FACE 1 IS THE FACE OPPOSITE NODE 1
    # THE SECOND THREE LINES ARE BASED ON THE MANUAL FROM SHAWN'S WEBSITE
    #--------------------------------------------------------------------------
    
    if data.GetCellType(cell) == vtk.VTK_TETRA:
        
      nodesInFaceOpposite1.SetId(0, nodesInCell.GetId(0))
      nodesInFaceOpposite1.SetId(1, nodesInCell.GetId(2))
      nodesInFaceOpposite1.SetId(2, nodesInCell.GetId(3))
      
      nodesInFaceOpposite2.SetId(0, nodesInCell.GetId(0))
      nodesInFaceOpposite2.SetId(1, nodesInCell.GetId(1))
      nodesInFaceOpposite2.SetId(2, nodesInCell.GetId(3))
      
      nodesInFaceOpposite3.SetId(0, nodesInCell.GetId(0))
      nodesInFaceOpposite3.SetId(1, nodesInCell.GetId(1))
      nodesInFaceOpposite3.SetId(2, nodesInCell.GetId(2))
      
      nodesInFaceOpposite4.SetId(0, nodesInCell.GetId(1))
      nodesInFaceOpposite4.SetId(1, nodesInCell.GetId(2))
      nodesInFaceOpposite4.SetId(2, nodesInCell.GetId(3))

    elif data.GetCellType(cell) == vtk.VTK_TRIANGLE:

      nodesInEdgeOpposite1.SetId(0, nodesInCell.GetId(0))
      nodesInEdgeOpposite1.SetId(1, nodesInCell.GetId(2))

      nodesInEdgeOpposite2.SetId(0, nodesInCell.GetId(0))
      nodesInEdgeOpposite2.SetId(1, nodesInCell.GetId(1))

      nodesInEdgeOpposite3.SetId(0, nodesInCell.GetId(1))
      nodesInEdgeOpposite3.SetId(1, nodesInCell.GetId(2))

    #--------------------------------------------------------------------------------
    # Find all cells using face (defined by node ID list) other than a specified cell
    # (this is a VTK Topological Inquiry function)
    #--------------------------------------------------------------------------------
    if data.GetCellType(cell) == vtk.VTK_TETRA:

      data.GetCellNeighbors(cell, nodesInFaceOpposite1, cellID1)
      data.GetCellNeighbors(cell, nodesInFaceOpposite2, cellID2)
      data.GetCellNeighbors(cell, nodesInFaceOpposite3, cellID3)
      data.GetCellNeighbors(cell, nodesInFaceOpposite4, cellID4)

    elif data.GetCellType(cell) == vtk.VTK_TRIANGLE:

      data.GetCellNeighbors(cell, nodesInEdgeOpposite1, cellID1)
      data.GetCellNeighbors(cell, nodesInEdgeOpposite2, cellID2)
      data.GetCellNeighbors(cell, nodesInEdgeOpposite3, cellID3)

    #---------------------------------------------------------------------------------
    # If a cell is internal, then the face should be topologically shared by two
    # cells, while if it is on boundary, some faces may be shared by only one 
    # using this idea, non-boundary neighbors and boundary neighbors are distinguished
    #---------------------------------------------------------------------------------
    if data.GetCellType(cell) == vtk.VTK_TETRA:

      if cellID1.GetNumberOfIds() == 1:
        cellAdjacency[0] = cellID1.GetId(0)
      elif cellID1.GetNumberOfIds() == 0:
        cellAdjacency[0] = -1

      if cellID2.GetNumberOfIds() == 1:
        cellAdjacency[1] = cellID2.GetId(0)
      elif cellID2.GetNumberOfIds() == 0:
        cellAdjacency[1] = -1

      if cellID3.GetNumberOfIds() == 1:
        cellAdjacency[2] = cellID3.GetId(0)
      elif cellID3.GetNumberOfIds() == 0:
        cellAdjacency[2] = -1

      if cellID4.GetNumberOfIds() == 1:
        cellAdjacency[3] = cellID4.GetId(0)
      elif cellID4.GetNumberOfIds() == 0:
        cellAdjacency[3] = -1

    elif data.GetCellType(cell) == vtk.VTK_TRIANGLE:

      if cellID1.GetNumberOfIds() == 1: 
        cellAdjacency[0] = cellID1.GetId(0)
      elif cellID1.GetNumberOfIds() == 0:
        cellAdjacency[0] = -1

      if cellID2.GetNumberOfIds() == 1:
        cellAdjacency[1] = cellID2.GetId(0)
      elif cellID2.GetNumberOfIds() == 0:
        cellAdjacency[1] = -1

      if cellID3.GetNumberOfIds() == 1:
        cellAdjacency[2] = cellID3.GetId(0)
      elif cellID3.GetNumberOfIds() == 0:
        cellAdjacency[2] = -1
    
    #-----------------------------------------------------------------------
    # Update this cell adjacency array into the overall mesh adjacency array
    #-----------------------------------------------------------------------
    adjacencyArray[cell, :] = cellAdjacency
  
  print("Finished Converting Adjacency\n")

  return coordinatesArray, connectivityArray, adjacencyArray

#------------------------------------------------------------------------------------------------------------
## @brief Function to convert the VTK velocity data files into binary format
#
# @note LT mesh velocity data format:
# - velocity data output name: <Prefix>.TSTAMP.bin
# - velocity data format: ts u0 v0 w0 u1 v1 w1 u(N-1) v(N-1) w(N-1)
# 
# @note The utility of this file should be such that the is called from within a loop
# over multiple velocity time-step data sets for a time-varying field
#
# @param[in] a_FileName              Name of the vtk file to be converted to binary
# @param[in] a_OutFilePrefixRoot     Prefix for the velocity output name as described in documentation
# @param[in] a_OutFileTimeIndex      Time index for time-stamp in the output file name in documentation
# @param[in] a_OutFileTimeStamp      Actual time-stamp value to be written into the output binary file
# @param[in] a_DataFiledName (opt)   The name of the data filed that is written as the interpolant field
# @param[in] a_LegacyDataType (opt)  If legacy vtk files are used, the data type vtu/vtp has to be specified
#
#------------------------------------------------------------------------------------------------------------
def convertVTKDataToBin(a_FileName, a_OutFilePrefixRoot, a_OutFileTimeIndex, a_OutFileTimeStamp, a_DataFieldName = 'velocity', a_LegacyDataType = 'none'):
    
  #-----------------------------------------
  # Read the input file and extract the data
  #-----------------------------------------
  data        = readVTK(a_FileName, a_Legacy=a_LegacyDataType)
  numNodes    = data.GetNumberOfPoints()
  dataValues  = data.GetPointData().GetArray(a_DataFieldName)
  
  velocityDataArrayToWrite  = np.zeros((numNodes,3), dtype=np.float64, order='F')
  timeStampArrayToWrite     = a_OutFileTimeStamp*np.ones(1)

  for p in range(numNodes):

    vel = np.asarray(dataValues.GetTuple(p))
    velocityDataArrayToWrite[p,:] = vel

  #---------------------------------------------------------------------------------
  # Open the file for velocity data in binary mode, and use the intrinsic 'tofile()'
  # method for numpy arrays to write the data into the file
  #---------------------------------------------------------------------------------
  dataFileName  = a_OutFilePrefixRoot +'.vel.' + str(a_OutFileTimeIndex) + '.bin'
  dataFileObj   = open(dataFileName, 'wb')
  timeStampArrayToWrite.tofile(dataFileObj)
  velocityDataArrayToWrite.tofile(dataFileObj)
  dataFileObj.close()

#-------------------------------------------------------------------------------------------------
## @brief Function to take in a vtk mesh-data file and generate the surface mesh-topology arrays
# in Python/Fortran readable formats. **2D MESHES ONLY**
#
# This function will extract the nodes and cells that form the boundaries of a 2D domain and
# return the associated connectivity and adjacency matrices.
#
# @param[in]   a_Filename               Name of the vtk file to be converted to binary
# @param[in]   a_InputLegacyDataType    If legacy vtk files are used, the data type vtu/vtp has to
#                                       be specified  
#
#--------------------------------------------------------------------------------------------------
def getSurfaceData(a_FileName, a_InputLegacyDataType = 'none'):

  #------------------
  # Get geometry data
  #------------------
  geoFilter = vtk.vtkGeometryFilter()
  geoFilter.SetInputData(a_Data)
  geoFilter.Update()

  #--------------------------
  # Extract domain boundaries
  #--------------------------
  featureEdges = vtk.vtkFeatureEdges()
  featureEdges.SetInputConnection(geo_filter.GetOutputPort())
  featureEdges.BoundaryEdgesOn()
  featureEdges.FeatureEdgesOff()
  featureEdges.Update()
  edgeData = featureEdges.GetOutput()

  numCells = edgeData.GetNumberOfCells()
  numNodes = edgeData.GetNumberOfPoints()

  #-------------------------------------------------------------------
  # Get connectivity data
  #
  # NOTE: Turning on ScalarConnectivity (not shown) will exclude holes 
  #       in the mesh
  #-------------------------------------------------------------------
  connectivity = vtk.vtkPolyDataConnectivityFilter()
  connectivity.SetInputConnection(featureEdges.GetOutputPort())
  connectivity.Update()
  connectivityData = connectivity.GetOutput()

  #----------------------------------------------
  # VTK List object for storing boundary node IDs
  #----------------------------------------------
  nodeList = vtk.vtkIdList()

  connectivityArray = np.zeros((numCells, 2), dtype=np.int32, order='F')

  for cellID in range(numCells):
    connectivityData.GetCellPoints(cellID, nodeList)
    connectivityArray[cellID,0] = nodeList.GetId(0)
    connectivityArray[cellID,1] = nodeList.GetId(1)

  coordinatesArray = np.zeros((numNodes, 3), dtype=np.float64, order='F')
  for nodeID in range(numNodes):
    coordinatesArray[p,:] = connectivityData.GetPoint(p)[0], connectivityData.GetPoint(p)[1], \
                            connectivityData.GetPoint(p)[2]

  return coordinatesArray, connectivityArray

#----------------------------------------------------------------------------------------------------------
## @brief Function to create an array of particle initial cell IDs to pass into Fortran
#
# Use instead of getCellIDs_VTK when seeding along a cartesian grid.
#
# @param[in]  a_Coordinates   numpy array of coordinates for every particle
# @param[in]  a_Locator       a vtk locator object
#
#----------------------------------------------------------------------------------------------------------
def getCellIDs_Cartesian(a_Coordinates, a_Locator):

  numParticles = a_Coordinates.shape[0]

  cellIDs = np.zeros((numParticles), dtype=np.int32, order='F')

  for p in range(numParticles):
    cellIDs[p]  = a_Locator.FindCell(a_Coordinates[p,:])

  return cellIDs

#--------------------------------------------------------------------------------------
## @brief Function to create an array of particle initial cell IDs to pass into Fortran
#
# Use instead of getCellIDs_Cartesian when seeding along a VTK points object.
#
# @param[in]  a_Data       a vtk reader object
# @param[in]  a_Locator    a vtk locator object
#
#--------------------------------------------------------------------------------------
def getCellIDs_VTK(a_Data, a_Locator):

  numPts = a_Data.GetNumberOfPoints()

  cellIDs = np.zeros((numPts), dtype=np.int32, order='F')

  for p in range(numPts):

    xyz         = a_Data.GetPoint(p)
    cellIDs[p]  = a_Locator.FindCell(xyz)

  return cellIDs

#---------------------------------------------------------------------------------
## @brief Function to create an array of particle coordinates to pass into Fortran
#
# Use for VTK point objects (e.g. seeding points along a clot boundary)
#
# @param[in]  a_Data       vtk reader object
# @param[out] coordinates  numpy array of particle coordinates
#
#----------------------------------------------------------------------------------
def getCoordinates(a_Data):

  numPts       = a_Data.GetNumberOfPoints()
  coordinates  = np.zeros((numPts,3), dtype=np.float64, order='F')

  for p in range(numPts):
    xyz               = a_Data.GetPoint(p)
    coordinates[p,:]  = xyz

  return coordinates
