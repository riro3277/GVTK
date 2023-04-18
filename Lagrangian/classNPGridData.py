import sys, os
import numpy as np
import vtk
#--------------------------------------------------------------
# FOR VTK IMPLEMENTATION ONLY
#----------------------------
# Class to encapsulate the data and methods for grid based data
#--------------------------------------------------------------
class GNPGridData(object):

    def __init__(self, a_NumDim=3, a_NumScalars=0, a_NumVectors=0):

        self.m_NumDim           = a_NumDim
        self.m_NumScalarData    = a_NumScalars
        self.m_NumVectorData    = a_NumScalars
        self.m_Coordinates  = None
        self.m_Connectivity = None
        self.m_Adjacency    = None

    @property
    def numNodes(self):
        return self.m_Coordinates.shape[0]

    @property
    def numCells(self):
        return self.m_Connectivity.shape[0]

    @property
    def numNodesPerCell(self):
        return self.m_Connectivity.shape[1]

    @property
    def coordinates(self):
        return self.m_Coordinates

    @coordinates.setter
    def coordinates(self, a_Coordinates):
        self.m_Coordinates = a_Coordinates

    @property
    def connectivity(self):
        return self.m_Connectivity

    @connectivity.setter
    def connectivity(self, a_Connectivity):
        self.m_Connectivity = a_Connectivity

    @property
    def adjacency(self):
        return self.m_Adjacency

    @adjacency.setter
    def adjacency(self, a_Adjacency):
        self.m_Adjacency = a_Adjacency

    def getGridDataFromFile(self, a_CoordinatesFile, a_ConnectivityFile, a_AdjacencyFile):

        if a_CoordinatesFile.endswith('.dat') or a_CoordinatesFile.endswith('.DAT'):
            self.m_Coordinates = np.loadtxt(a_CoordinatesFile)
        elif a_CoordinatesFile.endswith('.bin') or a_CoordinatesFile.endswith('.BIN'):
            sys.exit('GNPGridData: getGridDataFromFile binary file read not supported')

        if a_ConnectivityFile.endswith('.dat') or a_ConnectivityFile.endswith('.DAT'):
            self.m_Connectivity = np.loadtxt(a_ConnectivityFile)
        elif a_ConnectivityFile.endswith('.bin') or a_ConnectivityFile.endswith('.BIN'):
            sys.exit('GNPGridData: getGridDataFromFile binary file read not supported')

        if a_AdjacencyFile.endswith('.dat') or a_AdjacencyFile.endswith('.DAT'):
            self.m_Connectivity = np.loadtxt(a_AdjacencyFile)
        elif a_AdjacencyFile.endswith('.bin') or a_AdjacencyFile.endswith('.BIN'):
            sys.exit('GNPGridData: getGridDataFromFile binary file read not supported')

        print("Grid Initialized")

    def addVelocityToGrid(self, a_Array=None, a_File=None):
        print("Velocity Data Added Successfully")

    def addBoundaryMarkerToGrid(self, a_Array=None, a_File=None):
        print("Boundary Data Added Successfully")

    def gridLocate(self, a_X, a_CurrentCellID):
        print("Locating Particle in Cell")

    def interpolateGridData(self, a_X):
        print("Interpolating Grid Data")
