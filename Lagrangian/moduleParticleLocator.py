import sys, os
import vtk
import numpy as np

#----------------------------------------------------------------
# Base class for all indexed cell walk locators to be implemented 
# as part of this library
# NOTE: Use this to replace the other class constructors
#----------------------------------------------------------------
class CellWalkLocator(object):
  
  def __init__(self, a_Mesh=None, a_File=None):
    
    if a_Mesh is None and a_File is None:
      sys.exit('Error in instantiating TriangleLocator')
    elif a_Mesh is not None:
      self.m_Mesh = a_Mesh
    elif a_File is not None:
      if a_File.endswith('vtu'):
        reader = vtk.vtkXMLUnstructuredGridReader()
      elif a_File.endswith('vtp'):
        reader = vtk.vtkXMLPolyDatareader()
        
      reader.SetFileName(a_File)
      reader.Update()
      data = reader.GetOutput()
      self.m_Mesh = data

#----------------------------------------------------------------------------------------
# Triangle cell-walking algorithm, used from paper:
# Lohner,R. and Ambrosiano,J.(1990)."A Vectorized Particle-Tracer For Unstructured Grids"
# NOTE: For now we are not re-instantiating the derived class from the base class, but
# this needs to be changed in the near future.
#----------------------------------------------------------------------------------------
class TriangleLocator(object):
  
  def __init__(self, a_Mesh=None, a_File=None):
    
    if a_Mesh is None and a_File is None:
      sys.exit('Error in instantiating TriangleLocator')
    elif a_Mesh is not None:
      self.m_Mesh = a_Mesh
    elif a_File is not None:
      if a_File.endswith('vtu'):
        reader = vtk.vtkXMLUnstructuredGridReader()
      elif a_File.endswith('vtp'):
        reader = vtk.vtkXMLPolyDatareader()
        
      reader.SetFileName(a_File)
      reader.Update()
      data = reader.GetOutput()
      self.m_Mesh = data
      
  def locateParticle(self, a_X, a_ID):
    
    nodesInCell = vtk.vtkIdList()
    nodesInCell.SetNumberOfIds(3)
    
    nodesOpposite_1 = vtk.vtkIdList()   # nodes opposite face id 1
    nodesOpposite_2 = vtk.vtkIdList()   # nodes opposite face id 2
    nodesOpposite_3 = vtk.vtkIdList()   # nodes opposite face id 3

    nodesOpposite_1.SetNumberOfIds(2)
    nodesOpposite_2.SetNumberOfIds(2)
    nodesOpposite_3.SetNumberOfIds(2)

    Nbr_1       = vtk.vtkIdList()
    Nbr_2       = vtk.vtkIdList()
    Nbr_3       = vtk.vtkIdList()

    nowID       = a_ID
    isCheck     = True

    while isCheck == True:
      
      #--------------------------
      # element nodal coordinates
      #--------------------------
      self.a_Mesh.GetCellPoints(nowID, nodesInCell)
      X1  = self.m_Mesh.GetPoint(nodesInCell.GetId(0))    # A
      X2  = self.m_Mesh.GetPoint(nodesInCell.GetId(1))    # B 
      X3  = self.m_Mesh.GetPoint(nodesInCell.GetId(2))    # C
      
      #-------------------------------------------
      # extract the nodes for each edge of the cell
      #-------------------------------------------
      nodesOpposite_1.SetId(nodesInCell.GetId(1))
      nodesOpposite_1.SetId(nodesInCell.GetId(2))
      
      nodesOpposite_2.SetId(nodesInCell.GetId(0))
      nodesOpposite_2.SetId(nodesInCell.GetId(2))
      
      nodesOpposite_3.SetId(nodesInCell.GetId(0))
      nodesOpposite_3.SetId(nodesInCell.GetId(1))
      
      #-----------------------------------------------
      # evaluate shape functions for particle location
      # notation: x_ij = Xi - Xj etc,
      #-----------------------------------------------
      x_21    = X2[0] - X1[0]     # x_BA
      y_31    = X3[1] - X1[1]     # y_CA
      y_21    = X2[1] - X1[1]     # y_BA
      x_31    = X3[0] - X1[0]     # x_CA
      x_p1    = a_X[0] - X1[0]    # x_PA
      y_p1    = a_X[1] - X1[1]    # y_PA
      zeta    = (y_31*x_p1 - x_31*y_p1)/(x_21*y_31 - y_21*x_31)
      eta     = (-y_21*x_p1 + x_21*y_p1)/(x_21*y_31 - y_21*x_31)
      N1      = 1.0 - zeta - eta
      N2      = zeta
      N3      = eta

      isInside    = (N1 >= 0.0) and (N1 <= 1.0) and (N2 >= 0.0) and (N2 <= 1.0) and (N3 >= 0.0) and (N3 <= 1.0)

      if isInside == True:
        
        isCheck = False
        locID   = nowID
      
      else:
        
        N_min   = np.min(N1, N2, N3)
        
        if N_min == N1:
          self.m_Mesh.GetCellNeighbors(nowID, nodesOpposite_1, nbrID)
        elif N_min == N2:
          self.m_Mesh.GetCellNeighbors(nowID, nodesOpposite_2, nbrID)
        elif N_min == N3:
          self.m_Mesh.GetCellNeighbors(nowID, nodesOpposite_3, nbrID)
          
        if nbrID == -1:
          isCheck = False
          locID   = -1
        else:
          isCheck = True
          nowID = nbrID
        
    return locID

#-----------------------------------------------------------------------------
# Tetrahedral cell walking algorithm implemented from paper:
# Kenwright, D.N. and Lane, D. (1995) "Optimization of Time-Dependent Particle 
# Tracing Using Tetrahedral Decomposition"
#-----------------------------------------------------------------------------
class TetrahedronLocator(object):

  def __init__(self, a_Mesh=None, a_File=None):
    
    if a_Mesh is None and a_File is None:
      sys.exit('Error in instantiating TriangleLocator')
    elif a_Mesh is not None:
      self.m_Mesh = a_Mesh
    elif a_File is not None:
      if a_File.endswith('vtu'):
        reader = vtk.vtkXMLUnstructuredGridReader()
      elif a_File.endswith('vtp'):
        reader = vtk.vtkXMLPolyDatareader()
        
      reader.SetFileName(a_File)
      reader.Update()
      data = reader.GetOutput()
      self.m_Mesh = data


#----------------------------------------------------------------------------------------
# To Be Implemented:
#-------------------
# a neighbor-tree based recursive implementation
# 1. start: E0 = root cell where particles are
# 2. list(E)        = list of neighbors of cell E
# 3. inside(E)      = function returning if E contains point P
# 4. checknbrs(E)   = recursive check of which'the neighbor in the tree has E
# 5. algorithm of checknbrs(E):
#   a.  if P is not inside(E):
#   b.      for E in list(E):
#   c.          check = checknbrs(E)
#   d.          if check is True: break
#----------------------------------------------------------------------------------------
class NeighborListLocator:

    def __init__(self, a_Mesh):
        self.m_Mesh = a_Mesh

    def nbrlist(self, a_ID):
        
        nodesInCell = vtk.vtkIdList()
        nodesInCell.SetNumberOfIds(3)

        nodesOpposite_1 = vtk.vtkIdList()
        nodesOpposite_2 = vtk.vtkIdList()
        nodesOpposite_3 = vtk.vtkIdList()

        nodesOpposite_1.SetNumberOfIds(2)
        nodesOpposite_2.SetNumberOfIds(2)
        nodesOpposite_3.SetNumberOfIds(2)

        Nbr_1       = vtk.vtkIdList()
        Nbr_2       = vtk.vtkIdList()
        Nbr_3       = vtk.vtkIdList()

        self.a_Mesh.GetCellPoints(a_ID, nodesInCell)

        #-------------------------------------------
        # extract the nodesfor each edge of the cell
        #-------------------------------------------
        nodesOpposite_1.SetId(nodesInCell.GetId(1))
        nodesOpposite_1.SetId(nodesInCell.GetId(2))

        nodesOpposite_2.SetId(nodesInCell.GetId(0))
        nodesOpposite_2.SetId(nodesInCell.GetId(2))

        nodesOpposite_3.SetId(nodesInCell.GetId(0))
        nodesOpposite_3.SetId(nodesInCell.GetId(1))

        #---------------------------
        # extract the cell neighbors
        #---------------------------
        self.a_Mesh.GetCellNeighbors(nowID, nodesOpposite_1, Nbr_1)
        self.a_Mesh.GetCellNeighbors(nowID, nodesOpposite_2, Nbr_2)
        self.a_Mesh.GetCellNeighbors(nowID, nodesOpposite_3, Nbr_3)

        

    def locateParticle(self, a_X, a_ID):

        #-------------------------------------
        # find the number of cells in the mesh
        #-------------------------------------
        numCells = self.a_Mesh.getNumberOfCells()

        #--------------------------
        # initialize all vtk arrays
        #--------------------------
        nodesInCell = vtk.vtkIdList()
        nodesInCell.SetNumberOfIds(3)

        nodesOpposite_1 = vtk.vtkIdList()
        nodesOpposite_2 = vtk.vtkIdList()
        nodesOpposite_3 = vtk.vtkIdList()

        nodesOpposite_1.SetNumberOfIds(2)
        nodesOpposite_2.SetNumberOfIds(2)
        nodesOpposite_3.SetNumberOfIds(2)

        Nbr_1       = vtk.vtkIdList()
        Nbr_2       = vtk.vtkIdList()
        Nbr_3       = vtk.vtkIdList()

        NbrNodes_1  = vtk.vtkIdList()
        NbrNodes_2  = vtk.vtkIdList()
        NbrNodes_3  = vtk.vtkIdList()

        nowID       = a_ID
        isLocated   = False

        while isLocated == False:

            #-------------------------------------------
            # find the nodes in the current locator cell
            #-------------------------------------------
            self.a_Mesh.GetCellPoints(nowID, nodesInCell)

            #-------------------------------------------
            # extract the nodesfor each edge of the cell
            #-------------------------------------------
            nodesOpposite_1.SetId(nodesInCell.GetId(1))
            nodesOpposite_1.SetId(nodesInCell.GetId(2))

            nodesOpposite_2.SetId(nodesInCell.GetId(0))
            nodesOpposite_2.SetId(nodesInCell.GetId(2))

            nodesOpposite_3.SetId(nodesInCell.GetId(0))
            nodesOpposite_3.SetId(nodesInCell.GetId(1))

            #---------------------------
            # extract the cell neighbors
            #---------------------------
            self.a_Mesh.GetCellNeighbors(nowID, nodesOpposite_1, Nbr_1)
            self.a_Mesh.GetCellNeighbors(nowID, nodesOpposite_2, Nbr_2)
            self.a_Mesh.GetCellNeighbors(nowID, nodesOpposite_3, Nbr_3)

            if Nbr_1.GetId(0) != -1: 
                self.a_Mesh.GetCellPoints(Nbr_1, NbrNodes_1)
                check = pointInTriangle(self.a_Mesh.GetPoint(NbrNodes_1.GetId(0)),
                        self.a_Mesh.GetPoint(NbrNodes_1.GetId(1)),
                        self.a_Mesh.GetPoint(NbrNodes_1.GetId(2)))
                
            if Nbr_2.GetId(0) != -1:
                self.a_Mesh.GetCellPoints(Nbr_2, NbrNodes_2)
                check = pointInTriangle(self.a_Mesh.GetPoint(NbrNodes_2.GetId(0)),
                        self.a_Mesh.GetPoint(NbrNodes_2.GetId(1)),
                        self.a_Mesh.GetPoint(NbrNodes_2.GetId(2))) 

            if Nbr_3.GetId(0) != -1:
                self.a_mesh.GetCellPoints(Nbr_3, NbrNodes_3)
                check = pointInTriangle(self.a_Mesh.GetPoint(NbrNodes_3.GetId(0)),
                        self.a_Mesh.GetPoint(NbrNodes_3.GetId(1)),
                        self.a_Mesh.GetPoint(NbrNodes_3.GetId(2)))

            





        



