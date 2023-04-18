import sys, os
import vtk
import numpy as np
import numpy.random as nprand

#---------------------------------------------------------
# data class encapsulating the methods for a geometric bin
# that helps get a binned wall contact detection possible
#---------------------------------------------------------
class BoundaryBin(object):

    def __init__(self, a_IX=None, a_IY=None, a_IZ=None):
        
        self.m_IX   = a_IX
        self.m_IY   = a_IY
        self.m_IZ   = a_IZ
        self.m_List = []

    def addToBin(self, a_ID):

        self.m_List.append(a_ID)

    def isBinEmpty(self):

        if not self.m_List:
            return True
        else:
            return False
    
#------------------------------------------------------
# container class to hold an iterable list of mesh-bins
#------------------------------------------------------
class BoundaryBinList(object):

    def __init__(self, a_Nx=10, a_Ny=10, a_Nz=None):
        
        #
        # initialize the number of cells
        #
        self.m_Nx   = a_Nx
        self.m_Ny   = a_Ny
        self.m_Nz   = a_Nz
        
        #
        # based upon initialization, determine number of dimensions
        #
        if a_Nz is None:
            self.m_ND   = 2
            self.m_Bins = [None]*(self.m_Nx*self.m_Ny)
        else:
            self.m_ND   = 3
            self.m_Bins = [None]*(self.m_Nx*self.m_Ny*self.m_Nz)

    def setBoundingBox(self, a_BoundaryMesh):
        
        #
        # get bounding box from mesh
        #
        bbox = a_BoundaryMesh.GetBounds()

        self.m_Xmin = bbox[0]
        self.m_Xmax = bbox[1]
        self.m_Ymin = bbox[2]
        self.m_Ymax = bbox[3]

        if self.m_ND == 2:

            self.m_Zmin = None
            self.m_Zmax = None

        elif self.m_ND == 3:

            self.m_Zmin = bbox[4]
            self.m_Zmax = bbox[5]
        
        #
        # assign the number of points (for creating structured grid visual)    
        #
        self.m_Px   = self.m_Nx + 1
        self.m_Py   = self.m_Ny + 1

        if self.m_ND == 2:
            self.m_Pz = None
        elif self.m_ND == 3:
            self.m_Pz = self.m_Nz + 1

    def createBins(self):
        
        self.m_Bins = []

        if self.m_ND == 2:
            for j in range(self.m_Ny):
                for i in range(self.m_Nx):
                    B = BoundaryBin(a_IX=i,a_IY=j)
                    self.m_Bins.append(B)

        elif self.m_ND == 3:
            for k in range(self.m_Nz):
                for j in range(self.m_Ny):
                    for i in range(self.m_Nx):
                        B = BoundaryBin(a_IX=i,a_IY=j,a_IZ=k)
                        self.m_Bins.append(B)

    def populateBinsFromBoundaryMesh(self, a_BoundaryMesh):
        
        cellPointID = vtk.vtkIdList()

        for c in range(a_BoundaryMesh.GetNumberOfCells()):

            if a_BoundaryMesh.GetCellType(c) == vtk.VTK_LINE:

                a_BoundaryMesh.GetCellPoints(c, cellPointID)
                x0  = a_BoundaryMesh.GetPoint(cellPointID.GetId(0))
                x1  = a_BoundaryMesh.GetPoint(cellPointID.GetId(1))
                ix  = int(self.m_Nx*(0.5*(x0[0] + x1[0]) - self.m_Xmin)/(self.m_Xmax - self.m_Xmin))
                iy  = int(self.m_Ny*(0.5*(x0[1] + x1[1]) - self.m_Ymin)/(self.m_Ymax - self.m_Ymin))
                self.m_Bins[self.linearIndex(ix,iy)].addToBin(c)

    def populateBinsFromSurfaceMesh(self, a_SurfaceMesh):

        cellPointID = vtk.vtkIdList()

        for c in range(a_SurfaceMesh.GetNumberOfCells()):

            if a_SurfaceMesh.GetCellType(c) == vtk.VTK_TRIANGLE:

                a_SurfaceMesh.GetCellPoints(c, cellPointID)
                x0  = a_SurfaceMesh.GetPoint(cellPointID.GetId(0))
                x1  = a_SurfaceMesh.GetPoint(cellPointID.GetId(1))
                x2  = a_SurfaceMesh.GetPoint(cellPointID.GetId(2))
                ix  = int(self.m_Nx*((1.0/3.0)*(x0[0]+x1[0]+x2[0]) - self.m_Xmin)/(self.m_Xmax - self.m_Xmin))
                iy  = int(self.m_Ny*((1.0/3.0)*(x0[1]+x1[1]+x2[1]) - self.m_Ymin)/(self.m_Ymax - self.m_Ymin))
                iz  = int(self.m_Nz*((1.0/3.0)*(x0[2]+x1[2]+x2[2]) - self.m_Zmin)/(self.m_Zmax - self.m_Zmin))
                self.m_Bins[self.linearIndex(ix,iy,a_K=iz)].addToBin(c)

    def visualizeBinnedBoundary3D(self, a_File):

        binGrid = vtk.vtkStructuredGrid()
        dims    = np.array([self.m_Px,self.m_Py, self.m_Pz])
        binGrid.SetDimensions(dims)

        points  = vtk.vtkPoints()
        points.SetNumberOfPoints(dims[0]*dims[1]*dims[2])

        marker  = vtk.vtkDoubleArray()
        marker.SetName('boundary')
        marker.SetNumberOfComponents(1)
        marker.SetNumberOfTuples(self.m_Nx*self.m_Ny*self.m_Nz)
        
        deltaZ  = (self.m_Zmax - self.m_Zmin)/float(dims[2] - 1)
        deltaY  = (self.m_Ymax - self.m_Ymin)/float(dims[1] - 1)
        deltaX  = (self.m_Xmax - self.m_Xmin)/float(dims[0] - 1)

        for k in range(dims[2]):
            kOffset = k*dims[0]*dims[1]
            for j in range(dims[1]):
                jOffset = j*dims[0]
                for i in range(dims[0]):
                    offset  = i + jOffset + kOffset
                    x       = self.m_Xmin + i*deltaX
                    y       = self.m_Ymin + j*deltaY
                    z       = self.m_Zmin + k*deltaZ
                    points.InsertPoint(offset, [x,y,z])
        
        for k in range(self.m_Nz):
            kOffset = k*self.m_Nx*self.m_Ny
            for j in range(self.m_Ny):
                jOffset = j*self.m_Nx
                for i in range(self.m_Nx):
                    offset  = i + jOffset + kOffset
                    check   = self.m_Bins[offset].isBinEmpty()
                    if check == True:
                        marker.InsertTuple1(offset, 0.0)
                    else:
                        marker.InsertTuple1(offset, 1.0)

        binGrid.SetPoints(points)
        binGrid.GetCellData().AddArray(marker)
        
        writer = vtk.vtkStructuredGridWriter()
        writer.SetFileName(a_File)
        writer.SetInputData(binGrid)
        writer.Update()
        writer.Write()

    def linearIndex(self, a_I, a_J, a_K=None):

        if a_K is None:
            index = a_I + self.m_Ny*a_J
        else:
            index = a_I + self.m_Nx*(a_J + self.m_Ny*a_K)

        return index


if __name__=="__main__":

    meshFile    = '/Users/debanjan-shaddenlab/Box Sync/Brown-Langevin/LVAD-Model-Test/inputs-LVAD-Inc90-Azi45/surface-geometry.vtp'
    outFile     = '/Users/debanjan-shaddenlab/Box Sync/Brown-Langevin/LVAD-Model-Test/boundaryBin.vtk'
    reader      = vtk.vtkXMLPolyDataReader()
    reader.SetFileName(meshFile)
    reader.Update()
    meshData    = reader.GetOutput()

    bList   = BoundaryBinList(a_Nx=40, a_Ny=40, a_Nz=100)
    bList.setBoundingBox(meshData)
    bList.createBins()
    bList.populateBinsFromSurfaceMesh(meshData)
    bList.visualizeBinnedBoundary3D(outFile)

