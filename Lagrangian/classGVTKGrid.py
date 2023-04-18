import sys, os
import numpy as np
import vtk

#--------------------------------------------------------------
# FOR VTK IMPLEMENTATION ONLY
#----------------------------
# Class to encapsulate the data and methods for grid based data
#--------------------------------------------------------------
class GVTKGridData(object):
    """A Generalized VTK based class for holding grid/mesh based data

    Attributes:
    -----------
    m_vtkData : vtk-data
        The grid or mesh data (vtk unstructured grid or polydata object)
    m_NumCells : int
        The number of cells/elements in the grid
    m_NumPoints : int
        The number of points/nodes in the grid

    """

    def __init__(self):
        """Initialize the GVTKGridData object

        Parameters
        ----------
        none

        Returns
        -------
        none

        """

        self.m_vtkData  = None
        self.m_NumCells = None
        self.m_NumPoints = None

    @property
    def vtkData(self):
        """Accessor function (getter) for member attribute m_vtkData

        """

        return self.m_vtkData

    @vtkData.setter
    def vtkData(self, a_VTKData):
        """Accessor function (setter) for member attribute m_vtkData

        Parameters
        ----------
        a_VTKData : vtk-data
            The vtk unstructured grid or polydata object as input

        Returns
        -------
        none

        """

        self.m_vtkData      = a_VTKData
        self.m_NumCells     = a_VTKData.GetNumberOfCells()
        self.m_NumPoints    = a_VTKData.GetNumberOfPoints()

    @property
    def numCells(self):
        return self.m_NumCells

    @property
    def numNodes(self):
        return self.m_NumPoints

    def getDataFromFile(self, a_FileName, a_FileType='vtu', a_FoamBlockID=0):
        """Function to extract the grid data from a file and store it into a
        GVTKGridData object

        Parameters
        ----------
        a_FileName : string
            Name of the file from which the grid data object will be read

        a_FileType : string, optional
            Type of the file/extension (vtu by default)

        a_FoamBlockID : int, optional
            Parameter specific to OpenFOAM files

        Returns
        -------
        none

        """

        if a_FileName.endswith('vtk') and a_FileType == 'vtu':

            reader = vtk.vtkUnstructuredGridReader()
            reader.SetFileName(a_FileName)
            reader.Update()

            self.m_vtkData = reader.GetOutput()

        elif a_FileName.endswith('vtk') and a_FileType == 'vtp':

            reader = vtk.vtkPolyDataReader()
            reader.SetFileName(a_FileName)
            reader.Update()

            self.m_vtkData = reader.GetOutput()

        elif a_FileName.endswith('vtu'):

            reader = vtk.vtkXMLUnstructuredGridReader()
            reader.SetFileName(a_FileName)
            reader.Update()

            self.m_vtkData = reader.GetOutput()

        elif a_FileName.endswith('vtp'):

            reader = vtk.vtkXMLPolyDataReader()
            reader.SetFileName(a_FileName)
            reader.Update()

            self.m_vtkData = reader.GetOutput()

        elif a_FileName.endswith('foam'):

            reader = vtk.vtkOpenFOAMReader()
            reader.SetFileName(a_FileName)
            reader.SkipZeroTimeOn()
            reader.Update()
            foamData = reader.GetOutput()
            numBlock = foamData.GetNumberOfBlocks()
            blockData = foamData.GetBlock(a_FoamBlockID)

            self.m_vtkData = blockData

    def extractDataArray(self, a_DataName):
        return self.m_vtkData.GetPointData().GetArray(a_DataName)

    def renameDataArray(self, a_OldName, a_NewName):
        if self.m_vtkData.GetPointData().HasArray(a_OldName) == 1:
            self.m_vtkData.GetPointData().GetArray(a_OldName).SetName(a_NewName)

    def addDataArray(self, a_DataName, a_DataArray):
        if a_DataArray.GetName() != a_DataName:
            a_DataArray.SetName(a_DataName)
        self.m_vtkData.GetPointData().AddArray(a_DataArray)

    def addDataArrayFromFile(self, a_FileName, a_AddName, a_ReadName, a_FileType='vtu', a_FoamBlockID=0):

        if a_FileName.endswith('vtk'):

            reader = vtk.vtkUnstructuredGridReader()
            reader.SetFileName(a_FileName)
            reader.Update()
            da = reader.GetOutput().GetPointData().GetArray(a_ReadName)

        elif a_FileName.endswith('vtu'):

            reader = vtk.vtkXMLUnstructuredGridReader()
            reader.SetFileName(a_FileName)
            reader.Update()
            da = reader.GetOutput().GetPointData().GetArray(a_ReadName)

        elif a_FileName.endswith('foam'):

            reader = vtk.vtkOpenFOAMReader()
            reader.SetFileName(a_FileName)
            reader.SkipZeroTimeOn()
            reader.Update()
            foamData = reader.GetOutput()
            numBlock = foamData.GetNumberOfBlocks()
            blockData = foamData.GetBlock(a_FoamBlockID)
            da = blockData.GetPointData().GetArray(a_ReadName)

        da.SetName(a_AddName)
        self.m_vtkData.GetPointData().AddArray(da)

    def removeDataArray(self, a_DataName):
        if self.m_vtkData.GetPointData().HasArray(a_DataName) == 1:
            self.m_vtkData.GetPointData().RemoveArray(a_DataName)

    def buildLocator(self, a_LocatorType=None):

        if a_LocatorType is None:
            self.m_Locator = vtk.vtkCellTreeLocator()
        else:
            if a_LocatorType == 'oct':
                self.m_Locator = vtk.vtkCellLocator()
            elif a_LocatorType == 'tre':
                self.m_Locator = vtk.vtkCellTreeLocator()
            elif a_LocatorType == 'bsp':
                self.m_Locator = vtk.vtkModifiedBSPTree()
            else:
                sys.exit('ERROR: Incompatible locator type specified to GVTKGridData')

        self.m_Locator.SetDataSet(self.m_vtkData)
        self.m_Locator.BuildLocator()

    def gridInterpolateAveraging(self, a_X, a_DataName, a_GetGradient=False, a_GetStatus=False, a_DataDim=3, a_DataType='double'):
        ## TODO: DEPRECATED IMPLEMENTATION

        cell = self.m_Locator.FindCell(a_X)
        cellNodeIds = vtk.vtkIdList()
        dataArray = self.m_vtkData.GetPointData().GetArray(a_DataName)

        if a_DataDim == 1:
            if a_DataType == 'double' or 'float':
                interpData = 0.0
                if a_GetGradient == True:
                    interpGrad = np.zeros(3,dtype=np.float32)
            elif a_DataType == 'int':
                interpData = 0
                if a_GetGradient == True:
                    interpGrad = np.zeros(3,dtype=np.int32)
        elif a_DataDim == 2:
            if a_DataType == 'double' or a_DataType == 'float':
                interpData = np.zeros(2, dtype=np.float32)
                if a_GetGradient == True:
                    interpGrad = np.zeros((2,3), dtype=np.float32)
            elif a_DataType == 'int':
                interpData = np.zeros(2, dtype=np.int32)
                if a_GetGradient == True:
                    interpGrad = np.zeros((2,3), dtype=np.int32)
        elif a_DataDim == 3:
            if a_DataType == 'double' or a_DataType == 'float':
                interpData = np.zeros(3, dtype=np.float32)
                if a_GetGradient == True:
                    interpGrad = np.zeros((3,3), dtype=np.float32)
            elif a_DataType == 'int':
                interpData = np.zeros(3, dtype=np.int32)
                if a_GetGradient == True:
                    interpGrad = np.zeros((3,3), dtypa=np.float32)

        if cell != -1:

            self.m_vtkData.GetCellPoints(cell, cellNodeIds)
            numNodes = cellNodeIds.GetNumberOfIds()

            if a_DataDim == 1:

                for n in range(numNodes):
                    interpData = interpData + dataArray.GetTuple1(cellNodeIds.GetId(n))

                if a_DataType == 'double' or a_DataType == 'float':
                    interpData = (1.0/float(numNodes))*interpData
                else:
                    interpData = int((1.0/float(numNodes))*interpData)

            elif a_DataDim == 2:

                for n in range(numNodes):
                    interpData[0] = interpData[0] + dataArray.GetTuple2(cellNodeIds.GetId(n))[0]
                    interpData[1] = interpData[1] + dataArray.GetTuple2(cellNodeIds.GetId(n))[1]

                if a_DataType == 'double' or a_DataType == 'float':
                    interpData = (1.0/float(numNodes))*interpData
                else:
                    interpData = int((1.0/float(numNodes))*interpData)

            elif a_DataDim == 3:

                for n in range(numNodes):
                    interpData[0] = interpData[0] + dataArray.GetTuple3(cellNodeIds.GetId(n))[0]
                    interpData[1] = interpData[1] + dataArray.GetTuple3(cellNodeIds.GetId(n))[1]
                    interpData[2] = interpData[2] + dataArray.GetTuple3(cellNodeIds.GetId(n))[2]

                if a_DataType == 'double' or a_DataType == 'float':
                    interpData = (1.0/float(numNodes))*interpData
                else:
                    interpData = int((1.0/float(numNodes))*interpData)

        if a_GetGradient == False and a_GetStatus == False:
            return interpData
        elif a_GetGradient == True and a_GetStatus == False:
            return [interpData, interpGrad]
        elif a_GetGradient == True and a_GetStatus == False:
            return [interpData, interpGrad, cell]
        elif a_GetGradient == False and a_GetStatus == True:
            return [interpData, cell]

    def gridInterpolateNodalBasis(self, a_X, a_DataName, a_GetGradient=False, a_GetStatus=False, a_DataDim=3, a_DataType='double'):

        import moduleMath as MM

        cell = self.m_Locator.FindCell(a_X)
        cellNodeIds = vtk.vtkIdList()
        dataArray = self.m_vtkData.GetPointData().GetArray(a_DataName)

        if a_DataDim == 1:
            if a_DataType == 'double' or 'float':
                interpData = 0.0
                if a_GetGradient == True:
                    interpGrad = np.zeros(3,dtype=np.float32)
            elif a_DataType == 'int':
                interpData = 0
                if a_GetGradient == True:
                    interpGrad = np.zeros(3,dtype=np.int32)
        elif a_DataDim == 2:
            if a_DataType == 'double' or a_DataType == 'float':
                interpData = np.zeros(2, dtype=np.float32)
                if a_GetGradient == True:
                    interpGrad = np.zeros((2,3), dtype=np.float32)
            elif a_DataType == 'int':
                interpData = np.zeros(2, dtype=np.int32)
                if a_GetGradient == True:
                    interpGrad = np.zeros((2,3), dtype=np.int32)
        elif a_DataDim == 3:
            if a_DataType == 'double' or a_DataType == 'float':
                interpData = np.zeros(3, dtype=np.float32)
                if a_GetGradient == True:
                    interpGrad = np.zeros((3,3), dtype=np.float32)
            elif a_DataType == 'int':
                interpData = np.zeros(3, dtype=np.int32)
                if a_GetGradient == True:
                    interpGrad = np.zeros((3,3), dtypa=np.float32)

        if cell != -1:

            self.m_vtkData.GetCellPoints(cell, cellNodeIds)
            numNodes = cellNodeIds.GetNumberOfIds()

            if numNodes == 4:

                x1 = self.m_vtkData.GetPoint(cellNodeIds.GetId(0))
                x2 = self.m_vtkData.GetPoint(cellNodeIds.GetId(1))
                x3 = self.m_vtkData.GetPoint(cellNodeIds.GetId(2))
                x4 = self.m_vtkData.GetPoint(cellNodeIds.GetId(3))

            else:

                sys.exit('ERROR: GVTKGridData interpolation currently only for tetrahedral cells')

            if a_DataDim == 1:

                if numNodes == 4:

                    # Developer notes AS:  made d1, d2, d3, d4 take nodal data as a list to avoid error 
                    #'object of type float has no len() for getTetrahedralCellLocalToGlobal function'
                    
                    d1 = [dataArray.GetTuple1(cellNodeIds.GetId(0))]
                    d2 = [dataArray.GetTuple1(cellNodeIds.GetId(1))]
                    d3 = [dataArray.GetTuple1(cellNodeIds.GetId(2))]
                    d4 = [dataArray.GetTuple1(cellNodeIds.GetId(3))]

                    # print('-'*50)
                    # print('d1 = %s'% d1)
                    # print('d2 = %s'% d2)
                    # print('d3 = %s'% d3)
                    # print('d4 = %s'% d4)
                    # sys.exit()


                    eta         = MM.getTetrahedralCellGlobalToLocal(x1,x2,x3,x4,a_X)
                    interpData  = MM.getTetrahedralCellLocalToGlobal(d1,d2,d3,d4,eta)

                    if a_GetGradient == True:
                        interpGrad  = MM.getTetrahedralScalarGradient(d1,d2,d3,d4,x1,x2,x3,x4,a_X)

                else:

                    sys.exit('ERROR: GVTKGridData interpolation currently only for tetrahedral cells')

            elif a_DataDim == 2:

                if numNodes == 4:

                    d1 = dataArray.GetTuple2(cellNodeIds.GetId(0))
                    d2 = dataArray.GetTuple2(cellNodeIds.GetId(1))
                    d3 = dataArray.GetTuple2(cellNodeIds.GetId(2))
                    d4 = dataArray.GetTuple2(cellNodeIds.GetId(3))

                    eta         = MM.getTetrahedralCellGlobalToLocal(x1,x2,x3,x4,a_X)
                    interpData  = MM.getTetrahedralCellLocalToGlobal(d1,d2,d3,d4,eta)

                    if a_GetGradient == True:
                        interpGrad  = MM.getTetrahedralVectorGradient(d1,d2,d3,d4,x1,x2,x3,x4,a_X)

                else:

                    sys.exit('ERROR: GVTKGridData interpolation currently only for tetrahedral cells')

            elif a_DataDim == 3:

                if numNodes == 4:

                    d1 = dataArray.GetTuple3(cellNodeIds.GetId(0))
                    d2 = dataArray.GetTuple3(cellNodeIds.GetId(1))
                    d3 = dataArray.GetTuple3(cellNodeIds.GetId(2))
                    d4 = dataArray.GetTuple3(cellNodeIds.GetId(3))

                    eta         = MM.getTetrahedralCellGlobalToLocal(x1,x2,x3,x4,a_X)
                    interpData  = MM.getTetrahedralCellLocalToGlobal(d1,d2,d3,d4,eta)

                    if a_GetGradient == True:
                        interpGrad  = MM.getTetrahedralVectorGradient(d1,d2,d3,d4,x1,x2,x3,x4,a_X)

                else:

                    sys.exit('ERROR: GVTKGridData interpolation currently only for tetrahedral cells')
                    
            elif a_DataDim == 9:

                if numNodes == 4:

                    d1 = dataArray.GetTuple9(cellNodeIds.GetId(0))
                    d2 = dataArray.GetTuple9(cellNodeIds.GetId(1))
                    d3 = dataArray.GetTuple9(cellNodeIds.GetId(2))
                    d4 = dataArray.GetTuple9(cellNodeIds.GetId(3))

                    eta         = MM.getTetrahedralCellGlobalToLocal(x1,x2,x3,x4,a_X)
                    interpData  = MM.getTetrahedralCellLocalToGlobal(d1,d2,d3,d4,eta)

                    if a_GetGradient == True:
                        interpGrad  = MM.getTetrahedralVectorGradient(d1,d2,d3,d4,x1,x2,x3,x4,a_X)

                else:

                    sys.exit('ERROR: GVTKGridData interpolation currently only for tetrahedral cells')

        if a_GetGradient == False and a_GetStatus == False:
            return interpData
        elif a_GetGradient == True and a_GetStatus == False:
            # return [interpData, interpGrad]
            return interpGrad    #----------Edits for Gradient calculation return statement is done by Sreeparna----------#
        elif a_GetGradient == True and a_GetStatus == True:
            return [interpData, interpGrad, cell]
        elif a_GetGradient == False and a_GetStatus == True:
            return [interpData, cell]
