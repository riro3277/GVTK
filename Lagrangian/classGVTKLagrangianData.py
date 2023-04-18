import sys, os
import numpy as np
import vtk

from moduleFileIO import initializeTracerFromFile

#----------------------------------------------------------------------
# FOR VTK IMPLEMENTATION ONLY
#----------------------------
# Class to encapsulate the data and methods for particle based data
# this could be passive tracers, active tracers, and inertial particles
#----------------------------------------------------------------------
class GVTKLagrangianData(object):
    """A Generalized VTK based class for holding Lagrangian tracer based data

    This class is meant for a VTK based implementation only. The Lagrangian data
    could include passive tracers, active agents, and inertial particles. The VTK
    based implementations are meant for quick visualization and prototyping or
    proof of concept computations. Higher performance computations will require
    using the NumPy based implementations (TO BE IMPLEMENTED).

    Attributes
    ----------
    m_Points : vtk-points
        The set of coordinates for the Lagrangian tracers/particles
    m_NumParticles : int
        The current number of Lagrangian tracers/particles
    m_DataFields : list
        The set of data attrributes for each Lagrangian tracers/particles
        (this is a list of vtkDataArray objects)
    m_NumData : int
        The current number of data attributes for each Lagrangian tracers/particles

    """

    def __init__(self, a_File=None, a_VTKData=None, a_DataFields=None):
        """Initialize the GVTKLagrangianData object

        While all three initialization arguments are optional, it is mandatory to
        specify at least one argument to correctly initialize the object. The most
        convenient, and recommended, option is to initialize from a file.

        Parameters
        ----------
        a_FileName : str, optional
            Name of the vtk file (with path) from which tracers are initialized
            (this filename has to be for a vtk polydata file)
        a_VTKData : vtk-polydata, optional
            The vtk polydata object from which tracers are initialized
        a_DataFields : list, optional
            List of data fields to be added at initialization time

        Returns
        -------
        none

        """

        if a_File == None and a_VTKData != None:

            if a_VTKData.IsA('vtkPolyData') != 1:
                sys.exit('ERROR: GVTKLagrangianData not initialized as polydata')
            else:
                self.m_Points = a_VTKData.GetPoints()
                self.m_NumParticles = self.m_Points.GetNumberOfPoints()

                if a_DataFields != None:
                    self.m_DataFields = [a_VTKData.GetPointData().GetArray(f) for f in a_DataFields]
                    self.m_NumData = len(self.m_DataFields)
                else:
                    self.m_DataFields = []
                    self.m_NumData = 0

        elif a_File != None and a_VTKData == None:

            data = initializeTracerFromFile(a_File)

            self.m_Points = data.GetPoints()
            self.m_NumParticles = self.m_Points.GetNumberOfPoints()

            if a_DataFields != None:
                self.m_DataFields = [data.GetPointData().GetArray(f) for f in a_DataFields]
                self.m_NumData = len(self.m_DataFields)
            else:
                self.m_DataFields = []
                self.m_NumData = 0

        else:

            sys.exit('ERROR: GVTKLagrangianData needs file or vtkData object')

    @property
    def points(self):
        """Property function to get the underlying particle locations in form of
        a VTKPoints object

        """

        return self.m_Points

    @property
    def numParticles(self):
        """Property function to get the number of particles

        """

        return self.m_NumParticles

    @property
    def numData(self):
        """Property function to get the number of data arrays associated with
        the Lagrangian data

        """

        return self.m_NumData

    def getX(self, a_ID):
        """Retrieve the position coordinates of a specified particle/tracer

        Parameters
        ----------
        a_ID : int
            ID of the Lagrangian entity whose coordinates are sought

        Returns
        -------
        1X3 list : coordinates of the Lagrangian entity

        """

        return self.m_Points.GetPoint(a_ID)

    def setX(self, a_ID, a_X):
        """Set the position coordinates of a specified particle/tracer

        Parameters
        ----------
        a_ID : int
            ID of the Lagrangian entity whose coordinates are to be set
        a_X : 1X3 float
            Position of the Lagrangian entity

        Returns
        -------
        none

        """

        self.m_Points.SetPoint(a_ID, a_X[0], a_X[1], a_X[2])

    def addScalarData(self, a_ArrayName, a_InitValues, a_DType='double'):
        """Add a new scalar data array/field to the GVTKLagrangianData object

        The scalar data can be set in two different ways by settint the argument
        a_InitValues
        - a_InitValues set as a single scalar (int/float/double) will set the
        same fixed scalar for all entries of the array
        - a_InitValues set as an array of length m_NumParticles will set each
        entry of the scalar to its corresponding entry in the input array

        Parameters
        ----------
        a_ArrayName : str
            Name of the scalar array field that will be added to the object
        a_InitValues : mixed
            The actual array values to be added (see documentation for the two
            different choices for setting this input argument)
        a_DType : str, optional
            A string specifying the type of the data 'int'/'float'/'double'

        Returns
        -------
        none

        """

        for f in self.m_DataFields:
            if f.GetName() == a_ArrayName:
                self.m_DataFields.remove(f)

        if a_DType == 'double' or a_DType == 'float':
            da = vtk.vtkDoubleArray()
        elif a_DType == 'int':
            da = vtk.vtkIntArray()
        else:
            sys.exit('ERROR: GVTKLagrangianData addScalarData a_DType incorrectly specified')

        da.SetNumberOfTuples(self.m_NumParticles)
        da.SetNumberOfComponents(1)
        da.SetName(a_ArrayName)

        for p in range(self.m_NumParticles):
            if not hasattr(a_InitValues, '__len__'):
                da.SetTuple1(p, a_InitValues)
            elif len(a_InitValues) == self.m_NumParticles:
                da.SetTuple1(p, a_InitValues[p])
            else:
                sys.exit('ERROR: GVTKLagrangianData addScalarData needs either a scalar or array of size numParticles')

        self.m_DataFields.append(da)

    def setScalarData(self, a_SetValues, a_ArrayName=None, a_ArrayID=None, a_DataID=None, a_DType='double'):
        """Set the values of scalar data array/field in the GVTKLagrangianData object

        The scalar data can be set in two different ways by settint the argument
        a_SetValues
        - a_SetValues set as a single scalar (int/float/double) will set the
        same fixed scalar for all entries of the array
        - a_SetValues set as an array of length m_NumParticles will set each
        entry of the scalar to its corresponding entry in the input array

        Note
        ----
        Only a_ArrayName or a_ArrayID has to be set, not both.

        Parameters
        ----------
        a_InitValues : mixed
            The actual array values to be added (see documentation for the two
            different choices for setting this input argument)
        a_ArrayName : str, optional
            Name of the scalar array field that will be added to the object
        a_ArrayID : int, optional
            Position of the array in the m_DataFields list for the object
        a_DType : str, optional
            A string specifying the type of the data 'int'/'float'/'double'

        Returns
        -------
        none

        """

        if a_ArrayName == None and a_ArrayID == None:
            sys.exit('ERROR: GVTKLagrangianData setScalarData needs either an ArrayName or an ArrayID input')
        elif a_ArrayName != None and a_ArrayID != None:
            sys.exit('ERROR: GVTKLagrangianData setScalarData needs either an ArrayName or an ArrayID input')
        elif a_ArrayName != None and a_ArrayID == None:
            for i in range(len(self.m_DataFields)):
                tempField = self.m_DataFields[i]
                if tempField.GetName() == a_ArrayName:
                    da = tempField
        elif a_ArrayName == None and a_ArrayID != None:
            da = self.m_DataFields[a_ArrayID]

        if a_DataID == None:
            for p in range(self.m_NumParticles):
                if hasattr(a_SetValues, '__len__'):
                    if len(a_SetValues) != self.m_NumParticles:
                        sys.exit('When passing an array of scalar values, length of array must match number of particles')
                    da.SetTuple1(p, a_SetValues[p])
                else:
                    da.SetTuple1(p, a_SetValues)
        else:
            da.SetTuple1(a_DataID, a_SetValues)

    def getScalarData(self, a_ArrayName=None, a_ArrayID=None, a_DataID=None, a_DType='double'):
        """Get the values of scalar data array/field of the GVTKLagrangianData object

        Note
        ----
        - Only a_ArrayName or a_ArrayID has to be set, not both.
        - If a_DataID is set, then only that particle's data will be retrieved.
        - If a_DataID is not set, then the entire array will be retrived
        as a vtkArray object

        Parameters
        ----------
        a_ArrayName : str, optional
            Name of the scalar array field that will be added to the object
        a_ArrayID : int, optional
            Position of the array in the m_DataFields list for the object
        a_DataID : int, optional
            ID of the particle/tracer for which the data will be retrieved
        a_DType : str, optional
            A string specifying the type of the data 'int'/'float'/'double'

        Returns
        -------
        int/float/double scalar : if a_DataID is an int
        vtk-Array : if a_DataID is None

        """

        if a_ArrayName == None and a_ArrayID == None:
            sys.exit('ERROR: GVTKLagrangianData getScalarData needs either an ArrayName or an ArrayID input')
        elif a_ArrayName != None and a_ArrayID != None:
            sys.exit('ERROR: GVTKLagrangianData getScalarData needs either an ArrayName or an ArrayID input')
        elif a_ArrayName != None and a_ArrayID == None:
            for i in range(len(self.m_DataFields)):
                tempField = self.m_DataFields[i]
                if tempField.GetName() == a_ArrayName:
                    da = tempField
        elif a_ArrayName == None and a_ArrayID != None:
            da = self.m_DataFields[a_ArrayID]

        if a_DataID == None:
            return da
        else:
            return da.GetTuple1(a_DataID)

    def addVectorData(self, a_ArrayName, a_InitValues, a_DType='double'):
        """Add a new vector data array/field to the GVTKLagrangianData object

        The vector data can be set in two different ways by setting the argument
        a_InitValues
        - a_InitValues set as a single row vector (int/float/double) will set the
        same fixed scalar for all entries of the array
        - a_InitValues set as a 2D array of length m_NumParticles will set each
        entry of the scalar to its corresponding row entry in the input array

        Parameters
        ----------
        a_ArrayName : str
            Name of the scalar array field that will be added to the object
        a_InitValues : mixed
            The actual array values to be added (see documentation for the two
            different choices for setting this input argument)
        a_DType : str, optional
            A string specifying the type of the data 'int'/'float'/'double'

        Returns
        -------
        none

        """
        for f in self.m_DataFields:
            if f.GetName() == a_ArrayName:
                self.m_DataFields.remove(f)

        numTuples       = a_InitValues.shape[0]
        numComponents   = a_InitValues.shape[1]

        if a_DType == 'double' or a_DType == 'float':
            da = vtk.vtkDoubleArray()
        elif a_DType == 'int':
            da = vtk.vtkIntArray()
        else:
            sys.exit('ERROR: GVTKLagrangianData addVectorData a_DType incorrectly specified')

        da.SetNumberOfTuples(numTuples*3)
        da.SetNumberOfComponents(numComponents)
        da.SetName(a_ArrayName)

        for p in range(self.m_NumParticles):
            if numTuples == 1:
                if numComponents == 2:
                    da.SetTuple2(p, a_InitValues[0], a_InitValues[1])
                elif numComponents == 3:
                    da.SetTuple3(p, a_InitValues[0], a_InitValues[1], a_InitValues[2])
            elif numTuples == self.m_NumParticles:
                if numComponents == 2:
                    da.SetTuple2(p, a_InitValues[p,0], a_InitValues[p,1])
                elif numComponents == 3:
                    da.SetTuple3(p, a_InitValues[p,0], a_InitValues[p,1], a_InitValues[p,2])
            else:
                sys.exit('ERROR: GVTKLagrangianData addVectorData needs either a single vector or array of numParticles rows')

        self.m_DataFields.append(da)

    def setVectorData(self, a_SetValues, a_ArrayName=None, a_ArrayID=None, a_DataID=None, a_DType='double'):
        """Set the values of vector data array/field in the GVTKLagrangianData object

        The vector data can be set in two different ways by setting the argument
        a_SetValues
        - a_SetValues set as a single row vector (int/float/double) will set the
        same fixed scalar for all entries of the array
        - a_SetValues set as a 2D array of length m_NumParticles will set each
        entry of the scalar to its corresponding entry in the input array

        Note
        ----
        Only a_ArrayName or a_ArrayID has to be set, not both.

        Parameters
        ----------
        a_InitValues : mixed
            The actual array values to be added (see documentation for the two
            different choices for setting this input argument)
        a_ArrayName : str, optional
            Name of the scalar array field that will be added to the object
        a_ArrayID : int, optional
            Position of the array in the m_DataFields list for the object
        a_DType : str, optional
            A string specifying the type of the data 'int'/'float'/'double'

        Returns
        -------
        none

        """
        if a_DataID == None:
            numTuples     = a_SetValues.shape[0]
            numComponents = a_SetValues.shape[1]
        else:
            numTuples     = 1
            numComponents = a_SetValues.shape[0]

        if a_ArrayName == None and a_ArrayID == None:
            sys.exit('ERROR: GVTKLagrangianData setVectorData needs either an ArrayName or an ArrayID input')
        elif a_ArrayName != None and a_ArrayID != None:
            sys.exit('ERROR: GVTKLagrangianData setVectorData needs either an ArrayName or an ArrayID input')
        elif a_ArrayName != None and a_ArrayID == None:
            for i in range(len(self.m_DataFields)):
                tempField = self.m_DataFields[i]
                if tempField.GetName() == a_ArrayName:
                    da = tempField
        elif a_ArrayName == None and a_ArrayID != None:
            da = self.m_DataFields[a_ArrayID]

        if numTuples == 1:
            if numComponents == 2:
                da.SetTuple2(a_DataID, a_SetValues[0], a_SetValues[1])
            elif numComponents == 3:
                da.SetTuple3(a_DataID, a_SetValues[0], a_SetValues[1], a_SetValues[2])
        elif numTuples == self.m_NumParticles:
            for p in range(self.m_NumParticles):
                if numComponents == 2:
                    da.SetTuple2(p, a_SetValues[p,0], a_SetValues[p,1])
                elif numComponents == 3:
                    da.SetTuple3(p, a_SetValues[p,0], a_SetValues[p,1], a_SetValues[p,2])
        else:
            sys.exit('ERROR: GVTKLagrangianData setVectorData needs either a single vector or array of numParticles rows')

    def getVectorData(self, a_ArrayName=None, a_ArrayID=None, a_DataID=None, a_DType='doube'):
        """Get the values of scalar data array/field of the GVTKLagrangianData object

        Note
        ----
        - Only a_ArrayName or a_ArrayID has to be set, not both.
        - If a_DataID is set, then only that particle's data will be retrieved.
        - If a_DataID is not set, then the entire array will be retrived
        as a vtkArray object

        Parameters
        ----------
        a_ArrayName : str, optional
            Name of the scalar array field that will be added to the object
        a_ArrayID : int, optional
            Position of the array in the m_DataFields list for the object
        a_DataID : int, optional
            ID of the particle/tracer for which the data will be retrieved
        a_DType : str, optional
            A string specifying the type of the data 'int'/'float'/'double'

        Returns
        -------
        int/float/double vector : if a_DataID is an int
        vtk-Array : if a_DataID is None

        """

        if a_ArrayName == None and a_ArrayID == None:
            sys.exit('ERROR: GVTKLagrangianData getVectorData needs either an ArrayName or an ArrayID input')
        elif a_ArrayName != None and a_ArrayID != None:
            sys.exit('ERROR: GVTKLagrangianData getVectorData needs either an ArrayName or an ArrayID input')
        elif a_ArrayName != None and a_ArrayID == None:
            for i in range(len(self.m_DataFields)):
                tempField = self.m_DataFields[i]
                if tempField.GetName() == a_ArrayName:
                    da = tempField
        elif a_ArrayName == None and a_ArrayID != None:
            da = self.m_DataFields[a_ArrayID]

        numTuples       = da.GetNumberOfTuples()
        numComponents   = da.GetNumberOfComponents()

        if a_DataID == None:
            return da
        else:
            if numComponents == 2:
                return da.GetTuple2(a_DataID)
            elif numComponents == 3:
                return da.GetTuple3(a_DataID)

    def injectParticles(self, a_InputParticles):
        """Inject a specified number of particles with specified data fields
        into the existing GVTKLagrangianData object

        Note
        ----
        The injected data a_InputParticles can be in two forms
        - a vtkPoints object : in which case all compatible data fields are
        initialized to zero valued arrays
        - a vtkPolyData object : in which case data fields are imported from
        the polydata object directly

        Parameters
        ----------
        a_InputParticles : vtk-points or vtk-polydata
            The particle data to be injected into the GVTKLagrangianData object

        Returns
        -------
        none

        """

        #------------------------------------------------------------------------
        # check whether the input arguments are both of the same type,
        # and that they both have the same number of data arrays - exit otherwise
        #------------------------------------------------------------------------
        if a_InputParticles.IsA('vtkPoints') == 1:

            #------------------------------------------
            # get the number of points in each data set
            #------------------------------------------
            n_In    = self.m_NumParticles
            n_Add   = a_InputParticles.GetNumberOfPoints()

            #------------------------------------------------------------
            # add points from the injection dataset into original dataset
            #------------------------------------------------------------
            for p in range(n_Add):

                xyz = a_InputParticles.GetPoint(p)
                self.m_Points.InsertPoint(p + n_In, xyz)

            #-----------------------------------------------------------
            # finally, update the internal number of particles
            # note that for a vtkPoints initialization, the number of
            # data fields remain unchanged as only coordinates are added
            #-----------------------------------------------------------
            self.m_NumParticles = self.m_Points.GetNumberOfPoints()

        elif a_InputParticles.IsA('vtkPolyData') == 1:

            if ( self.m_NumData != a_InputParticles.GetPointData().GetNumberOfArrays() ):
                sys.exit("ERROR: GVTKLagrangianData.injectParticles - both arguments need same number of data")

            #------------------------------------------
            # get the number of points in each data set
            #------------------------------------------
            n_In    = self.m_NumParticles
            n_Add   = a_InputParticles.GetNumberOfPoints()

            #------------------------------------------------------------
            # add points from the injection dataset into original dataset
            #------------------------------------------------------------
            for p in range(n_Add):

                xyz = a_InputParticles.GetPoint(p)
                self.m_Points.InsertPoint(p + n_In, xyz)

            #-------------------------------------------------------------
            # now add tuples from each array in the injection dataset into
            # the corresponding array in the original dataset
            #-------------------------------------------------------------
            for a in range(self.m_NumData):

                arrName = self.m_DataFields[a].GetName()

                for p in range(n_Add):

                    daTuple = a_InputParticles.GetPointData().GetArray(arrName).GetTuple(p)
                    self.m_DataFields[a].InsertTuple(p + n_In, daTuple)

            #-------------------------------------------------------
            # finally update the number of particles as well as
            # the number of data fields at the end of the injection
            #-------------------------------------------------------
            self.m_NumParticles = self.m_Points.GetNumberOfPoints()
            self.m_NumData      = len(self.m_DataFields)

        else:

            sys.exit("ERROR: GVTKLagrangianData.injectParticles - input should be points or polyData")

    def writeData(self, a_FileName):
        """Write out the particle information stored in GVTKLagrangianData
        object into a specified vtkPolyData file

        Parameters
        ----------
        a_FileName : str
            Name of the file (with path) where the particle data is to be dumped

        Returns
        -------
        none

        """

        if a_FileName.endswith('vtp'):
            writer = vtk.vtkXMLPolyDataWriter()
        else:
            writer = vtk.vtkPolyDataWriter()

        writer.SetFileName(a_FileName)
        writer.SetFileTypeToBinary()

        wd = vtk.vtkPolyData()
        wd.SetPoints(self.m_Points)
        for a in self.m_DataFields:
            wd.GetPointData().AddArray(a)

        if vtk.VTK_MAJOR_VERSION <= 5.0:
            writer.SetInput(wd)
        else:
            writer.SetInputData(wd)

        writer.Update()
        writer.Write()
