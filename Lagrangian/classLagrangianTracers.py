#---------------------------------------------------------------------------------------------
# This is a module that encapsulates the physics and numerics of a light tracer integration
# calculation through a specified flow field
# 
# Author:       Debanjan Mukherjee
# Institution:  University of California, Berkeley
# Last Edit:    June 2018
#---------------------------------------------------------------------------------------------
import sys, os
import vtk 
import numpy as np

try:
    from tracerInputModule import *
    from fileIO import *
except ImportError:
    sys.exit('Could not import user defined modules')
    
class LagrangianTracer(object):

    def __init__(self, a_InputData, a_DebugMode=False):

      self.m_InputData  = a_InputData
      self.m_DebugMode  = a_DebugMode
      self.m_Tracers    = vtk.vtkPolyData()
    
    ### IMPLEMENTATION INCOMPLETE
    def addTracerScalarData(self, a_ArrayName, a_DType='DOUBLE'):

        if a_DType == 'DOUBLE':
            data = vtk.vtkDoubleArray()

        data.SetName(a_ArrayName)
        data.SetNumberOfTuples(self.m_Particles.GetNumberOfPoints())
        data.SetNumberOfComponents(1)
    
    ### IMPLEMENTATION INCOMPLETE
    def addTracerVectorData(self, a_ArrayName, a_DType='DOUBLE'):

        if a_DType == 'DOUBLE':
            data = vtk.vtkDoubleArray()

        data.SetName(a_ArrayName)
        data.SetNumberOfTuples(self.m_Particles.GetNumberOfPoints())
        data.SetNumberOfComponents(self.m_InputData.getProblemDimension())
    
    ### IMPLEMENTATION INCOMPLETE
    def runCompute(self):

        #----------------------------------------------------------
        # initialize tracer data from an external VTK polydata file
        #----------------------------------------------------------
        self.m_Tracers  = initializeTracerFromFile(self.m_InputData.getTracerInput(), a_IsPolyDataForm=True)
        tracerInject    = initializeTracerFromFile(self.m_InputData.getTracerInput(), a_IsPolyDataForm=True)

        #---------------------------------------------------
        # set up a tracer data object for writing into files
        #---------------------------------------------------
        tracerOutput    = vtk.vtkPolyData()
        tracerWriter    = vtk.vtkXMLPolyDataWriter()

        #---------------------------------------------------------------------------
        # generate the time synchronization map for flow data and integration points
        #---------------------------------------------------------------------------
        timeWindowDict = self.m_InputData.getDataTimeWindows()

        if not(self.m_InputData.isSteadyFlowData()):
            self.m_InputData.printDataTimeWindows()

        #-------------------
        # start time counter
        #-------------------
        timeIndex = 0

    def integrateForwardEuler(self):
        print "TO BE IMPLEMENTED"
        


