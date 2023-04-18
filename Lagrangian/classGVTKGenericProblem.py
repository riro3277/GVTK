import sys, os
import numpy as np
import vtk

#---------------------------------------------------------------------------
# FOR VTK IMPLEMENTATION ONLY:
#-----------------------------
# this defines a base problem class that serves as an abstract class
# that can be extended to create specific problem or physics implementations
#---------------------------------------------------------------------------
class GVTKProblem(object):

    def __init__(self, a_InputData, a_DebugMode=False):

        self.m_InputData    = a_InputData
        self.m_DebugMode    = a_DebugMode
        self.m_Mode         = None

    @property
    def lagrangianData(self):
        return self.m_LagrangianData

    @lagrangianData.setter
    def lagrangianData(self, a_LagrangianObject):
        self.m_LagrangianData = a_LagrangianObject

    @property
    def gridData(self):
        return self.m_GridData

    @gridData.setter
    def gridData(self, a_GridDataObject):
        self.m_gridData = a_GridDataObject
