#---------------------------------------------------------------------------------------------
# This is a module that encapsulates the physics and numerics of a light tracer integration
# calculation through a standard specified analytical flow field. 
#
# This module has primarily been included such that numerical analysis of theoretical 
# understanding of tracer dynamics with canonical flow fields can be enabled by the library.
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
    
class TraceAnalytical(object):

    def __init__(self, a_InputData, a_DebugMode=False):

        self.m_InputData    = a_InputData
        self.m_DebugMode    = a_DebugMode
        self.m_Tracers      = None
    
    def runCompute(self):

        #-------------------------------------------------------------------------
        # set-up analytical flow to be integrated from the analytical flow library
        #-------------------------------------------------------------------------
        flow = ABCFlow(a_A=np.sqrt(3.0), a_B=np.sqrt(2.0), a_C=1.0)

        #----------------------------------------------------------
        # initialize tracer data from an external VTK polydata file
        #----------------------------------------------------------
        self.m_Tracers  = initializeTracerFromFile(self.m_InputData.getTracerInput())
        tracerInject    = initializeTracerFromFile(self.m_InputData.getTracerInput())

        #---------------------------------------------------
        # set-up a tracer data object for writing into files
        #---------------------------------------------------
        tracerOutput    = vtk.vtkPolyData()
        tracerWriter    = vtk.vtkXMLPolyDataWriter()

        #-------------------
        # start time counter
        #-------------------
        timeIndex   = 0
        simTime     = self.m_InputData.getSimulationStartTime()

        while simTime <= self.m_InputData.getSimulationStopTime():

            print "Integrating from", timeIndex, " to ", timeIndex+1, "simTime", simTime

            #------------------------------------------------------------------
            # boolean condition that governs when new points are to be injected
            #------------------------------------------------------------------
            isInjectPoints = True

            #-----------------------------------------
            # inject tracers into the domain if needed
            #-----------------------------------------
            if isInjectPoints:
                self.m_Tracers  = injectPointsWithRandomPerturbations(self.m_Tracers, tracerInject)
                print "New particles injected"
                print "Now integrating", tracerPoints.GetNumberOfPoints(), "Particles"

            boundaryCondition = 1   # currently a dummy variable

            if self.m_InputData.getIntegrationScheme() == 'feu':
                self.integrateForwardEuler(flow, simTime, boundaryCondition)
            elif self.m_InputData.getIntegrationScheme() == 'rk4':
                self.integrateRK4(flow, simTime, boundaryCondition)
            
            #---------------------------------------------------------------------
            # at the end of appropriate number of time steps dump data into a file
            #---------------------------------------------------------------------
            writeTracerDataToFile(self.m_Tracers, self.m_InputData.getTracerOutputFile(a_ID=timeIndex))

            #----------------------------------------
            # update time indices and simulation time
            #----------------------------------------
            timeIndex   = timeIndex + 1
            simTime     = simTime + self.m_InputData.getIntegrationTimeStep() 

    #--------------------------------------------------------------------------------------
    # time integration routine for individual particles: Forward Euler Method modified such 
    # that this becomes a instanced method of the class instead of being a static method
    #--------------------------------------------------------------------------------------
    def integrateForwardEuler(self, a_FlowObject, a_T, a_BoundaryCondition):
        
        dT  = self.m_InputData.getIntegrationTimeStep()
        
        for p in range(self.m_Tracers.GetNumberOfPoints()):

            xyz = self.m_Tracers.GetPoint(p)
            uvw = a_FlowObject.eval(xyz[0], xyz[1], xyz[2], a_T)
            x   = xyz[0] + uvw[0]*dT
            y   = xyz[1] + uvw[1]*dT
            z   = xyz[2] + uvw[2]*dT
            
            #### AD-HOC REPLACE WITH BOUNDARY CONDITION FUNCTION
            #### NEED TO IMPROVE BOUNDARY CONDITION IMPLEMENTATION

            self.m_Tracers.SetPoint(p, xyz[0], xyz[1], xyz[2])
    
    #--------------------------------------------------------------------------------------
    # time integration routine for individual particles: 4-stage Runge Kutta Method modified 
    # such that this becomes a instanced method of the class instead of being a static method
    #--------------------------------------------------------------------------------------
    def integrateRK4(self, a_FlowObject, a_T, a_BoundaryCondition):

        dT  = self.m_InputData.getIntegrationTimeStep()

        for p in range(self.m_Tracers.GetNumberOfPoints()):
            
            xyz = a_TracerPoints.GetPoint(p)
            k1  = a_FlowObject.eval(xyz[0],             xyz[1],             xyz[2],             a_T)
            k2  = a_FlowObject.eval(xyz[0] + 0.5*dT*k1, xyz[1] + 0.5*dT*k1, xyz[2] + 0.5*dT*k1, a_T + 0.5*dT)
            k3  = a_FlowObject.eval(xyz[0] + 0.5*dT*k2, xyz[1] + 0.5*dT*k2, xyz[2] + 0.5*dT*k2, a_T + 0.5*dT)
            k4  = a_FlowObject.eval(xyz[0] + dT*k3,     xyz[1] + dT*k3,     xyz[2] + dT*k3,     a_T + dT)
            x   = xyz[0] + (dT/6.0)*(k1 + 2.0*k2 + 2.0*k3 + k4)
            y   = xyz[1] + (dT/6.0)*(k1 + 2.0*k2 + 2.0*k3 + k4)
            z   = xyz[2] + (dT/6.0)*(k1 + 2.0*k2 + 2.0*k3 + k4)

            self.m_Tracers.SetPoint(p, x, y, z)


