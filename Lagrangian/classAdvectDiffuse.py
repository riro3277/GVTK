#---------------------------------------------------------------------------------------------
# This is a module that encapsulates the physics and numerics of a Lagrangian tracer based 
# advection diffusion calculation
#
# The main objective here is to perform a drift-diffusion type calculation for (nearly) zero
# inertia particles/tracers, and the numerical integration is implemented using a simple
# Euler Maruyama integration scheme
# 
# Author:       Debanjan Mukherjee
# Institution:  University of California, Berkeley
# Last Edit:    June 2018
#---------------------------------------------------------------------------------------------
import sys, os
import vtk
import numpy as np
import numpy.random as nprand

class EulerMaruyama(object):
    
    def __init__(self, a_InputData, a_DebugMode=False):

        self.m_InputData    = a_InputData
        self.m_DebugMode    = a_Debugmode
        self.m_Particles    = None
        self.m_Diffusivity  = self.m_InputData.getTracerDiffusivity()

    def runCompute(self):
        
        #----------------------------------------------------------
        # initialize tracer data from an external VTK polydata file
        #----------------------------------------------------------
        self.m_Tracers  = initializeTracerFromFile(self.m_InputData.getTracerInput())
        tracerInject    = initializeTracerFromFile(self.m_InputData.getTracerInput())
        
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
        timeIndex   = 0
        simTime     = self.m_InputData.getSimulationStartTime()
        tWin_0      = timeWindowDict['T_Low'][0]
        tWin_1      = timeWindowDict['T_Up'][0]
        
        while simTime <= self.m_InputData.getSimulationStopTime():
            
            print "Integrating from", timeIndex, " to ", timeIndex + 1, "simTime", simTime
            
            #-----------------------------------------------------------------------
            # create a boolean condition that governs when new data files are loaded
            #-----------------------------------------------------------------------
            if self.m_InputData.isSteadyFlowData():
                
                isLoadFrame     = (simTime == self.m_InputData.getSimulationStartTime())
                isDataPointSync = True
            
            else:
                
                isLoadFrame = (simTime == self.m_InputData.getSimulationStartTime()) \
                        or (timeWindowDict['T_Low'][timeIndex] != tWin_0 or timeWindowDict['T_Up'][timeIndex] != tWin_1)
                        
                isDataPointSync = timeWindowDict['T_Low'][timeIndex] == timeWindowDict['T_Up'][timeIndex]
                
            #---------------------------------------------------------------------------
            # create a boolean condition that governs when new points are to be injected
            #---------------------------------------------------------------------------
            #isInjectPoints  = isLoadFrame    ### THIS NEEDS MORE POLISHING
            #isInjectPoints  = True
            
            if simTime <= 0.5: 
                isInjectPoints = True
            else:
                isInjectPoints = False
                
            #-------------------------------------------
            # a set of debug messages for program status
            #-------------------------------------------
            if isLoadFrame: print "Will Load Velocity Data"
                
            if isDataPointSync:
                if self.m_InputData.isSteadyFlowData():
                    if simTime == self.m_InputData.getSimulationStartTime():
                        print "Data And Integration Times Are Synced"
                else:
                    print "Data And Integration Times Are Synced"
                
            #-----------------------------------------
            # inject tracers into the domain if needed
            #-----------------------------------------
            if isInjectPoints:
                self.m_Tracers  = injectPoints(self.m_Tracers, tracerInject) ### THIS NEEDS MORE WORK
                #self.m_Tracers = injectPointsWithRandomPerturbations(self.m_Tracers, tracerInject) 
                
                print "New Particles Injected"
                print "Now Integrating", self.m_Tracers.GetNumberOfPoints(), "Particles"
                
            #-------------------------------------------------------
            # load data file based on the evaluated boolean variable
            #-------------------------------------------------------
            if isLoadFrame:
                
                if self.m_InputData.isSteadyFlowData():
                    
                    flowSingle      = self.m_InputData.getFlowDataFileName()
                    velocitySingle  = extractFlowDataFromFile(flowSingle, a_DataName=self.m_InputData.getVelDataName())
                    gridDataObject  = extractDataFromFile(flowSingle)
                
                else:
                    
                    tWin_0  = timeWindowDict['T_Low'][timeIndex]    # updating the time window
                    tWin_1  = timeWindowDict['T_Up'][timeIndex]     # updating the time window
                    
                    if isDataPointSync:
                        
                        flowSingle      = self.m_InputData.getFlowDataFileName(a_ID = timeWindowDict['ID_Low'][timeIndex])
                        velocitySingle  = extractFlowDataFromFile(flowSingle, a_DataName=self.m_InputData.getVelDataName())
                        gridDataObject  = extractDataFromFile(flowSingle)
                    
                    else:
                        
                        flowPlus        = self.m_InputData.getFlowDataFileName(a_ID = timeWindowDict['ID_Up'][timeIndex])
                        flowMinus       = self.m_InputData.getFlowDataFileName(a_ID = timeWindowDict['ID_Low'][timeIndex])
                        velocityPlus    = extractFlowDataFromFile(flowPlus,  a_DataName=self.m_InputData.getVelDataName())
                        velocityMinus   = extractFlowDataFromFile(flowMinus, a_DataName=self.m_InputData.getVelDataName())
                        gridDataObject  = extractDataFromFile(flowMinus)
                        
            #--------------------------------------------------------------
            # build a cell locator (if it is a fixed mesh, built only once)
            #--------------------------------------------------------------
            if simTime == self.m_InputData.getSimulationStartTime() and self.m_InputData.isFixedMesh():
                
                print "Building Cell Locator Maps"
                
                if isDataPointSync:
                    locatorObj  = createCellLocator(flowSingle)         # syntax for standard vtkLocator
                else:
                    locatorObj  = createCellLocator(flowMinus)          # syntax for standard vtkLocator
                    
            #--------------------------------------------------
            # now proceed with integration of each tracer point
            #--------------------------------------------------
            if isDataPointSync:
                velocities = [velocitySingle]
            else:
                velocities = [velocityMinus, velocityPlus]
                
            boundaryCondition = 1   # currently a dummy variable
            
            self.integrateEulerMaruyama(locatorObj, velocities, gridDataObject, simTime, [tWin_0,tWin_1], boundaryCondition)
                
            #--------------------------------------------------------------
            # at the end of appropriate number of steps dump data into file
            #--------------------------------------------------------------
            writeTracerDataToFile(self.m_Tracers, self.m_InputData.getTracerOutputFile(a_ID=timeIndex))
            
            #----------------------------------------
            # update time indices and simulation time
            #----------------------------------------
            timeIndex   = timeIndex + 1
            simTime     = simTime + self.m_InputData.getIntegrationTimeStep()

        #---------------------------------------------------------------------------------------------
        # first-order single step drift-diffusion integrator using a random walk Euler Maruyama scheme
        #---------------------------------------------------------------------------------------------
        def eulerMaruyamaIntegrator(self, a_Locator, a_Velocities, a_GridData, a_T, a_Window, a_BoundaryCondition):

            dT  = self.m_InputData.getIntegrationTimeStep()

            for p in range(self.m_Particles.GetNumberOfPoints()):

                xyz         = self.m_Particles.GetPoint(p)
                cell        = a_Locator.FindCell(xyz)
                cellPtIds   = vtk.vtkIdList()

                if len(a_Velocities) == 2:

                    if cell != -1:
                        
                        a_GridData.GetCellPoints(cell, cellPtIds)

                        velMinus_N1 = np.asarray(a_Velocities[0].GetTuple3(cellPtIds.GetId(0)))
                        velMinus_N2 = np.asarray(a_Velocities[0].GetTuple3(cellPtIds.GetId(1)))
                        velMinus_N3 = np.asarray(a_Velocities[0].GetTuple3(cellPtIds.GetId(2)))
                        velMinus    = (1.0/3.0)*(velMinus_N1 + velMinus_N2 + velMinus_N3)

                        velPlus_N1  = np.asarray(a_Velocities[1].GetTuple3(cellPtIds.GetId(0)))
                        velPlus_N2  = np.asarray(a_Velocities[1].GetTuple3(cellPtIds.GetId(1)))
                        velPlus_N3  = np.asarray(a_Velocities[1].GetTuple3(cellPtIds.GetId(2)))
                        velPlus     = (1.0/3.0)*(velPlus_N1 + velPlus_N2 + velPlus_N3)

                    else:

                        velMinus    = np.array([0.0, 0.0, 0.0])
                        velPlus     = np.array([0.0, 0.0, 0.0])

                    vel = velMinus + (a_T - a_Windo[0])*(velPlus - velMinus)/(a_Window[1] - a_Window[0])

                elif len(a_Velocities) == 1:

                    if cell != -1:

                        a_GridData.GetCellPoints(cell, cellPtIds)
                        
                        vel_N1 = np.asarray(a_Velocities[0].GetTuple3(cellPtIds.GetId(0)))
                        vel_N2 = np.asarray(a_Velocities[0].GetTuple3(cellPtIds.GetId(1)))
                        vel_N3 = np.asarray(a_Velocities[0].GetTuple3(cellPtIds.GetId(2)))
                        vel    = (1.0/3.0)*(vel_N1 + vel_N2 + vel_N3)

                    else:

                        vel     = np.array([0.0, 0.0, 0.0])

                xNew    = xyz[0] + vel[0]*dT + np.sqrt(2.0*self.m_Diffusivity*dT)*nprand.randn()
                yNew    = xyz[1] + vel[1]*dT + np.sqrt(2.0*self.m_Diffusivity*dT)*nprand.randn()

                if self.m_InputData.getProblemDimension() == 3:
                    zNew = xyz[2] + vel[2]*dT + np.sqrt(2.0*self.m_Diffusivity*dT)*nprand.randn()
                else:
                    zNew = xyz[2]

                #### AD-HOC REPLACE WITH BOUNDARY CONDITION FUNCTION
                #### NEED TO IMPROVE BOUNDARY CONDITION IMPLEMENTATION

                self.m_Particles.SetPoint(p, xNew, yNew, zNew)





