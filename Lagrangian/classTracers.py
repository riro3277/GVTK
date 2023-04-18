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
    from tracerInput import *
    from moduleFileIO import *
    from moduleVTKLocators import *
except ImportError:
    sys.exit('Could not import user defined modules')

class SimpleTracer(object):

    def __init__(self, a_InputData, a_DebugMode=False):

        self.m_InputData    = a_InputData
        self.m_DebugMode    = a_DebugMode
        self.m_Tracers      = None

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

            print("Integrating from", timeIndex, " to ", timeIndex + 1, "simTime", simTime)

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
            isInjectPoints  = False

            #if simTime <= 0.5:
            #    isInjectPoints = True
            #else:
            #    isInjectPoints = False

            #-------------------------------------------
            # a set of debug messages for program status
            #-------------------------------------------
            if isLoadFrame: print("Will Load Velocity Data")

            if isDataPointSync:
                if self.m_InputData.isSteadyFlowData():
                    if simTime == self.m_InputData.getSimulationStartTime():
                        print("Data And Integration Times Are Synced")
                else:
                    print("Data And Integration Times Are Synced")

            #-----------------------------------------
            # inject tracers into the domain if needed
            #-----------------------------------------
            if isInjectPoints:
                self.m_Tracers  = injectPoints(self.m_Tracers, tracerInject) ### THIS NEEDS MORE WORK
                #self.m_Tracers = injectPointsWithRandomPerturbations(self.m_Tracers, tracerInject)

                print("New Particles Injected")
                print("Now Integrating", self.m_Tracers.GetNumberOfPoints(), "Particles")

            #-------------------------------------------------------
            # load data file based on the evaluated boolean variable
            #-------------------------------------------------------
            if isLoadFrame:

                if self.m_InputData.isSteadyFlowData():

                    flowSingle      = self.m_InputData.getFlowDataFileName(a_Foam=True)
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

                print("Building Cell Locator Maps")

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

            if self.m_InputData.getIntegrationScheme() == 'feu':
                self.integrateForwardEuler(locatorObj, velocities, gridDataObject, simTime, [tWin_0,tWin_1], boundaryCondition)

            #--------------------------------------------------------------
            # at the end of appropriate number of steps dump data into file
            #--------------------------------------------------------------
            # print("Writing to file", self.m_InputData.getTracerOutputFile(a_ID1=0, a_ID2=timeIndex))
            if timeIndex%10 == 0:
                writeTracerDataToFile(self.m_Tracers, self.m_InputData.getTracerOutputFile(a_ID1=0, a_ID2=timeIndex))

            #----------------------------------------
            # update time indices and simulation time
            #----------------------------------------
            timeIndex   = timeIndex + 1
            simTime     = simTime + self.m_InputData.getIntegrationTimeStep()

    #--------------------------------------------------------------------------------------
    # time integration routine for individual particles: Forward Euler Method modified such
    # that this becomes a instanced method of the class instead of being a static method
    #--------------------------------------------------------------------------------------
    def integrateForwardEuler(self, a_Locator, a_Velocities, a_GridData, a_T, a_Window, a_BoundaryConditionn):

        for p in range(self.m_Tracers.GetNumberOfPoints()):

            xyz         = self.m_Tracers.GetPoint(p)
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

                    velMinus    = np.array([0.0,0.0,0.0])
                    velPlus     = np.array([0.0,0.0,0.0])

                velInterp   = velMinus + (a_T - a_Window[0])*(velPlus - velMinus)/(a_Window[1] - a_Window[0])
                xyz         = xyz + velInterp*self.m_InputData.getIntegrationTimeStep()

            elif len(a_Velocities) == 1:

                if cell != -1:

                    a_GridData.GetCellPoints(cell, cellPtIds)

                    vel_N1 = np.asarray(a_Velocities[0].GetTuple3(cellPtIds.GetId(0)))
                    vel_N2 = np.asarray(a_Velocities[0].GetTuple3(cellPtIds.GetId(1)))
                    vel_N3 = np.asarray(a_Velocities[0].GetTuple3(cellPtIds.GetId(2)))

                    vel    = (1.0/3.0)*(vel_N1 + vel_N2 + vel_N3)

                else:

                    vel     = np.array([0.0,0.0,0.0])

                xyz = xyz + vel*self.m_InputData.getIntegrationTimeStep()

                #### AD-HOC REPLACE WITH BOUNDARY CONDITION FUNCTION
                #### NEED TO IMPROVE BOUNDARY CONDITION IMPLEMENTATION

            self.m_Tracers.SetPoint(p, xyz[0], xyz[1], xyz[2])
