#---------------------------------------------------------------------------------------------
# This is a modified version of the simpleTracers code, which is designed based on 
# a simple code design philosophy that heavily makes use of the VTK library and numpy/scipy
# for all major geometry and compute operations. Tests on this code will help set-up a pilot
# study, as well as help transform this into a Fortran code equivalent for the SLIDES library
# and develop a plugin for particle tracer computation in SimVascular software suite.
# 
# Author:       Debanjan Mukherjee
# Institution:  University of California, Berkeley
# Last Edit:    02/19/2018
#---------------------------------------------------------------------------------------------
import sys, os
import vtk 
import numpy as np

from tracerInputModule import *

#-------------------------------------------------------------------------
# Function to inject points at specified location using a specified method
# NOTE: THIS IS THE DEPRECATED VERSION. DELETE THIS IN FUTURE VERSIONS.
#-------------------------------------------------------------------------
def injectPoints_OLD(a_Points, a_InputPoints, a_IsRandom=False):
    
    if a_IsRandom:
        randomY = np.random.uniform(0.0, 6.0, size=a_InputPoints.GetNumberOfPoints())
        
    for p in range(a_InputPoints.GetNumberOfPoints()):
        if a_IsRandom:
            a_Points.InsertNextPoint(0.0, randomY[p], 0.0)
        else:
            a_Points.InsertNextPoint(a_InputPoints.GetPoint(p))

    return a_Points

#-------------------------------------------------------------------------
# Function to inject points at specified location using a specified method
#-------------------------------------------------------------------------
def injectPoints(a_Points, a_InputPoints):
    
    oldNumPts   = a_Points.GetNumberOfPoints()
    addNumPts   = a_InputPoints.GetNumberOfPoints() 
        
    for p in range(addNumPts):
        addID = oldNumPts + p
        a_Points.InsertPoint(addID, a_InputPoints.GetPoint(p))

    return a_Points

#-------------------------------------------------------------------------
# Function to inject points at specific location from a chosen input file
# with initial cell locations seeded for faster/efficient particle-in-cell 
# computations. This is mandatory for any method for cell location other 
# than employing the standard VTK cell locator functionalities
#-------------------------------------------------------------------------
def injectPointsWithSeedLocations(a_Points, a_Locations, a_InputPoints, a_Locator):

    oldNumPts   = a_Points.GetNumberOfPoints()
    addNumPts   = a_InputPoints.GetNumberOfPoints()

    for p in range(addNumPts):
        addID



#---------------------------------------------------------------------------------
# Function to extract points from a specified filename, returns a vtkPoints object
#---------------------------------------------------------------------------------
def extractPointsFromFile(a_FileName):

    if a_FileName.endswith('vtk'):
        reader = vtk.vtkPolyDataReader()
    elif a_FileName.endsiwth('vtp'):
        reader = vtk.vtkXMLPolyDataReader()

    reader.SetFileName(a_FileName)
    reader.Update()
    points = reader.GetOutput().GetPoints()

    return points

#----------------------------------------------------------------------------------------
# Initialize all tracerData coordinates etc. from an input file for initial configuration
#----------------------------------------------------------------------------------------
def initializeTracerFromFile(a_FileName):

    tracerInput = extractPointsFromFile(a_FileName)
    tracerPoints = vtk.vtkPoints()
    tracerPoints.SetNumberOfPoints(tracerInput.GetNumberOfPoints())
    for p in range(tracerPoints.GetNumberOfPoints()):
        tracerPoints.SetPoint(p, tracerInput.GetPoint(p))

    return tracerPoints

#-----------------------------------------------------------------------------------
# Function to extract data from a specified filename, returns a vtkData(Set?) object
#-----------------------------------------------------------------------------------
def extractDataFromFile(a_FileName):

    if a_FileName.endswith('vtk'):
        reader = vtk.vtkUnstructuredGridReader()
    elif a_FileName.endswith('vtu'):
        reader = vtk.vtkXMLUnstructuredGridReader()

    reader.SetFileName(a_FileName)
    reader.Update()

    return reader.GetOutput()

#----------------------------------------------------------------------
# Function to write tracer data back to file
# NOTE: THIS IS THE DEPRECATED VERSION. DELETE THIS IN FUTURE VERSIONS.
#----------------------------------------------------------------------
def writeTracerDataToFile_OLD(a_TracerPoints, a_OutputFileName):

    if a_OutputFileName.endswith('vtp'):
        tracerWriter = vtk.vtkXMLPolyDataWriter()
    else:
        tracerWriter = vtk.vtkPolyDataWriter()

    tracerWriter.SetFileName(a_OutputFileName)
    
    tracerOutput = vtk.vtkPolyData()
    tracerOutput.SetPoints(a_TracerPoints)

    if vtk.VTK_MAJOR_VERSION <= 5.0:
        tracerWriter.SetInput(tracerOutput)
    else:
        tracerWriter.SetInputData(tracerOutput)

    tracerWriter.Update()
    tracerWriter.Write()

#-------------------------------------------------------------------
# Function to write tracer data back to file
# This new version can accept both vtkPoints and vtkPolyData objects
#-------------------------------------------------------------------
def writeTracerDataToFile(a_Tracers, a_OutputFileName):

    if a_OutputFileName.endswith('vtp'):
        tracerWriter = vtk.vtkXMLPolyDataWriter()
    else:
        tracerWriter = vtk.vtkPolyDataWriter()

    tracerWriter.SetFileName(a_OutputFileName)

    if a_TracerPoints.IsA('vtkPoints') == 1:
    
        tracerOutput = vtk.vtkPolyData()
        tracerOutput.SetPoints(a_Tracers)

        if vtk.VTK_MAJOR_VERSION <= 5.0:
            tracerWriter.SetInput(tracerOutput)
        else:
            tracerWriter.SetInputData(tracerOutput)

        tracerWriter.Update()
        tracerWriter.Write()

    elif a_TracerPoints.IsA('vtkPolyData') == 1:
        
        if vtk.VTK_MAJOR_VERSION <= 5.0:
            tracerWriter.SetInput(a_Tracers)
        else:
            tracerWriter.SetInputData(a_Tracers)

        tracerWriter.Update()
        tracerWriter.Write()

#----------------------------------------------------------------------------------------------
# Function to extract flow data from a specified filename, returns a vtkUnstructuredGrid object
#----------------------------------------------------------------------------------------------
def extractFlowDataFromFile(a_FileName, a_DataName=None):

    if a_FileName.endswith('vtk'):
        reader = vtk.vtkUnstructuredGridReader()
    elif a_FileName.endswith('vtu'):
        reader = vtk.vtkXMLUnstructuredGridReader()

    reader.SetFileName(a_FileName)
    reader.Update()

    if a_DataName is not None:
        return reader.GetOutput().GetPointData().GetArray(a_DataName)
    else:
        return reader.GetOutput()


#--------------------------------------------------------------
# create a cell locator object from input mesh data from a file
#--------------------------------------------------------------
def createCellLocator(a_FileName):

    if a_FileName.endswith('vtu'):
        reader = vtk.vtkXMLUnstructuredGridReader()
    elif a_FileName.endswith('vtp'):
        reader = vtk.vtkXMLPolyDataReader()
    
    reader.SetFileName(a_FileName)
    reader.Update()

    locator = vtk.vtkCellLocator()
    locator.SetDataSet(reader.GetOutput())
    locator.BuildLocator()
    
    return locator 

#------------------------------------------------------------------------
# time integration routine for individual particles: Forward Euler Method
#------------------------------------------------------------------------
def integrateForwardEuler(a_TracerPoints, a_Locator, a_Velocities, a_GridData, a_T, a_DT, a_Window, a_BoundaryCondition):
    
    for p in range(a_TracerPoints.GetNumberOfPoints()):
        
        xyz         = a_TracerPoints.GetPoint(p)
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
            xyz         = xyz + velInterp*a_DT

        elif len(a_Velocities) == 1:

            if cell != -1:
                a_GridData.GetCellPoints(cell, cellPtIds)

                vel_N1 = np.asarray(a_Velocities[0].GetTuple3(cellPtIds.GetId(0)))
                vel_N2 = np.asarray(a_Velocities[0].GetTuple3(cellPtIds.GetId(1)))
                vel_N3 = np.asarray(a_Velocities[0].GetTuple3(cellPtIds.GetId(2)))
                vel    = (1.0/3.0)*(vel_N1 + vel_N2 + vel_N3)
            else:
                vel     = np.array([0.0,0.0,0.0])

            xyz = xyz + vel*a_DT

        #### AD-HOC
        #### REPLACE WITH BOUNDARY CONDITION FUNCTION
            
        #if xyz[0] > 0.060: xyz[0] = 0.061
        
        #if xyz[0] < -0.060: xyz[0] = -0.061
        
        #if xyz[1] > 0.015: xyz[1] = 0.016
        
        #if xyz[1] < -0.015: xyz[1] = -0.016

        if xyz[0] > 45.0: xyz[0] = 45.01
        
        if xyz[0] < 0.0: xyz[0] = -0.01
        
        if xyz[1] > 6.0: xyz[1] = 6.01
        
        if xyz[1] < 0.0: xyz[1] = -0.01

        a_TracerPoints.SetPoint(p, xyz[0], xyz[1], xyz[2])

    return a_TracerPoints

#--------------------------------------------------------------------------
# time integration routine for individual particles: Runge Kutta Integrator
#--------------------------------------------------------------------------
def integrateRK4(a_TracerPoints, a_Locator, a_Velocities, a_DT, a_DataTimeWindow, a_BoundaryCondition):

    return 0

#------------------
# main program loop
#------------------
if __name__=="__main__":

    #----------------------------
    # parse command line argument
    #----------------------------
    if len(sys.argv) != 2:
        sys.exit("Need Input Filename As An Argument")

    inputFile = sys.argv[1].strip()

    #--------------------------
    # set-up problem input data
    #--------------------------
    inputData = SimInputs(inputFile)
    inputData.readInputFile()
    
    #----------------------------------------------------------
    # initialize tracer data from an external VTK polydata file
    #----------------------------------------------------------
    tracerPoints = initializeTracerFromFile(inputData.getTracerInput())
    tracerInject = initializeTracerFromFile(inputData.getTracerInput())

    #---------------------------------------------------
    # set up a tracer data object for writing into files
    #---------------------------------------------------
    tracerOutput    = vtk.vtkPolyData()
    tracerWriter    = vtk.vtkXMLPolyDataWriter()

    #---------------------------------------------------------------------------
    # generate the time synchronization map for flow data and integration points
    #---------------------------------------------------------------------------
    timeWindowDict = inputData.getDataTimeWindows()
    
    if not(inputData.isSteadyFlowData()):
        inputData.printDataTimeWindows()

    #-------------------
    # start time counter
    #-------------------
    timeIndex   = 0
    simTime     = inputData.getSimulationStartTime()
    tWin_0      = timeWindowDict['T_Low'][0]
    tWin_1      = timeWindowDict['T_Up'][0]

    while simTime <= inputData.getSimulationStopTime():

        print "Integrating from", timeIndex, " to ", timeIndex + 1, "simTime", simTime

        #-----------------------------------------------------------------------
        # create a boolean condition that governs when new data files are loaded
        #-----------------------------------------------------------------------
        if inputData.isSteadyFlowData():

            isLoadFrame     = (simTime == inputData.getSimulationStartTime())
            isDataPointSync = True

        else:
            
            isLoadFrame = (simTime == inputData.getSimulationStartTime()) \
                        or (timeWindowDict['T_Low'][timeIndex] != tWin_0 or timeWindowDict['T_Up'][timeIndex] != tWin_1)

            isDataPointSync = timeWindowDict['T_Low'][timeIndex] == timeWindowDict['T_Up'][timeIndex]
        
        #---------------------------------------------------------------------------
        # create a boolean condition that governs when new points are to be injected
        #---------------------------------------------------------------------------
        isInjectPoints = isLoadFrame    ### THIS NEEDS MORE POLISHING

        #-------------------------------------------
        # a set of debug messages for program status
        #-------------------------------------------
        if isLoadFrame: 
            print "Will Load Velocity Data"
        
        if isDataPointSync:
            if inputData.isSteadyFlowData():
                if simTime == inputData.getSimulationStartTime():
                    print "Data And Integration Times Are Synced"
            else:
                print "Data And Integration Times Are Synced"

        if isInjectPoints: 
            print "New Particles Injected"
            print "Now Integrating", tracerPoints.GetNumberOfPoints(), "Particles"

        #-----------------------------------------
        # inject tracers into the domain if needed
        #-----------------------------------------
        if isInjectPoints:
            print "injecting here"
            tracerPoints = injectPoints(tracerPoints, tracerInject) ### THIS NEEDS MORE WORK
        
        #-------------------------------------------------------
        # load data file based on the evaluated boolean variable
        #-------------------------------------------------------
        if isLoadFrame:
            
            if inputData.isSteadyFlowData():
                
                flowSingle      = inputData.getFlowDataFileName()
                velocitySingle  = extractFlowDataFromFile(flowSingle, a_DataName=inputData.getVelDataName())
                gridDataObject  = extractDataFromFile(flowSingle)

            else:

                tWin_0  = timeWindowDict['T_Low'][timeIndex]    # updating the time window
                tWin_1  = timeWindowDict['T_Up'][timeIndex]     # updating the time window

                if isDataPointSync:
                    flowSingle      = inputData.getFlowDataFileName(a_ID = timeWindowDict['ID_Low'][timeIndex])
                    velocitySingle  = extractFlowDataFromFile(flowSingle, a_DataName=inputData.getVelDataName())
                    gridDataObject  = extractDataFromFile(flowSingle)        
                else:
                    flowPlus        = inputData.getFlowDataFileName(a_ID = timeWindowDict['ID_Up'][timeIndex])
                    flowMinus       = inputData.getFlowDataFileName(a_ID = timeWindowDict['ID_Low'][timeIndex])
                    velocityPlus    = extractFlowDataFromFile(flowPlus,  a_DataName=inputData.getVelDataName())
                    velocityMinus   = extractFlowDataFromFile(flowMinus, a_DataName=inputData.getVelDataName())
                    gridDataObject  = extractDataFromFile(flowMinus)

        #--------------------------------------------------------------
        # build a cell locator (if it is a fixed mesh, built only once)
        #--------------------------------------------------------------
        if simTime == inputData.getSimulationStartTime() and inputData.isFixedMesh():

            print "Building Cell Locator Maps"

            if isDataPointSync:
                locatorObj  = createCellLocator(flowSingle)         # syntax for standard vtkLocator
                #locatorObj  = TriangleLocator(a_File=flowSingle)    # syntax for custom locators
            else:
                locatorObj  = createCellLocator(flowMinus)          # syntax for standard vtkLocator
                #locatorObj  = TriangleLocator(a_File=flowMinus)     # syntax for custom locators
        
        #--------------------------------------------------
        # now proceed with integration of each tracer point
        #--------------------------------------------------
        if isDataPointSync:
            velocities = [velocitySingle]
        else:
            velocities = [velocityMinus, velocityPlus]
        
        boundaryCondition = 1
        tracerPoints = integrateForwardEuler(tracerPoints, locatorObj, velocities, gridDataObject, simTime,
                inputData.getIntegrationTimeStep(), [tWin_0, tWin_1], boundaryCondition)
        
        #--------------------------------------------------------------
        # at the end of appropriate number of steps dump data into file
        #--------------------------------------------------------------
        writeTracerDataToFile(tracerPoints, inputData.getTracerOutputFile(a_ID=timeIndex))
        
        #----------------------------------------
        # update time indices and simulation time
        #----------------------------------------
        timeIndex   = timeIndex + 1
        simTime     = simTime + inputData.getIntegrationTimeStep()
