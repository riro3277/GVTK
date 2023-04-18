#----------------------------------------------------------------------------------------
# An implementation of the Euler Maruyama Method for solving a tracer advection-diffusion
# equation to create a Lagrangian Tracer Transport Model:
#
# Euler Maruayam Integration Scheme With Weiner Process:
# - dX = U(X,t)*dt + sqrt(2*D_0)*dB
# - U(X,t): velocity fields obtained from external CFD calculations
#
# NOTE: This version is using no vectorised code, mainly because it is trying to
# access the lists/data/coordinates using a VTK Object
#
# Author:       Debanjan Mukherjee
# Institution:  University of California, Berkeley
# Last Edited:  June 2018
#----------------------------------------------------------------------------------------
import sys, os
import vtk
import numpy as np
import numpy.random as nprand

from moduleTracerInput import *

#-----------------------------------
# BEGIN UTILITY FUNCTION DEFINITIONS
#-----------------------------------

#-------------------------------------------------------------------------
# Function to inject points at specified location using a specified method
# -------
# Params:
# -------
# a_Points:         vtkPoints object to which new points will be added
# a_InputPoints:    vtkPoints object that holds all injected points
# --------
# Returns:
# --------
# a_Points: a vtkPoints object with old and injected points
#-------------------------------------------------------------------------
def injectPoints(a_Points, a_InputPoints):

    oldNumPts   = a_Points.GetNumberOfPoints()
    addNumPts   = a_InputPoints.GetNumberOfPoints()

    for p in range(addNumPts):
        addID = oldNumPts + p
        a_Points.InsertPoint(addID, a_InputPoints.GetPoint(p))

    return a_Points

#---------------------------------------------------------------------------------
# Function to extract points from a specified filename, returns a vtkPoints object
# -------
# Params:
# -------
# a_FileName:   a vtk polydata file to read a collection of points/particles
# --------
# Returns:
# --------
# points:   a vtkPoints Object
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
# -------
# Params:
# -------
# a_FileName:   a vtk polydata file to read a collection of points/particles
#----------------------------------------------------------------------------------------
def initializeTracerFromFile(a_FileName):

    tracerInput = extractPointsFromFile(a_FileName)
    tracerPoints = vtk.vtkPoints()
    tracerPoints.SetNumberOfPoints(tracerInput.GetNumberOfPoints())
    for p in range(tracerPoints.GetNumberOfPoints()):
        tracerPoints.SetPoint(p, tracerInput.GetPoint(p))

    return tracerPoints

#-----------------------------------------------------------------------------------
# Function to extract mesh/grid data from a specified file. This is just a helper
# function to provide generic file read write capabilities.
# NOTE: MAKE THIS MORE GENERIC OR ELSE DEPRECATE THIS IN FUTURE VERSIONS
# -------
# Params:
# -------
# a_FileName:   a vtk unstructured grid file to read mesh/flow data from
# --------
# Returns:
# --------
# data: a vtkUnstructuredGrid data object
#-----------------------------------------------------------------------------------
def extractDataFromFile(a_FileName):

    if a_FileName.endswith('vtk'):
        reader = vtk.vtkUnstructuredGridReader()
    elif a_FileName.endswith('vtu'):
        reader = vtk.vtkXMLUnstructuredGridReader()

    reader.SetFileName(a_FileName)
    reader.Update()

    return reader.GetOutput()

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

#-------------------------------------------------------------------
# Function to write tracer data back to file
# This new version can accept both vtkPoints and vtkPolyData objects
# -------
# Params:
# -------
#-------------------------------------------------------------------
def writeTracerDataToFile(a_Tracers, a_OutputFileName):

    if a_OutputFileName.endswith('vtp'):
        tracerWriter = vtk.vtkXMLPolyDataWriter()
    else:
        tracerWriter = vtk.vtkPolyDataWriter()

    tracerWriter.SetFileName(a_OutputFileName)

    if a_Tracers.IsA('vtkPoints') == 1:

        tracerOutput = vtk.vtkPolyData()
        tracerOutput.SetPoints(a_Tracers)

        if vtk.VTK_MAJOR_VERSION <= 5.0:
            tracerWriter.SetInput(tracerOutput)
        else:
            tracerWriter.SetInputData(tracerOutput)

        tracerWriter.Update()
        tracerWriter.Write()

    elif a_Tracers.IsA('vtkPolyData') == 1:

        if vtk.VTK_MAJOR_VERSION <= 5.0:
            tracerWriter.SetInput(a_Tracers)
        else:
            tracerWriter.SetInputData(a_Tracers)

        tracerWriter.Update()
        tracerWriter.Write()

#--------------------------------------------------------------
# create a cell locator object from input mesh data from a file
#--------------------------------------------------------------
def createCellLocator(a_FileName, a_LocatorType=None):

    if a_FileName.endswith('vtu'):
        reader = vtk.vtkXMLUnstructuredGridReader()
    elif a_FileName.endswith('vtp'):
        reader = vtk.vtkXMLPolyDataReader()

    reader.SetFileName(a_FileName)
    reader.Update()

    if a_LocatorType is None:
        locator = vtk.vtkCellTreeLocator()
    else:
        if a_LocatorType == 'oct':
            locator = vtk.vtkCellLocator()
        elif a_LocatorType == 'tre':
            locator = vtk.vtkCellTreeLocator()
        elif a_LocatorType == 'bsp':
            locator = vtk.vtkModifiedBSPTree()

    locator.SetDataSet(reader.GetOutput())
    locator.BuildLocator()

    return locator

#-----------------------------------
# END UTILITY FUNCTION DEFINITIONS
#-----------------------------------

#-----------------------------------
# BEGIN INTEGRATOR DEFINITIONS
#-----------------------------------

#-----------------------------------------------------------------------------
# Function implementing a one-step Euler-Maruayam integrator for modeling the
# trajectory of Lagrangian tracers undergoing drift-diffusion processes
#-----------------------------------------------------------------------------
def eulerMaruyamaIntegrator(a_Particles, a_Locator, a_Velocities, a_GridData, a_T, a_DT, a_D0, a_Window, a_BoundaryCondition):

    for p in range(a_Particles.GetNumberOfPoints()):

        xyz         = a_Particles.GetPoint(p)
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

            vel   = velMinus + (a_T - a_Window[0])*(velPlus - velMinus)/(a_Window[1] - a_Window[0])

        elif len(a_Velocities) == 1:

            if cell != -1:

                a_GridData.GetCellPoints(cell, cellPtIds)

                vel_N1 = np.asarray(a_Velocities[0].GetTuple3(cellPtIds.GetId(0)))
                vel_N2 = np.asarray(a_Velocities[0].GetTuple3(cellPtIds.GetId(1)))
                vel_N3 = np.asarray(a_Velocities[0].GetTuple3(cellPtIds.GetId(2)))
                vel    = (1.0/3.0)*(vel_N1 + vel_N2 + vel_N3)

            else:
                vel     = np.array([0.0,0.0,0.0])

        xNew    = xyz[0] + vel[0]*a_DT + ((2.0*a_D0*a_DT)**0.5)*nprand.randn()
        yNew    = xyz[1] + vel[1]*a_DT + np.sqrt(2.0*a_D0*a_DT)*nprand.randn()

        #### AD-HOC
        #### Replace with a boundary condition function

        #if xNew > 0.060: xNew = 0.061

        #if xNew < -0.060: xNew = -0.061

        #if yNew > 0.015: yNew = 0.016

        #if yNew < -0.015: yNew = -0.016

        if xNew > 45.0: xNew = 45.01

        if xNew < 0.0:  xNew = -0.01

        if yNew > 6.0:  yNew = 6.01

        if yNew < 0.0:  yNew = -0.01

        #a_Particles.SetPoint(p, xyz[0], xyz[1], xyz[2])
        a_Particles.SetPoint(p, xNew, yNew, xyz[2])

    return a_Particles

#-----------------------------------
# END INTEGRATOR DEFINITIONS
#-----------------------------------

#-----------------------------------
# BEGIN COMPUTE SCRIPT DEFINITION
#-----------------------------------

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
                locatorObj  = createCellLocator(flowSingle)
            else:
                locatorObj  = createCellLocator(flowMinus)

        #--------------------------------------------------
        # now proceed with integration of each tracer point
        #--------------------------------------------------
        if isDataPointSync:
            velocities = [velocitySingle]
        else:
            velocities = [velocityMinus, velocityPlus]

        boundaryCondition = 1
        tracerPoints = eulerMaruyamaIntegrator(tracerPoints, locatorObj, velocities, gridDataObject, simTime,
                inputData.getIntegrationTimeStep(), inputData.getTracerDiffusivity() , [tWin_0, tWin_1], boundaryCondition)

        #--------------------------------------------------------------
        # at the end of appropriate number of steps dump data into file
        #--------------------------------------------------------------
        writeTracerDataToFile(tracerPoints, inputData.getTracerOutputFile(a_ID=timeIndex))

        #----------------------------------------
        # update time indices and simulation time
        #----------------------------------------
        timeIndex   = timeIndex + 1
        simTime     = simTime + inputData.getIntegrationTimeStep()

#-----------------------------------
# END COMPUTE SCRIPT DEFINITION
#-----------------------------------
