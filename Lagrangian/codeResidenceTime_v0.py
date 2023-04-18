#---------------------------------------------------------------------------------------------
# This is a single-file minimal code version for calculating residence time from an 
# unsteady/steady simulated flow data implemented based on a simple code design philosophy that
# heavily makes use of the VTK library and numpy/scipy for all geometry and compute operations. 
# 
# Additional detailed features to be included in future versions
#
# Author:       Debanjan Mukherjee
# Institution:  University of California, Berkeley
# Last Edit:    June 2018
#---------------------------------------------------------------------------------------------
import sys, os
import vtk 
import numpy as np

from moduleTracerInput import *

#---------------------------------------------------------------------------
# Function to inject polydata particles at a specified location using a file
#---------------------------------------------------------------------------
def injectParticles(a_Particles, a_InputParticles):
    
    if ( (a_Particles.IsA('vtkPolyData') != 1) or (a_InputParticles.IsA('vtkPolyData') != 1) ):
        sys.exit("Module: moduleInjection, Function: injectParticles - both arguments should be polyData")
    
    if ( a_Particles.GetPointData().GetNumberOfArrays() != a_InputParticles.GetPointData().GetNumberOfArrays() ):
        sys.exit("Module: moduleInjection, Function: injectParticles - both arguments need same number of data")

    n_In    = a_Particles.GetNumberOfPoints()
    n_Add   = a_InputParticles.GetNumberOfPoints()

    for p in range(n_Add):

        xyz = a_InputParticles.GetPoint(p)
        a_Particles.GetPoints().InsertPoint(p + n_In, xyz)

    for a in range(a_Particles.GetPointData().GetNumberOfArrays()):

        for p in range(n_Add):
            
            daTuple = a_InputParticles.GetPointData().GetArray(a).GetTuple(p)
            a_Particles.GetPointData().GetArray(a).InsertTuple(p + n_In, daTuple)

    return a_Particles

#----------------------------------------------------------------------------------------------
# This is a modified version of tracer initialization routine,  where only the tracer injection 
# coordinates are read in from the file, but all other data arrays are configured separately.
#----------------------------------------------------------------------------------------------
def configureTracersFromFile(a_FileName, a_StartGlobalIDFrom=None, a_ScalarNamesZero=None):

    #
    # read the polydata object from the file
    #
    if a_FileName.endswith('vtk'):
        reader = vtk.vtkPolyDataReader()
    elif a_FileName.endswith('vtp'):
        reader = vtk.vtkXMLPolyDataReader()

    reader.SetFileName(a_FileName)
    reader.Update()

    polyData    = reader.GetOutput()
    nP          = polyData.GetNumberOfPoints()

    #
    # check whether velocity data was provided in the file
    # if yes:   then that is the velocity array initialized
    # if no:    create a vtkArray that will hold the velocity
    #           and initialize velocities based on inputs
    #
    if not(polyData.GetPointData().HasArray('velocity')):
        velocity = vtk.vtkDoubleArray()
        velocity.SetArrayName('velocity')
        velocity.SetNumberOfTuples(nP)
        velocity.SetNumberOfComponents(3)
        for p in range(nP):
            velocity.InsertTuple(p, (0.0,0.0,0.0))

        polyData.GetPointData().AddArray(velocity)

    #
    # check whether a global id array has been included in polydata
    # if yes:   then that array has to be overwritten
    # if no:    then a new id array has to be created
    #
    if not(polyData.GetPointData().HasArray('id')):

        idArray = vtk.vtkTypeInt32Array()
        idArray.SetName('id')
        idArray.SetNumberOfTuples(nP)
        idArray.SetNumberOfComponents(1)
        for p in range(nP):
            idArray.InsertTuple(p, p + a_StartGlobalIDFrom)

        polyData.GetPointData().AddArray(idArray)

    else:

        idArray = polyData.GetPointData().GetArray('id')
        for p in range(nP):
            idArray.InsertTuple(p, p + a_StartGlobalIDFrom)

    #
    # check whether a status array has been included in polydata
    # if yes:   then read that array straight from file
    # if no:    then initialize all status to 1 (that is active)
    #
    if not(polyData.GetPointData().HasArray('status')):

        statusArray = vtk.vtkTypeInt32Array()
        statusArray.SetName('status')
        statusArray.SetNumberOfTuples(nP)
        statusArray.SetNumberOfComponents(1)
        for p in range(nP):
            statusArray.InsertTuple(p, 1)

        polyData.GetPointData().AddArray(statusArray)
    
    #
    # all the remaining scalar data arrays are simply created and 
    # instantiated to a constant vector of zeros
    #
    for array in a_ScalarnamesZero:

        scalar = vtk.vtkDoubleArray()
        scalar.SetName(array)
        scalar.SetNumberOfTuples(nP)
        scalar.SetNumberOfComponents(1)

        for p in range(nP):
            scalar.InsertTuple(p, 0.0)

        polyData.GetPointData().AddArray(scalar)

    return polyData

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

#------------------------------------------------------------------------------
# check whether the tracer is within the ROI for the residence time calculation
# NOTE:
# - for now the ROI is simply a bounding box or a vtkLocator
# - further advancements to defining complex ROI need to be incorporated
#------------------------------------------------------------------------------
def isTracerInROI(a_X, a_ROI):

    if a_ROI.IsA('vtkLocator'):
        
        if a_ROI.FindCell(a_X) != -1:
            return True
        else:
            return False
        
    else:
        
        if (a_X[0] >= a_ROI[0]) and (a_X[0] <= a_ROI[1]) \
            and (a_X[1] >= a_ROI[2]) and (a_X[1] <= a_ROI[3]) \
            and (a_X[2] >= a_ROI[4]) and (a_X[2] <= a_ROI[5]):
        
        return True
    
    else:
        
        return False

#------------------------------------------------------------------------
# time integration routine for individual particles: Forward Euler Method
#------------------------------------------------------------------------
def integrateForwardEuler(a_Tracers, a_MeshLocator, a_ROI, a_Velocities, a_GridData, a_T, a_DT, a_Window, a_BoundaryCondition):
    
    #
    # as the tracer object enters the integration routine the tracer status data is 
    # taken in from the previous time-step to be used in this. the status data is 
    # only updated after the numerical integration update is completed
    #

    cellPtIds   = vtk.vtkIdList()
    resTimeArr  = a_Tracers.GetPointData().GetArray('residence')

    for p in range(a_Tracers.GetNumberOfPoints()):
        
        #
        # update the Lagrangian tracer status based on whether particle is in ROI
        #
        xyz         = a_Tracers.GetPoint(p)
        resT        = resTimeArr.GetTuple(p)
        cell        = a_Locator.FindCell(xyz)
        indFunc     = 1 if ( (cell != -1) and (isTracerInROI(xyz, a_ROI)) ) else 0

        if (cell != -1):
            
            if len(a_Velocities) == 2:

                a_GridData.GetCellPoints(cell, cellPtIds)
                
                velMinus_N1 = np.asarray(a_Velocities[0].GetTuple3(cellPtIds.GetId(0)))
                velMinus_N2 = np.asarray(a_Velocities[0].GetTuple3(cellPtIds.GetId(1)))
                velMinus_N3 = np.asarray(a_Velocities[0].GetTuple3(cellPtIds.GetId(2)))
                velMinus    = (1.0/3.0)*(velMinus_N1 + velMinus_N2 + velMinus_N3)

                velPlus_N1  = np.asarray(a_Velocities[1].GetTuple3(cellPtIds.GetId(0)))
                velPlus_N2  = np.asarray(a_Velocities[1].GetTuple3(cellPtIds.GetId(1)))
                velPlus_N3  = np.asarray(a_Velocities[1].GetTuple3(cellPtIds.GetId(2)))
                velPlus     = (1.0/3.0)*(velPlus_N1 + velPlus_N2 + velPlus_N3) 
                
                velInterp   = velMinus + (a_T - a_Window[0])*(velPlus - velMinus)/(a_Window[1] - a_Window[0])
                xyz         = xyz + velInterp*a_DT
                
                if isTracerInROI(xyz, a_ROI): 
                    resT = resT + a_DT
            
            elif len(a_Velocities) == 1:
                
                a_GridData.GetCellPoints(cell, cellPtIds)

                vel_N1  = np.asarray(a_Velocities[0].GetTuple3(cellPtIds.GetId(0)))
                vel_N2  = np.asarray(a_Velocities[0].GetTuple3(cellPtIds.GetId(1)))
                vel_N3  = np.asarray(a_Velocities[0].GetTuple3(cellPtIds.GetId(2)))

                vel     = (1.0/3.0)*(vel_N1 + vel_N2 + vel_N3)
                xyz     = xyz + vel*a_DT
                
                if isTracerInROI(xyz, a_ROI): 
                    resT = resT + a_DT

            a_Tracers.SetPoint(p, xyz[0], xyz[1], xyz[2])
            a_Tracers.GetPointData().GetArray('residence').InsertTuple(p, resT)

        elif (cell == -1):

            a_Tracers.GetPointData().GetArray('status').InsertTuple(p, -1)

    return a_Tracers

#--------------------------------------------------------------------------
# time integration routine for individual particles: Runge Kutta Integrator
#--------------------------------------------------------------------------
def integrateRK4(a_TracerPoints, a_Locator, a_Velocities, a_GridData, a_T, a_DT, a_Window, a_BoundaryCondition):

    return 0

#------------------
# main program loop
#------------------
if __name__=="__main__":

    #
    # parse command line argument
    #
    if len(sys.argv) != 2:
        sys.exit("Need Input Filename As An Argument")

    inputFile = sys.argv[1].strip()

    #
    # set-up problem input data
    #
    inputData = SimInputs(inputFile)
    inputData.readInputFile()
    
    #
    # initialize tracer data in the ROI from an external VTK polydata file
    #
    tracers = configureTracersFromFile(inputData.getTracerInput(), a_StartGlobalIDFrom=0, a_ScalarNamesZero=['residence'])
    
    #
    # set up a tracer data object for writing into files
    #
    tracerWriter = vtk.vtkXMLPolyDataWriter()

    #
    # based on the type of ROI build a locator for the ROI
    #
    roiObject = getResidenceTimeROI(inputData.getROIFile())

    #
    # generate the time synchronization map for flow data and integration points
    #
    timeWindowDict = inputData.getDataTimeWindows()
    
    if not(inputData.isSteadyFlowData()): inputData.printDataTimeWindows()

    #
    # start time counter
    #
    timeIndex   = 0
    simTime     = inputData.getSimulationStartTime()
    tWin_0      = timeWindowDict['T_Low'][0]
    tWin_1      = timeWindowDict['T_Up'][0]

    # Note:
    # -----
    # if inputData.isResidenceTimeComputeModeMapped():
    #       only integrate one launch - no additional injections
    #       write only the mapped files and/or write intermediate files
    # if inputData.isResidenceTimeComputeModeStreaming():
    #       perform injected particles' integration
    #       write intermediate files 
    # for both modes: calculate statistics using scipy unless requested

    #
    # start time integration loop
    #
    while simTime <= inputData.getSimulationStopTime():

        print "Integrating from", timeIndex, " to ", timeIndex + 1, "simTime", simTime

        #
        # resolve particle injection, and inject particles if needed
        #
        injectParticles = inputData.isResidenceTimeComputeModeMapped() and 

        #
        # create a boolean condition that governs when new flow-data files are loaded
        #
        if inputData.isSteadyFlowData():

            isLoadFrame     = (simTime == inputData.getSimulationStartTime())
            isDataPointSync = True

        else:
            
            isLoadFrame = (simTime == inputData.getSimulationStartTime()) \
                        or (timeWindowDict['T_Low'][timeIndex] != tWin_0 or timeWindowDict['T_Up'][timeIndex] != tWin_1)

            isDataPointSync = timeWindowDict['T_Low'][timeIndex] == timeWindowDict['T_Up'][timeIndex]
        
        #
        # a set of debug messages for program status
        #
        if isLoadFrame: 
            print "Will Load Velocity Data"
        
        if isDataPointSync:
            if inputData.isSteadyFlowData():
                if simTime == inputData.getSimulationStartTime():
                    print "Data And Integration Times Are Synced"
            else:
                print "Data And Integration Times Are Synced"
        
        #
        # load data file based on the evaluated boolean variable
        #
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

        #
        # build a cell locator (if it is a fixed mesh, build only once)
        #
        if simTime == inputData.getSimulationStartTime() and inputData.isFixedMesh():

            print "Building Cell Locator Maps"

            if isDataPointSync:
                locatorObj  = createCellLocator(flowSingle)         # syntax for standard vtkLocator
            else:
                locatorObj  = createCellLocator(flowMinus)          # syntax for standard vtkLocator
        
        #
        # now proceed with integration and update of each tracer point
        #
        if isDataPointSync:
            velocities = [velocitySingle]
        else:
            velocities = [velocityMinus, velocityPlus]
        
        boundaryCondition = 1
        tracerPoints = integrateForwardEuler(tracerPoints, locatorObj, velocities, gridDataObject, simTime,
                inputData.getIntegrationTimeStep(), [tWin_0, tWin_1], boundaryCondition)
        
        #
        # at the end of appropriate number of steps dump data into file if streaming calculation
        #
        if inputData.isResidenceTimeComputeModeStreaming():
            writeTracerDataToFile(tracerPoints, inputData.getTracerOutputFile(a_ID=timeIndex))
        
        #
        # update time indices and simulation time
        #
        timeIndex   = timeIndex + 1
        simTime     = simTime + inputData.getIntegrationTimeStep()

    
