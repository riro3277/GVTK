#---------------------------------------------------------------------------------------------
# This is a modified version of the simpleTracers code, which is designed based on 
# a simple code design philosophy that heavily makes use of the VTK library and numpy/scipy
# for all major geometry and compute operations. 
# 
# The objective of this module is to:
# - simply perform an advection computation for tracers released in a vector field
# - develop a plugin for particle tracer computation in SimVascular software suite.
# 
# Additional detailed features to be included in future versions
#
# Author:       Debanjan Mukherjee
# Institution:  University of California, Berkeley
# Last Edit:    03/19/2018
#---------------------------------------------------------------------------------------------
import sys, os
import vtk 
import numpy as np

from tracerInputModule import *
from analyticalFlowLibrary import *

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

#------------------------------------------------------------------------
# A utility function to delete points from an ensemble based on a list of 
# point ids or index list provided to the function
# NOTE: To be finalized
#------------------------------------------------------------------------
def deletePoints(a_Points, a_IdsToDelete):

  newPoints = vtk.vtkPoints()
  numPoints = a_Points.GetNumberOfPoints()

  check = lambda x : (isinstance(a_IdsToDelete, int) and x == a_IdsToDelete) or (isinstance(a_IdsToDelete, list) and (x in a_IdsToDelete))
  
  for id in range(numPoints):
    if check(id) == False:
      xyz = a_Points.GetPoint(id)
      newPoints.InsertNextPoint(xyz)

  a_Points.ShallowCopy(newPoints)

  return a_Points

#-------------------------------------------------------------------------------------------
# This generates coordinates of points/tracers injected into a domain such that to the 
# original set of points (a_Points) a new set of points are added (a_InputPoints) but 
# prior to adding these new points are given a random perturbation of 1+a_PerturbationExtent
# You can always set the a_PerturbationExtent to 0.0 to make sure only the a_InputPoints 
# coordinates are directly added into the domain
#-------------------------------------------------------------------------------------------
def injectPointsWithRandomPerturbations(a_Points, a_InputPoints, a_PerturbationExtent=0.1):

    oldNumPts   = a_Points.GetNumberOfPoints()
    addNumPts   = a_InputPoints.GetNumberOfPoints()
    
    for p in range(addNumPts):
        addID   = oldNumPts + p
        xyz     = np.asarray(a_InputPoints.GetPoint(p))
        xyz[0]  = xyz[0]*(1.0 + np.random.uniform(low=-a_PerturbationExtent, high=a_PerturbationExtent))
        xyz[1]  = xyz[1]*(1.0 + np.random.uniform(low=-a_PerturbationExtent, high=a_PerturbationExtent))
        xyz[2]  = xyz[2]*(1.0 + np.random.uniform(low=-a_PerturbationExtent, high=a_PerturbationExtent))
        a_Points.InsertPoint(addID, xyz)

    return a_Points

#-------------------------------------------------------------------------
# Function to inject points at specific location from a chosen input file
# with initial cell locations seeded for faster/efficient particle-in-cell 
# computations. This is mandatory for any method for cell location other 
# than employing the standard VTK cell locator functionalities
#
# TO BE IMPLEMENTED
#-------------------------------------------------------------------------
def injectPointsWithSeedLocations(a_Points, a_Locations, a_InputPoints, a_Locator):

    oldNumPts   = a_Points.GetNumberOfPoints()
    addNumPts   = a_InputPoints.GetNumberOfPoints()


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

#------------------------------------------------------------------------
# time integration routine for individual particles: Forward Euler Method
#------------------------------------------------------------------------
def integrateForwardEuler(a_TracerPoints, a_FlowObject, a_T, a_DT, a_BoundaryCondition):
    
    for p in range(a_TracerPoints.GetNumberOfPoints()):
        
        xyz = a_TracerPoints.GetPoint(p)
        uvw = a_FlowObject.eval(xyz[0], xyz[1], xyz[2], a_T)
        x   = xyz[0] + uvw[0]*a_DT
        y   = xyz[1] + uvw[1]*a_DT
        z   = xyz[2] + uvw[2]*a_DT

        #### AD-HOC REPLACE WITH BOUNDARY CONDITION FUNCTION (TO BE IMPLEMENTED)

        a_TracerPoints.SetPoint(p, x, y, z)

    return a_TracerPoints

#--------------------------------------------------------------------------
# time integration routine for individual particles: Runge Kutta Integrator
#--------------------------------------------------------------------------
def integrateRK4(a_TracerPoints, a_FlowObject, a_T, a_DT, a_BoundaryCondition):

    for p in range(a_TracerPoints.GetNumberOfPoints()):

        xyz = a_TracerPoints.GetPoint(p)
        k1  = a_FlowObject.eval(xyz[0],                 xyz[1],                 xyz[2],                 a_T)
        k2  = a_FlowObject.eval(xyz[0] + 0.5*a_DT*k1,   xyz[1] + 0.5*a_DT*k1,   xyz[2] + 0.5*a_DT*k1,   a_T + 0.5*a_DT)
        k3  = a_FlowObject.eval(xyz[0] + 0.5*a_DT*k2,   xyz[1] + 0.5*a_DT*k2,   xyz[2] + 0.5*a_DT*k2,   a_T + 0.5*a_DT)
        k4  = a_FlowObject.eval(xyz[0] + a_DT*k3,       xyz[1] + a_DT*k3,       xyz[2] + a_DT*k3,       a_T + a_DT)
        x   = xyz[0] + (a_DT/6.0)*(k1 + 2.0*k2 + 2.0*k3 + k4)
        y   = xyz[1] + (a_DT/6.0)*(k1 + 2.0*k2 + 2.0*k3 + k4)
        z   = xyz[2] + (a_DT/6.0)*(k1 + 2.0*k2 + 2.0*k3 + k4)

        #### AD-HOC REPLACE WITH BOUNDARY CONDITION FUNCTION (TO BE IMPLEMENTED)

        a_TracerPoints.SetPoint(p, x, y, z)

    return a_TracerPoints

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
    
    #---------------------------------------------------------
    # set-up analytical flow to be integrated from the library
    #---------------------------------------------------------
    flow = ABCFlow(a_A=np.sqrt(3.0), a_B=np.sqrt(2.0), a_C=1.0)

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

    #-------------------
    # start time counter
    #-------------------
    timeIndex   = 0
    simTime     = inputData.getSimulationStartTime()

    while simTime <= inputData.getSimulationStopTime():

        print "Integrating from", timeIndex, " to ", timeIndex + 1, "simTime", simTime
        
        #---------------------------------------------------------------------------
        # create a boolean condition that governs when new points are to be injected
        #---------------------------------------------------------------------------
        isInjectPoints = True

        #-----------------------------------------
        # inject tracers into the domain if needed
        #-----------------------------------------
        if isInjectPoints:            
            tracerPoints = injectPointsWithRandomPerturbations(tracerPoints, tracerInject)
            print "New Particles Injected"
            print "Now Integrating", tracerPoints.GetNumberOfPoints(), "Particles"
        
        boundaryCondition   = 1
        tracerPoints        = integrateForwardEuler(tracerPoints, flow, simTime, inputData.getIntegrationTimeStep(), boundaryCondition)
        
        #--------------------------------------------------------------
        # at the end of appropriate number of steps dump data into file
        #--------------------------------------------------------------
        writeTracerDataToFile(tracerPoints, inputData.getTracerOutputFile(a_ID=timeIndex))
        
        #----------------------------------------
        # update time indices and simulation time
        #----------------------------------------
        timeIndex   = timeIndex + 1
        simTime     = simTime + inputData.getIntegrationTimeStep()
