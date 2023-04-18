#---------------------------------------------------------------------------------------------
# This is a single file implementation of a script to compute tracer pathlines from analytical
# flow field data, which is designed based on a simple code design philosophy that heavily 
# makes use of the VTK library and numpy/scipy for all major geometry and compute operations. 
# 
# The objective of this module is to:
# - perform an FTLE computation for analytical flow data/fields
# - provide a single script for future development and testing of FTLE related capabilities.
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
from analyticalFlowLibrary import *

#---------------------------------------------------------------------
# utility function that converts Cartesian indices into a linear index 
# to map into data arrays
#---------------------------------------------------------------------
def linearIndexing(a_Indices, a_Dims):

    if len(a_Indices) == 2:
        return a_Indices[0] + a_Dims[1]*a_Indices[1]
    elif len(a_Indices) == 3:
        return a_Indices[0] + a_Dims[0]*(a_Indices[1] + a_Dims[1]*a_Indices[2])

#---------------------------------------------------------------------------------
# utility function that convertsa a bounding box with a specified number of points
# along each axis into a vtkStructuredGridObject for performing FTLE calculations
#---------------------------------------------------------------------------------
def convertBoxToStructuredGrid(a_BBox, a_Dims):

    #outGrid = vtk.vtkStructuredGrid()
    #outGrid.SetDimensions(a_Dims)

    outGrid = vtk.vtkPolyData()

    points  = vtk.vtkPoints()
    points.SetNumberOfPoints(a_Dims[0]*a_Dims[1]*a_Dims[2])

    deltaX  = (a_BBox[0,1] - a_BBox[0,0])/float(a_Dims[0] - 1)
    deltaY  = (a_BBox[1,1] - a_BBox[1,0])/float(a_Dims[1] - 1)
    deltaZ  = (a_BBox[2,1] - a_BBox[2,0])/float(a_Dims[2] - 1)

    for k in xrange(a_Dims[2]):
        for j in xrange(a_Dims[1]):
            for i in xrange(a_Dims[0]):
                offset  = linearIndexing([i,j,k], a_Dims)
                x       = a_BBox[0,0] + i*deltaX
                y       = a_BBox[1,0] + j*deltaY
                z       = a_BBox[2,0] + k*deltaZ
                points.InsertPoint(offset, [x,y,z])

    outGrid.SetPoints(points)
    return outGrid

#--------------------------------------------------------------
# read cartesian grid for seeding tracers for FTLE calculations
#--------------------------------------------------------------
def readCartesianTracers(a_StructuredGridFile):

    if a_StructuredGridFile.endswith('vts'):
        reader = vtk.vtkXMLStructuredGridReader()
    elif a_StructuredGridFile.endswith('vtk'):
        reader = vtk.vtkStructuredGridReader()

    reader.SetFileName(a_StructuredGridFile)
    reader.Update()

    return reader.GetOutput()

#----------------------------------------------------------------------------------
# function that reads raw cartesian grid data, and performs a set of pre-processing
# operations on it for enabling flow-map and ftle computations
#----------------------------------------------------------------------------------
def generateParticleID(a_StructuredGrid):

    getID   = vtk.vtkIdFilter()
    getID.SetInputData(a_StructuredGrid)
    getID.PointIdsOn()
    getID.CellIdsOn()
    getID.SetIdsArrayName('pId')
    getID.Update()

    return getID.GetOutput()

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

#-----------------------------------------------------
# Function to write cartesian tracer data back to file
#-----------------------------------------------------
def writeCartesianDataToFile(a_StructuredGrid, a_OutputFileName):

    if a_OutputFileName.endswith('vts'):
        dataWriter = vtk.vtkXMLStructuredGridWriter()
    elif a_OutputFileName.endswith('vtk'):
        dataWriter = vtk.vtkStructuredGridWriter()

    dataWriter.SetFileName(a_OutputFileName)

    if vtk.VTK_MAJOR_VERSION <= 5.0:
        dataWriter.SetInput(a_StructuredGrid)
    else:
        dataWriter.SetInputData(a_StructuredGrid)

    dataWriter.Update()
    dataWriter.Write()

def appendToPath():
    print "Functionality to be implemented. Presently embedded in main code"

def writePathlineToPolyData(a_XPath, a_YPath, a_ZPath, a_PathlineFile):

    if a_XPath.ndim == 1:
        
        points = vtk.vtkPoints()
        points.SetNumberOfPoints(a_XPath.size)
    
        line = vtk.vtkPolyLine()
        line.GetPointIds().SetNumberOfIds(a_XPath.size)

        for p in xrange(a_XPath.size):
            points.InsertPoint(p, (a_XPath[p], a_YPath[p], a_ZPath[p]))
            line.GetPointIds().SetId(p,p)

        cells   = vtk.vtkCellArray()
        cells.InsertNextCell(line)

        polyData = vtk.vtkPolyData()
        polyData.SetPoints(points)
        polyData.SetLines(cells)

        if a_PathlineFile.endswith('vtk'):
            writer = vtk.vtkPolyDataWriter()
        elif a_PathlineFile.endswith('vtp'):
            writer = vtk.vtkXMLPolyDataWriter()

        writer.SetFileName(a_PathlineFile)

        if vtk.VTK_MAJOR_VERSION <= 5.0:
            writer.SetInput(polyData)
        else:
            writer.SetInputData(polyData)

        writer.Update()
        writer.Write()

    elif a_XPath.ndim == 2:

        points = vtk.vtkPoints()
        points.SetNumberOfPoints(a_XPath.shape[0]*a_XPath.shape[1])

        cells   = vtk.vtkCellArray()
        
        pCount  = 0
        
        for p in xrange(a_XPath.shape[0]):
            
            line = vtk.vtkPolyLine()
            line.GetPointIds().SetNumberOfIds(a_XPath.shape[1])

            for t in xrange(a_XPath.shape[1]):
                points.InsertPoint(pCount, (a_XPath[p,t], a_YPath[p,t], a_ZPath[p,t]))
                line.GetPointIds().SetId(t, pCount)
                pCount = pCount + 1

            cells.InsertNextCell(line)
        
        polyData = vtk.vtkPolyData()
        polyData.SetPoints(points)
        polyData.SetLines(cells)

        if a_PathlineFile.endswith('vtk'):
            writer = vtk.vtkPolyDataWriter()
        elif a_PathlineFile.endswith('vtp'):
            writer = vtk.vtkXMLPolyDataWriter()

        writer.SetFileName(a_PathlineFile)

        if vtk.VTK_MAJOR_VERSION <= 5.0:
            writer.SetInput(polyData)
        else:
            writer.SetInputData(polyData)

        writer.Update()
        writer.Write()

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
        
        if a_TracerPoints.IsA('vtkPolyData'):
            a_TracerPoints.GetPoints().SetPoint(p, x, y, z)

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

    #-----------------------------------------------------------------------
    # set-up cartesian grid based on problem inputs for seeding tracers
    #-----------------------------------------------------------------------
    box     = np.array([[0.0, 2.0*np.pi],[0.0, 2.0*np.pi],[0.0, 2.0*np.pi]])
    dims    = np.array([10,10,10])
    tracers = convertBoxToStructuredGrid(box, dims) 

    #--------------------------------------------------------------
    # configure the tracers on this cartesian grid for computations
    #--------------------------------------------------------------
    tracers = generateParticleID(tracers)

    #-------------------
    # start time counter
    #-------------------
    timeIndex   = 0
    startTime   = inputData.getSimulationStartTime()
    simTime     = startTime
    stopTime    = inputData.getSimulationStopTime()

    #-------------------------------------------------------------------------------
    # initialize lists for holding pathline coordinates.
    # this can be improved, for now we use list.append() efficiency over np.append()
    # NOTE: THIS ONLY WORKS FOR FIXED NUMBER OF INPUT PARTICLES
    #-------------------------------------------------------------------------------
    numRows = int(tracers.GetNumberOfPoints()/inputData.getPathlinePointSubsampleInterval())
    numCols = int((stopTime-startTime)/inputData.getIntegrationTimeStep())
    print numRows, numCols
    xPath   = np.zeros((numRows, numCols), dtype=np.float32)
    yPath   = np.zeros((numRows, numCols), dtype=np.float32)
    zPath   = np.zeros((numRows, numCols), dtype=np.float32)
    nPathT  = 0

    #--------------------------------------------
    # integrate trajectories and compute flow map
    #--------------------------------------------
    while simTime <= stopTime:

        print "Integrating from {0:10d} to {1:10d} simTime {2:10f}".format(timeIndex, timeIndex + 1, simTime)
        
        #----------------------------------------------------------------------
        # perform time integration update based on choice of integration scheme
        #----------------------------------------------------------------------
        boundaryCondition   = 1
        if inputData.getIntegrationScheme() == 'feu':
            tracers = integrateForwardEuler(tracers, flow, simTime, inputData.getIntegrationTimeStep(), boundaryCondition)
        
        #---------------------------------------------------------------
        # append points to pathlines based on input data on sub-sampling 
        #---------------------------------------------------------------
        if np.mod(timeIndex, inputData.getPathlineTimeSubsampleInterval()) == 0:
            
            nPathP  = 0
            
            for p in xrange(tracers.GetNumberOfPoints()):
                if np.mod(tracers.GetPointData().GetArray('pId').GetTuple1(p), inputData.getPathlinePointSubsampleInterval()) == 0:
                    xyz                     = tracers.GetPoint(p)
                    xPath[nPathP, nPathT]   = xyz[0]
                    yPath[nPathP, nPathT]   = xyz[1]
                    zPath[nPathP, nPathT]   = xyz[2]
                    nPathP                  = nPathP + 1
            
            nPathT = nPathT + 1

        #----------------------------------------
        # update time indices and simulation time
        #----------------------------------------
        timeIndex   = timeIndex + 1
        simTime     = simTime + inputData.getIntegrationTimeStep()

    #-----------------------------------------------------------------
    # at the end of total flow map integration, compute the FTLE field
    #-----------------------------------------------------------------
    print "Writing pathlines into a polydata file"
    writePathlineToPolyData(xPath, yPath, zPath, inputData.getPathlineFile())

    #for r in xrange(numRows):
    #    writePathlineToPolyData(xPath[r,:], yPath[r,:], zPath[r,:], inputData.getPathlineFile(a_ID=r))        
