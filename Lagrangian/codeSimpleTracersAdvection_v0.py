#---------------------------------------------------------------------------------------------
# This is a modified version of the simpleTracers code, which is designed based on
# a simple code design philosophy that heavily makes use of the VTK library and numpy/scipy
# for all major geometry and compute operations.
#
# The objective of this module is to:
# - simply perform an advection computation for tracers released in a vector field.
# - develop a plugin for particle tracer computation in SimVascular software suite.
#
# Additional detailed features to be included in future versions
#
# All range() converted into xrange()
#
# Author:       Debanjan Mukherjee
# Institution:  University of California, Berkeley
# Last Edit:    June 2018
#---------------------------------------------------------------------------------------------

#-----------------------------------
# BEGIN MODULE IMPORTS
#-----------------------------------

import sys, os
import vtk
import numpy as np

from moduleTracerInput import *

#-----------------------------------
# END MODULE IMPORTS
#-----------------------------------

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

    for p in xrange(addNumPts):
        addID = oldNumPts + p
        a_Points.InsertPoint(addID, a_InputPoints.GetPoint(p))

    return a_Points

def injectPointsRandomly(a_Points, a_N, a_XRandomize=None, a_YRandomize=None, a_ZRandomize=None):

    oldNumPts   = a_Points.GetNumberOfPoints()
    addNumPts   = a_N

    a_Points.ComputeBounds()
    boundBox    = a_Points.GetBounds()
    xMin        = boundBox[0]
    xMax        = boundBox[1]
    yMin        = boundBox[2]
    yMax        = boundBox[3]
    zMin        = boundBox[4]
    zMax        = boundBox[5]

    if a_XRandomize is not None:
        x   = np.random.uniform(low=xMin, high=xMax, size=addNumPts)
        for p in xrange(addNumPts):
            addID   = oldNumPts + p
            a_Points.InsertPoint(addID, [x[p], yMin, zMin])

    if a_YRandomize is not None:
        y   = np.random.uniform(low=yMin, high=yMax, size=addNumPts)
        for p in xrange(addNumPts):
            addID   = oldNumPts + p
            a_Points.InsertPoint(addID, [xMin, y[p], zMin])

    if a_ZRandomize is not None:
        z   = np.random.uniform(low=zMin, high=zMax, size=addNumPts)
        for p in xrange(addNumPts):
            addID   = oldNumPts + p
            a_Points.InsertPoint(addID, [xMin, yMin, z[p]])

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
    for p in xrange(tracerPoints.GetNumberOfPoints()):
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
# -------
# Params:
# -------
# a_FileName:   a vtk unstructured grid file to read data from
# a_DataName:   (optional) name of specific data field that is to be read in
# --------
# Returns:
# --------
# data: a vtkDataArray object (if a_DataName specified), a vtkUnstructuredGrid object otherwise
#----------------------------------------------------------------------------------------------
def extractFlowDataFromFile(a_FileName, a_DataName=None, a_FoamBlockID=0):

    if a_FileName.endswith('vtk'):

        reader = vtk.vtkUnstructuredGridReader()
        reader.SetFileName(a_FileName)
        reader.Update()

        if a_DataName is not None:
            return reader.GetOutput().GetPointData().GetArray(a_DataName)
        else:
            return reader.GetOutput()

    elif a_FileName.endswith('vtu'):

        reader = vtk.vtkXMLUnstructuredGridReader()
        reader.SetFileName(a_FileName)
        reader.Update()

        if a_DataName is not None:
            return reader.GetOutput().GetPointData().GetArray(a_DataName)
        else:
            return reader.GetOutput()

    elif a_FileName.endswith('foam'):

        reader = vtk.vtkOpenFOAMReader()
        reader.SetFileName(a_FileName)
        reader.Update()
        foamData = reader.GetOutput()
        numBlock = foamData.GetNumberOfBlocks()
        blockData = foamData.GetBlock(a_FoamBlockID)

        if a_DataName is not None:
            return blockData.GetPointData().GetArray(a_DataName)
        else:
            return blockData

#-------------------------------------------------------------------
# Function to write tracer data back to file
# This new version can accept both vtkPoints and vtkPolyData objects
# -------
# Params:
# -------
# a_Tracers:        vtkPoints or vtkPolyData object
# a_OutputFileName: polydata filename to dump data into
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

#--------------------------------------------------------------------------
# Function to create a cell locator object from input mesh data from a file
# -------
# Params:
# -------
# a_FileName:       file containing the mesh where particles are to be located
# a_LocatorType:    (optional) type of vtk locator being employed
# --------
# Returns:
# --------
# locator:  a vtkLocator object (not sure what this object data entails)
#--------------------------------------------------------------------------
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

#------------------------------------------------------------------------
# time integration routine for individual particles: Forward Euler Method
#------------------------------------------------------------------------
def integrateForwardEuler(a_TracerPoints, a_Locator, a_Velocities, a_GridData, a_T, a_DT, a_Window, a_BoundaryCondition):

    for p in xrange(a_TracerPoints.GetNumberOfPoints()):

        ## IMPLEMENTATION WITHOUT NUMPY ARRAY CALLS
        velMinus    = np.empty(3, dtype=np.float32)
        velPlus     = np.empty(3, dtype=np.float32)
        vel         = np.empty(3, dtype=np.float32)

        #
        # step 1: locate
        #
        xyz         = a_TracerPoints.GetPoint(p)
        cell        = a_Locator.FindCell(xyz)
        cellPtIds   = vtk.vtkIdList()

        #
        # interpolate and integrate (combined since one-step scheme)
        #
        if len(a_Velocities) == 2:

            if cell != -1:
                a_GridData.GetCellPoints(cell, cellPtIds)

                ## ORIGINAL IMPLEMENTATION:
                #velMinus_N1 = np.asarray(a_Velocities[0].GetTuple3(cellPtIds.GetId(0)))
                #velMinus_N2 = np.asarray(a_Velocities[0].GetTuple3(cellPtIds.GetId(1)))
                #velMinus_N3 = np.asarray(a_Velocities[0].GetTuple3(cellPtIds.GetId(2)))
                #velMinus    = (1.0/3.0)*(velMinus_N1 + velMinus_N2 + velMinus_N3)

                #velPlus_N1  = np.asarray(a_Velocities[1].GetTuple3(cellPtIds.GetId(0)))
                #velPlus_N2  = np.asarray(a_Velocities[1].GetTuple3(cellPtIds.GetId(1)))
                #velPlus_N3  = np.asarray(a_Velocities[1].GetTuple3(cellPtIds.GetId(2)))
                #velPlus     = (1.0/3.0)*(velPlus_N1 + velPlus_N2 + velPlus_N3)

                ## IMPLEMENTATION WITHOUT NUMPY ARRAY CALLS:
                velMinus_N1 = a_Velocities[0].GetTuple3(cellPtIds.GetId(0))
                velMinus_N2 = a_Velocities[0].GetTuple3(cellPtIds.GetId(1))
                velMinus_N3 = a_Velocities[0].GetTuple3(cellPtIds.GetId(2))
                velMinus[0] = (1.0/3.0)*(velMinus_N1[0] + velMinus_N2[0] + velMinus_N3[0])
                velMinus[1] = (1.0/3.0)*(velMinus_N1[1] + velMinus_N2[1] + velMinus_N3[1])
                velMinus[2] = (1.0/3.0)*(velMinus_N1[2] + velMinus_N2[2] + velMinus_N3[2])

            else:

                ## ORIGINAL IMPLEMENTATION:
                #velMinus    = np.array([0.0,0.0,0.0])
                #velPlus     = np.array([0.0,0.0,0.0])

                ## IMPLEMENTATION WITHOUT NUMPY ARRAY CALLS:
                velMinus[0] = 0.0
                velMinus[1] = 0.0
                velMinus[2] = 0.0
                velPlus[0]  = 0.0
                velPlus[1]  = 0.0
                velPlus[2]  = 0.0

            velInterp   = velMinus + (a_T - a_Window[0])*(velPlus - velMinus)/(a_Window[1] - a_Window[0])
            xyz         = xyz + velInterp*a_DT

        elif len(a_Velocities) == 1:

            if cell != -1:
                a_GridData.GetCellPoints(cell, cellPtIds)

                ## ORIGINAL IMPLEMENTATION:
                #vel_N1 = np.asarray(a_Velocities[0].GetTuple3(cellPtIds.GetId(0)))
                #vel_N2 = np.asarray(a_Velocities[0].GetTuple3(cellPtIds.GetId(1)))
                #vel_N3 = np.asarray(a_Velocities[0].GetTuple3(cellPtIds.GetId(2)))
                #vel    = (1.0/3.0)*(vel_N1 + vel_N2 + vel_N3)

                ## IMPLEMENTATION WITHOUT NUMPY ARRAY CALLS:
                vel_N1 = a_Velocities[0].GetTuple3(cellPtIds.GetId(0))
                vel_N2 = a_Velocities[0].GetTuple3(cellPtIds.GetId(1))
                vel_N3 = a_Velocities[0].GetTuple3(cellPtIds.GetId(2))
                vel[0]  = (1.0/3.0)*(vel_N1[0] + vel_N2[0] + vel_N3[0])
                vel[1]  = (1.0/3.0)*(vel_N1[1] + vel_N2[1] + vel_N3[1])
                vel[2]  = (1.0/3.0)*(vel_N1[2] + vel_N2[2] + vel_N3[2])

            else:

                ## ORIGINAL IMPLEMENTATION:
                #vel     = np.array([0.0,0.0,0.0])

                ## IMPLEMENTATION WITHOUT NUMPY ARRAY CALLS:
                vel[0]  = 0.0
                vel[1]  = 0.0
                vel[2]  = 0.0

            xyz = xyz + vel*a_DT

        #### AD-HOC REPLACE WITH BOUNDARY CONDITION FUNCTION

        #if xyz[0] > 0.060: xyz[0] = 0.061

        #if xyz[0] < -0.060: xyz[0] = -0.061

        #if xyz[1] > 0.015: xyz[1] = 0.016

        #if xyz[1] < -0.015: xyz[1] = -0.016

        if xyz[0] > 45.0: xyz[0] = 45.01 # This implementation does not respect walls vs outlets

        if xyz[0] < 0.0: xyz[0] = -0.01

        if xyz[1] > 6.0: xyz[1] = 6.01

        if xyz[1] < 0.0: xyz[1] = -0.01

        a_TracerPoints.SetPoint(p, xyz[0], xyz[1], xyz[2])

    return a_TracerPoints

#--------------------------------------------------------------------------
# time integration routine for individual particles: Runge Kutta Integrator
#
# implementation of this integrator follows the standard operation loop:
# locate -> interpolate in space -> interpolate in time -> integrate
#
# k1 = v(x0, t0)
# k2 = v(x0 + 0.5*dt*k1, t0+0.5*dt)
# k3 = v(x0 + 0.5*dt*k2, t0+0.5*dt)
# k4 = v(x0 + dt*k3,     t0 + dt)
# x  = x0 + (dt/6)*(k1 + 2*k2 + 2*k3 + k4)
# t* = t0 + 0.5*dt and t** = t0+dt should both be in the time window loaded
#
# Thus the only time when a_Velocities will have one entry is when the flow
# field being integrated is steady
#--------------------------------------------------------------------------
def integrateRK4(a_TracerPoints, a_Locator, a_Velocities, a_GridData, a_T, a_DT, a_Window, a_BoundaryCondition):

    checkTimeOut  = (((a_T + 0.5*a_DT) > a_Window[1]) or ((a_T + a_DT) > a_Window[1])) and len(a_Window) == 2

    if checkTimeOut:
        sys.exit("RK Integrators need all time interpolation points to be adjusted")

    for p in xrange(a_TracerPoints.GetNumberOfPoints()):

        xyz         = np.asarray(a_TracerPoints.GetPoint(p))
        cell        = a_Locator.FindCell(xyz)
        cellPtIds   = vtk.vtkIdList()

        #
        # stage 1 update: location-step
        #
        cell    = a_Locator.FindCell(xyz)

        #
        # stage 1 update: interpolation-step
        #
        if len(a_Velocities) == 1:
            #
            # steady flow field: only one velocity field loaded always
            #
            if cell != -1:

                a_GridData.GetCellPoints(cell, cellPtIds)
                vel_N1  = np.asarray(a_Velocities[0].GetTuple3(cellPtIds.GetId(0)))
                vel_N2  = np.asarray(a_velocities[0].GetTuple3(cellPtIds.GetId(1)))
                vel_N3  = np.asarray(a_Velocities[0].GetTuple3(cellPtIds.GetId(2)))
                vel     = (1.0/3.0)*(velMinus_N1 + velMinus_N2 + velMinus_N3)

            else:

                vel = np.zeros(3)

        else:
            #
            # unsteady flow field: interpolate in time as per loaded data frames
            #
            if cell != -1:

                a_GridData.GetCellPoints(cell, cellPtIds)

                velMinus_N1 = np.asarray(a_Velocities[0].GetTuple3(cellPtIds.GetId(0)))
                velMinus_N2 = np.asarray(a_Velocities[0].GetTuple3(cellPtIds.GetId(1)))
                velMinus_N3 = np.asarray(a_Velocities[0].GetTuple3(cellPtIds.GetId(2)))

                velPlus_N1  = np.asarray(a_Velocities[1].GetTuple3(cellPtIds.GetId(0)))
                velPlus_N2  = np.asarray(a_Velocities[1].GetTuple3(cellPtIds.GetId(1)))
                velPlus_N3  = np.asarray(a_Velocities[1].GetTuple3(cellPtIds.GetId(2)))

                velMinus    = (1.0/3.0)*(velMinus_N1 + velMinus_N2 + velMinus_N3)
                velPlus     = (1.0/3.0)*(velPlus_N1  + velPlus_N2  + velPlus_N3)

            else:

                velMinus    = np.zeros(3)
                velPlus     = np.zeros(3)

            vel = velMinus + (a_T - a_Window[0])*(velPlus - velMinus)/(a_Window[1] - a_Window[0])

        #
        # stage 1 update: integration-step
        #
        k1  = vel

        #
        # stage 2 update: location-step
        #
        cell    = a_Locator.FindCell(xyz + 0.5*k1*a_DT)

        #
        # stage 2 update: interpolation-step
        #
        if len(a_Velocities) == 1:
            #
            # steady flow field: only one velocity field loaded always
            #
            if cell!= -1:

                a_GridData.GetCellPoints(cell, cellPtIds)
                vel_N1  = np.asarray(a_Velocities[0].GetTuple3(cellPtIds.GetId(0)))
                vel_N2  = np.asarray(a_velocities[0].GetTuple3(cellPtIds.GetId(1)))
                vel_N3  = np.asarray(a_Velocities[0].GetTuple3(cellPtIds.GetId(2)))
                vel     = (1.0/3.0)*(velMinus_N1 + velMinus_N2 + velMinus_N3)

            else:

                vel = np.zeros(3)

        else:
            #
            # unsteady flow field: interpolate in time as per loaded data frames
            #
            if cell != -1:

                a_GridData.GetCellPoints(cell, cellPtIds)

                if (a_T + 0.5*a_DT >= a_window[0]) and (a_T + 0.5*a_DT <= a_Window[1]):

                    velMinus_N1 = np.asarray(a_Velocities[0].GetTuple3(cellPtIds.GetId(0)))
                    velMinus_N2 = np.asarray(a_Velocities[0].GetTuple3(cellPtIds.GetId(1)))
                    velMinus_N3 = np.asarray(a_Velocities[0].GetTuple3(cellPtIds.GetId(2)))

                    velPlus_N1  = np.asarray(a_Velocities[1].GetTuple3(cellPtIds.GetId(0)))
                    velPlus_N2  = np.asarray(a_Velocities[1].GetTuple3(cellPtIds.GetId(1)))
                    velPlus_N3  = np.asarray(a_Velocities[1].GetTuple3(cellPtIds.GetId(2)))

                    velMinus    = (1.0/3.0)*(velMinus_N1 + velMinus_N2 + velMinus_N3)
                    velPlus     = (1.0/3.0)*(velPlus_N1  + velPlus_N2  + velPlus_N3)

                    vel         = velMinus + (a_T - a_Window[0])*(velPlus - velMinus)/(a_Window[1] - a_Window[0])

                elif (a_T + 0.5*a_DT > a_Window[1]):

                    velMinus_N1 = np.asarray(a_Velocities[1].GetTuple3(cellPtIds.GetId(0)))
                    velMinus_N2 = np.asarray(a_Velocities[1].GetTuple3(cellPtIds.GetId(1)))
                    velMinus_N3 = np.asarray(a_Velocities[1].GetTuple3(cellPtIds.GetId(2)))

                    velPlus_N1  = np.asarray(a_Velocities[2].GetTuple3(cellPtIds.GetId(0)))
                    velPlus_N2  = np.asarray(a_Velocities[2].GetTuple3(cellPtIds.GetId(1)))
                    velPlus_N3  = np.asarray(a_Velocities[2].GetTuple3(cellPtIds.GetId(2)))

                    velMinus    = (1.0/3.0)*(velMinus_N1 + velMinus_N2 + velMinus_N3)
                    velPlus     = (1.0/3.0)*(velPlus_N1  + velPlus_N2  + velPlus_N3)

                    vel         = velMinus + (a_T - a_Window[1])*(velPlus - velMinus)/(a_Window[2] - a_Window[1])

            else:

                velMinus    = np.zeros(3)
                velPlus     = np.zeros(3)
                vel         = np.zeros(3)

        #
        # stage 2 update: integration-step
        #
        k2  = vel

        #
        # stage 3 update: location-step
        #
        cell = a_Locator.FindCell(xyz + 0.5*k2*a_DT)

        #
        # stage 3 update: interpolation-step
        #
        if len(a_Velocities) == 1:
            #
            # steady flow field: only one velocity field loaded always
            #
            if cell != -1:

                a_GridData.GetCellPoints(cell, cellPtIds)
                vel_N1  = np.asarray(a_Velocities[0].GetTuple3(cellPtIds.GetId(0)))
                vel_N2  = np.asarray(a_velocities[0].GetTuple3(cellPtIds.GetId(1)))
                vel_N3  = np.asarray(a_Velocities[0].GetTuple3(cellPtIds.GetId(2)))
                vel     = (1.0/3.0)*(velMinus_N1 + velMinus_N2 + velMinus_N3)

            else:

                vel = np.zeros(3)

        else:
            #
            # unsteady flow field: interpolate in time as per loaded data frames
            #
            if cell != -1:

                a_GridData.GetCellPoints(cell, cellPtIds)

                if (a_T + 0.5*a_DT >= a_window[0]) and (a_T + 0.5*a_DT <= a_Window[1]):

                    velMinus_N1 = np.asarray(a_Velocities[0].GetTuple3(cellPtIds.GetId(0)))
                    velMinus_N2 = np.asarray(a_Velocities[0].GetTuple3(cellPtIds.GetId(1)))
                    velMinus_N3 = np.asarray(a_Velocities[0].GetTuple3(cellPtIds.GetId(2)))

                    velPlus_N1  = np.asarray(a_Velocities[1].GetTuple3(cellPtIds.GetId(0)))
                    velPlus_N2  = np.asarray(a_Velocities[1].GetTuple3(cellPtIds.GetId(1)))
                    velPlus_N3  = np.asarray(a_Velocities[1].GetTuple3(cellPtIds.GetId(2)))

                    velMinus    = (1.0/3.0)*(velMinus_N1 + velMinus_N2 + velMinus_N3)
                    velPlus     = (1.0/3.0)*(velPlus_N1  + velPlus_N2  + velPlus_N3)

                    vel         = velMinus + (a_T - a_Window[0])*(velPlus - velMinus)/(a_Window[1] - a_Window[0])

                elif (a_T + 0.5*a_DT > a_Window[1]):

                    velMinus_N1 = np.asarray(a_Velocities[1].GetTuple3(cellPtIds.GetId(0)))
                    velMinus_N2 = np.asarray(a_Velocities[1].GetTuple3(cellPtIds.GetId(1)))
                    velMinus_N3 = np.asarray(a_Velocities[1].GetTuple3(cellPtIds.GetId(2)))

                    velPlus_N1  = np.asarray(a_Velocities[2].GetTuple3(cellPtIds.GetId(0)))
                    velPlus_N2  = np.asarray(a_Velocities[2].GetTuple3(cellPtIds.GetId(1)))
                    velPlus_N3  = np.asarray(a_Velocities[2].GetTuple3(cellPtIds.GetId(2)))

                    velMinus    = (1.0/3.0)*(velMinus_N1 + velMinus_N2 + velMinus_N3)
                    velPlus     = (1.0/3.0)*(velPlus_N1  + velPlus_N2  + velPlus_N3)

                    vel         = velMinus + (a_T - a_Window[1])*(velPlus - velMinus)/(a_Window[2] - a_Window[1])

            else:

                velMinus    = np.zeros(3)
                velPlus     = np.zeros(3)
                vel         = np.zeros(3)

        #
        # stage 3 update: integration-step
        #
        k3 = vel

        #
        # stage 4 update: location-step
        #
        cell = a_Locator.FindCell(xyz + k3*a_DT)

        #
        # stage 4 update: interpolation-step
        #
        if len(a_Velocities) == 1:
            #
            # steady flow field: only one velocity field loaded always
            #
            if cell != -1:

                a_GridData.GetCellPoints(cell, cellPtIds)
                vel_N1  = np.asarray(a_Velocities[0].GetTuple3(cellPtIds.GetId(0)))
                vel_N2  = np.asarray(a_velocities[0].GetTuple3(cellPtIds.GetId(1)))
                vel_N3  = np.asarray(a_Velocities[0].GetTuple3(cellPtIds.GetId(2)))
                vel     = (1.0/3.0)*(velMinus_N1 + velMinus_N2 + velMinus_N3)

            else:

                vel = np.zeros(3)
        else:
            #
            # unsteady flow field: interpolate in time as per loaded data frames
            #
            if cell != -1:

                a_GridData.GetCellPoints(cell, cellPtIds)

                if (a_T + a_DT >= a_window[0]) and (a_T + a_DT <= a_Window[1]):

                    velMinus_N1 = np.asarray(a_Velocities[0].GetTuple3(cellPtIds.GetId(0)))
                    velMinus_N2 = np.asarray(a_Velocities[0].GetTuple3(cellPtIds.GetId(1)))
                    velMinus_N3 = np.asarray(a_Velocities[0].GetTuple3(cellPtIds.GetId(2)))

                    velPlus_N1  = np.asarray(a_Velocities[1].GetTuple3(cellPtIds.GetId(0)))
                    velPlus_N2  = np.asarray(a_Velocities[1].GetTuple3(cellPtIds.GetId(1)))
                    velPlus_N3  = np.asarray(a_Velocities[1].GetTuple3(cellPtIds.GetId(2)))

                    velMinus    = (1.0/3.0)*(velMinus_N1 + velMinus_N2 + velMinus_N3)
                    velPlus     = (1.0/3.0)*(velPlus_N1  + velPlus_N2  + velPlus_N3)

                    vel         = velMinus + (a_T - a_Window[0])*(velPlus - velMinus)/(a_Window[1] - a_Window[0])

                elif (a_T + a_DT > a_Window[1]):

                    velMinus_N1 = np.asarray(a_Velocities[1].GetTuple3(cellPtIds.GetId(0)))
                    velMinus_N2 = np.asarray(a_Velocities[1].GetTuple3(cellPtIds.GetId(1)))
                    velMinus_N3 = np.asarray(a_Velocities[1].GetTuple3(cellPtIds.GetId(2)))

                    velPlus_N1  = np.asarray(a_Velocities[2].GetTuple3(cellPtIds.GetId(0)))
                    velPlus_N2  = np.asarray(a_Velocities[2].GetTuple3(cellPtIds.GetId(1)))
                    velPlus_N3  = np.asarray(a_Velocities[2].GetTuple3(cellPtIds.GetId(2)))

                    velMinus    = (1.0/3.0)*(velMinus_N1 + velMinus_N2 + velMinus_N3)
                    velPlus     = (1.0/3.0)*(velPlus_N1  + velPlus_N2  + velPlus_N3)

                    vel         = velMinus + (a_T - a_Window[1])*(velPlus - velMinus)/(a_Window[2] - a_Window[1])

            else:

                velMinus    = np.zeros(3)
                velPlus     = np.zeros(3)
                vel         = np.zeros(3)

        #
        # stage 4 update: integration-step
        #
        k4  = vel
        xyz = xyz + (a_DT/6.0)*(k1 + 2.0*k2 + 2.0*k3 + k4)

        #### AD-HOC REPLACE WITH BOUNDARY CONDITION FUNCTION (TO BE IMPLEMENTED)

        a_TracerPoints.SetPoint(p, xyz[0], xyz[1], xyz[2])

    return a_TracerPoints

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

    isDebug = True

    #--------------------------
    # set-up problem input data
    #--------------------------
    inputData = SimInputs(inputFile)
    inputData.readInputFile()

    #----------------------------------------------------------
    # initialize tracer data from an external VTK polydata file
    #----------------------------------------------------------
    tracerPoints    = initializeTracerFromFile(inputData.getTracerInput())
    tracerInject    = initializeTracerFromFile(inputData.getTracerInput())
    numInit         = tracerPoints.GetNumberOfPoints()

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

        print("Integrating {0:8d} tracers from {1:6d} to {2:6d} sim-time: {3:6f}".\
                format(tracerPoints.GetNumberOfPoints(),timeIndex,timeIndex + 1,simTime))

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
        #isInjectPoints = isLoadFrame    ### THIS NEEDS MORE POLISHING
        isInjectPoints  = simTime == inputData.getSimulationStartTime()

        #-------------------------------------------
        # a set of debug messages for program status
        #-------------------------------------------
        if isDebug:
            if isLoadFrame:
                print("Will Load Velocity Data")

            if isDataPointSync:
                if inputData.isSteadyFlowData():
                    if simTime == inputData.getSimulationStartTime():
                        print("Data And Integration Times Are Synced")
                else:
                    print("Data And Integration Times Are Synced")

            if isInjectPoints:
                print("New Particles Injected")
                print("Now Integrating", tracerPoints.GetNumberOfPoints(), "Particles")

        #-----------------------------------------
        # inject tracers into the domain if needed
        #-----------------------------------------
        if isInjectPoints:
            if isDebug: print("injecting here")

            tracerPoints = injectPoints(tracerPoints, tracerInject) ### THIS NEEDS MORE WORK
            #tracerPoints    = injectPointsRandomly(tracerPoints, numInit, a_YRandomize=1)

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

#-----------------------------------
# END COMPUTE SCRIPT DEFINITION
#-----------------------------------
