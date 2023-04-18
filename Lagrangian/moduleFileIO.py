#---------------------------------------------------------------------------------------------
# This is a module for handling all file IO operations for the library
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
# NEW VERSION
#----------------------------------------------------------------------------------------
def initializeTracerFromFile(a_FileName):

    if a_FileName.endswith('vtk'):
        reader = vtk.vtkPolyDataReader()
    elif a_FileName.endswith('vtp'):
        reader = vtk.vtkXMLPolyDataReader()

    reader.SetFileName(a_FileName)
    reader.Update()

    tracers = reader.GetOutput()

    return tracers

#---------------------------------------------------------------------------------------------
# This is a modified version of tracer initialization routine, where only the tracer injection
# coordinates are read in from the file, but all other data arrays are configured separately.
#
# This routine necessarily outputs only polydata format tracers
#
# -------
# Params:
# -------
# a_FileName            : polydata (vtp/vtk) filename to read in data from
# a_StartGlobalIdFrom   : the index from which global id for points need to be recomputed
#                         leave this to None if you want Id's to be computed from 0
# a_VelocityZero        : if set to True, this will initialize the velocity of the injected
#                         particles to be all 0.0 in each component
# a_VelocityFromFile    : if set to True, this will initialize the velocity of the injected
#                         particles to be all from the polydata file that has been read
# a_VelocityRandom      : if set to True, this will initialize the velocity of the injected
#                         to be all between the upper and lower bounds supplied here as tuple
# a_ScalarNamesZero     : the list of all scalar array names which will be initialized to 0.0
#---------------------------------------------------------------------------------------------
def placeTracersFromFile(a_FileName, a_StartGlobalIdFrom=None,
        a_VelocityZero=False, a_VelocityFromFile=True, a_VelocityRandom=False, a_ArrayNamesZero=None):

    if a_FileName.endswith('vtk'):
        reader = vtk.vtkPolyDataReader()
    elif a_FileName.endswith('vtp'):
        reader = vtk.vtkXMLPolyDataReader()

    reader.SetFileName(a_FileName)
    reader.Update()

    tracerPolyData = reader.GetOutput()

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
def extractGridDataFromFile(a_FileName, a_FoamBlockID=0):

    if a_FileName.endswith('vtk'):

        reader = vtk.vtkUnstructuredGridReader()
        reader.SetFileName(a_FileName)
        reader.Update()

        return reader.GetOutput()

    elif a_FileName.endswith('vtu'):

        reader = vtk.vtkXMLUnstructuredGridReader()
        reader.SetFileName(a_FileName)
        reader.Update()

        return reader.GetOutput()

    elif a_FileName.endswith('foam'):

        reader = vtk.vtkOpenFOAMReader()
        reader.SetFileName(a_FileName)
        reader.SkipZeroTimeOn()
        reader.Update()
        foamData = reader.GetOutput()
        numBlock = foamData.GetNumberOfBlocks()
        blockData = foamData.GetBlock(a_FoamBlockID)

        return blockData

#----------------------------------------------------------------------------------------------
# Function to extract flow data from a specified filename, returns a vtkUnstructuredGrid object
#----------------------------------------------------------------------------------------------
def extractDataArrayFromFile(a_FileName, a_DataName=None, a_FoamBlockID=0):

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
        reader.SkipZeroTimeOn()
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

        tracerWriter.SetFileTypeToBinary()
        tracerWriter.Update()
        tracerWriter.Write()

    elif a_Tracers.IsA('vtkPolyData') == 1:

        if vtk.VTK_MAJOR_VERSION <= 5.0:
            tracerWriter.SetInput(a_Tracers)
        else:
            tracerWriter.SetInputData(a_Tracers)

        tracerWriter.SetFileTypeToBinary()
        tracerWriter.Update()
        tracerWriter.Write()
