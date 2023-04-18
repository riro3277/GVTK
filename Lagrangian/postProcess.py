#-------------------------------------------------------------------------------------------------
# This is a collection of functions to perform post-processing computations for particle transport
# simulations based on an input file.
#
# Additional detailed features to be included in future versions.
#
# Author:       Zachariah Irwin
# Institution:  University of Colorado, Boulder
# Last Edit:    March 2020
#-------------------------------------------------------------------------------------------------
import sys, os

try:
  import numpy as np
except ImportError:
  print("Module Error! Numpy not installed.")
  sys.exit()

try:
  import vtk
except ImportError:
  print("Module Error! VTK module with Python bindings not installed.")
  sys.exit()

#---------------------------------------------------------------------
## @brief Utility function that converts Cartesian indices into a 
# linear index to map into data arrays
#
# @param[in]   a_Indices    i,j,k index
# @param[in]   a_Dims       Array number of points along x,y,z
#
#---------------------------------------------------------------------
def linearIndexing(a_Indices, a_Dims):

  if len(a_Indices) == 2:
    return a_Indices[0] + a_Dims[0] * a_Indices[1]
  elif len(a_Indices) == 3:
    return a_Indices[0] + a_Dims[0] * (a_Indices[1] + a_Dims[1] * a_Indices[2])

#---------------------------------------------------------------------------
## @brief Utility function to write Lagrangian numpy arrays into
# a vtkStructuredGrid file
#
# @param[in]  a_LagFlags        Flags denoting which fields to save to vtk
# @param[in]  a_VecFlags        Flags denoting which directions along which
#                               stretch/strain was computed
# @param[in]  a_RefCoordinates  Numpy array of grid coordinates in reference
#                               configuration
# @param[in]  a_FTLE_T          Numpy array of time-dependent FTLE field
# @param[in]  a_FTLE_NoT        Numpy array of time-independent FTLE field
# @param[in]  a_Stretch_1|2|3   Numpy array of stretch along directions 1|2|3
# @param[in]  a_Strain_1|2|3    Numpy array of strain along directions 1|2|3
# @param[in]  a_OutputFileName  Name of vtk output file
#----------------------------------------------------------------------------
def writeLagrangianFieldToVTK(a_RefCoordinates, a_RefNums, a_FTLE_T, a_FTLE_NoT, \
                              a_Stretch_1, a_Stretch_2, a_Stretch_3, \
                              a_Strain_1, a_Strain_2, a_Strain_3, a_LagFlags, a_VecFlags, a_OutputFile):
  
  xN = a_RefNums[0]
  yN = a_RefNums[1]
  zN = a_RefNums[2]
    
  vtkGrid = vtk.vtkStructuredGrid()
  vtkGrid.SetDimensions([xN, yN, zN])

  points  = vtk.vtkPoints()
  points.SetNumberOfPoints(xN*yN*zN)

  if a_LagFlags[0] == True:

    ftleArr_T = vtk.vtkDoubleArray()
    ftleArr_T.SetName('FTLE - T')
    ftleArr_T.SetNumberOfComponents(1)
    ftleArr_T.SetNumberOfTuples(xN*yN*zN)

    ftleArr_NoT = vtk.vtkDoubleArray()
    ftleArr_NoT.SetName('FTLE - No T')
    ftleArr_NoT.SetNumberOfComponents(1)
    ftleArr_NoT.SetNumberOfTuples(xN*yN*zN)

  if a_LagFlags[1] == True:

    if a_VecFlags[0] == True:

      stretchArr_1 = vtk.vtkDoubleArray()
      stretchArr_1.SetName('Stretch: 1')
      stretchArr_1.SetNumberOfComponents(1)
      stretchArr_1.SetNumberOfTuples(xN*yN*zN)

    if a_VecFlags[1] == True:

      stretchArr_2 = vtk.vtkDoubleArray()
      stretchArr_2.SetName('Stretch: 2')
      stretchArr_2.SetNumberOfComponents(1)
      stretchArr_2.SetNumberOfTuples(xN*yN*zN)

    if a_VecFlags[2] == True:

      stretchArr_3 = vtk.vtkDoubleArray()
      stretchArr_3.SetName('Stretch: 3')
      stretchArr_3.SetNumberOfComponents(1)
      stretchArr_3.SetNumberOfTuples(xN*yN*zN)

  if a_LagFlags[2] == True:

    if a_VecFlags[0] == True:

      strainArr_1 = vtk.vtkDoubleArray()
      strainArr_1.SetName('Strain: 1')
      strainArr_1.SetNumberOfComponents(1)
      strainArr_1.SetNumberOfTuples(xN*yN*zN)

    if a_VecFlags[1] == True:

      strainArr_2 = vtk.vtkDoubleArray()
      strainArr_2.SetName('Strain: 2')
      strainArr_2.SetNumberOfComponents(1)
      strainArr_2.SetNumberOfTuples(xN*yN*zN)

    if a_VecFlags[2] == True:

      strainArr_3 = vtk.vtkDoubleArray()
      strainArr_3.SetName('Strain: 3')
      strainArr_3.SetNumberOfComponents(1)
      strainArr_3.SetNumberOfTuples(xN*yN*zN)

  for p in range(a_RefCoordinates.shape[0]):
    xyz = a_RefCoordinates[p,:]
    points.InsertPoint(p, xyz)

    if a_LagFlags[0] == True:
      ftleArr_T.SetTuple1(p, a_FTLE_T[p])
      ftleArr_NoT.SetTuple1(p, a_FTLE_NoT[p])
    if a_LagFlags[1] == True:
      if a_VecFlags[0] == True:
        stretchArr_1.SetTuple1(p, a_Stretch_1[p])
      if a_VecFlags[1] == True:
        stretchArr_2.SetTuple1(p, a_Stretch_2[p])
      if a_VecFlags[2] == True:
        stretchArr_3.SetTuple1(p, a_Stretch_3[p])
    if a_LagFlags[2] == True:
      if a_VecFlags[0] == True:
        strainArr_1.SetTuple1(p, a_Strain_1[p])
      if a_VecFlags[1] == True:
        strainArr_2.SetTuple1(p, a_Strain_2[p])
      if a_VecFlags[2] == True:
        strainArr_3.SetTuple1(p, a_Strain_3[p])

  vtkGrid.SetPoints(points)
  if a_LagFlags[0] == True:
    vtkGrid.GetPointData().AddArray(ftleArr_T)
    vtkGrid.GetPointData().AddArray(ftleArr_NoT)
  if a_LagFlags[1] == True:
    if a_VecFlags[0] == True:
      vtkGrid.GetPointData().AddArray(stretchArr_1)
    if a_VecFlags[1] == True:
      vtkGrid.GetPointData().AddArray(stretchArr_2)
    if a_VecFlags[2] == True:
      vtkGrid.GetPointData().AddArray(stretchArr_3)
  if a_LagFlags[2] == True:
    if a_VecFlags[0] == True:
      vtkGrid.GetPointData().AddArray(strainArr_1)
    if a_VecFlags[1] == True:
      vtkGrid.GetPointData().AddArray(strainArr_2)
    if a_VecFlags[2] == True:
      vtkGrid.GetPointData().AddArray(strainArr_3)

  writer = vtk.vtkStructuredGridWriter()
  writer.SetFileName(a_OutputFile)
  writer.SetInputData(vtkGrid)
  writer.SetFileTypeToBinary()
  writer.Update()
  writer.Write()
