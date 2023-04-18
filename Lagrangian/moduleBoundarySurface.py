#------------------------------------------------------------------------------
# This module packages functionality for performing a naive check for particles
# or Lagrangian tracers interacting with a complex surface mesh triangulation
#
# The python based naive check is improved when f2py compiled modules are used 
#
# Author:       Debanjan Mukherjee
# Institution:  University of California, Berkeley
# Last Edited:  June 2018
#------------------------------------------------------------------------------
import sys, os
import vtk
import numpy as np

def readBoundarySurface(a_BoundaryFile):
    
    if a_BoundaryFile.endswith('vtp'):
        reader = vtk.vtkXMLPolyDataReader()
    elif a_BoundaryFile.endswith('vtk'):
        reader = vtk.vtkPolyDataReader()

    reader.SetFileName(a_BoundaryFile)
    reader.Update()

    return reader.GetOutput()

def createSurfaceNormals(a_SurfaceMesh):

    normalFilter = vtk.vtkPolyDataNormals()
    normalFilter.SetInputData(a_SurfaceMesh)
    normalFilter.ComputeCellNormalsOn()
    normalFiler.Update()

    return normalFilter.GetOutput()

def checkParticleSurfaceContact(a_SurfaceMesh, a_X, a_R=None, a_V=None):

    for t in xrange(a_SurfaceMesh.GetNumberOfCells()):
        

    return 0.0

    
