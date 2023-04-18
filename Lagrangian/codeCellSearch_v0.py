# ---------------------------------------------------------------------------------------------
# The objective of this module is to:
# - generate an array of cell data given a starting cell
# - expand this array to find a given point within an unstructured mesh
#
# Author:       Zachariah Irwin
# Institution:  University of Colorado, Boulder
# Last Edit:    March 2019
# ---------------------------------------------------------------------------------------------

# -----------------------------------
# BEGIN MODULE IMPORTS
# -----------------------------------
from __future__ import print_function
from codeBruteLocator_v1 import *
import sys, os
import vtk
import numpy as np
import time
# -----------------------------------
# END MODULE IMPORTS
# -----------------------------------

# -----------------------------------
# BEGIN UTILITY FUNCTION DEFINITIONS
# -----------------------------------

# --------------------------------------------------------------------------
# Function to build a list of neighboring Cell IDs of a provided Cell ID
# -------
# Params:
# -------
# a_CellID:                 vtk Cell ID
# a_MeshData:               vtkUnstructuredGrid object
# a_NumLevels:              integer corresponding to depth of a_CellIDList
# --------
# Returns:
# --------
# a_CellIDList:             list of neighboring Cell IDs
# --------------------------------------------------------------------------
def cellArrayBuilder(a_CellID, a_MeshData, a_NumLevels):

    currentLevel = 0
    a_CellIDList = []
    numberOfCells = a_MeshData.GetNumberOfCells()
    internalIDList = vtk.vtkIdList()
    for internalCellID in xrange(numberOfCells):
        internalIDList.InsertNextId(internalCellID)

    while currentLevel <= a_NumLevels:

        for cellID in a_CellIDList:

            cell           = a_MeshData.GetCell(cellID)
            cellPoints     = cell.GetPoints()
            localNeighbors = a_MeshData.GetCellNeighbors(cellID, cellPoints, internalIDList)

            a_CellIDList.append(localNeighbors)
            localNeighbors = set([id for id in localNeighbors if localNeighbors.count(id) > 1])
            localNeighbors = list(localNeighbors)


    return a_CellIDList
