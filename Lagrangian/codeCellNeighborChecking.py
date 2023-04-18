# ---------------------------------------------------------------------------------------------
# The objective of this module is to:
# - locate a point in an unstructured mesh using a simple geometrical check
# - compare this brute force calculation to vtk's CellLocator algorithm
#
# Author:       Debanjan Mukherjee, Zachariah Irwin
# Institution:  University of Colorado, Boulder
# Last Edit:    March 2019
# ---------------------------------------------------------------------------------------------

# -----------------------------------
# BEGIN MODULE IMPORTS
# -----------------------------------
from __future__ import print_function
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

def extractDataFromFile(a_FileName):
    if a_FileName.endswith('vtk'):
        reader = vtk.vtkUnstructuredGridReader()
    elif a_FileName.endswith('vtu'):
        reader = vtk.vtkXMLUnstructuredGridReader()

    reader.SetFileName(a_FileName)
    reader.Update()

    return reader.GetOutput()


# --------------------------------------------------------------------------
# Function to create a cell locator object from input mesh data from a file
# -------
# Params:
# -------
# a_MeshData:       data containing the mesh where particles are to be located
# a_LocatorType:    (optional) type of vtk locator being employed
# --------
# Returns:
# --------
# locator:  a vtkLocator object (not sure what this object data entails)
# --------------------------------------------------------------------------
def createCellLocator(a_MeshData, a_LocatorType=None):
    if a_LocatorType is None:
        locator = vtk.vtkCellTreeLocator()
    else:
        if a_LocatorType == 'oct':
            locator = vtk.vtkCellLocator()
        elif a_LocatorType == 'tre':
            locator = vtk.vtkCellTreeLocator()
        elif a_LocatorType == 'bsp':
            locator = vtk.vtkModifiedBSPTree()

    locator.SetDataSet(a_MeshData)
    locator.BuildLocator()

    return locator


# --------------------------------------------------------------------------
# Function to create a cell locator object from input mesh data from a file
# using a brute-force location technique
# -------
# Params:
# -------
# a_FileName:       file containing the mesh where particles are to be located
# a_Point:          coordinates of the point in the global configuration
# a_LocatorType:    (optional) type of vtk locator being employed
# --------
# Returns:
# --------
# cellID:           cell in which in the point is located
# --------------------------------------------------------------------------
def createBruteForceLocator(a_MeshData, a_Point, a_NumDim=3):
    locID = -1
    numberOfCells = a_MeshData.GetNumberOfCells()

    for cellID in xrange(numberOfCells):

        cell = a_MeshData.GetCell(cellID)
        cellType = cell.GetCellType()

        if cellType == vtk.VTK_TETRA and a_NumDim == 3:

            idList0 = vtk.vtkIdList()
            idList1 = vtk.vtkIdList()
            idList2 = vtk.vtkIdList()
            idList3 = vtk.vtkIdList()

            face0 = cell.GetFace(0)
            face0IDs = face0.GetPointIds()
            face1 = cell.GetFace(1)
            face1IDs = face1.GetPointIds()
            face2 = cell.GetFace(2)
            face2IDs = face2.GetPointIds()
            face3 = cell.GetFace(3)
            face3IDs = face3.GetPointIds()

            print(cell)

            a_MeshData.GetCellNeighbors(cellID, face0IDs, idList0)
            a_MeshData.GetCellNeighbors(cellID, face1IDs, idList1)
            a_MeshData.GetCellNeighbors(cellID, face2IDs, idList2)
            a_MeshData.GetCellNeighbors(cellID, face3IDs, idList3)

            parametricCoordinates = convertToReference(cell, a_Point)

            if checkInTetrahedron(parametricCoordinates) is True:
                locID = cellID

        elif cellType == vtk.VTK_TRIANGLE and a_NumDim == 2:

            parametricCoordinates = convertToReference(cell, a_Point)

            if checkInTriangle(parametricCoordinates) is True:
                locID = cellID

    return locID


# --------------------------------------------------------------------------
# Function to check whether or not a point is within a given cell
# -------
# Params:
# -------
# a_ParametricCoordinates:       the point's coordinates in the reference
#                                domain of a given cell
# --------
# Returns:
# --------
# boolean
# --------------------------------------------------------------------------
def checkInTetrahedron(a_ParametricCoordinates):
    xi = a_ParametricCoordinates[0]
    eta = a_ParametricCoordinates[1]
    zeta = a_ParametricCoordinates[2]

    if (xi >= 0) and (eta >= 0) and (zeta >= 0) and (1 - xi - eta - zeta >= 0):
        return True
    else:
        return False


def checkInTriangle(a_ParametricCoordinates):
    return 0


# --------------------------------------------------------------------------
# Function to convert a point to a cell's reference coordinates
# -------
# Params:
# -------
# a_Cell:                   vtk cell object
# a_Point:                  point coordinates in global configuration
# --------
# Returns:
# --------
# a_ParametricCoordinates:  point coordinates in reference configuration
# --------------------------------------------------------------------------
def convertToReference(a_Cell, a_Point):
    cellType = a_Cell.GetCellType()  # redundant variable assignment with current debugging setup

    if cellType == vtk.VTK_TETRA:
        zP = a_Point[2]
        yP = a_Point[1]
        xP = a_Point[0]

        cellPoints = a_Cell.GetPoints()

        xyz1 = cellPoints.GetPoint(0)
        xyz2 = cellPoints.GetPoint(1)
        xyz3 = cellPoints.GetPoint(2)
        xyz4 = cellPoints.GetPoint(3)

        z4 = xyz4[2]
        z3 = xyz3[2]
        z2 = xyz2[2]
        z1 = xyz1[2]

        y4 = xyz4[1]
        y3 = xyz3[1]
        y2 = xyz2[1]
        y1 = xyz1[1]

        x4 = xyz4[0]
        x3 = xyz3[0]
        x2 = xyz2[0]
        x1 = xyz1[0]

        a11 = (z4 - z1) * (y3 - y4) - (z3 - z4) * (y4 - y1)
        a21 = (z4 - z1) * (y1 - y2) - (z1 - z2) * (y4 - y1)
        a31 = (z2 - z3) * (y1 - y2) - (z1 - z2) * (y2 - y3)
        a12 = (x4 - x1) * (z3 - z4) - (x3 - x4) * (z4 - z1)
        a22 = (x4 - x1) * (z1 - z2) - (x1 - x2) * (z4 - z1)
        a32 = (x2 - x3) * (z1 - z2) - (x1 - x2) * (z2 - z3)
        a13 = (y4 - y1) * (x3 - x4) - (y3 - y4) * (x4 - x1)
        a23 = (y4 - y1) * (x1 - x2) - (y1 - y2) * (x4 - x1)
        a33 = (y2 - y3) * (x1 - x2) - (y1 - y2) * (x2 - x3)

        V = (x2 - x1) * ((y3 - y1) * (z4 - z1) - (z3 - z1) * (y4 - y1)) + \
            (x3 - x1) * ((y1 - y2) * (z4 - z1) - (z1 - z2) * (y4 - y1)) + \
            (x4 - x1) * ((y2 - y1) * (z3 - z1) - (z2 - z1) * (y3 - y1))

        aMatrix = np.matrix([[a11, a12, a13], [a21, a22, a23], [a31, a32, a33]])
        aMatrix = (1.0 / V) * aMatrix
        x = np.matrix([[xP - x1], [yP - y1], [zP - z1]])

        a_ParametricCoords = np.dot(aMatrix, x)

        return a_ParametricCoords


# -----------------------------------
# END UTILITY FUNCTION DEFINITIONS
# -----------------------------------

# -----------------------------------
# BEGIN COMPUTE SCRIPT DEFINITION
# -----------------------------------

if __name__ == "__main__":
    # meshFile = sys.argv[1]
    # -----------------------------------------------------------------
    # modifying user inputs for testing the brute force locator
    # 1. Input mesh file where location is to be performed
    # 2. The single point coordinates whose location is to be obtained
    # -----------------------------------------------------------------
    meshFile = 'h5.vtk'
    xyz = [0.1, 0.2, 0.3]

    meshData = extractDataFromFile(meshFile)

    # -----------------------------------------
    # find injected point using both locators
    # -----------------------------------------
    t_vtk_1 = time.time()
    locatorObj = createCellLocator(meshData, a_LocatorType='tre')
    cell = locatorObj.FindCell(xyz)
    t_vtk_2 = time.time()
    t_vtk_elapsed = t_vtk_2 - t_vtk_1
    print('Found point in cell #', cell, 'after %1.2f seconds using tree algorithm.' % t_vtk_elapsed)

    t_brute_1 = time.time()
    cell = createBruteForceLocator(meshData, xyz)
    t_brute_2 = time.time()
    t_brute_elapsed = t_brute_2 - t_brute_1
    print('Found point in cell #', cell, 'after %1.2f seconds using brute force algorithm.' % t_brute_elapsed)

# -----------------------------------
# END COMPUTE SCRIPT DEFINITION
# -----------------------------------
