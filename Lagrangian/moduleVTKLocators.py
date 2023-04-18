import sys, os
import vtk
import numpy as np

#--------------------------------------------------------------
# create a cell locator object from input mesh data from a file
#--------------------------------------------------------------
def createCellLocator(a_FileName, a_LocatorType=None, a_FoamBlockID=0):

    if a_FileName.endswith('vtu'):
        reader = vtk.vtkXMLUnstructuredGridReader()
        reader.SetFileName(a_FileName)
        reader.Update()
        gridData = reader.GetOutput()
    elif a_FileName.endswith('vtp'):
        reader = vtk.vtkXMLPolyDataReader()
        reader.SetFileName(a_FileName)
        reader.Update()
        gridData = reader.GetOutput()
    elif a_FileName.endswith('foam'):
        reader = vtk.vtkOpenFOAMReader()
        reader.SetFileName(a_FileName)
        reader.SkipZeroTimeOn()
        reader.Update()
        foamData = reader.GetOutput()
        numBlock = foamData.GetNumberOfBlocks()
        gridData = foamData.GetBlock(a_FoamBlockID)

    if a_LocatorType is None:
        locator = vtk.vtkCellTreeLocator()
    else:
        if a_LocatorType == 'oct':
            locator = vtk.vtkCellLocator()
        elif a_LocatorType == 'tre':
            locator = vtk.vtkCellTreeLocator()
        elif a_LocatorType == 'bsp':
            locator = vtk.vtkModifiedBSPTree()

    locator.SetDataSet(gridData)
    locator.BuildLocator()

    return locator

def createCellLocatorFromData(a_GridData, a_LocatorType=None):

    if (a_GridData.IsA('vtkUnstructuredGrid') != 1):
        sys.exit("Module: moduleVTKLocators, Function: createCellLocator - grid data should be unstructured grid")
    else:
        if a_LocatorType is None:
            locator = vtk.vtkCellTreeLocator()
        else:
            if a_LocatorType == 'oct':
                locator = vtk.vtkCellLocator()
            elif a_LocatorType == 'tre':
                locator = vtk.vtkCellTreeLocator()
            elif a_LocatorType == 'bsp':
                locator = vtk.vtkModifiedBSPTree()

        locator.SetDataSet(a_GridData)
        locator.BuildLocator()
        return locator
