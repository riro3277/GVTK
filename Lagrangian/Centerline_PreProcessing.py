import vtk
import numpy as np
import numpy as np
import sys
from vtk.util.numpy_support import vtk_to_numpy
from matplotlib import pyplot as plt
from vtk.numpy_interface import dataset_adapter as dsa
import functools as ft

sys.path.append('/mnt/d/GVTK/GVTK/Lagrangian/')
try:
    from classGVTKGrid import *
except:
    sys.exit("Import Failed")


def SDF(sdf_file, surface_file, out_File):
    grid = GVTKGridData()
    mesh = sdf_file
    surface = surface_file
    outFile = out_File
    # vtk.vtkPolydataDistanceField(mesh, surface, outFile)
    grid.getDataFromFile(sdf_file, a_FileType ='vtu')
    posP = [17,-276,-828]
    a_locator = grid.buildLocator(a_LocatorType = "oct")
    locatorObj = grid.m_Locator

    d   = grid.gridInterpolateNodalBasis(posP,'sdf', a_GetGradient=False, a_GetStatus=False, a_DataDim=1, a_DataType='double')
    print(d)
    return(grid, locatorObj)

pathtoCL = "/mnt/c/MTH/317_CL_3D.vtp"

reader = vtk.vtkXMLPolyDataReader()
reader.SetFileName(pathtoCL)
reader.Update()

data = reader.GetOutput()
numPts  = data.GetNumberOfPoints()

P_coords = reader.GetOutput().GetPoints().GetData()
coords = vtk_to_numpy(P_coords)
print(len(coords))
print(coords[0])
grid, locator = SDF("/mnt/c/MTH/SDF_317_aorta.vtu", "/mnt/c/MTH/317_AORTA_SURF.stl", "/mnt/c/MTH/OUT.vtu")
SDF_val = grid.gridInterpolateNodalBasis(coords[0],'sdf', a_GetGradient=False, a_GetStatus=False, a_DataDim=1, a_DataType='double')
print(SDF_val)

mesh = reader.GetOutput()
meshNew = dsa.WrapDataObject(mesh)
hsuArray = np.zeros((5,40,40))
array = vtk.vtkFloatArray()
array.SetNumberOfComponents(0)
array.SetName("Test")
for i in range(reader.GetNumberOfPoints()):
    array.InsertNextValue(0)
meshNew.GetPointData().AddArray(array)
writer = vtk.vtkXMLPolyDataWriter()
writer.SetFileName("/mnt/c/MTH/New.vtp")
writer.SetInputData(meshNew.VTKObject)
writer.Write()
