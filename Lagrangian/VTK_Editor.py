import vtk
import numpy as np
import numpy as np
import sys
from vtk.util.numpy_support import vtk_to_numpy
from matplotlib import pyplot as plt
from vtk.numpy_interface import dataset_adapter as dsa
import functools as ft
# noinspection PyUnresolvedReferences
import vtkmodules.vtkInteractionStyle
# noinspection PyUnresolvedReferences
import vtkmodules.vtkRenderingOpenGL2
from vtkmodules.vtkCommonColor import vtkNamedColors
from vtkmodules.vtkCommonCore import vtkPoints
from vtkmodules.vtkCommonDataModel import (
    vtkCellArray,
    vtkPolyData,
    vtkPolyLine
)
from vtkmodules.vtkRenderingCore import (
    vtkActor,
    vtkPolyDataMapper,
    vtkRenderWindow,
    vtkRenderWindowInteractor,
    vtkRenderer
)
import math
sys.path.append('/mnt/d/GVTK/GVTK/Lagrangian/')
try:
    from classGVTKGrid import *
except:
    sys.exit("Import Failed")

baseDir = '/mnt/c/MTH/317_HtoB/'
file_path = '/mnt/c/MTH/317_HtoB/P317_CL_Final.vtp'  # Replace with your file path
reader = vtk.vtkXMLPolyDataReader()
reader.SetFileName(file_path)
reader.Update()

centerline = reader.GetOutput()

# Find the longest line
max_length = 0.0
max_length_index = -1

for i in range(centerline.GetNumberOfLines()):
    line = centerline.GetCell(i)
    length = 0.0

    # Calculate the length of each line
    points = line.GetPoints()
    for j in range(points.GetNumberOfPoints() - 1):
        p1 = points.GetPoint(j)
        p2 = points.GetPoint(j + 1)
        length += vtk.vtkMath.Distance2BetweenPoints(p1, p2)

    # Check if this line is longer than previous max
    if length > max_length:
        max_length = length
        max_length_index = i
print(max_length_index)
# Create a new polydata to store modified centerline
new_centerline = vtk.vtkPolyData()
new_points = vtk.vtkPoints()
new_lines = vtk.vtkCellArray()

# Copy points and lines excluding the longest line
for i in range(centerline.GetNumberOfPoints()):
    new_points.InsertNextPoint(centerline.GetPoint(i))

for i in range(centerline.GetNumberOfLines()):
    if i != max_length_index:
        line = centerline.GetCell(i)
        new_lines.InsertNextCell(line)


new_centerline.SetPoints(new_points)
new_centerline.SetLines(new_lines)

for i in range(centerline.GetPointData().GetNumberOfArrays()):
    array = centerline.GetPointData().GetArray(i)
    new_centerline.GetPointData().AddArray(array)
# Create a writer to save the modified centerline
writer = vtk.vtkXMLPolyDataWriter()
writer.SetFileName(baseDir + 'edited_centerline.vtp')  # Specify the output file name
writer.SetInputData(new_centerline)
writer.Write()
