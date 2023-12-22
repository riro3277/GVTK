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
baseDir = "/mnt/c/MTH/317_HtoB/"
pathtoCL = baseDir + "SV_Python_CL.vtp"
pathtoSurf = baseDir + "P317_SURF.stl"
pathtoSDF = baseDir + "P317_SDF.vtu"
pathOutput = baseDir + "OUT.vtu"
pathFinal = baseDir + "SV_Python_CLv2.vtp"


def BFS_ID(polydata, start_point_id):
    visited = set()
    queue = [start_point_id]
    point_ids = vtk.vtkIntArray()
    point_ids.SetNumberOfComponents(1)
    point_ids.SetName("PointIDs")
    point_positions = polydata.GetPoints()

    new_id = 0
    while queue:
        current_point_id = queue.pop(0)
        if current_point_id not in visited:
            visited.add(current_point_id)
            point_ids.InsertNextValue(new_id)
            point_position = point_positions.GetPoint(current_point_id)
            print(f"New ID: {new_id}, Position: {point_position}")
            new_id += 1

            neighbors = vtk.vtkIdList()
            polydata.GetPointCells(current_point_id, neighbors)

            for i in range(neighbors.GetNumberOfIds()):
                cell_id = neighbors.GetId(i)
                cell = polydata.GetCell(cell_id)
                point_ids_of_cell = cell.GetPointIds()

                for j in range(point_ids_of_cell.GetNumberOfIds()):
                    neighbor_point_id = point_ids_of_cell.GetId(j)
                    if neighbor_point_id not in visited:
                        queue.append(neighbor_point_id)

    return point_ids

def find_nearest_points(tree, point, num_points=2):
    """
    Find nearest points to the given point using vtkKdTree.
    Returns the indices of the nearest points.
    """

    locator = vtk.vtkKdTree()
    locator.SetDataSet(tree)
    locator.BuildLocator()

    result = vtk.vtkIdList()
    locator.FindClosestNPoints(num_points, point, result)

    return [result.GetId(i) for i in range(result.GetNumberOfIds())]



def clean_data(data):
    point_data = data.GetPointData()
    for i in range(point_data.GetNumberOfArrays()):
        point_data.RemoveArray("Curvature")
        point_data.RemoveArray("FrenetNormal")
        point_data.RemoveArray("FrenetTangent")
        point_data.RemoveArray("FrenetBinormal")
        point_data.RemoveArray("Torsion")
        point_data.RemoveArray("Radius")
        # point_data.RemoveArray(i)

    # Remove all cell data arrays
    cell_data = data.GetCellData()
    for i in range(cell_data.GetNumberOfArrays()):
        # cell_data.RemoveArray(i)
        # cell_data.RemoveArray("Length")
        cell_data.RemoveArray("Tortuosity")


    # Remove all field data arrays
    field_data = data.GetFieldData()
    for i in range(field_data.GetNumberOfArrays()):
        field_data.RemoveArray(i)

def add_points_between(centerline, num_points_to_add):

    num_cells = centerline.GetNumberOfCells()

    # Define the number of points you want to add between the existing points
    num_points_to_add = 10  # Change this value according to your needs

    # Create a new polydata to store the interpolated points
    new_centerline = vtk.vtkPolyData()
    points = vtk.vtkPoints()
    new_centerline.SetPoints(points)

    # Create cells object to store lines connecting the points
    polylines = vtk.vtkCellArray()
    new_centerline.SetLines(polylines)
    for i in range(centerline.GetCellData().GetNumberOfArrays()):
        new_centerline.GetCellData().AddArray(centerline.GetCellData().GetArray(i))
    # Iterate through each cell in the centerline
    for cell_id in range(num_cells):
        cell = centerline.GetCell(cell_id)
        num_points_in_cell = cell.GetNumberOfPoints()

        # Retrieve existing points in the current cell
        existing_points = [cell.GetPoints().GetPoint(i) for i in range(num_points_in_cell)]

        # Interpolate points between the existing points
        interpolated_points = []
        for i in range(num_points_in_cell - 1):
            for j in range(num_points_to_add + 1):  # +1 to include the end point
                t = j / float(num_points_to_add + 1)
                new_point = [
                    existing_points[i][k] + t * (existing_points[i + 1][k] - existing_points[i][k])
                    for k in range(3)
                ]
                interpolated_points.append(new_point)
            if i == num_points_in_cell - 2:
                # print("here")
                last_point = existing_points[i+1]
                interpolated_points.append(last_point)

        # Add points and connect them using a polyline within the current cell
        polyline = vtk.vtkPolyLine()
        for new_point in interpolated_points:
            point_id = points.InsertNextPoint(new_point)
            polyline.GetPointIds().InsertNextId(point_id)
        polylines.InsertNextCell(polyline)
    # writer = vtk.vtkXMLPolyDataWriter()
    # writer.SetFileName(baseDir + "interpolated_centerline.vtp")  # Replace with your desired output file name
    # writer.SetInputData(new_centerline)
    # writer.Write()
    print("More Points Added")
    return new_centerline

def BranchID(polydata):
    cell_data = polydata.GetCellData()
    num_components = polydata.GetCellData().GetArray("Topology").GetNumberOfComponents()
    num_cells = polydata.GetNumberOfCells()
    bidArray = vtk.vtkDoubleArray()
    bidArray.SetNumberOfComponents(1)
    bidArray.SetNumberOfTuples(num_cells)
    bidArray.SetName("BranchId")
    magnitude_arr = []
    bid = 1
    bid_arr = np.zeros(num_cells)

    for i in range(num_cells-1):
        value1 = [polydata.GetCellData().GetArray("Topology").GetValue(i * num_components + j) for j in range(num_components)]
        value2 = [polydata.GetCellData().GetArray("Topology").GetValue((i+1) * num_components + j) for j in range(num_components)]
        if bid_arr[i] == 0:
            if value1[1] == value2[0]:
                bid_arr[i] = bid
                bid_arr[i+1] = bid
            else:
                bid_arr[i] = bid
            bid += 1
        else:
            bid = bid
        if i == num_cells-2:
            if bid_arr[i+1] == 0:
                bid_arr[i+1] = bid_arr[i] + 1

        bidArray.SetValue(i, bid_arr[i])
    polydata.GetCellData().AddArray(bidArray)

    cell_data_array = polydata.GetCellData().GetArray("BranchID")
    cell_to_point = vtk.vtkCellDataToPointData()
    cell_to_point.SetInputData(polydata)
    cell_to_point.PassCellDataOn()
    cell_to_point.Update()
    result_polydata = cell_to_point.GetOutput()

    # Retrieve the point data array
    polydata = result_polydata
    # writer = vtk.vtkXMLPolyDataWriter()
    # writer.SetFileName(baseDir + "branchid_values.vtp")  # Replace with your desired file name
    # writer.SetInputData(polydata)
    # writer.Write()
    return polydata

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
    # print(d)
    return(grid, locatorObj)



reader = vtk.vtkXMLPolyDataReader()
reader.SetFileName(pathtoCL)
reader.Update()

data = reader.GetOutput()
numPts  = data.GetNumberOfPoints()

# BFS_ID(data, 0)
# sys.exit()
# clean_data(data)
# add_points_between(data, 10)
#
# data = BranchID(data)
# data = connect_points(data)
# sys.exit()


# Write the modified data (without arrays) to a new VTP file
# writer = vtk.vtkXMLPolyDataWriter()
# writer.SetFileName(baseDir + "modified_file.vtp")  # Output file name
# writer.SetInputData(data)
# writer.Write()

P_coords = data.GetPoints().GetData()
coords = vtk_to_numpy(P_coords)



grid, locator = SDF(pathtoSDF, pathtoSurf, pathOutput)

#
# Create a vtkPoints object and store the points in it
points = vtkPoints()
# mesh = reader.GetOutput()
# meshNew = dsa.WrapDataObject(mesh)
array = vtk.vtkFloatArray()
array.SetNumberOfComponents(0)
array.SetName("MaximumInscribedSphereRadius")
for i in range(data.GetNumberOfPoints()):
    points.InsertNextPoint(coords[i])
    SDF_val = grid.gridInterpolateNodalBasis(coords[i],'sdf', a_GetGradient=False, a_GetStatus=False, a_DataDim=1, a_DataType='double')
    array.InsertNextValue(SDF_val)

# polyLine = vtkPolyLine()
# polyLine.GetPointIds().SetNumberOfIds(numPts)
# for i in range(0, numPts):
#     polyLine.GetPointIds().SetId(i, i)
# # Create a cell array to store the lines in and add the lines to it
# cells = vtkCellArray()
# cells.InsertNextCell(polyLine)
# # Create a polydata to store everything in
# polyData = vtkPolyData()
#
# # Add the points to the dataset
# polyData.SetPoints(points)

# Add the lines to the dataset
# data.SetLines(cells)



data.GetPointData().AddArray(array)
writer = vtk.vtkXMLPolyDataWriter()
writer.SetFileName(pathFinal)
writer.SetInputData(data)
writer.Write()
