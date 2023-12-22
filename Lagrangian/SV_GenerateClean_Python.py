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
pathtoCL = baseDir + "P317_CL.vtp"
pathtoSurf = baseDir + "P317_SURF.stl"
pathtoSDF = baseDir + "P317_SDF.vtu"
pathOutput = baseDir + "OUT.vtu"
pathFinal = baseDir + "P317_CL_Radv2.vtp"
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

def check_point_line_membership(centerline, point_id_to_check):
    num_cells = centerline.GetNumberOfCells()
    occurrences = 0

    # Iterate through each cell (line segment)
    for i in range(num_cells):
        cell = centerline.GetCell(i)
        point_ids = cell.GetPointIds()

        # Check if the point_id_to_check exists in the current cell
        if point_id_to_check in [point_ids.GetId(j) for j in range(point_ids.GetNumberOfIds())]:
            occurrences += 1

    return occurrences
def check_line_points(line_arr, point_pos):
    occ = 0
    for i in range(len(line_arr)):
        line_points = line_arr[i]
        if point_pos == line_points[0][0] or point_pos == line_points[0][1]:
            occ += 1
    return occ

def count_surrounding_points(centerline, point_id):
    num_points = centerline.GetNumberOfPoints()

    # Create an array to store neighboring point IDs
    neighbors = set()

    # Iterate through each cell (line segment)
    for i in range(centerline.GetNumberOfCells()):
        cell = centerline.GetCell(i)
        point_ids = cell.GetPointIds()

        # Check if the point_id is part of the current cell
        if point_id in [point_ids.GetId(j) for j in range(point_ids.GetNumberOfIds())]:
            # Add neighboring points of the current cell to the set
            for j in range(point_ids.GetNumberOfIds()):
                if point_ids.GetId(j) != point_id:
                    neighbors.add(point_ids.GetId(j))

    # Number of surrounding points (excluding the original point itself)
    surrounding_points_count = len(neighbors)
    return surrounding_points_count
def Generate_Clean(polydata):
    points = vtk.vtkPoints()
    lines = vtk.vtkCellArray()
    pointIds = vtk.vtkIdList()
    all_points = polydata.GetPoints()

    cleaner = vtk.vtkCleanPolyData()
    cleaner.SetInputData(polydata)
    cleaner.PointMergingOn()
    cleaner.Update()
    centerline = cleaner.GetOutput()
    n_centerlines = centerline.GetNumberOfCells()
    print("cells", n_centerlines)

    centId = vtk.vtkIntArray()
    centId.SetName("CenterlineID")
    centId.SetNumberOfValues(centerline.GetNumberOfPoints() * n_centerlines)
    centId.SetNumberOfComponents(n_centerlines)
    centId.Fill(0)

    points.InsertNextPoint(centerline.GetPoint(0))
    pointIds.InsertNextId(0)
    start_PIDs = []
    start_locs = []
    end_PIDs = []
    end_locs = []
    line_points_arr = []
    for c in range(n_centerlines):
        cell = centerline.GetCell(c)
        cell_points = cell.GetNumberOfPoints()
        all_cell_points = cell.GetPoints()
        # print("cell points", cell_points)
        for p in range(cell_points):
            i = cell.GetPointId(p)
            if p == cell_points-1:
                pos = all_cell_points.GetPoint(p)
                id = cell.GetPointId(p)
                end_PIDs.append((id,pos))
            if p == 0:
                pos = all_cell_points.GetPoint(p)
                id = cell.GetPointId(p)
                start_PIDs.append((id,pos))

            # print("PID", i, "Neighbord", count_surrounding_points(centerline, p))
            if pointIds.IsId(i) == -1 and i>0:
                points.InsertNextPoint(centerline.GetPoint(i))
                pointIds.InsertNextId(i)
                id_prev = pointIds.IsId(cell.GetPointId(p-1))
                id_this = pointIds.IsId(cell.GetPointId(p))
                pos2 = all_cell_points.GetPoint(p)
                pos1 = all_cell_points.GetPoint(p-1)
                line_points_arr.append([(pos1, pos2)])
                # print("id_prev", "id_this", pos2, pos1)
                line = vtk.vtkLine()
                line.GetPointIds().SetId(0, id_prev)
                line.GetPointIds().SetId(1, id_this)
                lines.InsertNextCell(line)
            else:
                centId.SetComponent(pointIds.IsId(i), c, 1)
    newdata = vtk.vtkPolyData()
    newdata.SetPoints(points)
    newdata.SetLines(lines)
    newdata.Modified()

    #Check end points of each branch to see if thereis either 0 line connections or >1. If there is 1, add a line
    # print(line_points_arr)
    checkdata = newdata
    checkpoints = vtk.vtkPoints()
    newpoints = newdata.GetPoints()
    newpoints_number = newdata.GetNumberOfPoints()
    checkout = vtk.vtkPolyData()
    all_cells = newdata.GetNumberOfCells()
    locator = vtk.vtkPointLocator()
    locator.SetDataSet(newdata)
    locator.BuildLocator()

    # Find the closest points to the target point using the locator
    closest_points = vtk.vtkIdList()
    lines = newdata.GetLines()

    # Create a vtkIdList to hold the IDs of the points in the line
    point_ids = vtk.vtkIdList()
    branch_end = []
    print("points", newpoints_number)
    missing_line = []
    for i in range(len(start_PIDs)):
        line_count = check_point_line_membership(newdata, start_PIDs[i][0])
        if line_count == 2:
            missing_line.append((start_PIDs[i][0], start_PIDs[i][1]))

    for c in range(lines.GetNumberOfCells()):
        cell = newdata.GetCell(c)
        if cell.GetCellType() == vtk.VTK_LINE:
            point_ids = cell.GetPointIds()
            all_cell_points = cell.GetPoints()
            line_points = []
            for i in range(point_ids.GetNumberOfIds()):

                pid = point_ids.GetId(i)
                point = points.GetPoint(pid)
                line_count = check_line_points(line_points_arr, point)
                # print(line_count)
                line_points.append(point)
                if line_count == 1:
                    branch_end.append((pid, point))
                    checkpoints.InsertNextPoint(point)


    for i in range(len(missing_line)):
        dist = []
        point = missing_line[i][1]
        missing_pid = missing_line[i][0]
        print(missing_pid)
        # sys.exit()
        for j in range(len(branch_end)):
            branch_end_point = branch_end[j][1]
            # print("point", point)
            # print("branch_end_point", branch_end_point)
            val = np.sqrt(((point[0]-branch_end_point[0])**2) + (((point[1]-branch_end_point[1]))**2) + (((point[2]-branch_end_point[2]))**2))
            if val != 0 or val != 0.0:
                dist.append(val)
            else:
                pass
        min_ele= dist[0]
        for i in range(1, len(dist)):
            if dist[i] < min_ele:
                min_ele = dist[i]
        min_pos = dist.index(min_ele)
        print("min_ele", min_ele)
        print("min_pos", min_pos)
        print("dist", dist)

        end_id = branch_end[min_pos][0]
        line = vtk.vtkLine()
        line.GetPointIds().SetId(0, missing_pid)
        line.GetPointIds().SetId(1, end_id)
        lines.InsertNextCell(line)

    newdata = vtk.vtkPolyData()
    newdata.SetPoints(points)
    newdata.SetLines(lines)
    newdata.Modified()

    del branch_end
    del start_PIDs
    del end_PIDs
    del line_points_arr






    checkout.SetPoints(checkpoints)

    writer = vtk.vtkXMLPolyDataWriter()
    writer.SetFileName(baseDir + "SV_Python_CL.vtp")  # Output file name
    writer.SetInputData(newdata)
    writer.Write()
    writer = vtk.vtkXMLPolyDataWriter()
    writer.SetFileName(baseDir + "Check_Ends.vtp")  # Output file name
    writer.SetInputData(checkout)
    writer.Write()
if __name__ == "__main__":
    reader = vtk.vtkXMLPolyDataReader()
    reader.SetFileName(pathtoCL)
    reader.Update()

    data = reader.GetOutput()
    data = add_points_between(data, 10)
    Generate_Clean(data)
