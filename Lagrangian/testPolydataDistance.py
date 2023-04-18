import numpy as np
import vtk
import sys

meshFile    = '/Users/debanjanmukherjee/Workfiles/Work/Atheroembolic-Renal-Disease-Simulations/VTU-TSteps/velocity_0.vtu'
surfaceFile = '/Users/debanjanmukherjee/Workfiles/Work/Atheroembolic-Renal-Disease-Simulations/surface.vtp'
outFile     = '/Users/debanjanmukherjee/Workfiles/Work/Atheroembolic-Renal-Disease-Simulations/sdf.vtu'

mReader = vtk.vtkXMLUnstructuredGridReader()
mReader.SetFileName(meshFile)
mReader.Update()
meshData = mReader.GetOutput()

sReader = vtk.vtkXMLPolyDataReader()
sReader.SetFileName(surfaceFile)
sReader.Update()
surfData = sReader.GetOutput()

numGridPoints = meshData.GetNumberOfPoints()

sdf = vtk.vtkImplicitPolyDataDistance()
sdf.SetInput(surfData)

distance = vtk.vtkDoubleArray()
distance.SetName('sdf')
distance.SetNumberOfTuples(numGridPoints)
distance.SetNumberOfComponents(1)

for p in range(numGridPoints):
    sys.stdout.write('Computing Point '+str(p)+' of a total '+str(numGridPoints)+'\r')
    sys.stdout.flush()
    xyz = meshData.GetPoint(p)
    distance.SetValue(p, sdf.EvaluateFunction(xyz))

meshData.GetPointData().AddArray(distance)

mWriter = vtk.vtkXMLUnstructuredGridWriter()
mWriter.SetFileName(outFile)
mWriter.SetInputData(meshData)
mWriter.SetDataModeToBinary()
mWriter.Update()
mWriter.Write()
