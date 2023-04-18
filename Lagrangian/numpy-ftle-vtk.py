import numpy as np
import sys, time, vtk
from numpy import linalg as lan

try:
    import classGVTKFTLE as M_FTLE
except ImportError:
    sys.exit('Could not import user define modules')

start = time.time()

ftleInput = M_FTLE.FTLEInputs(sys.argv[1])
ftleInput.readInputFile()

flowPath          = ftleInput.m_FlowDataDirectory
flowFile          = ftleInput.m_FlowDataFile
vtkOutputFile     = ftleInput.m_OutputFile
ftleChoice        = ftleInput.m_FTLEChoice
isScaledTime      = ftleInput.m_ScaledTime

reader  = vtk.vtkXMLUnstructuredGridReader()
reader.SetFileName(flowPath+flowFile)
reader.Update()

flowData    = reader.GetOutput()
velocity    = flowData.GetPointData().GetArray(ftleInput.m_FlowFieldName)

locator = vtk.vtkCellTreeLocator()
locator.SetDataSet(flowData)
locator.BuildLocator()

xMin    = ftleInput.m_xMin
xMax    = ftleInput.m_xMax
yMin    = ftleInput.m_yMin
yMax    = ftleInput.m_yMax
zMin    = ftleInput.m_zMin
zMax    = ftleInput.m_zMax
xN      = ftleInput.m_xBins
yN      = ftleInput.m_yBins
zN      = ftleInput.m_zBins
t0      = ftleInput.m_StartTime
dT      = ftleInput.m_dT
tEnd    = ftleInput.m_StopTime

gridX   = np.zeros((xN,yN,zN), dtype=np.float32)
gridY   = np.zeros((xN,yN,zN), dtype=np.float32)
gridZ   = np.zeros((xN,yN,zN), dtype=np.float32)
inGrid  = np.zeros((xN,yN,zN), dtype=np.int32) + 1
FTLE    = np.zeros((xN,yN,zN), dtype=np.float32) + 1.0

defGrad = np.zeros((3,3), dtype=np.float64)

dX      = (xMax - xMin)/float(xN - 1)
dY      = (yMax - yMin)/float(yN - 1)
dZ      = (zMax - zMin)/float(zN - 1)

for k in range(zN):
    for j in range(yN):
        for i in range(xN):

            gridX[i,j,k] = xMin + float(i)*dX
            gridY[i,j,k] = yMin + float(j)*dY
            gridZ[i,j,k] = zMin + float(k)*dZ

gridFTLE_X  = np.copy(gridX)
gridFTLE_Y  = np.copy(gridY)
gridFTLE_Z  = np.copy(gridZ)

t = t0
vel = np.zeros(3,dtype=np.float64)

while t < tEnd:

    if t < tEnd-dT:
        print("Time Step: {0:.4f} Out Of Total Time: {1:.2f}".format(t,tEnd), end='\r')
    else:
        print("Time Step: {0:.4f} Out Of Total Time: {1:.2f}".format(t,tEnd))

    for k in range(zN):
        for j in range(yN):
            for i in range(xN):
                point   = [gridX[i,j,k], gridY[i,j,k], gridZ[i,j,k]]
                cell    = locator.FindCell(point)

                cellPtIds = vtk.vtkIdList()

                if cell != -1:

                    flowData.GetCellPoints(cell, cellPtIds)

                    velNode1 = velocity.GetTuple3(cellPtIds.GetId(0))
                    velNode2 = velocity.GetTuple3(cellPtIds.GetId(1))
                    velNode3 = velocity.GetTuple3(cellPtIds.GetId(2))
                    velNode4 = velocity.GetTuple3(cellPtIds.GetId(3))

                    vel[0]  = (1.0/4.0)*(velNode1[0]+velNode2[0]+velNode3[0]+velNode4[0])
                    vel[1]  = (1.0/4.0)*(velNode1[1]+velNode2[1]+velNode3[1]+velNode4[1])
                    vel[2]  = (1.0/4.0)*(velNode1[2]+velNode2[2]+velNode3[2]+velNode4[2])

                elif cell == -1:

                    inGrid[i,j,k] = -1

                    vel[0]  = 0.0
                    vel[1]  = 0.0
                    vel[2]  = 0.0

                gridX[i,j,k] = point[0] + vel[0]*dT
                gridY[i,j,k] = point[1] + vel[1]*dT
                gridZ[i,j,k] = point[2] + vel[2]*dT

    t = t + dT

#------------------------------------------------------------------
# comuting deformation gradient, Cauchy-Green tensor and FTLE field
# ----------------
# general stencil:
#-----------------
# F[0,0] = (x_i_plus - x_i_minus)/(x0_i_plus - x0_i_minus)
# F[0,1] = (x_j_plus - x_j_minus)/(y0_j_plus - y0_j_minus)
# F[0,2] = (x_k_plus - x_k_minus)/(z0_k_plus - z0_k_minus)
# F[1,0] = (y_i_plus - y_i_minus)/(x0_i_plus - x0_i_minus)
# F[1,1] = (y_j_plus - y_j_minus)/(y0_j_plus - y0_j_minus)
# F[1,2] = (y_k_plus - y_k_minus)/(z0_k_plus - z0_k_minus)
# F[2,0] = (z_i_plus - z_i_minus)/(x0_i_plus - x0_i_minus)
# F[2,1] = (z_j_plus - z_j_minus)/(y0_j_plus - y0_j_minus)
# F[2,2] = (z_k_plus - z_k_minus)/(z0_k_plus - z0_k_minus)
#------------------------------------------------------------------
for k in range(zN):
    for j in range(yN):
        for i in range(xN):

            if i == 0:
                x_i_minus   = gridX[i,j,k]
                y_i_minus   = gridY[i,j,k]
                z_i_minus   = gridZ[i,j,k]
                x_i_plus    = gridX[i+1,j,k]
                y_i_plus    = gridY[i+1,j,k]
                z_i_plus    = gridZ[i+1,j,k]
                ddX         = 1.0
            elif i == (xN-1):
                x_i_minus   = gridX[i-1,j,k]
                y_i_minus   = gridY[i-1,j,k]
                z_i_minus   = gridZ[i-1,j,k]
                x_i_plus    = gridX[i,j,k]
                y_i_plus    = gridY[i,j,k]
                z_i_plus    = gridZ[i,j,k]
                ddX         = 1.0
            else:
                x_i_minus   = gridX[i-1,j,k]
                y_i_minus   = gridY[i-1,j,k]
                z_i_minus   = gridZ[i-1,j,k]
                x_i_plus    = gridX[i+1,j,k]
                y_i_plus    = gridY[i+1,j,k]
                z_i_plus    = gridZ[i+1,j,k]
                ddX         = 2.0

            if j == 0:
                x_j_minus   = gridX[i,j,k]
                y_j_minus   = gridY[i,j,k]
                z_j_minus   = gridZ[i,j,k]
                x_j_plus    = gridX[i,j+1,k]
                y_j_plus    = gridY[i,j+1,k]
                z_j_plus    = gridZ[i,j+1,k]
                ddY         = 1.0
            elif j == (yN-1):
                x_j_minus   = gridX[i,j-1,k]
                y_j_minus   = gridY[i,j-1,k]
                z_j_minus   = gridZ[i,j-1,k]
                x_j_plus    = gridX[i,j,k]
                y_j_plus    = gridY[i,j,k]
                z_j_plus    = gridZ[i,j,k]
                ddY         = 1.0
            else:
                x_j_minus   = gridX[i,j-1,k]
                y_j_minus   = gridY[i,j-1,k]
                z_j_minus   = gridZ[i,j-1,k]
                x_j_plus    = gridX[i,j+1,k]
                y_j_plus    = gridY[i,j+1,k]
                z_j_plus    = gridZ[i,j+1,k]
                ddY         = 2.0

            if k == 0:
                x_k_minus   = gridX[i,j,k]
                y_k_minus   = gridY[i,j,k]
                z_k_minus   = gridZ[i,j,k]
                x_k_plus    = gridX[i,j,k+1]
                y_k_plus    = gridY[i,j,k+1]
                z_k_plus    = gridZ[i,j,k+1]
                ddZ         = 1.0
            elif k == (zN-1):
                x_k_minus   = gridX[i,j,k-1]
                y_k_minus   = gridY[i,j,k-1]
                z_k_minus   = gridZ[i,j,k-1]
                x_k_plus    = gridX[i,j,k]
                y_k_plus    = gridY[i,j,k]
                z_k_plus    = gridZ[i,j,k]
                ddZ         = 1.0
            else:
                x_k_minus   = gridX[i,j,k-1]
                y_k_minus   = gridY[i,j,k-1]
                z_k_minus   = gridZ[i,j,k-1]
                x_k_plus    = gridX[i,j,k+1]
                y_k_plus    = gridY[i,j,k+1]
                z_k_plus    = gridZ[i,j,k+1]
                ddZ         = 2.0

            defGrad[0,0]    = (x_i_plus - x_i_minus)/(ddX*dX)
            defGrad[0,1]    = (x_j_plus - x_j_minus)/(ddY*dY)
            defGrad[0,2]    = (x_k_plus - x_k_minus)/(ddZ*dZ)
            defGrad[1,0]    = (y_i_minus - y_i_plus)/(ddX*dX)
            defGrad[1,1]    = (y_j_minus - y_j_plus)/(ddY*dY)
            defGrad[1,2]    = (y_k_minus - y_k_plus)/(ddZ*dZ)
            defGrad[2,0]    = (z_i_minus - z_i_plus)/(ddX*dX)
            defGrad[2,1]    = (z_j_minus - z_j_plus)/(ddY*dY)
            defGrad[2,2]    = (z_k_minus - z_k_plus)/(ddZ*dZ)

            cauchyGreen     = np.dot(defGrad.T, defGrad)
            eigVal, eigVec  = lan.eigh(cauchyGreen)
            lambdaVal       = np.max(eigVal)

            if ftleChoice == 'fwd':
                if lambdaVal > 1.0-1.0e-14:
                    if isScaledTime:
                        FTLE[i,j,k] = 0.5*np.log(lambdaVal)/(tEnd-t0)
                    else:
                        FTLE[i,j,k] = 0.5*np.log(lambdaVal)
                else:
                    FTLE[i,j,k]     = 0.0
            elif ftleChoice == 'bwd':
                if lambdaVal > 1.0 - 1.0e-14:
                    if isScaledTime:
                        FTLE[i,j]   = 0.5*np.log(lambdaVal)/(t0-tEnd)
                    else:
                        FTLE[i,j]   = 0.5*np.log(lambdaVal)
                else:
                    FTLE[i,j]       = 0.0

vtkGrid = vtk.vtkStructuredGrid()
vtkGrid.SetDimensions([xN,yN,zN])

points = vtk.vtkPoints()
points.SetNumberOfPoints(xN*yN*zN)

ftleArr = vtk.vtkDoubleArray()
ftleArr.SetName('ftle')
ftleArr.SetNumberOfComponents(1)
ftleArr.SetNumberOfTuples(xN*yN*zN)

offset = 0

for k in range(zN):
    for j in range(yN):
        for i in range(xN):
            x = gridFTLE_X[i,j,k]
            y = gridFTLE_Y[i,j,k]
            z = gridFTLE_Z[i,j,k]
            offset = i + xN*j + xN*yN*k
            points.InsertPoint(offset, [x,y,z])
            ftleArr.SetTuple1(offset, FTLE[i,j,k])

vtkGrid.SetPoints(points)
vtkGrid.GetPointData().AddArray(ftleArr)

writer  = vtk.vtkStructuredGridWriter()
writer.SetFileName(flowPath+vtkOutputFile)
writer.SetInputData(vtkGrid)
writer.Update()
writer.Write()

stop = time.time()

print('Total Calcualtion Runtime:', (stop-start)/3600, 'hours\n')
