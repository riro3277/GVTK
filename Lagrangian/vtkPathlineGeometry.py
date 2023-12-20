"""File: vtkCurveGeometryAnalysis

A module that implements discrete differential geometry operations of polyline objects.

Notes:
------
Function performs the discrete geometry calculations of curvature and torsion of the
3D curve presented as a polyline object and stores it into a polydata VTK file.
The calculation steps are outlined as follows:

- let the points of the polyline be x[k]
- get: ds[k]  = length(x[k+1]-x[k])
- get: T[k]   = (x[k+1]-x[k]) / ds[k]
- get: B[k]   = cross(T[k],T[k+1]) (vector product)
- get: N[k]   = cross(B[k],T[k])
- [ T,B,N defines the discrete Frenet Frame of the curve ]
- get: kappa[k] = length(T[k+1] - T[k]) / ds[k]       (curvature)
- get: tau[k]   = -dot((B[k+1] - B[k]),N[k]) / ds[k]  (torsion)
- [ with obvious restrictions on range of k ]

Author:       Debanjan Mukherjee
Institution:  University of California, Berkeley
Last edited:  June 2018

TODO:
    - Convert the read write functions into vtkGenIO function calls
    - Add a renderer for visualizing the vectors for the seret frenet frame

"""

#---------------------
# Begin module imports
#---------------------
import sys, os
import numpy as np
# import vtkGenIO as VTKIO

try:
    import vtk
    from vtk.util import numpy_support as vnp
except ImportError:
    sys.exit('VTK Modules Not Loaded!')

from geomdl import BSpline, utilities
from geomdl.fitting import interpolate_curve

#-------------------
# End module imports
#-------------------

def getPolyLineObject(a_InputVTKPolyData):

    if a_InputVTKPolyData.endswith('vtk'):
        reader = vtk.vtkPolyDataReader()
    elif a_InputVTKPolyData.endswith('vtp'):
        reader = vtk.vtkXMLPolyDataReader()

    reader.SetFileName(a_InputVTKPolyData)
    reader.Update()

    return reader.GetOutput()

def removeDuplicatePoints(a_PolyData):

    cleanFilter = vtk.vtkCleanPolyData()
    cleanFilter.SetInputData(a_PolyData)
    cleanFilter.Update()

    return cleanFilter.GetOutput()

def interpBSpline(a_PolyData, a_SplineDegree, a_Mode='interp'):
    '''
    Convert points on a_PolyData into the geomdl BSpline object
    '''

    # Get the points from a_PolyData
    numPoints = a_PolyData.GetNumberOfPoints()
    points = [list(a_PolyData.GetPoint(k)) for k in range(numPoints)]

    # Catch invalid input
    try:
        assert(a_Mode == 'ctrlpts' or a_Mode == 'interp')
    except AssertionError:
        sys.exit('interpolation mode must be either `interp` or `ctrlpts`.')

    # Run interpolation
    if a_Mode == 'ctrlpts':
        # Make points on a_PolyData control points for the bspline
        spline = BSpline.Curve()
        spline.degree = a_SplineDegree
        spline.ctrlpts = points
        spline.knotvector = utilities.generate_knot_vector(spline.degree, numPoints)
    elif a_Mode == 'interp':
        # Interpolate a spline through the points
        spline = interpolate_curve(points, a_SplineDegree)
    return spline

def convertBSplineToPolyData(a_Spline, a_NumResampledPoints=2000):

    # Get resampled points long the curve
    a_Spline.sample_size = a_NumResampledPoints
    a_Spline.evaluate()
    points = np.array(a_Spline.evalpts)

    # Create vtk points object
    vtkPoints = vtk.vtkPoints()
    vtkPoints.SetNumberOfPoints(a_NumResampledPoints)
    for i, p in enumerate(points):
        vtkPoints.SetPoint(i, p)

    # Create line connecting the points
    vtkLines = vtk.vtkCellArray()
    for i in range(a_NumResampledPoints-1):
        vtkIdList = vtk.vtkIdList()
        vtkIdList.InsertNextId(i)
        vtkIdList.InsertNextId(i+1)
        vtkLines.InsertNextCell(vtkIdList)

    # Write points and line into polydata and return
    polyData = vtk.vtkPolyData()
    polyData.SetPoints(vtkPoints)
    polyData.SetLines(vtkLines)
    return polyData

def computeCurvatureTorsion(a_PolyData):

    numPoints   = a_PolyData.GetNumberOfPoints()

    ds          = np.zeros(numPoints)

    kappa       = np.zeros(numPoints)
    tau         = np.zeros(numPoints)

    T           = np.zeros((numPoints,3))
    B           = np.zeros((numPoints,3))
    N           = np.zeros((numPoints,3))
    tang        = np.zeros(3)

    #
    # loop through timesteps/number of points to calculate
    # steps 1 and 2 from documentation above
    #
    for k in range(numPoints):

        #
        # unless we reach the last two points on the curve follow the steps in the documentation
        #
        if k != (numPoints-1):

            xyz0    = a_PolyData.GetPoint(k)
            xyz1    = a_PolyData.GetPoint(k+1)
            dist01  = ((xyz0[0] - xyz1[0])**2 + (xyz0[1] - xyz1[1])**2 + (xyz0[2] - xyz1[2])**2)**0.5

            tang[0] = (xyz1[0] - xyz0[0])/dist01
            tang[1] = (xyz1[1] - xyz0[1])/dist01
            tang[2] = (xyz1[2] - xyz0[2])/dist01

            ds[k]   = dist01
            T[k,:]  = tang
            #ds      = np.append(ds, dist01)
            #T       = np.append(T, tang, axis = 0)

        #
        # for the last point on the curve, set the values to be the same as the penultimate point
        #
        else:

            ds[k]   = ds[numPoints-2]
            T[k,:]  = T[numPoints-2,:]
            #ds = np.append(ds, ds[numPoints-2])
            #T  = np.append(T, T[numPoints-2,:], axis = 0)

    #
    # loop through timesteps to calculate steps 3 and 4 from documentation above
    #
    for k in range(numPoints):

        #
        # unless we reach the last two points on the curve follow the steps in the documentation
        #
        if k != (numPoints - 1):

            binorm  = np.cross(T[k,:], T[k+1,:])
            B[k, :] = binorm
            #B       = np.append(B, binorm, axis = 0)

            normal  = np.cross(B[k,:],T[k,:])
            N[k, :] = normal
            #N       = np.append(N, normal, axis = 0)

        #
        # for the last point on the curve, set the values to be the same as the penultimate point
        #
        else:

            B[k,:]  = B[numPoints-2,:]
            N[k,:]  = N[numPoints-2,:]

            #B = np.append(B, B[numPoints-2,:], axis = 0)
            #N = np.append(N, N[numPoints-2,:], axis = 0)

    #
    # loop through timesteps to calculate steps 5 and 6 from documentation above
    #
    for k in range(numPoints):

        #
        # unless we reach the last two points on the curve follow the steps in the documentation
        #
        if k != (numPoints - 1):

            diffT = ((T[k+1,0] - T[k,0])**2 + (T[k+1,1] - T[k,1])**2 + (T[k+1,2] - T[k,2])**2)**0.5

            dotBN = np.dot((B[k+1,:] - B[k,:]),N[k,:])

            kappa[k]    = diffT/ds[k]
            tau[k]      = -dotBN/ds[k]

            #kappa = np.append(kappa, diffT/ds[k])
            #tau   = np.append(tau, -dotBN/ds[k])

        #
        # for the last point on the curve, set the values to be the same as the penultimate point
        #
        else:

            kappa[k]    = kappa[numPoints-2]
            tau[k]      = tau[numPoints-2]

            #kappa = np.append(kappa, kappa[numPoints-2])
            #tau   = np.append(tau, tau[numPoints-2])

    pkappa  = vtk.vtkDoubleArray()
    pkappa.SetName('curvature')
    pkappa.SetNumberOfComponents(1)
    pkappa.SetNumberOfTuples(numPoints)

    ptau    = vtk.vtkDoubleArray()
    ptau.SetName('torsion')
    ptau.SetNumberOfComponents(1)
    ptau.SetNumberOfTuples(numPoints)

    tangents    = vtk.vtkDoubleArray()
    tangents.SetName('tangents')
    tangents.SetNumberOfComponents(3)
    tangents.SetNumberOfTuples(numPoints)

    normals     = vtk.vtkDoubleArray()
    normals.SetName('normals')
    normals.SetNumberOfComponents(3)
    normals.SetNumberOfTuples(numPoints)

    binormals   = vtk.vtkDoubleArray()
    binormals.SetName('binormals')
    binormals.SetNumberOfComponents(3)
    binormals.SetNumberOfTuples(numPoints)

    for k in range(numPoints):
        pkappa.SetTuple1(k, kappa[k])
        ptau.SetTuple1(k,   tau[k])
        tangents.SetTuple3(k,   T[k,0], T[k,1], T[k,2])
        normals.SetTuple3(k,    N[k,0], N[k,1], N[k,2])
        binormals.SetTuple3(k,  B[k,0], B[k,1], B[k,2])

    dataNew = vtk.vtkPolyData()
    dataNew.SetPoints(a_PolyData.GetPoints())
    dataNew.SetLines(a_PolyData.GetLines())
    dataNew.GetPointData().AddArray(pkappa)
    dataNew.GetPointData().AddArray(ptau)
    dataNew.GetPointData().AddArray(tangents)
    dataNew.GetPointData().AddArray(normals)
    dataNew.GetPointData().AddArray(binormals)

    return dataNew

def writePolyLine(a_PolyData, a_FileName):

    if a_FileName.endswith('vtk'):
        writer = vtk.vtkPolyDataWriter()
    elif a_FileName.endswith('vtp'):
        writer = vtk.vtkXMLPolyDataWriter()

    writer.SetFileName(a_FileName)
    writer.SetInputData(a_PolyData)
    writer.Update()
    writer.Write()

if __name__=="__main__":
# /Users/caoky/Downloads/vtkPathlineGeometry-new.py
    curveFile = '/mnt/c/MTH/317_HtoB/P317_CL.vtp'
    outFile = '/mnt/c/MTH/317_HtoB/P317_CL_BSpline.vtp'
    curve       = getPolyLineObject(curveFile)
    curve       = removeDuplicatePoints(curve)

    curveComp   = computeCurvatureTorsion(curve)
    writePolyLine(curveComp, outFile)


    # Interpolate into BSpline
    splineDegree = 2

    # Number of resampled points
    nResampledPoints = 5000

    # Make points on a_PolyData control points for the bspline
    # This is generally smoother, but we will get a slightly less accurate curve
    splineCurveCtrlPts = interpBSpline(curve, splineDegree, a_Mode='ctrlpts')
    splineCtrlPtsPolyData = convertBSplineToPolyData(splineCurveCtrlPts, nResampledPoints)
    # splineCtrlPtsComp = computeCurvatureTorsion(splineCtrlPtsPolyData)
    writePolyLine(splineCtrlPtsPolyData, outFile[:-4]+'_spline_ctrlpts.vtp')

    # Interpolate a spline through the points
    # This will give you a tighter fit, but you run the risk of overfitting
    # splineCurveInterp = interpBSpline(curve, splineDegree, a_Mode='interp')
    # splineInterpPolyData = convertBSplineToPolyData(splineCurveInterp, nResampledPoints)
    # # splineInterpComp = computeCurvatureTorsion(splineInterpPolyData)
    # writePolyLine(splineInterpPolyData, outFile[:-4]+"_spline_intrp.vtp")
