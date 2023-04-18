#----------------------------------------------------------------
# An implementation of the Euler Maruyama Method for solving a 
# tracer advection-diffusion equation to create a 
# Lagrangian Tracer Transport Model:
# dX = U(X,t)*dt + sqrt(2*D_0)*dB
# U(X,t): velocity fields obtained from external CFD calculations
#
# Required functionality:
# 1. Tracer intersection with spherical ensemble
# 2. Tracer intersection with superquadric ensemble
# 3. Tracer release from parametrized superquadric surfaces
# 4. Tracer intersection with Smoluchowski Radius 
#----------------------------------------------------------------
import sys, os
import numpy as np
import numpy.random as nprand
# MODIFY THIS INTO A BASE MODULE PROVIDING ESSENTIAL FUNCTIONALITIES
from tracerInputModule import *
from simpleTracersNew import *

import InOut

#------------------------------------------------------------------------------------
# for the collisional version, we need a separate array to track contact interactions
#------------------------------------------------------------------------------------
def initializeTracerTypesArray(a_N):
    
    tracerTypes = vtk.vtkIdList()
    tracerTypes.SetNumberOfIds(a_N)

    for p in range(a_N):
        tracerTypes.SetId(p,1)

    return tracerTypes

#-------------------------------------------------------------------------------------------
# This is an inefficient implementation of the OBB check for the collection of superquadrics
# The actual implementation has to be using f2py to make the calculations efficient
# This is a test module only
#-------------------------------------------------------------------------------------------
def ensembleOBBCheck(a_PNum, a_Dim, a_Matrix, a_Geom, a_X):

    if a_Dim == 2:
        x1 = a_X[0]
        x2 = a_X[1]
        x3 = 0.0
    elif a_Dim == 3:
        x1 = a_X[0]
        x2 = a_X[1]
        x3 = a_X[2]

    OBBInOut = -1

    for p in range(a_PNum):
        
        a1  = a_Geom[p,0]
        a2  = a_Geom[p,1]
        a3  = a_Geom[p,2]
        e1  = a_Geom[p,3]
        e2  = a_Geom[p,4]

        x1_BB0  = -0.5*a1
        x1_BB1  = 0.5*a1
        x2_BB0  = -0.5*a2
        x2_BB1  = 0.5*a2
        x3_BB0  = -0.5*a3
        x3_BB1  = 0.5*a3

        Nx  = a_Matrix[p,0]
        Ny  = a_Matrix[p,1]
        Nz  = a_Matrix[p,2]
        Ox  = a_Matrix[p,3]
        Oy  = a_Matrix[p,4]
        Oz  = a_Matrix[p,5]
        Ax  = a_Matrix[p,6]
        Ay  = a_Matrix[p,7]
        Az  = a_Matrix[p,8]
        Px  = a_Matrix[p,9]
        Py  = a_Matrix[p,10]
        Pz  = a_Matrix[p,11]

        x1_s  = Nx*x1 + Ny*x2 + Nz*x3 - (Px*Nx + Py*Ny + Pz*Nz)
        x2_s  = Ox*x1 + Oy*x2 + Oz*x3 - (Px*Ox + Py*Oy + Pz*Oz)
        x3_s  = Ax*x1 + Ay*x2 + Az*x3 - (Px*Ax + Py*Oy + Pz*Oz)
        
        numerator1  = Nx*x1 + Ny*x2 + Nz*x3 - Px*Nx - Py*Ny - Pz*Nz
        numerator2  = Ox*x1 + Oy*x2 + Oz*x3 - Px*Ox - Py*Oy - Pz*Oz
        numerator3  = Ax*x1 + Ay*x2 + Az*x3 - Px*Ax - Py*Ay - Pz*Az
        
        FInOut  = ( np.abs(numerator1/a1)**(2.0/e2) + np.abs(numerator2/a2)**(2.0/e2) )**(e2/e1) + np.abs(numerator3/a3)**(2.0/e1) 
        check   = (x1_s >= x1_BB0) and (x1_s <= x1_BB1) and (x2_s >= x2_BB0) and (x2_s <= x2_BB1) and (x3_s >= x3_BB0) and (x3_s <= x3_BB1)

        if FInOut < 1.0:
            OBBInOut = p
            break
    
    return OBBInOut

#-----------------------------------------------------------------------------------------
# Function that simply checks for whether a tracer is inside or outside a discrete element
#-----------------------------------------------------------------------------------------
def tracerDiscreteElementSimpleIntersection(a_Points, a_Types, a_Dim, a_Matrix, a_Geom):
    
    for i in range(a_Points.GetNumberOfPoints()):
    
        X = a_Points.GetPoint(i)

        if a_Dim == 2:
            x1 = X[0]
            x2 = X[1]
            x3 = 0.0
        elif a_Dim == 3:
            x1 = X[0]
            x2 = X[1]
            x3 = X[2]
        
        marker = 0

        for p in range(a_Matrix.shape[0]):
        
            a1  = a_Geom[p,0]
            a2  = a_Geom[p,1]
            a3  = a_Geom[p,2]
            e1  = a_Geom[p,3]
            e2  = a_Geom[p,4]
        
            Nx  = a_Matrix[p,0]
            Ny  = a_Matrix[p,1]
            Nz  = a_Matrix[p,2]
            Ox  = a_Matrix[p,3]
            Oy  = a_Matrix[p,4]
            Oz  = a_Matrix[p,5]
            Ax  = a_Matrix[p,6]
            Ay  = a_Matrix[p,7]
            Az  = a_Matrix[p,8]
        
            Px  = a_Matrix[p,9]
            Py  = a_Matrix[p,10]
            Pz  = a_Matrix[p,11]
        
            x1_s    = Nx*x1 + Ny*x2 + Nz*x3 - (Px*Nx + Py*Ny + Pz*Nz)
            x2_s    = Ox*x1 + Oy*x2 + Oz*x3 - (Px*Ox + Py*Oy + Pz*Oz)
            x3_s    = Ax*x1 + Ay*x2 + Az*x3 - (Px*Ax + Py*Oy + Pz*Oz)

            InOut   = ( abs(x1_s/a1)**(2.0/e2) + abs(x2_s/a2)**(2.0/e2) )**(e2/e1) + abs(x3_s/a3)**(2.0/e1)

            if InOut <= 1.0:
                marker = marker + 1
    
        if marker > 0:
            a_Types.SetId(i,0)
        else:
            a_Types.SetId(i,1)

    return [a_Points, a_Types]

#---------------------------------------------------------------------------------------------
# A more efficient implementation of the function to check whether a tracer particle is inside 
# a discrete element in an aggregate, using f2py wrapped fortran code for the contact check
#---------------------------------------------------------------------------------------------
def tracerDiscreteElementSimpleIntersectionWrapper(a_Points, a_Types, a_Dim, a_Geom, a_Matrix=None):

    for i in range(a_Points.GetNumberOfPoints()):

        X   = a_Points.GetPoint(i)
        Xp  = np.zeros(a_Dim)
        
        if a_Dim == 2:
            Xp[0] = X[0]
            Xp[1] = X[1]
        elif a_Dim == 3:
            Xp[0] = X[0]
            Xp[1] = X[1]
            Xp[2] = X[2]
        
        if a_Matrix is not None:
            value = InOut.generalisedsuperquadricensembleinout(a_Geom.shape[0], a_Dim, a_Matrix, a_Geom, Xp)
        else:
            value = InOut.ensembleinout(a_Geom.shape[0], a_Dim, a_Geom, Xp)

        if value == 1:
            a_Types.SetId(i,0)
        else:
            a_Types.SetId(i,1)

    return a_Types

#-------------------------------------------------------------------------
# Function to inject points at locations along the surface of discrete 
# elements in  a discrete element aggregate 
#-------------------------------------------------------------------------
def injectDiscreteElementTracers(a_Points, a_InputPoints, a_Dim, a_Matrix, a_Geom):

    return 0

#------------------------------------------------------
# Work in progress function for performing more precise 
# tracer-discrete element interaction calculations
#------------------------------------------------------
def tracerDiscreteElementInteraction(a_X, a_Dim, a_Matrix, a_Geom):

    if a_Dim == 2:
        x1 = a_X[0]
        x2 = a_X[1]
        x3 = 0.0
    elif a_Dim == 3:
        x1 = a_X[0]
        x2 = a_X[1]
        x3 = a_X[2]
    marker = 0

    for p in range(a_Matrix.shape[0]):
        
        a1  = a_Geom[p,0]
        a2  = a_Geom[p,1]
        a3  = a_Geom[p,2]
        e1  = a_Geom[p,3]
        e2  = a_Geom[p,4]
        
        Nx  = a_Matrix[p,0]
        Ny  = a_Matrix[p,1]
        Nz  = a_Matrix[p,2]
        Ox  = a_Matrix[p,3]
        Oy  = a_Matrix[p,4]
        Oz  = a_Matrix[p,5]
        Ax  = a_Matrix[p,6]
        Ay  = a_Matrix[p,7]
        Az  = a_Matrix[p,8]
        
        Px  = a_Matrix[p,9]
        Py  = a_Matrix[p,10]
        Pz  = a_Matrix[p,11]
        
        x1_s    = Nx*x1 + Ny*x2 + Nz*x3 - (Px*Nx + Py*Ny + Pz*Nz)
        x2_s    = Ox*x1 + Oy*x2 + Oz*x3 - (Px*Ox + Py*Oy + Pz*Oz)
        x3_s    = Ax*x1 + Ay*x2 + Az*x3 - (Px*Ax + Py*Oy + Pz*Oz)

    return nVec 

def readParticleVTKFile(a_ParticleVTK, a_SizeScalarName='scalars', a_Dim=2):

    try:
        import vtk
    except ImportError:
        print "Could Not Import Module vtk"

    try:
        if a_ParticleVTK.endswith('vtk'):
            reader = vtk.vtkPolyDataReader()
        elif a_ParticleVTK.endswith('vtp'):
            reader = vtk.vtkXMLPolyDataReader()

        reader.SetFileName(a_ParticleVTK)
        reader.Update()
    except IOError:
        print "Could Not Read File a_ParticleVTK"
        exit("Check file format and/or file extension!")

    data        = reader.GetOutput()
    numPoints   = data.GetNumberOfPoints()
    sizes       = data.GetPointData().GetArray(a_SizeScalarName)
    
    if a_Dim == 2:
        particleData = np.zeros((numPoints,3))
    elif a_Dim == 3:
        particleData = np.zeros((numPoints,4))

    for p in range(numPoints):

        xyzP = data.GetPoint(p)
        radP = sizes.GetValue(p)

        if a_Dim == 2:

            particleData[p,0] = xyzP[0]
            particleData[p,1] = xyzP[1]
            particleData[p,2] = radP

        elif a_Dim == 3:

            particleData[p,0] = xyzP[0]
            particleData[p,1] = xyzP[1] 
            particleData[p,2] = xyzP[2]
            particleData[p,3] = radP

    return particleData

#---------------------------------------------------------------
# NOTE: This version is using no vectorised code, mainly because 
# it is trying to access the lists using a VTK Object
#---------------------------------------------------------------
def eulerMaruyamaIntegrator(a_Particles, a_Indicators, a_Locator, a_Velocities, a_GridData, a_T, a_DT, a_D0, a_Window, a_BoundaryCondition):

    for p in range(a_Particles.GetNumberOfPoints()):

        if a_Indicators.GetId(p) == 1:
        
            xyz         = a_Particles.GetPoint(p)
            cell        = a_Locator.FindCell(xyz)
            cellPtIds   = vtk.vtkIdList()

            if len(a_Velocities) == 2:
            
                if cell != -1:

                    a_GridData.GetCellPoints(cell, cellPtIds)

                    velMinus_N1 = np.asarray(a_Velocities[0].GetTuple3(cellPtIds.GetId(0)))
                    velMinus_N2 = np.asarray(a_Velocities[0].GetTuple3(cellPtIds.GetId(1)))
                    velMinus_N3 = np.asarray(a_Velocities[0].GetTuple3(cellPtIds.GetId(2)))
                    velMinus    = (1.0/3.0)*(velMinus_N1 + velMinus_N2 + velMinus_N3)

                    velPlus_N1  = np.asarray(a_Velocities[1].GetTuple3(cellPtIds.GetId(0)))
                    velPlus_N2  = np.asarray(a_Velocities[1].GetTuple3(cellPtIds.GetId(1)))
                    velPlus_N3  = np.asarray(a_Velocities[1].GetTuple3(cellPtIds.GetId(2)))
                    velPlus     = (1.0/3.0)*(velPlus_N1 + velPlus_N2 + velPlus_N3)
                
                else:
                
                    velMinus    = np.array([0.0,0.0,0.0])
                    velPlus     = np.array([0.0,0.0,0.0])

                vel   = velMinus + (a_T - a_Window[0])*(velPlus - velMinus)/(a_Window[1] - a_Window[0])
            
            elif len(a_Velocities) == 1:

                if cell != -1:
                    a_GridData.GetCellPoints(cell, cellPtIds)

                    vel_N1 = np.asarray(a_Velocities[0].GetTuple3(cellPtIds.GetId(0)))
                    vel_N2 = np.asarray(a_Velocities[0].GetTuple3(cellPtIds.GetId(1)))
                    vel_N3 = np.asarray(a_Velocities[0].GetTuple3(cellPtIds.GetId(2)))
                    vel    = (1.0/3.0)*(vel_N1 + vel_N2 + vel_N3)
                
                else:
                    vel     = np.array([0.0,0.0,0.0])

            xNew    = xyz[0] + vel[0]*a_DT + ((2.0*a_D0*a_DT)**0.5)*nprand.randn()
            yNew    = xyz[1] + vel[1]*a_DT + np.sqrt(2.0*a_D0*a_DT)*nprand.randn()

        else:
            
            xyz     = a_Particles.GetPoint(p)
            xNew    = xyz[0]
            yNew    = xyz[1]
            
        #### AD-HOC
        #### Replace with a boundary condition function
            
        #if xNew > 0.060: xNew = 0.061
        
        #if xNew < -0.060: xNew = -0.061
        
        #if yNew > 0.015: yNew = 0.016
        
        #if yNew < -0.015: yNew = -0.016

        if xNew > 45.0: xNew = 45.01
        
        if xNew < 0.0: xNew = -0.01
        
        if yNew > 6.0: yNew = 6.01
        
        if yNew < 0.0: yNew = -0.01

        a_Particles.SetPoint(p, xNew, yNew, xyz[2])

    return a_Particles

if __name__=="__main__":

    #----------------------------
    # parse command line argument
    #----------------------------
    if len(sys.argv) != 2:
        sys.exit("Need Input Filename As An Argument")

    inputFile = sys.argv[1].strip()

    #--------------------------
    # set-up problem input data
    #--------------------------
    inputData = SimInputs(inputFile)
    inputData.readInputFile()
    
    particleEnsembleFile    = inputData.getDiscreteElementEnsemble()
    particleTransformFile   = inputData.getDiscreteElementTransforms()
    
    if particleEnsembleFile.endswith('vtk'):
        particleEnsemble    = readParticleVTKFile(particleEnsembleFile, a_SizeScalarName='radius')
    else:
        particleEnsemble    = np.loadtxt(particleEnsembleFile, skiprows=0)
    
    numParticles        = particleEnsemble.shape[0]

    if particleTransformFile is not None:
        transformData       = np.loadtxt(particleTransformFile, skiprows=1)
        particleGeometries  = particleEnsemble[:,6:]
    
    #----------------------------------------------------------
    # initialize tracer data from an external VTK polydata file
    #----------------------------------------------------------
    tracerPoints = initializeTracerFromFile(inputData.getTracerInput())
    tracerInject = initializeTracerFromFile(inputData.getTracerInput())

    #-----------------------------------------------------------
    # for collision check purposes, we create a contact ID array
    #-----------------------------------------------------------
    tracerTypes = initializeTracerTypesArray(tracerPoints.GetNumberOfPoints())

    #---------------------------------------------------
    # set up a tracer data object for writing into files
    #---------------------------------------------------
    tracerOutput    = vtk.vtkPolyData()
    tracerWriter    = vtk.vtkXMLPolyDataWriter()

    #---------------------------------------------------------------------------
    # generate the time synchronization map for flow data and integration points
    #---------------------------------------------------------------------------
    timeWindowDict = inputData.getDataTimeWindows()

    #-------------------
    # start time counter
    #-------------------
    timeIndex   = 0
    simTime     = inputData.getSimulationStartTime()
    tWin_0      = timeWindowDict['T_Low'][0]
    tWin_1      = timeWindowDict['T_Up'][0]

    if inputData.getDumpInterval() > 0:
        dumpCounter = 0

    while simTime <= inputData.getSimulationStopTime():

        print "Integrating from {0:6d} to {1:6d} simulation time {2:6.4f}".format(timeIndex, timeIndex + 1, simTime)

        #-----------------------------------------------------------------------
        # create a boolean condition that governs when new data files are loaded
        #-----------------------------------------------------------------------
        if inputData.isSteadyFlowData():

            isLoadFrame     = (simTime == inputData.getSimulationStartTime())
            isDataPointSync = True

        else:

            isLoadFrame = (simTime == inputData.getSimulationStartTime()) \
                        or (timeWindowDict['T_Low'][timeIndex] != tWin_0 or timeWindowDict['T_Up'][timeIndex] != tWin_1)
            
            isDataPointSync = timeWindowDict['T_Low'][timeIndex] == timeWindowDict['T_Up'][timeIndex]
        
        print "Here"

        #---------------------------------------------------------------------------
        # create a boolean condition that governs when new points are to be injected
        #---------------------------------------------------------------------------
        isInjectPoints = (timeIndex%10 == 0) #False  #isLoadFrame    ### THIS NEEDS MORE POLISHING

        #-------------------------------------------
        # a set of debug messages for program status
        #-------------------------------------------
        if isLoadFrame: 
            print "Will Load Velocity Data"

        if isDataPointSync:
            if inputData.isSteadyFlowData():
                if simTime == inputData.getSimulationStartTime():
                    print "Data And Integration Times Are Synced"
            else:
                print "Data And Integration Times Are Synced"

        if isInjectPoints: 
            print "New Particles Injected"
            print "Now Integrating", tracerPoints.GetNumberOfPoints(), "Particles"

        #-----------------------------------------
        # inject tracers into the domain if needed
        #-----------------------------------------
        if isInjectPoints:
            tracerPoints_New    = injectPoints(tracerPoints, tracerInject) ### THIS NEEDS MORE WORK
            tracerPoints        = tracerPoints_New
        
        print "Here 2"

        #-------------------------------------------------------
        # load data file based on the evaluated boolean variable
        #-------------------------------------------------------
        if isLoadFrame:

            if inputData.isSteadyFlowData():

                flowSingle      = inputData.getFlowDataFileName()
                velocitySingle  = extractFlowDataFromFile(flowSingle, a_DataName=inputData.getVelDataName())
                gridDataObject  = extractDataFromFile(flowSingle)

                print "Here 2.1"

            else:
                
                tWin_0  = timeWindowDict['T_Low'][timeIndex]    # updating the time window
                tWin_1  = timeWindowDict['T_Up'][timeIndex]     # updating the time window

                print "Times:", tWin_1, tWin_0

                if isDataPointSync:
                    print "Here 2.2"
                    flowSingle      = inputData.getFlowDataFileName(a_ID = timeWindowDict['ID_Low'][timeIndex])
                    print flowSingle
                    velocitySingle  = extractFlowDataFromFile(flowSingle, a_DataName=inputData.getVelDataName())
                    gridDataObject  = extractDataFromFile(flowSingle)        
                else:
                    print "Here 2.3"
                    flowPlus        = inputData.getFlowDataFileName(a_ID = timeWindowDict['ID_Up'][timeIndex])
                    flowMinus       = inputData.getFlowDataFileName(a_ID = timeWindowDict['ID_Low'][timeIndex])
                    print flowPlus
                    print flowMinus
                    print inputData.getVelDataName()
                    velocityPlus    = extractFlowDataFromFile(flowPlus,  a_DataName=inputData.getVelDataName())
                    print "Read plus"
                    velocityMinus   = extractFlowDataFromFile(flowMinus, a_DataName=inputData.getVelDataName())
                    print "Read minus"
                    gridDataObject  = extractDataFromFile(flowMinus)
        
        print "Here 3"

        #--------------------------------------------------------------
        # build a cell locator (if it is a fixed mesh, built only once)
        #--------------------------------------------------------------
        if simTime == inputData.getSimulationStartTime() and inputData.isFixedMesh():

            print "Building Cell Locator Maps"

            if isDataPointSync:
                locatorObj  = createCellLocator(flowSingle)
            else:
                locatorObj  = createCellLocator(flowMinus)
        
        #--------------------------------------------------
        # now proceed with integration of each tracer point
        #--------------------------------------------------
        if isDataPointSync:
            velocities = [velocitySingle]
        else:
            velocities = [velocityMinus, velocityPlus]
        
        print "Here pre-integration"
        
        boundaryCondition = 1
        tracerPoints = eulerMaruyamaIntegrator(tracerPoints, tracerTypes, locatorObj, velocities, gridDataObject, simTime,
                inputData.getIntegrationTimeStep(), 1.0e-5, [tWin_0, tWin_1], boundaryCondition)
        
        #--------------------------------------------------------------
        # update the tracer inidcators/types/ids based on whether they 
        # interacted with a discrete element or not
        #--------------------------------------------------------------
        #[tracerPoints, tracerTypes] = tracerDiscreteElementSimpleIntersection(tracerPoints, tracerTypes, 2, transformData, particleGeometries)
        if particleTransformFile is not None:
            tracerTypes = tracerDiscreteElementSimpleIntersectionWrapper(tracerPoints, tracerTypes, 2, particleGeometries, a_Matrix=transformData)
        else:
            tracerTypes = tracerDiscreteElementSimpleIntersectionWrapper(tracerPoints, tracerTypes, 2, particleEnsemble)
        
        #--------------------------------------------------------------
        # at the end of appropriate number of steps dump data into file
        #--------------------------------------------------------------
        if inputData.getDumpInterval() == 0:
            writeTracerDataToFile(tracerPoints, inputData.getTracerOutputFile(a_ID=timeIndex))
        else:
            if timeIndex%inputData.getDumpInterval() == 0:
                dumpCounter = dumpCounter + 1
                writeTracerDataToFile(tracerPoints, inputData.getTracerOutputFile(a_ID=dumpCounter))
        
        #----------------------------------------
        # update time indices and simulation time
        #----------------------------------------
        timeIndex   = timeIndex + 1
        simTime     = simTime + inputData.getIntegrationTimeStep()
