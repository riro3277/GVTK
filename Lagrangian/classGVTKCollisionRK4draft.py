#-------------------------------------------------------------------------------
# This module computes the embolus dynamics using a one way coupling scheme based
# on signed distance field algorithm for collision.

# edit version  : does not include the restart code

#
# Author:      Akshita Sahni
# Institution: University of Colorado, Boulder
# Last Edit:   Sept. 2021
#-------------------------------------------------------------------------------
import sys, os, vtk
import numpy as np

try:
    from classGVTKGrid import *
    from classGVTKLagrangianData import *
    from classGVTKGenericProblem import *
    from moduleVTKLocators import *
except:
    sys.exit('Could not import user defined modules')

class GVTKCollision(GVTKProblem):

    def __init__(self, a_InputData, a_DebugMode=True):

        super().__init__(a_InputData,a_DebugMode)
        self.m_LagrangianData = None
        self.m_GridData       = None

    @property
    def numParticles(self):
        return self.m_LagrangianData.numParticles

    @property
    def numCells(self):
        return self.m_GridData.numCells

    def runCompute(self):
 		#-----------------------------------------------------------------------
        # initialize lagrangian integration from an external VTK polydata file and
        # update the Lagrangian inertial object for the problem
        #-----------------------------------------------------------------------
        self.m_LagrangianData = GVTKLagrangianData(a_File=self.m_InputData.getTracerInput())
        #-----------------------------------------------------------------------
        # initialize blank grid data object as member of physics object
        # Note: The object gets instantiated here, but the actual data manipulation
        # happens at instances when flow or grid data files are read.
        #-----------------------------------------------------------------------
        self.m_GridData = GVTKGridData()

        #-----------------------------------------------------------------------
        # load particles from the same spot where particles were released
        #-----------------------------------------------------------------------
        if self.m_InputData.isInjectParticles():
            tracerInject = initializeTracerFromFile(self.m_InputData.getTracerInput())

        #---------------------------------------------------------------------------
        # generate the time synchronization map for flow data and integration points
        #---------------------------------------------------------------------------
        print('changed the time window version')
        timeWindowDict = self.m_InputData.getDataTimeWindows_v2()
        # timeWindowDict = self.m_InputData.getDataTimeWindows()

        #-----------------------------------------------------------------------------
        # display the time windows in case the integration involves unsteady flow data
        #-----------------------------------------------------------------------------
        print('print data time windows')
        if not(self.m_InputData.isSteadyFlowData()):
            self.m_InputData.printDataTimeWindows()

        #-----------------------------------------------------------------------
        # start time counter
        #-----------------------------------------------------------------------
        timeIndex = 0
        simTime   = self.m_InputData.getSimulationStartTime()
        stopTime  = self.m_InputData.getSimulationStopTime()
        dT        = self.m_InputData.getIntegrationTimeStep()
        tWin_0    = timeWindowDict['T_Low'][0] # t_low[0] usually 0
        tWin_1    = timeWindowDict['T_Up'][0] # t_up[0] usually 0
        isSingleStageIntegration = False # set true for feu and false for for rk4

        # get signed distance field
        ## TODO: move this to the input.dat file
        meshFile    = '/Volumes/easystore/New-LVAD-dataset-Post-processed-files/Case-ID-11-ConstantVAD-3rdcycle-post/volume-mesh-wSDF-igl.vtu'
        surfaceFile = '/Volumes/easystore/New-LVAD-dataset-Post-processed-files/Case-ID-11-ConstantVAD-3rdcycle-post/clean-Case-ID-11_LVAD-surface.stl'
        outFile     = '/Volumes/easystore/New-LVAD-dataset-Post-processed-files/Case-ID-11-ConstantVAD-3rdcycle-post/Case11-r0pt25mm/New-Case-ID-11_LVAD-sdf-grad.vtu'
        self.vtkPolydataDistanceField(meshFile, surfaceFile, outFile)

        # added gradient filter application on existing grid data
        self.m_GridData.getDataFromFile(outFile, a_FileType ='vtu')

        #-----------------------------------------------------------------------
        # start the simulation time loop
        #-----------------------------------------------------------------------
        while simTime <= stopTime:

            #--------------------------------------------------------------------
            # STEP 0: display a status message for the time step being integrated
            #--------------------------------------------------------------------
            if simTime != 0:
                print("Integrating from", timeIndex, "to", timeIndex + 1, "simTime", int(simTime), end='\r')
            elif round(simTime, 2) == stopTime-dT:
                print("Integrating from", timeIndex, "to", timeIndex + 1, "simTime", int(simTime))

            #-------------------------------------------------------------------------------
            # STEP 1: create a boolean condition that governs when new data files are loaded
            #-------------------------------------------------------------------------------
            if self.m_InputData.isSteadyFlowData():

                isLoadFrame     = (simTime == self.m_InputData.getSimulationStartTime())
                isDataPointSync = True

            else:

                isLoadFrame = (simTime == self.m_InputData.getSimulationStartTime()) \
                        or (timeWindowDict['T_Low'][timeIndex] != tWin_0 or timeWindowDict['T_Up'][timeIndex] != tWin_1)

                isDataPointSync = timeWindowDict['T_Low'][timeIndex] == timeWindowDict['T_Up'][timeIndex]

            #-------------------------------------------------------------------
            # a set of debug messages for program status
            #-------------------------------------------------------------------
            if isLoadFrame:
                print("Will Load Velocity Data")

            if isDataPointSync:
                if self.m_InputData.isSteadyFlowData():
                    if simTime == self.m_InputData.getSimulationStartTime():
                        print("Data And Integration Times Are Synced")
                else:
                    print("Data And Integration Times Are Synced")

            #-------------------------------------------------------------------
            # STEP 3: resolve injections and inject particles
            # for this simple tracer particle integration we are injecting particles
            # either always or never.
            # TODO: THIS NEEDS TO BE IMPROVED
            #-------------------------------------------------------------------
            isInjectPoints = False

            #-------------------------------------------------------------------
            # inject particles into the domain if needed
            # currently it injects every time-step
            #-------------------------------------------------------------------
            if isInjectPoints:
                tracers.injectParticles(tracerInject)
                print("New Particles Injected")
                print("Now Integrating", tracers.numParticles, "Particles")

            #---------------------------------------------------------------
            # STEP 4: load data file based on the evaluated boolean variable
            #---------------------------------------------------------------
            if isLoadFrame:

                if self.m_InputData.isSteadyFlowData():

                    flowSingle = self.m_InputData.getFlowDataFileName() # .vtu data
                    self.m_GridData.addDataArrayFromFile(flowSingle, 'v0',self.m_InputData.getVelDataName())

                else:

                    #-----------------------------------------------------------------------------------
                    # TODO: For the unsteady file read write, the data array updates have to be modified
                    # using the add and remove data array functionalities for the GVTKGridData object
                    #-----------------------------------------------------------------------------------
                    tWin_0  = timeWindowDict['T_Low'][timeIndex]    # updating the time window
                    tWin_1  = timeWindowDict['T_Up'][timeIndex]     # updating the time window

                    if simTime == self.m_InputData.getSimulationStartTime():

                        flowInitial = self.m_InputData.getFlowDataFileName(a_ID = timeWindowDict['ID_Low'][timeIndex])
                        self.m_GridData.addDataArrayFromFile(flowInitial, 'v0',self.m_InputData.getVelDataName())

                    elif isDataPointSync:

                        if isSingleStageIntegration == True:

                            self.m_GridData.removeDataArray('v0')
                            flowSingle = self.m_InputData.getFlowDataFileName(a_ID = timeWindowDict['ID_Low'][timeIndex])
                            self.m_GridData.addDataArrayFromFile(flowSingle, 'v0', self.m_InputData.getVelDataName())

                        else:

                            self.m_GridData.removeDataArray('v0')
                            self.m_GridData.removeDataArray('v1')
                            flowPlus  = self.m_InputData.getFlowDataFileName(a_ID = timeWindowDict['ID_Up'][timeIndex])
                            flowMinus = self.m_InputData.getFlowDataFileName(a_ID = timeWindowDict['ID_Low'][timeIndex])
                            self.m_GridData.addDataArrayFromFile(flowPlus, 'v1', self.m_InputData.getVelDataName())
                            self.m_GridData.addDataArrayFromFile(flowMinus, 'v0', self.m_InputData.getVelDataName())
                            print('Loaded a new velocity data')

                    else:

                        if isSingleStageIntegration == True:

                            self.m_GridData.removeDataArray('v0')
                            self.m_GridData.removeDataArray('v1')
                            flowPlus  = self.m_InputData.getFlowDataFileName(a_ID = timeWindowDict['ID_Up'][timeIndex])
                            flowMinus = self.m_InputData.getFlowDataFileName(a_ID = timeWindowDict['ID_Low'][timeIndex])
                            self.m_GridData.addDataArrayFromFile(flowPlus, 'v1', self.m_InputData.getVelDataName())
                            self.m_GridData.addDataArrayFromFile(flowMinus, 'v0', self.m_InputData.getVelDataName())
                            print('Loaded a new velocity data')

                        else:

                            if simTime < tWin_1 - dT:

                                self.m_GridData.removeDataArray('v0')
                                self.m_GridData.removeDataArray('v1')
                                self.m_GridData.removeDataArray('v2')
                                flowPlus  = self.m_InputData.getFlowDataFileName(a_ID = timeWindowDict['ID_Up'][timeIndex])
                                flowMinus = self.m_InputData.getFlowDataFileName(a_ID = timeWindowDict['ID_Low'][timeIndex])
                                self.m_GridData.addDataArrayFromFile(flowPlus, 'v1', self.m_InputData.getVelDataName())
                                self.m_GridData.addDataArrayFromFile(flowMinus, 'v0', self.m_InputData.getVelDataName())
                                print('Loaded a new velocity data')

                            elif simTime >= tWin_1 - dT:

                                self.m_GridData.removeDataArray('v0')
                                self.m_GridData.removeDataArray('v1')
                                self.m_GridData.removeDataArray('v2')
                                flowPlus  = self.m_InputData.getFlowDataFileName(a_ID = timeWindowDict['ID_Up'][timeIndex])
                                flowMinus = self.m_InputData.getFlowDataFileName(a_ID = timeWindowDict['ID_Low'][timeIndex])
                                flowNext  = self.m_InputData.getFlowDataFileName(a_ID = timeWindowDict['ID_Low'][timeIndex]+1)
                                self.m_GridData.addDataArrayFromFile(flowPlus, 'v1', self.m_InputData.getVelDataName())
                                self.m_GridData.addDataArrayFromFile(flowMinus, 'v0', self.m_InputData.getVelDataName())
                                self.m_GridData.addDataArrayFromFile(flowNext, 'v2', self.m_InputData.getVelDataName())
                                print('Loaded a new velocity data')

            #-----------------------------------------------------------------------------------------------
            # STEP 5: build a cell locator (if it is a fixed mesh, built only once) and set initial velocity
            #-----------------------------------------------------------------------------------------------
            if simTime == self.m_InputData.getSimulationStartTime() and self.m_InputData.isFixedMesh():
                print("Building Cell Locator Maps")
                if isDataPointSync:
                    locatorObj  = createCellLocatorFromData(self.m_GridData.m_vtkData, a_LocatorType = self.m_InputData.getLocatorType())# syntax for standard vtkLocator
                else:
                    locatorObj  = createCellLocatorFromData(self.m_GridData.m_vtkData, a_LocatorType = self.m_InputData.getLocatorType())# syntax for standard vtkLocator
                self.m_GridData.buildLocator(a_LocatorType = self.m_InputData.getLocatorType())

            #------------------------------------------
            # set initial velocity in m_LagrangianData
            #------------------------------------------
            #-----------------------------------
            # set initial velocity of particles
            #-----------------------------------
            InitialVelocity = np.zeros((self.numParticles, 3))
            InitialVelocity[:,0] = float(self.m_InputData.m_InitialVelocity['u'])
            InitialVelocity[:,1] = float(self.m_InputData.m_InitialVelocity['v'])
            InitialVelocity[:,2] = float(self.m_InputData.m_InitialVelocity['w'])

            if simTime == self.m_InputData.getSimulationStartTime():
                self.m_LagrangianData.addVectorData('Velocity', InitialVelocity)

            #----------------------------------------------------------------
            # STEP 7: now proceed with integration of each inertial particle
            #----------------------------------------------------------------
            boundaryCondition = 1   # currently a dummy variable

            if intScheme == feu:
                self.sdfCollision(simTime, tWin_0, tWin_1, boundaryCondition, locatorObj, a_DataSync=isDataPointSync, a_PolygonalCells=False)
            elif intScheme == rk4:
                self.sdfCollision_RK4(simTime, tWin_0, tWin_1, t_Next, boundaryCondition, locatorObj, a_DataSync=isDataPointSync, a_PolygonalCells=False)

            #----------------------------------------------------------------------
            # STEP 8: at the end of appropriate number of steps dump data into file
            #----------------------------------------------------------------------
            if timeIndex%self.m_InputData.m_WriteInterval == 0:
                self.m_LagrangianData.writeData(self.m_InputData.getTracerOutputFile(a_ID1=0, a_ID2=timeIndex))
            #------------------------------------------------
            # STEP 9: update time indices and simulation time
            #------------------------------------------------
            timeIndex = timeIndex + 1
            simTime   = simTime + dT

    ##TEST snippet
    ## add sdf rk4 here and check if the error is reproduced in the current version
    def sdfCollision_RK4(self, simTime, tWin_0, tWin_1, t_Next, boundaryCondition, a_Locator, a_DataSync=True, a_PolygonalCells=False):

        restitution     = 0.75
        flowVel   = self.m_GridData.extractDataArray('velocity')
        dist      = self.m_GridData.extractDataArray('sdf')
        normVec   = self.m_GridData.extractDataArray('sdf-grad')
        dT        = self.m_InputData.getIntegrationTimeStep()

        for p in range(self.numParticles):
            xyz_1 = self.m_LagrangianData.getX(p)
            v_i   = np.array(self.m_LagrangianData.getVectorData(a_ArrayName='Velocity', a_DataID=p))
            # print('getV = ', p, v_i)
            cell    = a_Locator.FindCell(xyz_1)
            cellPtIds = vtk.vtkIdList()
            tauP      = 4.0*self.m_InputData.getTracerDensity()*(self.m_InputData.m_TracerRadius**2)/(18.0*self.m_InputData.getFluidViscosity())

        #----------------------------
        # calculate average velocity
        #----------------------------
            if a_DataSync == True:

                if a_PolygonalCells == True:
                    vel_1, stat = self.m_GridData.gridInterpolateAveraging(xyz_1, 'v0', a_GetStatus=True)

                    if stat != -1:
                        self.m_GridData.m_vtkData.GetCellPoints(cell, cellPtIds)

                        dN1 = dist.GetTuple1(cellPtIds.GetId(0))
                        dN2 = dist.GetTuple1(cellPtIds.GetId(1))
                        dN3 = dist.GetTuple1(cellPtIds.GetId(2))
                        dN4 = dist.GetTuple1(cellPtIds.GetId(3))

                        gN1 = np.asarray(normVec.GetTuple3(cellPtIds.GetId(0)))
                        gN2 = np.asarray(normVec.GetTuple3(cellPtIds.GetId(1)))
                        gN3 = np.asarray(normVec.GetTuple3(cellPtIds.GetId(2)))
                        gN4 = np.asarray(normVec.GetTuple3(cellPtIds.GetId(3)))

                        d   = 0.25*(dN1 + dN2 + dN3 + dN4)
                        g   = 0.25*(gN1 + gN2 + gN3 + gN4)

                        # rk4 scheme
                        k_1   = v_i
                        el_1  = (1.0/tauP)*(vel_1-v_i)

                        xyz_2 = xyz_1 + 0.5*k_1*dT
                        vel_2 = self.m_GridData.gridInterpolateAveraging(xyz_2, 'v0')
                        k_2   = v_i + 0.5*el_1*dT
                        el_2  = (1.0/tauP)*(vel_2-v_i-0.5*el_1*dT)

                        xyz_3 = xyz_1 + 0.5*k_2*dT
                        vel_3 = self.m_GridData.gridInterpolateAveraging(xyz_3, 'v0')
                        k_3   = v_i + 0.5*el_2*dT
                        el_3  = (1.0/tauP)*(vel_3-v_i-0.5*el_2*dT)

                        xyz_4 = xyz_1 + k_3*dT
                        vel_4 = self.m_GridData.gridInterpolateAveraging(xyz_4, 'v0')
                        k_4   = v_i + el_3*dT
                        el_4  = (1.0/tauP)*(vel_4-v_i-el_3*dT)

                        # new particle position and velocity

                        xyz_new = xyz_1 + (dT/6)*(k_1  + 2*k_2  + 2*k_3  + k_4)
                        v_new   = v_i   + (dT/6)*(el_1 + 2*el_2 + 2*el_3 + el_4)

                        if d <= self.m_InputData.m_TracerRadius:
                            vn      = np.dot(v_new, g)
                            vt      = v_new - vn*g
                            v_new   = vt - restitution*g*vn

                    # else:
                    #     xyz_new = xyz_1
                    #     v_new   = np.zeros(self.m_InputData.m_SpaceDimension)

                else:
                    vel_1, stat = self.m_GridData.gridInterpolateNodalBasis(xyz_1, 'v0', a_GetStatus=True)

                    if stat != -1:
                        self.m_GridData.m_vtkData.GetCellPoints(cell, cellPtIds)

                        dN1 = dist.GetTuple1(cellPtIds.GetId(0))
                        dN2 = dist.GetTuple1(cellPtIds.GetId(1))
                        dN3 = dist.GetTuple1(cellPtIds.GetId(2))
                        dN4 = dist.GetTuple1(cellPtIds.GetId(3))

                        gN1 = np.asarray(normVec.GetTuple3(cellPtIds.GetId(0)))
                        gN2 = np.asarray(normVec.GetTuple3(cellPtIds.GetId(1)))
                        gN3 = np.asarray(normVec.GetTuple3(cellPtIds.GetId(2)))
                        gN4 = np.asarray(normVec.GetTuple3(cellPtIds.GetId(3)))

                        d   = 0.25*(dN1 + dN2 + dN3 + dN4)
                        g   = 0.25*(gN1 + gN2 + gN3 + gN4)

                        k_1   = v_i
                        el_1  = (1.0/tauP)*(vel_1-v_i)

                        xyz_2 = xyz_1 + 0.5*k_1*dT
                        vel_2 = self.m_GridData.gridInterpolateNodalBasis(xyz_2, 'v0')
                        k_2   = v_i + 0.5*el_1*dT
                        el_2  = (1.0/tauP)*(vel_2-v_i-0.5*el_1*dT)

                        xyz_3 = xyz_1 + 0.5*k_2*dT
                        vel_3 = self.m_GridData.gridInterpolateNodalBasis(xyz_3, 'v0')
                        k_3   = v_i + 0.5*el_2*dT
                        el_3  = (1.0/tauP)*(vel_3-v_i-0.5*el_2*dT)

                        xyz_4 = xyz_1 + k_3*dT
                        vel_4 = self.m_GridData.gridInterpolateNodalBasis(xyz_4, 'v0')
                        k_4   = v_i + el_3*dT
                        el_4  = (1.0/tauP)*(vel_4-v_i-el_3*dT)

                        # new particle position and velocity

                        xyz_new = xyz_1 + (dT/6)*(k_1  + 2*k_2  + 2*k_3  + k_4)
                        v_new   = v_i   + (dT/6)*(el_1 + 2*el_2 + 2*el_3 + el_4)

                        # else:
                        #     xyz_new = xyz_1
                        #     v_new   = np.zeros(self.m_InputData.m_SpaceDimension)

                        if d <= self.m_InputData.m_TracerRadius:
                            vn      = np.dot(v_new, g)
                            vt      = v_new - vn*g
                            v_new   = vt - restitution*g*vn

            else:
                t_1 = simTime
                t_2 = simTime + 0.5*dT
                t_3 = simTime + 0.5*dT
                t_4 = simTime + dT

                #---------
                # step 1
                #---------

                if a_PolygonalCells == True:
                    vel_1_L, stat = self.m_GridData.gridInterpolateAveraging(xyz_1, 'v0', a_GetStatus=True)
                    vel_1_U = self.m_GridData.gridInterpolateAveraging(xyz_1, 'v1')
                else:
                    vel_1_L, stat = self.m_GridData.gridInterpolateNodalBasis(xyz_1, 'v0', a_GetStatus=True)
                    vel_1_U = self.m_GridData.gridInterpolateNodalBasis(xyz_1, 'v1')

                #---------------------------
                # particles still suspended
                #---------------------------

                if stat!= -1:

                    self.m_GridData.m_vtkData.GetCellPoints(cell, cellPtIds)

                    dN1 = dist.GetTuple1(cellPtIds.GetId(0))
                    dN2 = dist.GetTuple1(cellPtIds.GetId(1))
                    dN3 = dist.GetTuple1(cellPtIds.GetId(2))
                    dN4 = dist.GetTuple1(cellPtIds.GetId(3))

                    gN1 = np.asarray(normVec.GetTuple3(cellPtIds.GetId(0)))
                    gN2 = np.asarray(normVec.GetTuple3(cellPtIds.GetId(1)))
                    gN3 = np.asarray(normVec.GetTuple3(cellPtIds.GetId(2)))
                    gN4 = np.asarray(normVec.GetTuple3(cellPtIds.GetId(3)))

                    d   = 0.25*(dN1 + dN2 + dN3 + dN4)
                    g   = 0.25*(gN1 + gN2 + gN3 + gN4)

                    vel_1 = vel_1_L + (simTime - tWin_0) * (vel_1_U - vel_1_L)/(tWin_1 - tWin_0)
                    k_1   = v_i
                    el_1  = (1.0/tauP)*(vel_1-v_i)

                    # if d <= self.m_InputData.m_TracerRadius:
                    #     vn      = np.dot(v_new, g)
                    #     vt      = v_new - vn*g
                    #     v_new   = vt - restitution*g*vn

                    # step 2
                    xyz_2 = xyz_1 + 0.50*k_1*dT

                    if t_2 <= tWin_1:
                        if a_PolygonalCells == True:
                            vel_2_L = self.m_GridData.gridInterpolateAveraging(xyz_2, 'v0')
                            vel_2_U = self.m_GridData.gridInterpolateAveraging(xyz_2, 'v1')
                        else:
                            vel_2_L = self.m_GridData.gridInterpolateNodalBasis(xyz_2, 'v0')
                            vel_2_U = self.m_GridData.gridInterpolateNodalBasis(xyz_2, 'v1')

                        vel_2 = vel_2_L + (simTime - tWin_0) * (vel_2_U - vel_2_L)/(tWin_1 - tWin_0)

                    elif t_2 > tWin_1:
                        if a_PolygonalCells == True:
                            vel_2_L = self.m_GridData.gridInterpolateAveraging(xyz_2, 'v1')
                            vel_2_U = self.m_GridData.gridInterpolateAveraging(xyz_2, 'v2')
                        else:
                            vel_2_L = self.m_GridData.gridInterpolateNodalBasis(xyz_2, 'v1')
                            vel_2_U = self.m_GridData.gridInterpolateNodalBasis(xyz_2, 'v2')

                        vel_2 = vel_2_L + (simTime - tWin_1) * (vel_2_U - vel_2_L)/(t_Next - tWin_1)

                    k_2   = v_i + 0.5*el_1*dT
                    el_2  = (1.0/tauP)*(vel_2-v_i-0.5*el_1*dT)


                    # step 3
                    xyz_3 = xyz_1 + 0.50*k_2*dT

                    if t_3 <= tWin_1:
                        if a_PolygonalCells == True:
                            vel_3_L = self.m_GridData.gridInterpolateAveraging(xyz_3, 'v0')
                            vel_3_U = self.m_GridData.gridInterpolateAveraging(xyz_3, 'v1')
                        else:
                            vel_3_L = self.m_GridData.gridInterpolateNodalBasis(xyz_3, 'v0')
                            vel_3_U = self.m_GridData.gridInterpolateNodalBasis(xyz_3, 'v1')

                        vel_3 = vel_3_L + (simTime - tWin_0) * (vel_3_U - vel_3_L)/(tWin_1 - tWin_0)

                    elif t_3 > tWin_1:
                        if a_PolygonalCells == True:
                            vel_3_L = self.m_GridData.gridInterpolateAveraging(xyz_3, 'v1')
                            vel_3_U = self.m_GridData.gridInterpolateAveraging(xyz_3, 'v2')
                        else:
                            vel_3_L = self.m_GridData.gridInterpolateNodalBasis(xyz_3, 'v1')
                            vel_3_U = self.m_GridData.gridInterpolateNodalBasis(xyz_3, 'v2')

                        vel_3 = vel_3_L + (simTime - tWin_1) * (vel_3_U - vel_3_L)/(t_Next - tWin_1)

                    k_3   = v_i + 0.5*el_2*dT
                    el_3  = (1.0/tauP)*(vel_3-v_i-0.5*el_2*dT)

                    # step 4
                    xyz_4 = xyz_1 + k_3*dT

                    if t_4 <= tWin_1:
                        if a_PolygonalCells == True:
                            vel_4_L = self.m_GridData.gridInterpolateAveraging(xyz_4, 'v0')
                            vel_4_U = self.m_GridData.gridInterpolateAveraging(xyz_4, 'v1')
                        else:
                            vel_4_L = self.m_GridData.gridInterpolateNodalBasis(xyz_4, 'v0')
                            vel_4_U = self.m_GridData.gridInterpolateNodalBasis(xyz_4, 'v1')

                        vel_4 = vel_4_L + (simTime - tWin_0) * (vel_4_U - vel_4_L)/(tWin_1 - tWin_0)

                    elif t_4 > tWin_1:
                        if a_PolygonalCells == True:
                            vel_4_L = self.m_GridData.gridInterpolateAveraging(xyz_4, 'v1')
                            vel_4_U = self.m_GridData.gridInterpolateAveraging(xyz_4, 'v2')
                        else:
                            vel_4_L = self.m_GridData.gridInterpolateNodalBasis(xyz_4, 'v1')
                            vel_4_U = self.m_GridData.gridInterpolateNodalBasis(xyz_4, 'v2')

                        vel_4 = vel_4_L + (simTime - tWin_1) * (vel_4_U - vel_4_L)/(t_Next - tWin_1)

                    k_4   = v_i + el_3*dT
                    el_4  = (1.0/tauP)*(vel_4-v_i-el_3*dT)

                    # new particle position and velocity

                    xyz_new = xyz_1 + (dT/6)*(k_1  + 2*k_2  + 2*k_3  + k_4)
                    v_new   = v_i   + (dT/6)*(el_1 + 2*el_2 + 2*el_3 + el_4)

                    if d <= self.m_InputData.m_TracerRadius:
                        vn      = np.dot(v_new, g)
                        vt      = v_new - vn*g
                        v_new   = vt - restitution*g*vn

        #--------------------------------
        # set new xyz for  point in loop
        #--------------------------------
            self.m_LagrangianData.setX(p, xyz_new)
            self.m_LagrangianData.setVectorData(v_new, a_ArrayName='Velocity', a_DataID=p)
            print('setV = ', p, v_new)
    ##TEST snippet

    def sdfCollision(self, simTime, tWin_0, tWin_1, boundaryCondition, a_Locator, a_DataSync=True, a_PolygonalCells=True):

        restitution = 0.75
        flowVel   = self.m_GridData.extractDataArray('velocity')
        dist      = self.m_GridData.extractDataArray('sdf')
        normVec   = self.m_GridData.extractDataArray('sdf-grad')
        dT        = self.m_InputData.getIntegrationTimeStep()

        #------------------------------------------
        # forward euler integration implementation
        #------------------------------------------
        for p in range(self.numParticles):
            posP  = self.m_LagrangianData.getX(p)
            v_i   = np.array(self.m_LagrangianData.getVectorData(a_ArrayName='Velocity', a_DataID=p))
            # print('getV = ', p, v_i)
            cell      = a_Locator.FindCell(posP)
            # cellPtIds = vtk.vtkIdList()
            tauP      = 4.0*self.m_InputData.getTracerDensity()*(self.m_InputData.m_TracerRadius**2)/(18.0*self.m_InputData.getFluidViscosity())

            if a_DataSync == True:

                if a_PolygonalCells == True:
                    v = self.m_GridData.gridInterpolateAveraging(posP, 'v0')
                else:
                    v = self.m_GridData.gridInterpolateNodalBasis(posP, 'v0')
            else:
                if a_PolygonalCells == True:
                    vPlus  = self.m_GridData.gridInterpolateAveraging(posP, 'v1')
                    vMinus = self.m_GridData.gridInterpolateAveraging(posP, 'v0')
                else:
                    vPlus  = self.m_GridData.gridInterpolateNodalBasis(posP, 'v1')
                    vMinus = self.m_GridData.gridInterpolateNodalBasis(posP, 'v0')

                v = vMinus + (simTime - tWin_0) * (vPlus - vMinus)/(tWin_1 - tWin_0)

            if cell !=-1:

                v_i   = v_i + (1.0/tauP)*dT*(v-v_i)
                posP = posP + v_i*dT

                cell      = a_Locator.FindCell(posP)
                cellPtIds = vtk.vtkIdList()

                if cell !=-1:

                    self.m_GridData.m_vtkData.GetCellPoints(cell, cellPtIds)

                    dN1 = dist.GetTuple1(cellPtIds.GetId(0))
                    dN2 = dist.GetTuple1(cellPtIds.GetId(1))
                    dN3 = dist.GetTuple1(cellPtIds.GetId(2))
                    dN4 = dist.GetTuple1(cellPtIds.GetId(3))

                    gN1 = np.asarray(normVec.GetTuple3(cellPtIds.GetId(0)))
                    gN2 = np.asarray(normVec.GetTuple3(cellPtIds.GetId(1)))
                    gN3 = np.asarray(normVec.GetTuple3(cellPtIds.GetId(2)))
                    gN4 = np.asarray(normVec.GetTuple3(cellPtIds.GetId(3)))

                    d   = 0.25*(dN1 + dN2 + dN3 + dN4)
                    g   = 0.25*(gN1 + gN2 + gN3 + gN4)

                    if d <= self.m_InputData.m_TracerRadius:
                        ## check collision with caps
                        ## Case-13 - 4 caps excluding cannula inlet
                        ## check separetly along cap 1, 2, 3, 4
                        ## create array of cap normal vectors, list of cap points
                        ## for n in range(numCaps) :
                        ##     d_cap[n] = g_p[n] . (posP[0]-xp[n], posP[1]-yp[n], posP[n]-zp[n])/ |g_p[n]|
                        ##     if d_cap <= self.m_InputData.m_TracerRadius:
                        ##         v_i = 0.0
                        ##     else :
                        ##          vn      = np.dot(v_i, g)
                        ##          vt      = v_i - vn*g
                        ##          v_i     = vt - restitution*g*vn
                        vn      = np.dot(v_i, g)
                        vt      = v_i - vn*g
                        v_i     = vt - restitution*g*vn

        #--------------------------------
        # set new xyz for  point in loop
        #--------------------------------
            self.m_LagrangianData.setX(p, posP)
            self.m_LagrangianData.setVectorData(v_i, a_ArrayName='Velocity', a_DataID=p)
            # print('setV = ', p, v_i)


    def vtkPolydataDistanceField(self, a_MeshFile, a_SurfaceFile, a_OutFile):

        mReader = vtk.vtkXMLUnstructuredGridReader()
        mReader.SetFileName(a_MeshFile)
        mReader.Update()
        meshData = mReader.GetOutput()

        # sReader = vtk.vtkXMLPolyDataReader()
        # sReader.SetFileName(a_SurfaceFile)
        # sReader.Update()
        # surfData = sReader.GetOutput()

        # numGridPoints = meshData.GetNumberOfPoints()

        # sdf = vtk.vtkImplicitPolyDataDistance()
        # sdf.SetInput(surfData)

        # distance = vtk.vtkDoubleArray()
        # distance.SetName('sdf')
        # distance.SetNumberOfTuples(numGridPoints)
        # distance.SetNumberOfComponents(1)

        # for p in range(numGridPoints):
        #     sys.stdout.write('Computing Point '+str(p)+' of a total '+str(numGridPoints)+'\r')
        #     sys.stdout.flush()
        #     xyz = meshData.GetPoint(p)
        #     distance.SetValue(p, sdf.EvaluateFunction(xyz))

        # meshData.GetPointData().AddArray(distance)

        gradientData = vtk.vtkGradientFilter()
        gradientData.SetInputData(meshData)
        gradientData.ComputeGradientOn()
        gradientData.ComputeQCriterionOff()
        gradientData.ComputeDivergenceOff()
        gradientData.ComputeVorticityOff()
        gradientData.SetInputArrayToProcess(0,0,0,vtk.vtkDataObject.FIELD_ASSOCIATION_POINTS, 'sdf')
        gradientData.SetResultArrayName('sdf-grad')
        gradientData.Update()
        gradientData = gradientData.GetOutput()

        mWriter = vtk.vtkXMLUnstructuredGridWriter()
        mWriter.SetFileName(a_OutFile)
        mWriter.SetInputData(gradientData)
        mWriter.SetDataModeToBinary()
        mWriter.Update()
        mWriter.Write()
