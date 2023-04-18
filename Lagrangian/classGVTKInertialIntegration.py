#-------------------------------------------------------------------------------
# This is a module that encapsulates the physics and numerics of an inertial
# model integration calculation through a specified flow field
#
# Author:      Joseph Wilson
# Institution: University of Colorado, Boulder
# Last Edit:   Sept 2021
#-------------------------------------------------------------------------------

import sys, os, vtk, warnings
import numpy as np
from decimal import Decimal

try:
    from classGVTKGrid import *
    from classGVTKLagrangianData import *
    from classGVTKGenericProblem import *
except:
    sys.exit('Could not import user defined modules')

class GVTKInertial(GVTKProblem):

    def __init__(self, a_InputData, a_DebugMode=False):

        super().__init__(a_InputData, a_DebugMode)
        self.m_LagrangianData = None
        self.m_GridData       = None

    @property
    def numParticles(self):
        return self.m_LagrangianData.numParticles

    @property
    def numCells(self):
        return self.m_GridData.numCells

    def getDtPrecision(self):
        self.m_dt_prec = len(str(self.m_dT).split('.')[1])
        self.m_QuantizeVal = '1.' + '0'*self.m_dt_prec

    def setGrav(self, a_Dir, a_Dim, a_Val):
        base_grav = np.ones(a_Dim) * a_Val

        if a_Dir[0] == '-':
            base_grav *= -1

        if a_Dim == 3:
            if a_Dir[1] == 'x':
                base_grav *= np.array([1, 0, 0])
            elif a_Dir[1] == 'y':
                base_grav *= np.array([0, 1, 0])
            else:
                base_grav *= np.array([0, 0, 1])
        elif a_Dim == 2:
            if a_Dir[1] == 'x':
                base_grav *= np.array([1, 0])
            else:
                base_grav *= np.array([0, 1])

        self.m_Grav = base_grav

    def checkTimeStep(self, dT):
        if self.m_isConstantSize:
            tauV = self.m_MomentumResponseTime[0]
        else:
            # currently based on the average particle radius given by user
            tauV = 2.0*self.m_InputData.m_TracerDensity*self.m_InputData.m_AverageRadius**2/(9.0*self.m_InputData.m_FluidViscosity)
        compare_dT = dT/tauV
        if compare_dT > 0.5:
            recommended_dT = 0.5*tauV
            warn_message = 'Integration time step is high relative to the momentum response time which could cause unstable results. Recommend resetting time step to no greater than ' + str(round(recommended_dT, self.m_dt_prec+2)) + 's.'
            warnings.warn(warn_message)
            print()

    def calculateAverageVelocity(self):#, a_T, a_T0, a_T1, a_DataSyncGrid=True):
        problemDimension  = self.m_InputData.getProblemDimension()
        numFlowDataPoints = self.m_GridData.m_vtkData.GetNumberOfPoints()
        dataArray_0       = self.m_GridData.m_vtkData.GetPointData().GetArray('v0')
        v                 = 0.0

        for p in range(numFlowDataPoints):
            if problemDimension == 1:
                v += dataArray_0.GetTuple1(p)
            elif problemDimension == 2:
                vel = dataArray_0.GetTuple2(p)
                v += np.sqrt(vel[0]**2 + vel[1]**2)
            else:
                vel = dataArray_0.GetTuple3(p)
                v += np.sqrt(vel[0]**2 + vel[1]**2 + vel[2]**2)

        self.m_AverageFlowVelocity = v/numFlowDataPoints

    def calculateMomentumResponseTime(self, a_ID):
        if not hasattr(self, 'm_MomentumResponseTime'):
            self.m_MomentumResponseTime = np.zeros(self.numParticles, dtype=np.float32)
        tauV = (self.m_InputData.getTracerDensity()*(2.0*self.m_LagrangianData.getScalarData(a_ArrayName='Radius', a_DataID=a_ID))**2)/(18.0*self.m_InputData.getFluidViscosity())
        self.m_MomentumResponseTime[a_ID] = tauV

    def calculateStokesNumber(self, a_ID):
        if not hasattr(self, 'm_StokesNumber'):
            self.m_StokesNumber = np.zeros(self.numParticles, dtype=np.float32)
        self.m_StokesNumber[a_ID] = 0.5*self.m_MomentumResponseTime[a_ID]*self.m_AverageFlowVelocity/self.m_LagrangianData.getScalarData(a_ArrayName='Radius', a_DataID=a_ID)

    def calculateSlipVelocityReynoldsNumber(self, vel, a_ID):
        if not hasattr(self, 'm_SlipVelocityReynoldsNumber'):
            self.m_SlipVelocityReynoldsNumber = np.zeros(self.numParticles, dtype=np.float32)

        v_p      = self.m_LagrangianData.getVectorData(a_ArrayName='Velocity', a_DataID=a_ID)
        vel_norm = np.sqrt(( vel[0]-v_p[0])**2 + (vel[1]-v_p[1])**2 + (vel[2]-v_p[2])**2)
        rho      = self.m_InputData.getTracerDensity()
        mu       = self.m_InputData.getFluidViscosity()
        d_p        = self.m_LagrangianData.getScalarData(a_ArrayName='Radius', a_DataID=a_ID)*2.0
        self.m_SlipVelocityReynoldsNumber[a_ID] = rho*vel_norm*d_p/mu

    def setVelocityDistribution(self):

        averageVel  = 20.0
        vtpFileName = self.m_InputData.m_TracerInput
        stlFileName = vtpFileName.split('.')[0]+'.stl'

        #--------------------------
        # calculate normal of cell
        #--------------------------
        def calculatePointNormal(vtkObject, cellId):
            point0 = np.array(stlData.GetCell(cellId).GetPoints().GetPoint(0))
            point1 = np.array(stlData.GetCell(cellId).GetPoints().GetPoint(1))
            point2 = np.array(stlData.GetCell(cellId).GetPoints().GetPoint(2))
            return np.cross(point0-point1, point0-point2)/np.linalg.norm(np.cross(point0-point1, point0-point2))

        #---------------------------
        # numpy array of vtp points
        #
        # already there in code
        #---------------------------
        readerVTP = vtk.vtkXMLPolyDataReader()
        readerVTP.SetFileName(vtpFileName)
        readerVTP.Update()
        vtpData = readerVTP.GetOutput()

        numVTPPoints = vtpData.GetNumberOfPoints()
        vtpPoints = np.zeros((numVTPPoints, 3), dtype=np.float32)

        for i in range(numVTPPoints):
            vtpPoints[i] = vtpData.GetPoint(i)

        #---------------------------
        # numpy array of stl points
        #---------------------------
        readerSTL = vtk.vtkSTLReader()
        readerSTL.SetFileName(stlFileName)
        readerSTL.Update()
        stlData = readerSTL.GetOutput()
        stlData.BuildLinks()

        numSTLPoints = stlData.GetNumberOfPoints()
        stlPoints = np.zeros((numSTLPoints, 3), dtype=np.float32)

        for i in range(numSTLPoints):
            stlPoints[i] = stlData.GetPoint(i)

        #----------------------------
        # create connectivity matrix
        #----------------------------
        numCells = stlData.GetNumberOfCells()
        connectivity = []
        for i in range(numSTLPoints):
            cellIDsUsed = []
            for j in range(numCells):
                if stlData.IsPointUsedByCell(i, j):
                    cellIDsUsed.append(j)
            connectivity.append(cellIDsUsed)

        #---------------------------------
        # calculate normal of cell in stl
        #---------------------------------
        normals = []
        for i in range(numSTLPoints):
            cellIDs = connectivity[i]
            pointNormal = calculatePointNormal(stlData, cellIDs[0])
            for j in range(1, len(cellIDs)):
                newRow = calculatePointNormal(stlData, i)
                pointNormal = np.vstack([pointNormal, newRow])
            pointNormal = np.mean(pointNormal, axis=0)
            normals.append(pointNormal)
        normals = np.array(normals)

        #------------------------------
        # create velocity distribution
        #------------------------------
        velMags = np.random.normal(loc=20.0, size=numVTPPoints)
        InitialVelocity = np.zeros((numVTPPoints, 3))
        for i in range(numVTPPoints):
            InitialVelocity[i] = velMags[i] * normals[i]

        self.m_LagrangianData.addVectorData('Velocity', InitialVelocity)

    def runCompute(self):

        #-----------------------------------------------------------------------
        # initialize inertial integration from an external VTK polydata file and
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
        timeWindowDict = self.m_InputData.getDataTimeWindows_v2()

        #-----------------------------------------------------------------------------
        # display the time windows in case the integration involves unsteady flow data
        #-----------------------------------------------------------------------------
        if not(self.m_InputData.isSteadyFlowData()):
            self.m_InputData.printDataTimeWindows()

        #-----------------------------------------------------
        # check for uniform paricle size or size distribution
        #-----------------------------------------------------
        self.m_isConstantSize = self.m_InputData.m_isConstantSize
        if self.m_isConstantSize:
            TracerRadius = self.m_InputData.m_TracerRadius
        else:
            TracerRadius = np.random.weibull(self.m_InputData.m_ShapeParam, self.m_LagrangianData.m_NumParticles) * self.m_InputData.m_AverageRadius
        self.m_LagrangianData.addScalarData('Radius', TracerRadius)

        #----------------------------------------------
        # calculate inverse momentum response times(s)
        #----------------------------------------------
        for i in range(self.numParticles):
            self.calculateMomentumResponseTime(i)
        self.m_LagrangianData.addScalarData('Momentum Response Time', self.m_MomentumResponseTime)

        #-------------------------------------------
        # set integration time step based on inputs
        #-------------------------------------------
        self.m_dT = self.m_InputData.m_Dt
        self.getDtPrecision()
        self.checkTimeStep(self.m_dT)

        #--------------------
        # start time counter
        #--------------------
        timeIndex = 0
        simTime   = self.m_InputData.getSimulationStartTime()
        stopTime  = self.m_InputData.getSimulationStopTime()
        startTime = simTime
        tWin_0    = timeWindowDict['T_Low'][0]
        tWin_1    = timeWindowDict['T_Up'][0]
        intScheme = self.m_InputData.getIntegrationScheme()
        numSteps  = stopTime/self.m_dT

        self.m_isSteadyFlow = self.m_InputData.isSteadyFlowData()
        if self.m_isSteadyFlow:
            tWin_2 = tWin_1
        else:
            tWin_2 = timeWindowDict['T_Up'][1]

        #--------------------------------------
        # check if OpenFoam data is being used
        #--------------------------------------
        isOpenFoam = False
        if self.m_InputData.m_FlowFieldName == 'U':
            flowFiles = os.listdir(self.m_InputData.m_FlowDirectory)
            for i in range(len(flowFiles)):
                check_start = flowFiles[i].startswith(self.m_InputData.m_FlowFileTag)
                check_end   = flowFiles[i].endswith('.foam')
                if check_start and check_end:
                    isOpenFoam = True

        #-------------
        # set gravity
        #-------------
        if self.m_InputData.m_GravEffects:
            self.setGrav(self.m_InputData.m_GravDir, self.m_InputData.m_SpaceDimension, self.m_InputData.m_GravVal)
        else:
            self.m_Grav = np.zeros(self.m_InputData.m_SpaceDimension)

        #-----------------------------------
        # set initial velocity of particles
        #-----------------------------------
        if True:
            self.setVelocityDistribution()
        else:
            InitialVelocity = np.zeros((self.numParticles, 3))
            InitialVelocity[:,0] = float(self.m_InputData.m_InitialVelocity['u'])
            InitialVelocity[:,1] = float(self.m_InputData.m_InitialVelocity['v'])
            InitialVelocity[:,2] = float(self.m_InputData.m_InitialVelocity['w'])
            self.m_LagrangianData.addVectorData('Velocity', InitialVelocity)

        '''
        implement new velocity code here.

        Need to account for:
            1. normal of closest surface (maybe passing stl corresponding to vtp file would be easiest)
            2. velocity distribution
                a. still need to do some research for this
                b. found ranges but no descriptions of curve (normal, rosin-rammler, etc)
                c. found no relation of initial velocity to size
        '''

        #------------------------------
        # check for resumed simulation
        #------------------------------
        self.m_ResumeSimulation = self.m_InputData.m_ResumeSimulation
        if self.m_ResumeSimulation:

            #--------------------------------------------------------------
            # get list of files already output and assign initial velocity
            #--------------------------------------------------------------
            dest_folder = self.m_InputData.m_RootPath + '/'.join(self.m_InputData.m_TracerOutputFile.split('/')[:-1]) + '/'
            prev_files = os.listdir(dest_folder)

            #----------------------------
            # find base file information
            #----------------------------
            base_file = '_'.join(prev_files[0].split('.')[0].split('_')[:-1])
            file_ext  = '.' + prev_files[0].split('.')[1]

            #--------------------------------------------
            # trim list of files down to just time index
            #--------------------------------------------
            for i in range(len(prev_files)):
                prev_files[i] = int(prev_files[i].split('_')[-1].split('.')[0])

            #-------------------------------------------------
            # update time index and update reader to max file
            #-------------------------------------------------
            timeIndex = max(prev_files)
            last_file = dest_folder + base_file + '_' + str(timeIndex) + file_ext
            reader    = vtk.vtkPolyDataReader()
            reader.SetFileName(last_file)
            reader.Update()

            #---------------------------
            # calculate current simTime
            #---------------------------
            simTime += self.m_InputData.m_WriteInterval * self.m_dT * timeIndex

            #-----------------------------------------------------------
            # get velocity, radius, and time data from last output file
            #-----------------------------------------------------------
            velStruct = reader.GetOutput().GetPointData().GetArray('Velocity')
            radStruct = reader.GetOutput().GetPointData().GetArray('Radius')
            MRTStruct = reader.GetOutput().GetPointData().GetArray('Momentum Response Time')

            for p in range(self.numParticles):

                point = reader.GetOutput().GetPoint(p)
                self.m_LagrangianData.setX(p, point)

                velocity = np.array(velStruct.GetTuple3(p))
                self.m_LagrangianData.setVectorData(velocity, a_ArrayName='Velocity', a_DataID=p)

                radius = radStruct.GetTuple1(p)
                self.m_LagrangianData.setScalarData(radius, a_ArrayName='Radius', a_DataID=p)

                self.m_MomentumResponseTime[p] = MRTStruct.GetTuple1(p)
                self.m_LagrangianData.setScalarData(self.m_MomentumResponseTime[p], a_ArrayName='Momentum Response Time', a_DataID=p)

            #--------------------
            # print resume point
            #--------------------
            print('Resuming simulation at', Decimal(str(simTime)).quantize(Decimal(self.m_QuantizeVal)), 's')

        #-----------------------------------------------------------------------
        # start the simulation time loop
        #-----------------------------------------------------------------------
        while simTime <= stopTime:

            #--------------------------------------------------------------------
            # STEP 0: display a status message for the time step being integrated
            # output initial position
            #--------------------------------------------------------------------
            if self.m_isSteadyFlow and simTime >= stopTime:
                print("Integrating from", timeIndex, "to", timeIndex + 1, "| simTime", Decimal(str(simTime)).quantize(Decimal(self.m_QuantizeVal)))
            elif not self.m_isSteadyFlow and (timeWindowDict['T_Low'][timeIndex] != tWin_0 or timeWindowDict['T_Up'][timeIndex] != tWin_1):
                print("Integrating from", timeIndex, "to", timeIndex + 1, "| simTime", Decimal(str(simTime)).quantize(Decimal(self.m_QuantizeVal)))
            elif simTime > 0 and not self.m_ResumeSimulation:
                print("Integrating from", timeIndex, "to", timeIndex + 1, "| simTime", Decimal(str(simTime)).quantize(Decimal(self.m_QuantizeVal)), end='\r')

            #-------------------------------------------------------------------------------
            # STEP 1: create a boolean condition that governs when new data files are loaded
            #-------------------------------------------------------------------------------
            if self.m_isSteadyFlow:

                isLoadFrame     = (simTime == startTime)
                isDataPointSync = True

                if self.m_ResumeSimulation:
                    isLoadFrame = True

            else:

                isLoadFrame = (simTime == startTime) \
                        or (timeWindowDict['T_Low'][timeIndex] != tWin_0 or timeWindowDict['T_Up'][timeIndex] != tWin_1)

                isDataPointSync = timeWindowDict['T_Low'][timeIndex] == timeWindowDict['T_Up'][timeIndex]

                if self.m_ResumeSimulation:
                    isLoadFrame = True

            #-------------------------------------------------------------------
            # a set of debug messages for program status
            #-------------------------------------------------------------------
            if isLoadFrame:
                print("Will Load Velocity Data")

            if isDataPointSync:
                if self.m_isSteadyFlow:
                    if simTime == startTime or self.m_ResumeSimulation:
                        print("Data And Integration Times Are Synced")
                else:
                    print("Data And Integration Times Are Synced")

            #-------------------------------------------------------------------
            # STEP 3: resolve injections and inject particles
            # for this inertial particle integration we are injecting particles
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

                isSingleStageIntegration = False

                if self.m_isSteadyFlow:

                    flowSingle = self.m_InputData.getFlowDataFileName(a_Foam=isOpenFoam)
                    self.m_GridData.getDataFromFile(flowSingle)
                    self.m_GridData.renameDataArray(self.m_InputData.getVelDataName(), 'v0')
                    # self.m_GridData.addDataArrayFromFile(boundaryFileName, 'dist', 'sdf')
                    # self.m_GridData.addVectorData()

                else:

                    #-----------------------------------------------------------------------------------
                    # TODO: For the unsteady file read write, the data array updates have to be modified
                    # using the add and remove data array functionalities for the GVTKGridData object
                    #-----------------------------------------------------------------------------------
                    tWin_0  = timeWindowDict['T_Low'][timeIndex]    # updating the time window
                    tWin_1  = timeWindowDict['T_Up'][timeIndex]     # updating the time window
                    if timeWindowDict['T_Up'][timeIndex] == timeWindowDict['T_Up'][-1]:
                        tWin_2 = tWin_1
                    else:
                        tWin_2 = timeWindowDict['T_Up'][timeIndex+1]

                    if simTime == startTime:

                        flowInitial = self.m_InputData.getFlowDataFileName(a_ID = timeWindowDict['ID_Low'][timeIndex])
                        self.m_GridData.getDataFromFile(flowInitial)
                        self.m_GridData.renameDataArray(self.m_InputData.getVelDataName(), 'v0')

                    elif self.m_ResumeSimulation:

                        if isDataPointSync:

                            flowInitial = self.m_InputData.getFlowDataFileName(a_ID = timeWindowDict['ID_Low'][timeIndex])
                            self.m_GridData.getDataFromFile(flowInitial)
                            self.m_GridData.renameDataArray(self.m_InputData.getVelDataName(), 'v0')

                        else:

                            if isSingleStageIntegration == True:

                                flowSingle = self.m_InputData.getFlowDataFileName(a_ID = timeWindowDict['ID_Low'][timeIndex])
                                self.m_GridData.addDataArrayFromFile(flowSingle, 'v0', self.m_InputData.getVelDataName())

                            if simTime < tWin_1 - dT:

                                flowPlus        = self.m_InputData.getFlowDataFileName(a_ID = timeWindowDict['ID_Up'][timeIndex])
                                flowMinus       = self.m_InputData.getFlowDataFileName(a_ID = timeWindowDict['ID_Low'][timeIndex])
                                self.m_GridData.addDataArrayFromFile(flowPlus, 'v1', self.m_InputData.getVelDataName())
                                self.m_GridData.addDataArrayFromFile(flowMinus, 'v0', self.m_InputData.getVelDataName())
                                print('Loaded a new velocity data')

                            elif simTime >= tWin_1 - dT:

                                flowPlus    = self.m_InputData.getFlowDataFileName(a_ID = timeWindowDict['ID_Up'][timeIndex])
                                flowMinus   = self.m_InputData.getFlowDataFileName(a_ID = timeWindowDict['ID_Low'][timeIndex])
                                flowNext    = self.m_InputData.getFlowDataFileName(a_ID = timeWindowDict['ID_Low'][timeIndex]+1)
                                self.m_GridData.addDataArrayFromFile(flowPlus, 'v1', self.m_InputData.getVelDataName())
                                self.m_GridData.addDataArrayFromFile(flowMinus, 'v0', self.m_InputData.getVelDataName())
                                self.m_GridData.addDataArrayFromFile(flowNext, 'v2', self.m_InputData.getVelDataName())
                                print('Loaded a new velocity data')

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

                            if (simTime + self.m_dT) < tWin_1:

                                self.m_GridData.removeDataArray('v0')
                                self.m_GridData.removeDataArray('v1')
                                self.m_GridData.removeDataArray('v2')
                                flowPlus  = self.m_InputData.getFlowDataFileName(a_ID = timeWindowDict['ID_Up'][timeIndex])
                                flowMinus = self.m_InputData.getFlowDataFileName(a_ID = timeWindowDict['ID_Low'][timeIndex])
                                self.m_GridData.addDataArrayFromFile(flowPlus, 'v1', self.m_InputData.getVelDataName())
                                self.m_GridData.addDataArrayFromFile(flowMinus, 'v0', self.m_InputData.getVelDataName())
                                print('Loaded a new velocity data')

                            elif (simTime + self.m_dT) >= tWin_1:

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
            if (simTime == startTime and self.m_InputData.isFixedMesh()) or (self.m_ResumeSimulation and self.m_InputData.isFixedMesh()):
                print("Building Cell Locator Maps")
                self.m_GridData.buildLocator()

            #-------------------------------------------------------------------
            # STEP 6: Output initial position. SCalculate average flow velocity,
            # Stokes Number, Slip Velocity Reynolds Number
            #-------------------------------------------------------------------
            if simTime == startTime or self.m_ResumeSimulation:
                self.calculateAverageVelocity()
                for i in range(self.numParticles):
                    self.calculateStokesNumber(i)
                    if True: #a_PolygonalCells:
                        vel_i = self.m_GridData.gridInterpolateAveraging(self.m_LagrangianData.getX(i), 'v0')
                    else:
                        vel_i = self.m_GridData.gridInterpolateNodalBasis(self.m_LagrangianData.getX(i), 'v0')
                    self.calculateSlipVelocityReynoldsNumber(vel_i, i)
                self.m_LagrangianData.addScalarData('Stokes Number', self.m_StokesNumber)
                self.m_LagrangianData.addScalarData('Slip Velocity Reynolds Number', self.m_SlipVelocityReynoldsNumber)
                self.m_LagrangianData.writeData(self.m_InputData.getTracerOutputFile(a_ID1=0, a_ID2=timeIndex))

                #-----------------------------------------------------------------
                # self.m_ResumeSimulation simulation now needs to be set to False
                #-----------------------------------------------------------------
                self.m_ResumeSimulation = False

            #----------------------------------------------------------------
            # STEP 7: now proceed with integration of each inertial particle
            #----------------------------------------------------------------
            boundaryCondition = 1   # currently a dummy variable
            if intScheme == 'feu':
                self.integrateForwardEuler(simTime, tWin_0, tWin_1, boundaryCondition, a_DataSync=isDataPointSync, a_PolygonalCells=True)
            elif intScheme == 'rk4':
                self.integrateRK4(simTime, tWin_0, tWin_1, tWin_2, boundaryCondition, a_DataSync=isDataPointSync, a_PolygonalCells=True)

            #------------------------------------------------
            # STEP 8: update time indices and simulation time
            #------------------------------------------------
            timeIndex += 1
            simTime    = round(simTime+self.m_dT, self.m_dt_prec)

            #----------------------------------------------------------------------
            # STEP 9: at the end of appropriate number of steps dump data into file
            #----------------------------------------------------------------------
            if timeIndex%self.m_InputData.m_WriteInterval == 0:
                self.m_LagrangianData.writeData(self.m_InputData.getTracerOutputFile(a_ID1=0, a_ID2=timeIndex))



    def integrateForwardEuler(self, simTime, tWin_0, tWin_1, boundaryCondition, a_DataSync=True, a_PolygonalCells=True, a_StoreXYZ=False):
        '''
        Forward Euler time integration routine for individual particles

        Parameters
        ----------
        a_T : float

        a_T0 : float

        a_T1 : float

        a_BoundaryCondition : dummy

        a_DataSync : bool, optional

        a_PolygonalCells : bool, optional

        Returns
        -------
        none
        '''

        #----------------------
        # simulation constants
        #----------------------
        dT      = self.m_dT
        buoy    = self.m_Grav*(1.0-self.m_InputData.getFluidDensity()/self.m_InputData.getTracerDensity())

        #----------------------------
        # calculate average velocity
        #----------------------------
        if not self.m_isSteadyFlow and a_DataSync:
            self.calculateAverageVelocity()

        for p in range(self.numParticles):

            tauV  = self.m_MomentumResponseTime[p]
            xyz_i = self.m_LagrangianData.getX(p)
            v_i   = np.array(self.m_LagrangianData.getVectorData(a_ArrayName='Velocity', a_DataID=p))

            if a_DataSync:
                if a_PolygonalCells == True:
                    vel, stat = self.m_GridData.gridInterpolateAveraging(xyz_i, 'v0', a_GetStatus=True)
                else:
                    vel, stat = self.m_GridData.gridInterpolateNodalBasis(xyz_i, 'v0', a_GetStatus=True)
            else:
                if a_PolygonalCells == True:
                    velPlus, stat  = self.m_GridData.gridInterpolateAveraging(xyz_i, 'v1', a_GetStatus=True)
                    velMinus, stat = self.m_GridData.gridInterpolateAveraging(xyz_i, 'v0', a_GetStatus=True)
                else:
                    velPlus, stat  = self.m_GridData.gridInterpolateNodalBasis(xyz_i, 'v1', a_GetStatus=True)
                    velMinus, stat = self.m_GridData.gridInterpolateNodalBasis(xyz_i, 'v0', a_GetStatus=True)

                vel = velMinus + (simTime - tWin_0) * (velPlus - velMinus)/(tWin_1 - tWin_0)

            #---------------------------
            # forward euler integration
            #---------------------------
            if stat != -1:
                v_new     = v_i + dT*((vel-v_i)/tauV + buoy)
                xyz_new   = xyz_i + v_new*dT

            else:
                v_new = np.zeros(self.m_InputData.m_SpaceDimension)
                xyz_new = xyz_i

            #-------------------------------
            # set new xyz for point in loop
            #-------------------------------
            self.m_LagrangianData.setX(p, xyz_new)
            self.m_LagrangianData.setVectorData(v_new, a_ArrayName='Velocity', a_DataID=p)

            #------------------------------------------------------------
            # calculate new stokes number and update in m_LagrangianData
            #------------------------------------------------------------
            self.calculateStokesNumber(p)
            self.m_LagrangianData.setScalarData(self.m_StokesNumber[p], a_ArrayName='Stokes Number', a_DataID=p)

            #------------------------------------------------------------------------
            # calculate slip velocity Reynolds number and update in m_LagrangianData
            #------------------------------------------------------------------------
            self.calculateSlipVelocityReynoldsNumber(vel, p)
            self.m_LagrangianData.setScalarData(self.m_SlipVelocityReynoldsNumber[p], a_ArrayName='Slip Velocity Reynolds Number', a_DataID=p)



    def integrateRK4(self, a_T, a_T0, a_T1, a_TNext, a_BoundaryCondition, a_DataSync=True, a_PolygonalCells=False):

        #----------------------
        # simulation constants
        #----------------------
        dT      = self.m_dT
        buoy    = self.m_Grav*(1.0-self.m_InputData.getFluidDensity()/self.m_InputData.getTracerDensity())

        #----------------------------
        # calculate average velocity
        #----------------------------
        if not self.m_isSteadyFlow and a_DataSync:
            self.calculateAverageVelocity()

        for p in range(self.numParticles):

            tauV  = self.m_MomentumResponseTime[p]
            xyz_1 = self.m_LagrangianData.getX(p)
            v_i   = np.array(self.m_LagrangianData.getVectorData(a_ArrayName='Velocity', a_DataID=p))

            if a_DataSync:

                if a_PolygonalCells:

                    vel_1, stat = self.m_GridData.gridInterpolateAveraging(xyz_1, 'v0', a_GetStatus=True)

                    if stat != -1:

                        k_1   = v_i
                        el_1  = (vel_1-v_i)/tauV + buoy

                        xyz_2 = xyz_1 + 0.5*k_1*dT
                        vel_2 = self.m_GridData.gridInterpolateAveraging(xyz_2, 'v0')
                        k_2   = v_i + 0.5*el_1*dT
                        el_2  = (vel_2-v_i-0.5*el_1*dT)/tauV + buoy

                        xyz_3 = xyz_1 + 0.5*k_2*dT
                        vel_3 = self.m_GridData.gridInterpolateAveraging(xyz_3, 'v0')
                        k_3   = v_i + 0.5*el_2*dT
                        el_3  = (vel_3-v_i-0.5*el_2*dT)/tauV + buoy

                        xyz_4 = xyz_1 + k_3*dT
                        vel_4 = self.m_GridData.gridInterpolateAveraging(xyz_4, 'v0')
                        k_4   = v_i + el_3*dT
                        el_4  = (vel_4-v_i-el_3*dT)/tauV + buoy

                        v_new   = v_i   + (dT/6)*(el_1 + 2*el_2 + 2*el_3 + el_4)
                        xyz_new = xyz_1 + (dT/6)*(k_1  + 2*k_2  + 2*k_3  + k_4)

                    else:
                        xyz_new = xyz_1
                        v_new   = np.zeros(self.m_InputData.m_SpaceDimension)

                else:

                    vel_1, stat = self.m_GridData.gridInterpolateNodalBasis(xyz_1, 'v0', a_GetStatus=True)

                    if stat != -1:

                        k_1   = v_i
                        el_1  = (vel_1-v_i)/tauV + buoy

                        xyz_2 = xyz_1 + 0.5*k_1*dT
                        vel_2 = self.m_GridData.gridInterpolateNodalBasis(xyz_2, 'v0')
                        k_2   = v_i + 0.5*el_1*dT
                        el_2  = (vel_2-v_i-0.5*el_1*dT)/tauV + buoy

                        xyz_3 = xyz_1 + 0.5*k_2*dT
                        vel_3 = self.m_GridData.gridInterpolateNodalBasis(xyz_3, 'v0')
                        k_3   = v_i + 0.5*el_2*dT
                        el_3  = (vel_3-v_i-0.5*el_2*dT)/tauV + buoy

                        xyz_4 = xyz_1 + k_3*dT
                        vel_4 = self.m_GridData.gridInterpolateNodalBasis(xyz_4, 'v0')
                        k_4   = v_i + el_3*dT
                        el_4  = (vel_4-v_i-el_3*dT)/tauV + buoy

                        v_new   = v_i   + (dT/6)*(el_1 + 2*el_2 + 2*el_3 + el_4)
                        xyz_new = xyz_1 + (dT/6)*(k_1  + 2*k_2  + 2*k_3  + k_4)

                    else:
                        xyz_new = xyz_1
                        v_new   = np.zeros(self.m_InputData.m_SpaceDimension)

            else:

                t_1 = a_T
                t_2 = round(a_T + 0.50*dT, self.m_dt_prec+1)
                t_3 = round(a_T + 0.50*dT, self.m_dt_prec+1)
                t_4 = round(a_T + dT, self.m_dt_prec)

                #--------
                # step 1
                #--------
                if a_PolygonalCells == True:
                    vel_1_L, stat = self.m_GridData.gridInterpolateAveraging(xyz_1, 'v0', a_GetStatus=True)
                    vel_1_U = self.m_GridData.gridInterpolateAveraging(xyz_1, 'v1')
                else:
                    vel_1_L, stat = self.m_GridData.gridInterpolateNodalBasis(xyz_1, 'v0', a_GetStatus=True)
                    vel_1_U = self.m_GridData.gridInterpolateNodalBasis(xyz_1, 'v1')

                #---------------------------
                # particles still suspended
                #---------------------------
                if stat != -1:

                    vel_1 = vel_1_L + (a_T - a_T0)*(vel_1_U - vel_1_L)/(a_T1 - a_T0)
                    k_1   = v_i
                    el_1  = (vel_1-v_i)/tauV + buoy

                    #--------
                    # step 2
                    #--------
                    xyz_2 = xyz_1 + 0.50*dT*k_1

                    if t_2 <= a_T1:
                        if a_PolygonalCells == True:
                            vel_2_L = self.m_GridData.gridInterpolateAveraging(xyz_2, 'v0')
                            vel_2_U = self.m_GridData.gridInterpolateAveraging(xyz_2, 'v1')
                        else:
                            vel_2_L = self.m_GridData.gridInterpolateNodalBasis(xyz_2, 'v0')
                            vel_2_U = self.m_GridData.gridInterpolateNodalBasis(xyz_2, 'v1')
                        vel_2 = vel_2_L + (a_T - a_T0)*(vel_2_U - vel_2_L)/(a_T1 - a_T0)

                    elif t_2 > a_T1:
                        if a_PolygonalCells == True:
                            vel_2_L = self.m_GridData.gridInterpolateAveraging(xyz_2, 'v1')
                            vel_2_U = self.m_GridData.gridInterpolateAveraging(xyz_2, 'v2')
                        else:
                            vel_2_L = self.m_GridData.gridInterpolateNodalBasis(xyz_2, 'v1')
                            vel_2_U = self.m_GridData.gridInterpolateNodalBasis(xyz_2, 'v2')
                        vel_2 = vel_2_L + (a_T - a_T1)*(vel_2_U - vel_2_L)/(a_TNext - a_T1)

                    k_2  = v_i + 0.50*dT*el_1
                    el_2 = (vel_2-k_2)/tauV + buoy

                    #--------
                    # step 3
                    #--------
                    xyz_3 = xyz_1 + 0.50*dT*k_2

                    if t_3 <= a_T1:
                        if a_PolygonalCells == True:
                            vel_3_L = self.m_GridData.gridInterpolateAveraging(xyz_3, 'v0')
                            vel_3_U = self.m_GridData.gridInterpolateAveraging(xyz_3, 'v1')
                        else:
                            vel_3_L = self.m_GridData.gridInterpolateNodalBasis(xyz_3, 'v0')
                            vel_3_U = self.m_GridData.gridInterpolateNodalBasis(xyz_3, 'v1')
                        vel_3 = vel_3_L + (a_T - a_T0)*(vel_3_U - vel_3_L)/(a_T1 - a_T0)

                    elif t_3 > a_T1:
                        if a_PolygonalCells == True:
                            vel_3_L = self.m_GridData.gridInterpolateAveraging(xyz_3, 'v1')
                            vel_3_U = self.m_GridData.gridInterpolateAveraging(xyz_3, 'v2')
                        else:
                            vel_3_L = self.m_GridData.gridInterpolateNodalBasis(xyz_3, 'v1')
                            vel_3_U = self.m_GridData.gridInterpolateNodalBasis(xyz_3, 'v2')
                        vel_3 = vel_3_L + (a_T - a_T1)*(vel_3_U - vel_3_L)/(a_TNext - a_T1)

                    k_3  = v_i + 0.50*dT*el_2
                    el_3 = (vel_3-k_3)/tauV + buoy

                    #--------
                    # step 4
                    #--------
                    xyz_4 = xyz_1 + dT*k_3

                    if t_4 <= a_T1:
                        if a_PolygonalCells == True:
                            vel_4_L = self.m_GridData.gridInterpolateAveraging(xyz_4, 'v0')
                            vel_4_U = self.m_GridData.gridInterpolateAveraging(xyz_4, 'v1')
                        else:
                            vel_4_L = self.m_GridData.gridInterpolateNodalBasis(xyz_4, 'v0')
                            vel_4_U = self.m_GridData.gridInterpolateNodalBasis(xyz_4, 'v1')
                        vel_4 = vel_4_L + (a_T - a_T0)*(vel_4_U - vel_4_L)/(a_T1 - a_T0)

                    elif t_4 > a_T1:
                        if a_PolygonalCells == True:
                            vel_4_L = self.m_GridData.gridInterpolateAveraging(xyz_4, 'v1')
                            vel_4_U = self.m_GridData.gridInterpolateAveraging(xyz_4, 'v2')
                        else:
                            vel_4_L = self.m_GridData.gridInterpolateNodalBasis(xyz_4, 'v1')
                            vel_4_U = self.m_GridData.gridInterpolateNodalBasis(xyz_4, 'v2')
                        vel_4 = vel_4_L + (a_T - a_T1)*(vel_4_U - vel_4_L)/(a_TNext - a_T1)

                    k_4  = v_i + dT*el_3
                    el_4 = (vel_4-k_4)/tauV + buoy

                    #---------------------------
                    # new velocity and position
                    #---------------------------
                    v_new   = v_i   + (dT/6)*(el_1 + 2*el_2 + 2*el_3 + el_4)
                    xyz_new = xyz_1 + (dT/6)*(k_1  + 2*k_2  + 2*k_3  + k_4)

                else:
                    xyz_new = xyz_1
                    v_new   = np.zeros(self.m_InputData.m_SpaceDimension)

            #---------------------------------------
            # set new xyz and v_p for point in loop
            #---------------------------------------
            self.m_LagrangianData.setX(p, xyz_new)
            self.m_LagrangianData.setVectorData(v_new, a_ArrayName='Velocity', a_DataID=p)

            #------------------------------------------------------------
            # calculate new stokes number and update in m_LagrangianData
            #------------------------------------------------------------
            self.calculateStokesNumber(p)
            self.m_LagrangianData.setScalarData(self.m_StokesNumber[p], a_ArrayName='Stokes Number', a_DataID=p)

            #------------------------------------------------------------------------
            # calculate slip velocity Reynolds number and update in m_LagrangianData
            #------------------------------------------------------------------------
            self.calculateSlipVelocityReynoldsNumber(vel_1, p)
            self.m_LagrangianData.setScalarData(self.m_SlipVelocityReynoldsNumber[p], a_ArrayName='Slip Velocity Reynolds Number', a_DataID=p)
