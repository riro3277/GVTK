#---------------------------------------------------------------------------------------------
# This is a module that encapsulates the physics and numerics of a light tracer integration
# calculation through a specified flow field
#
# Author:       Debanjan Mukherjee
# Institution:  University of California, Berkeley
# Last Edit:    June 2018
#---------------------------------------------------------------------------------------------
import sys, os
import vtk
import numpy as np
from decimal import Decimal

try:
    from classGVTKGrid import *
    from classGVTKLagrangianData import *
    from classGVTKGenericProblem import *
except ImportError:
    sys.exit('Could not import user defined modules')

class GVTKDriftDiffuse(GVTKProblem):

    def __init__(self, a_InputData, a_DebugMode=False):

        super().__init__(a_InputData, a_DebugMode)
        self.m_LagrangianData   = None
        self.m_GridData         = None

    @property
    def numParticles(self):
        return self.m_LagrangianData.numParticles

    @property
    def numCells(self):
        return self.m_GridData.numCells

    @property
    def numNodes(self):
        return self.m_GridData.numNodes

    def getDtPrecision(self, dT):
        self.m_dt_prec = len(str(dT).split('.')[1])
        self.m_QuantizeVal = '1.' + '0'*self.m_dt_prec

    def localTurbuleneParams(self, a_PolygonalCells):
        probDim = self.m_InputData.m_SpaceDimension
        self.m_EddyTime = np.zeros(self.numParticles)              # current time in eddy
        self.m_tInt     = np.ones(self.numParticles)               # eddy residence time
        self.m_pRand    = np.zeros((self.numParticles, probDim))   # random value for u'

        self.m_TrackedParticle = 1

        for p in range(self.numParticles):
            xyz_old = self.m_LagrangianData.getX(p)
            if a_PolygonalCells == True:
                omega, stat = self.m_GridData.gridInterpolateAveraging(xyz_old, 'omega0', a_DataDim=1, a_GetStatus=True)
            else:
                omega, stat = self.m_GridData.gridInterpolateNodalBasis(xyz_old, 'omega0', a_DataDim=1, a_GetStatus=True)

            if self.m_ResumeSimulation:

                if stat != -1:

                    self.m_tInt[p] *= 0.6/omega
                    for d in range(probDim):
                        self.m_pRand[p, d] = np.random.randn()

            else:

                if stat == -1:
                    sys.exit('Initial position of a particle is out of the domain. Edit particle file and try again.')

                self.m_tInt[p] *= 0.6/omega
                for d in range(probDim):
                    self.m_pRand[p, d] = np.random.randn()

    def deflection(self, a_xyz, a_V):

        wall = False

        #----------------------------------
        # out of bounds in the x-direction
        #----------------------------------
        if a_xyz[0] <= 0.0 or a_xyz[0] >= 7.5:
            a_V[0] *= -1
            wall = True

        #----------------------------------
        # out of bounds in the y-direction
        #----------------------------------
        elif a_xyz[1] <= 0.0 or a_xyz[1] >= 7.5:

            #-------------------
            # account for vents
            #-------------------
            if a_xyz[1] <= -0.15 or a_xyz[1] >= 7.65:
                a_V = a_V
            else:
                a_V[1] *= -1
                wall = True

        #----------------------------------
        # out of bounds in the z-direction
        #----------------------------------
        elif a_xyz[2] <= 0.0 or a_xyz[2] >= 3.0:
            a_V[2] *= -1
            wall = True

        #-----------------
        # return velocity
        #-----------------
        return a_V, wall

    def runCompute(self):

        #----------------------------------------------------------
        # initialize tracer data from an external VTK polydata file
        # and update the Lagrangian tracer object for the problem
        #----------------------------------------------------------
        self.m_LagrangianData = GVTKLagrangianData(a_File=self.m_InputData.getTracerInput())

        #--------------------------------------------------------------------------
        # initialize blank grid data object as member of physics object
        # Note: The object gets instantiated here, but the actual data manipulation
        # happens at instances when flow or grid data files are read.
        #--------------------------------------------------------------------------
        self.m_GridData = GVTKGridData()

        #-------------------------------------------------------------
        # load tracers from the same spot where particles are released
        #-------------------------------------------------------------
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

        #-----------------------------------------------------------------
        # set boolean to use constant D_eff or determine from data fields
        #-----------------------------------------------------------------
        self.m_LocalTurbulence = True

        #-------------------
        # start time counter
        #-------------------
        timeIndex = 0
        simTime   = self.m_InputData.getSimulationStartTime()
        stopTime  = self.m_InputData.getSimulationStopTime()
        dT        = self.m_InputData.getIntegrationTimeStep()
        numSteps  = stopTime/dT
        tWin_0    = timeWindowDict['T_Low'][0]
        tWin_1    = timeWindowDict['T_Up'][0]

        self.m_isSteadyFlow = self.m_InputData.isSteadyFlowData()
        self.getDtPrecision(dT)

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

            #---------------------------------------------------
            # update time index and get point data of last file
            #---------------------------------------------------
            timeIndex = max(prev_files)
            last_file = dest_folder + base_file + '_' + str(timeIndex) + file_ext
            reader    = vtk.vtkPolyDataReader()
            reader.SetFileName(last_file)
            reader.Update()
            for i in range(reader.GetOutput().GetNumberOfPoints()):
                point = reader.GetOutput().GetPoint(i)
                self.m_LagrangianData.setX(i, point)

            #------------------------------
            # update simTime and startTime
            #------------------------------
            simTime += self.m_InputData.m_WriteInterval * dT * timeIndex
            print('Resuming simulation at', simTime, 's')

        #-------------------------------
        # start the simulation time loop
        #-------------------------------
        while simTime <= stopTime:

            #--------------------------------------------------------------------
            # STEP 0: display a status message for the time step being integrated
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
            if self.m_InputData.isSteadyFlowData():

                isLoadFrame     = (simTime == self.m_InputData.getSimulationStartTime())
                isDataPointSync = True

                if self.m_ResumeSimulation:
                    isLoadFrame = True

            else:

                isLoadFrame = (simTime == self.m_InputData.getSimulationStartTime()) \
                        or (timeWindowDict['T_Low'][timeIndex] != tWin_0 or timeWindowDict['T_Up'][timeIndex] != tWin_1)

                isDataPointSync = timeWindowDict['T_Low'][timeIndex] == timeWindowDict['T_Up'][timeIndex]

                if self.m_ResumeSimulation:
                    isLoadFrame = True

            #-------------------------------------------
            # a set of debug messages for program status
            #-------------------------------------------
            if isLoadFrame:
                print("Will Load Velocity Data")

            if isDataPointSync:
                if self.m_InputData.isSteadyFlowData():
                    if simTime == self.m_InputData.getSimulationStartTime() or self.m_ResumeSimulation:
                        print("Data And Integration Times Are Synced")
                else:
                    print("Data And Integration Times Are Synced")

            #--------------------------------------------------------------
            # STEP 3: resolve injections and inject particles
            # for this simple tracer integration we are injecting particles
            # either always or never.
            # TODO: THIS NEEDS TO BE IMPROVED
            #--------------------------------------------------------------
            isInjectPoints = False

            #-----------------------------------------
            # inject tracers into the domain if needed
            # currently it injects every time-step
            #-----------------------------------------
            if isInjectPoints:
                tracers.injectParticles(tracerInject)
                print("New Particles Injected")
                print("Now Integrating", tracers.numParticles, "Particles")

            #---------------------------------------------------------------
            # STEP 4: load data file based on the evaluated boolean variable
            #---------------------------------------------------------------

            if isLoadFrame:

                isSingleStageIntegration = False

                if self.m_InputData.isSteadyFlowData():

                    flowSingle = self.m_InputData.getFlowDataFileName(a_Foam=isOpenFoam)
                    self.m_GridData.getDataFromFile(flowSingle)
                    self.m_GridData.renameDataArray(self.m_InputData.getVelDataName(), 'v0')

                    if self.m_LocalTurbulence:
                        self.m_GridData.renameDataArray('omega', 'omega0')
                        self.m_GridData.renameDataArray('k',   'k0')

                else:

                    #-----------------------------------------------------------------------------------
                    # TODO: For the unsteady file read write, the data array updates have to be modified
                    # using the add and remove data array functionalities for the GVTKGridData object
                    #-----------------------------------------------------------------------------------
                    tWin_0 = timeWindowDict['T_Low'][timeIndex]    # updating the time window
                    tWin_1 = timeWindowDict['T_Up'][timeIndex]     # updating the time window

                    if simTime == self.m_InputData.getSimulationStartTime():

                        flowInitial = self.m_InputData.getFlowDataFileName(a_ID = timeWindowDict['ID_Low'][timeIndex])
                        self.m_GridData.getDataFromFile(flowInitial)
                        self.m_GridData.renameDataArray(self.m_InputData.getVelDataName(), 'v0')

                        if self.m_LocalTurbulence:
                            self.m_GridData.renameDataArray('omega', 'omega0')
                            self.m_GridData.renameDataArray('k',   'k0')

                    elif self.m_ResumeSimulation:

                        if isDataPointSync:

                            flowInitial = self.m_InputData.getFlowDataFileName(a_ID = timeWindowDict['ID_Low'][timeIndex])
                            self.m_GridData.getDataFromFile(flowInitial)
                            self.m_GridData.renameDataArray(self.m_InputData.getVelDataName(), 'v0')

                            if self.m_LocalTurbulence:
                                self.m_GridData.renameDataArray('omega', 'omega0')
                                self.m_GridData.renameDataArray('k',   'k0')

                        else:

                            if isSingleStageIntegration == True:

                                flowSingle = self.m_InputData.getFlowDataFileName(a_ID = timeWindowDict['ID_Low'][timeIndex])
                                self.m_GridData.addDataArrayFromFile(flowSingle, 'v0', self.m_InputData.getVelDataName())

                                if self.m_LocalTurbulence:
                                    self.m_GridData.renameDataArray('omega', 'omega0')
                                    self.m_GridData.renameDataArray('k',   'k0')

                            if simTime < tWin_1 - dT:

                                flowPlus        = self.m_InputData.getFlowDataFileName(a_ID = timeWindowDict['ID_Up'][timeIndex])
                                flowMinus       = self.m_InputData.getFlowDataFileName(a_ID = timeWindowDict['ID_Low'][timeIndex])
                                self.m_GridData.addDataArrayFromFile(flowPlus, 'v1', self.m_InputData.getVelDataName())
                                self.m_GridData.addDataArrayFromFile(flowMinus, 'v0', self.m_InputData.getVelDataName())

                                if self.m_LocalTurbulence:
                                    self.m_GridData.renameDataArray(flowMinus, 'omega0', 'omega')
                                    self.m_GridData.renameDataArray(flowMinus, 'k0',   'k')
                                    self.m_GridData.renameDataArray(flowPlus,  'omega1', 'omega')
                                    self.m_GridData.renameDataArray(flowPlus,  'k1',   'k')

                                print('Loaded a new velocity data')

                            elif simTime >= tWin_1 - dT:

                                flowPlus    = self.m_InputData.getFlowDataFileName(a_ID = timeWindowDict['ID_Up'][timeIndex])
                                flowMinus   = self.m_InputData.getFlowDataFileName(a_ID = timeWindowDict['ID_Low'][timeIndex])
                                flowNext    = self.m_InputData.getFlowDataFileName(a_ID = timeWindowDict['ID_Low'][timeIndex]+1)
                                self.m_GridData.addDataArrayFromFile(flowPlus, 'v1', self.m_InputData.getVelDataName())
                                self.m_GridData.addDataArrayFromFile(flowMinus, 'v0', self.m_InputData.getVelDataName())
                                self.m_GridData.addDataArrayFromFile(flowNext, 'v2', self.m_InputData.getVelDataName())

                                if self.m_LocalTurbulence:
                                    self.m_GridData.renameDataArray(flowMinus, 'omega0', 'omega')
                                    self.m_GridData.renameDataArray(flowMinus, 'k0',   'k')
                                    self.m_GridData.renameDataArray(flowPlus,  'omega1', 'omega')
                                    self.m_GridData.renameDataArray(flowPlus,  'k1',   'k')
                                    self.m_GridData.renameDataArray(flowNext,  'omega2', 'omega')
                                    self.m_GridData.renameDataArray(flowNext,  'k2',   'k')

                                print('Loaded a new velocity data')

                    elif isDataPointSync:

                        if isSingleStageIntegration == True:

                            self.m_GridData.removeDataArray('v0')
                            flowSingle = self.m_InputData.getFlowDataFileName(a_ID = timeWindowDict['ID_Low'][timeIndex])
                            self.m_GridData.addDataArrayFromFile(flowSingle, 'v0', self.m_InputData.getVelDataName())

                            if self.m_LocalTurbulence:
                                self.m_GridData.removeDataArray('omega0')
                                self.m_GridData.removeDataArray('k0')
                                self.m_GridData.renameDataArray('omega', 'omega0')
                                self.m_GridData.renameDataArray('k',   'k0')

                        else:

                            self.m_GridData.removeDataArray('v0')
                            self.m_GridData.removeDataArray('v1')
                            flowPlus  = self.m_InputData.getFlowDataFileName(a_ID = timeWindowDict['ID_Up'][timeIndex])
                            flowMinus = self.m_InputData.getFlowDataFileName(a_ID = timeWindowDict['ID_Low'][timeIndex])
                            self.m_GridData.addDataArrayFromFile(flowPlus, 'v1', self.m_InputData.getVelDataName())
                            self.m_GridData.addDataArrayFromFile(flowMinus, 'v0', self.m_InputData.getVelDataName())

                            if self.m_LocalTurbulence:
                                self.m_GridData.removeDataArray('omega0')
                                self.m_GridData.removeDataArray('k0')
                                self.m_GridData.removeDataArray('omega1')
                                self.m_GridData.removeDataArray('k1')
                                self.m_GridData.renameDataArray(flowMinus, 'omega0', 'omega')
                                self.m_GridData.renameDataArray(flowMinus, 'k0',   'k')
                                self.m_GridData.renameDataArray(flowPlus,  'omega1', 'omega')
                                self.m_GridData.renameDataArray(flowPlus,  'k1',   'k')

                            print('Loaded a new velocity data')

                    else:

                        if isSingleStageIntegration == True:

                            self.m_GridData.removeDataArray('v0')
                            self.m_GridData.removeDataArray('v1')
                            flowPlus  = self.m_InputData.getFlowDataFileName(a_ID = timeWindowDict['ID_Up'][timeIndex])
                            flowMinus = self.m_InputData.getFlowDataFileName(a_ID = timeWindowDict['ID_Low'][timeIndex])
                            self.m_GridData.addDataArrayFromFile(flowPlus, 'v1', self.m_InputData.getVelDataName())
                            self.m_GridData.addDataArrayFromFile(flowMinus, 'v0', self.m_InputData.getVelDataName())

                            if self.m_LocalTurbulence:
                                self.m_GridData.removeDataArray('omega0')
                                self.m_GridData.removeDataArray('k0')
                                self.m_GridData.renameDataArray('omega', 'omega0')
                                self.m_GridData.renameDataArray('k',   'k0')

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

                                if self.m_LocalTurbulence:
                                    self.m_GridData.removeDataArray('omega0')
                                    self.m_GridData.removeDataArray('k0')
                                    self.m_GridData.removeDataArray('omega1')
                                    self.m_GridData.removeDataArray('k1')
                                    self.m_GridData.removeDataArray('omega2')
                                    self.m_GridData.removeDataArray('k2')
                                    self.m_GridData.renameDataArray(flowMinus, 'omega0', 'omega')
                                    self.m_GridData.renameDataArray(flowMinus, 'k0',   'k')
                                    self.m_GridData.renameDataArray(flowPlus,  'omega1', 'omega')
                                    self.m_GridData.renameDataArray(flowPlus,  'k1',   'k')
                                    self.m_GridData.renameDataArray(flowPlus,  'omega2', 'omega')
                                    self.m_GridData.renameDataArray(flowPlus,  'k2',   'k')

                                print('Loaded a new velocity data')

                            elif simTime >= tWin_1 - dT:

                                self.m_GridData.removeDataArray('v0')
                                self.m_GridData.removeDataArray('v1')
                                self.m_GridData.removeDataArray('v2')
                                flowPlus    = self.m_InputData.getFlowDataFileName(a_ID = timeWindowDict['ID_Up'][timeIndex])
                                flowMinus   = self.m_InputData.getFlowDataFileName(a_ID = timeWindowDict['ID_Low'][timeIndex])
                                flowNext    = self.m_InputData.getFlowDataFileName(a_ID = timeWindowDict['ID_Low'][timeIndex]+1)
                                self.m_GridData.addDataArrayFromFile(flowPlus, 'v1', self.m_InputData.getVelDataName())
                                self.m_GridData.addDataArrayFromFile(flowMinus, 'v0', self.m_InputData.getVelDataName())
                                self.m_GridData.addDataArrayFromFile(flowNext, 'v2', self.m_InputData.getVelDataName())

                                if self.m_LocalTurbulence:
                                    self.m_GridData.removeDataArray('omega0')
                                    self.m_GridData.removeDataArray('k0')
                                    self.m_GridData.removeDataArray('omega1')
                                    self.m_GridData.removeDataArray('k1')
                                    self.m_GridData.removeDataArray('omega2')
                                    self.m_GridData.removeDataArray('k2')
                                    self.m_GridData.renameDataArray(flowMinus, 'omega0', 'omega')
                                    self.m_GridData.renameDataArray(flowMinus, 'k0',   'k')
                                    self.m_GridData.renameDataArray(flowPlus,  'omega1', 'omega')
                                    self.m_GridData.renameDataArray(flowPlus,  'k1',   'k')
                                    self.m_GridData.renameDataArray(flowPlus,  'omega2', 'omega')
                                    self.m_GridData.renameDataArray(flowPlus,  'k2',   'k')

                                print('Loaded a new velocity data')

            #----------------------------------------------------------------------
            # STEP 5: build a cell locator (if it is a fixed mesh, built only once)60
            #----------------------------------------------------------------------
            if (simTime == self.m_InputData.getSimulationStartTime() and self.m_InputData.isFixedMesh()) or (self.m_ResumeSimulation and self.m_InputData.isFixedMesh()):
                print("Building Cell Locator Maps")
                self.m_GridData.buildLocator()

            #--------------------------------------
            # STEP 6: Inital simulation parameters
            #--------------------------------------
            if simTime == self.m_InputData.getSimulationStartTime() or self.m_ResumeSimulation:
                self.m_LagrangianData.writeData(self.m_InputData.getTracerOutputFile(a_ID1=0, a_ID2=timeIndex))
                if self.m_LocalTurbulence:
                    self.localTurbuleneParams(True)
                self.m_ResumeSimulation = False

                self.m_PrevVels = np.zeros((self.numParticles, 3))
                self.m_PrevXYZ  = np.zeros((self.numParticles, 3))

                for p in range(self.numParticles):
                    self.m_PrevXYZ[p] = self.m_LagrangianData.getX(p)
                    if True:
                        vel, stat = self.m_GridData.gridInterpolateAveraging(self.m_PrevXYZ[p], 'v0', a_GetStatus=True)
                        if self.m_LocalTurbulence:
                            omega = self.m_GridData.gridInterpolateAveraging(self.m_PrevXYZ[p], 'omega0', a_DataDim=1)
                            k = self.m_GridData.gridInterpolateAveraging(self.m_PrevXYZ[p], 'k0', a_DataDim=1)
                    else:
                        vel, stat = self.m_GridData.gridInterpolateNodalBasis(self.m_PrevXYZ[p], 'v0', a_GetStatus=True)
                        if self.m_LocalTurbulence:
                            omega = self.m_GridData.gridInterpolateNodalBasis(self.m_PrevXYZ[p], 'omega0', a_DataDim=1)
                            k = self.m_GridData.gridInterpolateNodalBasis(self.m_PrevXYZ[p], 'k0', a_DataDim=1)

            #------------------------------------------------------------
            # STEP 7: now proceed with integration of each tracer point
            # NOTE: Currently only the basic Euler-Maruyama Integrator is
            # programmed - ore integrators can be added in the future
            #------------------------------------------------------------
            boundaryCondition = 1   # currently a dummy variable

            self.integrateEulerMaruyama(simTime,tWin_0,tWin_1,boundaryCondition,a_DataSync=isDataPointSync, a_PolygonalCells=True)

            #------------------------------------------------
            # STEP 8: update time indices and simulation time
            #------------------------------------------------
            timeIndex += 1
            simTime    = round(simTime+dT, self.m_dt_prec)

            #----------------------------------------------------------------------
            # STEP 9: at the end of appropriate number of steps dump data into file
            #----------------------------------------------------------------------
            if timeIndex%self.m_InputData.m_WriteInterval == 0:
                self.m_LagrangianData.writeData(self.m_InputData.getTracerOutputFile(a_ID1=0, a_ID2=timeIndex))



    #--------------------------------------------------------------------------------------
    # time integration routine for individual particles: Forward Euler Method modified such
    # that this becomes a instanced method of the class instead of being a static method
    #--------------------------------------------------------------------------------------
    def integrateEulerMaruyama(self, a_T, a_T0, a_T1, a_BoundaryCondition, a_DataSync=True, a_PolygonalCells=False):

        dT = self.m_InputData.getIntegrationTimeStep()
        D0 = self.m_InputData.getTracerDiffusivity()

        xyz_new = np.zeros(3,dtype=np.float32)

        for p in range(self.numParticles):

            xyz_old  = self.m_LagrangianData.getX(p)
            xyz_prev = self.m_PrevXYZ[p]
            v_old    = self.m_PrevVels[p]

            if a_DataSync == True:
                if a_PolygonalCells == True:
                    vel, stat = self.m_GridData.gridInterpolateAveraging(xyz_old, 'v0', a_GetStatus=True)
                    if self.m_LocalTurbulence:
                        omega = self.m_GridData.gridInterpolateAveraging(xyz_old, 'omega0', a_DataDim=1)
                        k = self.m_GridData.gridInterpolateAveraging(xyz_old, 'k0', a_DataDim=1)
                else:
                    vel, stat = self.m_GridData.gridInterpolateNodalBasis(xyz_old, 'v0', a_GetStatus=True)
                    if self.m_LocalTurbulence:
                        omega = self.m_GridData.gridInterpolateNodalBasis(xyz_old, 'omega0', a_DataDim=1)
                        k = self.m_GridData.gridInterpolateNodalBasis(xyz_old, 'k0', a_DataDim=1)
            else:
                if a_PolygonalCells == True:
                    velPlus, stat = self.m_GridData.gridInterpolateAveraging(xyz_old, 'v1', a_GetStatus=True)
                    velMinus, stat = self.m_GridData.gridInterpolateAveraging(xyz_old, 'v0', a_GetStatus=True)
                    if self.m_LocalTurbulence:
                        omegaMinus = self.m_GridData.gridInterpolateAveraging(xyz_old, 'omega0', a_DataDim=1)
                        kMinus = self.m_GridData.gridInterpolateAveraging(xyz_old, 'k0', a_DataDim=1)
                        omegaPlus  = self.m_GridData.gridInterpolateAveraging(xyz_old, 'omega1', a_DataDim=1)
                        kPlus = self.m_GridData.gridInterpolateAveraging(xyz_old, 'k1', a_DataDim=1)
                else:
                    velPlus, stat = self.m_GridData.gridInterpolateNodalBasis(xyz_old, 'v1', a_GetStatus=True)
                    velMinus, stat = self.m_GridData.gridInterpolateNodalBasis(xyz_old, 'v0', a_GetStatus=True)
                    if self.m_LocalTurbulence:
                        omegaMinus = self.m_GridData.gridInterpolateNodalBasis(xyz_old, 'omega0', a_DataDim=1)
                        kMinus = self.m_GridData.gridInterpolateNodalBasis(xyz_old, 'k0', a_DataDim=1)
                        omegaPlus  = self.m_GridData.gridInterpolateNodalBasis(xyz_old, 'omega1', a_DataDim=1)
                        kPlus = self.m_GridData.gridInterpolateNodalBasis(xyz_old, 'k1', a_DataDim=1)

                vel = velMinus + (a_T - a_T0) * (velPlus - velMinus)/(a_T1 - a_T0)

                if self.m_LocalTurbulence:
                    omega = omegaMinus + (a_T - a_T0) * (omegaPlus - omegaMinus)/(a_T1 - a_T0)
                    k = kMinus + (a_T - a_T0) * (kPlus - kMinus)/(a_T1 - a_T0)

            #--------------------------------------------------------
            # reflect the particle back into the domain if on a wall
            #--------------------------------------------------------
            if stat == -1:
                vel, wall = self.deflection(xyz_old, v_old)
                if wall:
                    stat = 1
                    xyz_old = xyz_prev
                    if a_PolygonalCells == True:
                        vel, stat = self.m_GridData.gridInterpolateAveraging(xyz_old, 'v0', a_GetStatus=True)
                        if self.m_LocalTurbulence:
                            omega = self.m_GridData.gridInterpolateAveraging(xyz_old, 'omega0', a_DataDim=1)
                            k = self.m_GridData.gridInterpolateAveraging(xyz_old, 'k0', a_DataDim=1)
                    else:
                        vel, stat = self.m_GridData.gridInterpolateNodalBasis(xyz_old, 'v0', a_GetStatus=True)
                        if self.m_LocalTurbulence:
                            omega = self.m_GridData.gridInterpolateNodalBasis(xyz_old, 'omega0', a_DataDim=1)
                            k = self.m_GridData.gridInterpolateNodalBasis(xyz_old, 'k0', a_DataDim=1)

            if self.m_LocalTurbulence and stat != -1:
                for d in range(len(xyz_old)):
                    vel[d] += np.sqrt((2/3)*k) * self.m_pRand[p, d]
                self.m_EddyTime[p] += dT
                if self.m_EddyTime[p] >= self.m_tInt[p]:
                    self.m_tInt[p] = 0.6/omega
                    self.m_EddyTime[p] = 0.0
                    for d in range(len(xyz_old)):
                        self.m_pRand[p,d] = np.random.randn()

            for d in range(len(xyz_old)):
                if stat != -1:
                    xyz_new[d]  = xyz_old[d] + vel[d]*dT + np.sqrt(2.0 * D0 * dT)*np.random.randn()
                else:
                    xyz_new[d]  = xyz_old[d]

            self.m_LagrangianData.setX(p, xyz_new)

            self.m_PrevVels[p] = vel
            self.m_PrevXYZ[p]  = xyz_old
