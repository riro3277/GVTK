#-------------------------------------------------------------------------------
# This module computes the embolus dynamics using a one way coupling scheme based
# on signed distance field algorithm for collision and Maxey Riley Equations for
# fluid forces.
#
# Author:      Akshita Sahni
# Institution: University of Colorado, Boulder
# Last Edit:   Sept. 2022
#-------------------------------------------------------------------------------
import sys, os, vtk
import numpy as np
import time

from numpy import linalg as la

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
        timeWindowDict = self.m_InputData.getDataTimeWindows_v2()

        #-----------------------------------------------------------------------------
        # display the time windows in case the integration involves unsteady flow data
        #-----------------------------------------------------------------------------
        if not(self.m_InputData.isSteadyFlowData()):
            self.m_InputData.printDataTimeWindows()

        #-----------------------------------------------------------------------
        # start time counter
        #-----------------------------------------------------------------------
        timeIndex    = 0
        timeElapsed0 = 0.0
        timeElapsed1 = 0.0
        simTime   = self.m_InputData.getSimulationStartTime()
        stopTime  = self.m_InputData.getSimulationStopTime()
        dT        = self.m_InputData.getIntegrationTimeStep()
        tWin_0    = timeWindowDict['T_Low'][0]
        tWin_1    = timeWindowDict['T_Up'][0]
        t_Next    = timeWindowDict['T_Low'][0]
        data_DT   = self.m_InputData.getFlowDataTimeStep()
        isSingleStageIntegration = True # set true for feu and false for for rk4

        # get signed distance field
        meshFile    = self.m_InputData.getMeshdata()
        surfaceFile = self.m_InputData.getSurfacedata()
        outFile     = str(self.m_InputData.getOutdata())
        self.vtkPolydataDistanceField(meshFile, surfaceFile, outFile)

        # added gradient filter application on existing grid data
        self.m_GridData.getDataFromFile(outFile, a_FileType ='vtu')

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
            simTime += round(dT * timeIndex, 2)
            print('Resuming simulation at', simTime, 's')

        #-----------------------------------------------------------------------
        # start the simulation time loop
        #-----------------------------------------------------------------------
        while simTime <= stopTime:

            #--------------------------------------------------------------------
            # STEP 0: display a status message for the time step being integrated
            #--------------------------------------------------------------------
            # if simTime != 0:
            #     print("Integrating from", timeIndex, "to", timeIndex + 1, "simTime", int(simTime), end='\r')
            # elif round(simTime, 2) == stopTime-dT:
            #     print("Integrating from", timeIndex, "to", timeIndex + 1, "simTime", int(simTime))

            if self.m_InputData.isSteadyFlowData() and simTime >= stopTime:
                print("Integrating from", timeIndex, "to", timeIndex + 1, "simTime", float(simTime), end='\r')
            elif not self.m_InputData.isSteadyFlowData() and (timeWindowDict['T_Low'][timeIndex] != tWin_0 or timeWindowDict['T_Up'][timeIndex] != tWin_1):
                print("Integrating from", timeIndex, "to", timeIndex + 1, "simTime", float(simTime), end='\r')
            elif simTime > 0 and not self.m_ResumeSimulation:
                print("Integrating from", timeIndex, "to", timeIndex + 1, "simTime", float(simTime), end='\r')
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

            #-------------------------------------------------------------------
            # a set of debug messages for program status
            #-------------------------------------------------------------------
            if isLoadFrame:
                print("Will Load Velocity Data")

            if isDataPointSync:
                if self.m_InputData.isSteadyFlowData():
                    if simTime == self.m_InputData.getSimulationStartTime() or self.m_ResumeSimulation:
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
                    t_Next  = timeWindowDict['T_Low'][timeIndex] + 1    # updating the time window

                    if simTime == self.m_InputData.getSimulationStartTime():

                        flowInitial = self.m_InputData.getFlowDataFileName(a_ID = timeWindowDict['ID_Low'][timeIndex])
                        flowNext    = self.m_InputData.getFlowDataFileName(a_ID = timeWindowDict['ID_Low'][timeIndex] + 1)
                        self.m_GridData.addDataArrayFromFile(flowInitial, 'v0',self.m_InputData.getVelDataName())
                        self.m_GridData.addDataArrayFromFile(flowNext, 'v1+',self.m_InputData.getVelDataName())

                        gradVel = self.getVelocityGradient(flowInitial, self.m_InputData.getVelDataName())

                    elif self.m_ResumeSimulation:

                        if isDataPointSync:

                            if isSingleStageIntegration == True:

                                if simTime > self.m_InputData.getSimulationStartTime() and simTime < stopTime - data_DT:

                                    self.m_GridData.removeDataArray('v0')
                                    self.m_GridData.removeDataArray('v1-')
                                    self.m_GridData.removeDataArray('v1+')
                                    flowLive = self.m_InputData.getFlowDataFileName(a_ID = timeWindowDict['ID_Low'][timeIndex])
                                    if timeWindowDict['ID_Low'][timeIndex] != 0:
                                        flowPrev = self.m_InputData.getFlowDataFileName(a_ID = timeWindowDict['ID_Low'][timeIndex] - 1)
                                    else:
                                        flowPrev = self.m_InputData.getFlowDataFileName(a_ID = self.m_InputData.getFlowDataStopFileIndex())
                                    flowNext = self.m_InputData.getFlowDataFileName(a_ID = timeWindowDict['ID_Up'][timeIndex])
                                    self.m_GridData.addDataArrayFromFile(flowLive, 'v0', self.m_InputData.getVelDataName())
                                    self.m_GridData.addDataArrayFromFile(flowPrev, 'v1-', self.m_InputData.getVelDataName())
                                    self.m_GridData.addDataArrayFromFile(flowNext, 'v1+', self.m_InputData.getVelDataName())

                                elif simTime >= stopTime - data_DT :

                                    self.m_GridData.removeDataArray('v0')
                                    self.m_GridData.removeDataArray('v1-')
                                    flowLive = self.m_InputData.getFlowDataFileName(a_ID = timeWindowDict['ID_Low'][timeIndex])
                                    flowPrev = self.m_InputData.getFlowDataFileName(a_ID = timeWindowDict['ID_Low'][timeIndex] - 1)
                                    self.m_GridData.addDataArrayFromFile(flowLive, 'v0', self.m_InputData.getVelDataName())
                                    self.m_GridData.addDataArrayFromFile(flowPrev, 'v1-', self.m_InputData.getVelDataName())

                                gradVel = self.getVelocityGradient(flowLive, self.m_InputData.getVelDataName())

                            else:

                                self.m_GridData.removeDataArray('v0')
                                self.m_GridData.removeDataArray('v1+')
                                flowPlus  = self.m_InputData.getFlowDataFileName(a_ID = timeWindowDict['ID_Up'][timeIndex])
                                flowMinus = self.m_InputData.getFlowDataFileName(a_ID = timeWindowDict['ID_Low'][timeIndex])
                                self.m_GridData.addDataArrayFromFile(flowPlus, 'v1+', self.m_InputData.getVelDataName())
                                self.m_GridData.addDataArrayFromFile(flowMinus, 'v0', self.m_InputData.getVelDataName())
                                print('Loaded a new velocity data')

                        else:

                            if isSingleStageIntegration == True:

                                if simTime > 0 and simTime < data_DT:

                                    self.m_GridData.removeDataArray('v0')
                                    self.m_GridData.removeDataArray('v1+')
                                    flowNext  = self.m_InputData.getFlowDataFileName(a_ID = timeWindowDict['ID_Up'][timeIndex])
                                    flowInitial = self.m_InputData.getFlowDataFileName(a_ID = timeWindowDict['ID_Low'][timeIndex])
                                    self.m_GridData.addDataArrayFromFile(flowNext, 'v1+', self.m_InputData.getVelDataName())
                                    self.m_GridData.addDataArrayFromFile(flowInitial, 'v0', self.m_InputData.getVelDataName())

                                    gradVel = self.getVelocityGradient(flowInitial, self.m_InputData.getVelDataName())
                                    print('Loaded a new velocity data')

                                else:

                                    self.m_GridData.removeDataArray('v0')
                                    self.m_GridData.removeDataArray('v1-')
                                    self.m_GridData.removeDataArray('v1+')
                                    flowPrev = self.m_InputData.getFlowDataFileName(a_ID = timeWindowDict['ID_Low'][timeIndex])
                                    if timeWindowDict['ID_Low'][timeIndex] != 0:
                                        flowLast = self.m_InputData.getFlowDataFileName(a_ID = timeWindowDict['ID_Low'][timeIndex] - 1)
                                    else:
                                        flowLast = self.m_InputData.getFlowDataFileName(a_ID = self.m_InputData.getFlowDataStopFileIndex())
                                    flowNext = self.m_InputData.getFlowDataFileName(a_ID = timeWindowDict['ID_Up'][timeIndex])
                                    self.m_GridData.addDataArrayFromFile(flowPrev, 'v0', self.m_InputData.getVelDataName())
                                    self.m_GridData.addDataArrayFromFile(flowLast, 'v1-', self.m_InputData.getVelDataName())
                                    self.m_GridData.addDataArrayFromFile(flowNext, 'v1+', self.m_InputData.getVelDataName())

                                    gradVel = self.getVelocityGradient(flowPrev, self.m_InputData.getVelDataName())

                            else:

                                if simTime < tWin_1 - dT:

                                    self.m_GridData.removeDataArray('v0')
                                    self.m_GridData.removeDataArray('v1+')
                                    self.m_GridData.removeDataArray('v2')
                                    flowPlus  = self.m_InputData.getFlowDataFileName(a_ID = timeWindowDict['ID_Up'][timeIndex])
                                    flowMinus = self.m_InputData.getFlowDataFileName(a_ID = timeWindowDict['ID_Low'][timeIndex])
                                    self.m_GridData.addDataArrayFromFile(flowPlus, 'v1+', self.m_InputData.getVelDataName())
                                    self.m_GridData.addDataArrayFromFile(flowMinus, 'v0', self.m_InputData.getVelDataName())
                                    print('Loaded a new velocity data')

                                elif simTime >= tWin_1 - dT:

                                    self.m_GridData.removeDataArray('v0')
                                    self.m_GridData.removeDataArray('v1+')
                                    self.m_GridData.removeDataArray('v2')
                                    flowPlus  = self.m_InputData.getFlowDataFileName(a_ID = timeWindowDict['ID_Up'][timeIndex])
                                    flowMinus = self.m_InputData.getFlowDataFileName(a_ID = timeWindowDict['ID_Low'][timeIndex])
                                    flowNext  = self.m_InputData.getFlowDataFileName(a_ID = timeWindowDict['ID_Low'][timeIndex]+1)
                                    self.m_GridData.addDataArrayFromFile(flowPlus, 'v1+', self.m_InputData.getVelDataName())
                                    self.m_GridData.addDataArrayFromFile(flowMinus, 'v0', self.m_InputData.getVelDataName())
                                    self.m_GridData.addDataArrayFromFile(flowNext, 'v2', self.m_InputData.getVelDataName())
                                    print('Loaded a new velocity data')

                        if not self.m_InputData.isFixedMesh():
                            self.m_ResumeSimulation = False

                    elif isDataPointSync:

                        if isSingleStageIntegration == True:

                            if simTime > self.m_InputData.getSimulationStartTime() and simTime < stopTime - data_DT:

                                self.m_GridData.removeDataArray('v0')
                                self.m_GridData.removeDataArray('v1-')
                                self.m_GridData.removeDataArray('v1+')
                                flowLive = self.m_InputData.getFlowDataFileName(a_ID = timeWindowDict['ID_Low'][timeIndex])
                                if timeWindowDict['ID_Low'][timeIndex] != 0:
                                    flowPrev = self.m_InputData.getFlowDataFileName(a_ID = timeWindowDict['ID_Low'][timeIndex] - 1)
                                else:
                                    flowPrev = self.m_InputData.getFlowDataFileName(a_ID = self.m_InputData.getFlowDataStopFileIndex())
                                flowNext = self.m_InputData.getFlowDataFileName(a_ID = timeWindowDict['ID_Up'][timeIndex+1])
                                self.m_GridData.addDataArrayFromFile(flowLive, 'v0', self.m_InputData.getVelDataName())
                                self.m_GridData.addDataArrayFromFile(flowPrev, 'v1-', self.m_InputData.getVelDataName())
                                self.m_GridData.addDataArrayFromFile(flowNext, 'v1+', self.m_InputData.getVelDataName())

                            elif simTime >= stopTime - data_DT :

                                self.m_GridData.removeDataArray('v0')
                                self.m_GridData.removeDataArray('v1-')
                                flowLive = self.m_InputData.getFlowDataFileName(a_ID = timeWindowDict['ID_Low'][timeIndex])
                                flowPrev = self.m_InputData.getFlowDataFileName(a_ID = timeWindowDict['ID_Low'][timeIndex] - 1)
                                self.m_GridData.addDataArrayFromFile(flowLive, 'v0', self.m_InputData.getVelDataName())
                                self.m_GridData.addDataArrayFromFile(flowPrev, 'v1-', self.m_InputData.getVelDataName())

                            gradVel = self.getVelocityGradient(flowLive, self.m_InputData.getVelDataName())

                        else:

                            self.m_GridData.removeDataArray('v0')
                            self.m_GridData.removeDataArray('v1+')
                            flowPlus  = self.m_InputData.getFlowDataFileName(a_ID = timeWindowDict['ID_Up'][timeIndex])
                            flowMinus = self.m_InputData.getFlowDataFileName(a_ID = timeWindowDict['ID_Low'][timeIndex])
                            self.m_GridData.addDataArrayFromFile(flowPlus, 'v1+', self.m_InputData.getVelDataName())
                            self.m_GridData.addDataArrayFromFile(flowMinus, 'v0', self.m_InputData.getVelDataName())
                            print('Loaded a new velocity data')

                    else:

                        if isSingleStageIntegration == True:

                            if simTime > 0 and simTime < data_DT:

                                self.m_GridData.removeDataArray('v0')
                                self.m_GridData.removeDataArray('v1+')
                                flowNext  = self.m_InputData.getFlowDataFileName(a_ID = timeWindowDict['ID_Up'][timeIndex])
                                flowInitial = self.m_InputData.getFlowDataFileName(a_ID = timeWindowDict['ID_Low'][timeIndex])
                                self.m_GridData.addDataArrayFromFile(flowNext, 'v1+', self.m_InputData.getVelDataName())
                                self.m_GridData.addDataArrayFromFile(flowInitial, 'v0', self.m_InputData.getVelDataName())

                                gradVel = self.getVelocityGradient(flowInitial, self.m_InputData.getVelDataName())
                                print('Loaded a new velocity data')

                            else:

                                self.m_GridData.removeDataArray('v0')
                                self.m_GridData.removeDataArray('v1-')
                                self.m_GridData.removeDataArray('v1+')
                                flowPrev = self.m_InputData.getFlowDataFileName(a_ID = timeWindowDict['ID_Low'][timeIndex])
                                if timeWindowDict['ID_Low'][timeIndex] != 0:
                                    flowLast = self.m_InputData.getFlowDataFileName(a_ID = timeWindowDict['ID_Low'][timeIndex] - 1)
                                else:
                                    flowLast = self.m_InputData.getFlowDataFileName(a_ID = self.m_InputData.getFlowDataStopFileIndex())
                                flowNext = self.m_InputData.getFlowDataFileName(a_ID = timeWindowDict['ID_Up'][timeIndex])
                                self.m_GridData.addDataArrayFromFile(flowPrev, 'v0', self.m_InputData.getVelDataName())
                                self.m_GridData.addDataArrayFromFile(flowLast, 'v1-', self.m_InputData.getVelDataName())
                                self.m_GridData.addDataArrayFromFile(flowNext, 'v1+', self.m_InputData.getVelDataName())

                                gradVel = self.getVelocityGradient(flowPrev, self.m_InputData.getVelDataName())

                        else:

                            if simTime < tWin_1 - dT:

                                self.m_GridData.removeDataArray('v0')
                                self.m_GridData.removeDataArray('v1+')
                                self.m_GridData.removeDataArray('v2')
                                flowPlus  = self.m_InputData.getFlowDataFileName(a_ID = timeWindowDict['ID_Up'][timeIndex])
                                flowMinus = self.m_InputData.getFlowDataFileName(a_ID = timeWindowDict['ID_Low'][timeIndex])
                                self.m_GridData.addDataArrayFromFile(flowPlus, 'v1+', self.m_InputData.getVelDataName())
                                self.m_GridData.addDataArrayFromFile(flowMinus, 'v0', self.m_InputData.getVelDataName())
                                print('Loaded a new velocity data')

                            elif simTime >= tWin_1 - dT:

                                self.m_GridData.removeDataArray('v0')
                                self.m_GridData.removeDataArray('v1+')
                                self.m_GridData.removeDataArray('v2')
                                flowPlus  = self.m_InputData.getFlowDataFileName(a_ID = timeWindowDict['ID_Up'][timeIndex])
                                flowMinus = self.m_InputData.getFlowDataFileName(a_ID = timeWindowDict['ID_Low'][timeIndex])
                                flowNext  = self.m_InputData.getFlowDataFileName(a_ID = timeWindowDict['ID_Low'][timeIndex]+1)
                                self.m_GridData.addDataArrayFromFile(flowPlus, 'v1+', self.m_InputData.getVelDataName())
                                self.m_GridData.addDataArrayFromFile(flowMinus, 'v0', self.m_InputData.getVelDataName())
                                self.m_GridData.addDataArrayFromFile(flowNext, 'v2', self.m_InputData.getVelDataName())
                                print('Loaded a new velocity data')

            #-----------------------------------------------------------------------------------------------
            # STEP 5: build a cell locator (if it is a fixed mesh, built only once) and set initial velocity
            #-----------------------------------------------------------------------------------------------
            if simTime == self.m_InputData.getSimulationStartTime() and self.m_InputData.isFixedMesh() or  (self.m_ResumeSimulation and self.m_InputData.isFixedMesh()):
                print("Building Cell Locator Maps")
                t0 = time.process_time()
                if isDataPointSync:
                    locatorObj  = createCellLocatorFromData(self.m_GridData.m_vtkData, a_LocatorType = self.m_InputData.getLocatorType())         # syntax for standard vtkLocator
                else:
                    locatorObj  = createCellLocatorFromData(self.m_GridData.m_vtkData, a_LocatorType = self.m_InputData.getLocatorType())          # syntax for standard vtkLocator
                self.m_GridData.buildLocator(a_LocatorType = self.m_InputData.getLocatorType()) # check use
                t1 = time.process_time()
                timeElapsed0 = t1 - t0
                print('Elapsed time in building cell locator map: ', timeElapsed0)

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

            #Output initial positions of particles
            if simTime == self.m_InputData.getSimulationStartTime() or self.m_ResumeSimulation:
                self.m_LagrangianData.addVectorData('Velocity', InitialVelocity)
                self.m_ResumeSimulation = False
            #----------------------------------------------------------------
            # STEP 7: now proceed with integration of each inertial particle
            #----------------------------------------------------------------
            boundaryCondition = 1   # currently a dummy variable

            # measure time for particle collision
            t2 = time.process_time()
            intScheme = self.m_InputData.getIntegrationScheme()
            if intScheme == 'feu':
                self.sdfCollision(simTime, tWin_0, tWin_1, gradVel, boundaryCondition, locatorObj, a_DataSync=isDataPointSync, a_PolygonalCells=False)
            elif intScheme == 'rk4':
                self.sdfCollision_RK4(simTime, tWin_0, tWin_1, t_Next, gradVel, boundaryCondition, locatorObj, a_DataSync=isDataPointSync, a_PolygonalCells=False)
            t3 = time.process_time()

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
            timeElapsed1 += t3 - t2
            print('elapsed time in particle tracking/collision till last loop: ', timeElapsed1)
        print('Elapsed time in particle tracking/collision: ', timeElapsed1)

    def sdfCollision(self, simTime, tWin_0, tWin_1, gradVel, boundaryCondition, a_Locator, a_DataSync=True, a_PolygonalCells=False):

        dist            = self.m_GridData.extractDataArray('sdf')
        normVec         = self.m_GridData.extractDataArray('sdf-grad')
        gradients       = gradVel.GetPointData().GetArray('gradients')
        dT              = self.m_InputData.getIntegrationTimeStep()
        data_DT         = self.m_InputData.getFlowDataTimeStep()
        restitution     = 0.75
        c_1 = 0.1806
        c_2 = 0.6459
        c_3 = 0.4251
        c_4 = 6880.95
        d_1 = (3.0 * self.m_InputData.getFluidDensity())/(8.0 * self.m_InputData.m_TracerRadius)
        C_am    = 1.0
        accGrav    = np.zeros(3)
        accGrav[0] = 0.0
        accGrav[1] = 0.0
        accGrav[2] = -9.80665 * 1000.0 # mm/sec^2

        #------------------------------------------
        # forward euler integration implementation
        #------------------------------------------
        for p in range(self.numParticles):
            posP         = self.m_LagrangianData.getX(p) # particle position
            v_i          = np.array(self.m_LagrangianData.getVectorData(a_ArrayName='Velocity', a_DataID=p)) # particle velocity
            ## error with using the same name for cell assisgnment: FindCell argument 1: 'tuple' object does not support item assignment
            cell         = a_Locator.FindCell(list(posP))
            cellPtIds    = vtk.vtkIdList()

            delv_delT    = np.zeros(3)

            if simTime == self.m_InputData.getSimulationStartTime():
                if a_PolygonalCells == True:
                    v               = self.m_GridData.gridInterpolateAveraging(posP, 'v0')
                    vNext           = self.m_GridData.gridInterpolateAveraging(posP, 'v1+')

                    delv_delT[0]  = (vNext[0] - v[0])/data_DT
                    delv_delT[1]  = (vNext[1] - v[1])/data_DT
                    delv_delT[2]  = (vNext[2] - v[2])/data_DT

                else:
                    v           = self.m_GridData.gridInterpolateNodalBasis(posP, 'v0')
                    vNext       = self.m_GridData.gridInterpolateNodalBasis(posP, 'v1+')

                    delv_delT[0]  = (vNext[0] - v[0])/data_DT
                    delv_delT[1]  = (vNext[1] - v[1])/data_DT
                    delv_delT[2]  = (vNext[2] - v[2])/data_DT

            elif a_DataSync == True:
                if a_PolygonalCells == True:
                    v           = self.m_GridData.gridInterpolateAveraging(posP, 'v0')
                    vPrev       = self.m_GridData.gridInterpolateAveraging(posP, 'v1-')
                    vNext       = self.m_GridData.gridInterpolateAveraging(posP, 'v1+')

                    if simTime > self.m_InputData.getSimulationStartTime() and simTime < self.m_InputData.getSimulationStopTime() - data_DT:
                        delv_delT[0] = (0.5 * vNext[0] - 0.5 * vPrev[0])/data_DT
                        delv_delT[1] = (0.5 * vNext[1] - 0.5 * vPrev[1])/data_DT
                        delv_delT[2] = (0.5 * vNext[2] - 0.5 * vPrev[2])/data_DT

                    elif simTime >= self.m_InputData.getSimulationStopTime() - data_DT :
                        delv_delT[0] = (v[0] - vPrev[0])/data_DT
                        delv_delT[1] = (v[1] - vPrev[1])/data_DT
                        delv_delT[2] = (v[2] - vPrev[2])/data_DT

                else:
                    v           = self.m_GridData.gridInterpolateNodalBasis(posP, 'v0')
                    vPrev       = self.m_GridData.gridInterpolateNodalBasis(posP, 'v1-')
                    vNext       = self.m_GridData.gridInterpolateNodalBasis(posP, 'v1+')

                    if simTime > self.m_InputData.getSimulationStartTime() and simTime < self.m_InputData.getSimulationStopTime() - data_DT:
                        delv_delT[0] = (0.5 * vNext[0] - 0.5 * vPrev[0])/data_DT
                        delv_delT[1] = (0.5 * vNext[1] - 0.5 * vPrev[1])/data_DT
                        delv_delT[2] = (0.5 * vNext[2] - 0.5 * vPrev[2])/data_DT

                    elif simTime >= self.m_InputData.getSimulationStopTime() - data_DT:
                        delv_delT[0] = (v[0] - vPrev[0])/data_DT
                        delv_delT[1] = (v[1] - vPrev[1])/data_DT
                        delv_delT[2] = (v[2] - vPrev[2])/data_DT

            else:
                if a_PolygonalCells == True:
                    if simTime > 0 and simTime < data_DT:
                        vNext  = self.m_GridData.gridInterpolateAveraging(posP, 'v1+')
                        vPrev  = self.m_GridData.gridInterpolateAveraging(posP, 'v0')

                        delv_delT[0] = (vNext[0] - vPrev[0])/data_DT
                        delv_delT[1] = (vNext[1] - vPrev[1])/data_DT
                        delv_delT[2] = (vNext[2] - vPrev[2])/data_DT

                    else:
                        vNext      = self.m_GridData.gridInterpolateAveraging(posP, 'v1+')
                        vPrev      = self.m_GridData.gridInterpolateAveraging(posP, 'v0')
                        vLast      = self.m_GridData.gridInterpolateAveraging(posP, 'v1-')

                        delv_delT[0] = (0.5 * vNext[0] - 0.5 * vLast[0])/data_DT
                        delv_delT[1] = (0.5 * vNext[1] - 0.5 * vLast[1])/data_DT
                        delv_delT[2] = (0.5 * vNext[2] - 0.5 * vLast[2])/data_DT

                else:
                    if simTime > 0 and simTime < data_DT:
                        vNext      = self.m_GridData.gridInterpolateNodalBasis(posP, 'v1+')
                        vPrev      = self.m_GridData.gridInterpolateNodalBasis(posP, 'v0')

                        delv_delT[0] = (vNext[0] - vPrev[0])/data_DT
                        delv_delT[1] = (vNext[1] - vPrev[1])/data_DT
                        delv_delT[2] = (vNext[2] - vPrev[2])/data_DT

                    else:
                        vNext      = self.m_GridData.gridInterpolateNodalBasis(posP, 'v1+')
                        vPrev      = self.m_GridData.gridInterpolateNodalBasis(posP, 'v0')
                        vLast      = self.m_GridData.gridInterpolateNodalBasis(posP, 'v1-')

                        delv_delT[0] = (0.5 * vNext[0] - 0.5 * vLast[0])/data_DT
                        delv_delT[1] = (0.5 * vNext[1] - 0.5 * vLast[1])/data_DT
                        delv_delT[2] = (0.5 * vNext[2] - 0.5 * vLast[2])/data_DT

                v = vPrev + (simTime - tWin_0) * (vNext - vPrev)/(tWin_1 - tWin_0)

            if cell !=-1:
                self.m_GridData.m_vtkData.GetCellPoints(cell, cellPtIds)

                gradU1 = np.asarray(gradients.GetTuple9(cellPtIds.GetId(0)))
                gradU2 = np.asarray(gradients.GetTuple9(cellPtIds.GetId(1)))
                gradU3 = np.asarray(gradients.GetTuple9(cellPtIds.GetId(2)))
                gradU4 = np.asarray(gradients.GetTuple9(cellPtIds.GetId(3)))

                gradU = 0.25 * (gradU1   + gradU2   + gradU3   + gradU4)

                lhsNS    = np.zeros(3)
                lhsNS[0] = delv_delT[0] + v[0]*gradU[0] + v[1]*gradU[1] + v[2]*gradU[2]
                lhsNS[1] = delv_delT[1] + v[0]*gradU[3] + v[1]*gradU[4] + v[2]*gradU[5]
                lhsNS[2] = delv_delT[2] + v[0]*gradU[6] + v[1]*gradU[7] + v[2]*gradU[8]

                vor      = np.zeros(3)
                vor[0]   = gradU[7] - gradU[5]
                vor[1]   = gradU[2] - gradU[6]
                vor[2]   = gradU[3] - gradU[1]

                Re_p = (self.m_InputData.getTracerDensity() * la.norm(v - v_i) * 2 * self.m_InputData.m_TracerRadius)/self.m_InputData.getFluidViscosity()
                Re_g = self.m_InputData.getFluidDensity() * la.norm(gradU) * ((2.0 * self.m_InputData.m_TracerRadius)**2.0)/self.m_InputData.getFluidViscosity()
                alpha_LSA = 0.5 * (Re_g/(Re_p + 1.0e-16))
                CoeffLift_extra = (3.0/(4.0 * np.pi * self.m_InputData.m_TracerRadius)) * np.sqrt((self.m_InputData.getFluidDensity() * self.m_InputData.getFluidViscosity()) / (la.norm(vor) + 1.0e-16))

                if Re_p == 0:
                    coeffSL = 6.46 * CoeffLift_extra
                elif Re_p <= 40.0:
                    coeffSL = 6.46 * ((1.0 - 0.3314 * np.sqrt(alpha_LSA)) * np.exp(-Re_p/10.0) + 0.3314 * np.sqrt(alpha_LSA)) * CoeffLift_extra
                else:
                    coeffSL = 6.46 * 0.0524 * np.sqrt(0.5 * Re_g) * CoeffLift_extra

                C_D = (24.0/(Re_p + 1.0e-16)) * (1 + c_1*(Re_p + 1.0e-16)**c_2) + c_3/(1 + (c_4/(Re_p + 1.0e-16)))

                drag = d_1 * C_D * la.norm(v - v_i) * (v - v_i)

                shearGradLift = coeffSL * np.cross(v - v_i, vor)

                ambForce = self.m_InputData.getFluidDensity() * (1.0 + (C_am/2.0)) * lhsNS

                bodyForce = (self.m_InputData.getTracerDensity() - self.m_InputData.getFluidDensity()) * accGrav

                v_i = v_i + (1.0/(self.m_InputData.getTracerDensity() + self.m_InputData.getFluidDensity()*(C_am/2.0))) * dT * (drag + shearGradLift + ambForce + bodyForce)

                # dN1 = dist.GetTuple1(cellPtIds.GetId(0))
                # dN2 = dist.GetTuple1(cellPtIds.GetId(1))
                # dN3 = dist.GetTuple1(cellPtIds.GetId(2))
                # dN4 = dist.GetTuple1(cellPtIds.GetId(3))
                # d   = 0.25*(dN1 + dN2 + dN3 + dN4)

                gN1 = np.asarray(normVec.GetTuple3(cellPtIds.GetId(0)))
                gN2 = np.asarray(normVec.GetTuple3(cellPtIds.GetId(1)))
                gN3 = np.asarray(normVec.GetTuple3(cellPtIds.GetId(2)))
                gN4 = np.asarray(normVec.GetTuple3(cellPtIds.GetId(3)))

                d   = self.m_GridData.gridInterpolateNodalBasis(posP,'sdf', a_GetGradient=False, a_GetStatus=False, a_DataDim=1, a_DataType='double')
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
                    vn      = np.dot(v_i,g)/la.norm(g)

                    if vn < 0.0:
                        vt      = v_i - vn*g/la.norm(g)
                        v_i     = vt - restitution*vn*g/la.norm(g)

                posP = posP + v_i*dT

        #--------------------------------
        # set new xyz for  point in loop
        #--------------------------------
            self.m_LagrangianData.setX(p, posP)
            self.m_LagrangianData.setVectorData(v_i, a_ArrayName='Velocity', a_DataID=p)
            # print('setV = ', p, v_i)

    def sdfCollision_RK4(self, simTime, tWin_0, tWin_1, t_Next, gradVel, boundaryCondition, a_Locator, a_DataSync=True, a_PolygonalCells=False):
        dist            = self.m_GridData.extractDataArray('sdf')
        normVec         = self.m_GridData.extractDataArray('sdf-grad')
        gradients       = gradVel.GetPointData().GetArray('gradients')
        dT              = self.m_InputData.getIntegrationTimeStep()
        data_DT         = self.m_InputData.getFlowDataTimeStep()
        restitution     = 0.75
        c_1 = 0.1806
        c_2 = 0.6459
        c_3 = 0.4251
        c_4 = 6880.95
        d_1 = (3.0 * self.m_InputData.getFluidDensity())/(8.0 * self.m_InputData.m_TracerRadius)
        C_am    = 1.0
        accGrav    = np.zeros(3)
        accGrav[0] = 0.0
        accGrav[1] = 0.0
        accGrav[2] = -9.80665 * 1000.0 # mm/sec^2

        #-----------------------------------------------------------------------
        ##Fourth order Runge Kutta implementation
        #-----------------------------------------------------------------------
        for p in range(self.numParticles):
            posP_1       = self.m_LagrangianData.getX(p) # particle position
            v_i          = np.array(self.m_LagrangianData.getVectorData(a_ArrayName='Velocity', a_DataID=p)) # particle velocity
            cell         = a_Locator.FindCell(list(posP_1))
            cellPtIds    = vtk.vtkIdList()

            ##----------------------------------------------------------------------------------------------------
            ## When T_sim = T_start and integration time step dT is small in comparison to data time step data_DT:
            ## then time window remains the same, no need to update flow velocity (flow data from file) w.r.t time
            ## only update flow velocity w.r.t updated particle position
            ## this condition holds for only the first time loop where simtime = simulation start time of the particles
            ## V0 is the closest velocity to the particle
            ## V1+ is the next velocity to the V0
            ##----------------------------------------------------------------------------------------------------

            if simTime == self.m_InputData.getSimulationStartTime():

                    if cell !=-1:
                        self.m_GridData.m_vtkData.GetCellPoints(cell, cellPtIds)

                        if a_PolygonalCells == True:
                            #-------------v_1 is the first lagrangian data position mapped from Eulerian velocity field-----#
                            v_1           = self.m_GridData.gridInterpolateAveraging(posP_1, 'v0') 
                            vNext_1       = self.m_GridData.gridInterpolateAveraging(posP_1, 'v1+') 

                        else :
                            v_1           = self.m_GridData.gridInterpolateNodalBasis(posP_1, 'v0')
                            vNext_1       = self.m_GridData.gridInterpolateNodalBasis(posP_1, 'v1+')

                        ##----------------------------------------------------------------------------------------------------
                        ## Dv/Dt (for flow) = del_v/del_t + v.grad_v
                        ## Dv/Dt            = lhsNS (left hand side of Navier Stokes equation)
                        ## del_v/del_t      = delv_delT
                        ## grad_v           = gradU
                        ##
                        ## vor              = vorticity (curl u)
                        ##----------------------------------------------------------------------------------------------------
                        delv_delT     = np.zeros(3)
                        delv_delT[0]  = (vNext_1[0] - v_1[0])/data_DT
                        delv_delT[1]  = (vNext_1[1] - v_1[1])/data_DT
                        delv_delT[2]  = (vNext_1[2] - v_1[2])/data_DT

                        gradU1        = np.asarray(gradients.GetTuple9(cellPtIds.GetId(0)))
                        gradU2        = np.asarray(gradients.GetTuple9(cellPtIds.GetId(1)))
                        gradU3        = np.asarray(gradients.GetTuple9(cellPtIds.GetId(2)))
                        gradU4        = np.asarray(gradients.GetTuple9(cellPtIds.GetId(3)))

                        gradU         = 0.25 * (gradU1   + gradU2   + gradU3   + gradU4)

                        lhsNS         = np.zeros(3)
                        lhsNS[0]      = delv_delT[0] + v_1[0]*gradU[0] + v_1[1]*gradU[1] + v_1[2]*gradU[2]
                        lhsNS[1]      = delv_delT[1] + v_1[0]*gradU[3] + v_1[1]*gradU[4] + v_1[2]*gradU[5]
                        lhsNS[2]      = delv_delT[2] + v_1[0]*gradU[6] + v_1[1]*gradU[7] + v_1[2]*gradU[8]

                        vor           = np.zeros(3)
                        vor[0]        = gradU[7] - gradU[5]
                        vor[1]        = gradU[2] - gradU[6]
                        vor[2]        = gradU[3] - gradU[1]

                        Re_p          = (self.m_InputData.getTracerDensity() * la.norm(v_1 - v_i) * 2 * self.m_InputData.m_TracerRadius)/self.m_InputData.getFluidViscosity()
                        Re_g          = self.m_InputData.getFluidDensity() * la.norm(gradU) * ((2.0 * self.m_InputData.m_TracerRadius)**2.0)/self.m_InputData.getFluidViscosity()
                        alpha_LSA     = 0.5 * (Re_g/(Re_p + 1.0e-16))
                        C_D           = (24.0/(Re_p + 1.0e-16)) * (1 + c_1*(Re_p + 1.0e-16)**c_2) + c_3/(1 + (c_4/(Re_p + 1.0e-16)))
                        CoeffLift_extra = (3.0/(4.0 * np.pi * self.m_InputData.m_TracerRadius)) * np.sqrt((self.m_InputData.getFluidDensity() * self.m_InputData.getFluidViscosity()) / (la.norm(vor) + 1.0e-16))

                        if Re_p == 0:
                            coeffSL   = 6.46 * CoeffLift_extra
                        elif Re_p <= 40.0:
                            coeffSL   = 6.46 * ((1.0 - 0.3314 * np.sqrt(alpha_LSA)) * np.exp(-Re_p/10.0) + 0.3314 * np.sqrt(alpha_LSA)) * CoeffLift_extra
                        else:
                            coeffSL   = 6.46 * 0.0524 * np.sqrt(0.5 * Re_g) * CoeffLift_extra

                        drag          = d_1 * C_D * la.norm(v_1 - v_i) * (v_1 - v_i)
                        shearGradLift = coeffSL * np.cross(v_1 - v_i, vor)
                        ambForce      = self.m_InputData.getFluidDensity() * (1.0 + (C_am/2.0)) * lhsNS
                        bodyForce     = (self.m_InputData.getTracerDensity() - self.m_InputData.getFluidDensity()) * accGrav

                        ## rk4 scheme
                        k_1           = v_i
                        #--------el_1 is the total forces applied to the particle at posP_1---------#
                        el_1          = (1.0/(self.m_InputData.getTracerDensity() + self.m_InputData.getFluidDensity()*(C_am/2.0))) * (drag + shearGradLift + ambForce + bodyForce)

                        #---------updated particle position for Rk4 integration scheme-------#
                        posP_2        = posP_1 + 0.5*k_1*dT
                        print('Rk4 Particle position 2 (simtime=starttime)=', posP_2)
                        cell          = a_Locator.FindCell(list(posP_2))
                        cellPtIds     = vtk.vtkIdList()

                        if a_PolygonalCells == True:
                            #-------------v_2 is the second lagrangian data position updated using k_1. #
                            #             Mapping from Eulerian velocity field using PosP_2-----#
                            v_2           = self.m_GridData.gridInterpolateAveraging(posP_2, 'v0')
                            vNext_2       = self.m_GridData.gridInterpolateAveraging(posP_2, 'v1+')

                        else :
                            #-------------v_2 is the second lagrangian data position updated using k_1. #
                            #             Mapping from Eulerian velocity field using PosP_2-----#
                            v_2           = self.m_GridData.gridInterpolateNodalBasis(posP_2, 'v0')
                            vNext_2       = self.m_GridData.gridInterpolateNodalBasis(posP_2, 'v1+')

                        if cell !=-1:
                            self.m_GridData.m_vtkData.GetCellPoints(cell, cellPtIds)

                            delv_delT     = np.zeros(3)
                            delv_delT[0]  = (vNext_2[0] - v_2[0])/data_DT
                            delv_delT[1]  = (vNext_2[1] - v_2[1])/data_DT
                            delv_delT[2]  = (vNext_2[2] - v_2[2])/data_DT

                            gradU1        = np.asarray(gradients.GetTuple9(cellPtIds.GetId(0)))
                            gradU2        = np.asarray(gradients.GetTuple9(cellPtIds.GetId(1)))
                            gradU3        = np.asarray(gradients.GetTuple9(cellPtIds.GetId(2)))
                            gradU4        = np.asarray(gradients.GetTuple9(cellPtIds.GetId(3)))

                            gradU         = 0.25 * (gradU1   + gradU2   + gradU3   + gradU4)

                            lhsNS         = np.zeros(3)
                            lhsNS[0]      = delv_delT[0] + v_2[0]*gradU[0] + v_2[1]*gradU[1] + v_2[2]*gradU[2]
                            lhsNS[1]      = delv_delT[1] + v_2[0]*gradU[3] + v_2[1]*gradU[4] + v_2[2]*gradU[5]
                            lhsNS[2]      = delv_delT[2] + v_2[0]*gradU[6] + v_2[1]*gradU[7] + v_2[2]*gradU[8]

                            vor           = np.zeros(3)
                            vor[0]        = gradU[7] - gradU[5]
                            vor[1]        = gradU[2] - gradU[6]
                            vor[2]        = gradU[3] - gradU[1]

                            Re_p          = (self.m_InputData.getTracerDensity() * la.norm(v_2 - v_i - 0.5*el_1*dT) * 2 * self.m_InputData.m_TracerRadius)/self.m_InputData.getFluidViscosity()
                            Re_g          = self.m_InputData.getFluidDensity() * la.norm(gradU) * ((2.0 * self.m_InputData.m_TracerRadius)**2.0)/self.m_InputData.getFluidViscosity()
                            alpha_LSA     = 0.5 * (Re_g/(Re_p + 1.0e-16))
                            C_D           = (24.0/(Re_p + 1.0e-16)) * (1 + c_1*(Re_p + 1.0e-16)**c_2) + c_3/(1 + (c_4/(Re_p + 1.0e-16)))
                            CoeffLift_extra = (3.0/(4.0 * np.pi * self.m_InputData.m_TracerRadius)) * np.sqrt((self.m_InputData.getFluidDensity() * self.m_InputData.getFluidViscosity()) / (la.norm(vor) + 1.0e-16))

                            if Re_p == 0:
                                coeffSL   = 6.46 * CoeffLift_extra
                            elif Re_p <= 40.0:
                                coeffSL   = 6.46 * ((1.0 - 0.3314 * np.sqrt(alpha_LSA)) * np.exp(-Re_p/10.0) + 0.3314 * np.sqrt(alpha_LSA)) * CoeffLift_extra
                            else:
                                coeffSL   = 6.46 * 0.0524 * np.sqrt(0.5 * Re_g) * CoeffLift_extra

                            drag          = d_1 * C_D * la.norm(v_2 - v_i - 0.5*el_1*dT) * (v_2 - v_i - 0.5*el_1*dT)
                            shearGradLift = coeffSL * np.cross(v_2 - v_i - 0.5*el_1*dT, vor)
                            ambForce      = self.m_InputData.getFluidDensity() * (1.0 + (C_am/2.0)) * lhsNS
                            bodyForce     = (self.m_InputData.getTracerDensity() - self.m_InputData.getFluidDensity()) * accGrav

                            ## rk4 scheme
                            k_2           = v_i + 0.5*el_1*dT
                            #--------el_2 is the total forces applied to the particle at posP_2---------#
                            el_2          = (1.0/(self.m_InputData.getTracerDensity() + self.m_InputData.getFluidDensity()*(C_am/2.0))) * (drag + shearGradLift + ambForce + bodyForce)

                            #---------updated particle position for Rk4 integration scheme-------#
                            posP_3        = posP_1 + 0.5*k_2*dT
                            print('Rk4 Particle position 3 (simtime=starttime)=', posP_3)
                            cell          = a_Locator.FindCell(list(posP_3))
                            cellPtIds     = vtk.vtkIdList()

                            #-------------v_3 is the third lagrangian data position updated using k_2. #
                            #             Mapping from Eulerian velocity field using PosP_3-----#
                            if a_PolygonalCells == True:
                                v_3           = self.m_GridData.gridInterpolateAveraging(posP_3, 'v0')
                                vNext_3       = self.m_GridData.gridInterpolateAveraging(posP_3, 'v1+')

                            else :
                                v_3           = self.m_GridData.gridInterpolateNodalBasis(posP_3, 'v0')
                                vNext_3       = self.m_GridData.gridInterpolateNodalBasis(posP_3, 'v1+')

                            if cell !=-1:
                                self.m_GridData.m_vtkData.GetCellPoints(cell, cellPtIds)

                                # delv_delT     = np.zeros(3)
                                delv_delT[0]  = (vNext_3[0] - v_3[0])/data_DT
                                delv_delT[1]  = (vNext_3[1] - v_3[1])/data_DT
                                delv_delT[2]  = (vNext_3[2] - v_3[2])/data_DT

                                gradU1        = np.asarray(gradients.GetTuple9(cellPtIds.GetId(0)))
                                gradU2        = np.asarray(gradients.GetTuple9(cellPtIds.GetId(1)))
                                gradU3        = np.asarray(gradients.GetTuple9(cellPtIds.GetId(2)))
                                gradU4        = np.asarray(gradients.GetTuple9(cellPtIds.GetId(3)))

                                gradU         = 0.25 * (gradU1   + gradU2   + gradU3   + gradU4)

                                lhsNS         = np.zeros(3)
                                lhsNS[0]      = delv_delT[0] + v_3[0]*gradU[0] + v_3[1]*gradU[1] + v_3[2]*gradU[2]
                                lhsNS[1]      = delv_delT[1] + v_3[0]*gradU[3] + v_3[1]*gradU[4] + v_3[2]*gradU[5]
                                lhsNS[2]      = delv_delT[2] + v_3[0]*gradU[6] + v_3[1]*gradU[7] + v_3[2]*gradU[8]

                                vor           = np.zeros(3)
                                vor[0]        = gradU[7] - gradU[5]
                                vor[1]        = gradU[2] - gradU[6]
                                vor[2]        = gradU[3] - gradU[1]

                                Re_p          = (self.m_InputData.getTracerDensity() * la.norm(v_3 - v_i - 0.5*el_2*dT) * 2 * self.m_InputData.m_TracerRadius)/self.m_InputData.getFluidViscosity()
                                Re_g          = self.m_InputData.getFluidDensity() * la.norm(gradU) * ((2.0 * self.m_InputData.m_TracerRadius)**2.0)/self.m_InputData.getFluidViscosity()
                                alpha_LSA     = 0.5 * (Re_g/(Re_p + 1.0e-16))
                                C_D           = (24.0/(Re_p + 1.0e-16)) * (1 + c_1*(Re_p + 1.0e-16)**c_2) + c_3/(1 + (c_4/(Re_p + 1.0e-16)))
                                CoeffLift_extra = (3.0/(4.0 * np.pi * self.m_InputData.m_TracerRadius)) * np.sqrt((self.m_InputData.getFluidDensity() * self.m_InputData.getFluidViscosity()) / (la.norm(vor) + 1.0e-16))

                                if Re_p == 0:
                                    coeffSL     = 6.46 * CoeffLift_extra
                                elif Re_p <= 40.0:
                                    coeffSL     = 6.46 * ((1.0 - 0.3314 * np.sqrt(alpha_LSA)) * np.exp(-Re_p/10.0) + 0.3314 * np.sqrt(alpha_LSA)) * CoeffLift_extra
                                else:
                                    coeffSL     = 6.46 * 0.0524 * np.sqrt(0.5 * Re_g) * CoeffLift_extra

                                drag            = d_1 * C_D * la.norm(v_3 - v_i - 0.5*el_2*dT) * (v_3 - v_i - 0.5*el_2*dT)
                                shearGradLift   = coeffSL * np.cross(v_3 - v_i - 0.5*el_2*dT, vor)
                                ambForce        = self.m_InputData.getFluidDensity() * (1.0 + (C_am/2.0)) * lhsNS
                                bodyForce       = (self.m_InputData.getTracerDensity() - self.m_InputData.getFluidDensity()) * accGrav

                                ## rk4 scheme
                                k_3                 = v_i + 0.5*el_2*dT
                                #--------el_3 is the total forces applied to the particle at posP_3---------#
                                el_3                = (1.0/(self.m_InputData.getTracerDensity() + self.m_InputData.getFluidDensity()*(C_am/2.0))) * (drag + shearGradLift + ambForce + bodyForce)

                                #---------updated particle position for Rk4 integration scheme-------#
                                posP_4              = posP_1 + k_3*dT
                                print('Rk4 Particle position 4 (simtime=starttime)=', posP_4)
                                cell                = a_Locator.FindCell(list(posP_4))
                                cellPtIds           = vtk.vtkIdList()

                                 #-------------v_4 is the 4th lagrangian data position updated using k_3. #
                                 #             Mapping from Eulerian velocity field using PosP_4-----#
                                if a_PolygonalCells == True:
                                    v_4           = self.m_GridData.gridInterpolateAveraging(posP_4, 'v0')
                                    vNext_4       = self.m_GridData.gridInterpolateAveraging(posP_4, 'v1+')

                                else :
                                    v_4           = self.m_GridData.gridInterpolateNodalBasis(posP_4, 'v0')
                                    vNext_4       = self.m_GridData.gridInterpolateNodalBasis(posP_4, 'v1+')

                                if cell !=-1:
                                    self.m_GridData.m_vtkData.GetCellPoints(cell, cellPtIds)

                                    delv_delT     = np.zeros(3)
                                    delv_delT[0]  = (vNext_4[0] - v_4[0])/data_DT
                                    delv_delT[1]  = (vNext_4[1] - v_4[1])/data_DT
                                    delv_delT[2]  = (vNext_4[2] - v_4[2])/data_DT

                                    gradU1        = np.asarray(gradients.GetTuple9(cellPtIds.GetId(0)))
                                    gradU2        = np.asarray(gradients.GetTuple9(cellPtIds.GetId(1)))
                                    gradU3        = np.asarray(gradients.GetTuple9(cellPtIds.GetId(2)))
                                    gradU4        = np.asarray(gradients.GetTuple9(cellPtIds.GetId(3)))

                                    gradU         = 0.25 * (gradU1   + gradU2   + gradU3   + gradU4)

                                    lhsNS         = np.zeros(3)
                                    lhsNS[0]      = delv_delT[0] + v_4[0]*gradU[0] + v_4[1]*gradU[1] + v_4[2]*gradU[2]
                                    lhsNS[1]      = delv_delT[1] + v_4[0]*gradU[3] + v_4[1]*gradU[4] + v_4[2]*gradU[5]
                                    lhsNS[2]      = delv_delT[2] + v_4[0]*gradU[6] + v_4[1]*gradU[7] + v_4[2]*gradU[8]

                                    vor           = np.zeros(3)
                                    vor[0]        = gradU[7] - gradU[5]
                                    vor[1]        = gradU[2] - gradU[6]
                                    vor[2]        = gradU[3] - gradU[1]

                                    Re_p          = (self.m_InputData.getTracerDensity() * la.norm(v_4 - v_i - el_3*dT) * 2 * self.m_InputData.m_TracerRadius)/self.m_InputData.getFluidViscosity()
                                    Re_g          = self.m_InputData.getFluidDensity() * la.norm(gradU) * ((2.0 * self.m_InputData.m_TracerRadius)**2.0)/self.m_InputData.getFluidViscosity()
                                    alpha_LSA     = 0.5 * (Re_g/(Re_p + 1.0e-16))
                                    C_D           = (24.0/(Re_p + 1.0e-16)) * (1 + c_1*(Re_p + 1.0e-16)**c_2) + c_3/(1 + (c_4/(Re_p + 1.0e-16)))
                                    CoeffLift_extra = (3.0/(4.0 * np.pi * self.m_InputData.m_TracerRadius)) * np.sqrt((self.m_InputData.getFluidDensity() * self.m_InputData.getFluidViscosity()) / (la.norm(vor) + 1.0e-16))

                                    if Re_p == 0:
                                        coeffSL     = 6.46 * CoeffLift_extra
                                    elif Re_p <= 40.0:
                                        coeffSL     = 6.46 * ((1.0 - 0.3314 * np.sqrt(alpha_LSA)) * np.exp(-Re_p/10.0) + 0.3314 * np.sqrt(alpha_LSA)) * CoeffLift_extra
                                    else:
                                        coeffSL     = 6.46 * 0.0524 * np.sqrt(0.5 * Re_g) * CoeffLift_extra

                                    drag            = d_1 * C_D * la.norm(v_4 - v_i - el_3*dT) * (v_4 - v_i - el_3*dT)
                                    shearGradLift   = coeffSL * np.cross(v_4 - v_i - el_3*dT, vor)
                                    ambForce        = self.m_InputData.getFluidDensity() * (1.0 + (C_am/2.0)) * lhsNS
                                    bodyForce       = (self.m_InputData.getTracerDensity() - self.m_InputData.getFluidDensity()) * accGrav

                                    ## rk4 scheme
                                    k_4             = v_i + el_3*dT
                                    el_4            = (1.0/(self.m_InputData.getTracerDensity() + self.m_InputData.getFluidDensity()*(C_am/2.0))) * (drag + shearGradLift + ambForce + bodyForce)

                                    # new particle position and velocity after rk4 integration

                                    posP_new        = posP_1 + (dT/6)*(k_1  + 2*k_2  + 2*k_3  + k_4)
                                    v_new           = v_i   + (dT/6)*(el_1 + 2*el_2 + 2*el_3 + el_4)

                                    gN1 = np.asarray(normVec.GetTuple3(cellPtIds.GetId(0)))
                                    gN2 = np.asarray(normVec.GetTuple3(cellPtIds.GetId(1)))
                                    gN3 = np.asarray(normVec.GetTuple3(cellPtIds.GetId(2)))
                                    gN4 = np.asarray(normVec.GetTuple3(cellPtIds.GetId(3)))

                                    if a_PolygonalCells == True:
                                        d   = self.m_GridData.gridInterpolateAveraging(posP_new,'sdf', a_GetGradient=False, a_GetStatus=False, a_DataDim=1, a_DataType='double')
                                    else:
                                        d   = self.m_GridData.gridInterpolateNodalBasis(posP_new,'sdf', a_GetGradient=False, a_GetStatus=False, a_DataDim=1, a_DataType='double')
                                    g   = 0.25*(gN1 + gN2 + gN3 + gN4)

                                    if d <= self.m_InputData.m_TracerRadius:
                                        vn            = np.dot(v_new,g)/la.norm(g)

                                        if vn < 0.0:
                                            vt        = v_new - vn*g/la.norm(g)
                                            v_new     = vt - restitution*vn*g/la.norm(g)

                                    posP_new = posP_new + v_new*dT
                                    print('True Particle position (simtime=starttime)=', posP_new)

            ##----------------------------------------------------------------------------------------------------
            ## When DataSync holds and integration time step dT is small in comparison to data time step data_DT:
            ## then time window remains the same, no need to update flow velocity (flow data from file) w.r.t time
            ## only update flow velocity w.r.t updated particle position
            ## V0 is the closest velocity data to the particle
            ## V1+ is the next velocity data to the V0
            ## V1- is the velocity data previous to V0
            ##----------------------------------------------------------------------------------------------------

            elif a_DataSync == True:

                if cell !=-1:
                    self.m_GridData.m_vtkData.GetCellPoints(cell, cellPtIds)

                    if a_PolygonalCells == True:
                        v_1          = self.m_GridData.gridInterpolateAveraging(posP_1, 'v0')
                        vPrev_1      = self.m_GridData.gridInterpolateAveraging(posP_1, 'v1-')
                        vNext_1      = self.m_GridData.gridInterpolateAveraging(posP_1, 'v1+')
                    else :
                        v_1          = self.m_GridData.gridInterpolateNodalBasis(posP_1, 'v0')
                        vPrev_1      = self.m_GridData.gridInterpolateNodalBasis(posP_1, 'v1-')
                        vNext_1      = self.m_GridData.gridInterpolateNodalBasis(posP_1, 'v1+')

                    if simTime > self.m_InputData.getSimulationStartTime() and simTime < self.m_InputData.getSimulationStopTime() - data_DT:
                        delv_delT    = np.zeros(3)
                        delv_delT[0] = (0.5 * vNext_1[0] - 0.5 * vPrev_1[0])/data_DT
                        delv_delT[1] = (0.5 * vNext_1[1] - 0.5 * vPrev_1[1])/data_DT
                        delv_delT[2] = (0.5 * vNext_1[2] - 0.5 * vPrev_1[2])/data_DT

                    elif simTime >= self.m_InputData.getSimulationStopTime() - data_DT :
                        delv_delT[0] = (v_1[0] - vPrev_1[0])/data_DT
                        delv_delT[1] = (v_1[1] - vPrev_1[1])/data_DT
                        delv_delT[2] = (v_1[2] - vPrev_1[2])/data_DT

                    gradU1 = np.asarray(gradients.GetTuple9(cellPtIds.GetId(0)))
                    gradU2 = np.asarray(gradients.GetTuple9(cellPtIds.GetId(1)))
                    gradU3 = np.asarray(gradients.GetTuple9(cellPtIds.GetId(2)))
                    gradU4 = np.asarray(gradients.GetTuple9(cellPtIds.GetId(3)))

                    gradU  = 0.25 * (gradU1   + gradU2   + gradU3   + gradU4)

                    lhsNS    = np.zeros(3)
                    lhsNS[0] = delv_delT[0] + v_1[0]*gradU[0] + v_1[1]*gradU[1] + v_1[2]*gradU[2]
                    lhsNS[1] = delv_delT[1] + v_1[0]*gradU[3] + v_1[1]*gradU[4] + v_1[2]*gradU[5]
                    lhsNS[2] = delv_delT[2] + v_1[0]*gradU[6] + v_1[1]*gradU[7] + v_1[2]*gradU[8]

                    vor      = np.zeros(3)
                    vor[0]   = gradU[7] - gradU[5]
                    vor[1]   = gradU[2] - gradU[6]
                    vor[2]   = gradU[3] - gradU[1]

                    Re_p          = (self.m_InputData.getTracerDensity() * la.norm(v_1 - v_i) * 2 * self.m_InputData.m_TracerRadius)/self.m_InputData.getFluidViscosity()
                    Re_g          = self.m_InputData.getFluidDensity() * la.norm(gradU) * ((2.0 * self.m_InputData.m_TracerRadius)**2.0)/self.m_InputData.getFluidViscosity()
                    alpha_LSA     = 0.5 * (Re_g/(Re_p + 1.0e-16))
                    C_D           = (24.0/(Re_p + 1.0e-16)) * (1 + c_1*(Re_p + 1.0e-16)**c_2) + c_3/(1 + (c_4/(Re_p + 1.0e-16)))
                    CoeffLift_extra = (3.0/(4.0 * np.pi * self.m_InputData.m_TracerRadius)) * np.sqrt((self.m_InputData.getFluidDensity() * self.m_InputData.getFluidViscosity()) / (la.norm(vor) + 1.0e-16))

                    if Re_p == 0:
                        coeffSL   = 6.46 * CoeffLift_extra
                    elif Re_p <= 40.0:
                        coeffSL   = 6.46 * ((1.0 - 0.3314 * np.sqrt(alpha_LSA)) * np.exp(-Re_p/10.0) + 0.3314 * np.sqrt(alpha_LSA)) * CoeffLift_extra
                    else:
                        coeffSL   = 6.46 * 0.0524 * np.sqrt(0.5 * Re_g) * CoeffLift_extra

                    drag          = d_1 * C_D * la.norm(v_1 - v_i) * (v_1 - v_i)
                    shearGradLift = coeffSL * np.cross(v_1 - v_i, vor)
                    ambForce      = self.m_InputData.getFluidDensity() * (1.0 + (C_am/2.0)) * lhsNS
                    bodyForce     = (self.m_InputData.getTracerDensity() - self.m_InputData.getFluidDensity()) * accGrav

                    ## rk4 scheme
                    k_1           = v_i
                    el_1          = (1.0/(self.m_InputData.getTracerDensity() + self.m_InputData.getFluidDensity()*(C_am/2.0))) * (drag + shearGradLift + ambForce + bodyForce)

                    posP_2        = posP_1 + 0.5*k_1*dT
                    print('Rk4 Particle position 2 (datasync)=', posP_2)
                    cell          = a_Locator.FindCell(list(posP_2))
                    cellPtIds     = vtk.vtkIdList()

                    if cell !=-1:
                        self.m_GridData.m_vtkData.GetCellPoints(cell, cellPtIds)

                        if a_PolygonalCells == True:
                            v_2         = self.m_GridData.gridInterpolateAveraging(posP_2, 'v0')
                            vPrev_2     = self.m_GridData.gridInterpolateAveraging(posP_2, 'v1-')
                            vNext_2     = self.m_GridData.gridInterpolateAveraging(posP_2, 'v1+')
                        else :
                            v_2         = self.m_GridData.gridInterpolateNodalBasis(posP_2, 'v0')
                            vPrev_2     = self.m_GridData.gridInterpolateNodalBasis(posP_2, 'v1-')
                            vNext_2     = self.m_GridData.gridInterpolateNodalBasis(posP_2, 'v1+')

                        if simTime > self.m_InputData.getSimulationStartTime() and simTime < self.m_InputData.getSimulationStopTime() - data_DT:
                            delv_delT[0] = (0.5 * vNext_2[0] - 0.5 * vPrev_2[0])/data_DT
                            delv_delT[1] = (0.5 * vNext_2[1] - 0.5 * vPrev_2[1])/data_DT
                            delv_delT[2] = (0.5 * vNext_2[2] - 0.5 * vPrev_2[2])/data_DT

                        elif simTime >= self.m_InputData.getSimulationStopTime() - data_DT :
                            delv_delT[0] = (v_2[0] - vPrev_2[0])/data_DT
                            delv_delT[1] = (v_2[1] - vPrev_2[1])/data_DT
                            delv_delT[2] = (v_2[2] - vPrev_2[2])/data_DT

                        gradU1 = np.asarray(gradients.GetTuple9(cellPtIds.GetId(0)))
                        gradU2 = np.asarray(gradients.GetTuple9(cellPtIds.GetId(1)))
                        gradU3 = np.asarray(gradients.GetTuple9(cellPtIds.GetId(2)))
                        gradU4 = np.asarray(gradients.GetTuple9(cellPtIds.GetId(3)))

                        gradU  = 0.25 * (gradU1   + gradU2   + gradU3   + gradU4)

                        lhsNS    = np.zeros(3)
                        lhsNS[0] = delv_delT[0] + v_2[0]*gradU[0] + v_2[1]*gradU[1] + v_2[2]*gradU[2]
                        lhsNS[1] = delv_delT[1] + v_2[0]*gradU[3] + v_2[1]*gradU[4] + v_2[2]*gradU[5]
                        lhsNS[2] = delv_delT[2] + v_2[0]*gradU[6] + v_2[1]*gradU[7] + v_2[2]*gradU[8]

                        vor      = np.zeros(3)
                        vor[0]   = gradU[7] - gradU[5]
                        vor[1]   = gradU[2] - gradU[6]
                        vor[2]   = gradU[3] - gradU[1]

                        Re_p          = (self.m_InputData.getTracerDensity() * la.norm(v_2 - v_i - 0.5*k_1*dT) * 2 * self.m_InputData.m_TracerRadius)/self.m_InputData.getFluidViscosity()
                        Re_g          = self.m_InputData.getFluidDensity() * la.norm(gradU) * ((2.0 * self.m_InputData.m_TracerRadius)**2.0)/self.m_InputData.getFluidViscosity()
                        alpha_LSA     = 0.5 * (Re_g/(Re_p + 1.0e-16))
                        C_D           = (24.0/(Re_p + 1.0e-16)) * (1 + c_1*(Re_p + 1.0e-16)**c_2) + c_3/(1 + (c_4/(Re_p + 1.0e-16)))
                        CoeffLift_extra = (3.0/(4.0 * np.pi * self.m_InputData.m_TracerRadius)) * np.sqrt((self.m_InputData.getFluidDensity() * self.m_InputData.getFluidViscosity()) / (la.norm(vor) + 1.0e-16))

                        if Re_p == 0:
                            coeffSL   = 6.46 * CoeffLift_extra
                        elif Re_p <= 40.0:
                            coeffSL   = 6.46 * ((1.0 - 0.3314 * np.sqrt(alpha_LSA)) * np.exp(-Re_p/10.0) + 0.3314 * np.sqrt(alpha_LSA)) * CoeffLift_extra
                        else:
                            coeffSL   = 6.46 * 0.0524 * np.sqrt(0.5 * Re_g) * CoeffLift_extra

                        drag          = d_1 * C_D * la.norm(v_2 - v_i - 0.5*k_1*dT) * (v_2 - v_i - 0.5*k_1*dT)
                        shearGradLift = coeffSL * np.cross(v_2 - v_i - 0.5*k_1*dT, vor)
                        ambForce      = self.m_InputData.getFluidDensity() * (1.0 + (C_am/2.0)) * lhsNS
                        bodyForce     = (self.m_InputData.getTracerDensity() - self.m_InputData.getFluidDensity()) * accGrav

                        ## rk4 scheme
                        k_2           = v_i + 0.5*k_1*dT
                        el_2          = (1.0/(self.m_InputData.getTracerDensity() + self.m_InputData.getFluidDensity()*(C_am/2.0))) * (drag + shearGradLift + ambForce + bodyForce)

                        posP_3        = posP_1 + 0.5*k_2*dT
                        print('Rk4 Particle position 3 (datasync)=', posP_3)
                        cell          = a_Locator.FindCell(list(posP_3))
                        cellPtIds     = vtk.vtkIdList()

                        if cell !=-1:
                            self.m_GridData.m_vtkData.GetCellPoints(cell, cellPtIds)

                            if a_PolygonalCells == True:
                                v_3         = self.m_GridData.gridInterpolateAveraging(posP_3, 'v0')
                                vPrev_3     = self.m_GridData.gridInterpolateAveraging(posP_3, 'v1-')
                                vNext_3     = self.m_GridData.gridInterpolateAveraging(posP_3, 'v1+')
                            else :
                                v_3         = self.m_GridData.gridInterpolateNodalBasis(posP_3, 'v0')
                                vPrev_3     = self.m_GridData.gridInterpolateNodalBasis(posP_3, 'v1-')
                                vNext_3     = self.m_GridData.gridInterpolateNodalBasis(posP_3, 'v1+')

                            if simTime > self.m_InputData.getSimulationStartTime() and simTime < self.m_InputData.getSimulationStopTime() - data_DT:
                                delv_delT[0] = (0.5 * vNext_3[0] - 0.5 * vPrev_3[0])/data_DT
                                delv_delT[1] = (0.5 * vNext_3[1] - 0.5 * vPrev_3[1])/data_DT
                                delv_delT[2] = (0.5 * vNext_3[2] - 0.5 * vPrev_3[2])/data_DT

                            elif simTime >= self.m_InputData.getSimulationStopTime() - data_DT :
                                delv_delT[0] = (v_3[0] - vPrev_3[0])/data_DT
                                delv_delT[1] = (v_3[1] - vPrev_3[1])/data_DT
                                delv_delT[2] = (v_3[2] - vPrev_3[2])/data_DT

                            gradU1 = np.asarray(gradients.GetTuple9(cellPtIds.GetId(0)))
                            gradU2 = np.asarray(gradients.GetTuple9(cellPtIds.GetId(1)))
                            gradU3 = np.asarray(gradients.GetTuple9(cellPtIds.GetId(2)))
                            gradU4 = np.asarray(gradients.GetTuple9(cellPtIds.GetId(3)))

                            gradU  = 0.25 * (gradU1   + gradU2   + gradU3   + gradU4)

                            lhsNS    = np.zeros(3)
                            lhsNS[0] = delv_delT[0] + v_3[0]*gradU[0] + v_3[1]*gradU[1] + v_3[2]*gradU[2]
                            lhsNS[1] = delv_delT[1] + v_3[0]*gradU[3] + v_3[1]*gradU[4] + v_3[2]*gradU[5]
                            lhsNS[2] = delv_delT[2] + v_3[0]*gradU[6] + v_3[1]*gradU[7] + v_3[2]*gradU[8]

                            vor      = np.zeros(3)
                            vor[0]   = gradU[7] - gradU[5]
                            vor[1]   = gradU[2] - gradU[6]
                            vor[2]   = gradU[3] - gradU[1]

                            Re_p          = (self.m_InputData.getTracerDensity() * la.norm(v_3 - v_i - 0.5*k_2*dT) * 2 * self.m_InputData.m_TracerRadius)/self.m_InputData.getFluidViscosity()
                            Re_g          = self.m_InputData.getFluidDensity() * la.norm(gradU) * ((2.0 * self.m_InputData.m_TracerRadius)**2.0)/self.m_InputData.getFluidViscosity()
                            alpha_LSA     = 0.5 * (Re_g/(Re_p + 1.0e-16))
                            C_D           = (24.0/(Re_p + 1.0e-16)) * (1 + c_1*(Re_p + 1.0e-16)**c_2) + c_3/(1 + (c_4/(Re_p + 1.0e-16)))
                            CoeffLift_extra = (3.0/(4.0 * np.pi * self.m_InputData.m_TracerRadius)) * np.sqrt((self.m_InputData.getFluidDensity() * self.m_InputData.getFluidViscosity()) / (la.norm(vor) + 1.0e-16))

                            if Re_p == 0:
                                coeffSL   = 6.46 * CoeffLift_extra
                            elif Re_p <= 40.0:
                                coeffSL   = 6.46 * ((1.0 - 0.3314 * np.sqrt(alpha_LSA)) * np.exp(-Re_p/10.0) + 0.3314 * np.sqrt(alpha_LSA)) * CoeffLift_extra
                            else:
                                coeffSL   = 6.46 * 0.0524 * np.sqrt(0.5 * Re_g) * CoeffLift_extra

                            drag          = d_1 * C_D * la.norm(v_3 - v_i - 0.5*k_2*dT) * (v_3 - v_i - 0.5*k_2*dT)
                            shearGradLift = coeffSL * np.cross(v_3 - v_i - 0.5*k_2*dT, vor)
                            ambForce      = self.m_InputData.getFluidDensity() * (1.0 + (C_am/2.0)) * lhsNS
                            bodyForce     = (self.m_InputData.getTracerDensity() - self.m_InputData.getFluidDensity()) * accGrav

                            ## rk4 scheme
                            k_3           = v_i + 0.5*k_2*dT
                            el_3          = (1.0/(self.m_InputData.getTracerDensity() + self.m_InputData.getFluidDensity()*(C_am/2.0))) * (drag + shearGradLift + ambForce + bodyForce)

                            posP_4        = posP_1 + k_3*dT
                            print('Rk4 Particle position 4 (datasync)=', posP_4)
                            cell          = a_Locator.FindCell(list(posP_4))
                            cellPtIds     = vtk.vtkIdList()

                            if cell !=-1:
                                self.m_GridData.m_vtkData.GetCellPoints(cell, cellPtIds)

                                if a_PolygonalCells == True:
                                    v_4         = self.m_GridData.gridInterpolateAveraging(posP_4, 'v0')
                                    vPrev_4     = self.m_GridData.gridInterpolateAveraging(posP_4, 'v1-')
                                    vNext_4     = self.m_GridData.gridInterpolateAveraging(posP_4, 'v1+')
                                else :
                                    v_4         = self.m_GridData.gridInterpolateNodalBasis(posP_4, 'v0')
                                    vPrev_4     = self.m_GridData.gridInterpolateNodalBasis(posP_4, 'v1-')
                                    vNext_4     = self.m_GridData.gridInterpolateNodalBasis(posP_4, 'v1+')

                                if simTime > self.m_InputData.getSimulationStartTime() and simTime < self.m_InputData.getSimulationStopTime() - data_DT:
                                    delv_delT[0] = (0.5 * vNext_4[0] - 0.5 * vPrev_4[0])/data_DT
                                    delv_delT[1] = (0.5 * vNext_4[1] - 0.5 * vPrev_4[1])/data_DT
                                    delv_delT[2] = (0.5 * vNext_4[2] - 0.5 * vPrev_4[2])/data_DT

                                elif simTime >= self.m_InputData.getSimulationStopTime() - data_DT :
                                    delv_delT[0] = (v_4[0] - vPrev_4[0])/data_DT
                                    delv_delT[1] = (v_4[1] - vPrev_4[1])/data_DT
                                    delv_delT[2] = (v_4[2] - vPrev_4[2])/data_DT

                                gradU1 = np.asarray(gradients.GetTuple9(cellPtIds.GetId(0)))
                                gradU2 = np.asarray(gradients.GetTuple9(cellPtIds.GetId(1)))
                                gradU3 = np.asarray(gradients.GetTuple9(cellPtIds.GetId(2)))
                                gradU4 = np.asarray(gradients.GetTuple9(cellPtIds.GetId(3)))

                                gradU  = 0.25 * (gradU1   + gradU2   + gradU3   + gradU4)

                                lhsNS    = np.zeros(3)
                                lhsNS[0] = delv_delT[0] + v_4[0]*gradU[0] + v_4[1]*gradU[1] + v_4[2]*gradU[2]
                                lhsNS[1] = delv_delT[1] + v_4[0]*gradU[3] + v_4[1]*gradU[4] + v_4[2]*gradU[5]
                                lhsNS[2] = delv_delT[2] + v_4[0]*gradU[6] + v_4[1]*gradU[7] + v_4[2]*gradU[8]

                                vor      = np.zeros(3)
                                vor[0]   = gradU[7] - gradU[5]
                                vor[1]   = gradU[2] - gradU[6]
                                vor[2]   = gradU[3] - gradU[1]

                                Re_p          = (self.m_InputData.getTracerDensity() * la.norm(v_4 - v_i - k_3*dT) * 2 * self.m_InputData.m_TracerRadius)/self.m_InputData.getFluidViscosity()
                                Re_g          = self.m_InputData.getFluidDensity() * la.norm(gradU) * ((2.0 * self.m_InputData.m_TracerRadius)**2.0)/self.m_InputData.getFluidViscosity()
                                alpha_LSA     = 0.5 * (Re_g/(Re_p + 1.0e-16))
                                C_D           = (24.0/(Re_p + 1.0e-16)) * (1 + c_1*(Re_p + 1.0e-16)**c_2) + c_3/(1 + (c_4/(Re_p + 1.0e-16)))
                                CoeffLift_extra = (3.0/(4.0 * np.pi * self.m_InputData.m_TracerRadius)) * np.sqrt((self.m_InputData.getFluidDensity() * self.m_InputData.getFluidViscosity()) / (la.norm(vor) + 1.0e-16))

                                if Re_p == 0:
                                    coeffSL   = 6.46 * CoeffLift_extra
                                elif Re_p <= 40.0:
                                    coeffSL   = 6.46 * ((1.0 - 0.3314 * np.sqrt(alpha_LSA)) * np.exp(-Re_p/10.0) + 0.3314 * np.sqrt(alpha_LSA)) * CoeffLift_extra
                                else:
                                    coeffSL   = 6.46 * 0.0524 * np.sqrt(0.5 * Re_g) * CoeffLift_extra

                                drag          = d_1 * C_D * la.norm(v_4 - v_i - k_3*dT) * (v_4 - v_i - k_3*dT)
                                shearGradLift = coeffSL * np.cross(v_4 - v_i - k_3*dT, vor)
                                ambForce      = self.m_InputData.getFluidDensity() * (1.0 + (C_am/2.0)) * lhsNS
                                bodyForce     = (self.m_InputData.getTracerDensity() - self.m_InputData.getFluidDensity()) * accGrav

                                ## rk4 scheme
                                k_4           = v_i + k_3*dT
                                el_4          = (1.0/(self.m_InputData.getTracerDensity() + self.m_InputData.getFluidDensity()*(C_am/2.0))) * (drag + shearGradLift + ambForce + bodyForce)

                                # new particle position and velocity after rk4 integration

                                posP_new        = posP_1 + (dT/6)*(k_1  + 2*k_2  + 2*k_3  + k_4)
                                v_new           = v_i   + (dT/6)*(el_1 + 2*el_2 + 2*el_3 + el_4)

                                gN1 = np.asarray(normVec.GetTuple3(cellPtIds.GetId(0)))
                                gN2 = np.asarray(normVec.GetTuple3(cellPtIds.GetId(1)))
                                gN3 = np.asarray(normVec.GetTuple3(cellPtIds.GetId(2)))
                                gN4 = np.asarray(normVec.GetTuple3(cellPtIds.GetId(3)))

                                if a_PolygonalCells == True:
                                    d   = self.m_GridData.gridInterpolateAveraging(posP_new,'sdf', a_GetGradient=False, a_GetStatus=False, a_DataDim=1, a_DataType='double')
                                else:
                                    d   = self.m_GridData.gridInterpolateNodalBasis(posP_new,'sdf', a_GetGradient=False, a_GetStatus=False, a_DataDim=1, a_DataType='double')
                                g   = 0.25*(gN1 + gN2 + gN3 + gN4)

                                if d <= self.m_InputData.m_TracerRadius:
                                    vn            = np.dot(v_new,g)/la.norm(g)

                                    if vn < 0.0:
                                        vt        = v_new - vn*g/la.norm(g)
                                        v_new     = vt - restitution*vn*g/la.norm(g)

                                posP_new = posP_new + v_new*dT
                                print('True Particle position (datasync)=', posP_new)

            ##-----------------------------------------------------------------------------------------------------------
            ## When DataSync doesn't hold:
            ## then time window can change with each increment in time step, flow velocity needs to be updated w.r.t time
            ## and updated particle position
            ## V0 is the closest velocity data to the particle
            ## V1+ is the next velocity data to the V0
            ## V1- is the velocity data previous to V0
            ## V2 is the next velocity data to V1+, taken when Rk4 integration time is greater than timewindow (>tWin_1)
            ## t_2 is the second Rk4 integration time for the particle location 
            ## t_3 is the third Rk4 integration time for the particle location
            ## t_4 is the 4th Rk4 integration time for the particle location
            ##-----------------------------------------------------------------------------------------------------------

            else:
                simPrec = self.m_InputData.getSimPrec()
                t_2 = round(simTime + 0.5*dT, simPrec+1)
                t_3 = round(simTime + 0.5*dT, simPrec+1)
                t_4 = round(simTime + dT, simPrec+1)

                if cell !=-1:
                    self.m_GridData.m_vtkData.GetCellPoints(cell, cellPtIds)

                    if simTime > 0 and simTime < data_DT:

                        if a_PolygonalCells == True:
                            vNext_1  = self.m_GridData.gridInterpolateAveraging(posP_1, 'v1+')
                            vPrev_1  = self.m_GridData.gridInterpolateAveraging(posP_1, 'v0')
                        else:
                            vNext_1  = self.m_GridData.gridInterpolateNodalBasis(posP_1, 'v1+')
                            vPrev_1  = self.m_GridData.gridInterpolateNodalBasis(posP_1, 'v0')
                            
                        
                        delv_delT    = np.zeros(3)
                        delv_delT[0] = (vNext_1[0] - vPrev_1[0])/data_DT
                        delv_delT[1] = (vNext_1[1] - vPrev_1[1])/data_DT
                        delv_delT[2] = (vNext_1[2] - vPrev_1[2])/data_DT

                    #-------------------------------------------------------------------------------#
                    #  when the particle is positioned such as data sync is not true and 
                    #  greater than velocity delta time step i.e, the first velocity time step.
                    #-------------------------------------------------------------------------------#
                    else:

                        if a_PolygonalCells == True:
                            vNext_1      = self.m_GridData.gridInterpolateAveraging(posP_1, 'v1+')
                            vPrev_1      = self.m_GridData.gridInterpolateAveraging(posP_1, 'v0')
                            vLast_1      = self.m_GridData.gridInterpolateAveraging(posP_1, 'v1-')

                        else:
                            vNext_1      = self.m_GridData.gridInterpolateNodalBasis(posP_1, 'v1+')
                            vPrev_1      = self.m_GridData.gridInterpolateNodalBasis(posP_1, 'v0')
                            vLast_1      = self.m_GridData.gridInterpolateNodalBasis(posP_1, 'v1-')

                        delv_delT    = np.zeros(3)
                        delv_delT[0] = (0.5 * vNext_1[0] - 0.5 * vLast_1[0])/data_DT
                        delv_delT[1] = (0.5 * vNext_1[1] - 0.5 * vLast_1[1])/data_DT
                        delv_delT[2] = (0.5 * vNext_1[2] - 0.5 * vLast_1[2])/data_DT

                    v_1 = vPrev_1 + (simTime - tWin_0) * (vNext_1 - vPrev_1)/(tWin_1 - tWin_0)

                    gradU1 = np.asarray(gradients.GetTuple9(cellPtIds.GetId(0)))
                    gradU2 = np.asarray(gradients.GetTuple9(cellPtIds.GetId(1)))
                    gradU3 = np.asarray(gradients.GetTuple9(cellPtIds.GetId(2)))
                    gradU4 = np.asarray(gradients.GetTuple9(cellPtIds.GetId(3)))

                    gradU = 0.25 * (gradU1   + gradU2   + gradU3   + gradU4)

                    lhsNS    = np.zeros(3)
                    lhsNS[0] = delv_delT[0] + v_1[0]*gradU[0] + v_1[1]*gradU[1] + v_1[2]*gradU[2]
                    lhsNS[1] = delv_delT[1] + v_1[0]*gradU[3] + v_1[1]*gradU[4] + v_1[2]*gradU[5]
                    lhsNS[2] = delv_delT[2] + v_1[0]*gradU[6] + v_1[1]*gradU[7] + v_1[2]*gradU[8]

                    vor      = np.zeros(3)
                    vor[0]   = gradU[7] - gradU[5]
                    vor[1]   = gradU[2] - gradU[6]
                    vor[2]   = gradU[3] - gradU[1]

                    Re_p = (self.m_InputData.getTracerDensity() * la.norm(v_1 - v_i) * 2 * self.m_InputData.m_TracerRadius)/self.m_InputData.getFluidViscosity()
                    Re_g = self.m_InputData.getFluidDensity() * la.norm(gradU) * ((2.0 * self.m_InputData.m_TracerRadius)**2.0)/self.m_InputData.getFluidViscosity()
                    alpha_LSA = 0.5 * (Re_g/(Re_p + 1.0e-16))
                    C_D = (24.0/(Re_p + 1.0e-16)) * (1 + c_1*(Re_p + 1.0e-16)**c_2) + c_3/(1 + (c_4/(Re_p + 1.0e-16)))
                    CoeffLift_extra = (3.0/(4.0 * np.pi * self.m_InputData.m_TracerRadius)) * np.sqrt((self.m_InputData.getFluidDensity() * self.m_InputData.getFluidViscosity()) / (la.norm(vor) + 1.0e-16))

                    if Re_p == 0:
                        coeffSL = 6.46 * CoeffLift_extra
                    elif Re_p <= 40.0:
                        coeffSL = 6.46 * ((1.0 - 0.3314 * np.sqrt(alpha_LSA)) * np.exp(-Re_p/10.0) + 0.3314 * np.sqrt(alpha_LSA)) * CoeffLift_extra
                    else:
                        coeffSL = 6.46 * 0.0524 * np.sqrt(0.5 * Re_g) * CoeffLift_extra

                    drag = d_1 * C_D * la.norm(v_1 - v_i) * (v_1 - v_i)
                    shearGradLift = coeffSL * np.cross(v_1 - v_i, vor)
                    ambForce = self.m_InputData.getFluidDensity() * (1.0 + (C_am/2.0)) * lhsNS
                    bodyForce = (self.m_InputData.getTracerDensity() - self.m_InputData.getFluidDensity()) * accGrav

                    ## rk4 scheme step 1
                    k_1           = v_i
                    el_1          = (1.0/(self.m_InputData.getTracerDensity() + self.m_InputData.getFluidDensity()*(C_am/2.0))) * (drag + shearGradLift + ambForce + bodyForce)

                    posP_2        = posP_1 + 0.5*k_1*dT
                    print('Rk4 Particle position 2 (no datasync)=', posP_2)
                    cell          = a_Locator.FindCell(list(posP_2))
                    cellPtIds     = vtk.vtkIdList()

                    if cell !=-1:
                        self.m_GridData.m_vtkData.GetCellPoints(cell, cellPtIds)

                        if t_2 <= tWin_1:

                            if simTime > 0 and simTime < data_DT:

                                if a_PolygonalCells == True:
                                    vNext_2  = self.m_GridData.gridInterpolateAveraging(posP_2, 'v1+')
                                    vPrev_2  = self.m_GridData.gridInterpolateAveraging(posP_2, 'v0')
                                else:
                                    vNext_2  = self.m_GridData.gridInterpolateNodalBasis(posP_2, 'v1+')
                                    vPrev_2  = self.m_GridData.gridInterpolateNodalBasis(posP_2, 'v0')

                                delv_delT    = np.zeros(3)
                                delv_delT[0] = (vNext_2[0] - vPrev_2[0])/data_DT
                                delv_delT[1] = (vNext_2[1] - vPrev_2[1])/data_DT
                                delv_delT[2] = (vNext_2[2] - vPrev_2[2])/data_DT

                            else:

                                if a_PolygonalCells == True:
                                    vNext_2      = self.m_GridData.gridInterpolateAveraging(posP_2, 'v1+')
                                    vPrev_2      = self.m_GridData.gridInterpolateAveraging(posP_2, 'v0')
                                    vLast_2      = self.m_GridData.gridInterpolateAveraging(posP_2, 'v1-')

                                else:
                                    vNext_2      = self.m_GridData.gridInterpolateNodalBasis(posP_2, 'v1+')
                                    vPrev_2      = self.m_GridData.gridInterpolateNodalBasis(posP_2, 'v0')
                                    vLast_2      = self.m_GridData.gridInterpolateNodalBasis(posP_2, 'v1-')

                                delv_delT    = np.zeros(3)
                                delv_delT[0] = (0.5 * vNext_2[0] - 0.5 * vLast_2[0])/data_DT
                                delv_delT[1] = (0.5 * vNext_2[1] - 0.5 * vLast_2[1])/data_DT
                                delv_delT[2] = (0.5 * vNext_2[2] - 0.5 * vLast_2[2])/data_DT

                            v_2 = vPrev_2 + (simTime - tWin_0) * (vNext_2 - vPrev_2)/(tWin_1 - tWin_0)

                        elif t_2 > tWin_1:

                            if a_PolygonalCells == True:
                                vNext_2      = self.m_GridData.gridInterpolateAveraging(posP_2, 'v2')
                                vPrev_2      = self.m_GridData.gridInterpolateAveraging(posP_2, 'v1+')
                                vLast_2      = self.m_GridData.gridInterpolateAveraging(posP_2, 'v0')

                            else:
                                vNext_2      = self.m_GridData.gridInterpolateNodalBasis(posP_2, 'v2')
                                vPrev_2      = self.m_GridData.gridInterpolateNodalBasis(posP_2, 'v1+')
                                vLast_2      = self.m_GridData.gridInterpolateNodalBasis(posP_2, 'v0')

                            delv_delT    = np.zeros(3)
                            delv_delT[0] = (0.5 * vNext_2[0] - 0.5 * vLast_2[0])/data_DT
                            delv_delT[1] = (0.5 * vNext_2[1] - 0.5 * vLast_2[1])/data_DT
                            delv_delT[2] = (0.5 * vNext_2[2] - 0.5 * vLast_2[2])/data_DT

                            v_2 = vPrev_2 + (simTime - tWin_0) * (vNext_2 - vPrev_2)/(t_Next - tWin_1)

                        gradU1 = np.asarray(gradients.GetTuple9(cellPtIds.GetId(0)))
                        gradU2 = np.asarray(gradients.GetTuple9(cellPtIds.GetId(1)))
                        gradU3 = np.asarray(gradients.GetTuple9(cellPtIds.GetId(2)))
                        gradU4 = np.asarray(gradients.GetTuple9(cellPtIds.GetId(3)))

                        gradU = 0.25 * (gradU1   + gradU2   + gradU3   + gradU4)

                        lhsNS    = np.zeros(3)
                        lhsNS[0] = delv_delT[0] + v_2[0]*gradU[0] + v_2[1]*gradU[1] + v_2[2]*gradU[2]
                        lhsNS[1] = delv_delT[1] + v_2[0]*gradU[3] + v_2[1]*gradU[4] + v_2[2]*gradU[5]
                        lhsNS[2] = delv_delT[2] + v_2[0]*gradU[6] + v_2[1]*gradU[7] + v_2[2]*gradU[8]

                        vor      = np.zeros(3)
                        vor[0]   = gradU[7] - gradU[5]
                        vor[1]   = gradU[2] - gradU[6]
                        vor[2]   = gradU[3] - gradU[1]

                        Re_p = (self.m_InputData.getTracerDensity() * la.norm(v_2 - v_i - 0.5*el_1*dT) * 2 * self.m_InputData.m_TracerRadius)/self.m_InputData.getFluidViscosity()
                        Re_g = self.m_InputData.getFluidDensity() * la.norm(gradU) * ((2.0 * self.m_InputData.m_TracerRadius)**2.0)/self.m_InputData.getFluidViscosity()
                        alpha_LSA = 0.5 * (Re_g/(Re_p + 1.0e-16))
                        C_D = (24.0/(Re_p + 1.0e-16)) * (1 + c_1*(Re_p + 1.0e-16)**c_2) + c_3/(1 + (c_4/(Re_p + 1.0e-16)))
                        CoeffLift_extra = (3.0/(4.0 * np.pi * self.m_InputData.m_TracerRadius)) * np.sqrt((self.m_InputData.getFluidDensity() * self.m_InputData.getFluidViscosity()) / (la.norm(vor) + 1.0e-16))

                        if Re_p == 0:
                            coeffSL = 6.46 * CoeffLift_extra
                        elif Re_p <= 40.0:
                            coeffSL = 6.46 * ((1.0 - 0.3314 * np.sqrt(alpha_LSA)) * np.exp(-Re_p/10.0) + 0.3314 * np.sqrt(alpha_LSA)) * CoeffLift_extra
                        else:
                            coeffSL = 6.46 * 0.0524 * np.sqrt(0.5 * Re_g) * CoeffLift_extra

                        drag = d_1 * C_D * la.norm(v_2 - v_i  - 0.5*el_1*dT) * (v_2 - v_i  - 0.5*el_1*dT)
                        shearGradLift = coeffSL * np.cross(v_2 - v_i  - 0.5*el_1*dT, vor)
                        ambForce = self.m_InputData.getFluidDensity() * (1.0 + (C_am/2.0)) * lhsNS
                        bodyForce = (self.m_InputData.getTracerDensity() - self.m_InputData.getFluidDensity()) * accGrav

                        ## rk4 scheme step 2
                        k_2           = v_i + 0.5*el_1*dT
                        el_2          = (1.0/(self.m_InputData.getTracerDensity() + self.m_InputData.getFluidDensity()*(C_am/2.0))) * (drag + shearGradLift + ambForce + bodyForce)

                        posP_3        = posP_1 + 0.5*k_2*dT
                        print('Rk4 Particle position 3 (no datasync)=', posP_3)
                        cell          = a_Locator.FindCell(list(posP_3))
                        cellPtIds     = vtk.vtkIdList()

                        if cell !=-1:
                            self.m_GridData.m_vtkData.GetCellPoints(cell, cellPtIds)

                            if t_3 <= tWin_1:

                                if simTime > 0 and simTime < data_DT:

                                    if a_PolygonalCells == True:
                                        vNext_3  = self.m_GridData.gridInterpolateAveraging(posP_3, 'v1+')
                                        vPrev_3  = self.m_GridData.gridInterpolateAveraging(posP_3, 'v0')
                                    else:
                                        vNext_3  = self.m_GridData.gridInterpolateNodalBasis(posP_3, 'v1+')
                                        vPrev_3  = self.m_GridData.gridInterpolateNodalBasis(posP_3, 'v0')

                                    delv_delT    = np.zeros(3)
                                    delv_delT[0] = (vNext_3[0] - vPrev_3[0])/data_DT
                                    delv_delT[1] = (vNext_3[1] - vPrev_3[1])/data_DT
                                    delv_delT[2] = (vNext_3[2] - vPrev_3[2])/data_DT

                                else:

                                    if a_PolygonalCells == True:
                                        vNext_3      = self.m_GridData.gridInterpolateAveraging(posP_3, 'v1+')
                                        vPrev_3      = self.m_GridData.gridInterpolateAveraging(posP_3, 'v0')
                                        vLast_3      = self.m_GridData.gridInterpolateAveraging(posP_3, 'v1-')

                                    else:
                                        vNext_3      = self.m_GridData.gridInterpolateNodalBasis(posP_3, 'v1+')
                                        vPrev_3      = self.m_GridData.gridInterpolateNodalBasis(posP_3, 'v0')
                                        vLast_3      = self.m_GridData.gridInterpolateNodalBasis(posP_3, 'v1-')

                                    delv_delT    = np.zeros(3)
                                    delv_delT[0] = (0.5 * vNext_3[0] - 0.5 * vLast_3[0])/data_DT
                                    delv_delT[1] = (0.5 * vNext_3[1] - 0.5 * vLast_3[1])/data_DT
                                    delv_delT[2] = (0.5 * vNext_3[2] - 0.5 * vLast_3[2])/data_DT

                                v_3 = vPrev_3 + (simTime - tWin_0) * (vNext_3 - vPrev_3)/(tWin_1 - tWin_0)

                            elif t_3 > tWin_1:

                                if a_PolygonalCells == True:
                                    vNext_3      = self.m_GridData.gridInterpolateAveraging(posP_3, 'v2')
                                    vPrev_3      = self.m_GridData.gridInterpolateAveraging(posP_3, 'v1+')
                                    vLast_3      = self.m_GridData.gridInterpolateAveraging(posP_3, 'v0')

                                else:
                                    vNext_3      = self.m_GridData.gridInterpolateNodalBasis(posP_3, 'v2')
                                    vPrev_3      = self.m_GridData.gridInterpolateNodalBasis(posP_3, 'v1+')
                                    vLast_3      = self.m_GridData.gridInterpolateNodalBasis(posP_3, 'v0')

                                delv_delT    = np.zeros(3)
                                delv_delT[0] = (0.5 * vNext_3[0] - 0.5 * vLast_3[0])/data_DT
                                delv_delT[1] = (0.5 * vNext_3[1] - 0.5 * vLast_3[1])/data_DT
                                delv_delT[2] = (0.5 * vNext_3[2] - 0.5 * vLast_3[2])/data_DT

                                v_3 = vPrev_3 + (simTime - tWin_0) * (vNext_3 - vPrev_3)/(t_Next - tWin_1)

                            gradU1 = np.asarray(gradients.GetTuple9(cellPtIds.GetId(0)))
                            gradU2 = np.asarray(gradients.GetTuple9(cellPtIds.GetId(1)))
                            gradU3 = np.asarray(gradients.GetTuple9(cellPtIds.GetId(2)))
                            gradU4 = np.asarray(gradients.GetTuple9(cellPtIds.GetId(3)))

                            gradU = 0.25 * (gradU1   + gradU2   + gradU3   + gradU4)

                            lhsNS    = np.zeros(3)
                            lhsNS[0] = delv_delT[0] + v_3[0]*gradU[0] + v_3[1]*gradU[1] + v_3[2]*gradU[2]
                            lhsNS[1] = delv_delT[1] + v_3[0]*gradU[3] + v_3[1]*gradU[4] + v_3[2]*gradU[5]
                            lhsNS[2] = delv_delT[2] + v_3[0]*gradU[6] + v_3[1]*gradU[7] + v_3[2]*gradU[8]

                            vor      = np.zeros(3)
                            vor[0]   = gradU[7] - gradU[5]
                            vor[1]   = gradU[2] - gradU[6]
                            vor[2]   = gradU[3] - gradU[1]

                            Re_p = (self.m_InputData.getTracerDensity() * la.norm(v_3 - v_i - 0.5*el_2*dT) * 2 * self.m_InputData.m_TracerRadius)/self.m_InputData.getFluidViscosity()
                            Re_g = self.m_InputData.getFluidDensity() * la.norm(gradU) * ((2.0 * self.m_InputData.m_TracerRadius)**2.0)/self.m_InputData.getFluidViscosity()
                            alpha_LSA = 0.5 * (Re_g/(Re_p + 1.0e-16))
                            C_D = (24.0/(Re_p + 1.0e-16)) * (1 + c_1*(Re_p + 1.0e-16)**c_2) + c_3/(1 + (c_4/(Re_p + 1.0e-16)))
                            CoeffLift_extra = (3.0/(4.0 * np.pi * self.m_InputData.m_TracerRadius)) * np.sqrt((self.m_InputData.getFluidDensity() * self.m_InputData.getFluidViscosity()) / (la.norm(vor) + 1.0e-16))

                            if Re_p == 0:
                                coeffSL = 6.46 * CoeffLift_extra
                            elif Re_p <= 40.0:
                                coeffSL = 6.46 * ((1.0 - 0.3314 * np.sqrt(alpha_LSA)) * np.exp(-Re_p/10.0) + 0.3314 * np.sqrt(alpha_LSA)) * CoeffLift_extra
                            else:
                                coeffSL = 6.46 * 0.0524 * np.sqrt(0.5 * Re_g) * CoeffLift_extra

                            drag = d_1 * C_D * la.norm(v_3 - v_i  - 0.5*el_2*dT) * (v_3 - v_i  - 0.5*el_2*dT)
                            shearGradLift = coeffSL * np.cross(v_3 - v_i  - 0.5*el_2*dT, vor)
                            ambForce = self.m_InputData.getFluidDensity() * (1.0 + (C_am/2.0)) * lhsNS
                            bodyForce = (self.m_InputData.getTracerDensity() - self.m_InputData.getFluidDensity()) * accGrav

                            ## rk4 scheme step 3
                            k_3           = v_i + 0.5*el_2*dT
                            el_3          = (1.0/(self.m_InputData.getTracerDensity() + self.m_InputData.getFluidDensity()*(C_am/2.0))) * (drag + shearGradLift + ambForce + bodyForce)

                            posP_4        = posP_1 + k_3*dT
                            print('Rk4 Particle position 4 (no datasync)=', posP_4, simTime)
                            cell          = a_Locator.FindCell(list(posP_4))
                            cellPtIds     = vtk.vtkIdList()

                            if cell !=-1:
                                self.m_GridData.m_vtkData.GetCellPoints(cell, cellPtIds)

                                if t_4 <= tWin_1:

                                    if simTime > 0 and simTime < data_DT:

                                        if a_PolygonalCells == True:
                                            vNext_4  = self.m_GridData.gridInterpolateAveraging(posP_4, 'v1+')
                                            vPrev_4  = self.m_GridData.gridInterpolateAveraging(posP_4, 'v0')
                                        else:
                                            vNext_4  = self.m_GridData.gridInterpolateNodalBasis(posP_4, 'v1+')
                                            vPrev_4  = self.m_GridData.gridInterpolateNodalBasis(posP_4, 'v0')

                                        delv_delT    = np.zeros(3)
                                        delv_delT[0] = (vNext_4[0] - vPrev_4[0])/data_DT
                                        delv_delT[1] = (vNext_4[1] - vPrev_4[1])/data_DT
                                        delv_delT[2] = (vNext_4[2] - vPrev_4[2])/data_DT

                                    else:

                                        if a_PolygonalCells == True:
                                            vNext_4      = self.m_GridData.gridInterpolateAveraging(posP_4, 'v1+')
                                            vPrev_4      = self.m_GridData.gridInterpolateAveraging(posP_4, 'v0')
                                            vLast_4      = self.m_GridData.gridInterpolateAveraging(posP_4, 'v1-')

                                        else:
                                            vNext_4      = self.m_GridData.gridInterpolateNodalBasis(posP_4, 'v1+')
                                            vPrev_4      = self.m_GridData.gridInterpolateNodalBasis(posP_4, 'v0')
                                            vLast_4      = self.m_GridData.gridInterpolateNodalBasis(posP_4, 'v1-')

                                        delv_delT    = np.zeros(3)
                                        delv_delT[0] = (0.5 * vNext_4[0] - 0.5 * vLast_4[0])/data_DT
                                        delv_delT[1] = (0.5 * vNext_4[1] - 0.5 * vLast_4[1])/data_DT
                                        delv_delT[2] = (0.5 * vNext_4[2] - 0.5 * vLast_4[2])/data_DT

                                    v_4 = vPrev_4 + (simTime - tWin_0) * (vNext_4 - vPrev_4)/(tWin_1 - tWin_0)

                                elif t_4 > tWin_1:

                                    if a_PolygonalCells == True:
                                        vNext_4      = self.m_GridData.gridInterpolateAveraging(posP_4, 'v2')
                                        vPrev_4      = self.m_GridData.gridInterpolateAveraging(posP_4, 'v1+')
                                        vLast_4      = self.m_GridData.gridInterpolateAveraging(posP_4, 'v0')

                                    else:
                                        vNext_4      = self.m_GridData.gridInterpolateNodalBasis(posP_4, 'v2')
                                        vPrev_4      = self.m_GridData.gridInterpolateNodalBasis(posP_4, 'v1+')
                                        vLast_4      = self.m_GridData.gridInterpolateNodalBasis(posP_4, 'v0')

                                    delv_delT    = np.zeros(3)
                                    delv_delT[0] = (0.5 * vNext_4[0] - 0.5 * vLast_4[0])/data_DT
                                    delv_delT[1] = (0.5 * vNext_4[1] - 0.5 * vLast_4[1])/data_DT
                                    delv_delT[2] = (0.5 * vNext_4[2] - 0.5 * vLast_4[2])/data_DT

                                    v_4 = vPrev_4 + (simTime - tWin_0) * (vNext_4 - vPrev_4)/(t_Next - tWin_1)

                                gradU1 = np.asarray(gradients.GetTuple9(cellPtIds.GetId(0)))
                                gradU2 = np.asarray(gradients.GetTuple9(cellPtIds.GetId(1)))
                                gradU3 = np.asarray(gradients.GetTuple9(cellPtIds.GetId(2)))
                                gradU4 = np.asarray(gradients.GetTuple9(cellPtIds.GetId(3)))

                                gradU = 0.25 * (gradU1   + gradU2   + gradU3   + gradU4)

                                lhsNS    = np.zeros(3)
                                lhsNS[0] = delv_delT[0] + v_4[0]*gradU[0] + v_4[1]*gradU[1] + v_4[2]*gradU[2]
                                lhsNS[1] = delv_delT[1] + v_4[0]*gradU[3] + v_4[1]*gradU[4] + v_4[2]*gradU[5]
                                lhsNS[2] = delv_delT[2] + v_4[0]*gradU[6] + v_4[1]*gradU[7] + v_4[2]*gradU[8]

                                vor      = np.zeros(3)
                                vor[0]   = gradU[7] - gradU[5]
                                vor[1]   = gradU[2] - gradU[6]
                                vor[2]   = gradU[3] - gradU[1]

                                Re_p = (self.m_InputData.getTracerDensity() * la.norm(v_4 - v_i - 0.5*el_2*dT) * 2 * self.m_InputData.m_TracerRadius)/self.m_InputData.getFluidViscosity()
                                Re_g = self.m_InputData.getFluidDensity() * la.norm(gradU) * ((2.0 * self.m_InputData.m_TracerRadius)**2.0)/self.m_InputData.getFluidViscosity()
                                alpha_LSA = 0.5 * (Re_g/(Re_p + 1.0e-16))
                                C_D = (24.0/(Re_p + 1.0e-16)) * (1 + c_1*(Re_p + 1.0e-16)**c_2) + c_3/(1 + (c_4/(Re_p + 1.0e-16)))
                                CoeffLift_extra = (3.0/(4.0 * np.pi * self.m_InputData.m_TracerRadius)) * np.sqrt((self.m_InputData.getFluidDensity() * self.m_InputData.getFluidViscosity()) / (la.norm(vor) + 1.0e-16))

                                if Re_p == 0:
                                    coeffSL = 6.46 * CoeffLift_extra
                                elif Re_p <= 40.0:
                                    coeffSL = 6.46 * ((1.0 - 0.3314 * np.sqrt(alpha_LSA)) * np.exp(-Re_p/10.0) + 0.3314 * np.sqrt(alpha_LSA)) * CoeffLift_extra
                                else:
                                    coeffSL = 6.46 * 0.0524 * np.sqrt(0.5 * Re_g) * CoeffLift_extra

                                drag = d_1 * C_D * la.norm(v_4 - v_i  - 0.5*el_2*dT) * (v_4 - v_i  - 0.5*el_2*dT)
                                shearGradLift = coeffSL * np.cross(v_4 - v_i  - 0.5*el_2*dT, vor)
                                ambForce = self.m_InputData.getFluidDensity() * (1.0 + (C_am/2.0)) * lhsNS
                                bodyForce = (self.m_InputData.getTracerDensity() - self.m_InputData.getFluidDensity()) * accGrav

                                ## rk4 scheme step 4
                                k_4           = v_i + el_3*dT
                                el_4          = (1.0/(self.m_InputData.getTracerDensity() + self.m_InputData.getFluidDensity()*(C_am/2.0))) * (drag + shearGradLift + ambForce + bodyForce)

                                # new particle position and velocity after rk4 integration

                                posP_new        = posP_1 + (dT/6)*(k_1  + 2*k_2  + 2*k_3  + k_4)
                                v_new           = v_i   + (dT/6)*(el_1 + 2*el_2 + 2*el_3 + el_4)

                                gN1 = np.asarray(normVec.GetTuple3(cellPtIds.GetId(0)))
                                gN2 = np.asarray(normVec.GetTuple3(cellPtIds.GetId(1)))
                                gN3 = np.asarray(normVec.GetTuple3(cellPtIds.GetId(2)))
                                gN4 = np.asarray(normVec.GetTuple3(cellPtIds.GetId(3)))

                                if a_PolygonalCells == True:
                                    d   = self.m_GridData.gridInterpolateAveraging(posP_new,'sdf', a_GetGradient=False, a_GetStatus=False, a_DataDim=1, a_DataType='double')
                                else:
                                    d   = self.m_GridData.gridInterpolateNodalBasis(posP_new,'sdf', a_GetGradient=False, a_GetStatus=False, a_DataDim=1, a_DataType='double')
                                g   = 0.25*(gN1 + gN2 + gN3 + gN4)

                                if d <= self.m_InputData.m_TracerRadius:
                                    vn            = np.dot(v_new,g)/la.norm(g)

                                    if vn < 0.0:
                                        vt        = v_new - vn*g/la.norm(g)
                                        v_new     = vt - restitution*vn*g/la.norm(g)

                                posP_new = posP_new + v_new*dT
                                print('True Particle position (no datasync)=', posP_new)

                        #--------------------------------
                        # set new xyz for  point in loop
                        #--------------------------------
                            self.m_LagrangianData.setX(p, posP_new)
                            self.m_LagrangianData.setVectorData(v_new, a_ArrayName='Velocity', a_DataID=p)

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

    def getVelocityGradient(self, a_MeshFile, a_FieldName):

        mReader  = vtk.vtkXMLUnstructuredGridReader()
        mReader.SetFileName(a_MeshFile)
        mReader.Update()
        flowData = mReader.GetOutput()

        # compute gradient of v0
        gradientFilter = vtk.vtkGradientFilter()
        gradientFilter.ComputeDivergenceOff()
        gradientFilter.ComputeGradientOn()
        gradientFilter.ComputeQCriterionOff()
        gradientFilter.ComputeVorticityOff()
        gradientFilter.SetInputData(flowData)
        gradientFilter.SetResultArrayName('gradients')
        gradientFilter.SetInputArrayToProcess(0,0,0,vtk.vtkDataObject.FIELD_ASSOCIATION_POINTS, a_FieldName)
        gradientFilter.Update()
        velGradientData = gradientFilter.GetOutput()

        return velGradientData
