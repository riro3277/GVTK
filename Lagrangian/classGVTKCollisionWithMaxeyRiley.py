#-------------------------------------------------------------------------------
# This module computes the embolus dynamics using a one way coupling scheme based
# on signed distance field algorithm for collision and Maxey Riley Equations for
# fluid forces.
#
# Author:      Akshita Sahni, Sreeparna Majee, Ricardo Roopnarinsingh
# Institution: University of Colorado, Boulder
# Last Edit:   Nov. 2022
#-------------------------------------------------------------------------------
import sys, os, vtk
import numpy as np
import time

from numpy import linalg as la
from numba import jit
try:
    from classGVTKGrid import *
    from classGVTKLagrangianData import *
    from classGVTKGenericProblem import *
    from testModule import *
    #from moduleVTKLocators import *
except:
    sys.exit('Could not import user defined modules')

class GVTKCollision(GVTKProblem):

    def __init__(self, a_InputData, a_DebugMode=True):

        super().__init__(a_InputData,a_DebugMode)
        self.m_LagrangianData = None
        self.m_GridData       = None
        self.Mode = sys.argv[1]

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
        #timeWindowDict = self.m_InputData.getDataTimeWindows()

        #-----------------------------------------------------------------------------
        # display the time windows in case the integration involves unsteady flow data
        #-----------------------------------------------------------------------------
        if not(self.m_InputData.isSteadyFlowData()):
            self.m_InputData.printDataTimeWindows()

        #-----------------------------------------------------------------------
        # start time counter
        #-----------------------------------------------------------------------

        #Clears and creates Test Files for Regression Testings
        if self.Mode == 'test':
            self.Root_dir = self.m_InputData.getRootPath()
            testFilesDir = self.Root_dir + "/TestFiles/"
            for p in range(self.numParticles):
                if os.path.isfile(testFilesDir + 'TestFile' + str(p) + '.txt') == True:
                    with open(testFilesDir + 'TestFile' + str(p) + '.txt', 'w'):
                        pass

        timeIndex    = 0
        timeElapsed0 = 0.0
        timeElapsed1 = 0.0
        simPrec = self.m_InputData.getSimPrec()
        simTime   = round(self.m_InputData.getSimulationStartTime(), simPrec+1)
        stopTime  = self.m_InputData.getSimulationStopTime()
        dT        = self.m_InputData.getIntegrationTimeStep()
        tWin_0    = timeWindowDict['T_Low'][0]
        tWin_1    = timeWindowDict['T_Up'][0]
        data_DT   = self.m_InputData.getFlowDataTimeStep()
        isSingleStageIntegration = True # set true for feu and false for for rk4

        #-----------------------------------------------------------------------
        # get signed distance field
        #-----------------------------------------------------------------------
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
            print(dest_folder, prev_files)

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

                    #---------for new simulation from start--------#
                    if simTime == self.m_InputData.getSimulationStartTime():


                        flowInitial = self.m_InputData.getFlowDataFileName(a_ID = timeWindowDict['ID_Low'][timeIndex])
                        flowNext    = self.m_InputData.getFlowDataFileName(a_ID = timeWindowDict['ID_Low'][timeIndex] + 1)
                        self.m_GridData.addDataArrayFromFile(flowInitial, 'v0', self.m_InputData.getVelDataName())
                        self.m_GridData.addDataArrayFromFile(flowNext, 'v1+', self.m_InputData.getVelDataName())

                        gradVel = self.getVelocityGradient(flowInitial, self.m_InputData.getVelDataName())


                    #---------for resumed simulation----------#
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

                    #-------for lagrangian data is in sync with flow data------#
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

                    #-------for lagrangian data not in sync with flow data------#
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
                    self.m_GridData.buildLocator(a_LocatorType = self.m_InputData.getLocatorType())
                    locatorObj = self.m_GridData.m_Locator
                else:
                    self.m_GridData.buildLocator(a_LocatorType = self.m_InputData.getLocatorType())
                    locatorObj = self.m_GridData.m_Locator
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
            if self.m_InputData.m_InitialVelocityStatus == "True":
                if simTime == self.m_InputData.getSimulationStartTime() or self.m_ResumeSimulation:
                    self.m_LagrangianData.addVectorData('Velocity', InitialVelocity)
                    self.m_ResumeSimulation = False
            #----------------------------------------------------------------
            # STEP 7: now proceed with integration of each inertial particle
            #----------------------------------------------------------------
            boundaryCondition = 1   # currently a dummy variable

            # measure time for particle collision
            t2 = time.process_time()
            self.sdfCollision(simTime, tWin_0, tWin_1, gradVel, boundaryCondition, locatorObj, a_DataSync=isDataPointSync, a_PolygonalCells=False)
            t3 = time.process_time()

            #----------------------------------------------------------------------
            # STEP 8: at the end of appropriate number of steps dump data into file
            #----------------------------------------------------------------------
            if timeIndex%self.m_InputData.m_WriteInterval == 0:
                self.m_LagrangianData.writeData(self.m_InputData.getTracerOutputFile(a_ID1=0, a_ID2=timeIndex))
            # if self.Mode == "test":
            #      baseDir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
            #      sys.exit(print("Here", baseDir, self.Root_dir))
            #      case = Testing(baseDir, self.Root_dir)
            #      case.Compare()
            #      case.Documentation()


            #------------------------------------------------
            # STEP 9: update time indices and simulation time
            #------------------------------------------------
            timeIndex = timeIndex + 1
            simTime   = round(simTime + dT, simPrec+1)
            timeElapsed1 += t3 - t2
            print('elapsed time in particle tracking/collision till last loop: ', timeElapsed1)
        print('Elapsed time in particle tracking/collision: ', timeElapsed1)
        #Runs regression test module for doumentation and test to standard file comparison
        if self.Mode == "test":
             baseDir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
             case = Testing(self.Root_dir, baseDir)
             case.Compare()
             case.Documentation()
            #-------------------------

    #------collision function starts here--------#
    # @jit(nopython=True)
    def sdfCollision(self, simTime, tWin_0, tWin_1, gradVel, boundaryCondition, a_Locator, a_DataSync=True, a_PolygonalCells=False):

        dist            = self.m_GridData.extractDataArray('sdf')
        normVec         = self.m_GridData.extractDataArray('sdf-grad')
        gradients       = gradVel.GetPointData().GetArray('gradients')
        # print('gradients=', gradients)
        # sys.exit()
        dT              = self.m_InputData.getIntegrationTimeStep() #--------getting integration delta T value from input-------#
        data_DT         = self.m_InputData.getFlowDataTimeStep() #--------getting flow data delta T value from input-------#
        restitution     = self.m_InputData.getRestitution() #--------getting restitution value from input-------#

        #-------------------------------------------------------------------------------#
        #----------Calculating constants for each term in Maxey Riley Equation----------#
        #-------------------------------------------------------------------------------#
        c_1 = 0.1806
        c_2 = 0.6459
        c_3 = 0.4251
        c_4 = 6886.95
        d_1 = (3.0 * self.m_InputData.getFluidDensity())/(8.0 * self.m_InputData.m_TracerRadius)
        C_am    = 1.0
        accGrav = self.m_InputData.getGravityVector()
        # print('gravity=', accGrav)
        # sys.exit()

        #------------------------------------------
        # forward euler integration implementation
        #------------------------------------------
        TDrag = 0.0
        Tshear = 0.0
        Tadded = 0.0
        Tuforce = 0.0
        Tbody = 0.0
        c = 0.0

        for p in range(self.numParticles):
            posP         = self.m_LagrangianData.getX(p) # particle position

            if self.m_InputData.m_InitialVelocityStatus == "False":
                if simTime == self.m_InputData.getSimulationStartTime():
                        if a_PolygonalCells == True:
                            v           = self.m_GridData.gridInterpolateAveraging(posP, 'v0')
                        else:
                            v           = self.m_GridData.gridInterpolateNodalBasis(posP, 'v0')
                        InitialVelocity = np.zeros((self.numParticles, 3))
                        InitialVelocity[:,0] = float(v[0])
                        InitialVelocity[:,1] = float(v[1])
                        InitialVelocity[:,2] = float(v[2])

                        self.m_LagrangianData.addVectorData('Velocity', InitialVelocity)
            v_i          = np.array(self.m_LagrangianData.getVectorData(a_ArrayName='Velocity', a_DataID=p)) # particle velocity
            # print('velocityp=', v_i)
            # sys.exit()
            ## error with using the same name for cell assisgnment: FindCell argument 1: 'tuple' object does not support item assignment
            cell         = a_Locator.FindCell(list(posP))
            cellPtIds    = vtk.vtkIdList()

            delv_delT    = np.zeros(3)

            if simTime == self.m_InputData.getSimulationStartTime():
                if a_PolygonalCells == True:
                    v           = self.m_GridData.gridInterpolateAveraging(posP, 'v0')
                    vNext       = self.m_GridData.gridInterpolateAveraging(posP, 'v1+')

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

                #--------Edits for Gradient calculation based on Nodal Basis interpolation is done by Sreeparna----------#
                gradUU = self.m_GridData.gridInterpolateNodalBasis(posP,'v0', a_GetGradient=True, a_GetStatus=False, a_DataDim=3, a_DataType='double')
                gradU = np.array(gradUU.flat) #----Matrix is converted to arrays edited by Sreeparna----#

                #---------Calculate the material derivative of fluid velocity field-----------#
                lhsNS    = np.zeros(3)
                lhsNS[0] = delv_delT[0] + v[0]*gradU[0] + v[1]*gradU[1] + v[2]*gradU[2]
                lhsNS[1] = delv_delT[1] + v[0]*gradU[3] + v[1]*gradU[4] + v[2]*gradU[5]
                lhsNS[2] = delv_delT[2] + v[0]*gradU[6] + v[1]*gradU[7] + v[2]*gradU[8]

                #-------Calculate vorticity from the components of the velocity gradient tensor-------#
                vor      = np.zeros(3)
                vor[0]   = gradU[7] - gradU[5]
                vor[1]   = gradU[2] - gradU[6]
                vor[2]   = gradU[3] - gradU[1]

                vol_p = 4.0*3.14*(self.m_InputData.m_TracerRadius)**3.0/3.0
                #-------calculation of particle slip velocity Reynolds number--------#
                # Re_p = (self.m_InputData.getTracerDensity() * la.norm(v - v_i) * 2 * self.m_InputData.m_TracerRadius)/self.m_InputData.getFluidViscosity()
                Re_p = (self.m_InputData.getFluidDensity() * la.norm(v - v_i) * 2.0 * self.m_InputData.m_TracerRadius)/self.m_InputData.getFluidViscosity()

                #-------calculation of shear based Reynolds number--------#
                Re_g = self.m_InputData.getFluidDensity() * la.norm(gradU) * ((2.0 * self.m_InputData.m_TracerRadius)**2.0)/self.m_InputData.getFluidViscosity()

                alpha_LSA = 0.5 * (Re_g/(Re_p + 1.0e-20))
                # CoeffLift_extra = (3.0/(4.0 * np.pi * self.m_InputData.m_TracerRadius)) * np.sqrt((self.m_InputData.getFluidDensity() * self.m_InputData.getFluidViscosity()) / (la.norm(vor) + 1.0e-16))
                CoeffLift_extra = self.m_InputData.getFluidDensity()*(self.m_InputData.m_TracerRadius)**2.0/(vol_p)*np.sqrt((self.m_InputData.getFluidViscosity()/self.m_InputData.getFluidDensity())/(la.norm(vor) + 1.0e-20))

                if Re_p == 0:
                    coeffSL = 6.46 * CoeffLift_extra
                elif Re_p <= 40.0:
                    coeffSL = 6.46 * ((1.0 - 0.3314 * np.sqrt(alpha_LSA)) * np.exp(-Re_p/10.0) + 0.3314 * np.sqrt(alpha_LSA)) * CoeffLift_extra
                else:
                    # coeffSL = 6.46 * 0.0524 * np.sqrt(0.5 * Re_g) * CoeffLift_extra
                    coeffSL = 6.46 * 0.0524 * np.sqrt(alpha_LSA * Re_p) * CoeffLift_extra

                #-----------calculation of Drag coefficient----------#
                C_D = (24.0/(Re_p + 1.0e-20)) * (1 + c_1*(Re_p + 1.0e-20)**c_2) + c_3/(1 + (c_4/(Re_p + 1.0e-20)))

                #-----------calculation of Drag Force----------#
                drag = d_1 * C_D * la.norm(v - v_i) * (v - v_i)

                #-----------calculation of Shear Lift Force----------#
                shearGradLift = coeffSL * np.cross(v - v_i, vor)

                #--------------Calculation Of Added Mass Force----------#
                ambForce = self.m_InputData.getFluidDensity() * (1.0 + (C_am/2.0)) * lhsNS
                addedmass = self.m_InputData.getFluidDensity() * 0.5 * C_am * lhsNS
                undisturbedforce = self.m_InputData.getFluidDensity() * lhsNS

                #---------Calculation Of Buoyancy-----------#
                bodyForce = (self.m_InputData.getTracerDensity() - self.m_InputData.getFluidDensity()) * accGrav

                #--------Now sum up all the forces to calculate the total force on the particle----------#
                v_i = v_i + (1.0/(self.m_InputData.getTracerDensity() + self.m_InputData.getFluidDensity()*(C_am/2.0))) * dT * (drag + shearGradLift + ambForce + bodyForce)

                #------------Calculation of SDF Gradients----------#
                gN1 = np.asarray(normVec.GetTuple3(cellPtIds.GetId(0)))
                gN2 = np.asarray(normVec.GetTuple3(cellPtIds.GetId(1)))
                gN3 = np.asarray(normVec.GetTuple3(cellPtIds.GetId(2)))
                gN4 = np.asarray(normVec.GetTuple3(cellPtIds.GetId(3)))

                #---------calculate distance of particle position from wall-----------#
                d   = self.m_GridData.gridInterpolateNodalBasis(posP,'sdf', a_GetGradient=False, a_GetStatus=False, a_DataDim=1, a_DataType='double')
                g   = 0.25*(gN1 + gN2 + gN3 + gN4)

                #-----------checking for collision-----------#
                if d <= self.m_InputData.m_TracerRadius:
                #----calculate normal component of particle velocity----#
                    vn      = np.dot(v_i,g)/la.norm(g)
                #-------if the direction of the normal component of velocity is inside the domain------#
                    if vn < 0.0:
                #----calculate tangential component of particle velocity----#
                        vt      = v_i - vn*g/la.norm(g)
                #----calculate updated particle velocity after collision----#
                        v_i     = vt - restitution*vn*g/la.norm(g)


                #-------calculate updated particle position due to change in velocity-------#
                posP = posP + v_i*dT

        #--------------------------------
        # set new xyz for  point in loop
        #--------------------------------

        #--------------------------------
        # Regression Test File Creation
        #--------------------------------
            if self.Mode == 'test':

                #Root_dir is the cylindrical_flow folder
                Root_dir = self.Root_dir + "/TestFiles/"

                Standard_File = open(Root_dir+'TestFile' +  str(p) + '.txt', 'a')


                #Writes pos, vel, and particle forces to test files in cylinder_flow folder
                Standard_File.write('t:' + str(simTime) + ' Position:')
                for i in range(len(posP)):
                    Standard_File.write(str(posP[i]) + ",")
                Standard_File.write('\n')

                Standard_File.write('t:' + str(simTime) + ' Velocity:')
                for i in range(len(v_i)):
                    Standard_File.write(str(v_i[i]) + ",")
                Standard_File.write('\n')

                Standard_File.write('t:' + str(simTime) + ' Drag:')
                for i in range(len(drag)):
                    Standard_File.write(str(drag[i]) + ",")
                Standard_File.write('\n')

                Standard_File.write('t:' + str(simTime) + ' ShearGradLift:')
                for i in range(len(shearGradLift)):
                    Standard_File.write(str(shearGradLift[i]) + ",")
                Standard_File.write('\n')

                Standard_File.write('t:' + str(simTime) + ' BodyForce:')
                for i in range(len(bodyForce)):
                    Standard_File.write(str(bodyForce[i]) + ",")
                Standard_File.write('\n')

            self.m_LagrangianData.setX(p, posP)
            self.m_LagrangianData.setVectorData(v_i, a_ArrayName='Velocity', a_DataID=p)


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
