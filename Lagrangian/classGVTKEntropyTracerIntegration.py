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
import time
from decimal import Decimal

try:
    from classGVTKGrid import *
    from classGVTKLagrangianData import *
    from classGVTKGenericProblem import *
    from moduleFileIO import *
    from moduleVTKLocators import *
except ImportError:
    sys.exit('Could not import user defined modules')

class GVTKTracerIntegration(GVTKProblem):
    """A modular implementation of physics for passive tracer integration

    This is a derived class of the base class GVTKProblem. This follows standard
    definition of any instance of the base problem class - they must include a
    'runCompute' method, and a set of numerical integration methods.

    Attributes
    ----------
    m_InputData : SimInputs object
        (Base class attribute) Contains all input parameters for computation
    m_Debugode : bool
        (Base class attribute) Set True for verbose debug messages, False otherwise
    m_LagrangianData : GVTKLagrangianData
        Member lagrangian data container for the computation
    m_GridData : GVTKGridData
        Member vtk grid data container for the computation

    """

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

    def synchronizeDataUpdates(self, a_SimTime, a_TimeIndex, a_TimeDict, a_IsSingleStageIntegration):

        print("To Be Implemented")

    def runCompute(self):

        #----------------------------------------------------------
        # initialize tracer data from an external VTK polydata file
        # and update the Lagrangian tracer object for the problem
        #----------------------------------------------------------
        self.m_LagrangianDataArr = []
        self.Shannon_Entropy = self.m_InputData.getShannonEntropy()
        if self.Shannon_Entropy == "True":
            releases = self.m_InputData.getTracerInput()

            for i in range(0, len(releases)):
                self.m_LagrangianDataArr.append(GVTKLagrangianData(a_File = releases[i]))

        else:
            self.m_LagrangianDataArr.append(GVTKLagrangianData(a_File=self.m_InputData.getTracerInput()))


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

        #-------------------
        # start time counter
        #-------------------
        timeIndex   = 0
        simTime     = self.m_InputData.getSimulationStartTime()
        stopTime    = self.m_InputData.getSimulationStopTime()
        dT          = self.m_InputData.getIntegrationTimeStep()
        numSteps    = stopTime/dT
        tWin_0      = timeWindowDict['T_Low'][0]
        tWin_1      = timeWindowDict['T_Up'][0]


        self.m_isSteadyFlow = self.m_InputData.isSteadyFlowData()
        self.simPrec = self.m_InputData.getSimPrec()

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
                print("Integrating from", timeIndex, "to", timeIndex + 1, "| simTime")
            elif not self.m_isSteadyFlow and (timeWindowDict['T_Low'][timeIndex] != tWin_0 or timeWindowDict['T_Up'][timeIndex] != tWin_1):
                print("Integrating from", timeIndex, "to", timeIndex + 1, "| simTime")
            elif simTime > 0 and not self.m_ResumeSimulation:
                print("Integrating from", timeIndex, "to", timeIndex + 1, "| simTime", end='\r')

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

            #--------------------------------------------------------------------------------
            # STEP 3: resolve injections and inject particles
            # for simple tracer integration we are injecting particles either always or never
            # TODO: THIS NEEDS TO BE IMPROVED
            #--------------------------------------------------------------------------------
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

                    flowSingle      = self.m_InputData.getFlowDataFileName(a_Foam=isOpenFoam)
                    self.m_GridData.getDataFromFile(flowSingle)
                    self.m_GridData.renameDataArray(self.m_InputData.getVelDataName(), 'v0')
                    # self.m_GridData.addDataArrayFromFile(boundaryFileName, 'dist', 'sdf')

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

                    if simTime == self.m_InputData.getSimulationStartTime():

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

                        if not self.m_InputData.isFixedMesh():
                            self.m_ResumeSimulation = False

                    elif isDataPointSync:

                        if isSingleStageIntegration == True:

                            self.m_GridData.removeDataArray('v0')
                            flowSingle      = self.m_InputData.getFlowDataFileName(a_ID = timeWindowDict['ID_Low'][timeIndex])
                            self.m_GridData.addDataArrayFromFile(flowSingle, 'v0', self.m_InputData.getVelDataName())

                        else:

                            self.m_GridData.removeDataArray('v0')
                            self.m_GridData.removeDataArray('v1')
                            flowPlus    = self.m_InputData.getFlowDataFileName(a_ID = timeWindowDict['ID_Up'][timeIndex])
                            flowMinus   = self.m_InputData.getFlowDataFileName(a_ID = timeWindowDict['ID_Low'][timeIndex])
                            self.m_GridData.addDataArrayFromFile(flowPlus, 'v1', self.m_InputData.getVelDataName())
                            self.m_GridData.addDataArrayFromFile(flowMinus, 'v0', self.m_InputData.getVelDataName())
                            print('Loaded a new velocity data')

                    else:

                        if isSingleStageIntegration == True:

                            self.m_GridData.removeDataArray('v0')
                            self.m_GridData.removeDataArray('v1')
                            flowPlus    = self.m_InputData.getFlowDataFileName(a_ID = timeWindowDict['ID_Up'][timeIndex])
                            flowMinus   = self.m_InputData.getFlowDataFileName(a_ID = timeWindowDict['ID_Low'][timeIndex])
                            self.m_GridData.addDataArrayFromFile(flowPlus, 'v1', self.m_InputData.getVelDataName())
                            self.m_GridData.addDataArrayFromFile(flowMinus, 'v0', self.m_InputData.getVelDataName())
                            print('Loaded a new velocity data')

                        else:

                            if simTime < tWin_1 - dT:

                                self.m_GridData.removeDataArray('v0')
                                self.m_GridData.removeDataArray('v1')
                                self.m_GridData.removeDataArray('v2')
                                flowPlus        = self.m_InputData.getFlowDataFileName(a_ID = timeWindowDict['ID_Up'][timeIndex])
                                flowMinus       = self.m_InputData.getFlowDataFileName(a_ID = timeWindowDict['ID_Low'][timeIndex])
                                self.m_GridData.addDataArrayFromFile(flowPlus, 'v1', self.m_InputData.getVelDataName())
                                self.m_GridData.addDataArrayFromFile(flowMinus, 'v0', self.m_InputData.getVelDataName())
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
                                print('Loaded a new velocity data')

            #----------------------------------------------------------------------
            # STEP 5: build a cell locator (if it is a fixed mesh, built only once)
            #----------------------------------------------------------------------
            if (simTime == self.m_InputData.getSimulationStartTime() and self.m_InputData.isFixedMesh()) or (self.m_ResumeSimulation and self.m_InputData.isFixedMesh()):
                print("Building Cell Locator Maps")
                self.m_GridData.buildLocator()

            for i in range(0, len(self.m_LagrangianDataArr)):
                self.m_LagrangianData = self.m_LagrangianDataArr[i]
            #-----------------------------------------------
            # STEP 6: Output initial positions of particles
            #-----------------------------------------------
                if simTime == self.m_InputData.getSimulationStartTime() or self.m_ResumeSimulation:
                    self.m_LagrangianData.writeData(self.m_InputData.getTracerOutputFile(a_ID1=timeIndex)[i])
                    self.m_ResumeSimulation = False

                if self.m_ResumeSimulation:
                    for i in range(reader.GetOutput().GetNumberOfPoints()):
                        point = reader.GetOutput().GetPoint(i)
                        self.m_LagrangianData.setX(i, point)
        #----------------------------------------------------------
                # STEP 7: now proceed with integration of each tracer point
                #----------------------------------------------------------
                boundaryCondition = 1   # currently a dummy variable

                if self.m_InputData.getIntegrationScheme() == 'feu':
                    self.integrateForwardEuler(simTime,tWin_0,tWin_1,boundaryCondition,a_DataSync=isDataPointSync,a_PolygonalCells=True)
                elif self.m_InputData.getIntegrationScheme() == 'rk4':
                    self.integrateRK4(simTime,tWin_0,tWin_1,tWin_2,boundaryCondition,a_DataSync=isDataPointSync,a_PolygonalCells=True)

                #------------------------------------------------
                # STEP 8: update time indices and simulation time
                #------------------------------------------------
                timeIndex += 1
                simTime    = round(simTime+dT, self.simPrec)

                #----------------------------------------------------------------------
                # STEP 9: at the end of appropriate number of steps dump data into file
                #----------------------------------------------------------------------
                if timeIndex%self.m_InputData.m_WriteInterval == 0:
                    self.m_LagrangianData.writeData(self.m_InputData.getTracerOutputFile(a_ID1=timeIndex)[i])



    def integrateForwardEuler(self, a_T, a_T0, a_T1, a_BoundaryCondition, a_DataSync=True, a_PolygonalCells=False):
        """Forward Euler time integration routine for individual particles

        This implementation is designed in such a way that the integrator becomes
        an instanced method of the physics class itself (see VTK based OOP design notes).

        Notes
        -----
        For future attempts at using shared memory parallelism features for this
        library, the integrator routine is where that parallelism becomes key
        since all particles are integrated and updated without reference to each other.

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

        """

        for p in range(self.numParticles):

            xyz = self.m_LagrangianData.getX(p)

            if a_DataSync == True:
                if a_PolygonalCells == True:
                    vel = self.m_GridData.gridInterpolateAveraging(xyz, 'v0')
                else:
                    vel = self.m_GridData.gridInterpolateNodalBasis(xyz, 'v0')
            else:
                if a_PolygonalCells == True:
                    velPlus = self.m_GridData.gridInterpolateAveraging(xyz, 'v1')
                    velMinus = self.m_GridData.gridInterpolateAveraging(xyz, 'v0')
                else:
                    velPlus = self.m_GridData.gridInterpolateNodalBasis(xyz, 'v1')
                    velMinus = self.m_GridData.gridInterpolateNodalBasis(xyz, 'v0')

                vel = velMinus + (a_T - a_T0) * (velPlus - velMinus)/(a_T1 - a_T0)

            xyz = xyz + vel * self.m_InputData.getIntegrationTimeStep()

            #
            # NOTE: AD-HOC REPLACE WITH BOUNDARY CONDITION FUNCTION
            # NOTE: NEED TO IMPROVE BOUNDARY CONDITION IMPLEMENTATION
            #
            # if a_BoundaryCondition:
            #     print("To Be Implemented")

            self.m_LagrangianData.setX(p, xyz)

    #--------------------------------------------------------------------------
    # time integration routine for individual particles: Runge Kutta Integrator
    #
    # implementation of this integrator follows the standard operation loop:
    # locate -> interpolate in space -> interpolate in time -> integrate
    #
    # k1 = v(x0, t0)
    # k2 = v(x0 + 0.5*dt*k1, t0+0.5*dt)
    # k3 = v(x0 + 0.5*dt*k2, t0+0.5*dt)
    # k4 = v(x0 + dt*k3,     t0 + dt)
    # x  = x0 + (dt/6)*(k1 + 2*k2 + 2*k3 + k4)
    # t* = t0 + 0.5*dt and t** = t0+dt should both be in the time window loaded
    #--------------------------------------------------------------------------
    def integrateRK4(self, a_T, a_T0, a_T1, a_TNext, a_BoundaryCondition, a_Steady=True, a_DataSync=True, a_PolygonalCells=False):

        dT = self.m_InputData.getIntegrationTimeStep()

        for p in range(self.numParticles):

            xyz = self.m_LagrangianData.getX(p)

            t_1 = a_T
            t_2 = a_T + 0.50*dT
            t_3 = a_T + 0.50*dT
            t_4 = a_T + dT

            if a_DataSync == True:

                if a_PolygonalCells == True:

                    k_1     = self.m_GridData.gridInterpolateAveraging(xyz,'v0')
                    xyz_2   = xyz + 0.50*dT*k_1
                    k_2     = self.m_GridData.gridInterpolateAveraging(xyz_2,'v0')
                    xyz_3   = xyz + 0.50*dT*k_2
                    k_3     = self.m_GridData.gridInterpolateAveraging(xyz_3,'v0')
                    xyz_4   = xyz + dT*k_3
                    k_4     = self.m_GridData.gridInterpolateAveraging(xyz_4,'v0')

                else:

                    k_1     = self.m_GridData.gridInterpolateNodalBasis(xyz,'v0')
                    xyz_2   = xyz + 0.50*dT*k_1
                    k_2     = self.m_GridData.gridInterpolateNodalBasis(xyz_2,'v0')
                    xyz_3   = xyz + 0.50*dT*k_2
                    k_3     = self.m_GridData.gridInterpolateNodalBasis(xyz_3,'v0')
                    xyz_4   = xyz + dT*k_3
                    k_4     = self.m_GridData.gridInterpolateNodalBasis(xyz_4,'v0')

            else:

                if a_PolygonalCells == True:
                    vel_1_L = self.m_GridData.gridInterpolateAveraging(xyz, 'v0')
                    vel_1_U = self.m_GridData.gridInterpolateAveraging(xyz, 'v1')
                else:
                    vel_1_L = self.m_GridData.gridInterpolateNodalBasis(xyz, 'v0')
                    vel_1_U = self.m_GridData.gridInterpolateNodalBasis(xyz, 'v1')

                k_1     = vel_1_L + (a_T - a_T0)*(vel_1_U - vel_1_L)/(a_T1 - a_T0)

                if t_2 <= a_T1:

                    xyz_2   = xyz + 0.50*dT*k_1

                    if a_PolygonalCells == True:
                        vel_2_L = self.m_GridData.gridInterpolateAveraging(xyz_2, 'v0')
                        vel_2_U = self.m_GridData.gridInterpolateAveraging(xyz_2, 'v1')
                    else:
                        vel_2_L = self.m_GridData.gridInterpolateNodalBasis(xyz_2, 'v0')
                        vel_2_U = self.m_GridData.gridInterpolateNodalBasis(xyz_2, 'v1')

                    k_2     = vel_2_L + (a_T - a_T0)*(vel_2_U - vel_2_L)/(a_T1 - a_T0)

                elif t_2 > a_T1:

                    xyz_2   = xyz + 0.50*dT*k_1

                    if a_PolygonalCells == True:
                        vel_2_L = self.m_GridData.gridInterpolateAveraging(xyz_2, 'v1')
                        vel_2_U = self.m_GridData.gridInterpolateAveraging(xyz_2, 'v2')
                    else:
                        vel_2_L = self.m_GridData.gridInterpolateNodalBasis(xyz_2, 'v1')
                        vel_2_U = self.m_GridData.gridInterpolateNodalBasis(xyz_2, 'v2')

                    k_2     = vel_2_L + (a_T - a_T1)*(vel_2_U - vel_2_L)/(a_TNext - a_T1)

                if t_3 <= a_T1:

                    xyz_3   = xyz + 0.50*dT*k_2

                    if a_PolygonalCells == True:
                        vel_3_L = self.m_GridData.gridInterpolateAveraging(xyz_3, 'v0')
                        vel_3_U = self.m_GridData.gridInterpolateAveraging(xyz_3, 'v1')
                    else:
                        vel_3_L = self.m_GridData.gridInterpolateNodalBasis(xyz_3, 'v0')
                        vel_3_U = self.m_GridData.gridInterpolateNodalBasis(xyz_3, 'v1')

                    k_3     = vel_3_L + (a_T - a_T0)*(vel_3_U - vel_3_L)/(a_T1 - a_T0)

                elif t_3 > a_T1:

                    xyz_3   = xyz + 0.50*dT*k_2

                    if a_PolygonalCells == True:
                        vel_3_L = self.m_GridData.gridInterpolateAveraging(xyz_3, 'v1')
                        vel_3_U = self.m_GridData.gridInterpolateAveraging(xyz_3, 'v2')
                    else:
                        vel_3_L = self.m_GridData.gridInterpolateNodalBasis(xyz_3, 'v1')
                        vel_3_U = self.m_GridData.gridInterpolateNodalBasis(xyz_3, 'v2')

                    k_3     = vel_3_L + (a_T - a_T1)*(vel_3_U - vel_3_L)/(a_TNext - a_T1)

                if t_4 <= a_T1:

                    xyz_4   = xyz + dT*k_3

                    if a_PolygonalCells == True:
                        vel_4_L = self.m_GridData.gridInterpolateAveraging(xyz_4, 'v0')
                        vel_4_U = self.m_GridData.gridInterpolateAveraging(xyz_4, 'v1')
                    else:
                        vel_4_L = self.m_GridData.gridInterpolateNodalBasis(xyz_4, 'v0')
                        vel_4_U = self.m_GridData.gridInterpolateNodalBasis(xyz_4, 'v1')

                    k_4     = vel_4_L + (a_T - a_T0)*(vel_4_U - vel_4_L)/(a_T1 - a_T0)

                elif t_4 > a_T1:

                    xyz_4   = xyz + dT*k_3

                    if a_PolygonalCells == True:
                        vel_4_L = self.m_GridData.gridInterpolateAveraging(xyz_4, 'v1')
                        vel_4_U = self.m_GridData.gridInterpolateAveraging(xyz_4, 'v2')
                    else:
                        vel_4_L = self.m_GridData.gridInterpolateNodalBasis(xyz_4, 'v1')
                        vel_4_U = self.m_GridData.gridInterpolateNodalBasis(xyz_4, 'v2')

                    k_4     = vel_4_L + (a_T - a_T1)*(vel_4_U - vel_4_L)/(a_TNext - a_T1)

            xyz = xyz + (dT/6.0)*(k_1 + 2.0*k_2 + 2.0*k_3 + k_4)

            #
            # NOTE: AD-HOC REPLACE WITH BOUNDARY CONDITION FUNCTION
            # NOTE: NEED TO IMPROVE BOUNDARY CONDITION IMPLEMENTATION
            #
            # if a_BoundaryCondition:
            #     print("To Be Implemented")

            self.m_LagrangianData.setX(p, xyz)
