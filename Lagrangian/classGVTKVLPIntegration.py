#-------------------------------------------------------------------------------
# This is a module that encapsulates the physics and numerics of an inertial
# model integration calculation through a specified flow field
#
# Author:      Joseph Wilson
# Institution: University of Colorado, Boulder
# Last Edit:   Oct 2021
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

class GVTKViral(GVTKProblem):

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
            tauV = 2.0*self.m_InputData.m_TracerDensity*self.m_InputData.m_AverageRadius**2/(9.0*self.m_InputData.m_FluidViscosity)
        compare_dT = dT/tauV
        if compare_dT > 0.5:
            recommended_dT = 0.5*tauV
            warn_message = 'Integration time step is high relative to the momentum response time which could cause unstable results. Recommend resetting time step to no greater than ' + str(round(recommended_dT, self.m_dt_prec+2)) + 's.'
            warnings.warn(warn_message)
            print()

    def calculateAverageVelocity(self):
        problemDimension  = self.m_InputData.getProblemDimension()
        numFlowDataPoints = self.m_GridData.m_vtkData.GetNumberOfPoints()
        dataArray_0       = self.m_GridData.m_vtkData.GetPointData().GetArray('v0')
        v                 = 0.0

        for p in range(numFlowDataPoints):
            if problemDimension == 2:
                vel = dataArray_0.GetTuple2(p)
                v += np.sqrt(vel[0]**2 + vel[1]**2)
            else:
                vel = dataArray_0.GetTuple3(p)
                v += np.sqrt(vel[0]**2 + vel[1]**2 + vel[2]**2)

        self.m_AverageFlowVelocity = v/numFlowDataPoints

    def calculateAverageTemperature(self):
        numFlowDataPoints = self.m_GridData.m_vtkData.GetNumberOfPoints()
        dataArray         = self.m_GridData.m_vtkData.GetPointData().GetArray('T0')
        T                 = 0.0

        for p in range(numFlowDataPoints):
                T += dataArray.GetTuple1(p)

        self.m_AverageTemperature = T/numFlowDataPoints
        self.m_AverageSaturationPressure = self.calculateSaturationPressure(self.m_AverageTemperature, self.m_InputData.m_TempUnits, '(domain)')

    def calculateAveragePressure(self):
        numFlowDataPoints = self.m_GridData.m_vtkData.GetNumberOfPoints()
        dataArray         = self.m_GridData.m_vtkData.GetPointData().GetArray('T0')
        press             = 0.0

        for p in range(numFlowDataPoints):
                press += dataArray.GetTuple1(p)

        p_avg = (press/numFlowDataPoints) * self.m_PressureConversionFactor + self.m_AbsPressureConversion

        self.m_AveragePressure = p_avg

    def calculateMomentumResponseTime(self, a_R):
        d_p   = 2.0*a_R
        rho_p = self.m_InputData.getTracerDensity()
        mu_f  = self.m_InputData.getFluidViscosity()
        return rho_p*(d_p**2)/(18.0*mu_f)

    def calculateThermalResponseTime(self, a_R):
        d_p   = a_R*2.0
        rho_p = self.m_InputData.m_TracerDensity
        c_p   = self.m_InputData.m_SpecificHeat
        k     = self.m_FluidThermalConductivity
        return rho_p*c_p*(d_p**2)/(12.0*k)

    def calculateEvaporationFactor(self, a_R):
        d_p   = 2.0*a_R
        Sh    = 2.0
        rho_f = self.m_InputData.m_FluidDensity
        rho_p = self.m_InputData.m_TracerDensity
        D     = self.m_InputData.m_TracerDiffCoeff
        return 2.0*Sh*rho_f*D/(d_p*rho_p)

    def calculateSaturationPressure(self, a_T, a_T_Units, n):
        if a_T_Units == 'C':
            a_T += 273.1
        elif a_T_Units == 'F':
            a_T  = ((a_T-32)/1.8) + 273.1
        elif a_T_Units == 'R':
            a_T  = ((a_T-491.67)/1.8) + 273.1

        '''
        Below data pulled from:
        https://webbook.nist.gov/cgi/cbook.cgi?ID=C7732185&Mask=4&Type=ANTOINE&Plot=on
        '''
        if a_T<255.9:
            errorMessage = 'The temperature of particle '+str(n)+' is '+str(a_T)+a_T_Units+'. Cannot find partial pressures when particle temperature is <255.9 K'
            raise ValueError(errorMessage)
        elif a_T>=255.9 and a_T<274:
            A, B, C = 4.6543, 1432.264, -64.848
        elif a_T>=274 and a_T<304:
            A, B, C = 5.40221, 1838.675, -31.737
        elif a_T>=304 and a_T<333:
            A, B, C = 5.20389, 1733.926, -39.485
        elif a_T>=333 and a_T<363:
            A, B, C = 5.0768,  1659.793, -45.854
        elif a_T>=363 and a_T<379:
            A, B, C = 5.08354, 1663.125, -45.622
        elif a_T>=379 and a_T<=573:
            A, B, C = 3.55959,  643.748, -198.043
        else:
            errorMessage = 'The temperature of particle '+str(n)+' is '+str(a_T)+a_T_Units+'. Cannot find partial pressures when particle temperature is >573 K'
            raise ValueError(errorMessage)

        p_sat  = (10.0**(A-(B/(a_T+C)))) * 100000 #Pa

        if self.m_InputData.m_PressUnits == 'psi':
            p_sat *= 0.0001450377 #psi

        return p_sat

    def calculateMassFraction(self):
        rho_p       = self.m_InputData.m_TracerDensity
        relHum      = self.m_InputData.m_Humidity
        p_satAvg    = self.m_AverageSaturationPressure
        p_avg       = self.m_AveragePressure
        airMass     = self.m_InputData.m_FluidDensity * self.m_InputData.m_FluidVolume
        vaporMass   = airMass*0.622*relHum*p_satAvg/(p_avg-relHum*p_satAvg)
        dropletMass = 0

        for p in range(self.numParticles):
            dropletRad   = self.m_LagrangianData.getScalarData(a_ArrayName='Radius', a_DataID=p)
            dropletMass += (rho_p*4.0*np.pi/3.0)*(dropletRad**3)

        self.m_WaterMassFraction = (vaporMass+dropletMass)/(airMass+vaporMass+dropletMass)

    def calculatePressureConversionFactors(self):
        self.m_PressureConversionFactor = self.m_InputData.m_PressConversion

        if self.m_InputData.m_AbsRelPress == 'rel':
            if self.m_InputData.m_PressUnits == 'Pa':
                self.m_AbsPressureConversion = 101325
            elif self.m_InputData.m_PressUnits!='psi':
                self.m_AbsPressureConversion = 14.696
            else:
                sys.exit('Pressure units must be Pa or psi. If system units are not one of these, use a conversion factor.')
        else:
            self.m_AbsPressureConversion = 0.0

    def calculateFluidThermalConductivity(self, a_T, a_T_Units):
        if a_T_Units == 'C':
            a_T += 273.1
            if a_T<100:
                raise ValueError('No data available for thermal conductivity of air at temperatures less than 100K')
            elif a_T>600:
                raise ValueError('No data available for thermal conductivity of air at temperatures greater than 600K')
            kArray = np.array([9.34, 13.8, 18.1, 22.3, 26.3, 30.0, 33.8, 37.3, 40.7, 43.9, 46.9])*1e-6 #kW/m-K
            TArray = np.array([ 100,  150,  200,  250,  300,  350,  400,  450,  500,  550,  600]) #K

        elif a_T_Units == 'K':
            if a_T<100:
                raise ValueError('No data available for thermal conductivity of air at temperatures less than 100K')
            elif a_T>600:
                raise ValueError('No data available for thermal conductivity of air at temperatures greater than 600K')
            kArray = np.array([9.34, 13.8, 18.1, 22.3, 26.3, 30.0, 33.8, 37.3, 40.7, 43.9, 46.9])*1e-6 #kW/m-K
            TArray = np.array([ 100,  150,  200,  250,  300,  350,  400,  450,  500,  550,  600]) #K

        elif a_T_Units == 'F':
            a_T += 460.67
            if a_T<valLow:
                raise ValueError('No data available for thermal conductivity of air at temperatures less than 180°R')
            elif a_T>valUp:
                raise ValueError('No data available for thermal conductivity of air at temperatures greater than 1080°R')
            kArray = np.array([0.00548, 0.00818, 0.01069, 0.01303, 0.01524, 0.01734, 0.01933, 0.02124, 0.02308, 0.02486, 0.02658])/3600.0 #BTU/s-ft-R
            TArray = np.array([    180,     270,     360,     450,     540,     630,     720,     810,     900,     990,    1080]) #R

        elif a_T_Units == 'R':
            if a_T<valLow:
                raise ValueError('No data available for thermal conductivity of air at temperatures less than 180°R')
            elif a_T>valUp:
                raise ValueError('No data available for thermal conductivity of air at temperatures greater than 1080°R')
            kArray = np.array([0.00548, 0.00818, 0.01069, 0.01303, 0.01524, 0.01734, 0.01933, 0.02124, 0.02308, 0.02486, 0.02658])/3600.0 #BTU/s-ft-R
            TArray = np.array([    180,     270,     360,     450,     540,     630,     720,     810,     900,     990,    1080]) #R

        indUp  = np.searchsorted(TArray, a_T)
        indLow = indUp-1

        kLow   = kArray[indLow]
        kUp    = kArray[indUp]
        TLow   = TArray[indLow]
        TUp    = TArray[indUp]

        k_p    = kLow + (kUp-kLow)*(a_T-TLow)/(TUp-TLow)
        self.m_FluidThermalConductivity = k_p

    def calculateStokesNumber(self, a_ID):
        self.m_StokesNumber[a_ID] = 0.5*self.m_MomentumResponseTime[a_ID]*self.m_AverageFlowVelocity/self.m_LagrangianData.getScalarData(a_ArrayName='Radius', a_DataID=a_ID)

    def calculateSlipVelocityReynoldsNumber(self, vel, a_ID):
        v_p      = self.m_LagrangianData.getVectorData(a_ArrayName='Velocity', a_DataID=a_ID)
        vel_norm = np.sqrt(( vel[0]-v_p[0])**2 + (vel[1]-v_p[1])**2 + (vel[2]-v_p[2])**2)
        rho      = self.m_InputData.getTracerDensity()
        mu       = self.m_InputData.getFluidViscosity()
        d_p      = self.m_LagrangianData.getScalarData(a_ArrayName='Radius', a_DataID=a_ID)*2.0
        self.m_SlipVelocityReynoldsNumber[a_ID] = rho*vel_norm*d_p/mu



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
            TracerRadius = np.ones(self.numParticles)*self.m_InputData.m_TracerRadius
        else:
            TracerRadius = np.random.weibull(self.m_InputData.m_ShapeParam,\
             self.m_LagrangianData.m_NumParticles) * self.m_InputData.m_AverageRadius
        self.m_LagrangianData.addScalarData('Radius', TracerRadius)

        #---------------------------------------------
        # calculate equilibrium (final) particle size
        #---------------------------------------------
        self.m_TracerEquilibriumRadius = TracerRadius * (88.0/1000.0)**(1/3)

        #--------------------------------------
        # calculate momentum response times(s)
        #--------------------------------------
        self.m_MomentumResponseTime = self.m_InputData.getTracerDensity()*((2.0*TracerRadius)**2)/(18.0*self.m_InputData.getFluidViscosity())
        self.m_LagrangianData.addScalarData('Momentum Response Time', self.m_MomentumResponseTime)

        #-------------------------------------------------------------------
        # check integration time step based on inputs and warn if too large
        #-------------------------------------------------------------------
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

        #----------------------------------------------------
        # check for steady flow and set tWin_2 appropriately
        #----------------------------------------------------
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
                elif check_start and not check_end:
                    sys.exit('Inputs seem to indicate that OpenFOAM data is being used but no OpenFOAM file is found in the specified flow data directory.\
                    Double check the flow data directory given as well as the related inputs.')

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
        InitialVelocity      = np.zeros((self.numParticles, 3))
        InitialVelocity[:,0] = float(self.m_InputData.m_InitialVelocity['u'])
        InitialVelocity[:,1] = float(self.m_InputData.m_InitialVelocity['v'])
        InitialVelocity[:,2] = float(self.m_InputData.m_InitialVelocity['w'])
        self.m_LagrangianData.addVectorData('Velocity', InitialVelocity)

        #--------------------------------------
        # set initial temperature of particles
        #--------------------------------------
        InitialTemp    = np.zeros(self.numParticles)
        InitialTemp[:] = self.m_InputData.m_InitialTemp
        self.m_LagrangianData.addScalarData('Temperature', InitialTemp)

        #---------------------------------
        # calculate evaporation factor(s)
        #---------------------------------
        self.m_EvaporationFactor = 6.0*2.0*self.m_InputData.m_FluidDensity*\
         self.m_InputData.m_TracerDiffCoeff/(2.0*TracerRadius*self.m_InputData.m_TracerDensity)

        self.m_LagrangianData.addScalarData('Evaporation Factor', self.m_EvaporationFactor)

        #------------------------------------
        # create pressure conversion factors
        #------------------------------------
        self.calculatePressureConversionFactors()

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
            mrtStruct = reader.GetOutput().GetPointData().GetArray('Momentum Response Time')
            tmpStruct = reader.GetOutput().GetPointData().GetArray('Temperature')
            evfStruct = reader.GetOutput().GetPointData().GetArray('Evaporation Factor')

            for p in range(self.numParticles):

                point = reader.GetOutput().GetPoint(p)
                self.m_LagrangianData.setX(p, point)

                velocity = np.array(velStruct.GetTuple3(p))
                self.m_LagrangianData.setVectorData(velocity, a_ArrayName='Velocity', a_DataID=p)

                radius = radStruct.GetTuple1(p)
                self.m_LagrangianData.setScalarData(radius, a_ArrayName='Radius', a_DataID=p)

                self.m_MomentumResponseTime[p] = mrtStruct.GetTuple1(p)
                self.m_LagrangianData.setScalarData(self.m_MomentumResponseTime[p], a_ArrayName='Momentum Response Time', a_DataID=p)

                temp = tmpStruct.GetTuple1(p)
                self.m_LagrangianData.setScalarData(temp, a_ArrayName='Temperature', a_DataID=p)

                self.m_EvaporationFactor[p] = evfStruct.GetTuple1(p)
                self.m_LagrangianData.setScalarData(self.m_EvaporationFactor[p], a_ArrayName='Evaporation Factor', a_DataID=p)

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
                    self.m_GridData.renameDataArray(self.m_InputData.m_TemperatureFieldName, 'T0')
                    self.m_GridData.renameDataArray(self.m_InputData.m_PressureFieldName, 'p0')

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
                        self.m_GridData.renameDataArray(self.m_InputData.m_TemperatureFieldName, 'T0')
                        self.m_GridData.renameDataArray(self.m_InputData.m_PressureFieldName, 'p0')

                    elif self.m_ResumeSimulation:

                        if isDataPointSync:

                            flowInitial = self.m_InputData.getFlowDataFileName(a_ID = timeWindowDict['ID_Low'][timeIndex])
                            self.m_GridData.getDataFromFile(flowInitial)
                            self.m_GridData.renameDataArray(self.m_InputData.getVelDataName(), 'v0')
                            self.m_GridData.renameDataArray(self.m_InputData.m_TemperatureFieldName, 'T0')
                            self.m_GridData.renameDataArray(self.m_InputData.m_PressureFieldName, 'p0')

                        else:

                            if isSingleStageIntegration == True:

                                flowSingle = self.m_InputData.getFlowDataFileName(a_ID = timeWindowDict['ID_Low'][timeIndex])
                                self.m_GridData.addDataArrayFromFile(flowSingle, 'v0', self.m_InputData.getVelDataName())

                            if simTime < tWin_1 - dT:

                                flowPlus        = self.m_InputData.getFlowDataFileName(a_ID = timeWindowDict['ID_Up'][timeIndex])
                                flowMinus       = self.m_InputData.getFlowDataFileName(a_ID = timeWindowDict['ID_Low'][timeIndex])
                                self.m_GridData.addDataArrayFromFile(flowPlus, 'v1', self.m_InputData.getVelDataName())
                                self.m_GridData.addDataArrayFromFile(flowMinus, 'v0', self.m_InputData.getVelDataName())
                                self.m_GridData.addDataArrayFromFile(flowPlus, 'T1', self.m_InputData.m_TemperatureFieldName)
                                self.m_GridData.addDataArrayFromFile(flowMinus, 'T0', self.m_InputData.m_TemperatureFieldName)
                                self.m_GridData.addDataArrayFromFile(flowPlus, 'p0', self.m_InputData.m_PressureFieldName)
                                self.m_GridData.addDataArrayFromFile(flowMinus, 'p1', self.m_InputData.m_PressureFieldName)
                                print('Loaded a new velocity data')

                            elif simTime >= tWin_1 - dT:

                                flowPlus    = self.m_InputData.getFlowDataFileName(a_ID = timeWindowDict['ID_Up'][timeIndex])
                                flowMinus   = self.m_InputData.getFlowDataFileName(a_ID = timeWindowDict['ID_Low'][timeIndex])
                                flowNext    = self.m_InputData.getFlowDataFileName(a_ID = timeWindowDict['ID_Low'][timeIndex]+1)
                                self.m_GridData.addDataArrayFromFile(flowPlus, 'v1', self.m_InputData.getVelDataName())
                                self.m_GridData.addDataArrayFromFile(flowMinus, 'v0', self.m_InputData.getVelDataName())
                                self.m_GridData.addDataArrayFromFile(flowNext, 'v2', self.m_InputData.getVelDataName())
                                self.m_GridData.addDataArrayFromFile(flowPlus, 'T1', self.m_InputData.m_TemperatureFieldName)
                                self.m_GridData.addDataArrayFromFile(flowMinus, 'T0', self.m_InputData.m_TemperatureFieldName)
                                self.m_GridData.addDataArrayFromFile(flowNext, 'T2', self.m_InputData.m_TemperatureFieldName)
                                self.m_GridData.addDataArrayFromFile(flowPlus, 'p0', self.m_InputData.m_PressureFieldName)
                                self.m_GridData.addDataArrayFromFile(flowMinus, 'p1', self.m_InputData.m_PressureFieldName)
                                self.m_GridData.addDataArrayFromFile(flowNext, 'p2', self.m_InputData.m_PressureFieldName)
                                print('Loaded a new velocity data')

                    elif isDataPointSync:

                        if isSingleStageIntegration == True:

                            self.m_GridData.removeDataArray('v0')
                            flowSingle = self.m_InputData.getFlowDataFileName(a_ID = timeWindowDict['ID_Low'][timeIndex])
                            self.m_GridData.addDataArrayFromFile(flowSingle, 'v0', self.m_InputData.getVelDataName())

                        else:

                            self.m_GridData.removeDataArray('v0')
                            self.m_GridData.removeDataArray('v1')
                            self.m_GridData.removeDataArray('T0')
                            self.m_GridData.removeDataArray('T1')
                            self.m_GridData.removeDataArray('P0')
                            self.m_GridData.removeDataArray('P1')
                            flowPlus  = self.m_InputData.getFlowDataFileName(a_ID = timeWindowDict['ID_Up'][timeIndex])
                            flowMinus = self.m_InputData.getFlowDataFileName(a_ID = timeWindowDict['ID_Low'][timeIndex])
                            self.m_GridData.addDataArrayFromFile(flowPlus, 'v1', self.m_InputData.getVelDataName())
                            self.m_GridData.addDataArrayFromFile(flowMinus, 'v0', self.m_InputData.getVelDataName())
                            self.m_GridData.addDataArrayFromFile(flowPlus, 'T1', self.m_InputData.m_TemperatureFieldName)
                            self.m_GridData.addDataArrayFromFile(flowMinus, 'T0', self.m_InputData.m_TemperatureFieldName)
                            self.m_GridData.addDataArrayFromFile(flowPlus, 'p0', self.m_InputData.m_PressureFieldName)
                            self.m_GridData.addDataArrayFromFile(flowMinus, 'p1', self.m_InputData.m_PressureFieldName)
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
                                self.m_GridData.removeDataArray('T0')
                                self.m_GridData.removeDataArray('T1')
                                self.m_GridData.removeDataArray('T2')
                                self.m_GridData.removeDataArray('p0')
                                self.m_GridData.removeDataArray('p1')
                                self.m_GridData.removeDataArray('p2')
                                flowPlus  = self.m_InputData.getFlowDataFileName(a_ID = timeWindowDict['ID_Up'][timeIndex])
                                flowMinus = self.m_InputData.getFlowDataFileName(a_ID = timeWindowDict['ID_Low'][timeIndex])
                                self.m_GridData.addDataArrayFromFile(flowPlus, 'v1', self.m_InputData.getVelDataName())
                                self.m_GridData.addDataArrayFromFile(flowMinus, 'v0', self.m_InputData.getVelDataName())
                                self.m_GridData.addDataArrayFromFile(flowPlus, 'T1', self.m_InputData.m_TemperatureFieldName)
                                self.m_GridData.addDataArrayFromFile(flowMinus, 'T0', self.m_InputData.m_TemperatureFieldName)
                                self.m_GridData.addDataArrayFromFile(flowPlus, 'p0', self.m_InputData.m_PressureFieldName)
                                self.m_GridData.addDataArrayFromFile(flowMinus, 'p1', self.m_InputData.m_PressureFieldName)
                                print('Loaded a new velocity data')

                            elif (simTime + self.m_dT) >= tWin_1:

                                self.m_GridData.removeDataArray('v0')
                                self.m_GridData.removeDataArray('v1')
                                self.m_GridData.removeDataArray('v2')
                                self.m_GridData.removeDataArray('T0')
                                self.m_GridData.removeDataArray('T1')
                                self.m_GridData.removeDataArray('T2')
                                self.m_GridData.removeDataArray('p0')
                                self.m_GridData.removeDataArray('p1')
                                self.m_GridData.removeDataArray('p2')
                                flowPlus  = self.m_InputData.getFlowDataFileName(a_ID = timeWindowDict['ID_Up'][timeIndex])
                                flowMinus = self.m_InputData.getFlowDataFileName(a_ID = timeWindowDict['ID_Low'][timeIndex])
                                flowNext  = self.m_InputData.getFlowDataFileName(a_ID = timeWindowDict['ID_Low'][timeIndex]+1)
                                self.m_GridData.addDataArrayFromFile(flowPlus, 'v1', self.m_InputData.getVelDataName())
                                self.m_GridData.addDataArrayFromFile(flowMinus, 'v0', self.m_InputData.getVelDataName())
                                self.m_GridData.addDataArrayFromFile(flowNext, 'v2', self.m_InputData.getVelDataName())
                                self.m_GridData.addDataArrayFromFile(flowPlus, 'T1', self.m_InputData.m_TemperatureFieldName)
                                self.m_GridData.addDataArrayFromFile(flowMinus, 'T0', self.m_InputData.m_TemperatureFieldName)
                                self.m_GridData.addDataArrayFromFile(flowNext, 'T2', self.m_InputData.m_TemperatureFieldName)
                                self.m_GridData.addDataArrayFromFile(flowPlus, 'p0', self.m_InputData.m_PressureFieldName)
                                self.m_GridData.addDataArrayFromFile(flowMinus, 'p1', self.m_InputData.m_PressureFieldName)
                                self.m_GridData.addDataArrayFromFile(flowNext, 'p2', self.m_InputData.m_PressureFieldName)
                                print('Loaded a new velocity data')

            #-----------------------------------------------------------------------------------------------
            # STEP 5: build a cell locator (if it is a fixed mesh, built only once) and set initial velocity
            #-----------------------------------------------------------------------------------------------
            if (simTime == startTime and self.m_InputData.isFixedMesh()) or (self.m_ResumeSimulation and self.m_InputData.isFixedMesh()):
                print("Building Cell Locator Maps")
                self.m_GridData.buildLocator()

            #-------------------------------------------------------------------
            # STEP 6: Output initial position. Calculate average flow velocity,
            # Stokes Number, Slip Velocity Reynolds Number
            #-------------------------------------------------------------------
            if simTime == startTime or self.m_ResumeSimulation:

                #-----------------------------------------------------------------
                # calculate average values and the initial mass fraction of water
                #-----------------------------------------------------------------
                print('Calculating Average Field Values')
                self.calculateAverageVelocity()
                self.calculateAverageTemperature()
                self.calculateAveragePressure()
                self.calculateMassFraction()

                #-------------------------------------------------------------------
                # calculate thermal conductivity of the fluid based on average temp
                # calculate thermal response times of the particles
                #-------------------------------------------------------------------
                self.calculateFluidThermalConductivity(self.m_AverageTemperature, self.m_InputData.m_TempUnits)
                self.m_ThermalResponseTime = self.m_InputData.m_TracerDensity*self.m_InputData.m_SpecificHeat*((TracerRadius*2.0)**2)/(12.0*self.m_FluidThermalConductivity)
                self.m_LagrangianData.addScalarData('Thermal Response Time', self.m_ThermalResponseTime)

                #---------------------------------
                # add arrays to integration class
                #---------------------------------
                self.m_StokesNumber = np.zeros(self.numParticles, dtype=np.float32)
                self.m_SlipVelocityReynoldsNumber = np.zeros(self.numParticles, dtype=np.float32)

                #-----------------------------------------------------------
                # calculate stokes number and slip velocity reynolds number
                #-----------------------------------------------------------
                for p in range(self.numParticles):
                    self.calculateStokesNumber(p)
                    if True: #a_PolygonalCells:
                        vel_i = self.m_GridData.gridInterpolateAveraging(self.m_LagrangianData.getX(p), 'v0')
                    else:
                        vel_i = self.m_GridData.gridInterpolateNodalBasis(self.m_LagrangianData.getX(p), 'v0')
                    self.calculateSlipVelocityReynoldsNumber(vel_i, p)

                #------------------------------------------------------
                # add info to LagrangianData and output initial values
                #------------------------------------------------------
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

            #-------------------------------------------------
            # STEP 8: update time indices and simulation time
            #-------------------------------------------------
            timeIndex += 1
            simTime    = round(simTime+self.m_dT, self.m_dt_prec)

            #-----------------------------------------------------------------------
            # STEP 9: at the end of appropriate number of steps dump data into file
            #-----------------------------------------------------------------------
            if timeIndex%self.m_InputData.m_WriteInterval == 0:
                self.m_LagrangianData.writeData(self.m_InputData.getTracerOutputFile(a_ID1=0, a_ID2=timeIndex))

            #--------------------------------------------
            # STEP 10: calculate new water mass fraction
            #--------------------------------------------
            self.calculateMassFraction()



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
        dT       = self.m_dT
        buoy     = self.m_Grav*(1.0-self.m_InputData.getFluidDensity()/self.m_InputData.getTracerDensity())
        molWFrac = 18.0/29.0
        w_a      = self.m_WaterMassFraction

        #----------------------------
        # calculate average velocity
        #----------------------------
        if not self.m_isSteadyFlow and a_DataSync:
            self.calculateAverageVelocity()

        for p in range(self.numParticles):

            tauV  = self.m_MomentumResponseTime[p]
            tauT  = self.m_ThermalResponseTime[p]
            beta  = self.m_EvaporationFactor[p]
            xyz_i = self.m_LagrangianData.getX(p)
            v_i   = np.array(self.m_LagrangianData.getVectorData(a_ArrayName='Velocity', a_DataID=p))
            T_i   = self.m_LagrangianData.getScalarData(a_ArrayName='Temperature', a_DataID=p)
            D_i   = self.m_LagrangianData.getScalarData(a_ArrayName='Radius', a_DataID=p)*2.0
            R_eq  = self.m_TracerEquilibriumRadius[p]

            if a_DataSync:
                if a_PolygonalCells == True:
                    vel, stat = self.m_GridData.gridInterpolateAveraging(xyz_i, 'v0', a_GetStatus=True)
                    T         = self.m_GridData.gridInterpolateAveraging(xyz_i, 'T0', a_DataDim=1)
                    press     = self.m_GridData.gridInterpolateAveraging(xyz_i, 'p0', a_DataDim=1)*self.m_PressureConversionFactor + self.m_AbsPressureConversion
                else:
                    vel, stat = self.m_GridData.gridInterpolateNodalBasis(xyz_i, 'v0', a_GetStatus=True)
                    T         = self.m_GridData.gridInterpolateNodalBasis(xyz_i, 'T0', a_DataDim=1)
                    press     = self.m_GridData.gridInterpolateNodalBasis(xyz_i, 'p0', a_DataDim=1)*self.m_PressureConversionFactor + self.m_AbsPressureConversion
            else:
                if a_PolygonalCells == True:
                    velPlus, stat  = self.m_GridData.gridInterpolateAveraging(xyz_i, 'v1', a_GetStatus=True)
                    velMinus, stat = self.m_GridData.gridInterpolateAveraging(xyz_i, 'v0', a_GetStatus=True)
                    TPlus          = self.m_GridData.gridInterpolateAveraging(xyz_i, 'T1', a_DataDim=1)
                    TMinus         = self.m_GridData.gridInterpolateAveraging(xyz_i, 'T0', a_DataDim=1)
                    pPlus          = self.m_GridData.gridInterpolateAveraging(xyz_i, 'p1', a_DataDim=1)*self.m_PressureConversionFactor + self.m_AbsPressureConversion
                    pMinus         = self.m_GridData.gridInterpolateAveraging(xyz_i, 'p0', a_DataDim=1)*self.m_PressureConversionFactor + self.m_AbsPressureConversion
                else:
                    velPlus, stat  = self.m_GridData.gridInterpolateNodalBasis(xyz_i, 'v1', a_GetStatus=True)
                    velMinus, stat = self.m_GridData.gridInterpolateNodalBasis(xyz_i, 'v0', a_GetStatus=True)
                    TPlus          = self.m_GridData.gridInterpolateNodalBasis(xyz_i, 'T1', a_DataDim=1)
                    TMinus         = self.m_GridData.gridInterpolateNodalBasis(xyz_i, 'T0', a_DataDim=1)
                    pPlus          = self.m_GridData.gridInterpolateNodalBasis(xyz_i, 'p1', a_DataDim=1)*self.m_PressureConversionFactor + self.m_AbsPressureConversion
                    pMinus         = self.m_GridData.gridInterpolateNodalBasis(xyz_i, 'p0', a_DataDim=1)*self.m_PressureConversionFactor + self.m_AbsPressureConversion

                vel   = velMinus + (simTime - tWin_0) * (velPlus - velMinus)/(tWin_1 - tWin_0)
                T     = TMinus   + (simTime - tWin_0) * (TPlus   - TMinus)  /(tWin_1 - tWin_0)
                press = pMinus   + (simTime - tWin_0) * (pPlus   - pMinus)  /(tWin_1 - tWin_0)

            #-----------------------------------------------------
            # forward euler integration for velocity and position
            #-----------------------------------------------------
            if stat != -1:
                v_new   = v_i + dT*((vel-v_i)/tauV + buoy)
                xyz_new = xyz_i + v_new*dT
                T_new   = T_i + dT*(T-T_i)/tauT
                R_new   = 0.5*(D_i + dT*beta*(w_a-(molWFrac*self.calculateSaturationPressure(T_i, self.m_InputData.m_TempUnits, p)/press)))

                #--------------------------------------
                # check relation to equilibrium radius
                #--------------------------------------
                if R_new < R_eq:
                    R_new = R_eq

            else:
                v_new   = np.zeros(self.m_InputData.m_SpaceDimension)
                xyz_new = xyz_i
                T_new   = T_i
                R_new   = D_i*0.5

            #--------------------------------------------------
            # output new particle properties to LagrangianData
            #--------------------------------------------------
            self.m_LagrangianData.setX(p, xyz_new)
            self.m_LagrangianData.setVectorData(v_new, a_ArrayName='Velocity', a_DataID=p)
            self.m_LagrangianData.setScalarData(T_new, a_ArrayName='Temperature', a_DataID=p)
            self.m_LagrangianData.setScalarData(R_new, a_ArrayName='Radius', a_DataID=p)

            #---------------------------------------------
            # output new physics values to LagrangianData
            #---------------------------------------------
            self.m_MomentumResponseTime[p] = self.calculateMomentumResponseTime(R_new)
            self.m_ThermalResponseTime[p]  = self.calculateThermalResponseTime(R_new)
            self.m_EvaporationFactor[p]    = self.calculateEvaporationFactor(R_new)
            self.m_LagrangianData.setScalarData(self.m_MomentumResponseTime[p], a_ArrayName='Momentum Response Time', a_DataID=p)
            self.m_LagrangianData.setScalarData(self.m_ThermalResponseTime[p],  a_ArrayName='Thermal Response Time',  a_DataID=p)
            self.m_LagrangianData.setScalarData(self.m_EvaporationFactor[p],    a_ArrayName='Evaporation Factor',     a_DataID=p)

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
        dT        = self.m_dT
        buoy      = self.m_Grav*(1.0-self.m_InputData.getFluidDensity()/self.m_InputData.getTracerDensity())
        molWFrac  = 18.0/29.0
        w_a       = self.m_WaterMassFraction

        #----------------------------
        # calculate average velocity
        #----------------------------
        if not self.m_isSteadyFlow and a_DataSync:
            self.calculateAverageVelocity()

        for p in range(self.numParticles):

            xyz_1  = self.m_LagrangianData.getX(p)
            v_1    = np.array(self.m_LagrangianData.getVectorData(a_ArrayName='Velocity', a_DataID=p))
            T_p1   = self.m_LagrangianData.getScalarData(a_ArrayName='Temperature', a_DataID=p)
            D_1    = self.m_LagrangianData.getScalarData(a_ArrayName='Radius', a_DataID=p)*2.0
            tauV_1 = self.m_MomentumResponseTime[p]
            tauT_1 = self.m_ThermalResponseTime[p]
            beta_1 = self.m_EvaporationFactor[p]
            D_eq   = self.m_TracerEquilibriumRadius[p]*2.0

            if a_DataSync:

                if a_PolygonalCells:

                    vel_1, stat = self.m_GridData.gridInterpolateAveraging(xyz_1, 'v0', a_GetStatus=True)
                    T_air1      = self.m_GridData.gridInterpolateAveraging(xyz_1, 'T0', a_DataDim=1)
                    p_1         = self.m_GridData.gridInterpolateAveraging(xyz_1, 'p0', a_DataDim=1)*self.m_PressureConversionFactor + self.m_AbsPressureConversion

                    if stat != -1:

                        k_1     = v_1
                        el_1    = (vel_1-v_1)/tauV_1 + buoy
                        m_1     = (T_air1-T_p1)/tauT_1
                        n_1     = beta_1*(w_a-(molWFrac * self.calculateSaturationPressure(T_p1, self.m_InputData.m_TempUnits, p)/p_1))

                        xyz_2   = xyz_1 + 0.5*k_1*dT
                        v_2     = v_1   + 0.5*el_1*dT
                        T_p2    = T_p1  + 0.5*m_1*dT
                        D_2     = D_1   + 0.5*n_1*dT if D_1>D_eq else D_eq
                        tauV_2  = self.calculateMomentumResponseTime(D_2*0.5)
                        tauT_2  = self.calculateThermalResponseTime(D_2*0.5)
                        beta_2  = self.calculateEvaporationFactor(D_2*0.5)

                        vel_2   = self.m_GridData.gridInterpolateAveraging(xyz_2, 'v0')
                        T_air2  = self.m_GridData.gridInterpolateAveraging(xyz_2, 'T0', a_DataDim=1)
                        p_2     = self.m_GridData.gridInterpolateAveraging(xyz_2, 'p0', a_DataDim=1)*self.m_PressureConversionFactor + self.m_AbsPressureConversion

                        T_air2  = T_air2 if T_air2!=0 else T_air1

                        k_2     = v_2
                        el_2    = (vel_2-v_2)/tauV_2 + buoy
                        m_2     = (T_air2-T_p2)/tauT_2
                        n_2     = beta_2*(w_a-(molWFrac * self.calculateSaturationPressure(T_p2, self.m_InputData.m_TempUnits, p)/p_2))

                        xyz_3   = xyz_1 + 0.5*k_2*dT
                        v_3     = v_1   + 0.5*el_2*dT
                        T_p3    = T_p1  + 0.5*m_2*dT
                        D_3     = D_1   + 0.5*n_2*dT if D_1>D_eq else D_eq
                        tauV_3  = self.calculateMomentumResponseTime(D_3*0.5)
                        tauT_3  = self.calculateThermalResponseTime(D_3*0.5)
                        beta_3  = self.calculateEvaporationFactor(D_3*0.5)

                        vel_3   = self.m_GridData.gridInterpolateAveraging(xyz_3, 'v0')
                        T_air3  = self.m_GridData.gridInterpolateAveraging(xyz_3, 'T0', a_DataDim=1)
                        p_3     = self.m_GridData.gridInterpolateAveraging(xyz_3, 'p0', a_DataDim=1)*self.m_PressureConversionFactor + self.m_AbsPressureConversion

                        T_air3  = T_air3 if T_air3!=0 else T_air2

                        k_3     = v_3
                        el_3    = (vel_3-v_3)/tauV_3 + buoy
                        m_3     = (T_air3-T_p3)/tauT_3
                        n_3     = beta_3*(w_a-(molWFrac * self.calculateSaturationPressure(T_p3, self.m_InputData.m_TempUnits, p)/p_3))

                        xyz_4   = xyz_1 + k_3*dT
                        v_4     = v_1   + 0.5*el_3*dT
                        T_p4    = T_p1  + 0.5*m_3*dT
                        D_4     = D_1   + 0.5*n_3*dT if D_1>D_eq else D_eq
                        tauV_4  = self.calculateMomentumResponseTime(D_4*0.5)
                        tauT_4  = self.calculateThermalResponseTime(D_4*0.5)
                        beta_4  = self.calculateEvaporationFactor(D_4*0.5)

                        vel_4   = self.m_GridData.gridInterpolateAveraging(xyz_4, 'v0')
                        T_air4  = self.m_GridData.gridInterpolateAveraging(xyz_4, 'T0', a_DataDim=1)
                        p_4     = self.m_GridData.gridInterpolateAveraging(xyz_4, 'p0', a_DataDim=1)*self.m_PressureConversionFactor + self.m_AbsPressureConversion

                        k_4     = v_4
                        el_4    = (vel_4-v_4)/tauV_4 + buoy
                        m_4     = (T_air4-T_p4)/tauT_4
                        n_4     = beta_4*(w_a-(molWFrac * self.calculateSaturationPressure(T_p4, self.m_InputData.m_TempUnits, p)/p_4))

                        v_new   = v_1   + (dT/6)*(el_1 + 2*el_2 + 2*el_3 + el_4)
                        xyz_new = xyz_1 + (dT/6)*(k_1  + 2*k_2  + 2*k_3  + k_4)
                        T_new   = T_p1  + (dT/6)*(m_1  + 2*m_2  + 2*m_3  + m_4) if T_air4!=0 else T_p1
                        D_new   = D_1   + (dT/6)*(n_1  + 2*n_2  + 2*n_3  + n_4) if D_1>D_eq else D_eq

                    else:

                        xyz_new = xyz_1
                        v_new   = np.zeros(self.m_InputData.m_SpaceDimension)
                        T_new   = T_p1
                        D_new   = D_1

                else:

                    vel_1, stat = self.m_GridData.gridInterpolateNodalBasis(xyz_1, 'v0', a_GetStatus=True)
                    T_air1      = self.m_GridData.gridInterpolateNodalBasis(xyz_1, 'T0', a_DataDim=1)
                    p_1         = self.m_GridData.gridInterpolateNodalBasis(xyz_1, 'p0', a_DataDim=1)*self.m_PressureConversionFactor + self.m_AbsPressureConversion

                    if stat != -1:

                        k_1     = v_1
                        el_1    = (vel_1-v_1)/tauV_1 + buoy
                        m_1     = (T_air1-T_p1)/tauT_1
                        n_1     = beta_1*(w_a-(molWFrac * self.calculateSaturationPressure(T_p1, self.m_InputData.m_TempUnits, p)/p_1))

                        xyz_2   = xyz_1 + 0.5*k_1*dT
                        v_2     = v_1   + 0.5*el_1*dT
                        T_p2    = T_p1  + 0.5*m_1*dT
                        D_2     = D_1   + 0.5*n_1*dT if D_1>D_eq else D_eq
                        tauV_2  = self.calculateMomentumResponseTime(D_2*0.5)
                        tauT_2  = self.calculateThermalResponseTime(D_2*0.5)
                        beta_2  = self.calculateEvaporationFactor(D_2*0.5)

                        vel_2   = self.m_GridData.gridInterpolateNodalBasis(xyz_2, 'v0')
                        T_air2  = self.m_GridData.gridInterpolateNodalBasis(xyz_2, 'T0', a_DataDim=1)
                        p_2     = self.m_GridData.gridInterpolateNodalBasis(xyz_2, 'p0', a_DataDim=1)*self.m_PressureConversionFactor + self.m_AbsPressureConversion

                        T_air2  = T_air2 if T_air2!=0 else T_air1

                        k_2     = v_2
                        el_2    = (vel_2-v_2)/tauV_2 + buoy
                        m_2     = (T_air2-T_p2)/tauT_2
                        n_2     = beta_2*(w_a-(molWFrac * self.calculateSaturationPressure(T_p2, self.m_InputData.m_TempUnits, p)/p_2))

                        xyz_3   = xyz_1 + 0.5*k_2*dT
                        v_3     = v_1   + 0.5*el_2*dT
                        T_p3    = T_p1  + 0.5*m_2*dT
                        D_3     = D_1   + 0.5*n_2*dT if D_1>D_eq else D_eq
                        tauV_3  = self.calculateMomentumResponseTime(D_3*0.5)
                        tauT_3  = self.calculateThermalResponseTime(D_3*0.5)
                        beta_3  = self.calculateEvaporationFactor(D_3*0.5)

                        vel_3   = self.m_GridData.gridInterpolateNodalBasis(xyz_3, 'v0')
                        T_air3  = self.m_GridData.gridInterpolateNodalBasis(xyz_3, 'T0', a_DataDim=1)
                        p_3     = self.m_GridData.gridInterpolateNodalBasis(xyz_3, 'p0', a_DataDim=1)*self.m_PressureConversionFactor + self.m_AbsPressureConversion

                        T_air3  = T_air3 if T_air3!=0 else T_air2

                        k_3     = v_3
                        el_3    = (vel_3-v_3)/tauV_3 + buoy
                        m_3     = (T_air3-T_p3)/tauT_3
                        n_3     = beta_3*(w_a-(molWFrac * self.calculateSaturationPressure(T_p3, self.m_InputData.m_TempUnits, p)/p_3))

                        xyz_4   = xyz_1 + k_3*dT
                        v_4     = v_1   + 0.5*el_3*dT
                        T_p4    = T_p1  + 0.5*m_3*dT
                        D_4     = D_1   + 0.5*n_3*dT if D_1>D_eq else D_eq
                        tauV_4  = self.calculateMomentumResponseTime(D_4*0.5)
                        tauT_4  = self.calculateThermalResponseTime(D_4*0.5)
                        beta_4  = self.calculateEvaporationFactor(D_4*0.5)

                        vel_4   = self.m_GridData.gridInterpolateNodalBasis(xyz_4, 'v0')
                        T_air4  = self.m_GridData.gridInterpolateNodalBasis(xyz_4, 'T0', a_DataDim=1)
                        p_4     = self.m_GridData.gridInterpolateNodalBasis(xyz_4, 'p0', a_DataDim=1)*self.m_PressureConversionFactor + self.m_AbsPressureConversion

                        k_4     = v_4
                        el_4    = (vel_4-v_4)/tauV_4 + buoy
                        m_4     = (T_air4-T_p4)/tauT_4
                        n_4     = beta_4*(w_a-(molWFrac * self.calculateSaturationPressure(T_p4, self.m_InputData.m_TempUnits, p)/p_4))

                        v_new   = v_1   + (dT/6)*(el_1 + 2*el_2 + 2*el_3 + el_4)
                        xyz_new = xyz_1 + (dT/6)*(k_1  + 2*k_2  + 2*k_3  + k_4)
                        T_new   = T_p1  + (dT/6)*(m_1  + 2*m_2  + 2*m_3  + m_4) if T_air4!=0 else T_p1
                        D_new   = D_1   + (dT/6)*(n_1  + 2*n_2  + 2*n_3  + n_4) if D_1>D_eq else D_eq

                    else:

                        xyz_new = xyz_1
                        v_new   = np.zeros(self.m_InputData.m_SpaceDimension)
                        T_new   = T_p1
                        D_new   = D_1

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

            #------------------------------------------
            # set new xyz, v_p, T_p, and D_p for point
            #------------------------------------------
            self.m_LagrangianData.setX(p, xyz_new)
            self.m_LagrangianData.setVectorData(v_new, a_ArrayName='Velocity', a_DataID=p)
            self.m_LagrangianData.setScalarData(T_new, a_ArrayName='Temperature', a_DataID=p)
            self.m_LagrangianData.setScalarData(D_new*0.5, a_ArrayName='Radius', a_DataID=p)

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

            #-----------------------------------------------
            # update response times and evaporation factors
            #-----------------------------------------------
            self.m_MomentumResponseTime[p] = self.calculateMomentumResponseTime(D_new*0.5)
            self.m_ThermalResponseTime[p]  = self.calculateThermalResponseTime(D_new*0.5)
            self.m_EvaporationFactor[p]    = self.calculateEvaporationFactor(D_new*0.5)
            self.m_LagrangianData.setScalarData(self.m_MomentumResponseTime[p], a_ArrayName='Momentum Response Time', a_DataID=p)
            self.m_LagrangianData.setScalarData(self.m_ThermalResponseTime[p],  a_ArrayName='Thermal Response Time',  a_DataID=p)
            self.m_LagrangianData.setScalarData(self.m_EvaporationFactor[p],    a_ArrayName='Evaporation Factor',     a_DataID=p)
