#--------------------------------------------------------------------------------
# Module providing encapsulated input data handling capabilities to configure
# the particle dynamics simulations
#
# Author:       Debanjan Mukherjee
# Institution:  University of California, Berkeley
# Last Edited:  June 2018
#--------------------------------------------------------------------------------
import sys, os
import vtk
import numpy as np
from decimal import Decimal

#--------------------------------------------------------------------------------
# Problem input data structure, with all input data handled in a protected manner
#--------------------------------------------------------------------------------
class SimInputs:

    #--------------------------------------
    # initialize the class with member data
    #--------------------------------------
    def __init__(self, a_FileName):

        self.m_InputFile            = a_FileName
        temp                        = a_FileName.rfind('/')

        if temp != -1:
            self.m_RootPath   = a_FileName[0:temp]+'/'
        else:
            self.m_RootPath   = ''

        self.m_FlowGeometryPrimitives = {}
        self.m_LagrangianGeometryPrimitives = {}
        self.m_TimeSyncMapper       = {}
        self.m_LagrangianModuleParams = {}
        self.m_InitialVelocity = {}
        self.m_ResidenceModuleParams = {}
        self.m_ParticleInjectionParams = {}
        self.m_IDDVParams = {}

    #-----------------------------------------------------------------
    # read a formated ASCII input file and populate member data fields
    #-----------------------------------------------------------------
    def readInputFile(self):

        inputFileObj = open(self.m_InputFile)

        for line in inputFileObj:

            if not line.startswith('#'):

                lineDict = line.split('=')

                # checked
                if lineDict[0].strip() == 'Project directory':

                    tempRoot = lineDict[1].strip()
                    if tempRoot != 'default': self.m_RootPath = tempRoot

                # checked
                elif lineDict[0].strip() == 'Problem dimension':

                    self.m_SpaceDimension = int(lineDict[1].strip())

                # checked
                elif lineDict[0].strip() == 'Tracer input file':

                    self.m_TracerInput = self.m_RootPath + lineDict[1].strip()

                # checked
                elif lineDict[0].strip() == 'Discrete element input file':

                    tempList  = lineDict[1].strip().split()
                    self.m_DiscreteElementEnsembleFile  = self.m_RootPath + 'inputs/' + tempList[2].strip()
                    if tempList[6].strip() == 'none':
                        self.m_DiscreteElementTransformFile = None
                    else:
                        self.m_DiscreteElementTransformFile = self.m_RootPath + 'inputs/' + tempList[6].strip()

                # checked
                elif lineDict[0].strip() == 'Flow data directory':

                    self.m_FlowDirectory = self.m_RootPath + lineDict[1].strip()

                # checked
                elif lineDict[0].strip() == 'Flow data file tag':

                    self.m_FlowFileTag = lineDict[1].strip()

                #checked
                elif lineDict[0].strip() == 'SDF files':

                    tempList              = lineDict[1].strip().split()
                    self.m_Meshfile      = self.m_RootPath + tempList[2].strip()
                    self.m_Surfacefile = self.m_RootPath + tempList[6].strip()
                    self.m_OutputSDFGradfile   = self.m_RootPath + tempList[10].strip()

                #Pulling gravity vectort from input file
                elif lineDict[0].strip() == 'Gravity Vector':

                    tempList = lineDict[1].strip().split(';') #split into vectors

                    self.m_GravityVector = np.array([float(tempList[0].strip()), float(tempList[1].strip()), float(tempList[2].strip())])

                # checked
                elif lineDict[0].strip() == 'Restitution':

                    self.m_Restitution = float(lineDict[1].strip())

                # checked
                elif lineDict[0].strip() == 'Flow data field name':

                    self.m_FlowFieldName = lineDict[1].strip()

                # checked
                elif lineDict[0].strip() == 'Flow properties':

                    tempList              = lineDict[1].strip().split()
                    self.m_FlowModel      = int(tempList[2])
                    self.m_FluidViscosity = float(tempList[6])
                    self.m_FluidDensity   = float(tempList[10])

                # checked
                elif lineDict[0].strip() == 'Fixed mesh':

                    if lineDict[1].strip() == 'TRUE' or lineDict[1].strip() == 'True':
                        self.m_IsFixedMesh = True
                    else:
                        self.m_IsFixedMesh = False

                # checked
                elif lineDict[0].strip() == 'Particle in cell locator':

                    tempList = lineDict[1].strip().split()

                    if tempList[2].strip() == 'True':
                        self.m_LocatorType = 'oct'
                    elif tempList[6].strip() == 'True':
                        self.m_LocatorType = 'tre'
                    elif tempList[10].strip() == 'True':
                        self.m_LocatorType == 'bsp'

                # checked
                elif lineDict[0].strip() == 'Lagrangian tracer output file':

                    self.m_TracerOutputFile = lineDict[1].strip()

                # checked
                elif lineDict[0].strip() == 'Lagrangian field output file':

                    self.m_LagrangianOutputFile = lineDict[1].strip()

                elif lineDict[0].strip() == 'Write interval':

                    self.m_WriteInterval = int(lineDict[1].strip())

                # checked
                elif lineDict[0].strip() == 'Standard flow geometry attributes':

                    tempList = lineDict[1].strip().split()

                    self.m_FlowGeometryPrimitives['X0'] = None if tempList[2].strip()  == 'none' else float(tempList[2])
                    self.m_FlowGeometryPrimitives['Y0'] = None if tempList[4].strip()  == 'none' else float(tempList[4])
                    self.m_FlowGeometryPrimitives['Z0'] = None if tempList[6].strip()  == 'none' else float(tempList[6])
                    self.m_FlowGeometryPrimitives['DX'] = None if tempList[10].strip() == 'none' else float(tempList[10])
                    self.m_FlowGeometryPrimitives['DY'] = None if tempList[14].strip() == 'none' else float(tempList[14])
                    self.m_FlowGeometryPrimitives['DZ'] = None if tempList[18].strip() == 'none' else float(tempList[18])
                    # self.m_FlowGeometryPrimitives['R0'] = None if tempList[22].strip() == 'none' else float(tempList[22])

                # checked
                elif lineDict[0].strip() == 'Standard flow boundary tags':

                    tempList = lineDict[1].strip().split()

                    self.m_FlowGeometryPrimitives['B0'] = None if tempList[2].strip()  == 'none' else int(tempList[2].strip())
                    self.m_FlowGeometryPrimitives['B1'] = None if tempList[6].strip()  == 'none' else int(tempList[6].strip())
                    self.m_FlowGeometryPrimitives['B2'] = None if tempList[10].strip() == 'none' else int(tempList[10].strip())
                    self.m_FlowGeometryPrimitives['B3'] = None if tempList[14].strip() == 'none' else int(tempList[14].strip())
                    self.m_FlowGeometryPrimitives['B4'] = None if tempList[18].strip() == 'none' else int(tempList[18].strip())
                    self.m_FlowGeometryPrimitives['B5'] = None if tempList[22].strip() == 'none' else int(tempList[22].strip())

                # checked
                elif lineDict[0].strip() == 'Standard Lagrangian geometry attributes':

                    tempList = lineDict[1].strip().split()

                    self.m_LagrangianGeometryPrimitives['X0'] = None if tempList[2].strip()  == 'none' else float(tempList[2])
                    self.m_LagrangianGeometryPrimitives['Y0'] = None if tempList[4].strip()  == 'none' else float(tempList[4])
                    self.m_LagrangianGeometryPrimitives['Z0'] = None if tempList[6].strip()  == 'none' else float(tempList[6])
                    self.m_LagrangianGeometryPrimitives['X1'] = None if tempList[8].strip()  == 'none' else float(tempList[8])
                    self.m_LagrangianGeometryPrimitives['Y1'] = None if tempList[10].strip() == 'none' else float(tempList[10])
                    self.m_LagrangianGeometryPrimitives['Z1'] = None if tempList[12].strip() == 'none' else float(tempList[12])
                    self.m_LagrangianGeometryPrimitives['DX'] = None if tempList[16].strip() == 'none' else float(tempList[16])
                    self.m_LagrangianGeometryPrimitives['DY'] = None if tempList[20].strip() == 'none' else float(tempList[20])
                    self.m_LagrangianGeometryPrimitives['DZ'] = None if tempList[24].strip() == 'none' else float(tempList[24])

                # checked
                elif lineDict[0].strip() == 'Particle properties':

                    tempList                = lineDict[1].strip().split()
                    self.m_TracerDensity    = float(tempList[2])
                    self.m_TracerDiffCoeff  = float(tempList[6])

                # checked
                elif lineDict[0].strip() == 'Data file index':

                    tempList                = lineDict[1].strip().split()
                    self.m_DataIndexStart   = int(tempList[2])
                    self.m_DataIndexStop    = int(tempList[6])
                    self.m_DataIndexDelta   = int(tempList[10])
                    self.m_DataIndexCount   = int(tempList[14])

                # checked
                elif lineDict[0].strip() == 'Data file timing':

                    tempList                = lineDict[1].strip().split()
                    self.m_DataTimeStart    = float(tempList[2])
                    self.m_DataTimeStop     = float(tempList[6])
                    self.m_DataTimeDelta    = float(tempList[10])
                    self.m_DataPeriodic     = tempList[14].lower() == 'true'
                    self.m_dataDt_prec = len(tempList[10].split('.')[1]) #---Delta precision should be directly read from input file and split for calculation edited by Sreeparna---#

                # checked
                elif lineDict[0].strip() == 'Simulation timing':

                    tempList                = lineDict[1].strip().split()
                    self.m_SimTStart        = float(tempList[2])
                    self.m_SimTStop         = float(tempList[6])
                    self.m_Dt               = float(tempList[10])
                    self.m_Dt_prec = len(tempList[10].split('.')[1]) #---Delta precision should be directly read from input file and split for calculation edited by Sreeparna---#

                # checked
                elif lineDict[0].strip() == 'Integration setup':

                    tempList  = lineDict[1].strip().split()
                    self.m_IntegrationScheme  = tempList[2].strip()

                elif lineDict[0].strip() == 'Dump tracer output intervals':

                    self.m_DumpInterval     = int(lineDict[1].strip())

                # checked
                elif lineDict[0].strip() == 'Particle injection':

                    tempList = lineDict[1].strip().split()

                    self.m_IsInjectParticles = True if tempList[0] == 'True' else False
                    self.m_numberParticleInjections    = int(tempList[4])
                    self.m_startParticleInjections     = float(tempList[8])
                    self.m_intervalParticleInjections  = float(tempList[12])

                # checked
                elif lineDict[0].strip() == 'Particle pathlines':

                    tempList                    = lineDict[1].strip().split()
                    self.m_PathlineCompute      = True if (tempList[0] == 'True' or tempList[0] == 'TRUE') else False
                    self.m_PathSubsampleTime    = int(tempList[4])
                    self.m_PathSubsamplePoints  = int(tempList[6])
                    self.m_PathDDGCalculate     = True if (tempList[10] == 'True' or tempList[10] == 'TRUE') else False

                # checked
                elif lineDict[0].strip() == 'Residence time module setup':

                    tempList = lineDict[1].strip().split()

                    self.m_ResidenceModuleParams['Mode'] = tempList[2]

                    if len(tempList[6].split()) == 1:
                        self.m_ResidenceModuleParams['ROI'] = self.m_RootPath + 'inputs/' + tempList[6]
                    else:
                        self.m_ResidenceModuleParams['ROI'] = [float(x) for x in tempList[6].split()]

                    self.m_ResidenceModuleParams['InjFile'] = self.m_RootPath + 'inputs/' + tempList[10]
                    self.m_ResidenceModuleParams['InjPeriod'] = int(tempList[12])

                # checked
                elif lineDict[0].strip() == 'Lagrangian module setup':

                    tempList = lineDict[1].strip().split()
                    self.m_LagrangianModuleParams['isFWDFTLE'] = True if tempList[2].strip() == 'True' else False
                    self.m_LagrangianModuleParams['isBWDFTLE'] = True if tempList[6].strip() == 'True' else False
                    self.m_LagrangianModuleParams['isSretch'] = True if tempList[10].strip() == 'True' else False
                    self.m_LagrangianModuleParams['isStrain'] = True if tempList[14].strip() == 'True' else False

                # checked
                elif lineDict[0].strip() == 'Lagrangian kinematics unit vectors':

                    tempList = lineDict[1].strip().split()
                    self.m_LagrangianModuleParams['vector1'] = np.array([float(tempList[2].strip()), float(tempList[3].strip()), float(tempList[4].strip())])
                    self.m_LagrangianModuleParams['vector2'] = np.array([float(tempList[8].strip()), float(tempList[9].strip()), float(tempList[10].strip())])
                    self.m_LagrangianModuleParams['vector3'] = np.array([float(tempList[14].strip()), float(tempList[15].strip()), float(tempList[16].strip())])

                elif lineDict[0].strip() == 'Lagrangian initial velocity':

                    tempList = lineDict[1].strip().split(';')
                    self.m_InitialVelocity['u'] = tempList[0].strip()
                    self.m_InitialVelocity['v'] = tempList[1].strip()
                    self.m_InitialVelocity['w'] = tempList[2].strip()

                elif lineDict[0].strip() == 'Particle Size':

                    tempList = lineDict[1].strip().split(';')
                    self.m_isConstantSize = tempList[0].strip().split(':')[1].strip().lower() == 'true'
                    if self.m_isConstantSize:
                        self.m_TracerRadius   = float(tempList[1].strip().split(':')[1].strip())
                    else:
                        self.m_AverageRadius  = float(tempList[1].strip().split(':')[1].strip())
                        self.m_ShapeParam     = int(tempList[2].strip().split(':')[1].strip())

                elif lineDict[0].strip() == 'Gravitational Effects':

                    tempList = lineDict[1].strip().split(';')
                    self.m_GravEffects = tempList[0].strip().lower() == 'true'
                    if self.m_GravEffects:
                        self.m_GravDir  = tempList[1].strip().split(':')[1].strip().lower()
                        self.m_GravVal = float(tempList[2].strip().split(':')[1].strip())

                elif lineDict[0].strip() == 'Infectious Disease Model Additional Field Names':

                    tempList = lineDict[1].strip().split(';')
                    self.m_TemperatureFieldName = tempList[0].strip().split(':')[1].strip()
                    self.m_TempUnits            = tempList[1].strip().split(':')[1].strip()
                    self.m_PressureFieldName    = tempList[2].strip().split(':')[1].strip()
                    self.m_PressUnits           = tempList[3].strip().split(':')[1].strip()
                    self.m_AbsRelPress          = tempList[4].strip().split(':')[1].strip().lower()
                    self.m_PressConversion      = float(tempList[5].strip().split(':')[1].strip())

                elif lineDict[0].strip() == 'Infectious Disease Model Parameters':

                    tempList = lineDict[1].strip().split(';')
                    self.m_SpecificHeat   = float(tempList[0].strip().split(':')[1].strip())
                    self.m_InitialTemp    = float(tempList[1].strip().split(':')[1].strip())

                elif lineDict[0].strip() == 'Relative Humidity':

                    self.m_Humidity = float(lineDict[1].strip().split(':')[1].strip())

                elif lineDict[0].strip() == 'Fluid Domain Volume':

                    self.m_FluidVolume = float(lineDict[1].strip())

                elif lineDict[0].strip() == 'Resume simulation':

                    self.m_ResumeSimulation = lineDict[1].strip().lower() == 'true'

        inputFileObj.close()

    def getTracerInput(self):
        return self.m_TracerInput

    def getProblemDimension(self):
        return self.m_SpaceDimension

    def getDiscreteElementEnsemble(self):
        return self.m_DiscreteElementEnsembleFile

    def getDiscreteElementTransforms(self):
        return self.m_DiscreteElementTransformFile

    def getFlowDataFileName(self, a_ID=None, a_Legacy=False, a_Foam=False):

        if a_ID is not None:
            if a_Legacy == False:
                return self.m_FlowDirectory + self.m_FlowFileTag + str(a_ID) + ".vtu"
            else:
                return self.m_FlowDirectory + self.m_FlowFileTag + str(a_ID) + ".vtk"
        else:
            if a_Legacy == False and a_Foam == False:
                return self.m_FlowDirectory + self.m_FlowFileTag + ".vtu"
            elif a_Legacy == True and a_Foam == False:
                return self.m_FlowDirectory + self.m_FlowFileTag + ".vtk"
            elif a_Foam == True:
                return self.m_FlowDirectory + self.m_FlowFileTag + ".foam"

    #Added by SM and RR
    def getMeshdata(self):
        return self.m_Meshfile #sdf meshfile should be in vtu

    def getRootPath(self):
        return self.m_RootPath

    #Added by SM and RR
    def getSurfacedata(self):
        return self.m_Surfacefile #sdf surfacefile should be in stl

    #Added by SM and RR
    def getOutdata(self):
        return self.m_OutputSDFGradfile #sdf outfile should be in vtu

    #Added by SM and RR
    def getGravityVector(self):
        return self.m_GravityVector #units consistent with the velocity

    #Added by SM and RR
    def getRestitution(self):
        return self.m_Restitution # For restitution coefficient

    def getSimPrec(self):
        return self.m_Dt_prec # for simulation precision

## Added by AS
    def getFlowDataTimeStep(self):
        return self.m_DataTimeDelta

## Added by AS
    def getFlowDataStopFileIndex(self):
        return self.m_DataIndexStop

    def getSimulationStartTime(self):
        return self.m_SimTStart

    def getSimulationStopTime(self):
        return self.m_SimTStop

    def getIntegrationTimeStep(self):
        return self.m_Dt

    def getIntegrationScheme(self):
        return self.m_IntegrationScheme

    def getLocatorType(self):
        return self.m_LocatorType

    def getVelDataName(self):
        return self.m_FlowFieldName

    def getFluidDensity(self):
        return self.m_FluidDensity

    def getFluidViscosity(self):
        return self.m_FluidViscosity

    def getTracerDensity(self):
        return self.m_TracerDensity

    def getTracerDiffusivity(self):
        return self.m_TracerDiffCoeff

    def getTracerOutputFile(self, a_ID1=None, a_ID2=None):

      if (a_ID1 is not None) and (a_ID2 is not None):
          tempName    = self.m_TracerOutputFile.split('.')
          return self.m_RootPath + tempName[0] + '_Injection-' + str(a_ID1) + '_' + str(a_ID2)  + '.' + tempName[1]
      else:
          return self.m_RootPath + self.m_TracerOutputFile

    def getInitialVelocity(self):
        return self.m_InitialVelocity

    def getDumpInterval(self):
        return self.m_DumpInterval

    def getStandardDomainGeometryDefinition(self):
        return self.m_GeometryPrimitives

    def getPathlineTimeSubsampleInterval(self):
        return self.m_PathSubsampleTime

    def getPathlinePointSubsampleInterval(self):
        return self.m_PathSubsamplePoints

    def getPathlineFile(self, a_ID=None):
        if a_ID is None:
            tempName    = self.m_TracerOutputFile.split('.')
            return self.m_RootPath + tempName[0] + 'Path' + '.' + tempName[1]
        else:
            tempName    = self.m_TracerOutputFile.split('.')
            return self.m_RootPath + tempName[0] + 'Path' + str(a_ID) + '.' + tempName[1]

    def isDomainStandardGeometry(self):
        check = (self.m_FlowGeometryPrimitives['X0'] is not None) and \
                (self.m_FlowGeometryPrimitives['Y0'] is not None) and \
                (self.m_FlowGeometryPrimitives['Z0'] is not None) and \
                (self.m_FlowGeometryPrimitives['DX'] is not None) and \
                (self.m_FlowGeometryPrimitives['DY'] is not None) and \
                (self.m_FlowGeometryPrimitives['DZ'] is not None) and \
                (self.m_FlowGeometryPrimitives['R0'] is not None)

        return check

    def isDataLoopedPeriodic(self):
        return ( (self.m_SimTStop > self.m_DataTimeStop) and (self.m_DataPeriodic == True) )

    def isFixedMesh(self):
        return self.m_IsFixedMesh

    def isSteadyFlowData(self):
        return (self.m_DataIndexCount == 1)

    def isInjectParticles(self):
        return self.m_IsInjectParticles

    def isComputePathline(self):
        return self.m_PathlineCompute

    def isComputePathGeometry(self):
        return self.m_PathDDGCalculate

    #----------------------------------------------
    # BEGIN RESIDENCE TIME I/O FUNCTION DEFINITIONS
    #----------------------------------------------
    def isComputeResidenceTime(self):
        if self.m_ResidenceModuleParams['Mode'] == 'none':
            return False
        else:
            return True

    def isResidenceTimeComputeModeMapped(self):
        return self.m_ResidenceModuleParams['Mode'] == 'mapped'

    def isResidenceTimeComputeModeStreaming(self):
        return self.m_ResidenceModuleParams['Mode'] == 'stream'

    def isResidenceTimeROIFile(self):
        return isinstance(self.m_ResidenceModuleParams['ROI'], str)

    def isResidenceTimeROIBbox(self):
        return isinstance(self.m_ResidenceModuleParams['ROI'], list)

    def getResidenceTimeROIFileName(self):
        if isinstance(self.m_ResidenceModuleParams['ROI'], str):
            return self.m_ResidenceModuleParams['ROI']
        else:
            sys.exit('Error! Residence Time is configured to be a bounding box, no file name entered')

    def getResidenceTimeInjectionFile(self):
        return self.m_ResidenceModuleParams['InjFile']

    def getResidenceTimeInjectionInterval(self):
        return self.m_ResidenceModuleParams['InjPeriod']

    #----------------------------------------------
    # END RESIDENCE TIME I/O FUNCTION DEFINITIONS
    #----------------------------------------------

    #---------------------------------------------------------------------
    # For every tracer/particle trajectory integration timepoint, find the
    # corresponding time interval window for the data files
    #
    # NOTE:
    # - THIS IS AN EXTREMELY CRUCIAL FUNCTION/IMPLEMENTATION!!
    # - TAKE SPECIAL CARE WHILE MODIFYING THIS FOR FUTURE DEVELOPMENT!!
    #---------------------------------------------------------------------
    def getDataTimeWindows(self):
        print('data time windows version 1')
        if self.m_DataIndexCount == 1:

            t_Low     = np.array([0], dtype=np.float32)
            t_Up      = np.array([0], dtype=np.float32)
            ID_Low    = np.array([0], dtype=np.int32)
            ID_Up     = np.array([0], dtype=np.int32)
            timeIndex = np.array([0], dtype=np.int32)
            timeSync   = {'T_Low': t_Low, 'T_Up': t_Up, 'ID_Low': ID_Low, 'ID_Up': ID_Up, 'timeIndex': timeIndex}

        else:

            D0       = self.m_DataTimeStart
            D1       = self.m_DataTimeStop

            ID0      = self.m_DataIndexStart
            ID1      = self.m_DataIndexStop
            del_ID   = self.m_DataIndexDelta

            data_DT  = self.m_DataTimeDelta
            num_D    = self.m_DataIndexCount

            # simPrec  = len(str(self.m_Dt).split('.')[1])
            dataPrec = len(str(data_DT).split('.')[1])

            p_ID     = ID0 + np.arange(num_D)*del_ID     # array of data file indices
            p_DT     = D0 + np.arange(num_D)*data_DT     # array of data file times

            t0       = self.m_SimTStart
            t1       = self.m_SimTStop
            dT       = self.m_Dt

            num_T    = int((t1-t0)/dT)                   # number of integration time-instances
            p_simT   = t0 + np.arange(num_T+1)*dT        # array of all the integration times
            p_numT   = np.linspace(self.m_SimTStart, num_T, num_T+1).astype(int)

            t_Low    = np.floor(p_simT/data_DT)*data_DT  # the lower bound interval arrays
            t_Up     = np.ceil(p_simT/data_DT)*data_DT   # the upper bound interval arrays

            id_Low   = t_Low/data_DT                     # the lower bound data file index
            id_Up    = t_Up/data_DT                      # the upper bound data file index

            ## TEST IMPLEMENTATION FOR PERIODIC CASES
            if self.m_DataPeriodic == True:
                for i in range(id_Low.shape[0]):
                    if id_Low[i] >= num_D : id_Low[i] = np.mod(id_Low[i], num_D)

                for i in range(id_Up.shape[0]):
                    if id_Up[i] >= num_D : id_Up[i] = np.mod(id_Up[i], num_D)

            timeSync = {'T_Low': t_Low, 'T_Up': t_Up, 'ID_Low': p_ID[np.round(id_Low).astype(int)], 'ID_Up': p_ID[np.round(id_Up).astype(int)], 'timeIndex': p_numT}

        self.m_TimeSyncMapper = timeSync

        return timeSync

    #-----------------------------------------------------------------------
    # For every tracer/particle trajectory integration timepoint, find the
    # corresponding time interval window for the data files
    #
    # NOTE:
    # - THIS IS AN EXTREMELY CRUCIAL FUNCTION/IMPLEMENTATION!!
    # - TAKE SPECIAL CARE WHILE MODIFYING THIS FOR FUTURE DEVELOPMENT!!
    #
    # - This version was created by Joseph Wilson to correct issues stemming
    #   from the original version of the function. The original version
    #   is being investigated further and modified so it can be used
    #-----------------------------------------------------------------------
    def getDataTimeWindows_v2(self):
        print('data time windows version 2')
        if self.m_DataIndexCount == 1:

            t_Low     = np.array([0], dtype=np.float32)
            t_Up      = np.array([0], dtype=np.float32)
            ID_Low    = np.array([0], dtype=np.int32)
            ID_Up     = np.array([0], dtype=np.int32)
            timeIndex = np.array([0], dtype=np.int32)
            timeSync   = {'T_Low': t_Low, 'T_Up': t_Up, 'ID_Low': ID_Low, 'ID_Up': ID_Up, 'timeIndex': timeIndex}

        else:

            D0      = self.m_DataTimeStart
            D1      = self.m_DataTimeStop

            ID0     = self.m_DataIndexStart
            ID1     = self.m_DataIndexStop
            del_ID  = self.m_DataIndexDelta

            data_DT = self.m_DataTimeDelta
            num_D   = self.m_DataIndexCount

            simPrec  = self.m_Dt_prec # len(str(self.m_Dt).split('.')[1])
            dataPrec = self.m_dataDt_prec # len(str(data_DT).split('.')[1])
            # print('Dataprec = ', dataPrec)
            # sys.exit()

            num_T  = int((self.m_SimTStop-self.m_SimTStart)/self.m_Dt)+1
            p_numT = np.linspace(self.m_SimTStart, num_T, num_T+1).astype(int)

            id_Low = np.array([ID0])
            id_Up  = np.array([ID0])
            t_Low  = np.array([D0])
            t_Up   = np.array([D0])

            id_Low_track = ID0
            id_Up_track  = ID0+del_ID
            t_Low_track  = D0
            t_Up_track   = round(D0+data_DT, dataPrec)

            simTime   = round(self.m_SimTStart+self.m_Dt, simPrec)
            # print(simTime)
            # sys.exit()
            timeIndex = 1

            while simTime <= self.m_SimTStop:

                # print("Inside simtime loop", simTime, t_Low_track, t_Up_track)
                if simTime == t_Up_track:
                    t_Low_track   = np.round(t_Low_track+data_DT, dataPrec)
                    id_Low_track += del_ID
                elif simTime > t_Up_track:
                    t_Low_track   = np.round(t_Low_track+data_DT, dataPrec)
                    id_Low_track += del_ID
                    t_Up_track    = np.round(t_Low_track+data_DT, dataPrec)
                    id_Up_track  += del_ID

                id_Low = np.append(id_Low, id_Low_track)
                id_Up  = np.append(id_Up, id_Up_track)
                t_Low  = np.append(t_Low, t_Low_track)
                t_Up   = np.append(t_Up, t_Up_track)

                if simTime == t_Up_track:
                    t_Up_track   = np.round(t_Low_track+data_DT, dataPrec)
                    id_Up_track += del_ID

                simTime    = round(simTime+self.m_Dt, simPrec)
                timeIndex += 1

            ## TEST IMPLEMENTATION FOR PERIODIC CASES
            if self.m_DataPeriodic == True:
                for i in range(id_Low.shape[0]):
                    if id_Low[i] >= num_D : id_Low[i] = np.mod(id_Low[i], num_D)

                for i in range(id_Up.shape[0]):
                    if id_Up[i] >= num_D : id_Up[i] = np.mod(id_Up[i], num_D)

            timeSync = {'T_Low': t_Low, 'T_Up': t_Up, 'ID_Low': id_Low, 'ID_Up': id_Up, 'timeIndex': p_numT}

        self.m_TimeSyncMapper = timeSync

        return timeSync

    #-------------------------------------------------------------------------------
    # A utility function to print out the computed time windows from the input data.
    #
    # NOTE:
    # - may consider deprecating in future versions
    # - latest! - instead of deprecating have made timewindow an essential member
    #   data of the SimInputs class
    #-------------------------------------------------------------------------------
    def printDataTimeWindows(self):

        print("---------------------------------------------------------------------------------------")
        print("Printing all the time-window maps for loading data files in sync with integration steps")
        print("---------------------------------------------------------------------------------------")
        t_Low     = self.m_TimeSyncMapper['T_Low']
        t_Up      = self.m_TimeSyncMapper['T_Up']
        id_Low    = self.m_TimeSyncMapper['ID_Low']
        id_Up     = self.m_TimeSyncMapper['ID_Up']
        timeIndex = self.m_TimeSyncMapper['timeIndex']

        for i in range(len(t_Low)):
            print("T_Low: {0:10} T_Up: {1:10} ID_Low: {2:10} ID_Up: {3:10} timeIndex: {4:10}".format(t_Low[i], t_Up[i], id_Low[i], id_Up[i], timeIndex[i]))

        # ## TEST snippet
        # if self.m_TimeSyncMapper['ID_Low'].any() == self.m_TimeSyncMapper['ID_Up'].any():
        #     print(self.m_TimeSyncMapper['ID_Low'], self.m_TimeSyncMapper['ID_Up'])
        #     sys.exit()

        # import pandas as pd
        # # import openpyxl

        # # workbook = openpyxl.load_workbook("/Users/akshitasahni/Documents/PG/RMACC-job/dmgroup-scratch-dev_branch/Lagrangian-v2-Scratch/LVAD_collision_tests/versions/lagrangian-code/output.xlsx")
        # #
        # # for i in range(len(t_Low)):
        # #     df_ToExcel = pd.DataFrame({'T_Low': [t_Low[i]], 'T_Up': [t_Up[i]], 'ID_Low': [id_Low[i]], 'ID_Up': [id_Up[i]], 'timeindex': [timeIndex[i]]})
        # #     writer     = pd.ExcelWriter('/Users/akshitasahni/Documents/PG/RMACC-job/dmgroup-scratch-dev_branch/Lagrangian-v2-Scratch/LVAD_collision_tests/versions/lagrangian-code/output.xlsx', engine = 'openpyxl')
        # #     writer.book = workbook
        # #     df_ToExcel.to_excel(writer, sheet_name='Time window 1')
        # #     writer.save()
        # # writer.close()

        # out     = [[t_Low[i], t_Up[i], id_Low[i], id_Up[i], timeIndex[i]] for i in range(len(t_Low))]
        # dout    = pd.DataFrame(out)
        # outfile = '/Users/akshitasahni/Documents/PG/RMACC-job/dmgroup-scratch-dev_branch/Lagrangian-v2-Scratch/LVAD_collision_tests/versions/lagrangian-code/output.txt'
        # dout.to_csv(outfile, sep="\t", index=False)

        # print("Time sync data saved to excel file")
        # sys.exit()
        # ## TEST snippet
    #----------------------------------------------------------------------
    # a utility function that displays and prints all the configures input
    # variables, and allows teh user to verify the simulation configuration
    #----------------------------------------------------------------------
    def printAndVerify(self):

        print("The input file read for simulation was   : {:s}".format(self.m_InputFile))
        print("The number of spatial dimensions read    : {:d}".format(self.m_SpaceDimension))

        print("Tracer input file is located as follows  : {:s}".format(self.m_TracerInput))
        print("Tracer output is dumped into file        : {:s}".\
            format(self.m_RootPath+self.m_TracerOutputFile.split('.')[0]+'XX.'+self.m_TracerOutputFile.split('.')[1]))
        #print("Output files dumped once every N steps, N: {:d}\n".format(self.m_DumpInterval))


        print("Flow data file directory located at      : {:s}".format(self.m_FlowDirectory))
        print("Flow data file tag will be               : {:s}".format(self.m_FlowFileTag))
        print("Thus flow data files will be read as     : {:s}\n".format(self.m_FlowDirectory+self.m_FlowFileTag+"XX.vtk/vtu"))

        print("Flow data file index starts from         : {:d}".format(self.m_DataIndexStart))
        print("This is matched with starting time of    : {:f}\n".format(self.m_DataTimeStart))
        print("Flow data file index ends at             : {:d}".format(self.m_DataIndexStop))
        print("This is matched with end time of         : {:f}\n".format(self.m_DataTimeStop))
        print("Spacing between flow data file indices   : {:d}".format(self.m_DataIndexDelta))
        print("Spacing between flow data time instances : {:f}\n".format(self.m_DataTimeDelta))
        print("Number of flow data files saved          : {:d}\n".format(self.m_DataIndexCount))

        if self.m_DataIndexCount == 1:
            print("Looks like this configures a steady flow based simulation!\n")
        else:
            print("Looks like this configures an unsteady flow based simulation!\n")

        print("Discrete element ensemble file input     : {:s}".format(self.m_DiscreteElementEnsembleFile))
        print("Discrete element geometry transform file : {:s}\n".format(self.m_DiscreteElementTransformFile))

        if self.m_DiscreteElementEnsembleFile == 'none':
            print("Looks like simulation does not involve discrete elements!\n")

        print("Fluid viscosity                          : {:f}".format(self.m_FluidViscosity))
        print("Fluid density                            : {:f}".format(self.m_FluidDensity))
        print("Lagrangian tracer/particle density       : {:f}".format(self.m_TracerDensity))
        print("Lagrangian tracer/particle diffusion     : {:f}\n".format(self.m_TracerDiffCoeff))

        if self.isDomainStandardGeometry():

          print("Looks like simulation involves a standard computational domain geometry!\n")
          print("Standard domain anchor/center point      : {0:}, {1:}, {2:}".\
              format(self.m_FlowGeometryPrimitives['X0'], self.m_FlowGeometryPrimitives['Y0'], self.m_FlowGeometryPrimitives['Z0']))
          print("Standard domain spacing sizes            : {0:}, {1:}, {2:}".\
              format(self.m_FlowGeometryPrimitives['DX'], self.m_FlowGeometryPrimitives['DY'], self.m_FlowGeometryPrimitives['DZ']))
          print("Standard domain radius dimension         : {0:}".format(self.m_FlowGeometryPrimitives['R0']))

          print("Standard domain boundary condition tags  : {0:}, {1:}, {2:}, {3:}, {4:}, {5:}\n".\
              format( self.m_FlowGeometryPrimitives['B0'], self.m_FlowGeometryPrimitives['B1'],\
                      self.m_FlowGeometryPrimitives['B2'], self.m_FlowGeometryPrimitives['B3'],\
                      self.m_FlowGeometryPrimitives['B4'], self.m_FlowGeometryPrimitives['B5']))

        else:

            print("Looks like simulation involves a user-defined computational domain!\n")

        print("Numerical integration scheme to be used  : {:s}".format(self.m_IntegrationScheme))
        print("Integration start time                   : {:f}".format(self.m_SimTStart))
        print("Integration ends at                      : {:f}".format(self.m_SimTStop))
        print("Time-step for numerical integration      : {:f}".format(self.m_Dt))
        print("Choice for in-built VTK Locator          : {:s}\n".format(self.m_LocatorType))

        if self.m_LocatorType != 'TRE' or self.m_LocatorType != 'tre':
            print("Warning! Choice of vtkLocator is not the fastest!!\n")

        if self.isFixedMesh():
            print("Looks like this simulation is on a fixed-mesh!\n")
        else:
            print("Looks like this simulation is on a moving-mesh!\n")

        if self.m_PathlineCompute == True:
            print("Looks like this simulation computes particle pathlines!\n")

        print("Pathlines sub-sampled temporally every N : {:d}".format(self.m_PathSubsampleTime))
        print("Pathlines sub-sampled every P particles  : {:d}".format(self.m_PathSubsamplePoints))

        if self.m_PathDDGCalculate == True:
            print("Pathline geometry properties set to      : ON\n")
        else:
            print("Pathline geometry properties set to      : OFF\n")

        if self.m_ResidenceModuleParams:

            if self.isComputeResidenceTime() == True:
                print("Looks like this simulation computes particle residence time!\n")

                print("Residence time compute mode set to       : {:s}".format(self.m_ResidenceModuleParams['Mode']))

                if self.isResidenceTimeROIFile() == True:
                    print("Residence time R.O.I read from file      : {:s}".format(self.m_ResidenceModuleParams['ROI']))

                if self.isResidenceTimeROIBbox() == True:
                    print("Residence time R.O.I bbox X-extent       : {0:f}, {1:f}".format(self.m_ResidenceModuleParams['ROI'][0], self.m_ResidenceModuleParams['ROI'][1]))
                    print("Residence time R.O.I bbox Y-extent       : {0:f}, {1:f}".format(self.m_ResidenceModuleParams['ROI'][2], self.m_ResidenceModuleParams['ROI'][3]))
                    print("Residence time R.O.I bbox Z-extent       : {0:f}, {1:f}".format(self.m_ResidenceModuleParams['ROI'][4], self.m_ResidenceModuleParams['ROI'][5]))

                print("Tracer injection file for residence time : {:s}".format(self.m_ResidenceModuleParams['InjFile']))
                print("Injection frequency for residence time   : {:d}".format(self.m_ResidenceModuleParams['InjPeriod']))

#-------------------------------------------------------------------------------
# a simple read-print-verify mode is enabled with the input for a user/developer
# to test whether the input data structure is correctly configured
#-------------------------------------------------------------------------------
if __name__=="__main__":

    #
    # parse command-line argument
    #
    if len(sys.argv) != 2:
        sys.exit("Need Input Filename As An Argument")

    inputFile = sys.argv[1].strip()

    #
    # set-up problem input data
    #
    inputData = SimInputs(inputFile)
    inputData.readInputFile()
    inputData.getDataTimeWindows()

    #
    # print the input-data fields
    #
    os.system('clear')
    inputData.printAndVerify()
    inputData.printDataTimeWindows()
