#--------------------------------------------------------------------------------
# Module providing encapsulated input data handling capabilities to configure
# the particle dynamics simulations
#
# Authors:      Debanjan Mukherjee, Zachariah Irwin
# Institution:  University of Colorado, Boulder
# Last Edited:  February 2020
#--------------------------------------------------------------------------------
from __future__ import print_function
import sys, os
import vtk
import numpy as np

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

        self.m_FlowGeometryPrimitives          = {}
        self.m_LagrangianGeometryPrimitives    = {}
        self.m_TracerGeometryPrimitives        = {}
        self.m_TimeSyncMapper                  = {}
        self.m_LagrangianParams                = {}

    #-----------------------------------------------------------------
    # read a formated ASCII input file and populate member data fields
    #-----------------------------------------------------------------
    def readInputFile(self):

        inputFileObj = open(self.m_InputFile)

        for line in inputFileObj:

            if not line.startswith('#'):

                lineDict = line.split('=')

                if lineDict[0].strip() == 'Problem dimension':

                    self.m_SpaceDimension = int(lineDict[1].strip())

                elif lineDict[0].strip() == 'Project directory':

                    tempRoot = lineDict[1].strip()
                    if tempRoot != 'default': self.m_RootPath = tempRoot

                elif lineDict[0].strip() == 'Flow data directory':

                    self.m_FlowDirectory = self.m_RootPath + lineDict[1].strip()

                elif lineDict[0].strip() == 'Flow data file tag':

                    self.m_FlowFileTag = lineDict[1].strip()

                elif lineDict[0].strip() == 'Flow data field name':

                    self.m_FlowFieldName = lineDict[1].strip()

                elif lineDict[0].strip() == 'Tracer input file':

                    self.m_TracerInput = self.m_RootPath + lineDict[1].strip()

                elif lineDict[0].strip() == 'Discrete element input file':

                    tempList  = lineDict[1].strip().split()

                    self.m_DiscreteElementEnsembleFile  = self.m_RootPath + 'inputs/' + tempList[2].strip()

                    if tempList[6].strip() == 'none':
                        self.m_DiscreteElementTransformFile = None
                    else:
                        self.m_DiscreteElementTransformFile = self.m_RootPath + 'inputs/' + tempList[6].strip()

                elif lineDict[0].strip() == 'Lagrangian field output file':

                    self.m_LagrangianOutputFile = lineDict[1].strip()

                elif lineDict[0].strip() == 'Lagrangian tracer output file':

                    self.m_TracerOutputFile = lineDict[1].strip()

                elif lineDict[0].strip() == 'Data file index':

                    tempList                = lineDict[1].strip().split()
                    self.m_DataIndexStart   = int(tempList[2])
                    self.m_DataIndexStop    = int(tempList[6])
                    self.m_DataIndexDelta   = int(tempList[10])
                    self.m_DataIndexCount   = int(tempList[14])

                elif lineDict[0].strip() == 'Data file timing':

                    tempList                = lineDict[1].strip().split()
                    self.m_DataTimeStart    = float(tempList[2])
                    self.m_DataTimeStop     = float(tempList[6])
                    self.m_DataTimeDelta    = float(tempList[10])
                    self.m_DataPeriodic     = True if tempList[14] == 'True' or tempList[14] == 'TRUE' else False

                elif lineDict[0].strip() == 'Output file index':

                    tempList               = lineDict[1].strip().split()
                    self.m_OutIndexStart   = int(tempList[2])
                    self.m_OutIndexStop    = int(tempList[6])
                    self.m_OutIndexDelta   = int(tempList[10])

                elif lineDict[0].strip() == 'Output file timing':

                    tempList                = lineDict[1].strip().split()
                    self.m_OutTimeStart     = float(tempList[2])
                    self.m_OutTimeInterval  = float(tempList[6])

                elif lineDict[0].strip() == 'Simulation timing':

                    tempList                = lineDict[1].strip().split()
                    self.m_SimTStart        = float(tempList[2])
                    self.m_SimTStop         = float(tempList[6])
                    self.m_Dt               = float(tempList[10])

                elif lineDict[0].strip() == 'Integration setup':

                    tempList  = lineDict[1].strip().split()
                    self.m_IntegrationScheme  = tempList[2].strip()

                elif lineDict[0].strip() == 'Fixed mesh':

                    if lineDict[1].strip() == 'TRUE' or lineDict[1].strip() == 'True':
                        self.m_IsFixedMesh = True
                    else:
                        self.m_IsFixedMesh = False

                elif lineDict[0].strip() == 'Flow model':

                    self.m_FlowModel = int(lineDict[1].strip())

                elif lineDict[0].strip() == 'Lagrangian field calculations':

                    tempList = lineDict[1].strip().split()

                    if tempList[2].strip() == 'True':
                        self.m_LagrangianParams['FTLE'] = 1
                    elif tempList[2].strip() == 'False':
                        self.m_LagrangianParams['FTLE'] = 0
                    if tempList[6].strip() == 'True':
                        self.m_LagrangianParams['Stretch'] = 1
                    elif tempList[6].strip() == 'False':
                        self.m_LagrangianParams['Stretch'] = 0
                    if tempList[10].strip() == 'True':
                        self.m_LagrangianParams['Strain'] = 1
                    elif tempList[10].strip() == 'False':
                        self.m_LagrangianParams['Strain'] = 0

                elif lineDict[0].strip() == 'Unit vector 1':

                    tempList = lineDict[1].strip().split()

                    self.m_UnitVec1 = np.zeros((3), dtype=np.float64, order='F')
                    self.m_UnitVec1[0] = float(tempList[0])
                    self.m_UnitVec1[1] = float(tempList[2])
                    self.m_UnitVec1[2] = float(tempList[4])

                elif lineDict[0].strip() == 'Unit vector 2':

                    tempList = lineDict[1].strip().split()

                    self.m_UnitVec2 = np.zeros((3), dtype=np.float64, order='F')
                    self.m_UnitVec2[0] = float(tempList[0])
                    self.m_UnitVec2[1] = float(tempList[2])
                    self.m_UnitVec2[2] = float(tempList[4])

                elif lineDict[0].strip() == 'Unit vector 3':

                    tempList = lineDict[1].strip().split()

                    self.m_UnitVec3 = np.zeros((3), dtype=np.float64, order='F')
                    self.m_UnitVec3[0] = float(tempList[0])
                    self.m_UnitVec3[1] = float(tempList[2])
                    self.m_UnitVec3[2] = float(tempList[4])

                elif lineDict[0].strip() == 'Standard flow domain geometry attributes':

                    tempList = lineDict[1].strip().split()

                    self.m_FlowGeometryPrimitives['X0'] = None if tempList[2].strip() == 'none' else float(tempList[2])
                    self.m_FlowGeometryPrimitives['Y0'] = None if tempList[4].strip() == 'none' else float(tempList[4])
                    self.m_FlowGeometryPrimitives['Z0'] = None if tempList[6].strip() == 'none' else float(tempList[6])
                    self.m_FlowGeometryPrimitives['X1'] = None if tempList[8].strip() == 'none' else float(tempList[8])
                    self.m_FlowGeometryPrimitives['Y1'] = None if tempList[10].strip() == 'none' else float(tempList[10])
                    self.m_FlowGeometryPrimitives['Z1'] = None if tempList[12].strip() == 'none' else float(tempList[12])

                elif lineDict[0].strip() == 'Standard flow domain boundary condition tags':

                    tempList = lineDict[1].strip().split()

                    self.m_FlowGeometryPrimitives['B0'] = None if tempList[2].strip() == 'none' else int(tempList[2].strip())
                    self.m_FlowGeometryPrimitives['B1'] = None if tempList[6].strip() == 'none' else int(tempList[6].strip())
                    self.m_FlowGeometryPrimitives['B2'] = None if tempList[10].strip() == 'none' else int(tempList[10].strip())
                    self.m_FlowGeometryPrimitives['B3'] = None if tempList[14].strip() == 'none' else int(tempList[14].strip())
                    self.m_FlowGeometryPrimitives['B4'] = None if tempList[18].strip() == 'none' else int(tempList[18].strip())
                    self.m_FlowGeometryPrimitives['B5'] = None if tempList[22].strip() == 'none' else int(tempList[22].strip())

                elif lineDict[0].strip() == 'Standard Lagrangian field domain geometry attributes':

                    tempList = lineDict[1].strip().split()

                    self.m_LagrangianGeometryPrimitives['X0'] = None if tempList[2].strip() == 'none' else float(tempList[2])
                    self.m_LagrangianGeometryPrimitives['Y0'] = None if tempList[4].strip() == 'none' else float(tempList[4])
                    self.m_LagrangianGeometryPrimitives['Z0'] = None if tempList[6].strip() == 'none' else float(tempList[6])
                    self.m_LagrangianGeometryPrimitives['X1'] = None if tempList[8].strip() == 'none' else float(tempList[8])
                    self.m_LagrangianGeometryPrimitives['Y1'] = None if tempList[10].strip() == 'none' else float(tempList[10])
                    self.m_LagrangianGeometryPrimitives['Z1'] = None if tempList[12].strip() == 'none' else float(tempList[12])
                    self.m_LagrangianGeometryPrimitives['DX'] = None if tempList[16].strip() == 'none' else float(tempList[16])
                    self.m_LagrangianGeometryPrimitives['DY'] = None if tempList[20].strip() == 'none' else float(tempList[20])
                    self.m_LagrangianGeometryPrimitives['DZ'] = None if tempList[24].strip() == 'none' else float(tempList[24])

                elif lineDict[0].strip() == 'Particle injection domain geometry attributes':

                    tempList = lineDict[1].strip().split()

                    self.m_TracerGeometryPrimitives['X0'] = None if tempList[2].strip() == 'none' else float(tempList[2])
                    self.m_TracerGeometryPrimitives['Y0'] = None if tempList[4].strip() == 'none' else float(tempList[4])
                    self.m_TracerGeometryPrimitives['Z0'] = None if tempList[6].strip() == 'none' else float(tempList[6])
                    self.m_TracerGeometryPrimitives['X1'] = None if tempList[8].strip() == 'none' else float(tempList[8])
                    self.m_TracerGeometryPrimitives['Y1'] = None if tempList[10].strip() == 'none' else float(tempList[10])
                    self.m_TracerGeometryPrimitives['Z1'] = None if tempList[12].strip() == 'none' else float(tempList[12])
                    self.m_TracerGeometryPrimitives['DX'] = None if tempList[16].strip() == 'none' else float(tempList[16])
                    self.m_TracerGeometryPrimitives['DY'] = None if tempList[20].strip() == 'none' else float(tempList[20])
                    self.m_TracerGeometryPrimitives['DZ'] = None if tempList[24].strip() == 'none' else float(tempList[24])

                elif lineDict[0].strip() == 'Particle injection':

                    tempList = lineDict[1].strip().split()

                    if tempList[0] == 'False':
                        self.m_IsInjectParticles    = False
                    elif tempList[0] == 'True':
                        self.m_IsInjectParticles    = True

                    self.m_numberParticleInjections    = int(tempList[4])
                    self.m_startParticleInjections     = float(tempList[8])
                    self.m_intervalParticleInjections  = float(tempList[12])

                elif lineDict[0].strip() == 'Flow properties':

                    tempList              = lineDict[1].strip().split()
                    self.m_FluidDensity   = float(tempList[2])
                    self.m_FluidViscosity = float(tempList[6])

                elif lineDict[0].strip() == 'Particle properties':

                    tempList                = lineDict[1].strip().split()
                    self.m_TracerDensity    = float(tempList[2])
                    self.m_TracerDiffCoeff  = float(tempList[6])

                elif lineDict[0].strip() == 'Surface Normals':

                  tempList = lineDict[1].strip().split()

                  if tempList[2].strip() == 'True':
                    self.m_ConvertNormals = True
                  elif tempList[2].strip() == 'False':
                    self.m_ConvertNormals = False
                  if tempList[6].strip() == 'True':
                    self.m_FlipNormals = True
                  elif tempList[6].strip() == 'False':
                    self.m_FlipNormals = False
                  self.m_NormalScaling = float(tempList[10].strip())

                elif lineDict[0].strip() == 'Particle in cell locator':

                    tempList = lineDict[1].strip().split()

                    if tempList[2].strip() == 'True':
                        self.m_LocatorType = 'oct'
                    elif tempList[6].strip() == 'True':
                        self.m_LocatorType = 'tre'
                    elif tempList[10].strip() == 'True':
                        self.m_LocatorType == 'bsp'

                elif lineDict[0].strip() == 'Particle pathlines':

                    tempList                    = lineDict[1].strip().split()
                    self.m_PathlineCompute      = True if (tempList[0] == 'True' or tempList[0] == 'TRUE') else False
                    self.m_PathSubsampleTime    = int(tempList[4])
                    self.m_PathSubsamplePoints  = int(tempList[6])
                    self.m_PathDDGCalculate     = True if (tempList[10] == 'True' or tempList[10] == 'TRUE') else False

                elif lineDict[0].strip() == 'Residence time module setup':

                    tempList = lineDict[1].strip().split()

                    self.m_ResidenceTimeMode      = tempList[2]

                    if len(tempList[6].split()) == 1:
                        self.m_ResidenceTimeROI = self.m_RootPath + 'inputs/' + tempList[6]
                    else:
                        self.m_ResidenceTimeROI = [float(x) for x in tempList[6].split()]

                    self.m_ResidenceTimeInjFile   = self.m_RootPath + 'inputs/' + tempList[10]
                    self.m_ResidenceTimeInjPeriod = int(tempList[12])

        inputFileObj.close()

    def getTracerInput(self): return self.m_TracerInput

    def getProblemDimension(self): return self.m_SpaceDimension

    def getDiscreteElementEnsemble(self): return self.m_DiscreteElementEnsembleFile

    def getDiscreteElementTransforms(self): return self.m_DiscreteElementTransformFile

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

    def getSimulationStartTime(self): return self.m_SimTStart

    def getSimulationStopTime(self): return self.m_SimTStop

    def getIntegrationTimeStep(self): return self.m_Dt

    def getIntegrationScheme(self): return self.m_IntegrationScheme

    def getDataTimeDelta(self): return self.m_DataTimeDelta

    def getLocatorType(self): return self.m_LocatorType

    def getFlowModel(self): return self.m_FlowModel

    def getLagrangianParams(self): return self.m_LagrangianParams

    def getVelDataName(self): return self.m_FlowFieldName

    def getFluidDensity(self): return self.m_FluidDensity

    def getFluidViscosity(self): return self.m_FluidViscosity

    def getTracerDensity(self): return self.m_TracerDensity

    def getTracerDiffusivity(self): return self.m_TracerDiffCoeff

    def getLagrangianOutputFile(self, a_ID=None):

      if a_ID is not None:
          tempName    = self.m_LagrangianOutputFile.split('.')
          return self.m_RootPath + tempName[0] + str(a_ID) + '.' + tempName[1]
      else:
          return self.m_RootPath + self.m_LagrangianOutputFile

    def getTracerOutputFile(self, a_ID1=None, a_ID2=None):

      if (a_ID1 is not None) and (a_ID2 is not None):
          tempName    = self.m_TracerOutputFile.split('.')
          return self.m_RootPath + tempName[0] + 'Injection-' + str(a_ID1) + '_' + str(a_ID2)  + '.' + tempName[1]
      else:
          return self.m_RootPath + self.m_TracerOutputFile

    def getInjectionStartTime(self): return self.m_startParticleInjections

    def getInjectionInterval(self): return self.m_intervalParticleInjections

    def getNumberOfInjections(self): return self.m_numberParticleInjections

    def getStandardFlowDomainGeometryDefinition(self): return self.m_FlowGeometryPrimitives

    def getStandardLagrangianDomainGeometryDefinition(self): return self.m_LagrangianGeometryPrimitives

    def getStandardTracerDomainGeometryDefinition(self): return self.m_TracerGeometryPrimitives

    def getPathlineTimeSubsampleInterval(self): return self.m_PathSubsampleTime

    def getPathlinePointSubsampleInterval(self): return self.m_PathSubsamplePoints

    def getPathlineFile(self, a_ID=None):

        if a_ID is None:
            tempName    = self.m_TracerOutputFile.split('.')
            return self.m_RootPath + tempName[0] + 'Path' + '.' + tempName[1]
        else:
            tempName    = self.m_TracerOutputFile.split('.')
            return self.m_RootPath + tempName[0] + 'Path' + str(a_ID) + '.' + tempName[1]

    def getResidenceTimeROIFileName(self):

        if isinstance(self.m_ResidenceTimeROI, str):
            return self.m_ResidenceTimeROI
        else:
            sys.exit('Error! Residence Time is configured to be a bounding box, no file name entered')

    def getResidenceTimeInjectionFile(self): return self.m_ResidenceTimeInjFile

    def getResidenceTimeInjectionInterval(self): return self.m_ResidenceTimeInjPeriod

    def isDomainStandardGeometry(self):

        check = (self.m_FlowGeometryPrimitives['X0'] is not None) and \
                (self.m_FlowGeometryPrimitives['Y0'] is not None) and \
                (self.m_FlowGeometryPrimitives['Z0'] is not None) and \
                (self.m_LagrangianGeometryPrimitives['DX'] is not None) and \
                (self.m_LagrangianGeometryPrimitives['DY'] is not None) and \
                (self.m_LagrangianGeometryPrimitives['DZ'] is not None) and \
                (self.m_LagrangianGeometryPrimitives['R0'] is not None)

        return check

    def isDataLoopedPeriodic(self): return ( (self.m_SimTStop > self.m_DataTimeStop) and (self.m_DataPeriodic == True) )

    def isFixedMesh(self): return self.m_IsFixedMesh

    def isSteadyFlowData(self): return (self.m_DataIndexCount == 1)

    def isWriteTracers(self): return self.m_writeTracers

    def isComputeResidenceTime(self):
        if self.m_ResidenceTimeMode == 'none':
            return False
        else:
            return True

    def isComputePathline(self): return self.m_PathlineCompute

    def isComputePathGeometry(self): return self.m_PathDDGCalculate

    def isResidenceTimeComputeModeMapped(self): return self.m_ResidenceTimeMode == 'mapped'

    def isResidenceTimeComputeModeStreaming(self): return self.m_ResidenceTimeMode == 'stream'

    def isResidenceTimeROIFile(self): return isinstance(self.m_ResidenceTimeROI, str)

    def isResidenceTimeROIBbox(self): return isinstance(self.m_ResidenceTimeROI, list)

    #---------------------------------------------------------------------
    # For every tracer/particle trajectory integration timepoint, find the
    # corresponding time interval window for the data files
    #
    # NOTE:
    # - THIS IS AN EXTREMELY CRUCIAL FUNCTION/IMPLEMENTATION!!
    # - TAKE SPECIAL CARE WHILE MODIFYING THIS FOR FUTURE DEVELOPMENT!!
    #---------------------------------------------------------------------
    def getDataTimeWindows(self):

        if self.m_DataIndexCount == 1:

            t_Low   = np.array([0], dtype=np.float32)
            t_Up    = np.array([0], dtype=np.float32)
            ID_Low  = np.array([0], dtype=np.int32)
            ID_Up   = np.array([0], dtype=np.int32)
            timeSync = {'T_Low': t_Low, 'T_Up': t_Up, 'ID_Low': ID_Low, 'ID_Up': ID_Up}

        else:

            D0      = self.m_DataTimeStart
            D1      = self.m_DataTimeStop

            ID0     = self.m_DataIndexStart
            ID1     = self.m_DataIndexStop
            del_ID  = self.m_DataIndexDelta

            data_DT = self.m_DataTimeDelta
            num_D   = self.m_DataIndexCount

            p_ID    = ID0 + np.arange(num_D)*del_ID     # array of data file indices
            p_DT    = D0 + np.arange(num_D)*data_DT     # array of data file times

            t0      = self.m_SimTStart
            t1      = self.m_SimTStop
            dT      = self.m_Dt

            num_T   = int((t1-t0)/dT)                   # number of integration time-instances
            p_simT  = t0 + np.arange(num_T+1)*dT        # array of all the integration times

            t_Low   = np.floor(p_simT/data_DT)*data_DT           # the lower bound interval arrays
            t_Mid   = np.ceil(p_simT/data_DT)*data_DT            # the middle bound interval arrays
            t_Up    = np.ceil((p_simT + data_DT)/data_DT)*data_DT # the upper bound interval arrays

            id_Low   = t_Low/data_DT                     # the lower bound data file index
            id_Mid   = t_Mid/data_DT                     # the middle bound data file index
            id_Up    = t_Up/data_DT                      # the upper bound data file index

            ## TEST IMPLEMENTATION FOR PERIODIC CASES
            if self.m_DataPeriodic == True:
                for i in range(id_Low.shape[0]):
                    if id_Low[i] >= num_D : id_Low[i] = np.mod(id_Low[i], num_D)

                for i in range(id_Mid.shape[0]):
                    if id_Mid[i] >= num_D : id_Mid[i] = np.mod(id_Mid[i], num_D)

                for i in range(id_Up.shape[0]):
                    if id_Up[i] >= num_D : id_Up[i] = np.mod(id_Up[i], num_D)

            timeSync = {'T_Low': t_Low, 'T_Mid': t_Mid, 'T_Up': t_Up, 'ID_Low': p_ID[id_Low.astype(int)], 'ID_Mid' : p_ID[id_Mid.astype(int)], 'ID_Up': p_ID[id_Up.astype(int)]}

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

        #
        # SECTIONS OF OLD CODE. DEPRECATE IN FUTURE VERSIONS
        #
        # D0      = self.m_DataTimeStart
        # D1      = self.m_DataTimeStop
        # ID0     = self.m_DataIndexStart
        # ID1     = self.m_DataIndexStop
        # del_ID  = self.m_DataIndexDelta
        # data_DT = self.m_DataTimeDelta
        # num_D   = self.m_DataIndexCount
        # p_ID    = ID0 + np.arange(num_D)*del_ID
        # p_DT    = D0 + np.arange(num_D)*data_DT
        # t0      = self.m_SimTStart
        # t1      = self.m_SimTStop
        # dT      = self.m_Dt
        # num_T   = int((t1-t0)/dT)
        # p_simT  = t0 + np.arange(num_T+1)*dT
        # t_Low   = np.floor(p_simT/data_DT)*data_DT
        # t_Up    = np.ceil(p_simT/data_DT)*data_DT
        # id_Low  = t_Low/data_DT
        # id_Up   = t_Up/data_DT
        # id_Low  = id_Low.astype(int)
        # id_Up   = id_Up.astype(int)
        #
        # SECTIONS OF OLD CODE. DEPRECATE IN FUTURE VERSIONS
        #
        # print "Number of data files         : {:d}".format(self.m_DataIndexCount)
        # print "Number of integration times  : {:d}".format(num_T)
        # print "Length of data indices       : {:d}".format(len(p_ID))
        # print "Number of lower bounds       : {:d}".format(len(t_Low))
        # print "Number of upper bounds       : {:d}".format(len(t_Up))
        # print "Number of lower bound ids    : {:d}".format(len(id_Low))
        # print "Number of upper bound ids    : {:d}".format(len(id_Up))
        # print "Maximum file id in list      : {:d}".format(np.max(id_Low))
        # print "Minimum file id in list      : {:d}".format(np.max(id_Up))

        print("---------------------------------------------------------------------------------------")
        print("Printing all the time-window maps for loading data files in sync with integration steps")
        print("---------------------------------------------------------------------------------------")
        t_Low   = self.m_TimeSyncMapper['T_Low']
        t_Mid   = self.m_TimeSyncMapper['T_Mid']
        t_Up    = self.m_TimeSyncMapper['T_Up']
        id_Low  = self.m_TimeSyncMapper['ID_Low']
        id_Mid  = self.m_TimeSyncMapper['ID_Mid']
        id_Up   = self.m_TimeSyncMapper['ID_Up']

        for i in range(len(t_Low)):
            print("T_Low: {0:10} T_Mid: {1:10} T_Up: {2:10} ID_Low: {3:10} ID_Mid: {4:10} ID_Up: {5:10}".format(t_Low[i], t_Mid[i], t_Up[i], id_Low[i], id_Mid[i], id_Up[i]))

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
        print("Output files dumped once every N steps, N: {:d}\n".format(self.m_DumpInterval))


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

        print("Discrete element ensemble file input     : {!s}".format(self.m_DiscreteElementEnsembleFile))
        print("Discrete element geometry transform file : {!s}\n".format(self.m_DiscreteElementTransformFile))

        if self.m_DiscreteElementEnsembleFile == 'none': print("Looks like simulation does not invlve discrete elements!\n")

        print("Fluid viscosity                          : {:f}".format(self.m_FluidViscosity))
        print("Fluid density                            : {:f}".format(self.m_FluidDensity))
        print("Lagrangian tracer/particle density       : {:f}".format(self.m_TracerDensity))
        print("Lagrangian tracer/particle diffusion     : {:f}\n".format(self.m_TracerDiffCoeff))

        if self.isDomainStandardGeometry():

          print("Looks like simulation involves a standard computational domain geometry!\n")
          print("Standard domain anchor/center point      : {0:}, {1:}, {2:}".\
              format(self.m_FlowGeometryPrimitives['X0'], self.m_FlowGeometryPrimitives['Y0'], self.m_FlowGeometryPrimitives['Z0']))
          print("Standard domain end point                : {0:}, {1:}, {2:}".\
              format(self.m_FlowGeometryPrimitives['X1'], self.m_FlowGeometryPrimitives['Y1'], self.m_FlowGeometryPrimitives['Z1']))
          print("Standard domain spacing sizes            : {0:}, {1:}, {2:}".\
              format(self.m_LagrangianGeometryPrimitives['DX'], self.m_LagrangianGeometryPrimitives['DY'], self.m_LagrangianGeometryPrimitives['DZ']))
          print("Standard domain radius dimension         : {0:}".format(self.m_LagrangianGeometryPrimitives['R0']))

          print("Standard domain boundary condition tags  : {0:}, {1:}, {2:}, {3:}, {4:}, {5:}\n".\
              format( self.m_FlowGeometryPrimitives['B0'], self.m_FlowGeometryPrimitives['B1'],\
                      self.m_FlowGeometryPrimitives['B2'], self.m_FlowGeometryPrimitives['B3'],\
                      self.m_FlowGeometryPrimitives['B4'], self.m_FlowGeometryPrimitives['B5']))

        else:

            print("Looks like simulation involves a user-defined computational domain!\n")

        print("Numerical integration scheme to be used  : {!s}".format(self.m_IntegrationScheme))
        print("Integration start time                   : {:f}".format(self.m_SimTStart))
        print("Integration ends at                      : {:f}".format(self.m_SimTStop))
        print("Time-step for numerical integration      : {:f}".format(self.m_Dt))
        print("Choice for in-built VTK Locator          : {!s}\n".format(self.m_LocatorType))

        if self.m_LocatorType != 'TRE' or self.m_LocatorType != 'tre':
            print("Warning! Choice of vtkLocator is not the fastest!!\n")

        if self.isFixedMesh():
            print("Looks like this simulation is on a fixed-mesh!\n")
        else:
            print("Looks like this simulation is on a moving-mesh!\n")

        if self.m_PathlineCompute == True:
            print("Looks like this simulation computes particle pathlines!\n")

        print("Pathlines sub-sampled temporally every N : {:f}".format(self.m_PathSubsampleTime))
        print("Pathlines sub-sampled every P particles  : {:f}".format(self.m_PathSubsamplePoints))
        if self.m_PathDDGCalculate == True:
            print("Pathline geometry properties set to      : ON\n")
        else:
            print("Pathline geometry properties set to      : OFF\n")

        if self.isComputeResidenceTime() == True:
            print("Looks like this simulation computes particle residence time!\n")

        print("Residence time compute mode set to       : {!s}".format(self.m_ResidenceTimeMode))
        if self.isResidenceTimeROIFile() == True:
            print("Residence time R.O.I read from file      : {!s}".format(self.m_ResidenceTimeROI))

        if self.isResidenceTimeROIBbox() == True:
            print("Residence time R.O.I bbox X-extent       : {0:f}, {1:f}".format(self.m_ResidenceTimeROI[0], self.m_ResidenceTimeROI[1]))
            print("Residence time R.O.I bbox Y-extent       : {0:f}, {1:f}".format(self.m_ResidenceTimeROI[2], self.m_ResidenceTimeROI[3]))
            print("Residence time R.O.I bbox Z-extent       : {0:f}, {1:f}".format(self.m_ResidenceTimeROI[4], self.m_ResidenceTimeROI[5]))

        print("Tracer injection file for residence time : {:s}".format(self.m_ResidenceTimeInjFile))
        print("Injection frequency for residence time   : {:d}".format(self.m_ResidenceTimeInjPeriod))


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
    #inputData.printAndVerify()
    inputData.printDataTimeWindows()
