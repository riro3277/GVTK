import sys, vtk, os
import numpy as np
from scipy.spatial import ConvexHull



def arrayParser(str_arr, dimension, valType):
    '''
    parse input arrays when dimension and data type are specified

    inputs:
        str_arr (str): string array of the values
        dimension (int): dimension of the array, either 1 or 2
        valType: output data type
    outputs:
        parsed (nump array): parsed array
    '''
    lst_arr = str_arr.split(';')
    lst_arr = [x.strip(' ') for x in lst_arr]
    num_elems = len(lst_arr)
    num_arr   = np.zeros((num_elems, 3))
    for i in range(num_elems):
        vals = lst_arr[i][1:-1].split(',')
        vals = [x.strip(' ') for x in vals]
        num_arr[i,:] = np.fromiter(vals, dtype=valType)
    return num_arr



# def exposureRiskSphere(point, center, radii, num_exp_zones, viability):
#     '''
#     calculate exposure risk weight for a spherical exposure zone
#
#     inputs:
#         point (float array): cartesian coordinates of particle
#         center (float array): cartesian coordinates of center of head
#         radii (float array): array of exposure zone radii
#         num_exp_zones (int): number of exposure zones
#         viability (float): viability index
#     '''
#     #---------------------------------------------------------------
#     # magnitude of particle position vector relative to head center
#     #---------------------------------------------------------------
#     vel = point-center
#     r_p = np.sqrt(dot(vel,vel))
#
#     #----------------------------------------
#     # check if particle is in exposure zones
#     #----------------------------------------
#     zone_check = r_p <= radii[num_exp_zones-1]
#
#     #-------------------------
#     # calculate exposure risk
#     #-------------------------
#     if zone_check:
#         weights = np.linspace(100,0,num_exp_zones+1)
#         for i in range(num_exp_zones):
#             if i == 0:
#                 if r_p <= radii[0]:
#                     return weights[0]*viability
#             else:
#                 if r_p <= radii[i] and r_p > radii[i-1]:
#                     return weights[i]*viability
#     else:
#         return 0
#
#
#
# def exposureRiskCone(point, center, radii, zone_angle, num_exp_zones, viability, orientation):
#     '''
#     calculate exposure risk weight for a conical exposure zone
#
#     inputs:
#         point (float array): cartesian coordinates of particle
#         center (float array): cartesian coordinates of center of head
#         radii (float array): array of exposure zone radii
#         num_exp_zones (int): number of exposure zones
#         viability (float): viability index
#         orientation (float): cartesian vector of subject orientation
#         zone_angle (float): radian angle between orientation and side of exposure zones
#     '''
#     #---------------------------------------------------------------
#     # magnitude of particle position vector relative to head center
#     #---------------------------------------------------------------
#     vel = point-center
#     r_p = np.sqrt(dot(vel,vel))
#     zone_check = r_p <= radii[num_exp_zones-1]
#
#     #-----------------------------------
#     # check angle if zone_check is true
#     #-----------------------------------
#     if zone_check:
#         mag_p     = np.sqrt(np.dot(point,point))
#         mag_c     = np.sqrt(np.dot(center,center))
#         rel_angle = np.arccos(np.sqrt(np.dot(point,center)) / (mag_p * mag_c))
#
#         #-------------------------------------------------------
#         # if angle is within range then calculate exposure risk
#         #-------------------------------------------------------
#         if rel_angle <= zone_angle:
#             weights = np.linspace(100,0,num_exp_zones+1)
#             for i in range(num_exp_zones):
#                 if i == 0:
#                     if r_p <= radii[0]:
#                         return weights[0]*viability
#                 else:
#                     if r_p <= radii[i] and r_p > radii[i-1]:
#                         return weights[i]*viability
#     else:
#         return 0
#
#
#
# def exposureRiskCyl(point, center, radii, zone_angle, zone_height, num_exp_zones, viability, orientation):
#     '''
#     calculate exposure risk weight for a cylindrical section exposure zone
#
#     inputs:
#         point (float array): cartesian coordinates of particle
#         center (float array): cartesian coordinates of center of head
#         radii (float array): array of exposure zone radii
#         zone_angle (float): radian angle between orientation and side of exposure zones
#         zone_height (float):
#         num_exp_zones (int): number of exposure zones
#         viability (float): viability index
#         orientation (float): cartesian vector of subject orientation
#     '''
#     if point[2] >= center[2] and point[2] <= (center[2]+zone_height):
#         point  = point[0:2]
#         center = center[0:2]
#         vel = point - center
#         r_p = np.sqrt(dot(vel,vel))
#         if r_p <= radii[num_exp_zones-1]:
#             mag_p     = np.sqrt(np.dot(point,point))
#             mag_c     = np.sqrt(np.dot(center,center))
#             rel_angle = np.arccos(np.sqrt(np.dot(point,center)) / (mag_p * mag_c))
#             if rel_angle <= zone_angle:
#                 weights = np.linspace(100,0,num_exp_zones+1)
#                 for i in range(num_exp_zones):
#                     if i == 0:
#                         if r_p <= radii[0]:
#                             return weights[0]*viability
#                     else:
#                         if r_p <= radii[i] and r_p > radii[i-1]:
#                             return weights[i]*viability
#             else:
#                 return 0
#         else:
#             return 0
#     else:
#         return 0


def createZoneDictionary(a_ZoneDirectory, a_OccupantNames, a_ExposedIndices):

    #------------------------
    # zone sorting algorithm
    #------------------------
    def sortZones(a_FilesList, a_Dir):

        def calculateZoneVertices(a_STLFile, a_ZoneDirectory):
            reader = vtk.vtkSTLReader()
            reader.SetFileName(a_ZoneDirectory+a_STLFile)
            reader.Update()
            zone = reader.GetOutput()
            zoneNumPoints = zone.GetNumberOfPoints()
            zoneVertices = np.zeros((zoneNumPoints, 3))
            for i in range(zoneNumPoints):
                zoneVertices[i] = zone.GetPoint(i)
            return zoneVertices

        numFiles  = len(a_FilesList)
        sortCount = numFiles-1
        zoneVertices = []
        zoneVolumes  = []

        for i in range(numFiles):
            zoneVertices.append(calculateZoneVertices(a_FilesList[i], a_Dir))
            zoneVolumes.append(ConvexHull(zoneVertices[i]).volume)
            # tempVerts = calculateZoneVertices(a_FilesList[i], a_Dir)
            # zoneVertices.append("vertices of "+a_FilesList[i])
            # zoneVolumes.append(ConvexHull(tempVerts).volume)

        for i in range(numFiles-1):
            for j in range(sortCount):
                if zoneVolumes[j] > zoneVolumes[j+1]:
                    zoneVertices[j], zoneVertices[j+1] = zoneVertices[j+1], zoneVertices[j]
                    zoneVolumes[j],  zoneVolumes[j+1]  = zoneVolumes[j+1],  zoneVolumes[j]
            sortCount -= 1

        return zoneVolumes, zoneVertices

    #------------------------------------------
    # create list of names of exposed occupants
    #-------------------------------------------
    filesList = os.listdir(a_ZoneDirectory)
    exposedNames  = []
    for i in range(len(a_ExposedIndices)):
        exposedNames.append(a_OccupantNames[a_ExposedIndices[i]])

    #---------------------------------------------------------------
    # create dictionary of exposure zones for all exposed occupants
    #---------------------------------------------------------------
    numExpZones = int(len(filesList)/len(a_OccupantNames))
    zones = {}

    for i in range(len(exposedNames)):

        tempZones0 = []
        for j in range(len(filesList)):
            if exposedNames[i] in filesList[j]:
                tempZones0.append(filesList[j])

        zoneVolumes, zoneVertices = sortZones(tempZones0, a_ZoneDirectory)

        tempZones1 = {}
        for j in range(numExpZones):
            tempZones1[j] = {"volume": zoneVolumes[j], "vertices": zoneVertices[j]}

        zones[exposedNames[i]] = tempZones1

    return zones



def calculateExposureRisk(a_ZoneData, a_Point, a_Viability):

    def checkPointInZone(a_Volume0, a_Vertices0, a_Point):
        volume1 = ConvexHull(np.concatenate((a_Vertices0, [a_Point]))).volume
        return a_Volume0 == volume1

    numExpZones = len(a_ZoneData)
    weights     = np.linspace(100,0,numExpZones+1)

    for i in range(numExpZones):
        zoneFound = checkPointInZone(a_ZoneData[i]["volume"], a_ZoneData[i]["vertices"], a_Point)
        if zoneFound:
            return weights[i]*a_Viability

    if not zoneFound:
        return 0



class ExposureInputs:

    def __init__(self, FileName):
        self.m_InputFile = FileName
        temp = FileName.rfind('/')
        if temp != -1:
            self.m_RootPath = FileName[0:temp]+'/'
        else:
            self.m_RootPath = ''

        self.m_VariationalParameters = {}
        self.m_VariationalParametersEntry = {}
        self.m_HeaderNames = {}
        self.m_HeaderEntries = {}

    def readInputFile(self):
        inputFileObj = open(self.m_InputFile)
        for line in inputFileObj:

            if not line.startswith('#'):

                lineDict = line.split('=')

                if lineDict[0].strip() == 'Variational parameters':
                    tempList = lineDict[1].split(';')
                    self.m_VaryingParameters = tempList[0].strip() == 'True'

                    if self.m_VaryingParameters:
                        parameters = tempList[1:]
                        self.m_NumVariationalParameters = len(tempList[1:])
                        for i in range(self.m_NumVariationalParameters):
                            parameter = parameters[i].strip().split(':')
                            self.m_VariationalParameters[i] = parameter[0].strip()
                            self.m_VariationalParametersEntry[i] = parameter[1].strip()

                elif lineDict[0].strip() == 'Emitting subjects':
                    self.m_EmittingSubjects = lineDict[1].strip()

                    if self.m_VaryingParameters:
                        for i in range(self.m_NumVariationalParameters):
                            self.m_EmittingSubjects = self.m_EmittingSubjects.replace(self.m_VariationalParameters[i], self.m_VariationalParametersEntry[i])

                elif lineDict[0].strip() == 'Exposed subjects':
                    self.m_ExposedSubjects = lineDict[1].strip()

                    if self.m_VaryingParameters:
                        for i in range(self.m_NumVariationalParameters):
                            self.m_ExposedSubjects = self.m_ExposedSubjects.replace(self.m_VariationalParameters[i], self.m_VariationalParametersEntry[i])

                elif lineDict[0].strip() == 'Subject head centers':
                    self.m_HeadCenters = lineDict[1].strip()

                    if self.m_VaryingParameters:
                        for i in range(self.m_NumVariationalParameters):
                            self.m_HeadCenters = self.m_HeadCenters.replace(self.m_VariationalParameters[i], self.m_VariationalParametersEntry[i])

                elif lineDict[0].strip() == 'Subject orientations':

                    self.m_Orientations = lineDict[1].strip()

                    if self.m_VaryingParameters:
                        for i in range(self.m_NumVariationalParameters):
                            self.m_Orientations = self.m_Orientations.replace(self.m_VariationalParameters[i], self.m_VariationalParametersEntry[i])

                    for i in range(len(self.m_Orientations)):
                        orientation = self.m_Orientations[i].strip()

                        'check for orientations along primary axes'
                        if orientation == '+x':
                            self.m_Orientations[i] = '[1, 0, 0]'
                        elif orientation == '-x':
                            self.m_Orientations[i] = '[-1, 0, 0]'
                        elif orientation == '+y':
                            self.m_Orientations[i] = '[0, 1, 0]'
                        elif orientation == '-y':
                            self.m_Orientations[i] = '[0, -1, 0]'
                        elif orientation == '+z':
                            self.m_Orientations[i] = '[0, 0, 1]'
                        elif orientation == '-z':
                            self.m_Orientations[i] = '[0, 0, -1]'

                elif lineDict[0].strip() == 'Subject names':

                    self.m_SubjectNames = lineDict[1].strip().split(';')

                    for i in range(len(self.m_SubjectNames)):
                        self.m_SubjectNames[i] = self.m_SubjectNames[i].strip()

                elif lineDict[0].strip() == 'VTK results directory':
                    self.m_ResultsDirectory = lineDict[1].strip()

                    if self.m_VaryingParameters:
                        for i in range(self.m_NumVariationalParameters):
                            self.m_ResultsDirectory = self.m_ResultsDirectory.replace(self.m_VariationalParameters[i], self.m_VariationalParametersEntry[i])

                elif lineDict[0].strip() == 'File Step Size':
                    self.m_FileStepSize = lineDict[1].strip()

                    if self.m_VaryingParameters:
                        for i in range(self.m_NumVariationalParameters):
                            self.m_FileStepSize = self.m_FileStepSize.replace(self.m_VariationalParameters[i], self.m_VariationalParametersEntry[i])

                    self.m_FileStepSize = int(self.m_FileStepSize)

                elif lineDict[0].strip() == 'Viability':
                    tempList = lineDict[1].strip().split(';')
                    self.m_LinearViability = tempList[0].split(':')[1].strip() == 'True'
                    self.m_DecayRate       = tempList[1].split(':')[1].strip()
                    self.m_ViabilityDelta  = tempList[2].split(':')[1].strip()

                    if self.m_VaryingParameters:
                        for i in range(self.m_NumVariationalParameters):
                            self.m_DecayRate      = self.m_DecayRate.replace(self.m_VariationalParameters[i], self.m_VariationalParametersEntry[i])
                            self.m_ViabilityDelta = self.m_ViabilityDelta.replace(self.m_VariationalParameters[i], self.m_VariationalParametersEntry[i])

                    self.m_DecayRate      = float(self.m_DecayRate)
                    self.m_ViabilityDelta = float(self.m_ViabilityDelta)

                # elif lineDict[0].strip() == 'Exposure zones':
                #     tempList = lineDict[1].strip().split(';')
                #     self.m_ExposureZoneShape = tempList[0].split(':')[1].strip().lower()
                #     self.m_NumExposureZones  = tempList[1].split(':')[1].strip()
                #     self.m_Zone0Radius       = tempList[2].split(':')[1].strip()
                #     if len(tempList) >= 4:
                #         self.m_ZoneAngle = tempList[2].split(':')[1].strip()
                #     if len(tempList) == 5:
                #         self.m_ZoneHeight = tempList[3].split(':')[1].strip()
                #
                #     if self.m_VaryingParameters:
                #         for i in range(self.m_NumVariationalParameters):
                #             self.m_NumExposureZones = self.m_NumExposureZones.replace(self.m_VariationalParameters[i], self.m_VariationalParametersEntry[i])
                #             if len(tempList) == 4:
                #                 self.m_ZoneAngle = self.m_ZoneAngle.replace(self.m_VariationalParameters[i], self.m_VariationalParametersEntry[i])
                #
                #     self.m_NumExposureZones = int(self.m_NumExposureZones)
                #     if len(tempList) == 3:
                #         self.m_ZoneAngle = float(self.m_ZoneAngle)

                elif lineDict[0].strip() == 'Exposure zones directory':
                    self.m_ExposureZoneFolder = lineDict[1].strip()

                elif lineDict[0].strip() == 'Exposure risk filename':
                    self.m_ExpFileName = lineDict[1].strip()

                    if self.m_VaryingParameters:
                        for i in range(self.m_NumVariationalParameters):
                            self.m_ExpFileName= self.m_ExpFileName.replace(self.m_VariationalParameters[i], self.m_VariationalParametersEntry[i])

                elif lineDict[0].strip() == 'Exposure risk file case identifiers':
                    tempList = lineDict[1].strip().split(';')

                    for i in range(len(tempList)):
                        self.m_HeaderNames[i]   = tempList[i].strip().split(':')[0].strip()
                        self.m_HeaderEntries[i] = tempList[i].strip().split(':')[1].strip()

                        if self.m_VaryingParameters:
                            for j in range(self.m_NumVariationalParameters):
                                self.m_HeaderEntries[i] = self.m_HeaderEntries[i].replace(self.m_VariationalParameters[j], self.m_VariationalParametersEntry[j])

                elif lineDict[0].strip() == 'Grid animation filename':
                    self.m_AnimationFilename = lineDict[1].strip()

                    if self.m_VaryingParameters:
                        for i in range(self.m_NumVariationalParameters):
                            self.m_AnimationFilename = self.m_AnimationFilename.replace(self.m_VariationalParameters[i], self.m_VariationalParametersEntry[i])

                elif lineDict[0].strip() == 'Grid settings':
                    tempList = lineDict[1].strip().split(';')
                    self.m_GenerateHeatmaps = tempList[0].split(':')[1].strip() == 'True'

                    if self.m_GenerateHeatmaps:
                        x_tuple  = tempList[2].split(':')[1].strip().split(',')
                        y_tuple  = tempList[3].split(':')[1].strip().split(',')
                        self.m_GridResolution   = float(tempList[1].split(':')[1].strip())
                        self.m_GridXRange       = (float(x_tuple[0].strip()), float(x_tuple[1].strip()))
                        self.m_GridYRange       = (float(y_tuple[0].strip()), float(y_tuple[1].strip()))

                elif lineDict[0].strip() == 'Grid background':
                    tempList = lineDict[1].strip().split(';')
                    self.m_IncludeBackground = tempList[0].split(':')[1].strip() == 'True'

                    if self.m_IncludeBackground:
                        self.m_BackgroundImage   = tempList[1].split(':')[1].strip()

                        if self.m_VaryingParameters:
                            for i in range(self.m_NumVariationalParameters):
                                self.m_BackgroundImage = self.m_BackgroundImage.replace(self.m_VariationalParameters[i], self.m_VariationalParametersEntry[i])

        inputFileObj.close()
