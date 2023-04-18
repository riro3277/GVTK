import sys, os, vtk
import numpy as np

class FTLEInputs:

    def __init__(self, a_FileName):

        self.m_InputFile = a_FileName
        temp             = a_FileName.rfind('/')

        if temp != -1:
            self.m_RootPath = a_FileName[0:temp]+'/'
        else:
            self.m_RootPath = ''

    def readInputFile(self):

        inputFileObj = open(self.m_InputFile)

        for line in inputFileObj:

            if not line.startswith('#'):

                lineDict = line.split('=')

                if lineDict[0].strip() == 'Flow Data Directory':

                    self.m_FlowDataDirectory = lineDict[1].strip()

                elif lineDict[0].strip() == 'Flow Data File Tag':

                    self.m_FlowDataFile = lineDict[1].strip()

                elif lineDict[0].strip() == 'Flow data field name':

                    self.m_FlowFieldName = lineDict[1].strip()

                elif lineDict[0].strip() == 'FTLE File Name':

                    self.m_OutputFile = lineDict[1].strip()

                elif lineDict[0].strip() == 'FTLE Choice':

                    self.m_FTLEChoice = lineDict[1].strip()

                elif lineDict[0].strip() == 'Scaled Time':

                    self.m_ScaledTime = lineDict[1].strip().lower() == 'true'

                elif lineDict[0].strip() == 'X Bounds':

                    tempList = lineDict[1].strip().split(';')
                    self.m_xMin = float(tempList[0].strip())
                    self.m_xMax = float(tempList[1].strip())

                elif lineDict[0].strip() == 'Y Bounds':

                    tempList = lineDict[1].strip().split(';')
                    self.m_yMin = float(tempList[0].strip())
                    self.m_yMax = float(tempList[1].strip())

                elif lineDict[0].strip() == 'Z Bounds':

                    tempList = lineDict[1].strip().split(';')
                    self.m_zMin = float(tempList[0].strip())
                    self.m_zMax = float(tempList[1].strip())

                elif lineDict[0].strip() == 'Number of Bins':

                    tempList = lineDict[1].strip().split(';')
                    self.m_xBins = int(tempList[0].strip().split(':')[1].strip())
                    self.m_yBins = int(tempList[1].strip().split(':')[1].strip())
                    self.m_zBins = int(tempList[2].strip().split(':')[1].strip())

                elif lineDict[0].strip() == 'Simulation timing':

                    tempList = lineDict[1].strip().split(';')
                    self.m_StartTime = float(tempList[0].strip().split(':')[1].strip())
                    self.m_StopTime  = float(tempList[1].strip().split(':')[1].strip())
                    self.m_dT        = float(tempList[2].strip().split(':')[1].strip())

        inputFileObj.close()
