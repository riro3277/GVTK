# Standalone python script to get particle counts in specified bounded boxes within an anatomical vascular model - uses the collision code/embolus transport simulation results
import sys, os

try:
    import vtk
except ImportError:
    print("Could not import vtk in postProcessorFlowVC. All vtk dependent modules will not work")

try:
    import numpy as np
except ImportError:
    sys.exit("Could not import numpy in postProcessorFlowVC. This needs to be installed")

try:
    import matplotlib.pyplot as plt
except ImportError:
    print("Could not import matplotlib. All plotting commands will not work")


def getBoxCoordinateDictionaryFromFile(a_BoundingBoxFileName):

    fileObj     = open(a_BoundingBoxFileName, 'r')

    dictList    = []
    idV         = 0

    for line in fileObj:

        tempDict = {}

        lineList    = line.split()
        minVal      = np.asarray([float(lineList[3]), float(lineList[4]), float(lineList[5])])
        lengths     = np.asarray([float(lineList[7]), float(lineList[8]), float(lineList[9])])

        fileList    = lineList[1].split('.')

        tempDict['minVal']  = minVal
        tempDict['lengths'] = lengths
        tempDict['name']    = fileList[0]
        tempDict['index']   = idV

        # print('ID = ', idV)
        # print('minVal = ', minVal)
        # print('lengths = ', lengths)

        dictList.append(tempDict)

        idV = idV+1

    return dictList

def pointInBox(a_XYZ, a_min, a_Length):

    # print(a_min[0])
    # print(a_min[1])
    # print(a_min[2])
    # print(a_Length[0], a_Length[1], a_Length[2])

    if a_XYZ[0] >= a_min[0] and a_XYZ[0] <= a_min[0] + a_Length[0]:
       if a_XYZ[1] >= a_min[1] and a_XYZ[1] <= a_min[1] + a_Length[1]:
          if a_XYZ[2] >= a_min[2] and a_XYZ[2] <= a_min[2] + a_Length[2]:
            return True
          else:
            return False

#------------------------------------------------------------------------------------------------------------
# @brief Function to calculate the number of particles that enter a particular box ID
#
# @param[in] a_LastParticleVTK Full filename of the last particle vtk file, which has the final fate of all
# the particles
#
# @param[inout] a_BoxCoordinateDictionary Dictionary of Box IDs, that is modified by entering the number
# of particles that are exiting each vessel
#
# @param[in] a_OutputFile Full filename for the output file where the particle number distribution is written
#
# @param[out] outputDict A modified version of the dictionary generated for each vessel outlet
#
# @note The code adds a new field to the dictionary: 'num-particles' which is the number of particles that
# have exited from the vessel.
#
#------------------------------------------------------------------------------------------------------------
def countParticlesInBoxID(a_LastParticleVTK, a_BoxCoordinateDictionary, a_OutputFile=None):
    #
    # read the original polydata file
    #
    if a_LastParticleVTK.endswith('.vtp'):
      reader = vtk.vtkXMLPolyDataReader()

    elif a_LastParticleVTK.endswith('.vtk'):
      reader = vtk.vtkPolyDataReader()

    else:
      print("Unsupported File Extension")
      sys.exit()

    reader.SetFileName(a_LastParticleVTK)
    reader.Update()

    #
    # extract the data from this polydata file
    #
    data    = reader.GetOutput()
    numPts  = data.GetNumberOfPoints()

    for p in range(numPts):
      posP  = data.GetPoint(p)

      for vesselDict in a_BoxCoordinateDictionary:
        minVal  = vesselDict['minVal']
        lengths = vesselDict['lengths']

        if pointInBox(np.asarray(posP), minVal, lengths):
            if 'num-particles' in vesselDict:
                vesselDict['num-particles'] += 1
            else:
                vesselDict['num-particles'] = 1

    if not(a_OutputFile is None):
      fileObject = open(a_OutputFile,'w')

    outputDict = {}

    for vesselDict in a_BoxCoordinateDictionary:

        if 'num-particles' in vesselDict:
            strOut = "Vessel: {0:25} Particles: {1:10} Of {2:10}\n".\
                    format(vesselDict['name'], str(vesselDict['num-particles']), str(numPts))
            outputDict[vesselDict['name']] = vesselDict['num-particles']
        else:
            strOut = "Vessel: {0:25} Particles: {1:10} Of {2:10}\n".\
                    format(vesselDict['name'], str(0), str(numPts))
            outputDict[vesselDict['name']] = 0

        if not(a_OutputFile is None):
          fileObject.write(strOut)

    if not(a_OutputFile is None):
      fileObject.close()

    return outputDict

if __name__=="__main__":

    outDir = '/mnt/e/12_LVAD_velocityFiles/Case-23/'

    dictList = getBoxCoordinateDictionaryFromFile(outDir+'boxCoordinates.txt')
    lastVTK = '/mnt/e/12_LVAD_velocityFiles/Case-23/Output1/tracer_Injection-0_140000.vtk'

    countParticlesInBoxID(lastVTK, dictList, outDir+'box_particles.txt')


# #-------------------------------------------------------------------------------------------------
# ## @brief Function to run through all the outlet boxes in a specified directory and generate an
# # ASCII data file that has all the outlet center, and bounding lengths information enlisted.
# #
# # @param[in] a_BoxesDirectory   Directory, with full path, where all the box vtp files are stored
# # @param[in] a_OutputDirectory  Directory where the output file is to be dumped
# # @param[in] a_OutputFileName   Name of the output file created from this function
# # @param[in] a_NumBoxes         Total number of boxed/outlets for computation
# #
# #-------------------------------------------------------------------------------------------------
# def generateCentroidRadius(a_BoxesDirectory, a_OutputDirectory, a_OutputFileName, a_NumBoxes, a_VesselNames=None):

#     if not(a_OutputDirectory.endswith('/')):
#         a_OutputDirectory = a_OutputDirectory + '/'

#     BoxCount = 0

#     for fileName in os.listdir(a_BoxesDirectory):
#         if fileName.endswith('.vtk'):
#             BoxCount = BoxCount + 1

#     if BoxCount != a_NumBoxes:
#         sys.exit("Incompatible number of vtp files and number of boxes")

#     #
#     # generate output file in directory names Outlet-Bounds which is at the
#     # same level as the Boxes directory
#     #
#     fileObj     = open(a_OutputDirectory+a_OutputFileName, 'w')

#     #
#     # loop through the box vtk/vtp files, and write the data into the output file
#     #
#     idV = 0
#     for fileName in os.listdir(a_BoxesDirectory):

#         if ( fileName.endswith('.vtk') or fileName.endswith('.vtp') ):

#             returnData  = vtkGetAreaCentroid(a_BoxesDirectory+fileName, a_ReturnArea = 1, a_ReturnCentroid = 1,
#                                             a_ReturnBoundingBox = 0, a_ReturnBoundingRadius = 1,
#                                             a_ReturnEffectiveRadius = 0, a_IsAvailableNormals=False)
#             centroid    = returnData['centroid']
#             radius      = returnData['bound-radius']

#             if a_VesselNames is None:

#                 strOutput   = "FileName: {0:25} Centroid: {1:16} {2:16} {3:16} Bounding-Radius: {4:16}\n".\
#                         format(fileName, str(centroid[0]), str(centroid[1]), str(centroid[2]), str(radius))

#             else:

#                 strOutput   = "FileName: {0:25} Centroid: {1:16} {2:16} {3:16} Bounding-Radius: {4:16}\n".\
#                         format(a_VesselNames[idV], str(centroid[0]), str(centroid[1]), str(centroid[2]), str(radius))


#             fileObj.write(strOutput)

#             print(strOutput)

#             idV = idV + 1

#     fileObj.close()
