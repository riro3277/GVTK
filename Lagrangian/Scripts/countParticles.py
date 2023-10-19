# Standalone python script to get particle counts on end caps in an anatomical vascular model - uses the collision code/embolus transport simulation results

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

    
#-----------------------------------------------------------------------------------------
## @brief A utility function that takes a VTK file as input, and extracts the area, 
# centroid coordinates, and radius measures
#
# @param[in] a_FullInputFile Name (including path) of the VTK mesh/data file
# @param[in] a_ReturnArea   Set to 1 if area magnitude is to be returned
# @param[in] a_ReturnCentroid Set to 1 if centroid coordinates are to be returned
# @param[in] a_ReturnBoundingBox Set to 1 if face bounding box is to be returned
# @param[in] a_ReturnBoundingRadius Set to 1 if face bounding radius is to be returned
# @param[in] a_ReturnEffectiveRadius Set to 1 if face effective radius is to be returned
# @param[in] a_IsAvailableNormals Set True if surface mesh has normals pre-computed
# @param[in] a_UseVTKAreaCalculator Set True if in built VTK functions are to be used for
# computing face area and centroid
# @param[in] a_Verbose Set to True if you need detailed outputs of code execution status
#---------------------------------------------------------------------------------------
def vtkGetAreaCentroid(a_FullInputFile, \
                        a_ReturnArea            = 1, \
                        a_ReturnCentroid        = 1, \
                        a_ReturnBoundingBox     = 0, \
                        a_ReturnBoundingRadius  = 1,
                        a_ReturnEffectiveRadius = 1,\
                        a_IsAvailableNormals    = True,\
                        a_UseVTKAreaCalculator  = True,\
                        a_Verbose               = False):

    #
    # create an empty dictionary for the output variables
    # with easily identifiable keys
    #
    returnKeys = ['area', 'centroid', 'bbox', 'bound-radius', 'eff-radius']
    returnDict = dict.fromkeys(returnKeys)

    #
    # read the original polydata file
    #
    if a_FullInputFile.endswith('.vtp'):

      reader = vtk.vtkXMLPolyDataReader()
      reader.SetFileName(a_FullInputFile)
      reader.Update()

    elif a_FullInputFile.endswith('.vtk'):
      
      reader = vtk.vtkPolyDataReader()
      reader.SetFileName(a_FullInputFile)
      reader.ReadAllScalarsOn()
      reader.Update()


    #
    # extract the data from this polydata file
    #
    data    = reader.GetOutput()
    numPts  = data.GetNumberOfPoints()
    numElem = data.GetNumberOfPolys()

    #
    # perform a computation of point normals and obtaining them as a vtkDoubleArray()
    # or a vtkFloatArray()
    #
    if a_IsAvailableNormals:
        normalVec = data.GetPointData().GetNormals()
    else:
        normals = vtk.vtkPolyDataNormals()
        if vtk.VTK_MAJOR_VERSION <= 5.0:
            normals.SetInput(data)
        else:
            normals.SetInputData(data)
        normals.ComputePointNormalsOn()
        normals.ComputeCellNormalsOff()
        normals.Update()

        dataNew   = normals.GetOutput()

        normalVec = vtk.vtkDoubleArray()
        normalVec = dataNew.GetPointData().GetNormals()

    numVectors = normalVec.GetNumberOfTuples()
    numComp    = normalVec.GetNumberOfComponents()
    
    if a_Verbose:
        print('The following are the normal vector tuples')
        for i in range(numVectors):
            print(normalVec.GetTuple(i))

    #
    # now calculate the centroid
    #
    sumX = 0.0
    sumY = 0.0
    sumZ = 0.0

    for p in range(numPts):

        xyz = data.GetPoint(p)
    
        if p == 0:
            xP = np.array(xyz[0])
            yP = np.array(xyz[1])
            zP = np.array(xyz[2])
        else:
            xP = np.append(xP, xyz[0])
            yP = np.append(yP, xyz[1])
            zP = np.append(zP, xyz[2])

        sumX = sumX + xyz[0]
        sumY = sumY + xyz[1]
        sumZ = sumZ + xyz[2]

    centroid = np.array([sumX/float(numPts), sumY/float(numPts), sumZ/float(numPts)])
    
    if a_Verbose:
        print('Centroid coordinates are {0:5.8f}, {1:5.8f}, {2:5.8f}'.\
                format(centroid[0], centroid[1], centroid[2]))
    #
    # now calculate a bounding radius (or bounding box coordinates) for the inlet plane
    #
    for p in range(numPts):

        dist = np.sqrt((xP[p] - centroid[0])**2 + (yP[p] - centroid[1])**2 + (zP[p] - centroid[2])**2)

        if p == 0:
            xMax = xP[p]
            yMax = yP[p]
            zMax = zP[p]
            xMin = xP[p]
            yMin = yP[p]
            zMin = zP[p]
            radMax = dist
        else:

            if xP[p] > xMax:
                xMax = xP[p]

            if xP[p] < xMin:
                xMin = xP[p]

            if yP[p] > yMax:
                yMax = yP[p]

            if yP[p] < yMin:
                yMin = yP[p]

            if zP[p] > zMax:
                zMax = zP[p]

            if zP[p] < zMin:
                zMin = zP[p]

            if dist > radMax:
                radMax = dist

    if a_Verbose:
        print('The bounding box in x is [min={0:5.8f}, max={0:5.8f}]'.format(xMin, xMax))
        print('The bounding box in y is [min={0:5.8f}, max={0:5.8f}]'.format(yMin, yMax))
        print('The bounding box in z is [min={0:5.8f}, max={0:5.8f}]'.format(zMin, zMax))
        print('The bounding radius of the plane is {0:5.8f}'.format(radMax))

    if a_UseVTKAreaCalculator:

        for p in range(numElem):
            
            cell  = data.GetCell(p)
            a     = cell.ComputeArea()
            # print(cell)

            if p == 0:
                area = np.array(a)
            else:
                area = np.append(area, a)

        faceArea = np.sum(area)

    else:

        for p in range(numElem):

            numVerts = data.GetCell(p).GetNumberOfPoints()

            for n in range(numVerts):
        
                xyz = data.GetCell(p).GetPoints().GetPoint(n)

                if n == 0:
                    xF = np.array(xyz[0])
                    yF = np.array(xyz[1])
                    zF = np.array(xyz[2])
                else:
                    xF = np.append(xF, xyz[0])
                    yF = np.append(yF, xyz[1])
                    zF = np.append(zF, xyz[2])

            #
            # this part assumes that the faces are triangles and uses Heron's formula 
            # TODO: need to generalize this to any face type (potentially for working 
            # with Fluent solution data)
            #
            a       = np.sqrt((xF[0]-xF[1])**2 + (yF[0]-yF[1])**2 + (zF[0]-zF[1])**2)
            b       = np.sqrt((xF[1]-xF[2])**2 + (yF[1]-yF[2])**2 + (zF[1]-zF[2])**2)
            c       = np.sqrt((xF[2]-xF[0])**2 + (yF[2]-yF[0])**2 + (zF[2]-zF[0])**2)
            s       = 0.5*(a+b+c)

            if p == 0:
                area    = np.array(np.sqrt(s*(s-a)*(s-b)*(s-c)))
            else:
                area = np.append(area, np.sqrt(s*(s-a)*(s-b)*(s-c)))

        faceArea = np.sum(area)
    
    if a_Verbose:
        print('Total face area is {0:5.8f}'.format(faceArea))
        print('The effective radius of the plane is {0:5.8f}'.format(np.sqrt(faceArea/np.pi)))

    #
    # update the return dictionary based on the choices entered in the function
    ##
    if a_ReturnArea:
        returnDict['area'] = faceArea

    if a_ReturnCentroid:
        returnDict['centroid'] = centroid

    if a_ReturnBoundingBox:
        returnDict['bbox'] = [xMin, xMax, yMin, yMax, zMin, zMax]

    if a_ReturnBoundingRadius:
        returnDict['bound-radius'] = radMax

    if a_ReturnEffectiveRadius:
        returnDict['eff-radius'] = np.sqrt(faceArea/np.pi)


    return returnDict



#------------------------------------------------------------------------
## @brief Function to calculate the distance of a point to a sphere center
#
# @param[in] a_XYZ  The coordinates of the point (3D)
# @param[in] a_C    The coordinates of the sphere center (3D)
# @param[out] dist  The distance between point and center
#
#------------------------------------------------------------------------
def distancePointToSphere(a_XYZ, a_C):

    dist = np.sqrt((a_XYZ[0] - a_C[0])*(a_XYZ[0] - a_C[0]) \
            + (a_XYZ[1] - a_C[1])*(a_XYZ[1] - a_C[1]) \
            + (a_XYZ[2] - a_C[2])*(a_XYZ[2] - a_C[2]))

    return dist

#-----------------------------------------------------------------------
## @brief Function to check whether a point is inside/outside of a sphere
#
# @param[in] a_XYZ  The coordinates of the point (3D)
# @param[in] a_C    The coordinates of the center of the sphere (3D)
# @param[in] a_R    The radius of the sphere
# @param[out] bool  Boolean variable true if point is inside sphere 
#
#-----------------------------------------------------------------------
def pointInSphere(a_XYZ, a_C, a_R):

    if distancePointToSphere(a_XYZ, a_C) <= a_R:
        return True
    else:
        return False

#-------------------------------------------------------------------------------------------------
## @brief Function to run through all the outlet faces in a specified directory and generate an 
# ASCII data file that has all the outlet centroid, and bounding radius information enlisted.
#
# @param[in] a_FacesDirectory   Directory, with full path, where all the face vtp files are stored 
# @param[in] a_OutputDirectory  Directory where the output file is to be dumped
# @param[in] a_OutputFileName   Name of the output file created from this function
# @param[in] a_NumFaces         Total number of faces/outlets for computation
#
#-------------------------------------------------------------------------------------------------
def generateCentroidRadius(a_FacesDirectory, a_OutputDirectory, a_OutputFileName, a_NumFaces, a_VesselNames=None):
    
    if not(a_OutputDirectory.endswith('/')):
        a_OutputDirectory = a_OutputDirectory + '/'

    faceCount = 0
    
    for fileName in os.listdir(a_FacesDirectory):
        if fileName.endswith('.vtk'):
            faceCount = faceCount + 1

    if faceCount != a_NumFaces:
        sys.exit("Incompatible number of vtp files and number of faces")

    #
    # generate output file in directory names Outlet-Bounds which is at the 
    # same level as the Faces directory
    #
    fileObj     = open(a_OutputDirectory+a_OutputFileName, 'w')

    #
    # loop through the face vtk/vtp files, and write the data into the output file
    #
    idV = 0
    for fileName in os.listdir(a_FacesDirectory):

        if ( fileName.endswith('.vtk') or fileName.endswith('.vtp') ):

            returnData  = vtkGetAreaCentroid(a_FacesDirectory+fileName, a_ReturnArea = 1, a_ReturnCentroid = 1, 
                                            a_ReturnBoundingBox = 0, a_ReturnBoundingRadius = 1, 
                                            a_ReturnEffectiveRadius = 0, a_IsAvailableNormals=False)
            centroid    = returnData['centroid']
            radius      = returnData['bound-radius']

            if a_VesselNames is None:

                strOutput   = "FileName: {0:25} Centroid: {1:16} {2:16} {3:16} Bounding-Radius: {4:16}\n".\
                        format(fileName, str(centroid[0]), str(centroid[1]), str(centroid[2]), str(radius))

            else:
                
                strOutput   = "FileName: {0:25} Centroid: {1:16} {2:16} {3:16} Bounding-Radius: {4:16}\n".\
                        format(a_VesselNames[idV], str(centroid[0]), str(centroid[1]), str(centroid[2]), str(radius))


            fileObj.write(strOutput)

            print(strOutput)
            
            idV = idV + 1

    fileObj.close()



#--------------------------------------------------------------------------------------------------
## @brief Function to convert the vessel data from a file to a list of dictionaries that can
# be used for further numerical and post-processing calculations
# 
# @note Returned dictionary has 4 fields: center coordinates, sphere radius, name of vessel
# and id index of the vessel
#
# @param[in] a_CenterRadiusFileName Filename listing the centers and bounding sphere 
# radiuses for each outlet 
# @param[in] a_RadiusPadding (opt) A percentage/proportion of the radius added to the bounding 
# sphere in case it is originally not big enough to capture all particles 
# @param[in] a_VesselNames (opt) List of all actual vessel names for the vessel outlets
#
# @param[out] dictList List of dictionaries containing center and radius data of each vessel outlet
#
#--------------------------------------------------------------------------------------------------
def getCenterRadiusDictionaryFromFile(a_CenterRadiusFileName, a_RadiusPadding=None, a_VesselNames=None):

    fileObj     = open(a_CenterRadiusFileName, 'r')

    dictList    = []
    idV         = 0

    for line in fileObj:
        
        tempDict = {}

        lineList    = line.split()
        center      = np.asarray([float(lineList[3]), float(lineList[4]), float(lineList[5])])
        radius      = float(lineList[7])

        if a_RadiusPadding != None:
            radius = radius + a_RadiusPadding*radius

        fileList    = lineList[1].split('.')

        tempDict['center']  = center
        tempDict['radius']  = radius
        
        if a_VesselNames is None:
            tempDict['name']    = fileList[0]
        else:
            tempDict['name']    = a_VesselNames[idV]

        tempDict['index']   = idV

        dictList.append(tempDict)

        idV = idV+1
        
    return dictList



#------------------------------------------------------------------------------------------------------------
# @brief Function to calculate the number of particles that have exited through all the vessel outlets
#
# @param[in] a_LastParticleVTK Full filename of the last particle vtk file, which has the final fate of all 
# the particles
#
# @param[inout] a_CenterRadiusDictionary Dictionary of vessel types, that is modified by entering the number
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
def calculateParticlesExitingVessel(a_LastParticleVTK, a_CenterRadiusDictionary, a_OutputFile=None):
    
    #
    # create a reader based on the file name and file type (it is assumed that 
    # particle data files will all be of polydata type)
    #
    if a_LastParticleVTK.endswith('.vtp'):
        reader = vtk.vtkXMLPolyDataReader()
    elif a_LastParticleVTK.endswith('.vtk'):
        reader = vtk.vtkPolyDataReader()
    else:
        print("Unsupported File Extension")
        sys.exit()

    #
    # read and extract the coordinates for all the points
    #
    reader.SetFileName(a_LastParticleVTK)
    reader.Update()

    data      = reader.GetOutput()
    numPoints = data.GetNumberOfPoints()

    for p in range(numPoints):

        xyz = data.GetPoint(p)

        for vesselDict in a_CenterRadiusDictionary:

            center = vesselDict['center']
            radius = vesselDict['radius']

            if pointInSphere(np.asarray(xyz), center, radius):
                if 'num-particles' in vesselDict:
                    vesselDict['num-particles'] += 1
                else:
                    vesselDict['num-particles'] = 1

    if not(a_OutputFile is None):
      fileObject = open(a_OutputFile,'w')

    outputDict = {}

    for vesselDict in a_CenterRadiusDictionary:
        
        if 'num-particles' in vesselDict:
            strOut = "Vessel: {0:25} Particles: {1:10} Of {2:10}\n".\
                    format(vesselDict['name'], str(vesselDict['num-particles']), str(numPoints))
            outputDict[vesselDict['name']] = vesselDict['num-particles']
        else:
            strOut = "Vessel: {0:25} Particles: {1:10} Of {2:10}\n".\
                    format(vesselDict['name'], str(0), str(numPoints))
            outputDict[vesselDict['name']] = 0

        if not(a_OutputFile is None):
          fileObject.write(strOut)
    
    if not(a_OutputFile is None):
      fileObject.close()
      
    return outputDict

####################################################################################################################################
####################################################################################################################################
####################################################################################################################################

if __name__=="__main__":

    faces = '/mnt/e/12_LVAD_velocityFiles/Case-23/Caps/'
    facesFiles = os.listdir(faces)
    outDir = '/mnt/e/12_LVAD_velocityFiles/Case-23/'
    outFile = 'radii.txt'
    numFaces = len(facesFiles)

    generateCentroidRadius(faces, outDir, outFile, numFaces)
    
    dictList = getCenterRadiusDictionaryFromFile(outDir+'radii.txt')
    lastVTK = '/mnt/e/12_LVAD_velocityFiles/Case-23/Output1/tracer_Injection-0_140000.vtk'

    calculateParticlesExitingVessel(lastVTK, dictList, outDir+'particles.txt')
