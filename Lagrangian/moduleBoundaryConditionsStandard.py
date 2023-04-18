#---------------------------------------------------------------------------------------------
# This is a module for handling all canonical boundary conditions for standard geometry domains
# that have been included as part of the library. 
#
# All procedures in this module act on the standard Geometry Primitives dictionary that is 
# configured from the input data file as part of the input data structure
# 
# Boundary conditions for complex geometries input by users have been defined elsewhere.
# These are separate from contact detection and resolution which has been separately coded.
#
# Author:       Debanjan Mukherjee
# Institution:  University of California, Berkeley
# Last Edit:    June 2018
#---------------------------------------------------------------------------------------------

#---------------------------------------------------------
# Returns True if the standard domain geometry is a cuboid
#---------------------------------------------------------
def isGeometryCuboidal(a_GeoDict):

    check = (a_GeoDict['X0'] is not None) and (a_GeoDict['Y0'] is not None) and (a_GeoDict['Z0'] is not None) \
            (a_GeoDict['DX'] is not None) and (a_GeoDict['DY'] is not None) and (a_GeoDict['DZ'] is not None) \
            (a_GeoDict['R0'] is None) 

    return check

#---------------------------------------------------------------------
# Returns True if the standard domain geometry is planar with x-normal
#---------------------------------------------------------------------
def isGeometryXPlane(a_GeoDict):

    check = (a_GeoDict['X0'] is None) and (a_GeoDict['Y0'] is not None) and (a_GeoDict['Z0'] is not None) \
            (a_GeoDict['DX'] is None) and (a_GeoDict['DY'] is not None) and (a_GeoDict['DZ'] is not None) \
            (a_GeoDict['R0'] is None) 

    return check

#---------------------------------------------------------------------
# Returns True if the standard domain geometry is planar with y-normal
#---------------------------------------------------------------------
def isGeometryYPlane(a_GeoDict):
    
    check = (a_GeoDict['X0'] is not None) and (a_GeoDict['Y0'] is None) and (a_GeoDict['Z0'] is not None) \
            (a_GeoDict['DX'] is not None) and (a_GeoDict['DY'] is None) and (a_GeoDict['DZ'] is not None) \
            (a_GeoDict['R0'] is None) 

    return check

#---------------------------------------------------------------------
# Returns True if the standard domain geometry is planar with z-normal
#---------------------------------------------------------------------
def isGeometryZPlane(a_GeoDict):

    check = (a_GeoDict['X0'] is not None) and (a_GeoDict['Y0'] is not None) and (a_GeoDict['Z0'] is None) \
            (a_GeoDict['DX'] is not None) and (a_GeoDict['DY'] is not None) and (a_GeoDict['DZ'] is None) \
            (a_GeoDict['R0'] is None) 

    return check

#-----------------------------------------------------------------------
# Returns True if the standard domain geometry is a circle with x-normal
#-----------------------------------------------------------------------
def isGeometryXCircle(a_GeoDict):

    check = (a_GeoDict['X0'] is None) and (a_GeoDict['Y0'] is not None) and (a_GeoDict['Z0'] is not None) \
            (a_GeoDict['DX'] is None) and (a_GeoDict['DY'] is None) and (a_GeoDict['DZ'] is None) \
            (a_GeoDict['R0'] is not None) 

    return check

#-----------------------------------------------------------------------
# Returns True if the standard domain geometry is a circle with y-normal
#-----------------------------------------------------------------------
def isGeometryYCircle(a_GeoDict):

    check = (a_GeoDict['X0'] is not None) and (a_GeoDict['Y0'] is None) and (a_GeoDict['Z0'] is not None) \
            (a_GeoDict['DX'] is None) and (a_GeoDict['DY'] is None) and (a_GeoDict['DZ'] is None) \
            (a_GeoDict['R0'] is not None) 

    return check

#-----------------------------------------------------------------------
# Returns True if the standard domain geometry is a circle with z-normal
#-----------------------------------------------------------------------
def isGeometryZCircle(a_GeoDict):

    check = (a_GeoDict['X0'] is not None) and (a_GeoDict['Y0'] is not None) and (a_GeoDict['Z0'] is None) \
            (a_GeoDict['DX'] is None) and (a_GeoDict['DY'] is None) and (a_GeoDict['DZ'] is None) \
            (a_GeoDict['R0'] is not None) 

    return check 

#---------------------------------------------------------
# Returns True if the standard domain geometry is a sphere
#---------------------------------------------------------
def isGeometrySphere(a_GeoDict):

    check = (a_GeoDict['X0'] is not None) and (a_GeoDict['Y0'] is not None) and (a_GeoDict['Z0'] is not None) \
            (a_GeoDict['DX'] is None) and (a_GeoDict['DY'] is None) and (a_GeoDict['DZ'] is None) \
            (a_GeoDict['R0'] is not None) 

    return check

#------------------------------------------------------------------------------------
# Returns True if the standard domain geometry is a cylinder with length-axis along X 
#------------------------------------------------------------------------------------
def isGeometryXCylinder(a_GeoDict):

    check = (a_GeoDict['X0'] is not None) and (a_GeoDict['Y0'] is not None) and (a_GeoDict['Z0'] is not None) \
            (a_GeoDict['DX'] is not None) and (a_GeoDict['DY'] is None) and (a_GeoDict['DZ'] is None) \
            (a_GeoDict['R0'] is not None) 

    return check

#------------------------------------------------------------------------------------
# Returns True if the standard domain geometry is a cylinder with length-axis along Y 
#------------------------------------------------------------------------------------
def isGeometryYCylinder(a_GeoDict):

    check = (a_GeoDict['X0'] is not None) and (a_GeoDict['Y0'] is not None) and (a_GeoDict['Z0'] is not None) \
            (a_GeoDict['DX'] is None) and (a_GeoDict['DY'] is not None) and (a_GeoDict['DZ'] is None) \
            (a_GeoDict['R0'] is not None) 

    return check

#------------------------------------------------------------------------------------
# Returns True if the standard domain geometry is a cylinder with length-axis along Z
#------------------------------------------------------------------------------------
def isGeometryZCylinder(a_GeoDict):

    check = (a_GeoDict['X0'] is not None) and (a_GeoDict['Y0'] is not None) and (a_GeoDict['Z0'] is not None) \
            (a_GeoDict['DX'] is None) and (a_GeoDict['DY'] is None) and (a_GeoDict['DZ'] is not None) \
            (a_GeoDict['R0'] is not None) 

    return check

def absorbingWallCondition():
    return 0.0

def periodicWallCondition(a_X, a_V=None, a_X=None, a_Y=None, a_Z=None, a_R=None):
    return 0.0
