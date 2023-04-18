#--------------------------------------------------------------------------
# This module is a container for a variety of techniques to add or inject 
# Lagrangain particles into the computational domain at specified locations
#
# Author:       Debanjan Mukherjee
# Institution:  University of California, Berkeley
# Last Edit:    June 2018
#--------------------------------------------------------------------------
import sys, os
import vtk
import numpy as np

#--------------------------------------------------------------------------------
# Function to inject points at specified location using a specified method
# This function takea advantage of VTK's in-built dynamic allocation capabilities
# -------
# Params:
# -------
# a_Points:         vtkPoints object to which new points will be added
# a_InputPoints:    vtkPoints object that holds all injected points
# --------
# Returns:
# --------
# a_Points: a vtkPoints object with old and injected points
#--------------------------------------------------------------------------------
def injectPoints(a_Points, a_InputPoints):
    
    oldNumPts   = a_Points.GetNumberOfPoints()
    addNumPts   = a_InputPoints.GetNumberOfPoints() 
        
    for p in range(addNumPts):
        addID = oldNumPts + p
        a_Points.InsertPoint(addID, a_InputPoints.GetPoint(p))

    return a_Points

#---------------------------------------------------------------------------------
# A more elaborate implementation of injectPoints that can handle the injection
# of polydata representations of the Lagrangian particles.
#
# Key constraint: input polydata and copied polydata both need to have equal 
# number of arrays, and the same data array names, and if it does not have this
# code is going to exit with error
#
# This function takes advantage of VTK's in-built dynamic allocation capabilities
# 
# NOTE: 
# - this essentially is a re-implementation of the code that runs under-the-hood
#   for the vtkAppendPolyData() filter
# - for some reason vtkAppendPolyData() for this version does not seem to work
# - future implementations may consider re-implementing using this filter
#---------------------------------------------------------------------------------
def injectParticles(a_Particles, a_InputParticles):
    
    #
    # check whether the input arguments are both of the type polydata, 
    # and that they both have the same number of data arrays - exit otherwise
    #
    if ( (a_Particles.IsA('vtkPolyData') != 1) or (a_InputParticles.IsA('vtkPolyData') != 1) ):
        sys.exit("Module: moduleInjection, Function: injectParticles - both arguments should be polyData")
    
    if ( a_Particles.GetPointData().GetNumberOfArrays() != a_InputParticles.GetPointData().GetNumberOfArrays() ):
        sys.exit("Module: moduleInjection, Function: injectParticles - both arguments need same number of data")

    #
    # get the number of points in each data set
    #
    n_In    = a_Particles.GetNumberOfPoints()
    n_Add   = a_InputParticles.GetNumberOfPoints()

    #
    # add points from the injection dataset into original dataset
    #
    for p in range(n_Add):

        xyz = a_InputParticles.GetPoint(p)
        a_Particles.GetPoints().InsertPoint(p + n_In, xyz)

    #
    # now add tuples from each array in the injection dataset into
    # the corresponding array in the original dataset
    #
    for a in range(a_Particles.GetPointData().GetNumberOfArrays()):

        for p in range(n_Add):
            
            daTuple = a_InputParticles.GetPointData().GetArray(a).GetTuple(p)
            a_Particles.GetPointData().GetArray(a).InsertTuple(p + n_In, daTuple)

    return a_Particles

#---------------------------------------------------------------------------------
# An implementation where points are injected at the locations specified by a 
# set of seed points (given as a vtkPoints object), and each seed point is randomly
# perturbed around it's seed location by a uniform random radius factor
#---------------------------------------------------------------------------------
def injectPointsWithRandomPerturbations(a_Points, a_SeedPoints, a_PerturbationExtent=0.1, a_SeedPlaneNormal=None):
    
    oldNumPts   = a_Points.GetNumberOfPoints()
    addNumPts   = a_SeedPoints.GetNumberOfPoints()

    if a_SeedPlaneNormal is None:
        for p in range(addNumPts):
            xyz     = np.asarray(a_SeedPoints.GetPoint(p))
            xyz[0]  = xyz[0]*(1.0 + np.random.uniform(low=-a_PerturbationExtent, high=a_PerturbationExtent))
            xyz[1]  = xyz[1]*(1.0 + np.random.uniform(low=-a_PerturbationExtent, high=a_PerturbationExtent))
            xyz[2]  = xyz[2]*(1.0 + np.random.uniform(low=-a_PerturbationExtent, high=a_PerturbationExtent))
            a_Points.InsertPoint(oldNumPts + p, xyz)
    else:
        seedPlaneTangent    = np.asarray(a_SeedPoints.GetPoint(0)) - np.asarray(a_SeedPoints.GetPoint(1))
        seedPlaneTangent    = seedPlaneTangent/np.linalg.norm(seedPlaneTangent) 
        seedBiTangent       = np.cross(seedPlaneTangent, a_SeedPlaneNormal/np.linalg.norm(a_SeedPlaneNormal))
        seedBiTangent       = seedBiTangent/np.linalg.norm(seedBiTangent)
        
        for p in range(addNumPts):
            xyz     = np.asarray(a_SeedPoints.GetPoint(p))
            xyz     = xyz + \
                    np.random.uniform(low=-a_PerturbationExtent, hight=a_PerturbationExtent)*seedPlaneTangent + \
                    np.random.uniform(low=-a_PerturbationExtent, high=a_PerturbationExtent)*seedBiTangent
            a_Points.InsertPoint(oldNumPts + p, xyz)

#-----------------------------------------------------------------------------------
# An overloaded implementation of the 'injectPointsWithRandomPerturbations' function
# that operates on polydata representation of particle tracer data
# TO BE IMPLEMENTED
#-----------------------------------------------------------------------------------
def injectParticlesWithRandomPerturbations(a_Particles, a_SeedParticles, a_PerturbationExtent=0.1, a_SeedPlaneNormal=None):


def injectPointsSpecifiedConcentration(a_Points, a_Concentration, a_Location):
    return 0.0

