#-----------------------------------------------------------------------
# This module implements a range of geometry based checks for resolving 
# the contact/intersection between a sphere and a triangular element
#
# TODO: 
# - need to add a label for which voronoi region the point falls into
#
# Author:       Debanjan Mukherjee
# Last Edited:  June 2018
#-----------------------------------------------------------------------

#-----------------------------
# IMPLEMENTATION TO BE DECIDED
#-----------------------------
class Label(object):

    def __init__(self):
        self.m_R0 = 0

#-----------------------------
# IMPLEMENTATION TO BE DECIDED
#-----------------------------
class Triangle(object):
    
    def __init__(self, a_P0, a_P1, a_P2):

        self.m_P0   = a_P0
        self.m_P1   = a_P1
        self.m_P2   = a_P2
        self.m_U2   = None  # np.zeros(3)
        self.m_U1   = None  # np.zeros(3)
        self.m_U0   = None  # np.zeros(3)
        self.m_Q0   = None  # np.zeros(2)
        self.m_Q1   = None  # np.zeros(2)
        self.m_Q2   = None  # np.zeros(2)
        self.m_E0   = None  # np.zeros(2)
        self.m_E1   = None  # np.zeros(2)
        self.m_E2   = None  # np.zeros(2)
        self.m_N0   = None  # np.zeros(2)
        self.m_N1   = None  # np.zeros(2)
        self.m_N2   = None  # np.zeros(2)

#--------------------------------------------------------------------------------
# Finds the point that is closest to external point a_X on triangle ABC
# Implemenation based on computation described in "Real-Time Collision Detection"
#--------------------------------------------------------------------------------
def closestPointOnTriangleToPoint_1(a_X, a_A, a_B, a_C):
    
    ab  = a_B - a_A
    ac  = a_C - a_A
    bc  = a_C - a_B

    #
    # compute parametric position 's' for projection X' of X on AB
    # X' = A + s*AB, s = snom/(snom + sdenom)
    #
    snom    = np.dot(a_X - a_A, ab)
    sdenom  = np.dot(a_X - a_B, a_A - a_B)

    #
    # compute parametric position 't' for projection X' of X on AC
    # X' = X + A + t*AC, t = tnom/(tnom + tdenom)
    #
    tnom    = np.dot(a_X - a_A, ac)
    tdenom  = np.dot(a_X - a_C, a_A - a_C)

    #
    # early identification of a vertex region A
    #
    if (snom <= 0.0 and tnom <= 0.0):
        return a_A

    #
    # compute parametric position 'u' for projection X' of X on BC
    # X' = B + u*BC, u = unom/(unom + udenom)
    #
    unom    = np.dot(a_X - a_B, bc)
    udenom  = np.dot(a_X - a_C, a_B - a_C)

    #
    # early identification of vertex region B
    #
    if (sdenom <= 0.0 and unom <= 0.0):
        return a_B

    #
    # early identification of vertex region C
    #
    if (tdenom <= 0.0 and udenom <= 0.0):
        return a_C

    #
    # X is outside (or on) edge AB is the triple scalar product
    # [N XA XB] <= 0.0
    #
    n   = np.cross(a_B - a_A, a_C - a_A)
    vc  = np.dot(n, np.cross(a_A - a_X, a_B - a_X))

    #
    # If X is outside edge AB and within feature region of edge AB
    # then return the projection of X onto AB
    #
    if (vc <= 0.0 and snom >= 0.0 and sdenom >= 0.0):
        return a_A + (snom/(snom+sdenom))*ab

    #
    # X is outside (or on) edge BC if the triple scalar product
    # [N XB XC] <= 0.0
    #
    va  = np.dot(n, np.cross(a_B - a_X, a_C - a_X))

    #
    # If X is outside edge BC and within feature region of edge BC
    # then return the projection of X onto BC
    #
    if (va <= 0.0 and unom >= 0.0 and udenom >= 0.0):
        return a_B + (unom/(unom+udenom))*bc

    #
    # X is outside (or on) edge CA if the triple scalar product
    # [N, XC XA] <= 0.0
    #
    vb = np.dot(n, np.cross(a_C - a_P, a_A - a_P))

    #
    # If X is outside edge CA and within feature region of CA
    # then return the projection of X onto CA
    #
    if (vb <= 0.0 and tnom >= 0.0 and tdenom >= 0.0):
        return a_A + (tnom/(tnom+tdenom))*ac

    #
    # If none of the above conditions/features, then X must project 
    # inside face region. Compute closest point then using barycentric coordinates
    #
    u = va/(va+vb+vc)
    v = vb/(va+vb+vc)
    w = 1.0 - u - v
    return a_A*u + a_B*v + a_C*w

#----------------------------------------------------------------------------------
# Finds the point that is closest to external point a_X on triangle ABC.
# Implementation based on computation described in "Real-Time Collision Detection"
# This alternative implementation is more efficient as it uses only dot products.
#----------------------------------------------------------------------------------
def closestPointOnTriangleToPoint_2(a_X, a_A, a_B, a_C):

    #
    # First check if X is in the vertex Voronoi region outside vertex A
    #
    ab  = a_B - a_A
    ac  = a_C - a_A
    ax  = a_X - a_A
    
    d1  = np.dot(ab, ax)
    d2  = np.dot(ac, ax)

    if (d1 <= 0.0 and d2 <= 0.0):
        return a_A                  # barycentric coordinates (1,0,0)

    #
    # Then check if X is in the vertex Voronoi region outside vertex B
    #
    bx  = a_X - a_B
    d3  = np.dot(ab, bx)
    d4  = np.dot(ac, bx)

    if (d3 >= 0.0 and d4 <= d3):
        return a_B                  # barycentric coordinates (0,1,0)

    #
    # Then check if X is in edge Voronoi region of AB
    # and if so, return projection of X onto AB
    #
    vc  = d1*d4 - d3*d2
    if (vc <= 0.0 and d1 >= 0.0 and d3 <= 0.0):
        v = d1/(d1-d3)
        return a_A + v*ab           # barycentric coordinates (1-v,v,0)

    #
    # Then check if P is in the vertex Voronoi region outside vertex C
    #
    cx  = a_X - a_C
    d5  = np.dot(ab, cx)
    d6  = np.dot(ac, cx)
    if (d6 >= 0.0 and d5 <= d6):
        return a_C                  # barycentric coordinates (0,0,1)

    #
    # Then check if X is in edge Voronoi region of AC
    # and if so, return projection of X onto AC
    #
    vb  = d5*d2 - d1*d6
    if (vb <= 0.0 and d2 >= 0.0 and d6 <= 0.0):
        w = d2/(d2-d6)
        return a_A + w*ac           # barycentric coordinates (1-w,0,w)

    #
    # Then check if X is in edge Voronoi region of BC
    # and if so, return projection of X onto BC
    #
    va  = d3*d6 - d5*d4
    if (va <= 0.0 and (d4-d3) >= 0.0 and (d5-d6) >= 0.0):
        w = (d4 - d3)/((d4 - d3) + (d5 - d6))
        return a_B + w*(a_C - a_B)  # barycentric coordinates (0,1-w,w)

    #
    # If none of the above features/regions, then x is inside triangle
    # and closest point can be computed through its barycentric coordinates
    #
    denom   = 1.0/(va + vb + vc)
    v       = vb*denom
    w       = vc*denom
    return a_A + ab*v + ac*w

#---------------------------------------------------------------
# Function to check whether sphere intersects triangle
# Function returns True if intersection occurs, False otherwise
#---------------------------------------------------------------
def checkSphereTriangleIntersection(a_C, a_R, a_A, a_B, a_C):

    #
    # find point P on trianble ABC that is closest to sphere center C
    #
    p = closestPointOnTriangleToPoint_2(a_C, a_A, a_B, a_C)

    #
    # sphere and triangle intersect if the (squared) distance from 
    # sphere center to point P is less than the (squared) sphere radius
    #
    v = p - a_C
    
    return np.dot(v,v) <= a_R*a_R

if __name__=="__main__":

    #------------------------------------------------------------------------
    # Here we have provided a bunch of functionality to perform some
    # preliminary checks and tests on the functions programmed in this module
    #------------------------------------------------------------------------
    print "TESTS TO BE IMPLEMENTED"
    
