#-------------------------------------------------------------------------
# A python module to compute standard canonical analytical velocity fields
# Author: Debanjan Mukherjee
# TODO: 1. Add a visualization function for each class
# Last Revised: May 2018
#-------------------------------------------------------------------------
import sys, os
import numpy as np

#-----------------------------------------------------------------------
# A base class that can be used as a template for storing and retrieving 
# flow velocity and pressure fields for classical flows for theoretical 
# and analytical studies.
#-----------------------------------------------------------------------
class Flow(object):

  m_Params = []
  m_Name   = ''

  def __init__(self, a_Name, a_Params=None):

    self.m_Name = a_Name
    if a_Params is not None:
      self.m_Params = a_Params

  def getName(self):

    return self.m_Name
  
  def getParams(self):

    return self.m_Params


#---------------------------------------------------------------------------------
# Derived class that stores the unsteady version of Arnold-Beltrami-Childress flow 
#---------------------------------------------------------------------------------
class UnsteadyABC(Flow):

    def __init__(self, a_A=1.0, a_B=1.0, a_C=1.0):

        self.m_Name     = 'unsteadyabc'
        self.m_Params   = [a_A, a_B, a_C]

    def eval(self, a_X, a_Y, a_Z, a_T):

        v0  = ( self.m_Params[0] + 0.5*a_T*np.sin(np.pi*a_T) )*np.sin(a_Z) + self.m_Params[2]*np.cos(a_Y)
        v1  = self.m_Params[1]*np.sin(a_X) + ( self.m_Params[0] + 0.5*a_T*np.sin(np.pi*a_T) )*np.cos(a_Z)
        v2  = self.m_Params[2]*np.sin(a_Y) + self.m_Params[1]*np.cos(a_X)

        return np.array([v0, v1, v2])

#-------------------------------------------------------------------------------
# Derived class that stores the steady version of Arnold-Beltrami-Childress flow
#-------------------------------------------------------------------------------
class SteadyABC(Flow):
    
    def __init__(self, a_A=np.sqrt(3.0), a_B=np.sqrt(2.0), a_C=1.0):
        
        self.m_Name   = 'abcflow'
        self.m_Params = [a_A, a_B, a_C]
        
    def eval(self, a_X, a_Y, a_Z, a_T):
        
        v0  = self.m_Params[0]*np.sin(a_Z) + self.m_Params[2]*np.cos(a_Y)
        v1  = self.m_Params[1]*np.sin(a_X) + self.m_Params[0]*np.cos(a_Z)
        v2  = self.m_Params[2]*np.sin(a_Y) + self.m_Params[1]*np.cos(a_X)
        
        return np.array([v0, v1, v2])

#-------------------------------------------------------------------
# Derived class that stores the unsteady version of Double Gyre flow
#-------------------------------------------------------------------
class UnsteadyDoubleGyre(Flow):

    def __init__(self, a_Amp=0.1, a_Omega=0.20*np.pi, a_Epsilon=0.25):

        self.m_Name     = 'doublegyre'
        self.m_Params   = [a_Amp, a_Omega, a_Epsilon]

    def eval(self, a_X, a_Y, a_Z, a_T):
        
        a       = self.m_Params[2]*np.sin(self.m_Params[1]*a_T)
        b       = 1.0 - 2.0*self.m_Params[2]*np.sin(self.m_Params[1]*a_T)
        f       = a*a_X*a_X + b*a_X
        df_dx   = 2.0*a*a_X + b

        v0  = -np.pi*self.m_Params[0]*np.sin(np.pi*f)*np.cos(np.pi*a_Y)
        v1  = np.pi*self.m_Params[0]*np.cos(np.pi*f)*np.sin(np.pi*a_Y)*df_dx 
        v2  = 0.0

        return np.array([v0, v1, v2])

#-----------------------------------------------------------------
# Derived class that stores the steady version of Double Gyre flow
#-----------------------------------------------------------------
class SteadyDoubleGyre(Flow):

    def __init__(self, a_A=1.0):

        self.m_Name     = 'steadygyre'
        self.m_Params   = a_A

    def eval(self, a_X, a_Y, a_Z, a_T):

        v0  = -np.pi * np.sin(np.pi*a_X) * np.cos(np.pi*a_Y)
        v1  = np.pi * np.cos(np.pi*a_X) * np.sin(np.pi*a_Y)
        v2  = 0.0

        return np.array([v0, v1, v2])

#-----------------------------------------------------
# Derived class that stores the Lamb-Oseen vortex flow
#-----------------------------------------------------
class LambOseen(Flow):

    def __init__(self, a_Nu=0.001, a_Alpha=1.25643, a_U0=0.25):

        self.m_Name   = 'Lamb-Oseen'
        self.m_Params = [a_Nu, a_Alpha, a_U0]

    def eval(self, a_X, a_Y, a_Z, a_T):

        #m_R0       = np.sqrt(4.0*self.m_Params[0]*abs(a_T))
        m_R0       = 0.005
        a_D        = 0.29/2.0

        a_XY1 = [a_D,  a_D]
        a_XY2 = [-a_D, a_D]
        a_XY3 = [0.0,  a_D-np.sqrt((a_D*2.0)**2-a_D**2)]

        m_R1       = np.sqrt((a_X-a_XY1[0])**2 + (a_Y-a_XY1[1])**2)
        m_Theta1   = np.arctan2((a_X-a_XY1[0]), (a_Y-a_XY1[1]))
        m_R2       = np.sqrt((a_X-a_XY2[0])**2 + (a_Y-a_XY2[1])**2)
        m_Theta2   = np.arctan2((a_X-a_XY2[0]), (a_Y-a_XY2[1]))
        m_R3       = np.sqrt((a_X-a_XY3[0])**2 + (a_Y-a_XY3[1])**2)
        m_Theta3   = np.arctan2((a_X-a_XY3[0]), (a_Y-a_XY3[1]))

        m_ThetaDot1 = self.m_Params[2]*(1.0 + (0.5/self.m_Params[1]))*(m_R1/m_R0)*(1.0-\
                      np.exp(-self.m_Params[1]*(m_R0**2/m_R1**2), dtype=np.float64))
        m_ThetaDot2 = self.m_Params[2]*(1.0 + (0.5/self.m_Params[1]))*(m_R2/m_R0)*(1.0-\
                      np.exp(-self.m_Params[1]*(m_R0**2/m_R2**2),dtype=np.float64))
        m_ThetaDot3 = self.m_Params[2]*(1.0 + (0.5/self.m_Params[1]))*(m_R3/m_R0)*(1.0-\
                      np.exp(-self.m_Params[1]*(m_R0**2/m_R3**2),dtype=np.float64))

        v0 = -(m_ThetaDot1*np.cos(m_Theta1) + m_ThetaDot2*np.cos(m_Theta2) + m_ThetaDot3*np.cos(m_Theta3))
        v1 = m_ThetaDot1*np.sin(m_Theta1) + m_ThetaDot2*np.sin(m_Theta2) + m_ThetaDot3*np.sin(m_Theta3)
        v2 = 0.0

        return np.asarray([v0, v1, v2])

#--------------------------------------------------
# Derived class that stores the Bickley jet
# Parameters from https://doi.org/10.1063/1.3271342
#--------------------------------------------------
class Bickley(Flow):

    def __init__(self, a_U=62.66, a_L=1770.0*1e3, a_R=6371.0*1e3):

        self.m_Name = 'Bickley Jet'
        self.m_Params = [a_U, a_L, a_R]
        self.m_Kn = [2.0/a_R, 4.0/a_R, 6.0/a_R]
        self.m_Cn = [self.m_Cn[2]+((np.sqrt(5.0)-1.0)/2.0)*(self.m_Kn[1]/self.m_Kn[2])*(self.m_Cn[1]-self.m_Cn[2]),
                     0.205*a_U, 0.461*a_U]
        self.m_En = [0.075, 0.4, 0.3]

    def eval(self, a_X, a_Y, a_Z, a_T):

      v0 = self.m_Params[0]*(1/(np.cosh(a_Y/a_L)))**2 + \
      2**self.m_Params[0]*np.tanh(a_Y/self.m_Params[1])*((1/np.cosh((a_Y/self.m_Params[1])))**2)*\
      ((self.m_En[0]*(cos(self.m_Kn[0]*self.m_Cn[0]*a_T)*np.cos(self.m_Kn[0]*a_X)\
        + np.sin(self.m_Kn[0]*self.m_Cn[0]*a_T)*np.sin(self.m_Kn[0]*a_X)))\
      +(self.m_En[1]*(cos(self.m_Kn[1]*self.m_Cn[1]*a_T)*np.cos(self.m_Kn(2)*a_X)\
        + np.sin(self.m_Kn[1]*self.m_Cn[1]*a_T)*np.sin(self.m_Kn[1]*a_X)))\
      +(self.m_En[2]*(cos(self.m_Kn[2]*self.m_Cn[2]*a_T)*np.cos(self.m_Kn[2]*a_X)\
        + np.sin(self.m_Kn[2]*self.m_Cn[2]*a_T)*np.sin(self.m_Kn[2]*a_X))))

      v1 = self.m_Params[0]*self.m_Params[1]*(1/np.cosh((a_Y/self.m_Params[1]))**2)*\
      ((self.m_En[0]*(self.m_Kn[0]*np.cos(self.m_Kn[0]*a_X)*np.sin(self.m_Kn[0]*self.m_Cn[0]*a_T)\
        - np.cos(self.m_Kn[0]*self.m_Cn[0]*a_T)*self.m_Kn[0]*np.sin(self.m_Kn[0]*a_X)))\
      +(self.m_En[1]*(self.m_Kn[1]*np.cos(self.m_Kn[1]*a_X)*np.sin(self.m_Kn[1]*self.m_Cn[1]*a_T)\
        - np.cos(self.m_Kn[1]*self.m_Cn[1]*a_T)*self.m_Kn[1]*np.sin(self.m_Kn[1]*a_X)))\
      +(self.m_En[2]*(self.m_Kn[2]*np.cos(self.m_Kn[2]*a_X)*np.sin(self.m_Kn[2]*self.m_Cn[2]*a_T)\
        - np.cos(self.m_Kn[2]*self.m_Cn[2]*a_T)*self.m_Kn[2]*np.sin(self.m_Kn[2]*a_X))))

      v2 = 0.0

      return np.asarray([v0, v1, v2])
