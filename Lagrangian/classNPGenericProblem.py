import sys, os
import numpy as np
import vtk

class GNPProblem(object):

    def __init__(self):
        self.m_InputData = None
        self.m_Debugode = None
