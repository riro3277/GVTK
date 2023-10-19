import sys, os, vtk
import numpy as np
import time
import random

from numpy import linalg as la

try:
    from classGVTKGrid import *
    from classGVTKLagrangianData import *
    from classGVTKGenericProblem import *
    from testModule import *
except:
    sys.exit('Could not import user defined modules')


tracer_rootDir = "/mnt/d/lagrangiantest/PoiseuilleFlow/"
tracer_prefix = 'Seeds.vtk'
file = tracer_rootDir + tracer_prefix
outFile = tracer_rootDir + "AnalyticalForces.csv"
ParticleData = GVTKLagrangianData(a_File = file)
R = 1.0
r = 0.25
U = 100
InitVel = [0.0,0.0,0.0]
totalTime = 0.2
dT = 0.00005
steps = totalTime/dT
denF = 0.0011
denP = 0.0011
visF = 0.004
pos0 = np.zeros([ParticleData.numParticles, 3])
velp = np.zeros([ParticleData.numParticles, 3])
for p in range(ParticleData.numParticles):
    pos        = ParticleData.getX(p)
    #print(pos)
    pos0[p] = pos
    #print(pos0)
#print(pos0)
file = open(outFile, "w")
file.close()


for step in range(int(steps)):
#for step in range(1):
    simTime = step*dT
    restitution     = 0.75 #--------getting restitution value from input-------#

    #-------------------------------------------------------------------------------#
    #----------Calculating constants for each term in Maxey Riley Equation----------#
    #-------------------------------------------------------------------------------#
    c_1 = 0.1806
    c_2 = 0.6459
    c_3 = 0.4251
    c_4 = 6886.95
    d_1 = (3.0 * denF)/(8.0 * r)
    C_am    = 1.0
    accGrav = [0.0,0.0,-9806.65]
    # print('gravity=', accGrav)
    # sys.exit()

    #------------------------------------------
    # forward euler integration implementation
    #------------------------------------------
    TDrag = 0.0
    Tshear = 0.0
    Tadded = 0.0
    Tuforce = 0.0
    Tbody = 0.0
    c = 0.0

    for p in range(ParticleData.numParticles):
        # print("Pos0 for P", simTime, p, pos0[p])
        posP = pos0[p]
        #print("here",posP)
        # posP[0] = float(posP[0])
        # posP[1] = float(posP[1])
        # posP[2] = float(posP[2])


        if step == 0:
            v_i = InitVel
        else:
            v_i =velp[p]# particle velocity


        # print('velocityp=', v_i)
        # sys.exit()
        ## error with using the same name for cell assisgnment: FindCell argument 1: 'tuple' object does not support item assignment


        delv_delT    = np.zeros(3)




        v =[0.0, 0.0, U*(((R**2)-((posP[0]**2) + (posP[1]**2)))/(R**2))]
        # print("fluid velocity", v)

        duzdx = (-2.0*(U/(R**2)))*(posP[0])
        duzdy = (-2.0*(U/(R**2)))*(posP[1])
        gradU = [0.0,0.0,0.0,0.0,0.0,0.0,duzdx, duzdy, 0.0]
        delv_delT = [0.0,0.0,0.0]


        #--------Edits for Gradient calculation based on Nodal Basis interpolation is done by Sreeparna----------#
         #----Matrix is converted to arrays edited by Sreeparna----#
        # print("GRAD", gradU, posP)
        #---------Calculate the material derivative of fluid velocity field-----------#
        lhsNS    = np.zeros(3)
        lhsNS[0] = delv_delT[0] + v[0]*gradU[0] + v[1]*gradU[1] + v[2]*gradU[2]
        lhsNS[1] = delv_delT[1] + v[0]*gradU[3] + v[1]*gradU[4] + v[2]*gradU[5]
        lhsNS[2] = delv_delT[2] + v[0]*gradU[6] + v[1]*gradU[7] + v[2]*gradU[8]

        # print(gradU)

        #-------Calculate vorticity from the components of the velocity gradient tensor-------#
        vor      = np.zeros(3)
        vor[0]   = gradU[7] - gradU[5]
        vor[1]   = gradU[2] - gradU[6]
        vor[2]   = gradU[3] - gradU[1]

        vol_p = 4.0*3.1415926535898*((r)**3.0)/3.0
        #-------calculation of particle slip velocity Reynolds number--------#
        # Re_p = (self.m_InputData.getTracerDensity() * la.norm(v - v_i) * 2 * self.m_InputData.m_TracerRadius)/self.m_InputData.getFluidViscosity()
        Re_p = (denF * la.norm(v - v_i) * 2.0 * r)/visF

        #-------calculation of shear based Reynolds number--------#
        Re_g = denF * la.norm(gradU) * ((2.0 * r)**2.0)/visF

        alpha_LSA = 0.5 * (Re_g/(Re_p + 1.0e-20))
        # CoeffLift_extra = (3.0/(4.0 * np.pi * self.m_InputData.m_TracerRadius)) * np.sqrt((self.m_InputData.getFluidDensity() * self.m_InputData.getFluidViscosity()) / (la.norm(vor) + 1.0e-16))
        CoeffLift_extra = denF*(r)**2.0/(vol_p)*np.sqrt((visF/denF)/(la.norm(vor) + 1.0e-20))

        if Re_p == 0:
            coeffSL = 6.46 * CoeffLift_extra
        elif Re_p <= 40.0:
            coeffSL = 6.46 * ((1.0 - 0.3314 * np.sqrt(alpha_LSA)) * np.exp(-Re_p/10.0) + 0.3314 * np.sqrt(alpha_LSA)) * CoeffLift_extra
        else:
            # coeffSL = 6.46 * 0.0524 * np.sqrt(0.5 * Re_g) * CoeffLift_extra
            coeffSL = 6.46 * 0.0524 * np.sqrt(alpha_LSA * Re_p) * CoeffLift_extra

        #-----------calculation of Drag coefficient----------#
        C_D = (24.0/(Re_p + 1.0e-20)) * (1 + c_1*(Re_p + 1.0e-20)**c_2) + c_3/(1 + (c_4/(Re_p + 1.0e-20)))

        #-----------calculation of Drag Force----------#
        drag = d_1 * C_D * la.norm(v - v_i) * (v - v_i)

        #-----------calculation of Shear Lift Force----------#
        shearGradLift = coeffSL * np.cross(v - v_i, vor)

        #--------------Calculation Of Added Mass Force----------#
        ambForce = denF * (1.0 + (C_am/2.0)) * lhsNS
        # print("ambForce",ambForce)
        addedmass = denF * 0.5 * C_am * lhsNS
        undisturbedforce = denF * lhsNS

        #---------Calculation Of Buoyancy-----------#
        bodyForce = int((denP - denF)) * accGrav
        bodyForce = [1e-10,1e-10,1e-10]

        #--------Now sum up all the forces to calculate the total force on the particle----------#

        # print(drag, shearGradLift, ambForce, bodyForce)
        v_i = v_i + (1.0/(denP + denF*(C_am/2.0))) * dT * (drag + shearGradLift + ambForce)
        #print("here",v_i)
        # print(simTime, "Flow Vel", v, p)
        # print(simTime, "Particle Velocity", v_i, p)
        velp[p] = v_i

        #------------Calculation of SDF Gradients----------#
        velfmag = (v[0]**2)+(v[1]**2)+(v[2]**2)
        velpmag = (v_i[0]**2)+(v_i[1]**2)+(v_i[2]**2)
        dragp = drag*vol_p
        addmp= addedmass*vol_p
        undisturbedp = undisturbedforce*vol_p
        bodyforce = 0
        shearGradp = shearGradLift*vol_p

        dragmagp = (dragp[0]**2)+(dragp[1]**2)+(dragp[2]**2)
        addmmagp = (addmp[0]**2)+(addmp[1]**2) + (addmp[2]**2)
        undisturbedmagp = (undisturbedp[0]**2)+(undisturbedp[1]**2)+(undisturbedp[2]**2)
        shearGradmagp = (shearGradp[0]**2)+(shearGradp[1]**2)+(shearGradp[2]**2)
        #---------calculate distance of particle position from wall-----------#





        #-------calculate updated particle position due to change in velocity-------#
        posPnew = posP + v_i*dT

        posPnew = np.array(posPnew)
        #print("pos1",p, posP)
        pos0[p] = posPnew
        #print("new",p, posP)
        # print(pos0)


        print(simTime)
        out = open(outFile, 'a')
        temp = "p,{},t,{},Drag,{},Addm,{},Flow,{},Buoy,{},Lift,{}\n".format(p,round(simTime,5), dragmagp, addmmagp, undisturbedmagp,bodyforce, shearGradmagp)
        # temp = "p,{},t,{},Velf,{},Velp,{}\n".format(p,round(simTime,5), velfmag, velpmag)

        out.write(temp)
