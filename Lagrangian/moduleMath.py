
import numpy as np

def getTetrahedralCellGlobalToLocal(a_X1, a_X2, a_X3, a_X4, a_X):

    mat = np.zeros((3,3), dtype=np.float32)
    vec = np.zeros((3,1), dtype=np.float32)
    eta = np.zeros((3,1), dtype=np.float32)

    mat[0,0] = a_X1[0] - a_X4[0]
    mat[0,1] = a_X2[0] - a_X4[0]
    mat[0,2] = a_X3[0] - a_X4[0]

    mat[1,0] = a_X1[1] - a_X4[1]
    mat[1,1] = a_X2[1] - a_X4[1]
    mat[1,2] = a_X3[1] - a_X4[1]

    mat[2,0] = a_X1[2] - a_X4[2]
    mat[2,1] = a_X2[2] - a_X4[2]
    mat[2,2] = a_X3[2] - a_X4[2]
    
    # mat[0,0] = a_X1[0] - a_X4[0]
    # mat[0,1] = a_X1[1] - a_X4[1]
    # mat[0,2] = a_X1[2] - a_X4[2]

    # mat[1,0] = a_X2[0] - a_X4[0]
    # mat[1,1] = a_X2[1] - a_X4[1]
    # mat[1,2] = a_X2[2] - a_X4[2]

    # mat[2,0] = a_X3[0] - a_X4[0]
    # mat[2,1] = a_X3[1] - a_X4[1]
    # mat[2,2] = a_X3[2] - a_X4[2]

    vec[0,0] = a_X[0] - a_X4[0]
    vec[1,0] = a_X[1] - a_X4[1]
    vec[2,0] = a_X[2] - a_X4[2]

    eta = np.linalg.solve(mat,vec)

    return eta

def getTetrahedralCellLocalToGlobal(a_V1, a_V2, a_V3, a_V4, a_Eta):

    phi_1 = a_Eta[0]
    phi_2 = a_Eta[1]
    phi_3 = a_Eta[2]
    phi_4 = 1.0 - a_Eta[0] - a_Eta[1] - a_Eta[2]

# correction here len(a_V1 == 1) is incorrect
    if len(a_V1) == 1:
        val = a_V1*phi_1 + a_V2*phi_2 + a_V3*phi_3 + a_V4*phi_4
    else:
        val = np.zeros(len(a_V1), dtype=np.float32)
        for i in range(len(a_V1)):
            val[i] = a_V1[i]*phi_1 + a_V2[i]*phi_2 + a_V3[i]*phi_3 + a_V4[i]*phi_4

    return val

def getTetrahedralScalarGradient(a_V1, a_V2, a_V3, a_V4, a_X1, a_X2, a_X3, a_X4, a_X):

    vec = np.zeros(3, dtype=np.float32)

    vec[0] = (a_V1-a_V4)/(a_X1[0]-a_X4[0]) + (a_V2-a_V4)/(a_X2[0]-a_X4[0]) + (a_V3-a_V4)/(a_X3[0]-a_X4[0])
    vec[1] = (a_V1-a_V4)/(a_X1[1]-a_X4[1]) + (a_V2-a_V4)/(a_X2[1]-a_X4[1]) + (a_V3-a_V4)/(a_X3[1]-a_X4[1])
    vec[2] = (a_V1-a_V4)/(a_X1[2]-a_X4[2]) + (a_V2-a_V4)/(a_X2[2]-a_X4[2]) + (a_V3-a_V4)/(a_X3[2]-a_X4[2])

    return vec

# def getTetrahedralVectorGradient(a_V1, a_V2, a_V3, a_V4, a_X1, a_X2, a_X3, a_X4, a_X):


#     # mat = np.zeros((len(a_V1,3)), dtype=np.float32)
#     #for arg in args:
#     #    print(arg)
#     #sys.exit()
#     mat = np.zeros((len(a_V1),3), dtype=np.float32)
    
#     for r in range(len(a_V1)):
#         mat[r,0] = (a_V1[r]-a_V4[r])/(a_X1[0]-a_X4[0]+1e-16) + (a_V2[r]-a_V4[r])/(a_X2[0]-a_X4[0]+1e-16) + (a_V3[r]-a_V4[r])/(a_X3[0]-a_X4[0]+1e-16)
#         mat[r,1] = (a_V1[r]-a_V4[r])/(a_X1[1]-a_X4[1]+1e-16) + (a_V2[r]-a_V4[r])/(a_X2[1]-a_X4[1]+1e-16) + (a_V3[r]-a_V4[r])/(a_X3[1]-a_X4[1]+1e-16)
#         mat[r,2] = (a_V1[r]-a_V4[r])/(a_X1[2]-a_X4[2]+1e-16) + (a_V2[r]-a_V4[r])/(a_X2[2]-a_X4[2]+1e-16) + (a_V3[r]-a_V4[r])/(a_X3[2]-a_X4[2]+1e-16)
#     return mat

#---------Edits for Gradient calculation based on Matrix inversion is done by Sreeparna----------#
def getTetrahedralVectorGradient(a_V1, a_V2, a_V3, a_V4, a_X1, a_X2, a_X3, a_X4, a_X):


     mat = np.zeros((3,3), dtype=np.float32)
     vec = np.zeros((3,1), dtype=np.float32)
     eta = np.zeros((3,1), dtype=np.float32)
     F_inv = np.zeros((3,3), dtype=np.float32)
     eta_x = np.zeros((len(a_V1),3), dtype=np.float32)
     

     mat[0,0] = a_X1[0] - a_X4[0]
     mat[0,1] = a_X2[0] - a_X4[0]
     mat[0,2] = a_X3[0] - a_X4[0]

     mat[1,0] = a_X1[1] - a_X4[1]
     mat[1,1] = a_X2[1] - a_X4[1]
     mat[1,2] = a_X3[1] - a_X4[1]

     mat[2,0] = a_X1[2] - a_X4[2]
     mat[2,1] = a_X2[2] - a_X4[2]
     mat[2,2] = a_X3[2] - a_X4[2]
     
     F_inv = np.linalg.inv(mat)
     V_out = np.zeros( (3,3))
     J = np.vstack( (np.identity(3), np.array([[-1,-1,-1]])) )
     for i in range(3):
         Vi = [ a_V1[i], a_V2[i], a_V3[i], a_V4[i]]
         V_out[i, :] = np.dot(Vi, np.dot(J,F_inv)) 

     return V_out
     

         
     
     
