import numpy as np
from scipy.linalg import polar
import time

from elements.CElement import CElement

class CTetrahedron(CElement):

    def __init__(self, ID):
        self.SetID(ID)

    def SetVerticesCoordinatesAndGradient(self):

        self.verticesCoords    = np.zeros((4,3))
        self.verticesGradients = np.zeros((4,3))

        for i, vertex in enumerate(self.vertices):

            self.verticesCoords[i,:]    = vertex.GetCoordinates()
            self.verticesGradients[i,:] = vertex.GetGradient()

        # A = np.array([[x1, y1, z1, 1],
        #               [x2, y2, z2, 1],
        #               [x3, y3, z3, 1],
        #               [x4, y4, z4, 1]])

        self.A = np.hstack((self.verticesCoords, np.ones((4,1))))

    def ComputeVolume(self):
        self.volume = (1/6) * abs(np.linalg.det(self.A))
    
    def ComputePatchVolume(self):
        self.patchVolume = 0.0
        for elem in self.patchElements:
            self.patchVolume += elem.volume
    
    def ComputeGradient(self):
        
        # f_gradient_x = np.array([gradient_x_1, gradient_x_2, gradient_x_3, gradient_x_4])
        # f_gradient_y = np.array([gradient_y_1, gradient_y_2, gradient_y_3, gradient_y_4])
        # f_gradient_z = np.array([gradient_z_1, gradient_z_2, gradient_z_3, gradient_z_4])

        # coeff_gradient_x = np.linalg.solve(self.A, f_gradient_x)
        # coeff_gradient_y = np.linalg.solve(self.A, f_gradient_y)
        # coeff_gradient_z = np.linalg.solve(self.A, f_gradient_z)
        
        # self.gradient = np.array([coeff_gradient_x, coeff_gradient_y, coeff_gradient_z])
        self.gradient = np.transpose(np.linalg.solve(self.A, self.verticesGradients))

    def ComputeLocalErrorContribution(self):

        XA = (self.verticesCoords[0,0]+self.verticesCoords[1,0]+self.verticesCoords[2,0])/3 
        YA = (self.verticesCoords[0,1]+self.verticesCoords[1,1]+self.verticesCoords[2,1])/3 
        ZA = (self.verticesCoords[0,2]+self.verticesCoords[1,2]+self.verticesCoords[2,2])/3 

        XB = (self.verticesCoords[1,0]+self.verticesCoords[2,0]+self.verticesCoords[3,0])/3 
        YB = (self.verticesCoords[1,1]+self.verticesCoords[2,1]+self.verticesCoords[3,1])/3 
        ZB = (self.verticesCoords[1,2]+self.verticesCoords[2,2]+self.verticesCoords[3,2])/3 

        XC = (self.verticesCoords[2,0]+self.verticesCoords[3,0]+self.verticesCoords[0,0])/3 
        YC = (self.verticesCoords[2,1]+self.verticesCoords[3,1]+self.verticesCoords[0,1])/3 
        ZC = (self.verticesCoords[2,2]+self.verticesCoords[3,2]+self.verticesCoords[0,2])/3 

        XD = (self.verticesCoords[3,0]+self.verticesCoords[0,0]+self.verticesCoords[1,0])/3 
        YD = (self.verticesCoords[3,1]+self.verticesCoords[0,1]+self.verticesCoords[1,1])/3 
        ZD = (self.verticesCoords[3,2]+self.verticesCoords[0,2]+self.verticesCoords[1,2])/3 

        baryMat = np.transpose(np.array([[XA, YA, ZA, 1], 
                                         [XB, YB, ZB, 1],
                                         [XC, YC, ZC, 1],
                                         [XD, YD, ZD, 1]]))
        
        # gradEvalBary is a (3, 4) np.array storing the gradient 
        # evaluated at quadrature points along the column
        gradEvalBary = self.gradient @ baryMat

        self.localErrorContribution = self.volume/4 * np.dot(gradEvalBary, gradEvalBary.T)

    def ComputeLambdak(self):

        Mk = np.zeros((3,3))

        Mk[0,0] = 1/4 * np.sqrt(6) * (self.verticesCoords[1,0]-self.verticesCoords[0,0])
        Mk[1,0] = 1/4 * np.sqrt(6) * (self.verticesCoords[1,1]-self.verticesCoords[0,1])
        Mk[2,0] = 1/4 * np.sqrt(6) * (self.verticesCoords[1,2]-self.verticesCoords[0,2])

        Mk[0,1] = 1/4 * np.sqrt(2) * (2*self.verticesCoords[2,0]-self.verticesCoords[0,0]-self.verticesCoords[1,0])
        Mk[1,1] = 1/4 * np.sqrt(2) * (2*self.verticesCoords[2,1]-self.verticesCoords[0,1]-self.verticesCoords[1,1])
        Mk[2,1] = 1/4 * np.sqrt(2) * (2*self.verticesCoords[2,2]-self.verticesCoords[0,2]-self.verticesCoords[1,2])

        Mk[0,2] = 1/4 * (3*self.verticesCoords[3,0]-self.verticesCoords[0,0]-self.verticesCoords[1,0]-self.verticesCoords[2,0])
        Mk[1,2] = 1/4 * (3*self.verticesCoords[3,1]-self.verticesCoords[0,1]-self.verticesCoords[1,1]-self.verticesCoords[2,1])
        Mk[2,2] = 1/4 * (3*self.verticesCoords[3,2]-self.verticesCoords[0,2]-self.verticesCoords[1,2]-self.verticesCoords[2,2])

        self.lambda_k = np.linalg.svd(Mk, compute_uv=False)

    def ComputeMetric(self, toll=1.0, diam=1.0, card=1000):

        patchG = np.zeros((3,3))

        for elem in self.patchElements:
            factor = (elem.volume / self.patchVolume - 1)
            patchG += factor * factor * elem.localErrorContribution

        eigenValRefG, eigenVecRefG = np.linalg.eigh(patchG / self.patchVolume)

        referencePatchVolume = self.patchVolume / np.prod(self.lambda_k)

        valGmin = diam**(-3) * toll**2 / (3 * card * referencePatchVolume)

        valG1 = max(eigenValRefG[2], valGmin)
        valG2 = max(eigenValRefG[1], valGmin)
        valG3 = max(eigenValRefG[0], valGmin)

        factor = (toll**2 / (3 * card * referencePatchVolume))**(1/3) * (valG1 * valG2 * valG3)**(1/18)

        lambda_1 = factor * valG3**(-0.5)
        lambda_2 = factor * valG2**(-0.5)
        lambda_3 = factor * valG1**(-0.5)

        lambdaNewm2 = np.diag(np.array([lambda_1, lambda_2, lambda_3])**(-2))

        self.metric = eigenVecRefG @ lambdaNewm2 @ np.transpose(eigenVecRefG)

        return (valG2 == valGmin or valG1 == valGmin or valG3 == valGmin)
