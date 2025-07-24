import numpy as np
from scipy.linalg import polar, eig

from elements.CElement import CElement

class CTriangle(CElement):

    def __init__(self, ID):
        self.SetID(ID)

    def SetVerticesCoordinatesAndGradient(self):

        self.verticesCoords    = np.zeros((3,2))
        self.verticesGradients = np.zeros((3,2))

        for i, vertex in enumerate(self.vertices):

            self.verticesCoords[i,:]    = vertex.GetCoordinates()[:-1]
            self.verticesGradients[i,:] = vertex.GetGradient()

        # A = np.array([[x1, y1, 1],
        #               [x2, y2, 1],
        #               [x3, y3, 1]])

        self.A = np.hstack((self.verticesCoords, np.ones((3,1))))

    def ComputeVolume(self):
        self.volume = (1/2) * abs(np.linalg.det(self.A))
    
    def ComputePatchVolume(self):
        self.patchVolume = 0.0
        for elem in self.patchElements:
            self.patchVolume += elem.GetVolume()
    
    def ComputeGradient(self):
        
        # f_gradient_x = np.array([gradient_x_1, gradient_x_2, gradient_x_3])
        # f_gradient_y = np.array([gradient_y_1, gradient_y_2, gradient_y_3])

        # coeff_gradient_x = np.linalg.solve(self.A, f_gradient_x)
        # coeff_gradient_y = np.linalg.solve(self.A, f_gradient_y)

        # self.gradient = np.array([coeff_gradient_x, coeff_gradient_y])
        self.gradient = np.transpose(np.linalg.solve(self.A, self.verticesGradients))

    def ComputeLocalErrorContribution(self):

        XA = (self.verticesCoords[0,0]+self.verticesCoords[1,0])/2 
        YA = (self.verticesCoords[0,1]+self.verticesCoords[1,1])/2 

        XB = (self.verticesCoords[1,0]+self.verticesCoords[2,0])/2 
        YB = (self.verticesCoords[1,1]+self.verticesCoords[2,1])/2 

        XC = (self.verticesCoords[2,0]+self.verticesCoords[1,0])/2 
        YC = (self.verticesCoords[2,1]+self.verticesCoords[1,1])/2 

        baryMat = np.transpose(np.array([[XA, YA, 1], 
                                         [XB, YB, 1],
                                         [XC, YC, 1]]))
        
        # gradEvalBary is a (2, 3) np.array storing the gradient 
        # evaluated at quadrature points along the column
        gradEvalBary = self.gradient @ baryMat

        self.localErrorContribution = self.volume/3 * np.dot(gradEvalBary, gradEvalBary.T)

    def ComputeLambdak(self):

        Mk = np.zeros((2,2))
        tk = np.zeros(2)
        tk[0] = 1/3 * (self.verticesCoords[0,0]+self.verticesCoords[1,0]+self.verticesCoords[2,0])
        tk[1] = 1/3 * (self.verticesCoords[0,1]+self.verticesCoords[1,1]+self.verticesCoords[2,1])
        Mk[0,1] = 2 * tk[0] - self.verticesCoords[0,0] - self.verticesCoords[1,0]
        Mk[1,1] = 2 * tk[1] - self.verticesCoords[0,1] - self.verticesCoords[1,1]
        Mk[0,0] = 1/np.sqrt(3) * (2 * tk[0] - Mk[0,1] - 2 * self.verticesCoords[0,0])
        Mk[1,0] = 1/np.sqrt(3) * (2 * tk[1] - Mk[1,1] - 2 * self.verticesCoords[0,1])

        self.lambda_k = np.linalg.svd(Mk, compute_uv=False)

    def ComputeMetric(self, toll=1.0, diam=1.0, card=1000):

        patchG = np.zeros((2,2))

        for elem in self.patchElements:
            factor = (elem.volume / self.patchVolume - 1)
            patchG += factor * factor * elem.localErrorContribution

        eigenValRefG, eigenVecRefG = np.linalg.eigh(patchG / self.patchVolume)

        referencePatchVolume = self.patchVolume / np.prod(self.lambda_k)

        valGmin = diam**(-2) * toll**2 / (2 * card * referencePatchVolume)

        valG1 = max(eigenValRefG[1], valGmin)
        valG2 = max(eigenValRefG[0], valGmin)

        factor = (toll**2 / (2 * card * referencePatchVolume))**0.5

        lambda_1 = factor * valG2**(-0.5)
        lambda_2 = factor * valG1**(-0.5)

        lambdaNewm2 = np.diag(np.array([lambda_1, lambda_2])**(-2))

        self.metric = eigenVecRefG @ lambdaNewm2 @ np.transpose(eigenVecRefG)

        return valG2 == valGmin or valG1 == valGmin
