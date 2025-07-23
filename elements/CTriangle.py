import numpy as np
from scipy.linalg import polar, eig

from elements.CElement import CElement

class CTriangle(CElement):

    def __init__(self, ID):
        self.SetID(ID)

    def ComputeVolume(self):
        x1,y1 = self.vertices[0].GetCoordinates()
        x2,y2 = self.vertices[1].GetCoordinates()
        x3,y3 = self.vertices[2].GetCoordinates()
        self.volume = (1/2) * abs(x1*(y2 - y3) - x2*(y1 - y3) + x3*(y1 - y2))
    
    def ComputePatchVolume(self):
        self.patchVolume = 0.0
        for elem in self.patchElements:
            self.patchVolume += elem.GetVolume()
    
    def ComputeGradient(self):
        x1,y1 = self.vertices[0].GetCoordinates()
        x2,y2 = self.vertices[1].GetCoordinates()
        x3,y3 = self.vertices[2].GetCoordinates()

        gradient_x_1, gradient_y_1 = self.vertices[0].GetGradient()
        gradient_x_2, gradient_y_2 = self.vertices[1].GetGradient()
        gradient_x_3, gradient_y_3 = self.vertices[2].GetGradient()

        A = np.array([[x1, y1, 1],
                      [x2, y2, 1],
                      [x3, y3, 1]])
        
        f_gradient_x = np.array([gradient_x_1, gradient_x_2, gradient_x_3])
        f_gradient_y = np.array([gradient_y_1, gradient_y_2, gradient_y_3])

        coeff_gradient_x = np.linalg.solve(A, f_gradient_x)
        coeff_gradient_y = np.linalg.solve(A, f_gradient_y)

        self.gradient = np.array([coeff_gradient_x, coeff_gradient_y])

    def ComputeLambdak(self):
        x1,y1 = self.vertices[0].GetCoordinates()
        x2,y2 = self.vertices[1].GetCoordinates()
        x3,y3 = self.vertices[2].GetCoordinates()

        Mk = np.zeros((2,2))
        tk = np.zeros(2)
        tk[0] = 1/3 * (x1 + x2 + x3)
        tk[1] = 1/3 * (y1 + y2 + y3)
        Mk[0,1] = 2 * tk[0] - x1 - x2
        Mk[1,1] = 2 * tk[1] - y1 - y2
        Mk[0,0] = 1/np.sqrt(3) * (2 * tk[0] - Mk[0,1] - 2 * x1)
        Mk[1,0] = 1/np.sqrt(3) * (2 * tk[1] - Mk[1,1] - 2 * y1)

        _, Bk = polar(Mk, side='left')
        Lk, _ = np.linalg.eigh(Bk)
        self.lambda_k = [Lk[1], Lk[0]]

    def ComputeMetric(self, toll=1.0, diam=1.0, card=1000):

        patchG = np.zeros(3)

        for elem in self.patchElements:
            x1,y1 = elem.vertices[0].GetCoordinates()
            x2,y2 = elem.vertices[1].GetCoordinates()
            x3,y3 = elem.vertices[2].GetCoordinates()

            XA = (x1+x2)/2; YA = (y1+y2)/2
            XB = (x2+x3)/2; YB = (y2+y3)/2
            XC = (x3+x1)/2; YC = (y3+y1)/2

            E_xA, E_yA = ComputeError(elem, XA, YA, self.GetPatchVolume())
            E_xB, E_yB = ComputeError(elem, XB, YB, self.GetPatchVolume())
            E_xC, E_yC = ComputeError(elem, XC, YC, self.GetPatchVolume())

            patchG[0] += elem.GetVolume()/3 * (E_xA * E_xA + E_xB * E_xB + E_xC * E_xC)

            patchG[1] += elem.GetVolume()/3 * (E_xA * E_yA + E_xB * E_yB + E_xC * E_yC)

            patchG[2] += elem.GetVolume()/3 * (E_yA * E_yA + E_yB * E_yB + E_yC * E_yC)

        referencePatchG = np.zeros((2,2))
        referencePatchG[0,0] = patchG[0] 
        referencePatchG[0,1] = patchG[1] 
        referencePatchG[1,1] = patchG[2] 
        referencePatchG[1,0] = referencePatchG[0,1]
        referencePatchG /= self.GetPatchVolume()

        eigenValRefG, eigenVecRefG = np.linalg.eigh(referencePatchG)

        self.ComputeLambdak()
        referencePatchVolume = self.GetPatchVolume() / (self.lambda_k[0] * self.lambda_k[1])

        valGmin = diam**(-2) * toll**2 / (2 * card * referencePatchVolume)

        valG1 = max(eigenValRefG[1], valGmin)
        valG2 = max(eigenValRefG[0], valGmin)

        limited = 1 if valG2 == valGmin or valG1 == valGmin else 0

        factor = (toll**2 / (2 * card * referencePatchVolume))**0.5

        lambda_1 = factor * valG2**(-0.5)
        lambda_2 = factor * valG1**(-0.5)

        lambdaNewm2 = np.diag(np.array([lambda_1, lambda_2])**(-2))

        # vecG1 = eigenVecRefG[:,1]
        # vecG2 = eigenVecRefG[:,0]
        # RkNewT = np.hstack((vecG2[:,np.newaxis], vecG1[:,np.newaxis]))
        # RkNewT = eigenVecRefG    # equivalent form)

        self.metric = eigenVecRefG @ lambdaNewm2 @ np.transpose(eigenVecRefG)

        return limited


def ComputeError(elem, x, y, patchVolume):
    grad_coeffs = elem.GetGradient()
    gradient = np.dot(grad_coeffs, [x,y,1])
    factor   = (elem.GetVolume() / patchVolume - 1)
    return factor * gradient