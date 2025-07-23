import numpy as np
from scipy.linalg import polar
import time

from elements.CElement import CElement

class CTetrahedron(CElement):

    def __init__(self, ID):
        self.SetID(ID)

    def ComputeVolume(self):
        x1,y1,z1 = self.vertices[0].GetCoordinates()
        x2,y2,z2 = self.vertices[1].GetCoordinates()
        x3,y3,z3 = self.vertices[2].GetCoordinates()
        x4,y4,z4 = self.vertices[3].GetCoordinates()
        A = np.array([[x1, y1, z1, 1],
                    [x2, y2, z2, 1],
                    [x3, y3, z3, 1],
                    [x4, y4, z4, 1]])
        self.volume = (1/6) * abs(np.linalg.det(A))
    
    def ComputePatchVolume(self):
        self.patchVolume = 0.0
        for elem in self.patchElements:
            self.patchVolume += elem.GetVolume()
    
    def ComputeGradient(self):
        x1,y1,z1 = self.vertices[0].GetCoordinates()
        x2,y2,z2 = self.vertices[1].GetCoordinates()
        x3,y3,z3 = self.vertices[2].GetCoordinates()
        x4,y4,z4 = self.vertices[3].GetCoordinates()

        gradient_x_1, gradient_y_1, gradient_z_1 = self.vertices[0].GetGradient()
        gradient_x_2, gradient_y_2, gradient_z_2 = self.vertices[1].GetGradient()
        gradient_x_3, gradient_y_3, gradient_z_3 = self.vertices[2].GetGradient()
        gradient_x_4, gradient_y_4, gradient_z_4 = self.vertices[3].GetGradient()

        A = np.array([[x1, y1, z1, 1],
                    [x2, y2, z2, 1],
                    [x3, y3, z3, 1],
                    [x4, y4, z4, 1]])
        
        f_gradient_x = np.array([gradient_x_1, gradient_x_2, gradient_x_3, gradient_x_4])
        f_gradient_y = np.array([gradient_y_1, gradient_y_2, gradient_y_3, gradient_y_4])
        f_gradient_z = np.array([gradient_z_1, gradient_z_2, gradient_z_3, gradient_z_4])

        coeff_gradient_x = np.linalg.solve(A, f_gradient_x)
        coeff_gradient_y = np.linalg.solve(A, f_gradient_y)
        coeff_gradient_z = np.linalg.solve(A, f_gradient_z)
        
        self.gradient = np.array([coeff_gradient_x, coeff_gradient_y, coeff_gradient_z])

    def ComputeLambdak(self):
        x1,y1,z1 = self.vertices[0].GetCoordinates()
        x2,y2,z2 = self.vertices[1].GetCoordinates()
        x3,y3,z3 = self.vertices[2].GetCoordinates()
        x4,y4,z4 = self.vertices[3].GetCoordinates()

        Mk = np.zeros((3,3))
        tk = np.zeros(3)

        tk[0] = 1/4 * (x1 + x2 + x3 + x4)
        tk[1] = 1/4 * (y1 + y2 + y3 + y4)
        tk[2] = 1/4 * (z1 + z2 + z3 + z4)

        Mk[0,0] = 1/4 * np.sqrt(6) * (x2 - x1)
        Mk[1,0] = 1/4 * np.sqrt(6) * (y2 - y1)
        Mk[2,0] = 1/4 * np.sqrt(6) * (z2 - z1)

        Mk[0,1] = 1/4 * np.sqrt(2) * (2*x3 - x1 - x2)
        Mk[1,1] = 1/4 * np.sqrt(2) * (2*y3 - y1 - y2)
        Mk[2,1] = 1/4 * np.sqrt(2) * (2*z3 - z1 - z2)

        Mk[0,2] = 1/4 * (3*x4 - x1 - x2 - x3)
        Mk[1,2] = 1/4 * (3*y4 - y1 - y2 - y3)
        Mk[2,2] = 1/4 * (3*z4 - z1 - z2 - z3)

        _, Bk = polar(Mk, side='left')
        Lk, _ = np.linalg.eigh(Bk)
        self.lambda_k = [Lk[2], Lk[1], Lk[0]]

    def ComputeMetric(self, toll=1.0, diam=1.0, card=1000):

        patchG = np.zeros(6)

        for elem in self.patchElements:
            x1,y1,z1 = elem.vertices[0].GetCoordinates()
            x2,y2,z2 = elem.vertices[1].GetCoordinates()
            x3,y3,z3 = elem.vertices[2].GetCoordinates()
            x4,y4,z4 = elem.vertices[3].GetCoordinates()

            XA = (x1+x2+x3)/3; YA = (y1+y2+y3)/3; ZA = (z1+z2+z3)/3
            XB = (x2+x3+x4)/3; YB = (y2+y3+y4)/3; ZB = (z2+z3+z4)/3
            XC = (x3+x4+x1)/3; YC = (y3+y4+y1)/3; ZC = (z3+z4+z1)/3
            XD = (x4+x1+x2)/3; YD = (y4+y1+y2)/3; ZD = (z4+z1+z2)/3

            E_xA, E_yA, E_zA = ComputeError(elem, XA, YA, ZA, self.GetPatchVolume())
            E_xB, E_yB, E_zB = ComputeError(elem, XB, YB, ZB, self.GetPatchVolume())
            E_xC, E_yC, E_zC = ComputeError(elem, XC, YC, ZC, self.GetPatchVolume())
            E_xD, E_yD, E_zD = ComputeError(elem, XD, YD, ZD, self.GetPatchVolume())

            # G_xx
            patchG[0] += elem.GetVolume()/4 * (E_xA * E_xA + E_xB * E_xB + E_xC * E_xC + E_xD * E_xD)

            # G_xy
            patchG[1] += elem.GetVolume()/4 * (E_xA * E_yA + E_xB * E_yB + E_xC * E_yC + E_xD * E_yD)

            # G_yy
            patchG[2] += elem.GetVolume()/4 * (E_yA * E_yA + E_yB * E_yB + E_yC * E_yC + E_yD * E_yD)

            # G_xz
            patchG[3] += elem.GetVolume()/4 * (E_xA * E_zA + E_xB * E_zB + E_xC * E_zC + E_xD * E_zD)

            # G_yz
            patchG[4] += elem.GetVolume()/4 * (E_yA * E_zA + E_yB * E_zB + E_yC * E_zC + E_yD * E_zD)

            # G_zz
            patchG[5] += elem.GetVolume()/4 * (E_zA * E_zA + E_zB * E_zB + E_zC * E_zC + E_zD * E_zD)

        referencePatchG = np.zeros((3,3))
        referencePatchG[0,0] = patchG[0] 
        referencePatchG[0,1] = patchG[1]; referencePatchG[1,0] = referencePatchG[0,1] 
        referencePatchG[1,1] = patchG[2] 
        referencePatchG[0,2] = patchG[3]; referencePatchG[2,0] = referencePatchG[0,2] 
        referencePatchG[1,2] = patchG[4]; referencePatchG[2,1] = referencePatchG[1,2] 
        referencePatchG[2,2] = patchG[5] 
        referencePatchG /= self.GetPatchVolume()

        eigenValRefG, eigenVecRefG = np.linalg.eigh(referencePatchG)

        self.ComputeLambdak()
        referencePatchVolume = self.GetPatchVolume() / np.prod(self.lambda_k)

        valGmin = diam**(-3) * toll**2 / (3 * card * referencePatchVolume)

        valG1 = max(eigenValRefG[2], valGmin)
        valG2 = max(eigenValRefG[1], valGmin)
        valG3 = max(eigenValRefG[0], valGmin)

        limited = 1 if valG2 == valGmin or valG1 == valGmin or valG3 == valGmin else 0

        factor = (toll**2 / (3 * card * referencePatchVolume))**(1/3) * (valG1 * valG2 * valG3)**(1/18)

        lambda_1 = factor * valG3**(-0.5)
        lambda_2 = factor * valG2**(-0.5)
        lambda_3 = factor * valG1**(-0.5)

        lambdaNewm2 = np.diag(np.array([lambda_1, lambda_2, lambda_3])**(-2))

        # vecG1 = eigenVecRefG[:,2]
        # vecG2 = eigenVecRefG[:,1]
        # vecG3 = eigenVecRefG[:,0]
        # RkNewT = np.hstack((vecG3[:,np.newaxis], vecG2[:,np.newaxis], vecG1[:,np.newaxis]))
        # RkNewT = eigenVecRefG    # equivalent form)

        self.metric = eigenVecRefG @ lambdaNewm2 @ np.transpose(eigenVecRefG)

        return limited


def ComputeError(elem, x, y, z, patchVolume):
    grad_coeffs = elem.GetGradient()
    gradient = np.dot(grad_coeffs, [x,y,z,1])
    factor   = (elem.GetVolume() / patchVolume - 1)
    return factor * gradient