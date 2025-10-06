import numpy as np
from scipy.linalg import solve, eig, inv
class CElement():

    def __init__(self, ID):
        self.SetID(ID)

    def SetID(self, ID):
        if hasattr(self, 'ID'):
            raise KeyError('Already specified ID!')
        self.ID = int(ID)

    def SetSensors(self, sensors):
        self.sensors  = sensors
        self.nSensors = len(sensors)
        self.gradient = np.zeros((self.dim, self.dim+1, self.nSensors))     # remember, we store the coefficients of the gradient linear interpolation over the element
        self.metric   = np.zeros((self.dim, self.dim, self.nSensors))
        self.localAnisoError = np.zeros(self.nSensors)
        self.localErrorContribution = np.zeros((self.dim, self.dim, self.nSensors))

    def SetVerticesID(self, listID):
        self.verticesID = listID

    def SetVertices(self, mesh):
        self.vertices = []
        for id in self.verticesID:
            self.vertices.append(mesh.GetVertex(id))
    
    def SetPatchID(self):
        patchElementsID = []
        for vertex in self.vertices:
            patchElementsID += vertex.GetElementsNeighboursID()
        self.patchElementsID = list(set(patchElementsID))

    def SetPatchElements(self, mesh):
        self.patchElements = []
        for id in self.patchElementsID:
            self.patchElements.append(mesh.GetElement(id))

    def SetMetric(self, metric, iSensor=0):
        self.metric[:,:,iSensor] = metric

    def SetAnisotropicError(self, localAnisoError, iSensor=0):
        self.localAnisoError[iSensor] = localAnisoError

    def GetID(self):
        return self.ID

    def GetVerticesID(self):
        return self.verticesID
    
    def GetVertices(self):
        return self.vertices
    
    def GetPatchElements(self):
        return self.patchElements
    
    def GetVolume(self):
        return self.volume
    
    def GetPatchVolume(self):
        return self.patchVolume
        
    def GetGradient(self, iSensor=0):
        return self.gradient[:,:,iSensor]
    
    def GetMetric(self, iSensor=0):
        return self.metric[:,:,iSensor]
    
    def GetAnisotropicError(self, iSensor=0):
        return self.localAnisoError
    
    def GetFinalMetric(self):
        return self.finalMetric
    
    def IntersectMetric(self):

        if self.nSensors == 1:
            self.finalMetric = self.metric[:,:,0]
            return

        elif self.nSensors == 2:

            # N = M1^(-1) M2
            M1 = self.metric[:,:,0]
            M2 = self.metric[:,:,1]
            N = solve(M1, M2, overwrite_a=False, check_finite=False, assume_a='pos')

            # P = (e1 e2 e3) i.e., eigenvectors of N
            P = np.real(eig(N, left=False, right=True)[1]) 

            # eigenvalues of M1, M2 in basis P
            lam = np.diag(P.T @ M1 @ P)
            mu  = np.diag(P.T @ M2 @ P)

            # linear interpolation
            h1 = np.real(1 / np.sqrt(lam))
            h2 = np.real(1 / np.sqrt(mu ))
            t = 0.5
            H = ((1-t) * h1 + t * h2)**(-2)

            # min max 
            # H = np.amax(np.hstack((lam[:,np.newaxis], mu[:,np.newaxis])), axis=1)

            # self.finalMetric = P @ np.diag(H) @ P.T
            self.finalMetric = P @ np.diag(H) @ P.T

        else:
            raise ValueError('Metric intersection not implemented for more than 2 sensors at the moment.')
        
    

