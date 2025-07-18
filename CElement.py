import numpy as np
from scipy.linalg import polar, eig

class CElement():

    def __init__(self, ID):
        self.SetID(ID)

    def SetID(self, ID):
        if hasattr(self, 'ID'):
            raise KeyError('Already specified ID!')
        self.ID = int(ID)

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

    def GetID(self):
        return self.ID

    def GetVerticesID(self):
        return self.verticesID
    
    def GetVertices(self):
        return self.vertices
    
    def GetPatchElements(self):
        return self.patchElements
    
    def GetArea(self):
        if not hasattr(self, 'area'):
            self.ComputeArea()
        return self.area
    
    def GetPatchArea(self):
        return self.patchArea
        
    def GetGradient(self):
        if not hasattr(self, 'gradient'):
            self.ComputeGradient()
        return self.gradient
    
    def GetMetric(self):
        if not hasattr(self, 'metric'):
            self.ComputeMetric()
        return self.metric

    def ComputeArea(self):
        x1,y1 = self.vertices[0].GetCoordinates()
        x2,y2 = self.vertices[1].GetCoordinates()
        x3,y3 = self.vertices[2].GetCoordinates()
        self.area = (1/2) * abs(x1*(y2 - y3) - x2*(y1 - y3) + x3*(y1 - y2))
    
    def ComputePatchArea(self):
        self.patchArea = 0.0
        for elem in self.patchElements:
            self.patchArea += elem.GetArea()
    
    def ComputeGradient(self):
        x1,y1 = self.vertices[0].GetCoordinates()
        x2,y2 = self.vertices[1].GetCoordinates()
        x3,y3 = self.vertices[2].GetCoordinates()

        gradient_x_1,gradient_y_1 = self.vertices[0].GetGradient()
        gradient_x_2,gradient_y_2 = self.vertices[1].GetGradient()
        gradient_x_3,gradient_y_3 = self.vertices[2].GetGradient()

        A = np.array([[x1, y1, 1],
                      [x2, y2, 1],
                      [x3, y3, 1]])
        
        f_gradient_x = np.array([gradient_x_1, gradient_x_2, gradient_x_3])
        f_gradient_y = np.array([gradient_y_1, gradient_y_2, gradient_y_3])

        coeff_gradient_x = np.linalg.solve(A, f_gradient_x)
        coeff_gradient_y = np.linalg.solve(A, f_gradient_y)

        a_gradient_x, b_gradient_x, c_gradient_x = coeff_gradient_x
        a_gradient_y, b_gradient_y, c_gradient_y = coeff_gradient_y

        def gradient_x(x,y):
            return a_gradient_x * x + b_gradient_x * y + c_gradient_x
        def gradient_y(x,y):
            return a_gradient_y * x + b_gradient_y * y + c_gradient_y
        
        self.gradient = [gradient_x, gradient_y]

    def ComputePatchAveragedGradient(self):

        def patchAvgGradient_x(x,y):
            patchAvgGradient_x = 0.0
            borderCount = 0.0
            for elem in self.patchElements:
                inside, border = PointInTriangle([x,y], 
                                                    elem.GetVertices()[0].GetCoordinates(),
                                                    elem.GetVertices()[1].GetCoordinates(), 
                                                    elem.GetVertices()[2].GetCoordinates())
                if border:
                    elem_border = elem 
                    borderCount += 1
                patchAvgGradient_x += elem.GetArea() * elem.GetGradient()[0](x,y) * inside 
            if borderCount:
                patchAvgGradient_x = elem_border.GetArea() * elem_border.GetGradient()[0](x,y)
            return patchAvgGradient_x / self.GetPatchArea() 
        
        def patchAvgGradient_y(x,y):
            patchAvgGradient_y = 0.0
            borderCount = 0.0
            for elem in self.patchElements:
                inside, border = PointInTriangle([x,y], 
                                                    elem.GetVertices()[0].GetCoordinates(),
                                                    elem.GetVertices()[1].GetCoordinates(), 
                                                    elem.GetVertices()[2].GetCoordinates())
                if border:
                    elem_border = elem 
                    borderCount += 1
                patchAvgGradient_y += elem.GetArea() * elem.GetGradient()[1](x,y) * inside 
            if borderCount:
                patchAvgGradient_y = elem_border.GetArea() * elem_border.GetGradient()[1](x,y)
            return patchAvgGradient_y / self.GetPatchArea() 

        self.patchAvgGradient = [patchAvgGradient_x, patchAvgGradient_y]

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
        Lk, _ = eig(Bk)
        Lk = np.real(Lk)
        # r1_k = Rk[:,np.argmax(Lk)];
        lambda1_k = np.amax(Lk)
        # r2_k = Rk[:,np.argmin(Lk)];
        lambda2_k = np.amin(Lk)
        lambda_k = [lambda1_k,lambda2_k]
        self.lambda_k = lambda_k

    def ComputeMetric(self, params):

        patchG = np.zeros(3)
        for elem in self.patchElements:
            x1,y1 = elem.vertices[0].GetCoordinates()
            x2,y2 = elem.vertices[1].GetCoordinates()
            x3,y3 = elem.vertices[2].GetCoordinates()

            E_x = lambda x,y: self.patchAvgGradient[0](x,y) - elem.GetGradient()[0](x,y)
            E_y = lambda x,y: self.patchAvgGradient[1](x,y) - elem.GetGradient()[1](x,y)

            XA = (x1+x2)/2; YA = (y1+y2)/2
            XB = (x2+x3)/2; YB = (y2+y3)/2
            XC = (x3+x1)/2; YC = (y3+y1)/2

            patchG[0] += elem.GetArea()/3 * (
                E_x(XA,YA) * E_x(XA,YA) +
                E_x(XB,YB) * E_x(XB,YB) +
                E_x(XC,YC) * E_x(XC,YC) 
            )

            patchG[1] += elem.GetArea()/3 * (
                E_x(XA,YA) * E_y(XA,YA) +
                E_x(XB,YB) * E_y(XB,YB) +
                E_x(XC,YC) * E_y(XC,YC) 
            )

            patchG[2] += elem.GetArea()/3 * (
                E_y(XA,YA) * E_y(XA,YA) +
                E_y(XB,YB) * E_y(XB,YB) +
                E_y(XC,YC) * E_y(XC,YC) 
            )

        referencePatchG = np.zeros((2,2))
        referencePatchG[0,0] = patchG[0] / self.GetPatchArea()
        referencePatchG[0,1] = patchG[1] / self.GetPatchArea()
        referencePatchG[1,1] = patchG[2] / self.GetPatchArea()
        referencePatchG[1,0] = referencePatchG[0,1]

        eigenValRefG, eigenVecRefG = eig(referencePatchG)
        idx = eigenValRefG.argsort()[::-1]   
        eigenValRefG = eigenValRefG[idx]
        eigenVecRefG = eigenVecRefG[:,idx]

        valG1 = np.real(eigenValRefG[0]); vecG1 = eigenVecRefG[:,0]
        valG2 = np.real(eigenValRefG[1]); vecG2 = eigenVecRefG[:,1]

        toll = params['toll']
        card = params['card']

        self.ComputeLambdak()
        referencePatchArea = self.GetPatchArea() / (self.lambda_k[0] * self.lambda_k[1])

        valGmin = np.sqrt(2)**(-2) * toll**2 / (2 * card * referencePatchArea)

        valG1 = max(valG1, valGmin)
        valG2 = max(valG2, valGmin)

        factor = (toll**2 / (2 * card * referencePatchArea))**0.5

        lambda_1 = factor * valG2**(-0.5)
        lambda_2 = factor * valG1**(-0.5)

        lambdaNew = np.diag([lambda_1, lambda_2])
        RkNewT = np.hstack((vecG2[:,np.newaxis], vecG1[:,np.newaxis]))

        self.metric = RkNewT @ lambdaNew**(-2) @ RkNewT.T


def PointInTriangle(P, A, B, C):
    def dot(u, v): return u[0]*v[0] + u[1]*v[1]

    v0 = (C[0] - A[0], C[1] - A[1])
    v1 = (B[0] - A[0], B[1] - A[1])
    v2 = (P[0] - A[0], P[1] - A[1])

    d00 = dot(v0, v0)
    d01 = dot(v0, v1)
    d11 = dot(v1, v1)
    d20 = dot(v2, v0)
    d21 = dot(v2, v1)

    denom = d00 * d11 - d01 * d01
    v = (d11 * d20 - d01 * d21) / denom
    w = (d00 * d21 - d01 * d20) / denom
    u = 1 - v - w

    if (u > 0) and (v > 0) and (w > 0):
        inside = 1.0
    else:
        inside = 0.0

    if (u == 0) or (v == 0) or (w == 0):
        if (u < 0) or (v < 0) or (w < 0):
            border = False
        else:
            border = True
    else:
        border = False

    return inside, border







