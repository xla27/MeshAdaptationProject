import numpy as np

class CVertex():

    def __init__(self, ID, x, y, z=None):
        self.SetID(ID)
        if z:
            self.SetCoordinates(x, y, z=z)
        else:
            self.SetCoordinates(x, y)
    
    def SetID(self, ID):
        if hasattr(self, 'ID'):
            raise KeyError('Already specified ID!')
        self.ID = int(ID)
    
    def SetCoordinates(self, x, y, z=None):
        self.x = x
        self.y = y
        if z:
            self.z = z

    def SetSolution(self, solution):
        self.solution = solution

    def SetGradient(self, gradient):
        self.gradient = gradient

    def SetVerticesNeighboursID(self, vertexIDs):
        self.verticesNeighboursIDs = vertexIDs

    def SetElementsNeighboursID(self, elementsIDs):
        self.elementsNeighboursIDs = elementsIDs

    def GetID(self):
        return self.ID
    
    def GetCoordinates(self):
        if hasattr(self,'z'):
            return self.x, self.y, self.z
        else:   
            return self.x, self.y
    
    def GetSolution(self):
        return self.solution
    
    def GetGradient(self):
        return self.gradient
    
    def GetVerticesNeighboursID(self):
        return self.verticesNeighboursIDs 

    def GetElementsNeighboursID(self):
        return self.elementsNeighboursIDs 
    
    def GetMetric(self):
        return self.metric
    
    def ComputeMetric(self, mesh):
        elementsTotalVolume = 0.0
        elementsMetricSum = np.zeros((2,2)) if mesh.GetDim() == 2 else np.zeros((3,3))
        for id in self.elementsNeighboursIDs:
            element = mesh.GetElement(id)
            volume = element.GetVolume()
            metric = element.GetMetric()
            elementsTotalVolume += volume
            elementsMetricSum += volume * metric
        self.metric = elementsMetricSum / elementsTotalVolume






