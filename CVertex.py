import numpy as np

class CVertex():

    def __init__(self, ID, x, y):
        self.SetID(ID)
        self.SetCoordinates(x, y)
    
    def SetID(self, ID):
        if hasattr(self, 'ID'):
            raise KeyError('Already specified ID!')
        self.ID = int(ID)
    
    def SetCoordinates(self, x, y):
        self.x = x
        self.y = y

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
        elementsTotalArea = 0.0
        elementsMetricSum = np.zeros((2,2)) if mesh.dim == 2 else np.zeros((3,3))
        for id in self.elementsNeighboursIDs:
            element = mesh.GetElement(id)
            area = element.GetArea()
            metric = element.GetMetric()
            elementsTotalArea += area
            elementsMetricSum += area * metric
        self.metric = elementsMetricSum / elementsTotalArea






