import numpy as np

class CVertex():

    def __init__(self, ID, x, y, z=None):
        self.SetID(ID)
        if z is not None:
            self.SetCoordinates(x, y, z=z)
            self.dim = 3
        else:
            self.SetCoordinates(x, y, z=0.0)
            self.dim = 2
    
    def SetID(self, ID):
        if hasattr(self, 'ID'):
            raise KeyError('Already specified ID!')
        self.ID = int(ID)
    
    def SetCoordinates(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z

    def SetSolution(self, solution, iSensor=0):
        self.solution[iSensor] = solution

    def SetGradient(self, gradient, iSensor=0):
        self.gradient[:,iSensor] = gradient

    def SetVerticesNeighboursID(self, vertexIDs):
        self.verticesNeighboursIDs = vertexIDs

    def SetElementsNeighboursID(self, elementsIDs):
        self.elementsNeighboursIDs = elementsIDs

    def SetSensors(self, sensors):
        self.sensors = sensors
        self.nSensors = len(sensors)
        self.solution = np.zeros(self.nSensors)
        self.gradient = np.zeros((self.dim, self.nSensors))

    def GetID(self):
        return self.ID
    
    def GetCoordinates(self):
        return self.x, self.y, self.z  
    
    def GetSolution(self, iSensor=0):
        return self.solution[iSensor]
    
    def GetGradient(self, iSensor=0):
        return self.gradient[:,iSensor]
    
    def GetVerticesNeighboursID(self):
        return self.verticesNeighboursIDs 

    def GetElementsNeighboursID(self):
        return self.elementsNeighboursIDs 
    
    def GetMetric(self):
        return self.metric
    
    def ComputeMetric(self, mesh):
        elementsTotalVolume = 0.0
        elementsMetricSum = np.zeros((2,2)) if self.dim == 2 else np.zeros((3,3))
        for id in self.elementsNeighboursIDs:
            element = mesh.GetElement(id)
            volume = element.GetVolume()
            metric = element.GetFinalMetric()
            elementsTotalVolume += volume
            elementsMetricSum += volume * metric
        self.metric = elementsMetricSum / elementsTotalVolume






