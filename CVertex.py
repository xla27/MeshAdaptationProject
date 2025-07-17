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

    def SetVerticesNeighboursID(self, vertex_IDs):
        self.vertices_neighbours_IDs = vertex_IDs

    def SetElementsNeighboursID(self, elements_IDs):
        self.elements_neighbours_IDs = elements_IDs

    def GetID(self):
        return self.ID
    
    def GetCoordinates(self):
        return self.x, self.y
    
    def GetSolution(self):
        return self.solution
    
    def GetGradient(self):
        return self.gradient
    
    def GetVerticesNeighboursID(self):
        return self.vertices_neighbours_IDs 

    def GetElementsNeighboursID(self):
        return self.elements_neighbours_IDs 




