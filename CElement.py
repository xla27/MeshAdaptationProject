import numpy as np

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
            patchElementsID.append(vertex.GetElementsNeighboursID())
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
    
    def GetArea(self):
        return self.area
    
    def GetPatchArea(self):
        return self.patchArea
    
    def GetGradient(self):
        return self.gradient

    def ComputeArea(self):
        return NotImplementedError
    
    def ComputePatchArea(self):
        self.patchArea = 0.0
        for elem in self.patchElements:
            self.patchArea += elem.GetArea()
    
    def ComputeGradient(self):
        gradientsVertices = []
        for vertex in self.vertices:
            gradientsVertices.append(vertex.GetGradient())

        # Interpolazione
        #self.gradient = [gradient_x, gradient_y]
        return NotImplementedError

    def ComputePatchAveragedGradient(self):
        patchAvgGradient = 0.0
        for elem in self.patchElements:
            patchAvgGradient += elem.GetArea() * elem.GetGradient()
        self.patchAvgGradient = patchAvgGradient / self.GetPatchArea()






