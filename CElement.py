import numpy as np

class CElements():

    def __init__(self, ID):
        self.SetID(ID)

    def SetID(self, ID):
        if hasattr(self, 'ID'):
            raise KeyError('Already specified ID!')
        self.ID = int(ID)

    def SetVertices(self, vertex1, vertex2, vertex3):
        self.vertices = [vertex1, vertex2, vertex3]
    
    def SetPatchID(self):
        patch_elements_ID = []
        for vertex in self.vertices:
            patch_elements_ID.append(vertex.GetElementsNeighboursID())
        self.patch_elements_ID = list(set(patch_elements_ID))

    def SetPatchElements(self, mesh):
        self.patch_elements = []
        for id in self.patch_elements_ID:
            self.patch_elements.append(mesh.GetElement(id))
    
    def GetVertices(self):
        return self.vertices
    
    def GetArea(self):
        return self.area
    
    def GetPatchArea(self):
        return self.patch_area
    
    def GetGradient(self):
        return self.gradient

    def ComputeArea(self):
        return NotImplementedError
    
    def ComputePatchArea(self):
        self.patch_area = 0.0
        for elem in self.patch_elements:
            self.patch_area += elem.GetArea()
    
    def ComputeGradient(self):
        gradients_vertices = []
        for vertex in self.vertices:
            gradients_vertices.append(vertex.GetGradient())

        # Interpolazione
        #self.gradient = [gradient_x, gradient_y]
        return NotImplementedError

    def ComputePatchAveragedGradient(self):
        patch_avg_gradient = 0.0
        for elem in self.patch_elements:
            patch_avg_gradient += elem.GetArea() * elem.GetGradient()
        self.patch_avg_gradient = patch_avg_gradient / self.GetPatchArea()






