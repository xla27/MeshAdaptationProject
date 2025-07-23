
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
    
    def GetVolume(self):
        if not hasattr(self, 'volume'):
            self.ComputeVolume()
        return self.volume
    
    def GetPatchVolume(self):
        return self.patchVolume
        
    def GetGradient(self):
        if not hasattr(self, 'gradient'):
            self.ComputeGradient()
        return self.gradient
    
    def GetMetric(self):
        if not hasattr(self, 'metric'):
            self.ComputeMetric()
        return self.metric





