import copy
from utilities.input_output import *

class CMesh():
    """
    Class to convert .su2 to .mesh files and viceversa.
    Class to convert .csv or .dat SS2 solution files to .sol
    Works only with 2D-triangular and 3D-tetrahedral unstructured meshes.
    """

    def __init__(self, verbose=False):
        self.verbose = verbose
        return
    
    def SetDim(self, dim):
        """
        Setting the number of dimensions.
        """
        self.dim = dim

    def GetDim(self):
        """
        Returning the number of dimensions.
        """
        return self.dim
    
    def SetMeshDict(self, meshDict):
        """
        Setting the mesh dictionary once read either from .su2 or .mesh
        """
        self.meshDict = meshDict

    def GetMeshDict(self):
        """
        Returning the mesh dictionary.
        """
        return self.meshDict
    
    def SetMetricDict(self, metricDict):
        """
        Setting the dictionary for adaptation metric values.
        """
        self.metricDict = metricDict

    def GetMetricDict(self):
        """
        Returning the dictionary for adaptation metric values.
        """
        return self.metricDict
        
    def SetSU2MeditMarkersMap(self, su2MarkersList):
        """
        Constructing a unique between SU2 markers and Medit colors.
        """
        self.markersMap = []
        for meditTag, su2Tag in enumerate(su2MarkersList):
            # medit_tag starts from "1" since the tag "0" is left for the volume domain
            self.markersMap.append([str(meditTag+1), su2Tag])

    def GetSU2MeditMarkersMap(self):
        """
        Returning the map bewteen SU2 markers and Medit colors.
        """
        if self.markersMap:
            return self.markersMap
        else:
            raise ValueError('Su2-Medit markers map not set!')  

    def GetSU2Marker(self, meditTag):
        """
        Returning the SU2 markers correspondent to a Medit color.
        """
        for match in self.markersMap:
            if match[0] == meditTag:
                return match[1]

    def GetMeditMarker(self, su2Tag):
        """
        Returning the Medit color correspondent to a SU2 marker.
        """
        for match in self.markersMap:
            if match[1] == su2Tag:
                return match[0]

    def GetVertex(self, ID):
        return self.VertexElementlookup['Vertices'].get(ID)

    def GetElement(self, ID):
        return self.VertexElementlookup['Elements'].get(ID)

    def SetVertexElementLookup(self):
        keyElem = 'Triangles' if self.dim == 2 else 'Tetrahedra'
        self.VertexElementlookup = {
            'Vertices': {vert.GetID(): vert for vert in self.meshDict['Vertices']},
            'Elements': {elem.GetID(): elem for elem in self.meshDict[keyElem]}
        }

    def FinalizingDataStructure(self):
        '''
        For each vertex:
        - setting neighbouring elements
        - setting neighbouring points
        
        For each element
        - setting vertices through IDs
        - characterizng patch
        '''

        verticesID = [vert.GetID() for vert in self.meshDict['Vertices']]

        # helper to set the vertex eighbouring points
        verticesNeighbours = {vertID: set() for vertID in verticesID}
        # helper to set the vertex eighbouring elements
        elementsNeighbours = copy.deepcopy(verticesNeighbours)

        keyElem = 'Triangles' if self.dim == 2 else 'Tetrahedra'
        for element in self.meshDict[keyElem]:

            # setting the vertices (CVertex instances) belonging to each element
            element.SetVertices(self)
            eid = element.GetID()

            # adding the vertex neighbouring IDs iteratively looping on the elements
            vids = element.GetVerticesID()
            n = len(vids)
            for i in range(n):
                vi = vids[i]
                for j in range(i + 1, n):
                    vj = vids[j]
                    verticesNeighbours[vi].add(vj)
                    verticesNeighbours[vj].add(vi)
                    elementsNeighbours[vi].add(eid)
                    elementsNeighbours[vj].add(eid)

        for vert in self.meshDict['Vertices']:
            vert.SetVerticesNeighboursID(list(verticesNeighbours[vert.GetID()]))
            vert.SetElementsNeighboursID(list(elementsNeighbours[vert.GetID()]))

        # computing the mesh cardinality
        self.cardinality = len(self.meshDict[keyElem])
     
    def ReadMeshSU2(self, su2Filename):
        """ 
        Reads a .su2 mesh file and returns node coordinates, elements, and boundary markers in a dictionary data structure. 
        """
        meshDict = ReadSU2MeshASCII(self, su2Filename)

        return meshDict

    def ReadSolSU2(self, sensor, su2Filename):
        """
        Reads a .csv/.dat SU2 solution file to obtain the sensor and the gradient at each vertex. 
        """
        if '.dat' in su2Filename:
            ReadSU2RestartBinary(self, sensor, su2Filename)
        elif '.csv' in su2Filename:
            ReadSU2RestartASCII(self, sensor, su2Filename) 
    
    def ReadMeshMedit(self, meditFilename):
        """ 
        Reads a .mesh file and returns node coordinates, elements, and boundary markers  in a dictionary data structure. 
        """
        if meditFilename.endswith('.mesh'):
            meshDict = ReadMeditMeshASCII(self, meditFilename)
        elif meditFilename.endswith('.meshb'):
            meshDict = ReadMeditMeshBinary(self, meditFilename, verbose=self.verbose)

        return meshDict

    def WriteMeshSU2(self, su2Filename):
        """ 
        Writes a .su2 mesh file from given mesh data. 
        """
        WriteSU2MeshASCII(self, su2Filename)

        return
    
    def WriteMeshMedit(self, meditFilename):
        """ 
        Writes a .mesh mesh file from given mesh data. 
        """
        if meditFilename.endswith('.mesh'):
            WriteMeditMeshASCII(self, meditFilename)
        elif meditFilename.endswith('.meshb'):
            WriteMeditMeshBinary(self, meditFilename)
        
        return

    def WriteSolMedit(self, meditFilename):
        """ 
        Writes a .sol Medit file from given metric data. 
        """
        WriteMeditSolASCII(self, meditFilename)

        return
    
    def SU2ToMeditMesh(self, su2Filename, meditFilename):
        """
        Full mesh file conversion (reading-writing) from SU2 to Medit
        """
        self.ReadMeshSU2(su2Filename)
        self.WriteMeshMedit(meditFilename)
        if self.verbose:
            print(f"Converted {su2Filename} to {meditFilename}")

    def SU2ToMeditSol(self, su2Filename, meditFilename):
        """
        Full sol file conversion (reading-writing) from SU2 to Medit
        """
        self.ReadSolSU2(su2Filename)
        self.WriteSolMedit(meditFilename)
        if self.verbose:
            print(f"Converted {su2Filename} to {meditFilename}")

    def MeditToSU2Mesh(self, meditFilename, su2Filename):
        """
        Full mesh file conversion (reading-writing) from Medit to SU2
        """
        self.ReadMeshMedit(meditFilename)
        self.WriteMeshSU2(su2Filename)
        if self.verbose:
            print(f"Converted {meditFilename} to {su2Filename}")
   


