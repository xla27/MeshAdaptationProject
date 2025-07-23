import numpy as np
from scipy.spatial.distance import pdist
import copy

from elements import CTriangle, CTetrahedron
from CVertex import CVertex
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
        with open(su2Filename, "r") as f:
            lines = f.readlines()

        vertices = []
        elements = []
        boundaries = {}

        su2MarkersList = []
        
        i = 0
        while i < len(lines):
            line = lines[i].strip()

            if line.startswith("%"):
                i += 1
                continue

            if line.startswith("NDIME="):
                dim = int(line.split("=")[1].strip())
                self.SetDim(dim)

            elif line.startswith("NPOIN="):
                nVertices = int(line.split("=")[1].strip())
                if dim == 2:
                    for j in range(nVertices):
                        x, y = lines[i + 1 + j].split()[:dim]
                        vertex = CVertex(j, float(x), float(y))
                        vertices.append(vertex)
                elif dim == 3:
                    for j in range(nVertices):
                        x, y, z = lines[i + 1 + j].split()[:dim]
                        vertex = CVertex(j, float(x), float(y), z=float(z))
                        vertices.append(vertex)
                i += nVertices  # Move index past nodes

            elif line.startswith("NELEM="):
                nElements = int(line.split("=")[1].strip())
                for j in range(nElements):
                    elemType, *idVert = lines[i + 1 + j].split()[:(dim+2)]
                    element = CTriangle(j) if dim == 2 else CTetrahedron(j)
                    element.SetVerticesID([int(vert) for vert in idVert])
                    elements.append(element)

                i += nElements  # Move index past elements

            elif line.startswith("NMARK="):
                nMarkers = int(line.split("=")[1].strip())
                for _ in range(nMarkers):
                    i += 1
                    markerTag = lines[i].split("=")[1].strip()
                    su2MarkersList.append(markerTag)
                    i += 1
                    nFaces = int(lines[i].split("=")[1].strip())
                    boundaries[markerTag] = []
                    for j in range(nFaces):
                        faceType, *idVert = lines[i + 1 + j].split()
                        boundaries[markerTag].append([int(vert) for vert in idVert])
                    i += nFaces  # Move index past boundary elements

            i += 1

        meshDict = {'Dim': dim, 'Vertices': vertices}
        if int(elemType) == 5:
            meshDict['Triangles'] = elements
        if int(elemType) == 10:
            meshDict['Tetrahedra'] = elements

        # reordering markers in alphabetical orders
        #boundaries = {key: value for key, value in sorted(boundaries.items())}
        if int(faceType) == 3:
            meshDict['Edges'] =  boundaries
        if int(faceType) == 5:
            meshDict['Triangles'] =  boundaries

        self.SetMeshDict(meshDict)
        self.SetVertexElementLookup()
        self.SetSU2MeditMarkersMap(su2MarkersList)

        return meshDict

    def ReadSolSU2(self, sensor, su2Filename):
        """
        Reads a .csv/.dat SU2 solution file to obtain the sensor and the gradient at each vertex. 
        """
        if '.dat' in su2Filename:
            read_SU2_restart_binary(self, sensor, su2Filename)
        elif '.csv' in su2Filename:
            read_SU2_restart_ascii(self, sensor, su2Filename) 
    
    def ReadMeshMedit(self, meditFilename):
        """ 
        Reads a .mesh file and returns node coordinates, elements, and boundary markers  in a dictionary data structure. 
        """
        with open(meditFilename, "r") as f:
            lines = f.readlines()

        dim = None
        vertices = []
        elements = []
        boundaries = {}

        i = 0
        while i < len(lines):
            line = lines[i].strip()

            if line.startswith("Dimension"):
                dim = int(line.split()[1])
                self.SetDim(dim)

            elif line.startswith("Vertices"):
                nVert = int(lines[i+1])
                if dim == 2:
                    for j in range(nVert):
                        x, y, idDom = lines[i + 2 + j].split()
                        vertex = CVertex(j, float(x), float(y))
                        vertices.append(vertex)
                elif dim == 3:
                    for j in range(nVert):
                        x, y, z, idDom = lines[i + 2 + j].split()
                        vertex = CVertex(j, float(x), float(y), z=float(z))
                        vertices.append(vertex)
                    
                i += nVert  # Move index past nodes

            elif (line.startswith("Triangles") and dim == 2) or (line.startswith("Tetrahedra") and dim == 3):
                nElements = int(lines[i+1])
                elements = []
                for j in range(nElements):
                    elemData = list(map(int, lines[i + 2 + j].split()))
                    elemData = [elem-1 for elem in elemData]
                    element = CTriangle(j) if dim == 2 else CTetrahedron(j)
                    element.SetVerticesID([int(vert) for vert in elemData[:-1]]) # Last column is a region marker
                    elements.append(element)
                i += nElements  # Move index past elements             

            elif (line.startswith("Edges") and dim == 2) or (line.startswith("Triangles") and dim == 3):
                # These define boundary markers
                nFaces = int(lines[i+1])
                for j in range(nFaces):
                    faceData = list(map(int, lines[i + 2 + j].split()))
                    marker = str(faceData[-1])  # Last column is the boundary marker
                    faceData = [face - 1 for face in faceData[:-1]]
                    if marker not in boundaries.keys():
                        boundaries[marker] = []
                    boundaries[marker].append(faceData)  # Store only connectivity
                i += nFaces  # Move index past boundary elements

            i += 1

        meshDict = {'Dim': dim, 'Vertices': vertices}
        if dim == 2:
            meshDict['Triangles'] = elements
        if dim == 3:
            meshDict['Tetrahedra'] = elements

        # reordering markers in alphabetical orders
        if dim == 2:
            meshDict['Edges'] =  boundaries
        if dim == 3:
            meshDict['Triangles'] =  boundaries

        self.SetMeshDict(meshDict)

        return meshDict

    def WriteMeshSU2(self, su2Filename):
        """ 
        Writes a .su2 mesh file from given mesh data. 
        """
        mesh = self.GetMeshDict()
        dim = mesh['Dim']
        vertices = mesh["Vertices"]
        if dim == 2:
            elements = mesh['Triangles']
            elemType = 5
            boundaries = mesh['Edges']
            faceType = 3
        if dim == 3:
            elements = mesh['Tetrahedra']
            elemType = 10
            boundaries = mesh['Triangles']
            faceType = 5


        with open(su2Filename, "w") as f:
            f.write("NDIME= {}\n".format(dim))

            f.write("NELEM= {}\n".format(len(elements)))
            for i, elem in enumerate(elements):
                f.write("{} {} {}\n".format(elemType, " ".join(map(str, elem.GetVerticesID())), elem.GetID()))

            f.write("NPOIN= {}\n".format(len(vertices)))
            for i, vert in enumerate(vertices):
                f.write("{} {}\n".format(" ".join(map(str, vert.GetCoordinates())), vert.GetID()))

            f.write("NMARK= {}\n".format(len(boundaries)))
            for meditTag in boundaries.keys():
                f.write("MARKER_TAG= {}\n".format(self.GetSU2Marker(meditTag)))
                f.write("MARKER_ELEMS= {}\n".format(len(boundaries[meditTag])))
                for face in boundaries[meditTag]:
                    f.write("{} {}\n".format(faceType, " ".join(map(str, face))))

        return
    
    def WriteMeshMedit(self, meditFilename):
        """ 
        Writes a .mesh mesh file from given mesh data. 
        """
        mesh = self.GetMeshDict()
        dim = mesh['Dim']
        vertices = mesh['Vertices']
        if dim == 2:
            elements = mesh['Triangles']
            boundaries = mesh['Edges']
        if dim == 3:
            elements = mesh['Tetrahedra']
            boundaries = mesh['Triangles']

        with open(meditFilename, "w") as f:
            f.write("MeshVersionFormatted 2\n")
            f.write("Dimension {}\n".format(mesh['Dim']))
            
            # Write nodes
            f.write("\nVertices \n{}\n".format(len(mesh['Vertices'])))
            for vert in vertices:
                f.write(" ".join(map(str, vert.GetCoordinates())) + " 0\n")  # 0 is the default region ID
            
            # Write elements (assume triangles for 2D, tetrahedra for 3D)
            if dim == 2:
                f.write("\nTriangles \n{}\n".format(len(elements)))
            elif dim == 3:
                f.write("\nTetrahedra \n{}\n".format(len(elements)))
            
            for elem in elements:
                elemData = [el + 1 for el in elem.GetVerticesID()]
                f.write(" ".join(map(str, elemData)) + " 0\n")  # Last value is a region ID
            
            # Write boundary elements correctly
            if boundaries:
                if dim == 2:
                    f.write("\nEdges \n{}\n".format(sum(len(faces) for faces in boundaries.values())))
                elif dim == 3:
                    f.write("\nTriangles \n{}\n".format(sum(len(faces) for faces in boundaries.values())))
                
                for su2Tag in boundaries.keys():
                    for face in boundaries[su2Tag]:
                        face = [fa + 1 for fa in face]
                        f.write(" ".join(map(str, face)) + " {}\n".format(self.GetMeditMarker(su2Tag)))

            f.write("\nEnd\n")
        
        return

    def WriteSolMedit(self, meditFilename):
        """ 
        Writes a .sol Medit file from given metric data. 
        """

        dim = self.meshDict['Dim']
        nVert = len(self.meshDict['Vertices'])

        header = 'MeshVersionFormatted 2\nDimension %i\nSolAtVertices\n%i\n1 3\n' % (dim , nVert)
        footer = '\nEnd\n'

        metricDim = 3 if dim == 2 else 6
        solData = np.zeros((nVert, metricDim))

        for iVert, vertex in enumerate(self.meshDict['Vertices']):
            metric = vertex.GetMetric()

            solData[iVert, 0] = metric[0,0]
            solData[iVert, 1] = metric[0,1]
            solData[iVert, 2] = metric[1,1]

            if dim == 3:
                solData[iVert, 3] = metric[0,2]
                solData[iVert, 4] = metric[1,2]
                solData[iVert, 5] = metric[2,2]

        np.savetxt(meditFilename, solData, delimiter=' ', header=header, footer=footer, comments='', fmt='%1.5e')
        
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

    def WriteParamFile(self, configMmg, meshFilename):
        """
        Writing the .mmg2d/.mmg3d parameter file if required. 
        """
        paramRequired = isinstance(configMmg['hausd'], dict)
        if paramRequired:
            mesh = self.GetMeshDict()
            dim = mesh['Dim']
            if dim == 2:
                boundaries = mesh['Edges']
                elemType = 'Edges'
                mmgExt = '.mmg2d'
            if dim == 3:
                boundaries = mesh['Triangles']
                elemType = 'Triangles'
                mmgExt = '.mmg3d'

            paramFilename = meshFilename + mmgExt
            with open(paramFilename, 'w') as f:
                f.write('Parameters\n')
                f.write(str(len(boundaries.keys()))+'\n')
                f.write('\n')

                if len(boundaries.keys()) != len(configMmg['hausd'].keys()):
                    print('WARNING: Different number of markers between SU2 (%i) mesh and MMG parameters (%i). ' \
                                  'For unspecified markers, HAUSD = 0.01 is assumed.' %
                                   (len(boundaries.keys()), len(configMmg['hausd'].keys())))
                
                for su2Tag in configMmg['hausd'].keys():
                    f.write('%s %s %1.2e %1.2e %1.2e\n' % 
                            (self.GetMeditMarker(su2Tag), 
                            elemType, 
                            configMmg['hmin'], 
                            configMmg['hmax'], 
                            configMmg['hausd'][su2Tag]))       
        
        else:
            pass

        return      

# Reading helpers

CGNS_STRING_SIZE = 33  # Fixed string size per CGNS standard

def read_SU2_restart_binary(mesh, sensor, filename):
    """
    Read SU2 binary restart file and return fields and data array.

    Returns:
        fields (List[str]): Field names including "Point_ID".
        data (np.ndarray): Data array of shape (nPoints, nFields-1).

    Note that the Point_ID column is implicit in the ordering
    """

    restartFields = []  

    with open(filename, 'rb') as f:
        # Read 5 integers (magic number + metadata)
        header = np.fromfile(f, dtype=np.int32, count=5)
        if header.size != 5:
            raise RuntimeError("Error reading header from restart file.")
        
        magic_number, nFields, nPoints, _, _ = header

        # Check the magic number
        if magic_number != 535532:
            raise RuntimeError(f"{filename} is not a binary SU2 restart file.")

        # Read field names (each is CGNS_STRING_SIZE characters)
        for _ in range(nFields):
            name_bytes = f.read(CGNS_STRING_SIZE)
            name_str = name_bytes.decode('utf-8').strip('\x00').strip()
            restartFields.append(name_str)

        # Read restart data as a flat array of doubles
        data = np.fromfile(f, dtype=np.float64, count=nFields * nPoints)

        if data.size != nFields * nPoints:
            raise RuntimeError("Error reading restart data.")

        # Reshape to 2D: each row is a point, each column is a field
        data = data.reshape((nPoints, nFields))

    try:
        meshDict = mesh.GetMeshDict()
    except:
        raise ValueError('The mesh has not been read yet!')

    if 'z' in restartFields:
        meshDict['Dim'] = 3
        fieldsToRead = [sensor, 'Grad(Sensor)_x', 'Grad(Sensor)_y', 'Grad(Sensor)_z']
    else:
        meshDict['Dim'] = 2  
        fieldsToRead = [sensor, 'Grad(Sensor)_x', 'Grad(Sensor)_y']

    coords = np.zeros((nPoints, meshDict['Dim']))
    for iPoint in range(nPoints):
        vert = mesh.GetVertex(iPoint)
        solution = 0.0
        gradient = []
        for field in fieldsToRead:
            iField = restartFields.index(field)
            if field == sensor:
                solution += data[iPoint, iField]
            else:
                gradient.append(data[iPoint, iField])

        vert.SetSolution(solution)
        vert.SetGradient(gradient)

        coords[iPoint, 0] = data[iPoint, restartFields.index('x')]
        coords[iPoint, 1] = data[iPoint, restartFields.index('y')]
        if meshDict['Dim'] == 3:
            coords[iPoint, 2] = data[iPoint, restartFields.index('z')]

    mesh.diameter = np.amax(pdist(coords))


def read_SU2_restart_ascii(mesh, sensor, filename):
    """
    Read SU2 ASCII restart file and return fields and data array.

    Returns:
        fields (List[str]): Field names.
        data (np.ndarray): Data array of shape (nPoints, nFields).
    """

    # reading the first line to get the fields name
    try:
        with open(filename, "r") as f:
            line = f.readline()
    except:
        raise("The solution file must be in ASCII format!")

    restartFields = line.lstrip('"').rstrip('"\n')
    restartFields = restartFields.split('","')

    data = np.genfromtxt(filename, delimiter=',', skip_header=1, dtype=np.float64)
    nPoints = data.shape[0]

    try:
        meshDict = mesh.GetMeshDict()
    except:
        raise ValueError('The mesh has not been read yet!')

    if 'z' in restartFields:
        meshDict['Dim'] = 3
        fieldsToRead = [sensor, 'Grad(Sensor)_x', 'Grad(Sensor)_y', 'Grad(Sensor)_z']
    else:
        meshDict['Dim'] = 2  
        fieldsToRead = [sensor, 'Grad(Sensor)_x', 'Grad(Sensor)_y']

    coords = np.zeros((nPoints, meshDict['Dim']))
    for iPoint in range(nPoints):
        vert = mesh.GetVertex(iPoint)
        solution = 0.0
        gradient = []
        for field in fieldsToRead:
            iField = restartFields.index(field)
            if field == sensor:
                solution += data[iPoint, iField]
            else:
                gradient.append(data[iPoint, iField])

        vert.SetSolution(solution)
        vert.SetGradient(gradient)

        coords[iPoint, 0] = data[iPoint, restartFields.index('x')]
        coords[iPoint, 1] = data[iPoint, restartFields.index('y')]
        if meshDict['Dim'] == 3:
            coords[iPoint, 2] = data[iPoint, restartFields.index('z')]

    mesh.diameter = np.amax(pdist(coords))