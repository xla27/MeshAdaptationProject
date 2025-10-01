import sys, os, shutil
import subprocess
import time
import numpy as np
from struct import unpack, pack, calcsize, iter_unpack
import itertools

from geometry.elements import CTriangle, CTetrahedron
from geometry.CVertex import CVertex

# ------------------------------------------------------------
#  Mesh reading
# ------------------------------------------------------------

def ReadSU2MeshASCII(mesh, meshFilename):
    """ 
    Reads a .su2 mesh file and returns node coordinates, elements, and boundary markers in a dictionary data structure. 
    """
    with open(meshFilename, "r") as f:
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
            mesh.SetDim(dim)

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
                    faceType, *id_vert = lines[i + 1 + j].split()
                    boundaries[markerTag].append([int(vert) for vert in id_vert])
                i += nFaces  # Move index past boundary elements

        i += 1

    if int(elemType) != 5 and int(elemType) != 10:
        raise Exception('The .su2 mesh file containes volume elements different from triangles/tetrahedra')

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

    mesh.SetMeshDict(meshDict)
    mesh.SetVertexElementLookup()
    mesh.SetSU2MeditMarkersMap(su2MarkersList)

    return meshDict


def ReadMeditMeshASCII(mesh, meshFilename):

    with open(meshFilename, "r") as f:
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
            mesh.SetDim(dim)

        elif line.startswith("Vertices"):
            print(line)
            print(lines[i+1])
            nVertices = int(lines[i+1])
            if dim == 2:
                for j in range(nVertices):
                    vertex = CVertex(j, float(x), float(y))
                    vertices.append(vertex)
            elif dim == 3:
                for j in range(nVertices):
                    x, y, z, id_dom = lines[i + 2 + j].split()
                    vertex = CVertex(j, float(x), float(y), z=float(z))
                    vertices.append(vertex)                
            i += nVertices  # Move index past nodes

        elif (line.startswith("Triangles") and dim == 2) or (line.startswith("Tetrahedra") and dim == 3):
            nElements = int(lines[i+1])
            elements = [None]*nElements
            for j in range(nElements):
                elemData = list(map(int, lines[i + 2 + j].split()))
                elemData = [el-1 for el in elemData]
                element = CTriangle(j) if dim == 2 else CTetrahedron(j)
                element.SetVerticesID([int(vert) for vert in elemData[:-1]])  # Last column is a region marker
                elements.append(element)
            i += nElements  # Move index past elements             

        elif (line.startswith("Edges") and dim == 2) or (line.startswith("Triangles") and dim == 3):
            # These define boundary markers
            nFaces = int(lines[i+1])
            for j in range(nFaces):
                faceData = list(map(int, lines[i + 2 + j].split()))
                marker = str(faceData[-1])  # Last column is the boundary marker
                faceData = [fa - 1 for fa in faceData[:-1]]
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

    meshDict = CheckSurplusPoints(meshDict)

    mesh.SetMeshDict(meshDict)

    return meshDict


def ReadMeditMeshBinary(mesh, meshFilename, verbose=False):
    """Reload the mesh with meshFile (binary format).

    INPUT:
        meshFile :  the path of the .meshb file (binary format)
    """
    if not meshFilename.endswith(".meshb"):
        meshFilename += ".meshb"

    vertices = []
    elements = []
    boundaries = {}

    tic()
    f = open(meshFilename, "rb")

    code = read_int(f)
    if not code:
        raise Exception("Error in reading the binary file " + meshFilename)
    
    meshVersionFormatted = read_int(f)
    if not meshVersionFormatted in [1, 2]:
        raise Exception("Wrong MeshVersionFormatted. Should be 1 or 2.")
    
    if gmfKwdCod[read_int(f)] != 'GmfDimension':
        raise Exception("Error in reading the binary file " + meshFilename + ". "
                        "GmfDimension expected.")
    
    nextPos = read_int(f)
    if verbose: print(f"Next field position: {nextPos}.")
    dim = read_int(f)
    if verbose: print("Dimension " + str(dim))
    if not dim in [2, 3]:
        raise Exception("Error in reading the binary file " + meshFilename +
                        ". The dimension should be 2 or 3.")
    else:
        mesh.SetDim(dim)

    def readField(f, nfields, title):
        nextPos = read_int(f)
        if verbose: print(f"Next field position: {nextPos}.")
        n = read_int(f)
        if verbose: print(f"{n} {title}.")
        code = "i" * nfields * n
        nbytes = 4 * nfields * n
        field = np.asarray(unpack(code, f.read(nbytes)))
        return (n, field.reshape((n, nfields)))

    while True:
        try:
            kwdCod = read_int(f)
        except:
            if verbose: print("Warning: end of file without the END keyword.")
            break
        if verbose: print("Reading "+gmfKwdCod[kwdCod])

        # Vertices
        if gmfKwdCod[kwdCod] == 'GmfVertices':

            nextPos = read_int(f)
            if verbose: print(f"Next field position: {nextPos}.")

            nvert = read_int(f)
            if verbose: print(f"{nvert} Vertices.")

            if meshVersionFormatted == 1:
                # Float precision
                code = "="+("f"*dim + "i") * nvert
                nbytes = calcsize(code)
            elif meshVersionFormatted == 2:
                # Double precision
                code = "="+("d"*dim + "i") * nvert
                nbytes = calcsize(code)

            verticesArray = np.asarray(unpack(code, f.read(nbytes)))
            verticesList = verticesArray.reshape((nvert, dim+1)).tolist()
            if dim == 2:
                for j, vert in enumerate(verticesList):
                    vertex = CVertex(j, vert[0], vert[1])
                    vertices.append(vertex)
            else:
                for j, vert in enumerate(verticesList):
                    vertex = CVertex(j, vert[0], vert[1], z=vert[2])
                    vertices.append(vertex)

        # Edges
        elif gmfKwdCod[kwdCod] == 'GmfEdges':
            nEdges, edges = readField(f, 3, gmfKwdCod[kwdCod][3:])
            edges.tolist()

            if dim == 2:
                for edg in edges:
                    if str(edg[-1]) not in boundaries.keys():
                        boundaries[str(edg[-1])] = []
                    boundaries[str(edg[-1])].append([ed - 1 for ed in edg[:-1]])

        # Triangles
        elif gmfKwdCod[kwdCod] == 'GmfTriangles':
            nTrias, triangles = readField(f, 4, gmfKwdCod[kwdCod][3:])
            triangles.tolist()

            if dim == 2: # triangles as elements
                for j, tria in enumerate(triangles):
                    element = CTriangle(j)
                    element.SetVerticesID([int(tr-1) for tr in tria[:-1]])
                    elements.append(element)

            elif dim == 3: # triangles as boundaries
                for tria in triangles:
                    if str(tria[-1]) not in boundaries.keys():
                        boundaries[str(tria[-1])] = []
                    boundaries[str(tria[-1])].append([tr-1 for tr in tria[:-1]])

        # Tetrahedra
        elif gmfKwdCod[kwdCod] == 'GmfTetrahedra':
            nTetra, tetrahedra = readField(f, 5, gmfKwdCod[kwdCod][3:])
            tetrahedra.tolist()
            for j, tetra in enumerate(tetrahedra):
                element = CTetrahedron(j)
                element.SetVerticesID([int(te-1) for te in tetra[:-1]])
                elements.append(element)

        # Fields that need to be read but not saved
        elif (gmfKwdCod[kwdCod] == 'GmfCorners'          or 
              gmfKwdCod[kwdCod] == 'GmfRequiredVertices' or
              gmfKwdCod[kwdCod] == 'GmfRequiredEdges'    or
              gmfKwdCod[kwdCod] == 'GmfRidges'           or
              gmfKwdCod[kwdCod] == 'GmfRequiredTriangles'):
            _, _ = readField(f, 1, gmfKwdCod[kwdCod][3:])

        elif (gmfKwdCod[kwdCod] == 'GmfNormalAtVertices' or 
              gmfKwdCod[kwdCod] == 'GmfTangentAtVertices'):
            _, _ = readField(f, 2, gmfKwdCod[kwdCod][3:])


        elif (gmfKwdCod[kwdCod] == 'GmfNormals' or
              gmfKwdCod[kwdCod] == 'GmfTangents'):
            nextPos = read_int(f)
            if verbose: print(f"Next field position: {nextPos}.")
            nvn = read_int(f)
            if verbose: print(f"{nvn} Normals.")
            if meshVersionFormatted == 1:
                code = "f"*(3*nvn)
                nbytes = calcsize(code)
            else:
                code = "d"*(3*nvn)
                nbytes = calcsize(code)
            _ = np.asarray(unpack(code, f.read(nbytes)))

        # End
        elif gmfKwdCod[kwdCod] == 'GmfEnd':
            if verbose: print("End of mesh.")
            break

        else:
            raise KeyError("Error, field "+gmfKwdCod[kwdCod]+" not supported.")
        
    f.close()
    if verbose: print("Read " + meshFilename + " in "+toc()+".")

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

    meshDict = CheckSurplusPoints(meshDict)

    mesh.SetMeshDict(meshDict)
    
    return meshDict


def CheckSurplusPoints(meshDict):
    """
    Routine to check and eliminate for surplus points in medit mesh generation
    """

    dim = meshDict['Dim']

    vertices = meshDict['Vertices']
    nVertices = len(vertices)

    if dim == 2:
        elements   = meshDict['Triangles']
        boundaries = meshDict['Edges']
    elif dim == 3:
        elements   = meshDict['Tetrahedra']
        boundaries = meshDict['Triangles']

    print("Checking for surplus points")

    PointIDs = np.array(range(nVertices), dtype=int)
    elements = np.array([elem.GetVerticesID() for elem in elements], dtype=int)
    isIn = np.isin(PointIDs, elements)
    whichAreSurplusPoints = np.where(isIn == False)[0]
    
    if (len(whichAreSurplusPoints) > 0):
        print("There are surplus points in mesh file for unknown reasons.")
        # print("Surplus points:", whichAreSurplusPoints)
        # print("Deleting it and then fix connectivity...")
        # SubtractIDs = np.zeros((nVertices, ), dtype=int)
        # for iPoint in whichAreSurplusPoints:
        #     IDsOver = np.where(PointIDs > iPoint)[0]
        #     SubtractIDs[IDsOver] += 1
        
        # # Now I can just fix the connectivity
        # elements -= SubtractIDs[elements]
        # for marker in boundaries.keys():
        #     markerData = np.array(boundaries[marker], dtype=int)
        #     markerData -= SubtractIDs[markerData]
        #     boundaries[marker] = markerData

        # whichAreSurplusPointsReverse = np.sort(whichAreSurplusPoints)[::-1]
        # for iPoint in whichAreSurplusPointsReverse:
        #     # Remove from the list of vertices but in reverse order
        #     vertices.pop(iPoint)
        raise KeyError(f'This function needs to be reimplemented!')

    meshDict['Vertices'] = vertices
    # if dim == 2:
    #     meshDict['Triangles'] = elements
    # if dim == 3:
    #     meshDict['Tetrahedra'] = elements

    # reordering markers in alphabetical orders
    if dim == 2:
        meshDict['Edges'] = boundaries
    if dim == 3:
        meshDict['Triangles'] =  boundaries

    return meshDict


# ------------------------------------------------------------
#  Mesh writing
# ------------------------------------------------------------

def WriteSU2MeshASCII(mesh, meshFilename):
    """ 
    Writes a .su2 mesh file from given mesh data. 
    """
    meshDict = mesh.GetMeshDict()
    dim = meshDict['Dim']
    vertices = meshDict["Vertices"]
    if dim == 2:
        elements = meshDict['Triangles']
        elemType = 5
        boundaries = meshDict['Edges']
        faceType = 3
    if dim == 3:
        elements = meshDict['Tetrahedra']
        elemType = 10
        boundaries = meshDict['Triangles']
        faceType = 5


    with open(meshFilename, "w") as f:
        f.write("NDIME= {}\n".format(dim))

        f.write("NELEM= {}\n".format(len(elements)))
        for i, elem in enumerate(elements):
            f.write("{} {} {}\n".format(elemType, " ".join(map(str, elem.GetVerticesID())), elem.GetID()))

        f.write("NPOIN= {}\n".format(len(vertices)))
        if dim == 2:
            for i, vert in enumerate(vertices):
                f.write("{} {}\n".format(" ".join(map(str, vert.GetCoordinates()[:-1])), vert.GetID()))
        elif dim == 3:
            for i, vert in enumerate(vertices):
                f.write("{} {}\n".format(" ".join(map(str, vert.GetCoordinates())), vert.GetID()))

        f.write("NMARK= {}\n".format(len(boundaries)))
        for meditTag in boundaries.keys():
            f.write("MARKER_TAG= {}\n".format(mesh.GetSU2Marker(meditTag)))
            f.write("MARKER_ELEMS= {}\n".format(len(boundaries[meditTag])))
            for face in boundaries[meditTag]:
                f.write("{} {}\n".format(faceType, " ".join(map(str, face))))

    return


def WriteMeditMeshASCII(mesh, meshFilename):
    """ 
    Writes a .mesh mesh file from given mesh data. 
    """
    meshDict = mesh.GetMeshDict()
    dim = meshDict['Dim']
    vertices = meshDict['Vertices']
    if dim == 2:
        elements = meshDict['Triangles']
        boundaries = meshDict['Edges']
    if dim == 3:
        elements = meshDict['Tetrahedra']
        boundaries = meshDict['Triangles']

    with open(meshFilename, "w") as f:
        f.write("MeshVersionFormatted 2\n")
        f.write("Dimension {}\n".format(meshDict['Dim']))
        
        # Write nodes
        f.write("\nVertices \n{}\n".format(len(mesh['Vertices'])))
        if dim == 2:
            for vert in vertices:
                f.write(" ".join(map(str, vert.GetCoordinates()[:-1])) + " 0\n")  # 0 is the default region ID
        elif dim == 3:
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
                    f.write(" ".join(map(str, face)) + " {}\n".format(mesh.GetMeditMarker(su2Tag)))

        f.write("\nEnd\n")
    
    return
    

def WriteMeditMeshBinary(mesh, meshFilename, verbose=False):
    """Save a mesh in the INRIA binary file format"""
    
    if not meshFilename.endswith(".meshb"):
        meshFilename += ".meshb"

    meshDict = mesh.GetMeshDict()
    dim = meshDict['Dim']
    vertices = meshDict['Vertices']
    if dim == 2:
        elements = meshDict['Triangles']
        boundaries = meshDict['Edges']
    elif dim == 3:
        elements = meshDict['Tetrahedra']
        boundaries = meshDict['Triangles']

    tic()
    f = open(meshFilename, "wb")
    f.write(pack("i", 1))  # Write code
    f.write(pack("i", 2))  # Write MeshVersionFormatted 2

    f.write(pack("i", indicesGmf['GmfDimension']))
    f.write(pack("i", 20))  # NextPos
    if not dim in [2, 3]:
        raise Exception("Error, the mesh dimension is not 2 or 3")
    f.write(pack("i", dim))

    # Vertices
    if verbose: print(f"Write {len(vertices)} Vertices.")
    f.write(pack("i", indicesGmf['GmfVertices']))
    code = "=" + ("d" * dim + "i") * len(vertices)
    nextPos = f.tell() + calcsize("ii") + calcsize(code)
    if verbose: print(f"Next position: {nextPos}")
    f.write(pack("i", nextPos))
    verticesList = []
    if dim == 2:
        for vert in vertices:
            x, y, _ = vert.GetCoordinates()
            verticesList.append([x, y, 0]) # appending the region ID
    elif dim == 3:
        for vert in vertices:
            x, y, z = vert.GetCoordinates()
            verticesList.append([x, y, z, 0]) # appending the region ID
    f.write(pack("i", len(verticesList)))
    data = list(itertools.chain.from_iterable(verticesList))
    f.write(pack(code, *data))

    # Edges
    if dim == 2:
        nEdges = sum(len(faces) for faces in boundaries.values())
        if verbose: print(f"Write {nEdges} "+"Edges.")
        f.write(pack("i", indicesGmf['GmfEdges']))
        code = "i" * 3 * nEdges

        nextPos = f.tell()+calcsize(code)+calcsize("ii")
        if verbose: print(f"Next position: {nextPos}")
        f.write(pack("i", nextPos))  # NulPos
        f.write(pack("i", nEdges))
        edges = []
        for su2Tag in boundaries.keys():
            for edg in boundaries[su2Tag]:
                edges.append([ed+1 for ed in edg] + [int(mesh.GetMeditMarker(su2Tag))])
        data = list(itertools.chain.from_iterable(edges))
        f.write(pack(code, *data))

    # Triangles
    if dim == 2: # triangles as elements
        nTrias = len(elements)
    elif dim == 3: # triangles as boundaries
        nTrias = sum(len(faces) for faces in boundaries.values())

    if verbose: print(f"Write {nTrias} "+"Triangles.")
    f.write(pack("i", indicesGmf['GmfTriangles']))
    code = "i" * 4 * nTrias

    nextPos = f.tell()+calcsize(code)+calcsize("ii")
    if verbose: print(f"Next position: {nextPos}")
    f.write(pack("i", nextPos))  # NulPos
    f.write(pack("i", nTrias))
    triangles = []
    if dim == 2:
        for elem in elements:
            triangles.append([el + 1 for el in elem.GetVerticesID()] + [0])
    elif dim == 3:
        for su2Tag in boundaries.keys():
            for tria in boundaries[su2Tag]:
                triangles.append([tr + 1 for tr in tria] +[int(mesh.GetMeditMarker(su2Tag))])   

    data = list(itertools.chain.from_iterable(triangles))
    f.write(pack(code, *data))

    # Tetrahedra
    if dim == 3:
        nTetras = len(elements)

        if verbose: print(f"Write {nTetras} "+"Tetrahedra.")
        f.write(pack("i", indicesGmf['GmfTetrahedra']))
        code = "i" * 5 * nTetras

        nextPos = f.tell()+calcsize(code)+calcsize("ii")
        if verbose: print(f"Next position: {nextPos}")
        f.write(pack("i", nextPos))  # NulPos
        f.write(pack("i", nTetras))

        tetrahedra = []
        for elem in elements:
            tetrahedra.append([el + 1 for el in elem.GetVerticesID()] + [0])

        data = list(itertools.chain.from_iterable(tetrahedra))
        f.write(pack(code, *data))

    # End
    f.write(pack("i", indicesGmf['GmfEnd']))
    nextPos = f.tell()+calcsize("i")
    # Final size
    f.write(pack("i",nextPos))
    f.close()

    if verbose: print("Wrote "+meshFilename+" in "+toc()+".")

    return

# ------------------------------------------------------------
#  SU2 Restart reading
# ------------------------------------------------------------

CGNS_STRING_SIZE = 33  # Fixed string size per CGNS standard

def ReadSU2RestartBinary(mesh, sensor, filename):
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
        
        magicNumber, nFields, nPoints, _, _ = header

        # Check the magic number
        if magicNumber != 535532:
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

    mesh.diameter = ApproximateDiameter(coords)


def ReadSU2RestartASCII(mesh, sensor, filename):
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

    mesh.diameter = ApproximateDiameter(coords)


def ApproximateDiameter(points, start_idx=None):
    """
    Farthest Point Sampling (FPS) for 2D/3D points.
    
    Parameters
    ----------
    points : ndarray of shape (M, dim)
        Input points in 3D space.
    start_idx : int or None
        Index of starting point. If None, picks a random point.
    
    Returns
    -------
    p1_idx, p2_idx, max_dist
        Indices of two farthest points (approximation) and the distance between them.
    """
    M = points.shape[0]

    # Step 1: choose starting point
    if start_idx is None:
        start_idx = np.random.randint(M)
    
    # Step 2: find farthest point from start
    dists = np.linalg.norm(points - points[start_idx,:], axis=1)
    p2_idx = np.argmax(dists)

    # Step 3: find farthest point from p2
    dists = np.linalg.norm(points - points[p2_idx,:], axis=1)
    p3_idx = np.argmax(dists)

    # Step 4: output
    max_dist = np.linalg.norm(points[p3_idx] - points[p2_idx])
    return max_dist


# ------------------------------------------------------------
#  Metric Prescription
# ------------------------------------------------------------

def WriteMeditSolASCII(mesh, meditFilename):
    """ 
    Writes a .sol Medit file from given metric data. 
    """
    dim = mesh.meshDict['Dim']
    nVert = len(mesh.meshDict['Vertices'])

    header = 'MeshVersionFormatted 2\nDimension %i\nSolAtVertices\n%i\n1 3\n' % (dim , nVert)
    footer = '\nEnd\n'

    metricDim = 3 if dim == 2 else 6
    solData = np.zeros((nVert, metricDim))

    for iVert, vertex in enumerate(mesh.meshDict['Vertices']):
        metric = vertex.GetMetric()

        solData[iVert, 0] = metric[0,0]
        solData[iVert, 1] = metric[0,1]
        solData[iVert, 2] = metric[1,1]

        if dim == 3:
            solData[iVert, 3] = metric[0,2]
            solData[iVert, 4] = metric[1,2]
            solData[iVert, 5] = metric[2,2]

    np.savetxt(meditFilename, solData, delimiter=' ', header=header, footer=footer, comments='', fmt='%1.5e')

# ------------------------------------------------------------
#  MMG Parameter file
# ------------------------------------------------------------

def WriteParamFile(mesh, configMmg, meshFilename):
    """
    Writing the .mmg2d/.mmg3d parameter file if required. 
    """
    paramRequired = isinstance(configMmg['hausd'], dict)
    if paramRequired:
        meshDict = mesh.GetMeshDict()
        dim = meshDict['Dim']
        if dim == 2:
            boundaries = meshDict['Edges']
            elemType = 'Edges'
            mmgExt = '.mmg2d'
        if dim == 3:
            boundaries = meshDict['Triangles']
            elemType = 'Triangles'
            mmgExt = '.mmg3d'

        if meshFilename.endswith('.mesh'):
            paramFilename = meshFilename.replace('.mesh', mmgExt)
        elif meshFilename.endswith('.meshb'):
            paramFilename = meshFilename.replace('.meshb', mmgExt)

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
                        (mesh.GetMeditMarker(su2Tag), 
                        elemType, 
                        configMmg['hmin'], 
                        configMmg['hmax'], 
                        configMmg['hausd'][su2Tag]))       
    
    else:
        pass
    
# ------------------------------------------------------------
#  Other utilities for I/O
# ------------------------------------------------------------

def read_int(f):
    return unpack("i", f.read(4))[0]


def next_line(f):
    while True:
        line = f.readline().strip()
        if line:
            break
    return line


def readField(f, nfields, title, verbose):
    nextPos = read_int(f)
    print(f"Next field position: {nextPos}.")
    n = read_int(f)
    print(f"{n} {title}.")
    code = "i"*nfields*n
    nbytes = 4*nfields*n
    field = np.asarray(unpack(code, f.read(nbytes)))
    return (n, field.reshape((n, nfields)))


gmfKwdCod = ['GmfReserved1',
            'GmfVersionFormatted',
            'GmfReserved2',
            'GmfDimension',
            'GmfVertices',
            'GmfEdges',
            'GmfTriangles',
            'GmfQuadrilaterals',
            'GmfTetrahedra',
            'GmfPentahedra',
            'GmfHexahedra',
            'GmfReserved3',
            'GmfReserved4',
            'GmfCorners',
            'GmfRidges',
            'GmfRequiredVertices',
            'GmfRequiredEdges',
            'GmfRequiredTriangles',
            'GmfRequiredQuadrilaterals',
            'GmfTangentAtEdgeVertices',
            'GmfNormalAtVertices',
            'GmfNormalAtTriangleVertices',
            'GmfNormalAtQuadrilateralVertices',
            'GmfAngleOfCornerBound',
            'GmfReserved5',
            'GmfReserved6',
            'GmfReserved7',
            'GmfReserved8',
            'GmfReserved9',
            'GmfReserved10',
            'GmfReserved11',
            'GmfReserved12',
            'GmfReserved13',
            'GmfReserved14',
            'GmfReserved15',
            'GmfReserved16',
            'GmfReserved17',
            'GmfReserved18',
            'GmfReserved19',
            'GmfReserved20',
            'GmfReserved21',
            'GmfReserved22',
            'GmfReserved23',
            'GmfReserved24',
            'GmfReserved25',
            'GmfReserved26',
            'GmfReserved27',
            'GmfReserved28',
            'GmfReserved29',
            'GmfReserved30',
            'GmfBoundingBox',
            'GmfReserved31',
            'GmfReserved32',
            'GmfReserved33',
            'GmfEnd',
            'GmfReserved34',
            'GmfReserved35',
            'GmfReserved36',
            'GmfReserved37',
            'GmfTangents',
            'GmfNormals',
            'GmfTangentAtVertices',
            'GmfSolAtVertices',
            'GmfSolAtEdges',
            'GmfSolAtTriangles',
            'GmfSolAtQuadrilaterals',
            'GmfSolAtTetrahedra',
            'GmfSolAtPentahedra',
            'GmfSolAtHexahedra',
            'GmfDSolAtVertices',
            'GmfISolAtVertices',
            'GmfISolAtEdges',
            'GmfISolAtTriangles',
            'GmfISolAtQuadrilaterals',
            'GmfISolAtTetrahedra',
            'GmfISolAtPentahedra',
            'GmfISolAtHexahedra',
            'GmfIterations',
            'GmfTime',
            'GmfReserved38']

indicesGmf = dict([(value, i) for i, value in enumerate(gmfKwdCod)])


tclock = dict()


def tic(ref=0):
    global tclock
    tclock[ref] = time.time()


def toc(ref=0):
    global tclock
    return format(time.time()-tclock[ref], "0.2f")+"s"