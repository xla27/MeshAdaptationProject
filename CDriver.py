# Struttura codice

# lettura mesh .su2
# 1) creazione di Cvertex per ogni vertex contente coordinate, punti neighbours e elementi neighbours 
# 2) creazione di Celement per ogni elemento contenente ID e coordinate vertici, ID elementi patch
# 3) conversione in file .mesh

# lettura soluzione dal restart .csv
# 1) loading ad ogni CVertex di valore sensore e valore gradiente sensore (GG)

# main loop sugli elementi 
# 1) interpolazione dai vertici alla cella
# 2) calcolo della G
# 3) calcolo e autodecomposizione della jacobiana
# 4) nuovi autovalori e autovettori
# 5) calcolo metrica

# loop sui vertices per i calcolo della media nodewise
# scrittura file .sol con metrica nodewise
import numpy as np

from CMesh import CMesh
from CElement import CElement
from CVertex import CVertex

class CDriver():

    def __init__(self, parameters):
        return
    
    def Initialize(self, meshFilename, solFilename):

        self.mesh = CMesh()

        # reading mesh .su2
        self.mesh.ReadMeshSU2(meshFilename)

        # reading the solution
        self.mesh.ReadSolSU2(solFilename)

        # finalizing data structure
        self.mesh.FinalizingDataStructure()

    def Run(self):

        # mesh data structure
        meshDict = self.mesh.GetMeshDict()

        dim = self.mesh.GetDim()

        keyElem = 'Triangles' if dim == 2 else 'Tetrahedra'

        # computing the element-wise metric
        for element in meshDict[keyElem]:

            # compute element area
            element.ComputeArea()

            # computing the gradient on the element
            element.ComputeGradient()

            # creating the element patch
            element.SetPatchID()
            element.SetPatchElements(self.mesh)

            # computing the patch area
            element.ComputePatchArea()

            # computing the average patch gradient
            element.ComputePatchAveragedGradient()

            element.ComputeMetric()

        # computing the vertex-wise metric
        for vertex in meshDict['Vertices']:
            
            vertex.ComputeMetric(self.mesh)

    def Finalize(self, meditFilename, solFilename):

        # writing the mesh in medit format
        self.mesh.WriteMeshMedit(meditFilename)

        # writing the solution in sol format
        self.mesh.WriteSolMedit(solFilename)



    


