import numpy as np

from CMesh import CMesh
import time

class CDriver():

    def __init__(self, sensor, meshFilename, solFilename, params):
        self.sensor = sensor
        self.meshFilename = meshFilename
        self.solFilename = solFilename
        self.params = params
        return
    
    def ReadSU2(self):

        self.mesh = CMesh()

        # reading mesh .su2
        self.mesh.ReadMeshSU2(self.meshFilename)

        # reading the solution
        self.mesh.ReadSolSU2(self.sensor, self.solFilename)

        # setting the mesh diameter as param
        self.params['diam'] = self.mesh.diameter

        # finalizing data structure
        self.mesh.FinalizingDataStructure()
        self.params['card'] = self.mesh.cardinality

    def ComputeMetric(self):

        # mesh data structure
        meshDict = self.mesh.GetMeshDict()

        dim = self.mesh.GetDim()

        keyElem = 'Triangles' if dim == 2 else 'Tetrahedra'

        print('\tStart metric computation.')

        # computing the element-wise metric
        limited = 0
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

            # computing the element-wise metric
            limited += element.ComputeMetric(toll=self.params['toll'], 
                                            diam=self.params['diam'], 
                                            card=self.params['card'])
            
        print('\t\tAspect ratio limited by gmin in %i out of %i elements' %
              (limited, len(meshDict[keyElem])))

        # computing the vertex-wise metric
        for vertex in meshDict['Vertices']:

            vertex.ComputeMetric(self.mesh)

        print('\tEnd metric computation.')

    def WriteMedit(self, meditFilename, solFilename):

        # writing the mesh in medit format
        self.mesh.WriteMeshMedit(meditFilename)

        # writing the solution in sol format
        self.mesh.WriteSolMedit(solFilename)

        # write mmg parameters file
        self.mesh.WriteParamFile(self.params, meditFilename)

    def ReadMedit(self, meditFilename):
        self.mesh.ReadMeshMedit(meditFilename)

    def WriteSU2(self, su2Filename):

        self.mesh.WriteMeshSU2(su2Filename)




    


