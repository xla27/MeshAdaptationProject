import os, sys, shutil
import time, copy
import numpy as np
from SU2.run.interface import CFD as SU2_CFD
from SU2.io import redirect
from geometry.CMesh import CMesh
from utilities.input_output import WriteParamFile


class CDriver():

    def __init__(self, sensors, meshFilename, solFilename, params):
        self.sensors = sensors
        self.nSensors = len(sensors)
        self.meshFilename = meshFilename
        self.solFilename = solFilename
        self.params = params
        return
    
    def ComputeSU2Gradients(self, configCfd):

        # performing a single CFD iterations to compute the adaptation sensors gradients
        konfigCfd = copy.deepcopy(configCfd)

        konfigCfd['ITER'] = 1
        konfigCfd['SOLUTION_FILENAME'] = konfigCfd['RESTART_FILENAME']
        ext = '.dat' if konfigCfd['RESTART_FILENAME'].endswith('dat') else '.csv'

        for iSensor, sensor in enumerate(self.sensors):
            konfigCfd['ADAP_SENSOR'] = sensor
            konfigCfd['RESTART_FILENAME'] = f'restart_{sensor}{ext}'
            konfigCfd['CONV_FILENAME'] = f'history_{sensor}'

            with redirect.output(f'su2_{sensor}.out'): SU2_CFD(konfigCfd)

            os.remove(konfigCfd['CONV_FILENAME']+'.csv')

        return
    
    def ReadSU2(self, configCfd):

        self.mesh = CMesh()

        # reading mesh .su2
        self.mesh.ReadMeshSU2(self.meshFilename)

        # finalizing data structure
        self.mesh.FinalizingDataStructure(self.sensors)
        self.params['card'] = self.mesh.cardinality

        # reading the solution
        for iSensor, sensor in enumerate(self.sensors):
            ext = '.dat' if configCfd['RESTART_FILENAME'].endswith('dat') else '.csv'
            solFilename = f'restart_{sensor}{ext}'
            self.mesh.ReadSolSU2(sensor, iSensor, solFilename)

        # setting the mesh diameter as param
        self.params['diam'] = self.mesh.diameter

    def ComputeMetricAndAnisoError(self):

        # mesh data structure
        meshDict = self.mesh.GetMeshDict()

        dim = self.mesh.GetDim()

        keyElem = 'Triangles' if dim == 2 else 'Tetrahedra'

        # print('\tStart metric computation.')

        time_total_init = time.time()

        # element-wise operations on single elements
        for element in meshDict[keyElem]:

            # assigning vertices coordinates and gradients 
            element.SetVerticesCoordinatesAndGradient()

            # computing lambdak
            element.ComputeLambdak(computeRk=True)

            # compute element volume
            element.ComputeVolume()

            for iSensor in range(self.nSensors):

                # computing the gradient on the element
                element.ComputeGradient(iSensor=iSensor)

                # computing local error contribution
                element.ComputeLocalErrorContribution(iSensor=iSensor)


        # element-wise operations on patches
        limitedElements  = np.zeros(self.nSensors, dtype=int)
        globalAnisoError = np.zeros(self.nSensors)
        for element in meshDict[keyElem]:

            # creating the element patch
            element.SetPatchID()
            element.SetPatchElements(self.mesh)

            # computing the patch volume
            element.ComputePatchVolume()

            # computing the element-wise metric
            for iSensor, sensor in enumerate(self.sensors):
                limited, localAnisoError = element.ComputeMetricAndAnisoError(iSensor=iSensor,
                                                                            toll=self.params['toll'][iSensor], 
                                                                            diam=self.params['diam'], 
                                                                            card=self.params['card'])        
                limitedElements[iSensor]  += limited
                globalAnisoError[iSensor] += localAnisoError

            # intersect the metric
            element.IntersectMetric()

        # computing the vertex-wise metric
        for vertex in meshDict['Vertices']:

            vertex.ComputeMetric(self.mesh)

        # print('\tEnd metric computation.')

        return np.sqrt(globalAnisoError), limitedElements, time.time()-time_total_init

    def WriteMedit(self, meditFilename, solFilename):

        # writing the mesh in medit format
        self.mesh.WriteMeshMedit(meditFilename)

        # writing the solution in sol format
        self.mesh.WriteSolMedit(solFilename)

        # write mmg parameters file
        WriteParamFile(self.mesh, self.params, meditFilename)

    def ReadMedit(self, meditFilename):
        self.mesh.ReadMeshMedit(meditFilename)

    def WriteSU2(self, su2Filename):

        self.mesh.WriteMeshSU2(su2Filename)




    


