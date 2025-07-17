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
    
    def ReadSU2(self, meshFilename, solFilename):

        mesh = CMesh()

        # reading mesh .su2
        mesh.ReadMeshSU2(meshFilename)

        #