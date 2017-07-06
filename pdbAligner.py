import math
import numpy as np
import os
import sys

class PDBAligner:

    def transform( self, transformation, modAtoms ):
        translation = np.matrix( [transformation[0:3]]*len( modAtoms ) )
        rotation = transformation[3:6]

        rotationX = np.matrix( [[1.0, 0.0, 0.0], [0.0, math.cos( rotation[0] ), -math.sin( rotation[0] )], [0.0, math.sin( rotation[0] ), math.cos( rotation[0] )]] )
        rotationY = np.matrix( [[math.cos( rotation[1] ), 0.0, math.sin( rotation[1] )], [0.0, 1.0, 0.0], [-math.sin( rotation[1] ), 0.0, math.cos( rotation[1] )]] )
        rotationZ = np.matrix( [[math.cos( rotation[2] ), -math.sin( rotation[2] ), 0.0], [math.sin( rotation[2] ), math.cos( rotation[2] ), 0.0], [0.0, 0.0, 1.0]] )
        rotationXYZ = rotationZ * rotationY * rotationX

        transformedAtoms = np.matrix(copy.deepcopy( modAtoms ) )
        transformedAtoms = transformedAtoms + translation
        transformedAtoms = transformedAtoms * rotationXYZ.transpose()
       
        transformedAtoms = np.matrix.tolist( transformedAtoms )
        return transformedAtoms

    def calcRMSD( self, refAtoms, modAtoms ):
        #transformedAtoms = self.align( transformation, modAtoms )
        distanceSum = 0.0

        for coord in zip( refAtoms, modAtoms ):
            distanceSum += ( coord[0][0] - coord[1][0] )**2
            distanceSum += ( coord[0][1] - coord[1][1] )**2 
            distanceSum += ( coord[0][2] - coord[1][2] )**2

        score = math.sqrt( distanceSum / float( len( refAtoms ) ) )
        return score