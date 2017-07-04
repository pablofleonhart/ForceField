import math
import numpy as np
import os
import sys

class PDBAligner:

    def transform( self, modPosAtoms, translation, rotation ):
        tempX = np.zeros( ( len( modPosAtoms ), 3 ) )
        tempY = np.zeros( ( len( modPosAtoms ), 3 ) )
        solution = np.zeros( ( len( modPosAtoms ), 3 ) )

        rotationX = [ [1.0, 0.0, 0.0], [0.0, math.cos( rotation[0] ), -math.sin( rotation[0] )], [0.0, math.sin( rotation[0] ), math.cos( rotation[0] )] ]
        rotationY = [ [math.cos( rotation[1] ), 0.0, math.sin( rotation[1] )], [0.0, 1.0, 0.0], [-math.sin( rotation[1] ), 0.0, math.cos( rotation[1] )] ]
        rotationZ = [ [math.cos( rotation[2] ), -math.sin( rotation[2] ), 0.0], [math.sin( rotation[2] ), math.cos( rotation[2] ), 0.0], [0.0, 0.0, 1.0] ]

        # rotation
        for i in range( len( modPosAtoms ) ):
            for j in range( 3 ):
                for k in range( 3 ):
                    tempX[i][j] += modPosAtoms[i][k] * rotationX[k][j]

            for j in range( 3 ):
                for k in range( 3 ):
                    tempY[i][j] += tempX[i][k] * rotationY[k][j]

            for j in range( 3 ):
                for k in range( 3 ):
                    solution[i][j] += tempY[i][k] * rotationZ[k][j]

        # translation
        for i in range( len( modPosAtoms ) ):
            for j in range( 3 ):
                solution[i][j] += translation[j]

        return solution

    def calcRMSD( self, reference, solution ):
        sumDistance = 0

        for i in range( len( reference ) ):
            sumDistance += math.pow( reference[i][0] - solution[i][0], 2 )
            sumDistance += math.pow( reference[i][1] - solution[i][1], 2 )
            sumDistance += math.pow( reference[i][2] - solution[i][2], 2 )

        return math.sqrt( sumDistance/2.0 )