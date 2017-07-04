import os
import sys
import shutil
import numpy as np
import random
import multiprocessing
import datetime
import math
from scipy.stats import norm
import copy
from pdbAligner import PDBAligner
from aminoPhiPsi import AminoPhiPsi
import rmsd

def evals( acor, c ):
    return acor.evaluator( c )

class ACOR:
    pdbPattern = "{:6s}{:5d} {:^4s}{:1s}{:3s} {:1s}{:4d}{:1s} {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}"
    NHC_ATOMS = ("N", "H", "1H", "H1", "2H", "H2", "3H", "H3", "CA")
    generations = []
    values = []
    mod = []
    experimental = None
    modified = None
    # maximization or minimization problem
    maximize = False

    # variables
    parameters = None
    # number of variables
    numVar = 0
    # size of solution archive
    sizeSolutions = 200
    # number of ants
    numAnts = 150
    # parameter self.q
    q = 0.0001
    # standard deviation
    qk = q * sizeSolutions
    # parameter self.xi (like pheromone evaporation)
    xi = 0.85
    # maximum iterations
    maxIterations = 100
    # bounds of variables
    upperBound = []
    lowerBound = []

    def __init__( self, exp, mod, variables, maximization, iterations ):
        self.experimental = exp
        self.modified = mod
        self.parameters = variables
        self.numVar = len( variables )
        self.maximize = maximization
        self.maxIterations = iterations
        self.upperBound = [1] * self.numVar
        self.lowerBound = [0] * self.numVar

    def rotate_to( self, ang, atoms, aminoAcids, posAtoms ):
        angles = [math.pi]*len( ang )
        angles[0] = math.pi*2
        angles[len( ang ) - 1] = math.pi*2
        n_aa = len( aminoAcids )
        for i in xrange(n_aa):
            if i + min( aminoAcids) <= max(aminoAcids):
                #ROTATE PHI
                #print atoms, aminoAcids
                n_i = zip(atoms, aminoAcids).index((" N  ", i + min(aminoAcids)))   
                ca_i = zip(atoms, aminoAcids).index((" CA ", i + min(aminoAcids)))
                current_angles = angles
                #print current_angles
                dphi = math.atan2(math.sin(ang[2*i] - current_angles[2*i]), math.cos(ang[2*i] - current_angles[2*i]))
                #print "dphi", degrees( dphi )
                n_pos = posAtoms[n_i]
                ca_pos = posAtoms[ca_i]                
                ia = 0
                for atom in zip(atoms, aminoAcids):
                    if (i > 0) and (atom[1] > i + min(aminoAcids) or (atom[1] == i + min(aminoAcids) and (atom[0].strip() not in self.NHC_ATOMS))): 
                        posAtoms[ia] = self.rotate_atom_around_bond(dphi, posAtoms[ia], n_pos, ca_pos)
                        #print(atom[0], atom[1])   
                    ia += 1        
                #ROTATE PSI    
                c_i  = zip(atoms, aminoAcids).index((" C  ",  i + min(aminoAcids)))  
                ca_i = zip(atoms, aminoAcids).index((" CA ", i + min(aminoAcids)))
                current_angles = angles
                #print current_angles
                dpsi = math.atan2(math.sin(ang[2*i+1] - current_angles[2*i+1]), math.cos(ang[2*i+1] - current_angles[2*i+1]))              
                c_pos = posAtoms[c_i] 
                ca_pos = posAtoms[ca_i]
                ia = 0
                for atom in zip(atoms, aminoAcids):
                    if (i+min(aminoAcids) < max(aminoAcids)) and (atom[1] > i+min(aminoAcids) or (atom[1] == i+min(aminoAcids) and (atom[0]==" O  "))): 
                        posAtoms[ia] = self.rotate_atom_around_bond(dpsi, posAtoms[ia], ca_pos, c_pos)
                        #print(atom[0], atom[1])          
                ia += 1
                
    def normalize( self, v ):
        norm = np.linalg.norm( v )
        if norm == 0: 
            return v
        return v/norm  

    def rotate_atom_around_bond( self, theta, atom_pos, bond_start, bond_end ):
        #https://en.wikipedia.org/wiki/Rodrigues%27_rotation_formula
        v = np.array( atom_pos ) - np.array( bond_start )
        ant = np.array( bond_end ) - np.array( bond_start )
        ant = self.normalize( ant )
        rot_pos = v * np.cos( theta ) + ( np.cross( ant, v ) ) * np.sin( theta ) + ant * ( np.dot( ant, v ) ) * ( 1.0 - np.cos( theta ) )
        return list( rot_pos + np.array( bond_start ) )

    def calcKabschRMSD( self, mod ):
        P = np.array( self.experimental.posAtoms )
        self.q = np.array( mod )
        P -= rmsd.centroid( P )
        self.q -= rmsd.centroid( self.q )
        result = rmsd.kabsch_rmsd( P, self.q )
        return result

    def evaluator( self, x ):
        rotation = [ ( 2 * math.pi * i ) - math.pi for i in x ]
        rotation = np.hstack( ( 0.0, rotation, 0.0 ) )
        print "r0",rotation
        [rotation[i:i+2] for i in range( 0, len( rotation ), 2 )]
        print "r",rotation

        mod = copy.deepcopy( self.modified.posAtoms )

        app = AminoPhiPsi( "1L2Y-P.pdb" )
        app.adjustPhiPsi( rotation )

        fitness = self.calcKabschRMSD( mod )
        return fitness

    def multiprocessEvaluator( self, x ):
        nprocs = 4
        pool = multiprocessing.Pool( processes = nprocs )
        results = [pool.apply_async( evals, [self, c] ) for c in x]
        pool.close()
        pool.join()
        
        return [r.get() for r in results]

    def initialize( self ):
        return np.random.uniform( low = 0, high = 1, size = ( self.sizeSolutions, self.numVar ) )

    def evolve( self ):
        self.generations = []
        self.values = []

        # initilize matrices
        solutions = np.zeros( ( self.sizeSolutions, self.numVar ) )
        mFitness = np.zeros( ( self.sizeSolutions, 1 ) )

        print '-----------------------------------------'
        print 'Starting initilization of solution matrix'
        print '-----------------------------------------'

        initialSolution = self.initialize()
        vFitness = self.multiprocessEvaluator( initialSolution )

        for i in range( len( vFitness ) ):
            mFitness[i] = vFitness[i]

        solutions = np.hstack( ( initialSolution, mFitness ) )
        # sort according to fitness (last column)
        solutions = sorted( solutions, key = lambda row: row[-1], reverse = self.maximize )
        solutions = np.array( solutions )

        # initilize weight array with pdf function
        weights = np.zeros( ( self.sizeSolutions ) )
        for i in range( self.sizeSolutions ):
            weights[i] = ( 1/( self.qk * math.sqrt( 2 * math.pi ) ) ) * math.exp( -math.pow( i, 2 )/( 2 * math.pow( self.q, 2 ) * math.pow( self.sizeSolutions, 2 ) ) )

        # initialize variables
        iterations = 1
        best_par = []
        best_obj = []
        best_sol = []
        best_res = []
        worst_obj = []
        best_par.append(solutions[0][:self.numVar])
        best_obj.append(solutions[0][-1])
        best_sol.append(solutions[0][:])
        best_res.append(solutions[0][self.numVar:-1])
        worst_obj.append(solutions[-1][-1])

        p = weights/sum( weights )
        stop = 0

        while iterations <= self.maxIterations:
            print '-----------------------------------------'
            print 'Iteration', iterations
            print '-----------------------------------------'

            Stemp = np.zeros( ( self.numAnts, self.numVar ) )

            # for each ant..
            for ant in range( self.numAnts ):
                # ..it's choosed a solution
                cs = np.random.random_sample()
                #print cs
                total = 0
                for z in xrange( self.sizeSolutions-1, -1, -1 ):
                    total += p[z]
                    if cs <= total:
                        sol = z
                        break

                # for each variable..
                for i in range( self.numVar ):
                    # ..it's calc standard deviation of 'sol' solution
                    sigma = 0
                    for y in xrange( self.sizeSolutions ):
                        sigma += abs( solutions[y][i] - solutions[sol][i] )/( self.sizeSolutions-1 )

                    #print "sigma", i, "=", sigma

                    # calc 'i' value with gaussian functions
                    x = np.random.random_sample()
                    gi = weights[sol]*math.exp( -math.pow( x - solutions[sol][i], 2 ) / ( 2*math.pow( sigma, 2 ) ) )* (1/( sigma*math.pow( 2*math.pi, 2 ) ))

                    #print gi
                    #Stemp[ant][i] = sigma[i] * x + solutions[selection][i]
                    Stemp[ant][i] = gi
                    if Stemp[ant][i] > self.upperBound[i]:
                        Stemp[ant][i] = self.upperBound[i]
                    elif Stemp[ant][i] < self.lowerBound[i]:
                        Stemp[ant][i] = self.lowerBound[i]

            vFitness = self.multiprocessEvaluator( Stemp )
            mFitness = np.zeros( ( self.numAnts, 1 ) )

            for i in range( len( vFitness ) ):
                mFitness[i] = vFitness[i]

            # add responses and "fitness" column to solution
            Ssample = np.hstack( ( Stemp, mFitness ) )

            # add new solutions in the solutions table
            Solution_temp = np.vstack( ( solutions, Ssample ) )

            # sort according to "fitness"
            Solution_temp = sorted( Solution_temp, key = lambda row: row[-1], reverse = self.maximize )
            Solution_temp = np.array( Solution_temp )

            # keep best solutions
            solutions = Solution_temp[:self.sizeSolutions][:]

            #print "solutions", solutions
            # keep best after each iteration
            best_par.append(solutions[0][:self.numVar])
            best_obj.append(solutions[0][-1])
            best_res.append(solutions[0][self.numVar:-1])
            best_sol.append(solutions[0][:])
            worst_obj.append(solutions[-1][-1])

            #print "Best individual:", self.parameters
            #print best_sol[0][0:len(self.parameters)]
            print "Fitness:", solutions[0][:][8]
            self.generations.append( iterations )
            self.values.append( solutions[0][:][8] )

            iterations += 1

        best_sol = sorted( best_sol, key=lambda row: row[-1], reverse = self.maximize )

        print "Best individual:", self.parameters
        print best_sol[0][0:len(self.parameters)]
        print "Fitness:"
        print best_sol[0][-1]

        print self.generations
        print self.values

        rotation = [ (2*math.pi*i)-math.pi for i in best_sol[0][0:len(self.parameters)] ]
        rotation.append( 0.0 )

        rt = []
        rt.append( 0.0 )
        for i in xrange( len( rotation ) ):
            rt.append( rotation[i] )

        mod = copy.deepcopy( self.modified.posAtoms )

        self.rotate_to( rt, self.modified.atoms, self.modified.aminoAcids, mod )

        pdbNew = open( "1PLX-vFitness.pdb", "weights" )
        countTotal = 1
        acid = 0
        aa = None
        for z in range( 0, len( self.modified.atoms ) ):
            if self.modified.aminoAcids[z] != aa:
                aa = self.modified.aminoAcids[z]
                acid += 1
            pdbNew.write( self.pdbPattern.format( "ATOM", countTotal, str( self.modified.atoms[z] ), " ", str( self.modified.aAcids[z] ), " ", \
                          acid, " ", float( mod[z][0] ), float( mod[z][1] ), float( mod[z][2] ), float( 1.00 ), float( 0.0 ) ) + "\n" )

            countTotal += 1

        pdbNew.write( "TER\n" )
        pdbNew.close()