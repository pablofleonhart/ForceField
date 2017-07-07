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
from funcEnergy import EnergyFunction

def evals( acor, c ):
    return acor.evaluator( c )

class ACOR:
    pdbPattern = "{:6s}{:5d} {:^4s}{:1s}{:3s} {:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}{:2s}"
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

    def __init__( self, exp, mod, variables, maximization, iterations ):
        self.experimental = exp
        self.modified = mod
        self.parameters = variables
        self.numVar = len( variables )
        self.maximize = maximization
        self.maxIterations = iterations
        self.upperBound = [1] * self.numVar
        self.lowerBound = [0] * self.numVar
        self.generations = []
        self.energies = []
        self.rmsds = []
        self.mod = []

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
        rotation = [rotation[i:i+2] for i in range( 0, len( rotation ), 2 )]

        app = AminoPhiPsi( "1L2Y-P.pdb" )        
        app.pdb.adjustAtoms( self.experimental.atoms, self.experimental.aminoAcids )
        #print app.pdb.posAtoms
        app.adjustPhiPsi( rotation )

        '''print app.pdb.atoms
        fitness = self.calcKabschRMSD( app.pdb.posAtoms )'''
        fe = EnergyFunction( app.pdb )
        fitness = fe.getEnergy()
        #print fitness
        return fitness

    def multiprocessEvaluator( self, x ):
        nprocs = 8
        pool = multiprocessing.Pool( processes = nprocs )
        results = [pool.apply_async( evals, [self, c] ) for c in x]
        pool.close()
        pool.join()
        
        return [r.get() for r in results]

    def initialize( self ):
        return np.random.uniform( low = 0, high = 1, size = ( self.sizeSolutions, self.numVar ) )

    def evolve( self ):
        self.generations = []
        self.energies = []
        self.rmsds = []
        step = 20

        solutions = np.zeros( ( self.sizeSolutions, self.numVar ) )
        mFitness = np.zeros( ( self.sizeSolutions, 1 ) )

        print '-----------------------------------------'
        print 'Starting initilization of solution matrix'

        initialSolution = self.initialize()
        vFitness = self.multiprocessEvaluator( initialSolution )

        for i in range( len( vFitness ) ):
            mFitness[i] = vFitness[i]

        solutions = np.hstack( ( initialSolution, mFitness ) )
        solutions = sorted( solutions, key = lambda row: row[-1], reverse = self.maximize )
        solutions = np.array( solutions )

        weights = np.zeros( ( self.sizeSolutions ) )
        for i in range( self.sizeSolutions ):
            weights[i] = ( 1/( self.qk * math.sqrt( 2 * math.pi ) ) ) * math.exp( -math.pow( i, 2 )/( 2 * math.pow( self.q, 2 ) * math.pow( self.sizeSolutions, 2 ) ) )

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

            Stemp = np.zeros( ( self.numAnts, self.numVar ) )

            # for each ant..
            for ant in range( self.numAnts ):
                # ..it's choosed a solution
                cs = np.random.random_sample()
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

                    # calc 'i' value with gaussian functions
                    x = np.random.random_sample()
                    gi = weights[sol]*math.exp( -math.pow( x - solutions[sol][i], 2 ) / ( 2*math.pow( sigma, 2 ) ) )* (1/( sigma*math.pow( 2*math.pi, 2 ) ))

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

            # keep best after each iteration
            best_par.append(solutions[0][:self.numVar])
            best_obj.append(solutions[0][-1])
            best_res.append(solutions[0][self.numVar:-1])
            best_sol.append(solutions[0][:])
            worst_obj.append(solutions[-1][-1])

            #print "Best individual:", self.parameters
            #print best_sol[0][0:len( self.parameters )]
            print "Fitness:", solutions[0][:][len( self.parameters )]
            self.generations.append( iterations )
            self.energies.append( solutions[0][:][len( self.parameters )] )            

            rotation = [ (2*math.pi*i)-math.pi for i in best_sol[0][0:len( self.parameters )] ]
            rotation.append( 0.0 )

            rt = []
            rt.append( 0.0 )
            for i in xrange( len( rotation ) ):
                rt.append( rotation[i] )

            rotation = [rt[i:i+2] for i in range( 0, len( rt ), 2 )]

            app = AminoPhiPsi( "1L2Y-P.pdb" )
            app.pdb.adjustAtoms( self.experimental.atoms, self.experimental.aminoAcids )
            #print rotation
            app.adjustPhiPsi( rotation )
            rm = self.calcKabschRMSD( app.pdb.posAtoms )
            print "RMSD", rm
            self.rmsds.append( rm )

            if iterations%step == 0:                
                mod = app.pdb.posAtoms

                pdbNew = open( "1L2Y-F" + str( iterations ) + ".pdb", "w" )
                countTotal = 1
                acid = 0
                aa = None
                for z in range( 0, len( self.modified.atoms ) ):
                    if self.modified.aminoAcids[z] != aa:
                        aa = self.modified.aminoAcids[z]
                        acid += 1
                    pdbNew.write( self.pdbPattern.format( "ATOM", countTotal, str( self.modified.atoms[z] ), " ", str( self.modified.aAcids[z] ), " ", \
                                  acid, " ", float( mod[z][0] ), float( mod[z][1] ), float( mod[z][2] ), float( 1.00 ), float( 0.0 ), "", "" ) + "\n" )

                    countTotal += 1

                pdbNew.write( "TER\n" )
                pdbNew.close() 

            iterations += 1

        best_sol = sorted( best_sol, key=lambda row: row[-1], reverse = self.maximize )

        print "Best individual:", self.parameters
        print best_sol[0][0:len( self.parameters )]
        print "Fitness:", best_sol[0][-1]

        print self.generations
        print self.energies
        print self.rmsds

        ff = open( "outputs.txt", "w" )

        for i in xrange( len( self.generations ) ):
            ff.write( str( self.generations[i] ) + " " + str( self.energies[i] ) + " " + str( self.rmsds[i] ) + "\n" )
        
        ff.close()

        rotation = [ (2*math.pi*i)-math.pi for i in best_sol[0][0:len( self.parameters )] ]
        rotation.append( 0.0 )

        rt = []
        rt.append( 0.0 )
        for i in xrange( len( rotation ) ):
            rt.append( rotation[i] )

        rotation = [rt[i:i+2] for i in range( 0, len( rt ), 2 )]

        app = AminoPhiPsi( "1L2Y-P.pdb" )
        app.pdb.adjustAtoms( self.experimental.atoms, self.experimental.aminoAcids )
        #print rotation
        app.adjustPhiPsi( rotation )
        mod = app.pdb.posAtoms

        pdbNew = open( "1L2Y-F.pdb", "w" )
        countTotal = 1
        acid = 0
        aa = None
        for z in range( 0, len( self.modified.atoms ) ):
            if self.modified.aminoAcids[z] != aa:
                aa = self.modified.aminoAcids[z]
                acid += 1
            pdbNew.write( self.pdbPattern.format( "ATOM", countTotal, str( self.modified.atoms[z] ), " ", str( self.modified.aAcids[z] ), " ", \
                          acid, " ", float( mod[z][0] ), float( mod[z][1] ), float( mod[z][2] ), float( 1.00 ), float( 0.0 ), "", "" ) + "\n" )

            countTotal += 1

        pdbNew.write( "TER\n" )
        pdbNew.close()        