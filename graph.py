import math
import numpy as np
import copy
import sys
import time
import matplotlib.pyplot as plt
import rmsd

file = open( "outputs.txt", "r" )
content = file.readlines()
file.close()

scores_over_time = []
rmsds_over_time = []

for value in content:
        value = value.split()
        #print value[1], value[2]
        scores_over_time.append( value[1] )
        rmsds_over_time.append( value[2] )

#print content

fig, ax1 = plt.subplots()
ax1.plot(scores_over_time, 'b-')
ax1.set_title('Energy x RMSD' )
ax1.set_xlabel('Iterations')
# Make the y-axis label, ticks and tick labels match the line color.
ax1.set_ylabel('Energy', color='r')
ax1.tick_params('y', colors='r')
ax2 = ax1.twinx()
ax2.plot(rmsds_over_time, 'r-')
ax2.set_ylabel('RMSD', color='b')
ax2.tick_params('y', colors='b')
fig = plt.gcf()
print("###############################")
fig.savefig( "rmsd_energy.png", dpi=100, bbox_inches='tight')