#!/usr/bin/env python
import re
import numpy as np
import meshtools as mt
import sys
from pyvtk import *
import sys
sys.path.append('../../../globals/python_timedata')
import timedata as td

#
# script to define a timedata a dirac-like forcing term
# intergrated with respect to the p1-galerkin funcition.
# At this stage it only fix the value of the closest node 
# to the value of the dirac function.
#
# Usage:
#  python <grid> <dirac_coord> <dirac_value> <rhs_integrated>
# where :
# grid           (in ) : ASCII triangulation
# dirac_coord    (in ): steady_state time data with dirac coordinates
# dirac_value    (in ): time data with dirac values
# rhs_integrated (out): time data with integrated rhs
#

def dirac2rhs(dirac_values, indeces,rhs):
    rhs[:,:]=0.0
    for i in range(len(indeces)):
        rhs[indeces[i],:] = dirac_values[i,:]
    return rhs;
    

def get_nodes_indeces(dirac_coord, coord):
    #
    # find indeces of closest point
    #
    indeces=[]
    distances=np.zeros(len(coord))
    for i in range(len(dirac_coord)):
        for j in range(len(distances)):
            distances[j] = np.linalg.norm(coord[j,:]-dirac_coord[i,0:2])
        inode=np.argmin(distances)
        indeces.append(inode)
    return indeces;


grid_file=sys.argv[1]
dirac_coord_file=sys.argv[2]
dirac_value_file=sys.argv[3]
rhs_file=sys.argv[4]


#
# read grid coordinate
#
coord, _ , _ = mt.read_grid(grid_file)

#
# read dirac coordinates and get list of the clostest nodes
# 
dirac_coord  = td.read_steady_timedata(dirac_coord_file)
indeces = get_nodes_indeces(dirac_coord, coord)

##### define data variable and outfile
rhs=np.zeros([len(coord),1])
print rhs_file
file_out=open(rhs_file,'w')
td.write_head(file_out,rhs)


#
# read all input lines
#
infile = open(dirac_value_file, 'r')
input_lines = infile.readlines()
infile.close()

#
# read dimensions and initialiaze data
#
ndim,ndata=td.read_head(input_lines)
dirac_data=np.zeros([ndim,ndata])


#
# read sequentially all portion of lines 
#

# initialize counter
iline=1
portion=0
EOF=False

# read first portion and get data
time,nnz,EOF,dirac_data=td.read_portion(input_lines,iline,dirac_data)
iline=iline+2+nnz
portion=portion+1

# set the rhs
rhs=dirac2rhs(dirac_data, indeces,rhs)
td.write_portion(file_out,time,rhs)


# start reading cycle from second portion
while (time < 1e+30):
        time,nnz,EOF,dirac_data=td.read_portion(input_lines,iline,dirac_data)
        print(time)
        iline=iline+2+nnz
        portion=portion+1
        if ( EOF ):
           break

        # set the rhs
        rhs=dirac2rhs(dirac_data, indeces,rhs)
        td.write_portion(file_out,time,rhs)
        if ( time  >= 1.0e30): 
             file_out.write("time " + str(1.e30)+"\n")
             break
             
 
file_out.close()            

             


    
