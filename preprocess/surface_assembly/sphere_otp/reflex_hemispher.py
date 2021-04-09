#!pyhton
#
# project the point in gridin to the surface defined by
# x^2+y^+z^2=1 
# all point -0.8<z<8 are projected to remain on the same parallel
# remaing points are projecte by scaling the coordinates
 


# usage pyhton "gridin" "gridout" x y z

import sys
sys.path.append('../')
import meshtools_surface as mt
import numpy as np
import math

gridin=sys.argv[1]
gridout=sys.argv[2]

coord_north, triang_north, flag_north=mt.read_grid(gridin)

coord_equator = coord_north[abs(coord_north[:,2])<1e-12]
coord_south = coord_north[abs(coord_north[:,2])>1e-12]
coord_south[:,2] = -coord_south[:,2] 
print (len(coord_south))
print (len(coord_north))
print (len(coord_equator))
new_coord=np.concatenate((coord_north, coord_south),axis=0)

nnode=len(coord_north)
ntria=len(triang_north)

trasform=range(len(new_coord))
nnew=0
for i in range(len(coord_north)):
    if ( abs(coord_north[i][2]) > 1e-12 ) :
        trasform[i]=nnode+nnew
        nnew=nnew+1
        #print i, trasform[i]



print ('ntria',ntria)
copy=triang_north

triang_south = np.zeros(np.shape(triang_north),dtype=int)
for itria in range(ntria):
    for iloc in range(3):
        inode=triang_north[itria][iloc]
        #print inode, trasform[inode]
        #print triang_south[itria][iloc], triang_north[itria][iloc],copy[itria][iloc]
        triang_south[itria][iloc]=trasform[inode]
        #print triang_south[itria][iloc], triang_north[itria][iloc],copy[itria][iloc]

new_triang=np.concatenate((copy, triang_south),axis=0)

print (len(new_triang))
print (len(new_coord))

mt.write_grid(new_coord,new_triang,gridout,'dat')
#mt.write_grid(new_coord,new_triang,'sphere_test.vtk','vtk')
