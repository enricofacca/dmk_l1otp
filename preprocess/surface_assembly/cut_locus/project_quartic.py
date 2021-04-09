#!pyhton
#
# project the point in gridin to the surface defined by
# x^4+y^4+z^4=1 
#

# usage pyhton "gridin" "gridout" x y z

import sys
sys.path.append('../')
import meshtools_surface as mt
import numpy as np

gridin=sys.argv[1]
gridout=sys.argv[2]

coord, triang, flag=mt.read_grid(gridin)

newcoord=coord
for inode in range(coord.shape[0]):
    p=coord[inode][:]
    a=p[0]**4+p[1]**4+p[2]**4
    newcoord[inode][:]=(1.0/a)**(0.25) * coord[inode][:]
    
mt.write_grid(newcoord,triang,gridout,'dat')
mt.write_grid(newcoord,triang,'quartic.vtk','vtk')
