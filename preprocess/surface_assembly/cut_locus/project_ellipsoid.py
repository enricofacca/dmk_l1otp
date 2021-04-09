#!pyhton
#
# project the point in gridin to the surface defined by
# (x/a)^2+(y/b)^+(z/c)^2=1 
#


# usage pyhton "gridin" "gridout" x y z

import sys
sys.path.append('../')
import meshtools_surface as mt
import numpy as np

a=0.2
b=0.6
c=1.0


gridin=sys.argv[1]
gridout=sys.argv[2]

coord, triang, flag=mt.read_grid(gridin)

newcoord=coord
for inode in range(coord.shape[0]):
    p=coord[inode][:]
    x=(p[0]/a)**2+(p[1]/b)**2+(p[2]/c)**2
    newcoord[inode][:]=(1.0/x)**(0.5) * coord[inode][:]
    
mt.write_grid(newcoord,triang,gridout,'dat')
mt.write_grid(newcoord,triang,'ellipsoid.vtk','vtk')
