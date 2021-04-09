#!pyhton
#
# replace the closest point to p=[x,y,z] wiht x,y,z 
#

# usage pyhton "gridin" "gridout" x y z

import sys
sys.path.append('../')
import meshtools_surface as mt
import numpy as np

gridin=sys.argv[1]
gridout=sys.argv[2]
point=sys.argv[3:]  # in form x y z

print point
pdirac=np.asarray(point,dtype=np.float)
coord, triang, flag=mt.read_grid(gridin)

newcoord=coord
imin=0
dist_min=1.0e30
for inode in range(coord.shape[0]):
    p=coord[inode][:]
    dist=((p[0]-pdirac[0])**2 +
          (p[1]-pdirac[1])**2 +
          (p[2]-pdirac[2])**2)
    if (dist<dist_min) : 
        imin=inode
    
print(imin+1)
#newcoord[imin][:]=pdirac[:]

#mt.write_grid(newcoord,triang,gridout,'dat')
#mt.write_grid(newcoord,triang,'dirac.vtk','vtk')
