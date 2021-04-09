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

coord, triang, flag=mt.read_grid(gridin)

newcoord=coord
for inode in range(coord.shape[0]):
    p=coord[inode][:]
    if ( (p[2]>-0.9) and (p[2]<0.9)):
         l=np.sqrt(( 1.0 - p[2]**2 )/( p[0]**2+p[1]**2 ))
         if ( math.isnan(l) ):
             print inode, 'inside',l
         
         newcoord[inode][0] = l * coord[inode][0]
         newcoord[inode][1] = l * coord[inode][1]
         newcoord[inode][2] = coord[inode][2]
    else:
         l=(p[0])**2+(p[1])**2+(p[2])**2
         if ( math.isnan(l) ):
             print inode, 'inside',l
         newcoord[inode][:]=(1.0/l)**(0.5) * coord[inode][:]

    #print newcoord[inode][0]**2+ newcoord[inode][1]**2+newcoord[inode][2]**2
         

mt.write_grid(newcoord,triang,gridout,'dat')
mt.write_grid(newcoord,triang,'sphere_prj.vtk','vtk')
