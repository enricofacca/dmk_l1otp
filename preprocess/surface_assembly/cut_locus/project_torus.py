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

R=2.0
r=1.0

gridin=sys.argv[1]
gridout=sys.argv[2]

coord, triang, flag=mt.read_grid(gridin)

newcoord=coord
for inode in range(coord.shape[0]):
    x=coord[inode][0]
    y=coord[inode][1]
    z=coord[inode][2]
    #print(inode,'node')
    #print(x,y,z)
    rm=np.sqrt(x**2+y**2) 
    #print ( 'err', (rm - R)**2 +z**2-r**2 )
    phi=np.arctan2(y,x)
    theta = np.arctan2(z,(rm-R))
       

    
    #print('angles',phi*180/(np.pi),theta*180/(np.pi))
    new_x=(R+r*np.cos(theta))*np.cos(phi)
    new_y=(R+r*np.cos(theta))*np.sin(phi)
    new_z=r*np.sin(theta)
    #print(inode,'new node')
    #print([new_x,new_y,new_z])
    #print('differnce')
    #print (np.linalg.norm([x-new_x,y-new_y,z-new_z]))
    rm=np.sqrt(new_x**2+new_y**2) 
    #print ( 'new err', (rm - R)**2 +new_z**2-r**2 )
    #print ( ' ')
    newcoord[inode][:]=[new_x,new_y,new_z]

    
    
mt.write_grid(newcoord,triang,gridout,'dat')
mt.write_grid(newcoord,triang,'ellipsoid.vtk','vtk')
