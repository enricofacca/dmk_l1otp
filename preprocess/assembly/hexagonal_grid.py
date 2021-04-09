#!/usr/bin/env python
import numpy as np
import sys
import meshtools as mt

def hexagonal_grid(xcenter=0.0,ycenter=0.0,edge_size=1.0):    
    coord =np.zeros([7,3])
    triang=np.zeros([6,3],dtype=int)

    inode=0
    itria=0
    print(xcenter,ycenter)
    r=float(edge_size)
    x=float(xcenter)
    y=float(ycenter)
    coord[0,:]=[x,y,0]
    for i in range(6): # 0,1, ndivy-1
        inode = inode+1
        coord[inode,:] = [x+r*np.cos(i*np.pi/3),y+r*np.sin(i*np.pi/3),0]
        print(coord[inode,:])
        triang[itria,:] = np.sort([0,np.mod(i+1,6)+1,i+1])
        itria=itria+1

    return triang,coord;



# print str(fileout)[-3]+'.vtk'
#mt.write_grid(coord,triang,str(fileout)[:,-3]+'.vtk','vtk')



if __name__ == "__main__":
    if len(sys.argv) > 1 :
        topol,coord=hexagonal_grid(*sys.argv[2:])
        print(sys.argv[-1])
        mt.write_grid(coord,topol,sys.argv[1],'dat')
    else:
        raise SystemExit("usage:  python hexagonal_grid.py  <fileout> [xcenter ycenter edge]")
    


