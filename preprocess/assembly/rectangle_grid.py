#!/usr/bin/env python
import numpy as np
import sys
import meshtools as mt

def rectangle_grid(ndivx,ndivy,L=1.0):    
    # grid data
    print (ndivx, ndivy,L)
    ndivx=int(ndivx)
    ndivy=int(ndivy)
    
    ntria=2*ndivx*ndivy
    nnode=(ndivx+1)*(ndivy+1)
    coord=np.zeros([nnode,3])
    triang=np.zeros([ntria,3],dtype=int)
        
    # aux. data
    ntriax=2*ndivx
    nnodex=ndivx+1
    lenx = L / ndivx
    leny = lenx

    inode=0
    itria=0
    for iy in range(ndivy): # 0,1, ndivy-1
        for ix in range(ndivx): # 0,1, ndivx-1
            coord[inode,:] = (ix*lenx,iy*leny,0.0)
            inode = inode+1
            # ------
            # |    |
            # |\   |
            # | \  |
            # |  \ |
            # |   \|
            # ------ 
            # sw triangle 
            n1= iy * nnodex + ix 
            n2= n1 + 1
            n3= n1 + nnodex 
            triang[itria,:]=(n1,n2,n3)
            
            #print 'itria=',itria, 'n=', n1,n2,n3
            

            itria=itria+1
            # copy
            n1old=n1
            n2old=n2
            n3old=n3

            # ne triangle 
            n1= n2old
            n2= n3old + 1
            n3= n3old

            #print 'itria=',itria, 'n=', n1,n2,n3

            triang[itria,:]=(n1,n2,n3)

            itria=itria+1

        # add last point in x direction
        coord[inode,:] = ((ix+1)*lenx,iy*leny,0.0)
        inode = inode+1
   

    # add line of top points 
    for ix in range(ndivx+1): # 0,1, ndivx
        coord[inode,:] = (ix*lenx,(iy+1)*leny,0.0)
        inode = inode+1

    return triang,coord;



# print str(fileout)[-3]+'.vtk'
#mt.write_grid(coord,triang,str(fileout)[:,-3]+'.vtk','vtk')



if __name__ == "__main__":
    if len(sys.argv) > 1 :
        topol,coord=rectangle_grid(*sys.argv[1:3])
        print(sys.argv[-1])
        mt.write_grid(coord,topol,sys.argv[-1],'dat')
    else:
        raise SystemExit("usage:  python square_grid.py ndivx ndivy <fileout>")
    


