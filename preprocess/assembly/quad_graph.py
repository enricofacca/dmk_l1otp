#!/usr/bin/env python
import numpy as np
import sys
import meshtools as mt

def rectangle_grid(ndivx,ndivy,fileout):    
    # grid data
    print ndivx, ndivy
    ndivx=int(ndivx)
    ndivy=int(ndivy)
    
    ntria=ndivx*ndivy
    nnode=(ndivx+1)*(ndivy+1)
    coord=np.zeros([nnode,2])
    triang=np.zeros([ntria,2],dtype=int)
        
    # aux. data
    ntriax=2*ndivx
    nnodex=ndivx+1
    lenx = 1.0 / ndivx
    leny = lenx

    inode=0
    itria=0
    for iy in range(ndivy): # 0,1, ndivy-1
        for ix in range(ndivx): # 0,1, ndivx-1
            coord[inode,:] = (ix*lenx,iy*leny)
            inode = inode+1
            # 
            # |   
            # |   
            # ------ 
            
            # lower edge 
            n1= iy * nnodex + ix 
            n2= n1 + 1
            triang[itria,:]=(n1,n2)
            
            #print 'itria=',itria, 'n=', n1,n2,n3
            
            # vertical edge
            itria=itria+1
            # copy
            n1old=n1
            n2old=n2
            n3old=n3

            # ne triangle 
            n1= n1
            n2= n1 + ndivx
            
            #print 'itria=',itria, 'n=', n1,n2,n3

            triang[itria,:]=(n1,n1+ndivx)

            itria=itria+1

        # add last point in x direction
        coord[inode,:] = ((ix+1)*lenx,iy*leny)
        inode = inode+1
        triang[itria,:]=(n1,n1+ndivx)
                    
   

    # add line of top points 
    for ix in range(ndivx+1): # 0,1, ndivx
        coord[inode,:] = (ix*lenx,(iy+1)*leny)
        inode = inode+1

    mt.write_grid(coord,triang,fileout,'dat')
    # print str(fileout)[-3]+'.vtk'
    #mt.write_grid(coord,triang,str(fileout)[:,-3]+'.vtk','vtk')



if __name__ == "__main__":
    if len(sys.argv) > 1 :
        rectangle_grid(*sys.argv[1:])
    else:
        raise SystemExit("usage:  python square_grid.py ndivx ndivy <fileout>")
    


