#!/usr/bin/env python
import numpy as np
import sys
import meshtools as mt
import timecell as timecell


def forcing(flag_rhs,fgraph=None,frhs=None):
    if ( 'help' in flag_rhs ) :
        print ('Use: "make_graph_forcing list" to print options')
        print ('Use: "make_graph_forcing flag graph forcing"')
        print ('flag    = availble flag')
        print ('graph   = graph data')
        print ('forcing = timedata with forcing term')
        return
    
    print_description=False
    if ( 'list' in flag_rhs):
        print_description=True

    if ( not print_description):
        coord,topol,flags =mt.read_grid(fgraph)
        nnode=len(coord)
        ncell=len(topol)
        rhs=np.zeros([nnode,1])

    print(flag_rhs,'  :')
        
    if (flag_rhs=='rect_cnst'):
        
        if (print_description):
            print(flag_rhs,'  :')
            print('Piecewise constant forcing')
            print('Source==cnst  on 0.125<x<0.375')
            print('Source==-cnst on 0.125<x<0.375')
            print('cnst scale wiht number of nodes')
            return
        else:
            ndiv=int(np.sqrt(nnode)-1)
            print (ndiv)
            value=2.0/(ndiv/4+1)


            rhs=np.zeros([nnode,1])
            region=np.zeros([nnode,1])
            for inode in range(nnode):
                x=coord[inode,0]
                y=coord[inode,1]

                
                if ( (x>=1.0/8.0 and x<=3.0/8.0) and
                     (y>=1.0/4.0 and y<=3.0/4.0) ):
                    rhs[inode,0]=value
                if ( (x>=5.0/8.0 and x<=7.0/8.0) and
                     (y>=1.0/4.0 and y<=3.0/4.0) ):
                    rhs[inode,0]=-value
                if ( ( y > 3.0/4.0-5/ndiv) and y<0.9):
                    region[inode,0]=1
                
    if (flag_rhs=='eikonal0500'):
        if (print_description):
            print('Forcing FOr Eikonal solution for node x,y=(0.5,0)')
            print('Source==1  at the closer node to (0.5,0)')
            print('Sink==-1/(nnode-1) to balance')
            return
        else:
            rhs[:]=-1.0/(nnode-1)
            inode=mt.Inode(coord,[0.5,0.0,0.0])
            rhs[inode]=1.0

            print((coord[inode]))

    if (not print_description):
        # open file and write dimensions
        file_out=open(frhs, 'w')
        time=0.0
        timecell.write2file(file_out,time,True,rhs)
        steady=True
        timecell.write2file(file_out,time,False,rhs,steady)
        file_out.close()

        file_out=open('region.dat', 'w')
        time=0.0
        timecell.write2file(file_out,time,True,region)
        steady=True
        timecell.write2file(file_out,time,False,region,steady)
        file_out.close()


if __name__ == "__main__":
    if len(sys.argv) > 1 :
        forcing(*sys.argv[1:])
    else:
        raise SystemExit("usage:  python graph rhs_flag rhs_file")
