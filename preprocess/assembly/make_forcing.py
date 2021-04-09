# -*- coding: utf-8 -*-
#!/usr/bin/env python
import numpy as np
import sys
import os

#
# geometry
#
import meshtools as mt
import example_grid as ex_grid

#
# inputs
#
import example_forcing as ex_forcing
import common as common
import timecell as timecell

def define_forcing(flag_forcing,extra_info=None):
    #####################################################################
    # optimal transport density definition
    #####################################################################
    if ( flag_forcing == 'cnst_center'):
        def forcing(time,coord,flag):
            x=coord[0]
            y=coord[1]
            z=coord[2]
            value=0.0
            steady=True
            
            fvalue=0.0
            steady=True
            if ( (x >= 0.4) and (x<=0.6) and 
                 (y >= 0.4) and (y<=0.6) ) :
                fvalue=1.0
            else:
                fvalue=-1.0
            return fvalue,steady;
            
        return forcing; 
        
    if ( flag_forcing == 'eikonal'):
        def forcing(time,coord,flag):
            x=coord[0]
            y=coord[1]
            z=coord[2]
            value=0.0
            steady=True
            
            fvalue=0.0
            steady=True
            if ( (x >= 0.5-1e-10) and (x<=0.5+1e-10) and 
                 (y >= 0.5-1e-10) and (y<=0.5+1e-10) ) :
                fvalue=1.0
            else:
                fvalue=-1.0
            return fvalue,steady;
            
        return forcing;        
       

if __name__ == "__main__":
    if len(sys.argv) > 1:
        flag_forcing=sys.argv[1] 
        flag_node_baricenter=sys.argv[2] 
        fin_grid=sys.argv[3]
        fout_forcing=sys.argv[4]
        
        ############################################
        # reads flags and other values
        ############################################
       
        
       

        forcing_functions=[]
        forcing_functions.append(define_forcing(flag_forcing))
        
        if( len(forcing_functions) >= 1):
            # reads coord topol and flags 
            coord, topol, flags = mt.read_grid(fin_grid)
            if ( coord.shape[1] == 2 ):
                zcoord = np.zeros([coord.shape[0],1])
                coord=np.append(coord, zcoord, axis=1)
            nnode=len(coord)
            ncell=len(topol)
            bar_cell=mt.make_bar(coord,topol)


            # definitions 
            if ( flag_node_baricenter == 'node'):
                forcing_evalute = np.zeros([nnode,1])
            else:
                forcing_evalute = np.zeros([ncell,1])
            
            
                
            steady=False
            time=0.0
            
            print(coord.shape)
            print(bar_cell.shape)

            # open file and write dimensions
            file_out=open(fout_forcing, 'w')
            timecell.write2file(file_out,time,True,forcing_evalute)
            if ( flag_node_baricenter == 'node'):
                forcing_evalute, steady = timecell.build(
                    forcing_evalute, steady, 
                    forcing_functions,time,coord,flags)
            else:
                forcing_evalute, steady = timecell.build(
                    forcing_evalute, steady, 
                    forcing_functions,time,bar_cell,flags)

            # write 2 file
            timecell.write2file(file_out,time,False,forcing_evalute,steady)
            
            file_out.close()
                
        else:
            print('Forcing not defined')
        
    else:
        raise SystemExit("usage:  python make_optimal_tdens.py 'input ctrl file' 'grid file' 'output file'"  )
