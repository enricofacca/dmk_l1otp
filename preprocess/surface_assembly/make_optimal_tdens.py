# -*- coding: utf-8 -*-
#!/usr/bin/env python
import numpy as np
import sys
import os

#
# geometry
#
import meshtools_surface as mt
import example_grid as ex_grid

#
# inputs
#
import example_forcing as ex_forcing
import common as common
import timecell as timecell

def define_optimal_tdens(flag_grid, flag_source, flag_sink,
                         extra_grid, extra_source, extra_sink):
    #####################################################################
    # optimal transport density definition
    #####################################################################
    #
    # exact solution for dmk_surface paper
    #
    print flag_source, flag_sink
    if ( flag_source == 'north_band' and
         flag_sink   == 'south_band' ) :
        def  optdens(time,coord,flag_cell):
            x=coord[0]
            y=coord[1]
            z=coord[2]
            
            #phi=np.atan2(y,x)
            value=0.0
            if ( (x>=0.0) and (x<=1.0 ) and 
                 (y>=0.0) and (y<=1.0 ) ):
                r=np.sqrt(x**2+y**2+z**2)
                theta=np.arccos(z)
                    
                if ( z > np.sqrt(3)/2.0):
                    value=0.0
                    
                if ( (z <= np.sqrt(3)/2.0) and (z>=0.5) ):
                    value=(np.sqrt(3)/2.0 - z ) / np.sin(theta)
                
                if ( (z < 0.5) and (z>-0.5) ):
                    value=(np.sqrt(3)/2.0 - 0.5) / np.sin(theta)

                if ( (z <= -0.5) and (z>-np.sqrt(3)/2.0) ):
                    value=(np.sqrt(3)/2.0 + z) / np.sin(theta)
                
                if ( z<=-np.sqrt(3)/2.0 ):
                    value=0.0
            if (value<0.0):
                value=0.0
            steady=True
            return value,steady;
        
    return optdens;

if __name__ == "__main__":
    if len(sys.argv) > 1:
        fin_ctrl=sys.argv[1]  
        fin_grid=sys.argv[2]
        fout_optdens=sys.argv[3]
        
        ############################################
        # reads flags and other values
        ############################################
        ctrl = open(fin_ctrl, "r")
        ctrl_lines = ctrl.readlines()
        ctrl.close()
        
        # flags grid
        flag_grid=str(common.read_column('flag_grid',0,ctrl_lines))
        extra_grid=str(common.read_column('flag_grid',1,ctrl_lines))
        
        # flag_source
        flag_source=str(common.read_column('flag_source',0,ctrl_lines))
        extra_source=str(common.read_column('flag_source',1,ctrl_lines))
        
        #flag sink
        flag_sink=str(common.read_column('flag_sink',0,ctrl_lines))
        extra_sink=str(common.read_column('flag_sink',1,ctrl_lines))

        opt_tdens_functions=[]
        opt_tdens_functions.append(define_optimal_tdens(
            flag_grid, flag_source, flag_sink, 
            extra_grid, extra_source, extra_sink))
        
        if( len(opt_tdens_functions) >= 1):
            # reads coord topol and flags 
            coord, topol, flags = mt.read_grid(fin_grid)
            if ( coord.shape[1] == 2 ):
                zcoord = np.zeros([coord.shape[0],1])
                coord=np.append(coord, zcoord, axis=1)
            ncell=len(topol)
            bar_cell=mt.make_bar(coord,topol)


            # definitions 
            optdens_cell = np.zeros([ncell,1])
            steady=False
            time=0.0
            
            # open file and write dimensions
            file_out=open(fout_optdens, 'w')
            timecell.write2file(file_out,time,True,optdens_cell)
            
            optdens_cell, steady = timecell.build(
                optdens_cell, steady, 
                opt_tdens_functions,time,bar_cell,flags)
            # write 2 file
            timecell.write2file(file_out,time,False,optdens_cell,steady)
            
            file_out.close()
                
        else:
            print('Optimal transport density not defined')
        
    else:
        raise SystemExit("usage:  python make_optimal_tdens.py 'input ctrl file' 'grid file' 'output file'"  )
