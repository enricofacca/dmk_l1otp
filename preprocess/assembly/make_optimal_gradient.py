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

def define_optimal_gradient(flag_grid, flag_source, flag_sink,
                         extra_grid, extra_source, extra_sink,):
    #####################################################################
    # optimal transport density definition
    #####################################################################
    if ( # ( (flag_gridd='rect_cnst') or  (flag_source='rect_cnst_aligned') ) and 
            (flag_source=='rect_cnst') and 
            (flag_sink=='rect_cnst') ):
        def optimal_gradient_x(time,coord,flag):
            x=coord[0]
            y=coord[1]
            z=coord[2]
            value=-1.0
            steady=True
            return value,steady;
        def optimal_gradient_y(time,coord,flag):
            value=0.0
            steady=True
            return value,steady;

        def optimal_gradient_z(time,coord,flag):
            value=0.0
            steady=True
            return value,steady;

        return [optimal_gradient_x, optimal_gradient_y, optimal_gradient_z];        

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

        opt_grads_functions=[]
        opt_grads_functions.append(define_optimal_gradient(
            flag_grid, flag_source, flag_sink, 
            extra_grid, extra_source, extra_sink))

        print opt_grads_functions
        
        if( len(opt_grads_functions) >= 1):
            # reads coord topol and flags 
            coord, topol, flags = mt.read_grid(fin_grid)
            if ( coord.shape[1] == 2 ):
                zcoord = np.zeros([coord.shape[0],1])
                coord=np.append(coord, zcoord, axis=1)
            ncell=len(topol)
            bar_cell=mt.make_bar(coord,topol)


            # definitions 
            optgrad_cell = np.zeros([ncell,3])
            steady=False
            time=0.0
            
            # open file and write dimensions
            file_out=open(fout_optdens, 'w')
            timecell.write2file(file_out,time,True,optgrad_cell)
            
            for j in range(3):
                optgrad_cell[:,j], steady = timecell.build(
                    optgrad_cell[:,j], steady, 
                    [opt_grads_functions[0][j]],time,bar_cell,flags)
                
            # write 2 file
            timecell.write2file(file_out,time,False,optgrad_cell,steady)
            
            file_out.close()
                
        else:
            print('Optimal transport density not defined')
        
    else:
        raise SystemExit("usage:  python make_optimal_gradient.py 'input ctrl file' 'grid file' 'output file'"  )
