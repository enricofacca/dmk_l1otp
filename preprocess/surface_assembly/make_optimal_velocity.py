# -*- coding: utf-8 -*-
#!/usr/bin/env python
import numpy as np
import sys
import os
import sys
# insert at 1, 0 is the script path (or '' in REPL)
sys.path.insert(1, '../../../globals/python_timedata/')
import timedata as td



#
# geometry
#
import meshtools_surface as mt

#
# inputs
#
import common as common
import timecell as timecell

def define_optimal_vel(flag_source, flag_sink,
                         extra_source, extra_sink):
    #####################################################################
    # optimal transport density definition
    #####################################################################
    #
    # exact solution for dmk_surface paper
    #
    print flag_source, flag_sink
    if ( flag_source == 'north_band' and
         flag_sink   == 'south_band' ) :
        def  opt_vel(coord):
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
            #v1=np.array([x,y,z])
            #v2=np.array([x,-y,0])
            gx=-x*z/np.sqrt(1-z**2)
            gy=-y*z/np.sqrt(1-z**2)
            gz=np.sqrt(1-z**2)
            gradient=np.array([gx, gy, gz])
            gradient=gradient/np.linalg.norm(gradient)
            value_x=-value*gradient[0]
            value_y=-value*gradient[1]
            value_z=-value*gradient[2]
            return value_x,value_y,value_z;
        
    return opt_vel;

if __name__ == "__main__":
    if len(sys.argv) > 1:
        fin_ctrl=sys.argv[1]  
        fin_grid=sys.argv[2]
        fout_opt_vel=sys.argv[3]
        
        ############################################
        # reads flags and other values
        ############################################
        ctrl = open(fin_ctrl, "r")
        ctrl_lines = ctrl.readlines()
        ctrl.close()
        
        # flags grid
        #flag_grid=str(common.read_column('flag_grid',0,ctrl_lines))
        #extra_grid=str(common.read_column('flag_grid',1,ctrl_lines))
        
        # flag_source
        flag_source=str(common.read_column('flag_source',0,ctrl_lines))
        extra_source=str(common.read_column('flag_source',1,ctrl_lines))

        
        #flag sink
        flag_sink=str(common.read_column('flag_sink',0,ctrl_lines))
        extra_sink=str(common.read_column('flag_sink',1,ctrl_lines))
        


        opt_vel_functions=[]
        opt_vel_functions.append(define_optimal_vel(
            flag_source, flag_sink, 
            extra_source, extra_sink))
        
        print(len(opt_vel_functions),flag_source, flag_sink)
        
        if( len(opt_vel_functions) >= 1):
            # reads coord topol and flags 
            coord, topol, flags = mt.read_grid(fin_grid)
            if ( coord.shape[1] == 2 ):
                zcoord = np.zeros([coord.shape[0],1])
                coord=np.append(coord, zcoord, axis=1)
            ncell=len(topol)
            bar_cell=mt.make_bar(coord,topol)


            # eval definitions 
            opt_vel_cell = np.zeros([ncell,3])
            for icell in range(len(bar_cell)):
                #print (opt_vel_cell[icell,:])
                #print (bar_cell[icell,:])
                opt_vel_cell[icell,:] = opt_vel_functions[0](bar_cell[icell,:])
            
            # write 2 file
            td.write_steady_timedata(fout_opt_vel,opt_vel_cell)
                            
        else:
            print('Optimal transport density not defined')
        
    else:
        raise SystemExit("usage:  python make_optimal_vel.py 'input ctrl file' 'grid file' 'output file'"  )
