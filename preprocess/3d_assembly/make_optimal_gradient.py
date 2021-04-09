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
    # optimal gradient definition
    #####################################################################
    if (  (flag_grid=='sphere_plaplacian') and 
          (flag_source=='sphere_plaplacian') and 
          (flag_sink=='sphere_plaplacian') ):
        
        # radius_l := large radius
        # radius_m := medium radius
        # radius_s := small radius
        radius_l = 1.0
        radius_m = 2.0/3.0
        radius_s = 1.0/3.0 
        label_sink = -1
        label_middle = 0 
        label_source = 1
        center_domain = [0,0,0]
        
        # the process of building the mesh modifies the coordinates of the vetices
        # of the three isospheres, so we need to evaluate the new three radius
        # new_radius_s = 0.316418075535
        # new_radius_m = 0.64956704691
        # new_radius_l = 0.974350570365
        # update radius
        # radius_s = new_radius_s
        # radius_m = new_radius_m
        # radius_l = new_radius_l
        
        # cnst_value=1.0
        volume_source = (4.0/3.0)*np.pi*radius_s**3
        volume_sink   = (4.0/3.0)*np.pi*radius_l**3 - (4.0/3.0)*np.pi*radius_m**3
        cnst_value = 1.0 / volume_source
        cnst_sink_value = 1.0 / volume_sink
        # cnst_sink_value = (volume_source*cnst_value)/volume_sink
        # cnst_value = (volume_sink*cnst_sink_value)/volume_source
        c1 = cnst_value
        c2 = - cnst_sink_value


        #
        # read pflux  from input.ctrl 
        # absolute path must be written in extra_source
        #
        
        ctrl = open(extra_source, "r")
        ctrl_lines = ctrl.readlines()
        pflux=float(common.read_column('flag_pflux',0,ctrl_lines))
        ctrl.close()

        plapl=(2.0-pflux)/(1.0-pflux)

        # define optimal gradient
        def optimal_gradient_x(time,coord,flag_cell):
            
            x=coord[0]
            y=coord[1]
            z=coord[2]

            zeta_fun_value=0.0
            r_dist = np.sqrt(x**2+y**2+z**2)
            if (flag_cell > 0):
                zeta_fun_value = - c1*r_dist/3.0
                zeta_fun_value = np.abs(zeta_fun_value)
                # if (r_dist > 0.30):
                #     print r_dist,zeta_fun_value
            elif (flag_cell == 0):
                zeta_fun_value =  - c1*(radius_s)**3/(r_dist**2*3.0)
                zeta_fun_value = np.abs(zeta_fun_value)
                # if (r_dist > 0.64):
                #     print r_dist,zeta_fun_value
            else:
                zeta_fun_value = -1.0/r_dist**2*(
                    c1*radius_s**3/3.0+c2*(r_dist**3-radius_m**3)/3.0)
                zeta_fun_value = np.abs(zeta_fun_value)
            
            value = (zeta_fun_value/r_dist**2)**(1.0/(plapl-1.0)) * x /r_dist
            
            steady=True
            return value,steady;

        # define optimal gradient
        def optimal_gradient_y(time,coord,flag_cell):
            
            x=coord[0]
            y=coord[1]
            z=coord[2]

            zeta_fun_value=0.0
            r_dist = np.sqrt(x**2+y**2+z**2)
            if (flag_cell > 0):
                zeta_fun_value = - c1*r_dist/3.0
                zeta_fun_value = np.abs(zeta_fun_value)
                # if (r_dist > 0.30):
                #     print r_dist,zeta_fun_value
            elif (flag_cell == 0):
                zeta_fun_value =  - c1*(radius_s)**3/(r_dist**2*3.0)
                zeta_fun_value = np.abs(zeta_fun_value)
                # if (r_dist > 0.64):
                #     print r_dist,zeta_fun_value
            else:
                zeta_fun_value = -1.0/r_dist**2*(
                    c1*radius_s**3/3.0+c2*(r_dist**3-radius_m**3)/3.0)
                zeta_fun_value = np.abs(zeta_fun_value)
            
            value = (zeta_fun_value/r_dist**2)**(1.0/(plapl-1.0)) * y /r_dist
            
            steady=True
            return value,steady;

        # define optimal gradient
        def optimal_gradient_z(time,coord,flag_cell):
            
            x=coord[0]
            y=coord[1]
            z=coord[2]

            zeta_fun_value=0.0
            r_dist = np.sqrt(x**2+y**2+z**2)
            if (flag_cell > 0):
                zeta_fun_value = - c1*r_dist/3.0
                zeta_fun_value = np.abs(zeta_fun_value)
                # if (r_dist > 0.30):
                #     print r_dist,zeta_fun_value
            elif (flag_cell == 0):
                zeta_fun_value =  - c1*(radius_s)**3/(r_dist**2*3.0)
                zeta_fun_value = np.abs(zeta_fun_value)
                # if (r_dist > 0.64):
                #     print r_dist,zeta_fun_value
            else:
                zeta_fun_value = -1.0/r_dist**2*(
                    c1*radius_s**3/3.0+c2*(r_dist**3-radius_m**3)/3.0)
                zeta_fun_value = np.abs(zeta_fun_value)
            
            value = (zeta_fun_value/r_dist**2)**(1.0/(plapl-1.0)) * z /r_dist
            
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
