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

def define_optimal_tdens(flag_grid, flag_source, flag_sink,
                         extra_grid, extra_source, extra_sink,):
    #####################################################################
    # optimal transport density definition
    #####################################################################
    print((flag_source=='rect_cnst') and 
            (flag_sink=='rect_cnst'))
    print (flag_source, flag_sink)
    if (    (flag_source=='rect_cnst') and 
            (flag_sink=='rect_cnst') ):
        def optimal_tdens(time,coord,flag):
            x=coord[0]
            y=coord[1]
            z=coord[2]
            value=0.0
            steady=True
            
            # supportrs info
            base=np.array([0.125,0.25])
            width=0.25
            height=0.5
            shift=0.25
            forcing_value=2.0
            if ( (x >  base[0]        ) and 
                 (x <= base[0]+width  ) and 
                 (y >  base[1]        ) and 
                 (y <= base[1]+height ) )  :
                value=(x-base[0])*forcing_value
            if ( (x >  base[0]+width        ) and 
                 (x <= base[0]+width+shift  ) and 
                 (y >  base[1]              ) and 
                 (y <= base[1]+height       ) )  :
                value=(width)*forcing_value
            if ( (x >  base[0]+  width+shift ) and 
                 (x <= base[0]+2*width+shift ) and 
                 (y >  base[1]               ) and 
                 (y <= base[1]+height        ) )  :
                value=(base[0]+2*width+shift-x)*forcing_value
            return value,steady;
    #
    # optiamal tdens given in 
    # Branching structures emerging from a continuous
    # optimal transport model
    # Enrico Facca , Franco Cardin , Mario Putti
    # \Tdens^*=|1/r int_{0}^{r} t F(t) dt|^{p-2}/{p-1}
    if ((flag_grid  =='plaplacian') and 
        (flag_source=='plaplacian') and 
        (flag_sink  =='plaplacian') ):
        print(extra_source)
        def optimal_tdens(time,coord,flag):
            x=coord[0]
            y=coord[1]
            z=coord[2]
            value=0.0
            steady=True
            plapl=float(extra_source)
            r=np.sqrt(x**2+y**2)
            if ( r<1.0/3.0):
                value=r/2.0 # =r**2/(2r)
            if ( r>1.0/3.0) and ( r<2.0/3.0) :
                value=((1.0/3.0)**2/2.0)/r
            if ( r>2.0/3.0) :
                value=((1.0/3.0)**2/2.0)/r-(1.0/5.0)*(r**2-(2.0/3.0)**2)/(2.0*r)
            if (plapl<1e+15):
                value=abs(value)**((plapl-2)/(plapl-1))
            return value,steady;

                
    return optimal_tdens;        

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
