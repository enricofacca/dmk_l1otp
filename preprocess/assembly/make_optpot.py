# -*- coding: utf-8 -*-
#!/usr/bin/env python
import numpy as np
import sys
import os

#
# geometry
#
#sys.path.insert(0,'./geometry2d')
import meshtools as mt
import example_grid as ex_grid
from pyvtk import *

#
# inputs
#

import example_forcing as ex_forcing
import common as common
import cell as cell


############################################
# reads flags and other values
ctrl = open(f_ctrl, "r")
ctrl_lines = ctrl.readlines()
ctrl.close()

# flags grid
flag_grid=str(common.read_column('flag_grid',0,ctrl_lines))
extra_grid=str(common.read_column('flag_grid',1,ctrl_lines))
ndiv=int(common.read_column('ndiv',0,ctrl_lines))
nref=int(common.read_column('nref',0,ctrl_lines))

# flag_source
flag_source=str(common.read_column('flag_source',0,ctrl_lines))
extra_source=str(common.read_column('flag_source',1,ctrl_lines))

#flag sink
flag_sink=str(common.read_column('flag_sink',0,ctrl_lines))
extra_sink=str(common.read_column('flag_sink',1,ctrl_lines))

#####################################################################
# OPTPOT optimal potential
#####################################################################
optpot_functions=[]
optpot_functions.append(ex_forcing.known_optpot(flag_grid,flag_source,flag_sink))
# fill data
if( len(optpot_functions) == 1):
    # reads coord topol and flags 
    coord, topol, flags = mt.read_grid(extra_grid)
    if ( coord.shape[1] == 2 ):
        zcoord = np.zeros([coord.shape[0],1])
        coord=np.append(coord, zcoord, axis=1)

    optpot_node = np.zeros([nnode,1])
    for inode in range(len(coord)):
        optpot_node[inode]=optpot_functions[0](coord[inode],flags[inode])
        
    # write to file
    cell.write2file(filename_optdens,optdens_tria)
else:
    print('Optimal Potenial not defined')
    

