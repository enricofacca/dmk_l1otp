#setup 
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import network_impainting_tools as NetImp
from ExampleDerivedTypes import ( Dmkcontrols,
                                  Abstractgeometry,
                                  Dmkinputsdata,
                                  example_monge_kantorovich_spectral,
                                  otpdmk,data2grids,dmkp1p0_steady_data,
                                  build_subgrid_rhs,
                                  Tdenspotentialsystem)

import sys
sys.path.append('../../globals/python_timedata')
import timedata as td
sys.path.append('../preprocess/2d_assembly/')
import timecell as timecell
import make_optimal_tdens as make_opttdens
import meshtools as mt
sys.path.append('../preprocess/')


def assembly_optimal_tdens(bar_cell,flags,
                     flag_grid,flag_source,flag_sink,
                     extra_source,extra_sink):
    ncell=bar_cell.shape[0]    
    opt_tdens_functions=[]
    opt_tdens_functions.append(
        make_opttdens.define_optimal_tdens(
            flag_grid, flag_source, flag_sink, 
            ' ', ' ', ' '))
    optdens_cell = np.zeros(ncell)
    steady=False
    time=0.0
    optdens_cell, steady = timecell.build(
        optdens_cell, steady, 
        opt_tdens_functions,time,bar_cell,flags)

    return optdens_cell;

