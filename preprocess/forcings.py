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
import example_grid as ex_grid
import example_forcing as ex_forcing
import timecell as timecell
import make_optimal_tdens as make_opttdens
import meshtools as mt
sys.path.append('../preprocess/')
import grids as grids



def assembly_forcing(bar_cell,flags,
                     flag_source,flag_sink,
                     extra_source,extra_sink):
    ncell=bar_cell.shape[0]
    print(ncell)
    
    # build
    sources=[];
    sinks=[];
    dirac_sources=[];
    dirac_sinks=[];
    source_tria=np.zeros([ncell,1])
    sink_tria=np.zeros([ncell,1])
    steady_source=True
    steady_sink=True
    ex_forcing.example_source(sources,dirac_sources,
                              str(flag_source),str(extra_source))
    ex_forcing.example_sink(sinks,dirac_sinks,
                            str(flag_sink),str(extra_sink))
    print(len(sources),len(sinks))
    source_tria, steady_source=ex_forcing.make_source(
        source_tria,steady_source,
        sources,
        0,flags,bar_cell)
    print(min(source_tria),max(source_tria))
    

    sink_tria, steady_sink=ex_forcing.make_source(
        sink_tria,steady_sink,
        sinks,
        0,flags,bar_cell)
    print(min(sink_tria),max(sink_tria))
    forcing=source_tria-sink_tria

    return forcing.reshape([ncell,]);

