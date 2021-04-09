#setup 
import numpy as np
import sys
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import network_impainting_tools as NetImp
sys.path.append('../new_python_interface/')
from ExampleDerivedTypes import ( Abstractgeometry,
                                  data2grids,dmkp1p0_steady_data)

sys.path.append('../globals/python_timedata')
import timedata as td
sys.path.append('2d_assembly/')
import example_grid as ex_grid
import example_forcing as ex_forcing
import timecell as timecell
import make_optimal_tdens as make_opttdens
import meshtools as mt

def assembly_grid(flag_grid,length,extra_grid):
    (points,
     vertices,
     coordinates,
     elements,
     element_attributes )= ex_grid.example_grid(flag_grid,length,extra_grid)
    coord=np.array(coordinates)
    if ( coord.shape[1] == 2 ):
        zcoord = np.zeros([coord.shape[0],1])
        coord=np.append(coord, zcoord, axis=1)
    topol=np.array(elements)
    try:
        print (len(element_attributes))
        flags=np.array(element_attributes)
    except NameError:
        flags=np.int_(range(len(topol)))

    return topol,coord,flags;

def py2fortran(topol,coord):
    topolT=topol.transpose()
    coordT=coord.transpose()
    topolT=topolT+1

    grid=Abstractgeometry.abs_simplex_mesh()
    subgrid=Abstractgeometry.abs_simplex_mesh()
    data2grids(6, 3, len(coord), len(topol), coordT, topolT, grid,subgrid)
    return grid,subgrid;
