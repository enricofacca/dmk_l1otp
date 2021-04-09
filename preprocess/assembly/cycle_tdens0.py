# -*- coding: utf-8 -*-
#!/usr/bin/env python
import numpy as np
import sys
import meshtools as mt
import cell as cell


def write_tdens0(path_grid,dir_out,function_list): 
    #
    # read grid
    #
    coord, topol, flags = mt.read_grid(path_grid)
    if ( coord.shape[1] == 2 ):
        zcoord = np.zeros([coord.shape[0],1])
        coord=np.append(coord, zcoord, axis=1)
    ncell=len(topol)

    size_cell=mt.make_size(coord,topol)
    length = np.sqrt( 2.0* np.sum(size_cell) / float( len(size_cell) ) )
    bar_cell=mt.make_bar(coord,topol)

    
    #
    # init array
    #
    tdens0_cell=np.zeros([ncell,1])

    iexample=0
    print len(function_list)
    for funct in function_list:
        iexample=iexample+1
        filename_tdens0=dir_out+'/tdens0_'+str(iexample)+'.dat'
        for icell in range(ncell):
            tdens0_cell[icell]= funct(
                bar_cell[icell][:],
                flags[icell])
        cell.write2file(filename_tdens0,tdens0_cell)

path_grid=sys.argv[1]
dir_out=sys.argv[2]

function_list=[]
function_list.append(cell.example('1.0'))
function_list.append(cell.example('par0505'))

write_tdens0(path_grid,dir_out,function_list)

