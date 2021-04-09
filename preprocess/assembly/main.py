# -*- coding: utf-8 -*-
#!/usr/bin/env python

#######################################################
# This program takes 
# f_ctrl        : inputs file wiht flag/controls/ path 
#                 See 'inputs.ctrl' in ../script_controls/
# dir_inputs    : folder where store all inputs data
# dir_inputs    : folder where store data in vtk format
#                
#######################################################

import numpy as np
import subprocess
import sys
import os
import argparse
from shutil import copyfile


path_td= os.path.abspath(
    os.path.normpath(
        os.path.dirname(os.path.realpath(__file__))+
        '../../../../../globals/python/timedata'))
print(path_td)
sys.path.append(path_td)
import timedata as td

#
# geometry
#
#sys.path.insert(0,'./geometry2d')
import meshtools as mt
import example_grid as ex_grid
import example_grid3d as ex_grid3d



#
# inputs
#

import example_forcing as ex_forcing
import make_dirichlet as dirichlet
import make_neumann as neumann
import make_optimal_tdens as make_optdens

import common as common
import timecell as timecell
import cell as cell
import scalars as scalars

# DESCRIPTION
# write in "pathfile" the time-varying data
# defined by "functions" on the "evaluation_points".
# INPUTS/OUTPUT
# ndata     : number of evalution points 
# dimdata   : codomain dimension of functions
# tzero     : initial time, time-stepping is given by funcition common.next_time(time)
# functions : list of functions in the form
#     def fun(t,x,y,z,flag):
#         return value(real value) , steady(logical)
# evulation_points : numpy array with xyz coordinates of ndata points
# flags_points     : numpy array with integer flag for each point(ignored most of the time)
# path_file        : srting with output file path
#
def writetimedata(ndata,dimdata,tzero,functions,evaluation_points,flags_points,pathfile):
    # definitions 
    data=np.zeros([ncell,dimdata])
    steady=False
    time=tzero

    # open file and write dimensions
    file_out=open(pathfile, 'w')
    timecell.write2file(file_out,time,True,data)

    #cycle in time
    while (time <= tmax) and (not steady ):
        # evalute 
        data, steady = timecell.build(
        data, steady, 
        functions,time,evaluation_points,flags_points)
    # write 2 file
    timecell.write2file(file_out,time,False,data,steady)
    
    # next time
    time=common.next_time(time)
    file_out.close()

    return



def assembly_inputs(f_ctrl,dir_inputs,dir_inputs_vtk):
    if ( dir_inputs == ''):
        dir_inputs='./'
    if ( dir_inputs_vtk == ''):
        dir_inputs_vtk='./'

    # set grid destination
    file_grid=dir_inputs+'/'+'grid.dat'
    file_grid_vtk=dir_inputs_vtk+'/'+'grid'
    file_poly_vtk=dir_inputs_vtk+'/'+'poly'

    # set inputs destinations
    filename_forcing=dir_inputs+'/'+'forcing.dat'
    filename_source=dir_inputs+'/'+'source.dat'
    filename_sink=dir_inputs+'/'+'sink.dat'
    filename_dirac_forcing=dir_inputs+'/'+'dirac_forcing.dat'
    filename_dirac_source=dir_inputs+'/'+'dirac_source.dat'
    filename_dirac_sink=dir_inputs+'/'+'dirac_sink.dat'
    filename_dirichlet=dir_inputs+'/'+'dirichlet.dat'
    filename_neumann=dir_inputs+'/'+'neumann.dat'
    

    # cell time data (changing in time)
    filename_kappa=dir_inputs+'/'+'kappa.dat'
    filename_lambda=dir_inputs+'/'+'lambda.dat'


    # scalar time data (changing in time)
    filename_pflux=dir_inputs+'/'+'pflux.dat'
    filename_pode=dir_inputs+'/'+'pode.dat'
    filename_pmass=dir_inputs+'/'+'pmass.dat'
    filename_decay=dir_inputs+'/'+'decay.dat'


    # cell data (fix in time)
    filename_measure_source=dir_inputs+'/'+'measure_source.dat'
    filename_measure_sink=dir_inputs+'/'+'measure_sink.dat'
    filename_tdens0=dir_inputs+'/'+'tdens0.dat'
    filename_optdens=dir_inputs+'/optdens.dat'

    #vtk
    filename_inputs_vtk=dir_inputs_vtk+'/'+'inputs'


    ############################################
    # reads flags and other values
    ctrl = open(f_ctrl, "r")
    ctrl_lines = ctrl.readlines()

    # flags equation
    flag_equation=str(common.read_column('flag_equation',0,ctrl_lines))
    if (flag_equation==''):
        flag_equation='poisson'
    elif ( not (flag_equation !='poisson' or flag_equation !='saintvenant') ):
        print ('Only equation poisson and saintvenant' )
        sys.exit()

   

    flag_domain=str(common.read_column('flag_domain',0,ctrl_lines))
    if ( flag_domain =='2d'):
        nnodeincell = 3
        ambient_dim = 2
        cell_dim = 2
        print('Assembler for 2d')
    elif ( flag_domain =='3d'):
        nnodeincell = 4
        ambient_dim = 3
        cell_dim = 3
        print('Assembler for 3d')
    elif ( flag_domain =='3d'):
        nnodeincell = 3
        ambient_dim = 3
        cell_dim = 2
        print('Assembler for Surface')
    else:
        print('Flag domain ',flag_domain, ' not supported')
        print('Use 2d , 3d or surface')

    if (flag_equation == 'poisson'):
        dimdata = 1
    if (flag_equation == 'saintvenant'):
        dimdata = ambient_dim
    print('Assembler for ', flag_equation)
    
    
    
    # flags grid
    flag_grid=str(common.read_column('flag_grid',0,ctrl_lines))
    extra_grid=str(common.read_column('flag_grid',1,ctrl_lines))
    ndiv=int(common.read_column('ndiv',0,ctrl_lines))
    nref=int(common.read_column('nref',0,ctrl_lines))

    # flag_inputs
    flag_source=str(common.read_column('flag_source',0,ctrl_lines))
    extra_source=str(common.read_column('flag_source',1,ctrl_lines))

    flag_sink=str(common.read_column('flag_sink',0,ctrl_lines))
    extra_sink=str(common.read_column('flag_sink',1,ctrl_lines))


    flag_dirichlet=str(common.read_column('flag_dirichlet',0,ctrl_lines))
    extra_dirichlet=str(common.read_column('flag_dirichlet',1,ctrl_lines))

    flag_neumann=str(common.read_column('flag_neumann',0,ctrl_lines))
    extra_neumann=str(common.read_column('flag_neumann',1,ctrl_lines))

    # grid dependent quantities
    flag_tdens0=str(common.read_column('flag_tdens0',0,ctrl_lines))
    extra_tdens0=str(common.read_column('flag_tdens0',1,ctrl_lines))

    flag_kappa=str(common.read_column('flag_kappa',0,ctrl_lines))
    extra_kappa=str(common.read_column('flag_kappa',1,ctrl_lines))

    flag_lambda=str(common.read_column('flag_lambda',0,ctrl_lines))
    extra_lambda=str(common.read_column('flag_lambda',1,ctrl_lines))

    # scalar quantities
    flag_pflux=str(common.read_column('flag_pflux',0,ctrl_lines))
    extra_pflux=str(common.read_column('flag_pflux',1,ctrl_lines))

    flag_pode=str(common.read_column('flag_pode',0,ctrl_lines))
    extra_pode=str(common.read_column('flag_pode',1,ctrl_lines))
    
    flag_pmass=str(common.read_column('flag_pmass',0,ctrl_lines))
    extra_pmass=str(common.read_column('flag_pmass',1,ctrl_lines))
    
    flag_decay=str(common.read_column('flag_decay',0,ctrl_lines))
    extra_decay=str(common.read_column('flag_decay',1,ctrl_lines))
    
    
    # time interval
    tzero=float(common.read_column('tzero',0,ctrl_lines))
    tmax=float(common.read_column('tmax',0,ctrl_lines))
    
    
    
    sources=[];
    sinks=[];
    dirac_sources=[];
    dirac_sinks=[];
    kappas=[];
    tdens0_functions=[];
    optdens_functions=[];
    dirac_points=[];
    
    
    length=0.0
    
    ###############################
    # build mesh
    if ( flag_grid != 'use'):
        length=1.0/float(ndiv)
        print(ambient_dim)
        if (ambient_dim == 2) :
            points, vertices, coordinates,elements,element_attributes = ex_grid.example_grid(flag_grid,length,extra_grid)
        if (ambient_dim == 3):
            points, vertices, coordinates,elements,element_attributes = ex_grid3d.example_grid(flag_grid,length,extra_grid)
        
        coord=np.array(coordinates)
        if ( coord.shape[1] == 2 ):
            zcoord = np.zeros([coord.shape[0],1])
            coord=np.append(coord, zcoord, axis=1)
        
        topol=np.array(elements)
        #print (element_attributes)
        try:
            flags=np.array(element_attributes)
        except NameError:
            flags=np.int_(range(len(topol)))
        
    
        
       # mt.write_poly(points,vertices,str(file_poly_vtk),'vtk')
    else:
        # reads new coord topol and flags of refined grid
        coord, topol, flags = mt.read_grid(extra_grid)
        if ( coord.shape[1] == 2 ):
            zcoord = np.zeros([coord.shape[0],1])
            coord=np.append(coord, zcoord, axis=1)
    
        
    
    
    #
    # call external program to refine mesh
    #
    if (nref > 0):
        file_grid0=dir_inputs+'/'+'grid0.dat'
        #file_grid0_vtk=dir_inputs_vtk+'/'+'grid0'
        mt.write_grid(coord,topol,file_grid0,'dat',flags)
        #mt.write_grid(coord,topol,file_grid0_vtk,'vtk',flags)
        path_nref = os.path.abspath(os.path.normpath(
                    os.path.dirname(os.path.realpath(__file__))+
                    '/uniform_refinement/code/uniform_refinement.out'))
        cmd = path_nref + \
            ' '+str(file_grid0)+' '+str(file_grid)+' '+str(nref)
        subprocess.call(cmd,shell=True)
        # reads new coord topol and flags of refined grid
        coord, topol, flags = mt.read_grid(file_grid)
        if ( coord.shape[1] == 2 ):
            zcoord = np.zeros((coord.shape[0],1))
            coord = np.append(coord, zcoord, axis=1)
    
    
    mt.write_grid(coord,topol,file_grid,'dat',flags)
    #mt.write_grid(coord,topol,str(file_grid_vtk),'vtk',flags)

    #
    # computing basic geomtry info of the mesh
    #
    ncell=len(topol)
    nnode=len(coord)
    bar_cell=mt.make_bar(coord,topol)
    size_cell=mt.make_size(coord,topol)
    
    if (cell_dim ==2):
        all_edges=mt.FindEdges(points,topol)
        nedge = len(all_edges)
    if (cell_dim == 3 ):
        # this is wrong and should count the boundary faces
        nedge = 4*ncell 
        
    ###############################
    # prepare sources ans sink functions if not defined
    # append sources or dirac nodes   
    if ( ( len(sources) == 0) and (len(dirac_sources) == 0) ) :
        ex_forcing.example_source(sources,dirac_sources,
                                  str(flag_source),str(extra_source),topol,coord)
    
    # append sinks or driac nodes
    if ( ( len(sinks) == 0) and (len(dirac_sinks) == 0) ) :
        ex_forcing.example_sink(sinks,dirac_sinks,
                            str(flag_sink),str(extra_sink),topol,coord)
    
    
    ##########################
    # forcing init.
    source_tria=np.zeros([ncell,dimdata])
    sink_tria=np.zeros([ncell,dimdata])
    forcing_tria=np.zeros([ncell,dimdata])
    dirac_forcing_nodes=np.zeros([nnode,dimdata])
    dirac_source_nodes=np.zeros([nnode,dimdata])
    dirac_sink_nodes=np.zeros([nnode,dimdata])

    #out_source_tria=np.zeros([ncell*dimdata,1])
    #out_sink_tria=np.zeros([ncell*dimdata,1])
    #out_forcing_tria=np.zeros([ncell*dimdata,1])
    #out_dirac_forcing_nodes=np.zeros([nnode*dimdata,1])
    #out_dirac_source_nodes=np.zeros([nnode*dimdata,1])
    #out_dirac_sink_nodes=np.zeros([nnode*dimdata,1])
    
    
    steady_source=False
    steady_sink=False
    steady_forcing=False
    dirac_source_steady=False
    dirac_sink_steady=False
    dirac_forcing_steady=False

    
    if ( True ) :
        time=tzero
       
    
        # open file and write dimensions
        file_forcing=open(filename_forcing, 'w')
        file_source=open(filename_source, 'w')
        file_sink=open(filename_sink, 'w')
    
        timecell.write2file(file_forcing,time,True,
                            #out_forcing_tria)
                            forcing_tria)
        timecell.write2file(file_source,time,True,
                            #out_source_tria)
                            source_tria)
        timecell.write2file(file_sink,time,True,
                            #out_forcing_tria)
                            sink_tria)
    
        # open file and write dimensions
        file_dirac_forcing=open(filename_dirac_forcing, 'w')
        file_dirac_source=open(filename_dirac_source, 'w')
        file_dirac_sink=open(filename_dirac_sink, 'w')
    
        timecell.write2file(file_dirac_forcing,time,True,
                            #out_dirac_forcing_nodes)
                            dirac_forcing_nodes)
        timecell.write2file(file_dirac_source,time,True,
                            #out_dirac_source_nodes)
                            dirac_source_nodes)
        timecell.write2file(file_dirac_sink,time,True,
                            #out_dirac_sink_nodes)
                            dirac_sink_nodes)
            
        #cycle in time
        while ( (time <= tmax) and 
                ( (not steady_forcing ) or (not steady_dirac_forcing) )
            ):
            #make source
            source_tria, steady_source = timecell.build(
                source_tria,steady_source,
                sources,
                time,bar_cell,flags)

            
            
            dirac_source_nodes, dirac_source_steady= ex_forcing.make_dirac_source(
                dirac_source_nodes,dirac_source_steady,dirac_sources,time,
                str(extra_source),coord)
            
            if ( str(flag_source) == 'use') :
                cmd='cp '+ str(extra_source) +' '+str(filename_source)
                subprocess.call(cmd,shell=True)
                steady_source=True
    
    
            #make sink
            sink_tria, steady_sink = timecell.build(
                sink_tria,steady_sink,
                sinks,
                time,bar_cell,flags)
            
    
            dirac_sink_nodes, dirac_sink_steady = ex_forcing.make_dirac_sink(
                dirac_sink_nodes,dirac_sink_steady,dirac_sinks,time,
                str(extra_sink),coord)
    
            if ( str(flag_source) == 'use') :
                cmd='cp '+ str(extra_source) +' '+str(filename_source)
                subprocess.call(cmd,shell=True)
                steady_sink=True
    
    
            steady_forcing=(
                steady_source and 
                steady_sink)
            steady_dirac_forcing=( 
                dirac_source_steady and 
                dirac_sink_steady)
    
    
            if ( flag_dirichlet =='no') :
                print( 'correction' )
                source_tria,sink_tria, dirac_source_nodes,dirac_sink_nodes=ex_forcing.correction(
                    topol,
                    source_tria,sink_tria,
                    dirac_source_nodes,dirac_sink_nodes,
                    size_cell)
    
                
    
            forcing_tria = source_tria-sink_tria
            dirac_forcing_nodes = dirac_source_nodes-dirac_sink_nodes

            #out_forcing_tria[:,0] = forcing_tria.flatten('C')
            #out_source_tria[:,0]  = source_tria.flatten('C')
            #out_sink_tria[:,0]    = sink_tria.flatten('C')

            timecell.write2file( 
                file_forcing,time,False,
                #out_forcing_tria,steady_forcing)
                forcing_tria,steady_forcing)
            timecell.write2file(
                file_source,time,False,
                #out_source_tria,steady_source)
                source_tria,steady_source)
            timecell.write2file(
                file_sink,time,False,
                #out_sink_tria,steady_sink)
                sink_tria,steady_sink)

            #out_dirac_forcing_nodes[:,0]=dirac_forcing_nodes.flatten('F')
            #out_dirac_source_nodes[:,0]=dirac_source_nodes.flatten('F')
            #out_dirac_sink_nodes[:,0]=dirac_sink_nodes.flatten('F')
            
            timecell.write2file(
                file_dirac_forcing,time,False,
                #out_dirac_forcing_nodes,steady_dirac_forcing)
                dirac_forcing_nodes,steady_dirac_forcing)
            timecell.write2file(
                file_dirac_source,time,False,
                #out_dirac_source_nodes,dirac_source_steady)
                dirac_source_nodes,dirac_source_steady)
            timecell.write2file(
                file_dirac_sink,time,False,
                #out_dirac_sink_nodes,dirac_sink_steady)
                dirac_sink_nodes,dirac_sink_steady)
        
            time=common.next_time(time)
            
            ### close files
            file_forcing.close()
            file_source.close()
            file_sink.close()
            file_dirac_forcing.close()
            file_dirac_source.close()
            file_dirac_sink.close()
    
            #####################
            # write last measure
            #out=np.zeros([ncell,4])
            #out[:,0:3]=bar_cell[:,:]
            
            #out[:,3]=source_tria[:][0]*size_cell[:]
            #cell.write2file(filename_measure_source,out)
        
            #out[:,3]=sink_tria[:][0]*size_cell[:]
            #cell.write2file(filename_measure_sink,out)
            
    if ( str(flag_source) == 'use') :
        cmd='cp '+ str(extra_source) +' '+str(filename_source)
        subprocess.call(cmd,shell=True)
    
    if ( str(flag_sink) == 'use') :    
        cmd='cp '+ str(extra_sink) +' '+str(filename_sink)
        subprocess.call(cmd,shell=True)
    
    
    
    
    
    
    # #################################################################
    # kappa
    if ( flag_kappa == 'use'):
        #
        # copy kappa
        #
        copyfile(extra_kappa,filename_kappa)
    else:
        if ( not kappas):
            kappas=[]
            kappas.append(timecell.example(flag_kappa,extra_kappa))
            
        # definitions 
        kappa_tria=np.zeros([ncell,1])
        steady=False
        time=tzero
        
        # open file and write dimensions
        file_out=open(filename_kappa, 'w')
        timecell.write2file(file_out,time,True,kappa_tria)
        
        #cycle in time
        while (time <= tmax) and (not steady ):
            # evalute 
            kappa_tria, steady = timecell.build(
                kappa_tria, steady, 
                kappas,time,bar_cell,flags)
            # write 2 file
            timecell.write2file(file_out,time,False,kappa_tria,steady)
        
            # next time
            time=common.next_time(time)
            file_out.close()
    
    #####################################################################
    # lambda (lift for tdens)
    #####################################################################
    # Default : lambda=0.0, no file is written.
    if ( flag_lambda == 'use'):
        #
        # copy lambda
        #
        copyfile(extra_lambda,filename_lambda)
    elif ( float(flag_lambda) == 0.0 ):
        #
        # set to zero
        #
        td.write_zero_timedata(filename_lambda,1,ncell)
        
    else:
        lambda_functions=[]
        lambda_functions.append(timecell.example(flag_lambda,extra_lambda))
        writetimedata(ncell,1,tzero,lambda_functions,bar_cell,flags,filename_lambda)
    
    
    #####################################################################
    # tdens0
    #####################################################################
    if ( flag_tdens0 == 'use'):
        #
        # copy lambda
        #
        copyfile(extra_tdens0,filename_tdens0)
    else:
        if ( not tdens0_functions):
            tdens0_functions=[]
            tdens0_functions.append(timecell.example(flag_tdens0,extra_tdens0))
    
            # definitions 
            tdens0_tria=np.zeros([ncell,1])
            steady=False
            time=tzero
            
            # open file and write dimensions
            file_out=open(filename_tdens0, 'w')
            timecell.write2file(file_out,time,True,tdens0_tria)
        
            #cycle in time
            while (time <= tmax) and (not steady ):
                # evalute 
                tdens0_tria, steady = timecell.build(
                    tdens0_tria, steady, 
                    tdens0_functions,time,bar_cell,flags)
                # write 2 file
                timecell.write2file(file_out,time,False,tdens0_tria,steady)
                
                # next time
                time=common.next_time(time)
            file_out.close()
    
    
    #####################################################################
    # Dirichlet
    #####################################################################
    # Write Zero file if there is no Dirichlet boundary condition
    # Otherwise use examples define in define_dirichlet_node
    # Only fixed in time are supported now
    if ( flag_dirichlet == 'no'):
        td.write_zero_timedata(filename_dirichlet,2*dimdata,nnode)
    else:
        dir_out=np.zeros([len(coord),2*dimdata]);
        isdirdof,nodedir,divalue=dirichlet.define_dirichlet_node(coord,topol,
                                                        flag_dirichlet,
                                                        extra_dirichlet,
                                                        flag_equation)

        if ( flag_equation == 'poisson' ):
            dir_out[nodedir[:],0]=1.0
            dir_out[nodedir[:],1]=divalue[:]
            td.write_steady_timedata(filename_dirichlet,dir_out)
        elif ( flag_equation == 'saintvenant' ):
            dir_out[nodedir[:,0],0]=isdirdof[:,0]
            dir_out[nodedir[:,0],1]=divalue[:,0]
            dir_out[nodedir[:,1],2]=isdirdof[:,1]
            dir_out[nodedir[:,1],3]=divalue[:,1]
            td.write_steady_timedata(filename_dirichlet,dir_out)
    

    #####################################################################
    # Neumann
    #####################################################################
    # Write Zero file if there is no Neumann boundary condition
    # Otherwise use examples define in define_neumann_edge
    # Only fixed in time are supported now
    #
    # 2d case : innz_edge n1 n2 value    
    # 3d case : innz_face n1 n2 n3 value 
    # with
    # poisson: value = value 
    # siant-venant: value = value x, value y
    nneumann=(nnodeincell-1)+dimdata
    print('nnodeincell,dimdata',nnodeincell,dimdata)
    if ( flag_neumann == 'no'):
        td.write_zero_timedata(filename_neumann,nneumann,nedge)
    else:
        neu_out=np.zeros([nedge,nneumann]);
        neuedge,neunode,neuvalue=neumann.define_neumann_edge(flag_neumann,
                                                                extra_neumann,
                                                                coord,topol,
                                                                all_edges)
        if ( flag_equation == 'saintvenant' ):
            neu_out[0:neuedge,0:2]=neunode[:,:]
            neu_out[0:neuedge,2:4]=neuvalue[:,:]
            td.write_steady_timedata(filename_neumann,neu_out)
        else:
            print ('Only Saint-Venant equations' )
            sys.exit()
    
    
    #####################################################################
    # OPTDENS optimal transpost density
    #####################################################################
    if ( not optdens_functions):
        optdens_functions=[]
    try:
        if (len(optdens_functions) == 0):
            optdens_functions.append(make_optdens.define_optimal_tdens(
                flag_grid, flag_source, flag_sink, 
                extra_grid, extra_source, extra_sink ))
            
        optdens_tria=np.zeros([ncell,1])
        steady=False
        time=tzero
        
        # open file and write dimensions
        print(filename_optdens)
        file_out=open(filename_optdens, 'w')
        timecell.write2file(file_out,time,True,optdens_tria)
        
        #cycle in time
        while (time <= tmax) and (not steady ):
            # evalute
            optdens_tria, steady = timecell.build(
                optdens_tria, steady, 
                optdens_functions,time,bar_cell,flags)
            # write 2 file
            timecell.write2file(file_out,time,False,optdens_tria,steady)
            
            # next time
            time=common.next_time(time)
        file_out.close()
    except:
        print('Do optimal tdens defined')
        
    
    
    
    #######################################################################
    ##########################
    # pflux
    scalars.build_and_write(flag_pflux,extra_pflux,
                        tzero,tmax,
                        filename_pflux)
    
    # pode
    scalars.build_and_write(flag_pode,extra_pode,
                        tzero,tmax,
                        filename_pode)
    
    # pmass
    scalars.build_and_write(flag_pmass,extra_pmass,
                        tzero,tmax,
                        filename_pmass)
    
    # decay
    scalars.build_and_write(flag_decay,extra_decay,
                        tzero,tmax,
                        filename_decay)

    return;







if __name__ == "__main__":
    if len(sys.argv) > 1:
        f_ctrl=sys.argv[1]
        dir_inputs=sys.argv[2]
        dir_inputs_vtk=sys.argv[3]

        assembly_inputs(f_ctrl,dir_inputs,dir_inputs_vtk)
    else:
        raise SystemExit("usage: python assembly_inputs.py 'input ctrl file' 'directory_dat' 'directory_vtk'"  )
