#!/usr/bin/env python
import numpy as np
import common
import meshtools_surface as mt


# assembly procedure 
def make_source(source_cell,steady,sources,time,flags,bar_cell):
    f_cell=0.0
    steady=True
    source_cell[:]=0.0
    for icell in range(len(bar_cell)):
        for ifun in range(len(sources)):
            f_cell,steady_cell=sources[ifun](
                time,
                bar_cell[icell][0],
                bar_cell[icell][1],
                bar_cell[icell][2],
                flags[icell])
            source_cell[icell][:]+=float(f_cell)
            steady=steady and steady_cell
    return source_cell,steady;

def make_sink(sink_cell,steady,sinks,time,flags,bar_cell):
    f_cell=0.0
    steady=True
    sink_cell[:]=0
    for icell in range(len(bar_cell)):
        for ifun in range(len(sinks)):
            f_cell,steady_cell=sinks[ifun](
                time,
                bar_cell[icell,0],
                bar_cell[icell][1],
                bar_cell[icell][2],
                flags[icell])
            sink_cell[icell][:]+=float(f_cell)
            steady=steady and steady_cell
    return sink_cell,steady;


# assembly procedure 
def make_dirac_source(source_nodes,steady,dirac_sources,time,path_file,coord):
    steady_global=True
    source_nodes[:]=0
    for fdirac in dirac_sources:
        source_coord,values, steady=fdirac(time,path_file)
        steady_global = steady_global and steady
        print len(coord)
        for idirac in range(len(source_coord)):
            print source_coord[idirac]
            inode,dist=mt.FindClosestNode(
                range(len(coord)),coord,source_coord[idirac])
            source_nodes[inode[0]][:] += values[idirac]
            
    return source_nodes, steady_global;

# assembly procedure 
def make_dirac_sink(sink_nodes,steady,dirac_sinks,time,path_file,coord):
    steady_global=True
    sink_nodes[:]=0
    for fdirac in dirac_sinks:
        sink_coord=[]
        sink_values=[]
        sink_coord,values, steady=fdirac(time,path_file)
        steady_global = steady_global and steady
        for idirac in range(len(sink_coord)):
            inode,dist=mt.FindClosestNode(
                range(len(coord)),coord,sink_coord[idirac])
            sink_nodes[inode[0]][:] += values[idirac]
    return sink_nodes, steady_global;



def correction(topol,
               source_cell,sink_cell, 
               dirac_source_values, dirac_sink_values,
               size_cell, tol=1.0e-10):

    for icell in range(len(sink_cell)):
        if (source_cell[icell,0] * sink_cell[icell,0] != 0.0):
            print icell,source_cell[icell,0], sink_cell[icell,0]
            sys.exit('source and sink have no disjoint support')

    rhs = assembly_rhs(topol,
                       source_cell-sink_cell,
                       dirac_source_values-dirac_sink_values,size_cell)
    rhs_source = assembly_rhs(topol,
                       source_cell,
                       dirac_source_values,size_cell)
    rhs_sink  = assembly_rhs(topol,
                       sink_cell,
                       dirac_sink_values,size_cell)

    plus =np.sum(rhs_source)
    minus=np.sum(rhs_sink)
    imbalance =np.sum(rhs_source-rhs_sink)
        
    if ( abs(imbalance) >= tol):
        print 'imbalance=',imbalance, sum(rhs,all>0.0),sum(rhs,all<0.0)
        correct_factor = plus/minus 
        print 'correction factor =',correct_factor
        sink_cell=sink_cell*correct_factor
        #dirac_sink_values=dirac_sink_values*correct_factor
        rhs_sink  = assembly_rhs(topol,
                       sink_cell,
                       dirac_sink_values,size_cell)
        print 'imbalance=',np.sum(rhs_source-rhs_sink)


    return source_cell,sink_cell, dirac_source_values, dirac_sink_values;

def assembly_rhs(topol,lebesgue, dirac,cell_size):
    rhs=np.zeros(dirac.shape[0])
    nnodeincell=topol.shape[1]
    for icell in range(len(lebesgue)):
        for node in topol[icell,:]:
            rhs[node]=rhs[node] + lebesgue[icell][0] * cell_size[icell] / nnodeincell

    
    for node in range(len(dirac)):
        if (dirac[node][0] != 0.0):
            rhs[node] = dirac[node][0]
    return (rhs)
        



def normalize(source_cell,sink_cell, 
              dirac_source_values, dirac_sink_values,
              size_cell):
    fplus=np.dot(source_cell,size_cell)+sum(dirac_source_values)
    fminus=np.dot(sink_cell,size_cell)+sum(dirac_sink_values)
    source_cell=source_cell/fplus
    dirac_source_values=dirac_source_values/fplus
    sink_cell=sink_cell/fminus
    dirac_sink_values=dirac_sink_values/fminus
    return source_cell,sink_cell, dirac_source_values, dirac_sink_values;

# build 
def example_source(sources,dirac_sources,flag_source,extra_info):
    #  personal
    if (flag_source == 'cnst'):
        def source(time,x,y,z,flag_mesh):
            fvalue=0.0
            steady=True
            if (flag_mesh > 0 ):
                fvalue=2.0
            return fvalue, steady;
        sources.append(source)

    #  experiment for dmk of surface
    if (flag_source == 'north_band'):
        def source(time,x,y,z,flag_mesh):
            fvalue=0.0
            steady=True
            if ( (x >= 0.0) and (x<=1.0 ) and 
                 (y >= 0.0) and (y<=1.0 ) and
                 (z >= 0.5) and (z<=np.sqrt(3)/2.0) ) :
                fvalue=1.0
            return fvalue, steady;
        sources.append(source)

    #  personal
    if (flag_source == 'uniform'):
        def source(time,x,y,z,flag_mesh):
            fvalue=0.0
            steady=True
            fvalue=1.0
            return fvalue, steady;
        sources.append(source)

    
    #  sin
    if (flag_source == 'sin'):
        def source(time,x,y,z,flag_mesh):
            fvalue=0.0
            if (flag_mesh > 0 ):
                fvalue=np.sin(2.0*y)
            steady=True
            return fvalue,steady;
        sources.append(source)
            
    #  rectangles constat
    if (flag_source == 'cnst_rect'):
        def source(time,x,y,z,flag_mesh):
            fvalue=0.0
            steady=True
            if ( (x >= 1.0/8.0) and (x<=3.0/8.0) and 
                 (y >= 1.0/4.0) and (y<=3.0/4.0) ) :
                fvalue=2.0
            return fvalue,steady;
        sources.append(source)

    #  y branch problem
    if (flag_source == 'ybranch'):
        # read points
        info_path=common.remove_comments(extra_info,'!')
        def dirac_source(t,info_path):
            coord_points=[]
            values=[]
            steady=True
            input_file=open(info_path,'r')
            input_lines = input_file.readlines()
            for line in input_lines[0:3]:
                coord_dirac = [float(w) for w in line.split()[0:2]]
                branch_value = float(line.split()[2])
                if ( branch_value > 0.0 ):
                    values.append(branch_value)
                    coord_points.append(coord_dirac)                    
            input_file.close()
            return coord_points, values,steady;
        dirac_sources.append(dirac_source)
        
    #  y branch problem
    if (flag_source == 'dirac'):
        # read points
        info_path=common.remove_comments(extra_info,'!')
        def dirac_source(t,info_path):
            coord_points=[]
            values=[]
            steady=True
            input_file=open(info_path,'r')
            input_lines = input_file.readlines()
            for line in input_lines[0:]:
                coord = [float(w) for w in line.split()[0:3]]
                value = float(line.split()[3])
                # check for positive values only
                if ( value > 0.0 ):
                    values.append(value)
                    coord_points.append(coord)                    
            input_file.close()
            return coord_points, values,steady;
        dirac_sources.append(dirac_source)


    return sources, dirac_sources;

# build standard tdens_0  given id_tdens0
def example_sink(sinks,dirac_sinks,flag_sink,extra_info):    

    print flag_sink
    #  constant examples
    if (flag_sink == 'cnst'):
        def sink(time,x,y,z,flag_mesh):
            fvalue=0.0
            steady=True
            if (flag_mesh < 0 ):
                fvalue=2.0
            return fvalue,steady;
        sinks.append(sink)

    #  personal
    if (flag_sink == 'uniform'):
        
        def sink(time,x,y,z,flag_mesh):
            fvalue=0.0
            steady=True
            fvalue=1.0
            return fvalue, steady;
        sinks.append(sink)

        #  experiment for dmk of surface
    if (flag_sink == 'south_band'):
        def source(time,x,y,z,flag_mesh):
            fvalue=0.0
            steady=True
            r=np.sqrt(x**2+y**2+z**2)
            if ( (x >= 0.0  ) and (x<=1.0 ) and 
                 (y >= 0.0  ) and (y<=1.0 ) and
                 (z >= -np.sqrt(3)/2.0) and (z<=-0.5) ) :
                fvalue=1.0
            return fvalue, steady;
        sinks.append(source)



     
    #  y branch problem
    if (flag_sink == 'ybranch'):
        # read points
        file_path=common.remove_comments(extra_info,'!')
        def dirac_sink(time, file_path):
            coord_points=[]
            values=[]
            input_file  = open(file_path,'r')
            input_lines = input_file.readlines()
            for line in input_lines[0:3]:
                coord_dirac = [float(w) for w in line.split()[0:2]]
                branch_value = float(line.split()[2])
                if ( branch_value < 0.0 ):
                    values.append(branch_value)
                    coord_points.append(coord_dirac)
            input_file.close()
            steady=True

            return coord_points, values, steady;
        dirac_sinks.append(dirac_sink)
        
    #  y branch problem
    if (flag_sink == 'dirac'):
        # read points
        info_path=common.remove_comments(extra_info,'!')
        def dirac_sink(t,info_path):
            coord_points=[]
            values=[]
            steady=True
            input_file=open(info_path,'r')
            input_lines = input_file.readlines()
            for line in input_lines[0:]:
                coord = [float(w) for w in line.split()[0:2]]
                value = float(line.split()[2])
                # check for positive values only
                if ( value > 0.0 ):
                    values.append(value)
                    coord_points.append(coord)                    
            input_file.close()
            return coord_points, values,steady;
        dirac_sinks.append(dirac_sink)


    return sinks, dirac_sinks;    

    
