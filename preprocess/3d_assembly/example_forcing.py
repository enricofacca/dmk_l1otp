#!/usr/bin/env python
import numpy as np
import common
import meshtools as mt


# assembly procedure 
def make_source(source_cell,steady,sources,time,flags,bar_cell):
    f_cell=0.0
    steady=True
    source_cell[:]=0.0
    for icell in range(len(bar_cell)):
        for ifun in range(len(sources)):
            f_cell,steady_cell=sources[ifun](
                time,
                bar_cell[icell][:],
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
                bar_cell[icell][:],
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
        for idirac in range(len(source_coord)):
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
    if (flag_source == 'rect_cnst'):
        def source(time,coord,flag_mesh):
            fvalue=0.0
            x=coord[0]
            y=coord[1]
            z=coord[2]
            steady=True
            if (flag_mesh > 0 ):
                fvalue=2.0
            return fvalue, steady;
        sources.append(source)

   
    if (flag_source == 'convergence'):
        def source(time,coord,flag_mesh):
            fvalue=0.0
            x=coord[0]
            y=coord[1]
            z=coord[2]
            fvalue=3.0*np.pi*np.cos(np.pi*x)*np.cos(np.pi*y)*np.cos(np.pi*z)
            steady=True
            if (fvalue < 0.0 ):
                fvalue=0.0
            return fvalue, steady;
        sources.append(source)
            
    if (flag_source == 'tosi'):
        def source(time,coord,flag_mesh):
            fvalue=0.0
            x=coord[0]
            y=coord[1]
            z=coord[2]
            fvalue=3.0*(x**2+y**2+z**2)
            steady=True
            if (fvalue < 0.0 ):
                fvalue=0.0
            return fvalue, steady;
        sources.append(source)

    if (flag_source == 'sphere_plaplacian'):
        # define sphere radius
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
        
        # cnst_value=1.0
        volume_source = (4.0/3.0)*np.pi*radius_s**3
        volume_sink   = (4.0/3.0)*np.pi*radius_l**3 - (4.0/3.0)*np.pi*radius_m**3
        cnst_value = 1.0 / volume_source
        cnst_sink_value = 1.0 / volume_sink
        # cnst_sink_value = (volume_source*cnst_value)/volume_sink
        # cnst_value = (volume_sink*cnst_sink_value)/volume_source
        c1 = cnst_value
        c2 = - cnst_sink_value
        
        print 'Sink/Source - 1/19:',np.abs(c2/c1+1.0/19.0)
        
        def source(time,coord,flag_mesh):
            fvalue=0.0
            x=coord[0]
            y=coord[1]
            z=coord[2]
            steady=True
            r_dist = np.sqrt(x**2+y**2+z**2)
            if (flag_mesh > 0):
                fvalue=cnst_value
            return fvalue, steady;
        sources.append(source)


    return sources, dirac_sources;

# build standard tdens_0  given id_tdens0
def example_sink(sinks,dirac_sinks,flag_sink,extra_info):    
    #  personal
    if (flag_sink == 'rect_cnst'):
        def sink(time,coord,flag_mesh):
            fvalue=0.0
            x=coord[0]
            y=coord[1]
            z=coord[2]
            fvalue=0.0
            steady=True
            if (flag_mesh < 0 ):
                fvalue=2.0
            return fvalue,steady;
        sinks.append(sink)
    
    if (flag_sink == 'convergence'):
        def sink(time,coord,flag_mesh):
            fvalue=0.0
            x=coord[0]
            y=coord[1]
            z=coord[2]
            fvalue=3.0*np.pi*np.cos(np.pi*x)*np.cos(np.pi*y)*np.cos(np.pi*z)
            steady=True
            if (fvalue > 0.0 ):
                fvalue=0.0
            return fvalue, steady;
        sinks.append(sink)


    if (flag_sink == 'sphere_plaplacian'):
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
        print 'Sink/Source - 1/19:',np.abs(c2/c1+1.0/19.0)
        

        def sink(time,coord,flag_mesh):
            fvalue=0.0
            x=coord[0]
            y=coord[1]
            z=coord[2]
            steady=True
            r_dist = np.sqrt(x**2+y**2+z**2)
            if (flag_mesh < 0):
                fvalue=cnst_sink_value
            return fvalue, steady;
        sinks.append(sink)
        


    return sinks, dirac_sinks;              


    
