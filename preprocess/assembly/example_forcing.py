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
                bar_cell[icell,:],
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
                bar_cell[icell,:],
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
            source_nodes[inode[0],:][:] += values[idirac]
            
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
            sink_cell[icell,0]=0.0
            #print icell,source_cell[icell,0], sink_cell[icell,0]
            #sys.exit('source and sink have no disjoint support')

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
        print ('imbalance=',imbalance)
        correct_factor = plus/minus 
        print ('correction factor =',correct_factor)
        sink_cell=sink_cell*correct_factor
        #dirac_sink_values=dirac_sink_values*correct_factor
        rhs_sink  = assembly_rhs(topol,
                       sink_cell,
                       dirac_sink_values,size_cell)
        print ('imbalance=',np.sum(rhs_source-rhs_sink))


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
def example_source(sources,dirac_sources,flag_source,extra_info,grid_topology=None,grid_coordinates=None):
    #  personal
    if (flag_source == 'cnst'):
        def source(time,coord,flag_mesh):
            fvalue=0.0
            steady=True
            if (flag_mesh > 0 ):
                fvalue=2.0
            return fvalue, steady;
        sources.append(source)

    #  personal
    if (flag_source == 'uniform'):
        def source(time,coord,flag_mesh):
            fvalue=0.0
            steady=True
            fvalue=1.0
            return fvalue, steady;
        sources.append(source)

     #  personal
    if (flag_source == 'zero'):
        def source(time,coord,flag_mesh):
            fvalue=0.0
            steady=True
            fvalue=0.0
            return fvalue, steady;
        sources.append(source)


    if (flag_source == 'plaplacian'):
        # Source term used for testing convergence for
        # plaplacian equation in  
        # Branching Structures Emerging From A
        #Continuous Optimal Transport Model
        # Enrico Facca, Franco Cardin, Mario Putti
        def source(time,coord,flag_mesh):
            fvalue=0.0
            steady=True
            if (np.sqrt(x**2+y**2)< 1.0/3.0):
                fvalue=1
            return fvalue, steady;
        sources.append(source)

    if (flag_source == 'plaplacian_dirichlet'):
        # Source term used for testing convergence for
        # plaplacian equation in  
        # Branching Structures Emerging From A
        #Continuous Optimal Transport Model
        # Enrico Facca, Franco Cardin, Mario Putti
        def source(time,coord,flag_mesh):
            fvalue=1.0
            steady=True
            return fvalue, steady;
        sources.append(source)
        
    
    #  A parallel fast sweeping method for the Eikonal equation
    #  Miles Detrixhea, Frederic Giboua,b, Chohong Min
    if (flag_source == 'example1fmt'):
        def dirac_source(t,info_path):
            coord_points=[[0.0,0.0,0.0]]
            values=[1.0]
            steady=True
            return coord_points, values,steady;
        dirac_sources.append(dirac_source)

    if (flag_source == 'one2two'):
        def dirac_source(t,info_path):
            coord_points=[[0.5,0.1,0.0]]
            values=[1.0]
            steady=True
            return coord_points, values,steady;
        dirac_sources.append(dirac_source)


    # An accurate discontinuous Galerkin method for solving point-source
    # Eikonal equation in 2-D heterogeneous anisotropic media
    # Bouteiller,1 M. Benjemaa,2 L. Metivier1,3 and J. Virieux1
    if (flag_source == 'bouteiller18_tc1'):
        def dirac_source(t,info_path):
            coord_points=[[2.0,2.0,0.0]]
            values=[1.0]
            steady=True
            return coord_points, values,steady;
        dirac_sources.append(dirac_source)

    if (flag_source == 'three_forces_0'):
        # set point at vertices of equilateral triangles
        coord = [
            [np.cos(np.pi/2          ),np.sin(np.pi/2          )],
            [np.cos(np.pi/2+2*np.pi/3),np.sin(np.pi/2+2*np.pi/3)],
            [np.cos(np.pi/2+4*np.pi/3),np.sin(np.pi/2+4*np.pi/3)],
        ]        
        coord=np.asarray(coord)
        values=np.where(coord>0,coord,0)
        coord[:,:]=0.3*coord[:,:]
        coord[:,:]=coord+[0.5,0.5]

        coord=np.hstack((coord,np.zeros([3,1])))

        def dirac_source(t,info_path):
            steady=True
            return coord, values,steady;
        dirac_sources.append(dirac_source)

    if (flag_source == 'three_forces_1'):
        # set point at vertices of equilateral triangles
        coord = [
            [np.cos(3*np.pi/4),np.sin(3*np.pi/4)],
            [np.cos(5*np.pi/4),np.sin(5*np.pi/4)],
            [np.cos(7*np.pi/4),np.sin(7*np.pi/4)],
        ]        
        coord=np.asarray(coord)

        # Esempio 1
        values=[
            [+np.sqrt(2)/4,+np.sqrt(2)/4],
            [-np.sqrt(2)/2,-np.sqrt(2)/2],
            [+np.sqrt(2)/4,+np.sqrt(2)/4],
        ]
        values=np.asarray(values)
        values=np.where(values>0,values,0)

        coord[:,:]=0.3*coord[:,:]
        coord[:,:]=coord+[0.5,0.5]
        coord=np.hstack((coord,np.zeros([3,1])))

        def dirac_source(t,info_path):
            steady=True
            return coord, values,steady;
        dirac_sources.append(dirac_source)

    if (flag_source == 'three_forces_2'):
        # set point at vertices of equilateral triangles
        coord = [
            [np.cos(3*np.pi/4),np.sin(3*np.pi/4)],
            [np.cos(5*np.pi/4),np.sin(5*np.pi/4)],
            [np.cos(7*np.pi/4),np.sin(7*np.pi/4)],
        ]        
        coord=np.asarray(coord)

        # Esempio 2
        values=[
            [-0.5,-0.5],
            [+0.0,+1.0],
            [+0.5,-0.5],
        ]
        values=np.asarray(values)
        values=np.where(values>0,values,0)

        coord[:,:]=0.3*coord[:,:]
        coord[:,:]=coord+[0.5,0.5]
        coord=np.hstack((coord,np.zeros([3,1])))

        def dirac_source(t,info_path):
            steady=True
            return coord, values,steady;
        dirac_sources.append(dirac_source)

    if (flag_source == 'three_forces_3'):
        # set point at vertices of equilateral triangles
        coord = [
            [np.cos(3*np.pi/4),np.sin(3*np.pi/4)],
            [np.cos(5*np.pi/4),np.sin(5*np.pi/4)],
            [np.cos(7*np.pi/4),np.sin(7*np.pi/4)],
        ]        
        coord=np.asarray(coord)

        # Esempio 3
        values=[
            [-0.5,+0.5*(4-np.sqrt(3))],
            [-0.5*(np.sqrt(3)-1),-0.5*(3-np.sqrt(3))],
            [+0.5*np.sqrt(3),-0.5],
        ]
        values=np.asarray(values)
        values=np.where(values>0,values,0)

        coord[:,:]=0.3*coord[:,:]
        coord[:,:]=coord+[0.5,0.5]
        coord=np.hstack((coord,np.zeros([3,1])))

        def dirac_source(t,info_path):
            steady=True
            return coord, values,steady;
        dirac_sources.append(dirac_source)

    if (flag_source == 'three_forces_4'):
        # set point at vertices of equilateral triangles
        coord = [
            [np.cos(3*np.pi/4),np.sin(3*np.pi/4)],
            [np.cos(5*np.pi/4),np.sin(5*np.pi/4)],
            [np.cos(7*np.pi/4),np.sin(7*np.pi/4)],
        ]        
        coord=np.asarray(coord)

        # Esempio 4
        values=[
            [+1*np.sin(np.pi/12),-1*np.cos(np.pi/12)],
            [-0.5*np.sqrt(2),0.5*np.sqrt(2)],
            [0.5*np.sqrt(2)*(np.sqrt(3)-1)*np.cos(np.pi/6),
             0.5*np.sqrt(2)*(np.sqrt(3)-1)*np.sin(np.pi/6)],
        ]
        values=np.asarray(values)
        values=np.where(values>0,values,0)

        coord[:,:]=0.3*coord[:,:]
        coord[:,:]=coord+[0.5,0.5]
        coord=np.hstack((coord,np.zeros([3,1])))

        def dirac_source(t,info_path):
            steady=True
            return coord, values,steady;
        dirac_sources.append(dirac_source)

    if (flag_source == 'three_forces_5'):
        # set point at vertices of equilateral triangles
        coord = [
            [np.cos(3*np.pi/4),np.sin(3*np.pi/4)],
            [np.cos(5*np.pi/4),np.sin(5*np.pi/4)],
            [np.cos(7*np.pi/4),np.sin(7*np.pi/4)],
        ]        
        coord=np.asarray(coord)

        # Esempio 5
        values=[
            [-0.5*np.sqrt(3),0.5*np.sqrt(3)*(2-np.sqrt(3))],
            [0.5*(np.sqrt(3)-1),0.5*np.sqrt(3)*(np.sqrt(3)-1)],
            [0.5,-0.5*np.sqrt(3)],
        ]
        values=np.asarray(values)
        values=np.where(values>0,values,0)

        coord[:,:]=0.3*coord[:,:]
        coord[:,:]=coord+[0.5,0.5]
        coord=np.hstack((coord,np.zeros([3,1])))

        def dirac_source(t,info_path):
            steady=True
            return coord, values,steady;
        dirac_sources.append(dirac_source)

    if (flag_source == 'three_forces_6'):
        # set point at vertices of equilateral triangles
        coord = [
            [np.cos(3*np.pi/4),np.sin(3*np.pi/4)],
            [np.cos(5*np.pi/4),np.sin(5*np.pi/4)],
            [np.cos(7*np.pi/4),np.sin(7*np.pi/4)],
        ]        
        coord=np.asarray(coord)

        # Esempio 6
        values=[
            [0.5*np.sqrt(3),+0.5-1/np.sqrt(3)],
            [-0.5*(np.sqrt(3)+1),-(3+np.sqrt(3))/6],
            [0.5,0.5*np.sqrt(3)],
        ]
        values=np.asarray(values)
        values=np.where(values>0,values,0)

        coord[:,:]=0.3*coord[:,:]
        coord[:,:]=coord+[0.5,0.5]
        coord=np.hstack((coord,np.zeros([3,1])))

        def dirac_source(t,info_path):
            steady=True
            return coord, values,steady;
        dirac_sources.append(dirac_source)

    if (flag_source == 'example3'):

        def dirac_source(t,info_path):

            coord_source = []
            values = []

            for i in range(len(grid_coordinates)):
                x = grid_coordinates[i,0]
                y = grid_coordinates[i,1]
                if abs(y-0.4) < 1e-12 :
                    if abs(x-1.2) < 1e-12 :
                        coord_source.append([x,y,0])
                        values.append([0,-1])
            values = np.asarray(values)
            values = np.where(values>0,values,0)        
            steady=True
            return coord_source, values,steady;
        dirac_sources.append(dirac_source)

    if (flag_source == 'example3_semi_inside'):

        def dirac_source(t,info_path):

            coord_source = []
            values = []

            for i in range(len(grid_coordinates)):
                x = grid_coordinates[i,0]
                y = grid_coordinates[i,1]
                if abs(y-0.4) < 1e-12 :
                    if abs(x-1.2) < 1e-12 :
                        coord_source.append([x,y,0])
                        values.append([0,-1])
            values = np.asarray(values)
            values = np.where(values>0,values,0)        
            steady=True
            return coord_source, values,steady;
        dirac_sources.append(dirac_source)

    if (flag_source == 'example3_inside'):

        def dirac_source(t,info_path):

            coord_source = []
            values = []

            for i in range(len(grid_coordinates)):
                x = grid_coordinates[i,0]
                y = grid_coordinates[i,1]
                if abs(y-1.2) < 1e-12 :
                    if abs(x-1.5) < 1e-12 :
                        coord_source.append([x,y,0])
                        values.append([0,-1])
            values = np.asarray(values)
            values = np.where(values>0,values,0)        
            steady=True
            return coord_source, values,steady;
        dirac_sources.append(dirac_source)

    if (flag_source == 'beam'):

        def dirac_source(t,info_path):

            coord_source = []
            values = []

            for i in range(len(grid_coordinates)):
                x = grid_coordinates[i,0]
                y = grid_coordinates[i,1]
                if abs(y-1.2) < 1e-12 :
                    if abs(x-1.5) < 1e-12 :
                        coord_source.append([x,y,0])
                        values.append([0,-1])
            values = np.asarray(values)
            values = np.where(values>0,values,0)        
            steady=True
            return coord_source, values,steady;
        dirac_sources.append(dirac_source)

    if (flag_source == 'beam_shear'):

        def dirac_source(t,info_path):

            coord_source = []
            values = []

            for i in range(len(grid_coordinates)):
                x = grid_coordinates[i,0]
                y = grid_coordinates[i,1]
                if abs(x-2) < 1e-12 :
                    if (y > - 1e-12) and (y < 0.3 + 1e-12) :
                        coord_source.append([x,y,0])
                        values.append([0,-0.01])
            values = np.asarray(values)
            values = np.where(values>0,values,0)        
            steady=True
            return coord_source, values,steady;
        dirac_sources.append(dirac_source)

    if (flag_source == 'beam_shear_inside'):

        def dirac_source(t,info_path):

            coord_source = []
            values = []

            for i in range(len(grid_coordinates)):
                x = grid_coordinates[i,0]
                y = grid_coordinates[i,1]
                if abs(x-2.5) < 1e-12 :
                    if (y > 0.35 - 1e-12) and (y < 0.65 + 1e-12) :
                        coord_source.append([x,y,0])
                        values.append([0,-1])
            values = np.asarray(values)
            values = np.where(values>0,values,0)        
            steady=True
            return coord_source, values,steady;
        dirac_sources.append(dirac_source)

    if (flag_source == 'example2_concentrated_load'):

        def dirac_source(t,info_path):

            coord_source = []
            values = []

            for i in range(len(grid_coordinates)):
                x = grid_coordinates[i,0]
                y = grid_coordinates[i,1]
                if abs(x-2.0) < 1e-12 :
                    if abs(y-0.15) < 1e-12 :
                        coord_source.append([x,y,0])
                        values.append([0,-1])
            values = np.asarray(values)
            values = np.where(values>0,values,0)        
            steady=True
            return coord_source, values,steady;
        dirac_sources.append(dirac_source)

    if (flag_source == 'example2_concentrated_load_inside'):

        def dirac_source(t,info_path):

            coord_source = []
            values = []

            for i in range(len(grid_coordinates)):
                x = grid_coordinates[i,0]
                y = grid_coordinates[i,1]
                if abs(x-2.5) < 1e-12 :
                    if abs(y-0.5) < 1e-12 :
                        coord_source.append([x,y,0])
                        values.append([0,-1])
            values = np.asarray(values)
            values = np.where(values>0,values,0)        
            steady=True
            return coord_source, values,steady;
        dirac_sources.append(dirac_source)

    if (flag_source == 'example1_concentrated_load_compression'):

        def dirac_source(t,info_path):

            coord_source = []
            values = []

            for i in range(len(grid_coordinates)):
                x = grid_coordinates[i,0]
                y = grid_coordinates[i,1]
                if abs(y-2.0) < 1e-12 :
                    if abs(x-0.15) < 1e-12 :
                        coord_source.append([x,y,0])
                        values.append([0,-1])
            values = np.asarray(values)
            values = np.where(values>0,values,0)        
            steady=True
            return coord_source, values,steady;
        dirac_sources.append(dirac_source)

    if (flag_source == 'example1_concentrated_load_traction'):

        def dirac_source(t,info_path):

            coord_source = []
            values = []

            for i in range(len(grid_coordinates)):
                x = grid_coordinates[i,0]
                y = grid_coordinates[i,1]
                if abs(y-2.0) < 1e-12 :
                    if abs(x-0.15) < 1e-12 :
                        coord_source.append([x,y,0])
                        values.append([0,+1])
            values = np.asarray(values)
            values = np.where(values>0,values,0)        
            steady=True
            return coord_source, values,steady;
        dirac_sources.append(dirac_source)

    if (flag_source == 'example1_concentrated_load_compression_inside'):

        def dirac_source(t,info_path):

            coord_source = []
            values = []

            for i in range(len(grid_coordinates)):
                x = grid_coordinates[i,0]
                y = grid_coordinates[i,1]
                if abs(y-2.5) < 1e-12 :
                    if abs(x-0.5) < 1e-12 :
                        coord_source.append([x,y,0])
                        values.append([0,-1])
            values = np.asarray(values)
            values = np.where(values>0,values,0)        
            steady=True
            return coord_source, values,steady;
        dirac_sources.append(dirac_source)

    if (flag_source == 'example1_concentrated_load_traction_inside'):

        def dirac_source(t,info_path):

            coord_source = []
            values = []

            for i in range(len(grid_coordinates)):
                x = grid_coordinates[i,0]
                y = grid_coordinates[i,1]
                if abs(y-2.5) < 1e-12 :
                    if abs(x-0.5) < 1e-12 :
                        coord_source.append([x,y,0])
                        values.append([0,+1])
            values = np.asarray(values)
            values = np.where(values>0,values,0)        
            steady=True
            return coord_source, values,steady;
        dirac_sources.append(dirac_source)

    if (flag_source == 'column_traction'):

        def dirac_source(t,info_path):

            coord_source = []
            values = []

            for i in range(len(grid_coordinates)):
                x = grid_coordinates[i,0]
                y = grid_coordinates[i,1]
                if abs(y-2.0) < 1e-12 :
                    if (x > - 1e-12) and (x < 0.3 + 1e-12) :
                        coord_source.append([x,y,0])
                        values.append([0,1])
            values = np.asarray(values)
            values = np.where(values>0,values,0)        
            steady=True
            return coord_source, values,steady;
        dirac_sources.append(dirac_source)

    if (flag_source == 'column_traction_inside'):

        def dirac_source(t,info_path):

            coord_source = []
            values = []

            for i in range(len(grid_coordinates)):
                x = grid_coordinates[i,0]
                y = grid_coordinates[i,1]
                if abs(y-2.5) < 1e-12 :
                    if (x > 0.35 - 1e-12) and (x < 0.65 + 1e-12) :
                        coord_source.append([x,y,0])
                        values.append([0,1])
            values = np.asarray(values)
            values = np.where(values>0,values,0)        
            steady=True
            return coord_source, values,steady;
        dirac_sources.append(dirac_source)

    if (flag_source == 'column_compression'):

        def dirac_source(t,info_path):

            coord_source = []
            values = []

            for i in range(len(grid_coordinates)):
                x = grid_coordinates[i,0]
                y = grid_coordinates[i,1]
                if abs(y-2.0) < 1e-12 :
                    if (x > - 1e-12) and (x < 0.3 + 1e-12) :
                        coord_source.append([x,y,0])
                        values.append([0,-1])
            values = np.asarray(values)
            values = np.where(values>0,values,0)        
            steady=True
            return coord_source, values,steady;
        dirac_sources.append(dirac_source)

    if (flag_source == 'column_compression_inside'):

        def dirac_source(t,info_path):

            coord_source = []
            values = []

            for i in range(len(grid_coordinates)):
                x = grid_coordinates[i,0]
                y = grid_coordinates[i,1]
                if abs(y-2.5) < 1e-12 :
                    if (x > 0.35 - 1e-12) and (x < 0.65 + 1e-12) :
                        coord_source.append([x,y,0])
                        values.append([0,-1])
            values = np.asarray(values)
            values = np.where(values>0,values,0)        
            steady=True
            return coord_source, values,steady;
        dirac_sources.append(dirac_source)

    if (flag_source == 'example_internet1'):

        def dirac_source(t,info_path):

            coord_source = []
            values = []

            for i in range(len(grid_coordinates)):
                x = grid_coordinates[i,0]
                y = grid_coordinates[i,1]
                if abs(y-0.0) < 1e-12 :
                    if abs(x-1.5) < 1e-12 :
                        coord_source.append([x,y,0])
                        values.append([0,-1])
            values = np.asarray(values)
            values = np.where(values>0,values,0)        
            steady=True
            return coord_source, values,steady;
        dirac_sources.append(dirac_source)

    if (flag_source == 'example_internet1_semi_inside'):

        def dirac_source(t,info_path):

            coord_source = []
            values = []

            for i in range(len(grid_coordinates)):
                x = grid_coordinates[i,0]
                y = grid_coordinates[i,1]
                if abs(y-0.0) < 1e-12 :
                    if abs(x-1.5) < 1e-12 :
                        coord_source.append([x,y,0])
                        values.append([0,-1])
            values = np.asarray(values)
            values = np.where(values>0,values,0)        
            steady=True
            return coord_source, values,steady;
        dirac_sources.append(dirac_source)

    if (flag_source == 'example_internet2'):

        def dirac_source(t,info_path):

            coord_source = []
            values = []

            for i in range(len(grid_coordinates)):
                x = grid_coordinates[i,0]
                y = grid_coordinates[i,1]
                if abs(y-0.375) < 1e-12 :
                    if abs(x-1.5) < 1e-12 :
                        coord_source.append([x,y,0])
                        values.append([0,-1])
            values = np.asarray(values)
            values = np.where(values>0,values,0)        
            steady=True
            return coord_source, values,steady;
        dirac_sources.append(dirac_source)

    if (flag_source == 'example_internet2_semi_inside'):

        def dirac_source(t,info_path):

            coord_source = []
            values = []

            for i in range(len(grid_coordinates)):
                x = grid_coordinates[i,0]
                y = grid_coordinates[i,1]
                if abs(y-0.625) < 1e-12 :
                    if abs(x-1.5) < 1e-12 :
                        coord_source.append([x,y,0])
                        values.append([0,-1])
            values = np.asarray(values)
            values = np.where(values>0,values,0)        
            steady=True
            return coord_source, values,steady;
        dirac_sources.append(dirac_source)

    if (flag_source == 'example_book'):

        def dirac_source(t,info_path):

            coord_source = []
            values = []

            for i in range(len(grid_coordinates)):
                x = grid_coordinates[i,0]
                y = grid_coordinates[i,1]
                if abs(y-0.5) < 1e-12 :
                    if abs(x-0.5) < 1e-12 :
                        coord_source.append([x,y,0])
                        values.append([0,-1])
            values = np.asarray(values)
            values = np.where(values>0,values,0)        
            steady=True
            return coord_source, values,steady;
        dirac_sources.append(dirac_source)


    #  prigozhin
    if (flag_source == 'prigozhin'):
        def source(time,coord,flag_mesh):
            fvalue=0.0
            steady=True
            if (flag_mesh == 1 ):
                fvalue=2.0
            return fvalue, steady;
        sources.append(source)


    #  sin
    if (flag_source == 'sin'):
        def source(time,coord,flag_mesh):
            x=coord[0]
            y=coord[1]

            fvalue=0.0
            if (flag_mesh > 0 ):
                fvalue=np.sin(2.0*y)
            steady=True
            return fvalue,steady;
        sources.append(source)
            
    #  rectangles constat
    if (flag_source == 'rect_cnst'):
        def source(time,coord,flag_mesh):
            x=coord[0]
            y=coord[1]
            fvalue=0.0
            steady=True
            if ( (x >= 1.0/8.0) and (x<=3.0/8.0) and 
                 (y >= 1.0/4.0) and (y<=3.0/4.0) ) :
                fvalue=2.0
            return fvalue,steady;
        sources.append(source)

    #  rectangles constat
    if (flag_source == 'cnst_center'):
        def source(time,coord,flag_mesh):
            x=coord[0]
            y=coord[1]

            fvalue=0.0
            steady=True
            if ( (x >= 0.4) and (x<=0.6) and 
                 (y >= 0.4) and (y<=0.6) ) :
                fvalue=1.0
            return fvalue,steady;
        sources.append(source)
    
    #  rectangles cont
    if (flag_source == 'rect_continuous'):
        def source(time,coord,flag_mesh):
            x=coord[0]
            y=coord[1]

            fvalue=0.0
            steady=True
            if ( (x >= 1.0/8.0) and (x<=3.0/8.0) and 
                 (y >= 1.0/4.0) and (y<=3.0/4.0) ) :
                fvalue=(np.sin(4*np.pi*(x-1.0/8.0)) * 
                        np.sin(2*np.pi*(y-1.0/4.0)) ) 
            return fvalue,steady;
        sources.append(source)

    x17 ,  x18 , y17, y18, y19, y20, y21, y22, x19, x20= 0.1, 0.2, 0.1, 0.250, 0.475, 0.625, 0.750, 0.9, 0.45, 0.55
    if (flag_source == '5rcm'):
        def source(time,coord,flag_mesh):
            x=coord[0]
            y=coord[1]

            fvalue=0.0
            steady=True
            if ((x > x17 and x < x18 and y > y17 and y < y18) 
                or (x > x17 and x < x18 and y > y19 and y < y20) 
                or (x > x17 and x < x18 and y > y21 and y < y22) 
                or (x > x19 and x < x20 and y > y19 and y < y20)):
                fvalue=2.0
            return fvalue, steady;
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
            #print(len(input_lines),'source')
            for line in input_lines[0:]:
                coord = [float(w) for w in line.split()[0:3]]
                value = float(line.split()[3])
                #print(coord,value)
                # check for positive values only
                if ( value > 0.0 ):
                    values.append(value)
                    coord_points.append(coord)                    
            input_file.close()
            return coord_points, values,steady;
        dirac_sources.append(dirac_source)


    return sources, dirac_sources;

# build standard tdens_0  given id_tdens0
def example_sink(sinks,dirac_sinks,flag_sink,extra_info,grid_topology=None,grid_coordinates=None): 
    #  constant examples
    if (flag_sink == 'cnst'):
        def sink(time,coord,flag_mesh):
            fvalue=0.0
            steady=True
            if (flag_mesh < 0 ):
                fvalue=2.0
            return fvalue,steady;
        sinks.append(sink)

    #  personal
    if (flag_sink == 'uniform'):
        def sink(time,coord,flag_mesh):
            fvalue=0.0
            steady=True
            fvalue=1.0
            return fvalue, steady;
        sinks.append(sink)

    if (flag_sink== 'plaplacian'):
        # Source term used for testing convergence for
        # plaplacian equation in  
        # Branching Structures Emerging From A
        #Continuous Optimal Transport Model
        # Enrico Facca, Franco Cardin, Mario Putti
        def sink(time,coord,flag_mesh):
            x=coord[0]
            y=coord[1]


            fvalue=0.0
            steady=True
            if (np.sqrt(x**2+y**2)> 2.0/3.0):
                fvalue=1.0/5.0
            return fvalue, steady;
        sinks.append(sink)

    if (flag_sink== 'plaplacian_dirichlet'):
        # Source term used for testing convergence for
        # plaplacian equation in  
        # Branching Structures Emerging From A
        #Continuous Optimal Transport Model
        # Enrico Facca, Franco Cardin, Mario Putti
        def sink(time,coord,flag_mesh):
            fvalue=0.0
            steady=True
            return fvalue, steady;
        sinks.append(sink)

    if (flag_sink == 'one2two'):
        def dirac_sink(t,info_path):
            coord_points=[[0.4,0.9,0.0],[0.6,0.9,0.0]]
            values=[0.5,0.5]
            steady=True
            return coord_points, values,steady;
        dirac_sinks.append(dirac_sink)

    if (flag_sink == 'three_forces_0'):
        # set point at vertices of equilateral triangles
        coord = [
            [np.cos(np.pi/2          ),np.sin(np.pi/2          )],
            [np.cos(np.pi/2+2*np.pi/3),np.sin(np.pi/2+2*np.pi/3)],
            [np.cos(np.pi/2+4*np.pi/3),np.sin(np.pi/2+4*np.pi/3)],
        ]
        
        coord=np.asarray(coord)
        values=-np.where(coord<0,coord,0)
        coord[:,:]=0.3*coord[:,:]
        coord[:,:]=coord+[0.5,0.5]
        coord=np.hstack((coord,np.zeros([3,1])))
        print(coord)
        def dirac_sink(t,info_path):
            steady=True
            return coord, values,steady;
        dirac_sinks.append(dirac_sink)

    if (flag_sink == 'three_forces_1'):
        # set point at vertices of equilateral triangles
        coord = [
            [np.cos(3*np.pi/4),np.sin(3*np.pi/4)],
            [np.cos(5*np.pi/4),np.sin(5*np.pi/4)],
            [np.cos(7*np.pi/4),np.sin(7*np.pi/4)],
        ]        
        
        coord=np.asarray(coord)

        # Esempio 1
        values=[
            [+np.sqrt(2)/4,+np.sqrt(2)/4],
            [-np.sqrt(2)/2,-np.sqrt(2)/2],
            [+np.sqrt(2)/4,+np.sqrt(2)/4],
        ]
        values=np.asarray(values)
        values=-np.where(values<0,values,0)

        coord[:,:]=0.3*coord[:,:]
        coord[:,:]=coord+[0.5,0.5]
        coord=np.hstack((coord,np.zeros([3,1])))
        print(coord)
        def dirac_sink(t,info_path):
            steady=True
            return coord, values,steady;
        dirac_sinks.append(dirac_sink)

    if (flag_sink == 'three_forces_2'):
        # set point at vertices of equilateral triangles
        coord = [
            [np.cos(3*np.pi/4),np.sin(3*np.pi/4)],
            [np.cos(5*np.pi/4),np.sin(5*np.pi/4)],
            [np.cos(7*np.pi/4),np.sin(7*np.pi/4)],
        ]        
        
        coord=np.asarray(coord)

        # Esempio 2
        values=[
            [-0.5,-0.5],
            [+0.0,+1.0],
            [+0.5,-0.5],
        ]
        values=np.asarray(values)
        values=-np.where(values<0,values,0)

        coord[:,:]=0.3*coord[:,:]
        coord[:,:]=coord+[0.5,0.5]
        coord=np.hstack((coord,np.zeros([3,1])))
        print(coord)
        def dirac_sink(t,info_path):
            steady=True
            return coord, values,steady;
        dirac_sinks.append(dirac_sink)

    if (flag_sink == 'three_forces_3'):
        # set point at vertices of equilateral triangles
        coord = [
            [np.cos(3*np.pi/4),np.sin(3*np.pi/4)],
            [np.cos(5*np.pi/4),np.sin(5*np.pi/4)],
            [np.cos(7*np.pi/4),np.sin(7*np.pi/4)],
        ]        
        
        coord=np.asarray(coord)

        # Esempio 3
        values=[
            [-0.5,+0.5*(4-np.sqrt(3))],
            [-0.5*(np.sqrt(3)-1),-0.5*(3-np.sqrt(3))],
            [+0.5*np.sqrt(3),-0.5],
        ]
        values=np.asarray(values)
        values=-np.where(values<0,values,0)

        coord[:,:]=0.3*coord[:,:]
        coord[:,:]=coord+[0.5,0.5]
        coord=np.hstack((coord,np.zeros([3,1])))
        print(coord)
        def dirac_sink(t,info_path):
            steady=True
            return coord, values,steady;
        dirac_sinks.append(dirac_sink)

    if (flag_sink == 'three_forces_4'):
        # set point at vertices of equilateral triangles
        coord = [
            [np.cos(3*np.pi/4),np.sin(3*np.pi/4)],
            [np.cos(5*np.pi/4),np.sin(5*np.pi/4)],
            [np.cos(7*np.pi/4),np.sin(7*np.pi/4)],
        ]        
        
        coord=np.asarray(coord)

        # Esempio 4
        values=[
            [+1*np.sin(np.pi/12),-1*np.cos(np.pi/12)],
            [-1*np.sin(np.pi/4), +1*np.cos(np.pi/4) ],
            [0.5*np.sqrt(2)*(np.sqrt(3)-1)*np.cos(np.pi/6),
             0.5*np.sqrt(2)*(np.sqrt(3)-1)*np.sin(np.pi/6)],
        ]
        values=np.asarray(values)
        values=-np.where(values<0,values,0)

        coord[:,:]=0.3*coord[:,:]
        coord[:,:]=coord+[0.5,0.5]
        coord=np.hstack((coord,np.zeros([3,1])))
        print(coord)
        def dirac_sink(t,info_path):
            steady=True
            return coord, values,steady;
        dirac_sinks.append(dirac_sink)

    if (flag_sink == 'three_forces_5'):
        # set point at vertices of equilateral triangles
        coord = [
            [np.cos(3*np.pi/4),np.sin(3*np.pi/4)],
            [np.cos(5*np.pi/4),np.sin(5*np.pi/4)],
            [np.cos(7*np.pi/4),np.sin(7*np.pi/4)],
        ]        
        
        coord=np.asarray(coord)

        # Esempio 5
        values=[
            [-0.5*np.sqrt(3),0.5*np.sqrt(3)*(2-np.sqrt(3))],
            [0.5*(np.sqrt(3)-1),0.5*np.sqrt(3)*(np.sqrt(3)-1)],
            [0.5,-0.5*np.sqrt(3)],
        ]
        values=np.asarray(values)
        values=-np.where(values<0,values,0)

        coord[:,:]=0.3*coord[:,:]
        coord[:,:]=coord+[0.5,0.5]
        coord=np.hstack((coord,np.zeros([3,1])))
        print(coord)
        def dirac_sink(t,info_path):
            steady=True
            return coord, values,steady;
        dirac_sinks.append(dirac_sink)

    if (flag_sink == 'three_forces_6'):
        # set point at vertices of equilateral triangles
        coord = [
            [np.cos(3*np.pi/4),np.sin(3*np.pi/4)],
            [np.cos(5*np.pi/4),np.sin(5*np.pi/4)],
            [np.cos(7*np.pi/4),np.sin(7*np.pi/4)],
        ]        
        
        coord=np.asarray(coord)

        # Esempio 6
        values=[
            [0.5*np.sqrt(3),+0.5-1/np.sqrt(3)],
            [-0.5*(np.sqrt(3)+1),-(3+np.sqrt(3))/6],
            [0.5,0.5*np.sqrt(3)],
        ]
        values=np.asarray(values)
        values=-np.where(values<0,values,0)

        coord[:,:]=0.3*coord[:,:]
        coord[:,:]=coord+[0.5,0.5]
        coord=np.hstack((coord,np.zeros([3,1])))
        print(coord)
        def dirac_sink(t,info_path):
            steady=True
            return coord, values,steady;
        dirac_sinks.append(dirac_sink)

    x11 , x12 , y11, y12 = 0.750, 0.9, 0.475, 0.625
    if (flag_sink == '5rcm'):
        def sink(time,coord,flag_mesh):
            x=coord[0]
            y=coord[1]

            fvalue=0.0
            steady=True
            if (x > x11 and x < x12 and y > y11 and y < y12):
                fvalue=2.0
            return fvalue, steady;
        sinks.append(sink)

    if (flag_sink == 'example3'):

        def dirac_sink(t,info_path):

            coord_source = []
            values = []

            print('numero di nodi', len(grid_coordinates))
            for i in range(len(grid_coordinates)):
                x = grid_coordinates[i,0]
                y = grid_coordinates[i,1]
                if abs(y-0.4) < 1e-12 :
                    if abs(x-1.2) < 1e-12 :
                        print(x)
                        print(y)
                        coord_source.append([x,y,0])
                        values.append([0,-1])
            values = np.asarray(values)
            values = -np.where(values<0,values,0)        
            steady=True
            return coord_source, values,steady;
        dirac_sinks.append(dirac_sink)

    if (flag_sink == 'example3_semi_inside'):

        def dirac_sink(t,info_path):

            coord_source = []
            values = []

            print('numero di nodi', len(grid_coordinates))
            for i in range(len(grid_coordinates)):
                x = grid_coordinates[i,0]
                y = grid_coordinates[i,1]
                if abs(y-0.4) < 1e-12 :
                    if abs(x-1.2) < 1e-12 :
                        print(x)
                        print(y)
                        coord_source.append([x,y,0])
                        values.append([0,-1])
            values = np.asarray(values)
            values = -np.where(values<0,values,0)        
            steady=True
            return coord_source, values,steady;
        dirac_sinks.append(dirac_sink)

    if (flag_sink == 'example3_inside'):

        def dirac_sink(t,info_path):

            coord_source = []
            values = []

            print('numero di nodi', len(grid_coordinates))
            for i in range(len(grid_coordinates)):
                x = grid_coordinates[i,0]
                y = grid_coordinates[i,1]
                if abs(y-1.2) < 1e-12 :
                    if abs(x-1.5) < 1e-12 :
                        print(x)
                        print(y)
                        coord_source.append([x,y,0])
                        values.append([0,-1])
            values = np.asarray(values)
            values = -np.where(values<0,values,0)        
            steady=True
            return coord_source, values,steady;
        dirac_sinks.append(dirac_sink)

    if (flag_sink == 'beam'):

        def dirac_sink(t,info_path):

            coord_source = []
            values = []

            print('numero di nodi', len(grid_coordinates))
            for i in range(len(grid_coordinates)):
                x = grid_coordinates[i,0]
                y = grid_coordinates[i,1]
                if abs(y-1.2) < 1e-12 :
                    if abs(x-1.5) < 1e-12 :
                        print(x)
                        print(y)
                        coord_source.append([x,y,0])
                        values.append([0,-1])
            values = np.asarray(values)
            values = -np.where(values<0,values,0)        
            steady=True
            return coord_source, values,steady;
        dirac_sinks.append(dirac_sink)

    if (flag_sink == 'beam_shear'):

        def dirac_sink(t,info_path):

            coord_source = []
            values = []

            print('numero di nodi', len(grid_coordinates))
            for i in range(len(grid_coordinates)):
                x = grid_coordinates[i,0]
                y = grid_coordinates[i,1]
                if abs(x-2) < 1e-12 :
                    if (y >  - 1e-12) and (y < 0.3 + 1e-12) :
                        print(i)
                        coord_source.append([x,y,0])
                        values.append([0,-0.01])
            values = np.asarray(values)
            values = -np.where(values<0,values,0)        
            steady=True
            return coord_source, values,steady;
        dirac_sinks.append(dirac_sink)

    if (flag_sink == 'beam_shear_inside'):

        def dirac_sink(t,info_path):

            coord_source = []
            values = []

            print('numero di nodi', len(grid_coordinates))
            for i in range(len(grid_coordinates)):
                x = grid_coordinates[i,0]
                y = grid_coordinates[i,1]
                if abs(x-2.5) < 1e-12 :
                    if (y > 0.35 - 1e-12) and (y < 0.65 + 1e-12) :
                        print(i)
                        coord_source.append([x,y,0])
                        values.append([0,-1])
            values = np.asarray(values)
            values = -np.where(values<0,values,0)        
            steady=True
            return coord_source, values,steady;
        dirac_sinks.append(dirac_sink)

    if (flag_sink == 'example2_concentrated_load'):

        def dirac_sink(t,info_path):

            coord_source = []
            values = []

            print('numero di nodi', len(grid_coordinates))
            for i in range(len(grid_coordinates)):
                x = grid_coordinates[i,0]
                y = grid_coordinates[i,1]
                if abs(x-2.0) < 1e-12 :
                    if abs(y-0.15) < 1e-12 :
                        print(x)
                        print(y)
                        print(i)
                        coord_source.append([x,y,0])
                        values.append([0,-1])
            values = np.asarray(values)
            values = -np.where(values<0,values,0)        
            steady=True
            return coord_source, values,steady;
        dirac_sinks.append(dirac_sink)

    if (flag_sink == 'example2_concentrated_load_inside'):

        def dirac_sink(t,info_path):

            coord_source = []
            values = []

            print('numero di nodi', len(grid_coordinates))
            for i in range(len(grid_coordinates)):
                x = grid_coordinates[i,0]
                y = grid_coordinates[i,1]
                if abs(x-2.5) < 1e-12 :
                    if abs(y-0.5) < 1e-12 :
                        print(x)
                        print(y)
                        print(i)
                        coord_source.append([x,y,0])
                        values.append([0,-1])
            values = np.asarray(values)
            values = -np.where(values<0,values,0)        
            steady=True
            return coord_source, values,steady;
        dirac_sinks.append(dirac_sink)

    if (flag_sink == 'example1_concentrated_load_compression'):

        def dirac_sink(t,info_path):

            coord_source = []
            values = []

            print('numero di nodi', len(grid_coordinates))
            for i in range(len(grid_coordinates)):
                x = grid_coordinates[i,0]
                y = grid_coordinates[i,1]
                if abs(y-2.0) < 1e-12 :
                    if abs(x-0.15) < 1e-12 :
                        print(i)
                        coord_source.append([x,y,0])
                        values.append([0,-1])
            values = np.asarray(values)
            values = -np.where(values<0,values,0)        
            steady=True
            return coord_source, values,steady;
        dirac_sinks.append(dirac_sink)

    if (flag_sink == 'example1_concentrated_load_traction'):

        def dirac_sink(t,info_path):

            coord_source = []
            values = []

            print('numero di nodi', len(grid_coordinates))
            for i in range(len(grid_coordinates)):
                x = grid_coordinates[i,0]
                y = grid_coordinates[i,1]
                if abs(y-2.0) < 1e-12 :
                    if abs(x-0.15) < 1e-12 :
                        print(i)
                        coord_source.append([x,y,0])
                        values.append([0,+1])
            values = np.asarray(values)
            values = -np.where(values<0,values,0)        
            steady=True
            return coord_source, values,steady;
        dirac_sinks.append(dirac_sink)

    if (flag_sink == 'example1_concentrated_load_compression_inside'):

        def dirac_sink(t,info_path):

            coord_source = []
            values = []

            print('numero di nodi', len(grid_coordinates))
            for i in range(len(grid_coordinates)):
                x = grid_coordinates[i,0]
                y = grid_coordinates[i,1]
                if abs(y-2.5) < 1e-12 :
                    if abs(x-0.5) < 1e-12 :
                        print(i)
                        coord_source.append([x,y,0])
                        values.append([0,-1])
            values = np.asarray(values)
            values = -np.where(values<0,values,0)        
            steady=True
            return coord_source, values,steady;
        dirac_sinks.append(dirac_sink)

    if (flag_sink == 'example1_concentrated_load_traction_inside'):

        def dirac_sink(t,info_path):

            coord_source = []
            values = []

            print('numero di nodi', len(grid_coordinates))
            for i in range(len(grid_coordinates)):
                x = grid_coordinates[i,0]
                y = grid_coordinates[i,1]
                if abs(y-2.5) < 1e-12 :
                    if abs(x-0.5) < 1e-12 :
                        print(i)
                        coord_source.append([x,y,0])
                        values.append([0,+1])
            values = np.asarray(values)
            values = -np.where(values<0,values,0)        
            steady=True
            return coord_source, values,steady;
        dirac_sinks.append(dirac_sink)

    if (flag_sink == 'column_traction'):

        def dirac_sink(t,info_path):

            coord_source = []
            values = []

            print('numero di nodi', len(grid_coordinates))
            for i in range(len(grid_coordinates)):
                x = grid_coordinates[i,0]
                y = grid_coordinates[i,1]
                if abs(y-2.0) < 1e-12 :
                    if (x >  - 1e-12) and (x < 0.3 + 1e-12) :
                        print(i)
                        coord_source.append([x,y,0])
                        values.append([0,1])
            values = np.asarray(values)
            values = -np.where(values<0,values,0)        
            steady=True
            return coord_source, values,steady;
        dirac_sinks.append(dirac_sink)

    if (flag_sink == 'column_traction_inside'):

        def dirac_sink(t,info_path):

            coord_source = []
            values = []

            print('numero di nodi', len(grid_coordinates))
            for i in range(len(grid_coordinates)):
                x = grid_coordinates[i,0]
                y = grid_coordinates[i,1]
                if abs(y-2.5) < 1e-12 :
                    if (x > 0.35 - 1e-12) and (x < 0.65 + 1e-12) :
                        print(i)
                        coord_source.append([x,y,0])
                        values.append([0,1])
            values = np.asarray(values)
            values = -np.where(values<0,values,0)        
            steady=True
            return coord_source, values,steady;
        dirac_sinks.append(dirac_sink)

    if (flag_sink == 'column_compression'):

        def dirac_sink(t,info_path):

            coord_source = []
            values = []

            print('numero di nodi', len(grid_coordinates))
            for i in range(len(grid_coordinates)):
                x = grid_coordinates[i,0]
                y = grid_coordinates[i,1]
                if abs(y-2.0) < 1e-12 :
                    if (x >  - 1e-12) and (x < 0.3 + 1e-12) :
                        print(i)
                        coord_source.append([x,y,0])
                        values.append([0,-1])
            values = np.asarray(values)
            values = -np.where(values<0,values,0)        
            steady=True
            return coord_source, values,steady;
        dirac_sinks.append(dirac_sink)

    if (flag_sink == 'column_compression_inside'):

        def dirac_sink(t,info_path):

            coord_source = []
            values = []

            print('numero di nodi', len(grid_coordinates))
            for i in range(len(grid_coordinates)):
                x = grid_coordinates[i,0]
                y = grid_coordinates[i,1]
                if abs(y-2.5) < 1e-12 :
                    if (x > 0.35 - 1e-12) and (x < 0.65 + 1e-12) :
                        print(i)
                        coord_source.append([x,y,0])
                        values.append([0,-1])
            values = np.asarray(values)
            values = -np.where(values<0,values,0)        
            steady=True
            return coord_source, values,steady;
        dirac_sinks.append(dirac_sink)

    if (flag_sink == 'example_internet1'):

        def dirac_sink(t,info_path):

            coord_source = []
            values = []

            print('numero di nodi', len(grid_coordinates))
            for i in range(len(grid_coordinates)):
                x = grid_coordinates[i,0]
                y = grid_coordinates[i,1]
                if abs(y-0.0) < 1e-12 :
                    if abs(x-1.5) < 1e-12 :
                        print(i)
                        coord_source.append([x,y,0])
                        values.append([0,-1])
            values = np.asarray(values)
            values = -np.where(values<0,values,0)        
            steady=True
            return coord_source, values,steady;
        dirac_sinks.append(dirac_sink)

    if (flag_sink == 'example_internet1_semi_inside'):

        def dirac_sink(t,info_path):

            coord_source = []
            values = []

            print('numero di nodi', len(grid_coordinates))
            for i in range(len(grid_coordinates)):
                x = grid_coordinates[i,0]
                y = grid_coordinates[i,1]
                if abs(y-0.0) < 1e-12 :
                    if abs(x-1.5) < 1e-12 :
                        print(i)
                        coord_source.append([x,y,0])
                        values.append([0,-1])
            values = np.asarray(values)
            values = -np.where(values<0,values,0)        
            steady=True
            return coord_source, values,steady;
        dirac_sinks.append(dirac_sink)

    if (flag_sink == 'example_internet2'):

        def dirac_sink(t,info_path):

            coord_source = []
            values = []

            print('numero di nodi', len(grid_coordinates))
            for i in range(len(grid_coordinates)):
                x = grid_coordinates[i,0]
                y = grid_coordinates[i,1]
                if abs(y-0.375) < 1e-12 :
                    if abs(x-1.5) < 1e-12 :
                        print(i)
                        coord_source.append([x,y,0])
                        values.append([0,-1])
            values = np.asarray(values)
            values = -np.where(values<0,values,0)        
            steady=True
            return coord_source, values,steady;
        dirac_sinks.append(dirac_sink)

    if (flag_sink == 'example_internet2_semi_inside'):

        def dirac_sink(t,info_path):

            coord_source = []
            values = []

            print('numero di nodi', len(grid_coordinates))
            for i in range(len(grid_coordinates)):
                x = grid_coordinates[i,0]
                y = grid_coordinates[i,1]
                if abs(y-0.625) < 1e-12 :
                    if abs(x-1.5) < 1e-12 :
                        print(i)
                        coord_source.append([x,y,0])
                        values.append([0,-1])
            values = np.asarray(values)
            values = -np.where(values<0,values,0)        
            steady=True
            return coord_source, values,steady;
        dirac_sinks.append(dirac_sink)

    if (flag_sink == 'example_book'):

        def dirac_sink(t,info_path):

            coord_source = []
            values = []

            print('numero di nodi', len(grid_coordinates))
            for i in range(len(grid_coordinates)):
                x = grid_coordinates[i,0]
                y = grid_coordinates[i,1]
                if abs(y-0.5) < 1e-12 :
                    if abs(x-0.5) < 1e-12 :
                        print(i)
                        coord_source.append([x,y,0])
                        values.append([0,-1])
            values = np.asarray(values)
            values = -np.where(values<0,values,0)        
            steady=True
            return coord_source, values,steady;
        dirac_sinks.append(dirac_sink)

    #  A parallel fast sweeping method for the Eikonal equation
    #  Miles Detrixhea, Frederic Giboua,b, Chohong Min
    if (flag_sink == 'example1fmt'):
        def sink(time,coord,flag_mesh):
            fvalue=1.0
            steady=True
            return fvalue, steady;
        sinks.append(sink)

    # An accurate discontinuous Galerkin method for solving point-source
    # Eikonal equation in 2-D heterogeneous anisotropic media
    # Bouteiller,1 M. Benjemaa,2 L. Metivier1,3 and J. Virieux1
    if (flag_sink == 'bouteiller18_tc1'):
        def sink(t,coord,info_path):
            fvalue=1.0
            steady=True
            return fvalue, steady;
        sinks.append(sink)


    #  personal
    if (flag_sink == 'corner'):
        def sink(time,coord,flag_mesh):
            x=coord[0]
            y=coord[1]
            
            fvalue=0.0
            steady=True
            if ( x**2+(y-1.0)**2>0.005 ):
                fvalue=1.0
            return fvalue, steady;
        sinks.append(sink)

    #  rectangles constat
    if (flag_sink == 'cnst_center'):
        def source(time,coord,flag_mesh):
            x=coord[0]
            y=coord[1]

            fvalue=0.0
            steady=True
            if ( (x < 0.4) and (x>0.6) and 
                 (y < 0.4) and (y>0.6) ) :
                fvalue=1.0
            return fvalue,steady;
        sources.append(source)



    #  prigozhin
    if (flag_sink == 'prigozhin'):
        def sink(time,coord,flag_mesh):
            fvalue=0.0
            steady=True
            if (flag_mesh == -1 ):
                fvalue=2.0
            return fvalue, steady;
        sinks.append(sink)


    #  sin
    if (flag_sink == 'sin'):
        def sink(time,coord,flag_mesh):
            x=coord[0]
            y=coord[1]

            fvalue=0.0
            if (flag_mesh < 0 ):
                fvalue=np.sin(2.0*x)
            steady=True
            return fvalue,steady;
        sinks.append(sink)


    # classical example
    if (flag_sink == 'rect_cnst'):
        def sink(time,coord,flag_mesh):
            x=coord[0]
            y=coord[1]

            fvalue=0.0
            steady=True
            if ( (x >= 5.0/8.0) and 
                 (x <= 7.0/8.0) and 
                 (y >= 1.0/4.0) and 
                 (y<=3.0/4.0)     ) :
                fvalue=2.0
            return fvalue, steady;
        sinks.append(sink)
    
    #  rectangles continouos
    if (flag_sink == 'rect_continuous'):
        def sink(time,coord,flag_mesh):
            x=coord[0]
            y=coord[1]

            fvalue=0.0
            steady=True
            if ( (x >= 5.0/8.0) and (x<=7.0/8.0) and 
                 (y >= 1.0/4.0) and (y<=3.0/4.0) ) :
                fvalue=(np.sin(4*np.pi*(x-5.0/8.0)) * 
                        np.sin(2*np.pi*(y-1.0/4.0)) ) 
            return fvalue,steady;
        sinks.append(sink)

     
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
                coord_dirac = [float(w) for w in line.split()[0:3]]
                branch_value = float(line.split()[3])
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
                coord = [float(w) for w in line.split()[0:3]]
                value = float(line.split()[3])
                #print(coord,value)
                # check for positive values only
                if ( value > 0.0 ):
                    values.append(value)
                    coord_points.append(coord)                    
            input_file.close()
            return coord_points, values,steady;
        dirac_sinks.append(dirac_sink)


    return sinks, dirac_sinks;

# build 
def example_source_3d(sources,dirac_sources,flag_source,extra_info):
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
        
        print ('Sink/Source - 1/19:',np.abs(c2/c1+1.0/19.0))
        
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
def example_sink_3d(sinks,dirac_sinks,flag_sink,extra_info):    
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
        print ('Sink/Source - 1/19:',np.abs(c2/c1+1.0/19.0))
        

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



def known_optpot(flag_grid,flag_source,flag_sink):
    if ( flag_grid=='example1fmt' and 
         flag_source=='example1fmt' and 
         flag_sink=='example1fmt'):
        def optpot(coord,flag):
            x=coord[0]
            y=coord[1]
            z=coord[2]
            if (x>0 and y>0):
                f=x**2+y**2+3*x*y
            else:
                f=x**2+y**2
            return f;
    
    if ( flag_grid=='bouteiller18_tc1' and 
         flag_source=='bouteiller18_tc1' and 
         flag_sink=='bouteiller18_tc1'):
        def optpot(coord,flag):
            x=coord[0]
            y=coord[1]
            z=coord[2]
            #
            Szero=2
            gzero=(0.0,0.5)
            gzeronorm=np.sqrt(gzero[0]**2+gzero[1]**2)
            xzero=(2.0,2.0)
            S=1/(1/Szero + np.scal(gzero,coord[0:1]-xzero))
            
            f=1/Gzeronorm * np.arccosh(
                1+0.5*S*Szero*Gzeronrm**2*((x-xzero[0])**2+(y-xzer0[1])**2))
            return f;
    
    return optpot;
