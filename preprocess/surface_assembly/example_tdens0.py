#!/usr/bin/env python
import numpy as np
import jw_meshtools as mt
import common
import branch

ntype=3
nparameter=5

#define basic tdens0 functions
def default_tdens0_func(x,y,flag_func,a):
    #constant
    if (flag_func == 'cnst'):
        tdens0_value=a[0]
    #parabola
    if (flag_func == 'par'):    
        tdens0_value=a[0]+a[1]*(x-a[2])**2+a[3]*(y-a[4])**2
    #sin
    if (flag_func == 'sin'):    
        tdens0_value=a[0]+a[1]*np.sin(x*a[2])+a[3]*np.sin(y*a[4])
    #print tdens0_value
    return tdens0_value;


def build_tdens0(tdens0_functions,triang,bar_tria,flags): 
    tdens0_tria=np.zeros([len(bar_tria)])
    for itria in range(len(bar_tria)):
        for ifunc in range(len(tdens0_functions)):
            tvalue=tdens0_functions[ifunc](
                bar_tria[itria][0],
                bar_tria[itria][1],
                flags[itria])
            tdens0_tria[itria]+=tvalue
    return tdens0_tria


# build standard tdens_0  given id_tdens0
def example_tdens0(flag_tdens0,length,extra_info):
    #  personal
    if (flag_tdens0 == 'one'):
        # cnst one on whole domain 
        def id_region(x,y,flag):
            coefficients=[1.0, 0.0 , 0.0, 0.0 , 0.0 ]
            return default_tdens0_func(x,y,'cnst',coefficients)

    # parabola sen
    if (flag_tdens0 == 'par0505'):
        def id_region(x,y,flag):
            # whole dmain 
            func='par'
            coefficients=[0.01, 4.0 , 0.5, 4.0 , 0.5 ]
            return default_tdens0_func(x,y,'par',coefficients)
    # parabola
    if (flag_tdens0 == 'sin'):
        def id_region(x,y,flag):
            # whole dmain 
            func='sin'
            coefficients=[1.0, 2.0 , 0.5, 2.0 , 0.5 ]
            return default_tdens0_func(x,y,'sin',coefficients)
    # sin
    if (flag_tdens0 == 'par_eastwest'):
        def id_region(x,y,flag):
            if ( x <= 0.5 ):
                # west
                coefficients=[0.01, 2.0 , 0.5, .0 , 0.0 ]
                fvalue=default_tdens0_func(x,y,'par',coefficients)
            if ( x > 0.5 ):
                # east
                coefficients=[1.0, 2.0 , 0.5, 2.0 , 0.5 ]
                fvalue=default_tdens0_func(x,y,'par',coefficients)
            return fvalue;
    
    # y_branch
    if (flag_tdens0 == 'ybranch'):
        
        info_path=common.remove_comments(extra_info,'!')
        coord_points=[]
        values=[]
        input_file=open(info_path,'r')
        input_lines = input_file.readlines()
        for line in input_lines[0:3]:
            coord_dirac = [float(w) for w in line.split()[0:2]]
            branch_value = float(line.split()[2])
            values.append(branch_value)
            coord_points.append(coord_dirac)                    
            
            
        alpha=float(common.remove_comments(input_lines[4],'!'))
            
        # Build the optimal solution
        # add branching point (if exists)
        
        branch_points,branch_vertices = branch.y_branch(
            coord_points,values,alpha)
        input_file.close()

        def id_region(x,y,flag):
            point=[x,y]
            fvalue = 1.0e-10
            for segments in branch_vertices :
                distance=mt.distance_segment_point(
                    branch_points[segments[0]],
                    branch_points[segments[1]],
                    point)
                if ( distance < length / 1.5 ):
                    fvalue = 1.0e10
                    
            return fvalue;
    
    return id_region;

# build standard tdens_0  given id_tdens0
def build_by_seed(flag_tdens0, coord, triang, bar_tria, seed):
    #  personal
    def id_region(x,y):
        if ( x < 0.5 ) and ( y < 5 ):
            # west
            func='par'
            coefficients=[1.0, 2.0*seed , 0.5, 2.0*seed , 0.0 ]
        if ( x > 0.5 ) and ( y < 5 ):
            # east
            func='par'
            coefficients=[1.0, 2.0 , 0.5, 2.0 , 0.5 ]
        return func,coefficients;

    tdens0_tria=build_tdens0(id_region,coord,triang,bar_tria);
    return tdens0_tria;

