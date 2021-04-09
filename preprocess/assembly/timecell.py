#!/usr/bin/env python
import numpy as np
import common as common
from random import seed
from random import random
#
# Programs handling creation of time data defined on
# cell grid. Contains
#
# build: 
# inputs  : list of functions, coordinate of baricenter, time 
# outputs : values of the function at baricenters, logical for steady state
#
# write: write time data in proper format
#
# example : create example of functions give flags
# 

#
# evaluate cells functions at given time 
# it supports list of functions that are summed.
#
def build(cell_values,steady,
          functions,time,bar_cell,flags):
    steady=True
    for icell in range(len(cell_values)):
        for ifun in range(len(functions)):
            value,steady_fun = functions[ifun](
                time,
                bar_cell[icell,:],
                flags[icell])

            cell_values[icell,:]+=value
            steady = steady and steady_fun
    return cell_values,steady;

def eval(functions,time,bar_cell,flags):
    # get dimensions
    value,steady_fun = functions[ifun](
        time,
        bar_cell[0,:],flags[0])
    nout=len(value)
    ncell=len(topol)
    if (nout == 1) :
        cell_values=np.zeros(ncell)
    else:
        cell_values=np.zeros([nout,ncell])
        
    steady=True
    #
    # eval function at centroid
    # 
    if (len(functions)>1):
        for icell in range(len(cell_values)):
            for ifun in range(len(functions)):
                value,steady_fun = functions[ifun](
                    time,
                    bar_cell[icell,:],
                    flags[icell])
                cell_values[icell,:]+=value
                steady = steady and steady_fun
    else:
        for icell in range(len(cell_values)):
            value,steady_fun = functions(
                time,
                bar_cell[icell,:],
                flags[icell])
            cell_values[icell,:]+=value
            steady = steady and steady_fun
                
    return cell_values,steady;


#
# Write time data to file
#
def write2file(file_out,time,head,data,steady=False):
    if ( head ) :
        file_out.write(
            str(data.shape[1]) + " " + 
            str(data.shape[0]) + " !  dim ndata" 
            +"\n")
    else:
        nrm_data=np.zeros(len(data))
        for i in range(len(data)):
            nrm_data[i]=np.linalg.norm(data[:][i])
        ndata=(abs(nrm_data) != 0.0).sum()
        #print(ndata)
        file_out.write("time    "+ str(time)+" \n")
        file_out.write(str(ndata)+" \n")
        for i in range(len(data)):
            if ( np.sum(np.abs(data[i][:])) != 0.0):
                file_out.write(str(i+1)+" " + 
                               " ".join(map(str,data[i][:])) +"\n")
        if( steady ):
            file_out.write("time " + str(1.e30)+"\n")


#
# Defition of timedata example
# Given flag and extra_info return a function with
#
# function inputs  are time, coord,falg_mesh 
# function outputs are function_value and steady
#
# Use to define also varaible that do not vary in time
# like tdens0, optdens, etc. We impose steaddy for all time.
# We use this format to unify formats of input data
# 

def example(flag,extra_info=None,flag_mesh=None,length=None):
    # if is const contained in flag_kappa
    try :
        cnst=float(flag)
        def funct(time,coord,flag_mesh):
            value=cnst
            steady=True
            return cnst, steady;
    #
    # user defined example
    # 
    except ValueError:
        #
        #  parabola centered at origin
        #
        if (flag == 'par'):
            def funct(time,coord,flag_mesh):
                x=coord[0]
                y=coord[1]
                value=1.0+x**2+y**2
                steady=True
                return value,steady;
        #
        #  parabola centered at origin
        #
        if (flag == 'parabola_1'):
            def funct(time,coord,flag_mesh):
                x=coord[0]
                y=coord[1]
                value=1.0+x**2+y**2
                steady=True
                return value,steady;

        if (flag == 'parabola_2'):
            def funct(time,coord,flag_mesh):
                x=coord[0]
                y=coord[1]
                x_0 = 1
                y_0 = 1
                a = 1
                value = a*((x-x_0)**2 + (y-y_0)**2)
                steady=True
                return value,steady;

        if (flag == 'parabola_3'):
            def funct(time,coord,flag_mesh):
                x=coord[0]
                y=coord[1]
                #print('x',x)
                #print('y',y)
                x_0 = 1
                y_0 = 1
                a = 1
                value = -a*((x-x_0)**3 + (y-y_0)**3)
                #print('funct',value)
                steady=True
                return value,steady;

        if (flag == 'parabola_4'):
            def funct(time,coord,flag_mesh):
                x=coord[0]
                y=coord[1]
                x_0 = 1
                y_0 = 1
                a = 1
                value = a*((x-x_0)**4 + (y-y_0)**4)
            
                steady=True
                return value,steady;

        if (flag == 'parabola_5'):
            def funct(time,coord,flag_mesh):
                x=coord[0]
                y=coord[1]
                x_0 = 1
                y_0 = 1
                a = 1
                value = -a*((x-x_0)**5 + (y-y_0)**5)
            
                steady=True
                return value,steady;

        if (flag == 'parabola_6'):
            def funct(time,coord,flag_mesh):
                x=coord[0]
                y=coord[1]
                x_0 = 1
                y_0 = 1
                a = 1
                value = a*((x-x_0)**6 + (y-y_0)**6)
            
                steady=True
                return value,steady;

        if (flag == 'parabola_7'):
            def funct(time,coord,flag_mesh):
                x=coord[0]
                y=coord[1]
                x_0 = 1
                y_0 = 1
                a = 1
                value = -a*((x-x_0)**7 + (y-y_0)**7)
            
                steady=True
                return value,steady;

        if (flag == 'parabola_8'):
            def funct(time,coord,flag_mesh):
                x=coord[0]
                y=coord[1]
                x_0 = 1
                y_0 = 1
                a = 1
                value = a*((x-x_0)**8 + (y-y_0)**8)
            
                steady=True
                return value,steady;                

        if (flag == 'parabola_9'):
            def funct(time,coord,flag_mesh):
                x=coord[0]
                y=coord[1]
                x_0 = 1
                y_0 = 1
                a = 1
                value = -a*((x-x_0)**9 + (y-y_0)**9)
            
                steady=True
                return value,steady;

        if (flag == 'parabola_10'):
            def funct(time,coord,flag_mesh):
                x=coord[0]
                y=coord[1]
                x_0 = 1
                y_0 = 1
                a = 1
                value = a*((x-x_0)**10 + (y-y_0)**10)
            
                steady=True
                return value,steady;

        if (flag == 'parabola_11'):
            def funct(time,coord,flag_mesh):
                x=coord[0]
                y=coord[1]
                x_0 = 1
                y_0 = 1
                a = 1
                value = a*((x-x_0)**2 - (y-y_0)**5)
            
                steady=True
                return value,steady;                

        if ( flag == 'delta_1' ):
            def funct(time, coord, flag_mesh):
                x=coord[0]
                y=coord[1]
                a = 1

                if ((x>= 0.5 and x <= 0.6) and (y>= 0.5 and y<=0.6)):
                    value = 1.0
                else:
                    value = 0.00001
                steady= True
                return value,steady;

        if ( flag == 'delta_2' ):
            def funct(time, coord, flag_mesh):
                x=coord[0]
                y=coord[1]
                a = 1

                if ((x>= 0.1 and x <= 0.2) and (y>= 0.1 and y<=0.2)):
                    value = 1.0
                else:
                    value = 0.000000001
                steady= True
                return value,steady;

        if ( flag == 'delta_3' ):
            def funct(time, coord, flag_mesh):
                x=coord[0]
                y=coord[1]
                a = 1

                if ((x>= 0.8 and x <= 0.9) and (y>= 0.8 and y<=0.9)):
                    value = 1.0
                else:
                    value = 0.000000001
                steady= True
                return value,steady; 

        if ( flag == 'multiple_deltas' ):
            def funct(time, coord, flag_mesh):
                x=coord[0]
                y=coord[1]
                a = 1

                if ((x>= 0.1 and x <= 0.2) and (y>= 0.8 and y<=0.9)):
                    value = 1.0
                elif ((x>= 0.5 and x <= 0.6) and (y>= 0.5 and y<=0.6)):
                    value = 1.0                
                elif ((x>= 0.8 and x <= 0.9) and (y>= 0.8 and y<=0.9)):
                    value = 1.0
                elif ((x>= 0.01 and x <= 0.1) and (y>= 0.4 and y<=0.5)):
                    value = 1.0                 
                elif ((x>= 0.01 and x <= 0.1) and (y>= 0.1 and y<=0.2)):
                    value = 1.0                 
                elif ((x>= 0.5 and x <= 0.6) and (y>= 0.8 and y<=0.9)):
                    value = 1.0                 
                elif ((x>= 0.5 and x <= 0.6) and (y>= 0.1 and y<=0.2)):
                    value = 1.0                 
                elif ((x>= 0.90 and x <= 0.99) and (y>= 0.5 and y<=0.6)):
                    value = 1.0                 
                elif ((x>= 0.8 and x <= 0.9) and (y>= 0.1 and y<=0.2)):
                    value = 1.0                 

                else:
                    value = 0.0000000001
                steady= True
                return value,steady; 



        if ( flag == 'delta_4' ):
            def funct(time, coord, flag_mesh):
                x=coord[0]
                y=coord[1]
                a = 1

                if ((x>= 0.1 and x <= 0.2) and (y>= 0.8 and y<=0.9)):
                    value = 1.0
                else:
                    value = 0.000000001
                steady= True
                return value,steady;                              


        if ( flag == 'delta_5' ):
            def funct(time, coord, flag_mesh):
                x=coord[0]
                y=coord[1]
                a = 1

                if ((x>= 0.10 and x <= 0.12) and (y>= 0.4 and y<=0.5)):
                    value = 1.0
                else:
                    value = 0.000000001
                steady= True
                return value,steady;                              


        if ( flag == 'delta_6' ):
            def funct(time, coord, flag_mesh):
                x=coord[0]
                y=coord[1]
                a = 1

                if ((x>= 0.5 and x <= 0.6) and (y>= 0.8 and y<=0.9)):
                    value = 1.0
                else:
                    value = 0.000000001
                steady= True
                return value,steady;                              

        if ( flag == 'delta_7' ):
            def funct(time, coord, flag_mesh):
                x=coord[0]
                y=coord[1]
                a = 1

                if ((x>= 0.5 and x <= 0.6) and (y>= 0.1 and y<=0.2)):
                    value = 1.0
                else:
                    value = 0.000000001
                steady= True
                return value,steady;                              



        if ( flag == 'delta_8' ):
            def funct(time, coord, flag_mesh):
                x=coord[0]
                y=coord[1]
                a = 1

                if ((x>= 0.9 and x <= 0.99) and (y>= 0.4 and y<=0.5)):
                    value = 1.0
                else:
                    value = 0.000000001
                steady= True
                return value,steady;                              


        if ( flag == 'delta_9' ):
            def funct(time, coord, flag_mesh):
                x=coord[0]
                y=coord[1]
                a = 1

                if ((x>= 0.8 and x <= 0.9) and (y>= 0.1 and y<=0.2)):
                    value = 1.0
                else:
                    value = 0.000000001
                steady= True
                return value,steady;

        if ( flag == 'delta_10' ):
            def funct(time, coord, flag_mesh):
                x=coord[0]
                y=coord[1]
                a = 1

                if ((x>= 0.8 and x <= 0.9) and (y>= 0.1 and y<=0.2)):
                    value = 10.0
                else:
                    value = 0.000000001
                steady= True
                return value,steady;                                              


        if ( flag == 'delta_11' ):
            def funct(time, coord, flag_mesh):
                x=coord[0]
                y=coord[1]
                a = 1

                if ((x>= 0.1 and x <= 0.2) and (y>= 0.1 and y<=0.2)):
                    value = 5.0
                else:
                    value = 0.000000001
                steady= True
                return value,steady;

        if ( flag == 'delta_12' ):
            def funct(time, coord, flag_mesh):
                x=coord[0]
                y=coord[1]
                a = 1

                if ((x>= 0.4 and x <= 0.5) and (y>= 0.8 and y<=0.9)):
                    value = 10.0
                else:
                    value = 0.000000001
                steady= True
                return value,steady;

        if (flag == 'circle'):
            def funct(time,coord,flag_mesh):
                x=coord[0]
                y=coord[1]
                a=1
                value=1.0+(np.sin(a*np.pi*x))**2 + (np.cos(a*np.pi*y))**2
                steady=True
                return value,steady;


        if ( flag == 'multiple_deltas_2' ):
            def funct(time, coord, flag_mesh):
                x=coord[0]
                y=coord[1]
                a = 1

                if ((x>= 0.1 and x <= 0.2) and (y>= 0.8 and y<=0.9)):
                    value = 1.0
                elif ((x>= 0.5 and x <= 0.6) and (y>= 0.5 and y<=0.6)):
                    value = 5.0                
                elif ((x>= 0.8 and x <= 0.9) and (y>= 0.8 and y<=0.9)):
                    value = 1.0
                elif ((x>= 0.01 and x <= 0.1) and (y>= 0.4 and y<=0.5)):
                    value = 1.0                 
                elif ((x>= 0.01 and x <= 0.1) and (y>= 0.1 and y<=0.2)):
                    value = 1.0                 
                elif ((x>= 0.5 and x <= 0.6) and (y>= 0.8 and y<=0.9)):
                    value = 5.0                 
                elif ((x>= 0.5 and x <= 0.6) and (y>= 0.1 and y<=0.2)):
                    value = 5.0                 
                elif ((x>= 0.90 and x <= 0.99) and (y>= 0.5 and y<=0.6)):
                    value = 1.0                 
                elif ((x>= 0.8 and x <= 0.9) and (y>= 0.1 and y<=0.2)):
                    value = 1.0                 

                else:
                    value = 0.0000000001
                steady= True
                return value,steady; 


        if ( flag == 'multiple_deltas_3' ):
            def funct(time, coord, flag_mesh):
                x=coord[0]
                y=coord[1]
                a = 1

                if ((x>= 0.5 and x <= 0.6) and (y>= 0.8 and y<=0.9)):
                    value = 1.0
                elif ((x>= 0.5 and x <= 0.6) and (y>= 0.1 and y<=0.2)):
                    value = 1.0                               

                else:
                    value = 0.0000000001
                steady= True
                return value,steady; 



        if ( flag == 'circle_2' ):
            def funct(time,coord,flag_mesh):
                x=coord[0]
                y=coord[1]
                a=1
                value=1.0+(np.cos(a*np.pi*x))**2 + (np.sin(a*np.pi*y))**2
                steady=True
                return value,steady;

        if ( flag == 'Uniform' ):
            def funct(time,coord,flag_mesh):
                x=coord[0]
                y=coord[1]
                a=1
                value=1.0
                steady=True
                return value,steady;

                

        #
        # Kappa field in example 1 in fast marching tecnique
        #       
        if (flag == 'example1fmt'):
            def funct(time,coord,flag_mesh):
                x=coord[0]
                y=coord[1]
                if (x>0 and y>0):
                    value=np.sqrt(13*x**2+13*y**2+24*x*y)
                else:
                    value=2.0*np.sqrt(x**2+y**2)
                steady=True
                return value,steady;
        
        #
        # Kappa field in  test case 1 in 
        # An accurate discontinuous Galerkin method for solving point-source
        # Eikonal equation in 2-D heterogeneous anisotropic media
        # Bouteiller,1 M. Benjemaa,2 L. Metivier1,3 and J. Virieux1
        if (flag == 'bouteiller18_tc1'):
            def funct(time,coord,flag_mesh):
                x=coord[0]
                y=coord[1]
                z=coord[2]
                value=1.0/(1.0+0.5*y)
                steady=True
                return value,steady;
        
                #
        # Kappa field in  test case 2 in 
        # An accurate discontinuous Galerkin method for solving point-source
        # Eikonal equation in 2-D heterogeneous anisotropic media
        # Bouteiller,1 M. Benjemaa,2 L. Metivier1,3 and J. Virieux1
        if (flag == 'bouteiller18_tc2'):
            def funct(time,coord,flag_mesh):
                x=coord[0]
                y=coord[1]
                z=coord[2]
                value=1.0
                if (x-1.0)**2+(y-1.5)**2<0.25:
                    value=4
                steady=True
                return value,steady;

        #
        #  parabola centered at 0,0
        #
        if (flag == 'par'):
            def funct(time,coord,flag_mesh):
                x=coord[0]
                y=coord[1]
                value=1.0+x**2+y**2
                steady=True
                return value,steady;

        #
        #  parabola centered at 0,0
        #
        if (flag == 'prigozhin'):
            value_ellipse = float(extra_info)
            #print( value_ellipse)
            def funct(time,coord,flag_mesh):
                if (flag_mesh == 2):
                    value=value_ellipse
                else:
                    value=1.0
            
                steady=True
                return value,steady;

        #
        #  parabola centered at origin
        #
        if (flag == 'random'):
            info_path=common.remove_comments(extra_info,'!')
            input_file=open(info_path,'r')
            input_lines = input_file.readlines()
            min=float(common.remove_comments(input_lines[0],'!'))
            max=float(common.remove_comments(input_lines[1],'!'))
            iseed=float(common.remove_comments(input_lines[2],'!'))
        def funct(time,coord,flag_mesh):
                if (iseed != 0 ):
                    seed(iseed+flag_mesh)
                else:
                    seed()
                value=min+(max-min)*random()
                steady=True
                return value,steady;
            

        #
        # parabola centered at (0.5, 0.5)
        # 
        if (flag == 'par0505'):
            def funct(time,coord,flag_mesh):
                x=coord[0]
                y=coord[1]
                value=0.1+(x-0.5)**2+(y-0.5)**2
                steady=True
                return value,steady;

        #
        # oscilating data  
        # 
        if (flag == 'sin'):
            def funct(time,coord,flag_mesh):
                x=coord[0]
                y=coord[1]
                value=2.0+sin(8*np.pi*x)+sin(8*np.pi*y)
                steady=True
                return value,steady;

            
        # double parabola
        if (flag == 'par_eastwest'):
            def funct(coord,flag_mesh):
                x=coord[0]
                y=coord[1]
                # west
                if ( x <= 0.5 ):
                    value=0.01+2*(x-0.5)**2+y**2;
                # east
                if ( x > 0.5 ):
                    value=0.01+2*(x-0.5)**2+(y-0.5)**2;
                steady=True
                return value,steady;

        

    
        # y_branch
        if (flag == 'ybranch'):
            # read source-sink point
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
        
            def funct(coord,flag_mesh):
                x=coord[0]
                y=coord[1]
                value = 1.0e-10
                for segments in branch_vertices :
                    distance=mt.distance_segment_point(
                        branch_points[segments[0]],
                        branch_points[segments[1]],
                        point)
                    if ( distance < length / 1.5 ):
                        value = 1.0e10
                    
                steady=True
                return value,steady;
    
    return funct;



