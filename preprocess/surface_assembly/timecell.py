#!/usr/bin/env python
import numpy as np

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
# evaluate function at certain time 
#
def build(cell_values,steady,
          functions,time,bar_cell,flags):
    steady=True
    for icell in range(len(cell_values)):
        for ifun in range(len(functions)):
            value,steady_fun = functions[ifun](
                time,
                bar_cell[icell][:],
                flags[icell])
            cell_values[icell]+=value
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
        ndata=(abs(data) != 0.0).sum()
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
        #  parabola centered at origin
        if (flag == 'par'):
            def funct(time,coord,flag_mesh):
                x=coord[0]
                y=coord[1]
                value=1.0+x**2+y**2
                steady=True
                return value,steady;
        #  parabola centered at origin
        
        if (flag == 'virieux'):
            print flag
            def funct(time,coord,flag_mesh):
                x=coord[0]
                y=coord[1]
                z=coord[2]
                value=1.0
                if ( (x-0.5)**2+(y-0.5)**2+(z-0.5)**2<0.10):
                    value=4.0
                steady=True
                return value,steady;
        else:
            print 'flag =' + flag + 'not found'
    
    return funct;




