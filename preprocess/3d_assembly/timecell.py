#!/usr/bin/env python
import numpy as np

#
# prgrams to build functions
# taking real value on grid elements
#


#
# evaluate function at certain time 
#
def build(values,steady,
          functions,time,coord_eval,flags=None):
    if ( flags.all() == None ) :
        steady=True
        values[:]=0.0
        for i in range(len(values)):
            for ifun in range(len(functions)):
                value,steady_fun = functions[ifun](
                    time,
                    coord_eval[i][:],
                    i)
                values[i]+=value
                steady = steady and steady_fun
    else: 
        steady=True
        values[:]=0.0
        for i in range(len(values)):
            for ifun in range(len(functions)):
                value,steady_fun = functions[ifun](
                    time,
                    coord_eval[i][:],
                    flags[i])
                values[i]+=value
                steady = steady and steady_fun
        


    return values,steady;


# build examples of functions from
#
# function_scale [0,+\infty[ \mapsto \REAL^{\Dim}
# 
# where \Dim is the Cell number of the grid
# 
# flag_scalar (string) : define the example
# if flag_scalar is float number 
# function(t)=float(flag_scalar)
# other case aare defined by the flag used.
#
# 
# The string "extra_info" it may contain the path to 
# a file contiang additional information 
# defing the function ( constants, range of time, etc. )
#


# build standard tdens_0  given id_tdens0
def example(flag,extra_info):
    # if is const contained in flag_kappa
    try :
        cnst=float(flag)
        def func(time,coord,flag_mesh):
            value=cnst
            steady=True
            return cnst, steady;
    #
    # user defined example
    # 
    except ValueError:
        #  parabola centered at origin
        if (flag == 'par'):
            def func(time,coord,flag_mesh):
                x=coord[0]
                y=coord[1]
                value=1.0+x**2+y**2
                steady=True
                return value,steady;
        else:
            print 'flag =' + flag + 'not found'
    
    return func;


def write2file(file_out,time,head,data,steady=False):
    if ( head ) :
        file_out.write(
            str(data.shape[1]) + " " + 
            str(data.shape[0]) + " !  dim ndata" 
            +"\n")
    else:
        temp=np.zeros(data.shape[0])
        for i in range(data.shape[0]):
            temp[i]=np.linalg.norm(data[i,:])
        ndata=(abs(temp) != 0.0).sum()
        file_out.write("time    "+ str(time)+" \n")
        file_out.write(str(ndata)+" \n")
        for i in range(len(data)):
            if ( np.sum(np.abs(data[i][:])) != 0.0):
                file_out.write(str(i+1)+" " + 
                               " ".join(map(str,data[i][:])) +"\n")
        if( steady ):
            file_out.write("time " + str(1.e30)+"\n")


