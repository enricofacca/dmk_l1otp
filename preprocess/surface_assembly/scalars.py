#!/usr/bin/env python
import numpy as np
import common as common

#
# prgrams to build functions
# taking real value, possibly in several dimension, that 
#


#
# evaluate scalar function at certain time 
#
def make_scalars(function,time):
    scalars, steady = function(time)
    return scalars, steady;


# build standard function from
#
# function_scale [0,+\infty[ \mapsto \REAL^{\Dim}
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

def example(flag_scalar,extra_info):
    # if is const contained in falg_pflux0
    try :
        scalar0=float(flag_scalar)
        def function_scalar(time):
            scalar_value=scalar0
            steady=True
            return [scalar_value], steady;
    #
    # user defined example
    # 
    except ValueError:
        # scalar oscillating
        if (flag_scalar == 'sin'):
            def function_scalar(time):
                scalar_value=(2.0+np.sin(time))
                steady=False
                return [scalar_value], steady;
        # scalar oscillating
        if (flag_scalar == 'oneandhalf'):
            def function_scalar(time):
                scalar_value1=1.0
                scalar_value2=0.5
                steady=True
                return [scalar_value1,scalar_value2], steady;

        else:
            print 'flag ='+ flag_scalar + 'not found'
        
    return function_scalar;


##########################
def build_and_write(flag,extra,tzero,tmax,fname):
    if (flag!='no') or (flag==''):
        # define functions
        function=example(flag,extra)

        # open file and write dimensions
        file_out=open(fname, 'w')
        # sample to get dimension and write head
        time = tzero
        values, steady = function(time)
        write_time_inputs(file_out,time,True,values,steady)
        
        #cycle in time
        time=tzero
        steady=False
        while (time <= tmax) and (not steady ):
            values, steady = function(time)
            write_time_inputs(file_out,time,False,values,steady)
            time=common.next_time(time)
        
        #close 
        file_out.close()

def write_time_inputs(file_out,time,head,scalars,steady=False):
    if ( head ) :
        file_out.write(
            str(len(scalars)) + " " + 
            str(1)            + " !  dim " 
            +"\n")
    else:
        nscalars=1
        file_out.write("time    "+ str(time)+" \n")
        file_out.write(str(1)+" \n")
        file_out.write(str(1)+" " + " ".join(map(str,scalars)) +"\n")
        if( steady ):
            file_out.write("time " + str(1.e30)+"\n")
