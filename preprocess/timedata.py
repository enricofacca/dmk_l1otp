#!/usr/bin/env python
import numpy as np

#
# prgrams to build functions
# taking real value on grid elements
#


#
# evaluate list of functions at given time and
# on the coordinate passed
# E.g. inputs:
#       functions=[f1,f2,f3]
#       time=0.0
#       coord_eval=[[0,0,0], [1,1,1]]
#     outputs:
#       values[0]=f1(0.0,[0,0,0])+f2(0.0,[0,0,0])+f3(0.0,[0,0,0])
#       values[1]=f1(0.0,[1,1,1])+f2(0.0,[1,1,1])+f3(0.0,[1,1,1])
#       steady=defined by f1,f2,f3
def build(values,steady,
          functions,time,coord_eval,flags=None):
    #
    #
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

#
#
# write head in fiel containg time data
#
# E.g.
# 1 1  !  ndim ndata" 
def writehead2file(file_out,data):
    file_out.write(
        str(data.shape[1]) + " " + 
        str(data.shape[0]) + " !  ndim ndata" 
        +"\n")
    
#
# write data( with dimension ndim x ndata) to open fileout 
# in format:
#
# time t
# nnz 
# i1 d_{1,i1} d_{2,i1} ... d_{ndim,i1}  
# i2 d_{1,i2}    d_{2,i2} ... d_{ndim,i2}  
# ..
# innz d_{1,nnz} d_{2,nnz} ... d_{ndim,nnz}
# 
# If steady_state is true, it adds a talign string with
#
# time  1.0e30 
# 
def write2file(file_out,time,data,steady=False):
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


def read_steady_timedata(filename):
   fin=open(str(filename), 'r')
   lines = fin.readlines()
   dimdata=1
   ndata=int(lines[0].split()[1])
   data=np.zeros(ndata)   
   ninputs= int(lines[2])
   for i in range(ninputs):
      line=lines[3+i]
      idata=int(line.split()[0])-1
      data[idata] = float(line.split()[1])
   fin.close()
   return data;
