#!/usr/bin/env python
import re
import numpy as np
import jw_meshtools as mt
from meshpy.tools import uniform_refine_triangles
import numpy.linalg as la
import sys
from pyvtk import *
import meshpy.triangle as triangle

def make_kappa(flags,bar_tria,kappa_func):
    kappa_tria=np.empty([len(bar_tria)])
    for itria in range(len(bar_tria)):
        kappa_tria[itria]= kappa_func(flags[itria],bar_tria[itria][0],bar_tria[itria][1])
    return kappa_tria;


def make_optdens(bar_tria, optdens):
    optdens_tria=np.empty([len(bar_tria)])
    for itria in range(len(bar_tria)):
        optdens_tria[itria]=optdens(bar_tria[itria][0],bar_tria[itria][1])#,\
                             #base,width,height,shift,sink,source)
    return optdens_tria;


def write_measure(measure,bar_tria,area_tria,filename):
    nmeasure=(abs(measure) > 0.0).sum()
    file_out = open(filename, 'w')
    file_out.write(str(nmeasure)+" \n")
    for itria in range(len(measure)):
        if ( measure[itria] > 0.0):
            file_out.write(str(itria+1)+" "+
                           str(bar_tria[itria][0]) + " " + 
                           str(bar_tria[itria][1]) + " " + 
                           str(measure[itria]*area_tria[itria]) +"\n")
    file_out.close()
    
def write_time_inputs(file_out,time,data=None,steady=False,head=False):
    if ( head ) :
        file_out.write(
            str(data.ndim)+ " " + str(len(data)) + " !   dimdata, ndata" +"\n")
    else:
        if (len(data) != 1):
            ndata=(abs(data) > 0.0).sum()
        else:
            ndata=1
        file_out.write("time    "+ str(time)+" \n")
        file_out.write(str(ndata)+" \n")
        for i in range(len(data)):
            if ( abs(data[i]) > 0.0):
                file_out.write(str(i+1)+" " + str(data[i]) +"\n")
        if( steady ):
            file_out.write("time " + str(1.e30)+"\n")

def write_data(data,filename):
    file_out = open(filename, 'w')
    file_out.write(str(data.ndim)+" "+str(len(data))+" !dim number data \n")
    for itria in range(len(data)):
        file_out.write(str(data[itria]) +"\n")
    file_out.close()



        
        
    
    
        
    

