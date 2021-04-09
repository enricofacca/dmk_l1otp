#!/usr/bin/env python
import re
import numpy as np
import jw_meshtools as mt
from meshpy.tools import uniform_refine_triangles
import numpy.linalg as la
import sys
from pyvtk import *
import meshpy.triangle as triangle
import example_grid as ex_grid
import make_grid as make

################################################
# domain definition

def assembly(example_name, length, aligned=None):

    # default 
    def optdens(x,y):
            return 1.0e30;

    if (example_name=='one2twocircles') :
        points,vertices,mesh=ex_grid.one2twocircles(length)
        def sink(x,y):
            fvalue=2.0
            return fvalue;
        def source(x,y):
            fvalue=2.0
            return fvalue;
       
        def kappa(flag,x,y):
            kappa_value=1.0
            if (flag != 0) or (flag != 1) or (flag != -1):
                fvalue=float(flag)
                return kappa_value;


    
    if (example_name=='rectangles_cost'):    
        # build geomtry
        points,vertices,mesh,base,width,height,shift=ex_grid.rectangles(length,aligned)
      
        #build inputs
        #define sink and source
        def sink(x,y):
            fvalue=2.0
            return fvalue;
        def source(x,y):
            fvalue=2.0
            return fvalue;

        def kappa(flag,x,y):
            kappa_value=1.0
            if (flag != 0) or (flag != 1) or (flag != -1):
                fvalue=float(flag)
            return kappa_value;

        def optdens(x,y):#,base,width,height,shift,sink,source):
            value=0.0
            if (x >= base[0]) and (x <= base[0]+width) and \
               (y >= base[1]) and (y <= base[1]+height):
                value=(x-base[0])*source(x,y)
            if (x >= base[0]+shift) and (x <= base[0]+width+shift) and \
               (y >= base[1]) and (y <= base[1]+height):
                value=(base[0]+width+shift-x)*sink(x,y)
            if (x >= base[0]+width) and (x <= base[0]+shift) and \
               (y >= base[1]) and (y <= base[1]+height):
                value=(width)*source(x,y)
            return value;

         
    if (example_name=='rectangles_cont'):   
        points,vertices,mesh,base,width,height,shift=ex_grid.rectangles(length,aligned)

        #build inputs
        #define sink and source
        def sink(x,y):
            fvalue=2.0
            return fvalue;
        def source(x,y):
            fvalue=2.0
            return fvalue;

        def kappa(flag,x,y):
            kappa_value=1.0
            if (flag != 0) or (flag != 1) or (flag != -1):
                fvalue=float(flag)
            return kappa_value;

        def optdens(x,y,base,width,height,shift,sink,source):
            value=0.0
            if (x >= base[0]) and (x <= base[0]+width) and \
               (y >= base[1]) and (y <= base[1]+height):
                value=(x-base[0])*source(x,y)
            if (x >= base[0]+shift) and (x <= base[0]+width+shift) and \
               (y >= base[1]) and (y <= base[1]+height):
                value=(base[0]+width+shift-x)*sink(x,y)
            if (x >= base[0]+width) and (x <= base[0]+shift) and \
               (y >= base[1]) and (y <= base[1]+height):
                value=(width)*source(x,y)
            return value;
                
    return points, vertices, mesh, sink, source, kappa, optdens;




