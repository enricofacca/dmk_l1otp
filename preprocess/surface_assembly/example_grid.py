#!/usr/bin/env python
import re
import numpy as np
import meshtools_surface as mt
import sys
from pyvtk import *
import meshpy.triangle as make_triangle
import common 
import branch


def example_grid(flag_grid,length,extra_info=None,minimum_angle=27.0):
    sources=[];
    sinks=[];
    kappas=[];
    tdens0_functions=[];
    optdens_functions=[];
    dirac_source_points=[];
    dirac_sinks_points=[];
    dirac_source_values=[];
    dirac_sink_values=[];

    if (flag_grid == 'unitsquare'):
        # domain
        print flag_grid
        points,vertices=mt.MyRectangle([0,0],[1,1])
        
        #set info mesh    
        info = make_triangle.MeshInfo()
        info.set_points(points)
        info.set_facets(vertices)
        info.regions.resize(1)

        #build mesh
        mesh = make_triangle.build(info,
                                   attributes=True, 
                                   max_volume=(length**2/2.0))

    
    if (flag_grid == 'rectangles'):
        # domain
        print flag_grid
        points,vertices=mt.MyRectangle([0,0],[1,1])
    
        # supportrs info
        base=np.array([0.125,0.25])
        width=0.25
        height=0.5
        shift=0.5
        forcing_value=2.0


        # Support \Source
        sw_source=np.array(base)
        ne_source=sw_source+(width,height)
        psource=(sw_source+ne_source)/2.0
        p1,v1=mt.MyRectangle(sw_source,ne_source)
        se_source=p1[1]
        nw_source=p1[3]
        label_source=1
    
        # Support \Sink
        sw_sink=sw_source+(shift,0.0)
        ne_sink=sw_sink+(width,height)
        psink=(sw_sink+ne_sink)/2.0
        p2,v2=mt.MyRectangle(sw_sink,ne_sink)
        se_sink=p2[1]
        nw_sink=p2[3]
        label_sink=-1
    

        # union   
        points,vertices=mt.AddCurves(points,vertices,p1,v1)
        points,vertices=mt.AddCurves(points,vertices,p2,v2)

    

        #set info mesh    
        info = make_triangle.MeshInfo()
        info.set_points(points)
        info.set_facets(vertices)
        info.regions.resize(2)
        info.regions[0] = (list(psource) + [label_source, 1])
        info.regions[1] = (list(psink) + [label_sink,   1]) 

        #build mesh
        mesh = make_triangle.build(info,
                                   attributes=True, 
                                   max_volume=(length**2/2.0))

    if (flag_grid == 'rect_cnst'):
        # domain
        print flag_grid
        points,vertices=mt.MyRectangle([0,0],[1,1])
    
        # supportrs info
        base=np.array([0.125,0.25])
        width=0.25
        height=0.5
        shift=0.5
        forcing_value=2.0


        # Support \Source
        sw_source=np.array(base)
        ne_source=sw_source+(width,height)
        psource=(sw_source+ne_source)/2.0
        p1,v1=mt.MyRectangle(sw_source,ne_source)
        se_source=p1[1]
        nw_source=p1[3]
        label_source=1
    
        # Support \Sink
        sw_sink=sw_source+(shift,0.0)
        ne_sink=sw_sink+(width,height)
        psink=(sw_sink+ne_sink)/2.0
        p2,v2=mt.MyRectangle(sw_sink,ne_sink)
        se_sink=p2[1]
        nw_sink=p2[3]
        label_sink=-1
    

        # union   
        points,vertices=mt.AddCurves(points,vertices,p1,v1)
        points,vertices=mt.AddCurves(points,vertices,p2,v2)

    

        #set info mesh    
        info = make_triangle.MeshInfo()
        info.set_points(points)
        info.set_facets(vertices)
        info.regions.resize(2)
        info.regions[0] = (list(psource) + [label_source, 1])
        info.regions[1] = (list(psink) + [label_sink,   1]) 

        #build mesh
        mesh = make_triangle.build(info,
                                   attributes=True, 
                                   max_volume=(length**2/2.0))


        def source(time,x,y,flag_mesh):
            fvalue=0.0
            steady=True
            if (flag_mesh > 0 ):
                fvalue=forcing_value
            return fvalue, steady;
        sources.append(source)

        #  constant examples
        def sink(time,x,y,flag_mesh):
             fvalue=0.0
             steady=True
             if (flag_mesh < 0 ):
                 fvalue=forcing_value
                 return fvalue,steady;
        sinks.append(sink)


        def optdens(x,y,z,flag):
            value=0.0
            if ( (x >  base[0]        ) and 
                 (x <= base[0]+width  ) and 
                 (y >  base[1]        ) and 
                 (y <= base[1]+height ) )  :
                value=(x-base[0])*forcing_value
            if ( (x >  base[0]+width        ) and 
                 (x <= base[0]+width+shift  ) and 
                 (y >  base[1]              ) and 
                 (y <= base[1]+height       ) )  :
                value=(width)*forcing_value
            if ( (x >  base[0]+  width+shift ) and 
                 (x <= base[0]+2*width+shift ) and 
                 (y >  base[1]               ) and 
                 (y <= base[1]+height        ) )  :
                value=(base[0]+2*width+shift-x)*forcing_value
            return value;
        optdens_functions.append(optdens)


        

    if (flag_grid == 'rectangles_aligned'):
        # domain
        points,vertices=mt.MyRectangle([0,0],[1,1])
    
        # supportrs info
        base=np.array([0.125,0.25])
        width=0.25
        height=0.5
        shift=0.5


        # Support \Source
        sw_source=np.array(base)
        ne_source=sw_source+(width,height)
        psource=(sw_source+ne_source)/2.0
        p1,v1=mt.MyRectangle(sw_source,ne_source)
        se_source=p1[1]
        nw_source=p1[3]
        label_source=1
    
        # Support \Sink
        sw_sink=sw_source+(shift,0.0)
        ne_sink=sw_sink+(width,height)
        psink=(sw_sink+ne_sink)/2.0
        p2,v2=mt.MyRectangle(sw_sink,ne_sink)
        se_sink=p2[1]
        nw_sink=p2[3]
        label_sink=-1
    

        # union   
        points,vertices=mt.AddCurves(points,vertices,p1,v1)
        points,vertices=mt.AddCurves(points,vertices,p2,v2)

        
        p3,v3=mt.MyLineSegments(se_source,sw_sink)
        p4,v4=mt.MyLineSegments(ne_source,nw_sink)
        # union 
        points,vertices=mt.AddCurves(points,vertices,p3,v3)
        points,vertices=mt.AddCurves(points,vertices,p4,v4)
   

        #set info mesh    
        info = make_triangle.MeshInfo()
        info.set_points(points)
        info.set_facets(vertices)
        info.regions.resize(2)
        info.regions[0] = (list(psource) + [label_source, 1])
        info.regions[1] = (list(psink) + [label_sink,   1]) 

        #build mesh
        mesh = make_triangle.build(info, 
                                   attributes=True, 
                                   max_volume=(length**2/2.0))

    if (flag_grid == 'rect_cnst_aligned'):
        # domain
        points,vertices=mt.MyRectangle([0,0],[1,1])
    
        # supportrs info
        base=np.array([0.125,0.25])
        width=0.25
        height=0.5
        shift=0.5
        forcing_value = 2.0


        # Support \Source
        sw_source=np.array(base)
        ne_source=sw_source+(width,height)
        psource=(sw_source+ne_source)/2.0
        p1,v1=mt.MyRectangle(sw_source,ne_source)
        se_source=p1[1]
        nw_source=p1[3]
        label_source=1
    
        # Support \Sink
        sw_sink=sw_source+(shift,0.0)
        ne_sink=sw_sink+(width,height)
        psink=(sw_sink+ne_sink)/2.0
        p2,v2=mt.MyRectangle(sw_sink,ne_sink)
        se_sink=p2[1]
        nw_sink=p2[3]
        label_sink=-1
    

        # union   
        points,vertices=mt.AddCurves(points,vertices,p1,v1)
        points,vertices=mt.AddCurves(points,vertices,p2,v2)

        
        p3,v3=mt.MyLineSegments(se_source,sw_sink)
        p4,v4=mt.MyLineSegments(ne_source,nw_sink)
        # union 
        points,vertices=mt.AddCurves(points,vertices,p3,v3)
        points,vertices=mt.AddCurves(points,vertices,p4,v4)
   

        #set info mesh    
        info = make_triangle.MeshInfo()
        info.set_points(points)
        info.set_facets(vertices)
        info.regions.resize(2)
        info.regions[0] = (list(psource) + [label_source, 1])
        info.regions[1] = (list(psink) + [label_sink,   1]) 

        #build mesh
        mesh = make_triangle.build(info, 
                                   attributes=True, 
                                   max_volume=(length**2/2.0))


        def source(time,x,y,z,flag_mesh):
            fvalue=0.0
            steady=True
            if (flag_mesh > 0 ):
                fvalue=forcing_value
            return fvalue, steady;
        sources.append(source)

        #  constant examples
        def sink(time,x,y,z,flag_mesh):
            fvalue=0.0
            steady=True
            if (flag_mesh < 0 ):
                fvalue=forcing_value
            return fvalue,steady;
        sinks.append(sink)


        def optdens(coord,flag):
            x=coord[0]
            y=coord[1]
            z=coord[2]
            value=0.0
            #print x,y,z, base, width, height,forcing_value
            if ( (x >  base[0]        ) and 
                 (x <= base[0]+width  ) and 
                 (y >  base[1]        ) and 
                 (y <= base[1]+height ) )  :
                value=(x-base[0])*forcing_value
            if ( (x >  base[0]+width        ) and 
                 (x <= base[0]+shift        ) and 
                 (y >  base[1]              ) and 
                 (y <= base[1]+height       ) )  :
                value=(width)*forcing_value
            if ( (x >  base[0]+     +shift ) and 
                 (x <= base[0]+width+shift ) and 
                 (y >  base[1]               ) and 
                 (y <= base[1]+height        ) )  :
                value=(base[0]+width+shift-x)*forcing_value
            return value;
        optdens_functions.append(optdens)


   
    if (flag_grid == 'prigozhin'):
        # domain
        c0=[0.5,0.5]
        l0=0
        points,vertices=mt.MyRectangle([0,0],[1,1])
        
        #source
        c1=[0.3,0.65]
        r1=0.18
        l1=1
        p1, v1 = mt.CircleSegments(c1,r1,edge_length=length)
        
        #sink
        c2=[0.8,0.4]
        r2=0.3
        l2=-1
        p2, v2 = mt.CircleSegments(c2,r2,edge_length=length)
        for point in p2:
            point[0] = c2[0] + ( point[0] - c2[0] ) / 3.0

        #union
        points,vertices=mt.AddCurves(points,vertices,p1,v1)
        points,vertices=mt.AddCurves(points,vertices,p2,v2)
        
        #set info mesh
        info = make_triangle.MeshInfo()
        info.set_points(points)
        info.set_facets(vertices)
        info.regions.resize(3)
        info.regions[0] = (c0 + [l0, 1])
        info.regions[1] = (c1 + [l1, 1]) 
        info.regions[2] = (c2 + [l2, 1]) 
    
        # mesh
        mesh = make_triangle.build(info, 
                              attributes=True, 
                              max_volume=(length**2/2.0),
                              min_angle=minimum_angle)


 
    if (flag_grid == 'one2twocircles'):
        # domain
        points,vertices=mt.MyRectangle([0,0],[1,1])
        
        #source
        c1=[0.25,0.5]
        r1=0.2
        l1=1
        p1, v1 = mt.CircleSegments(c1,r1,edge_length=length)
        
        #sink
        c2=[0.75,0.35]
        r2=0.2
        l2=-1
        p2, v2 = mt.CircleSegments(c2,r2,edge_length=length)


        c3=[0.85,0.85]
        r3=0.1
        l3=-1
        p3, v3 = mt.CircleSegments(c3,r3,edge_length=length)

        #union
        points,vertices=mt.AddCurves(points,vertices,p1,v1)
        points,vertices=mt.AddCurves(points,vertices,p2,v2)
        points,vertices=mt.AddCurves(points,vertices,p3,v3)


        #set info mesh
        info = make_triangle.MeshInfo()
        info.set_points(points)
        info.set_facets(vertices)
        info.regions.resize(3)
        info.regions[0] = (c1 + [l1, 1])
        info.regions[1] = (c2 + [l2, 1]) 
        info.regions[2] = (c3 + [l3, 1]) 
    
        # mesh
        mesh = make_triangle.build(info, 
                              attributes=True, 
                              max_volume=(length**2/2.0),
                              min_angle=minimum_angle)



    #return points, vertices, mesh;#, base, width,height,shift;

    if (flag_grid == 'one2fourcircles'):
        # domain
        points,vertices=mt.MyRectangle([0,0],[1,1])
        
        #source
        c1=[0.5,0.5]
        r1=0.25
        l1=1
        p1, v1 = mt.CircleSegments(c1,r1,edge_length=length)
        def source(t,x,y,z,flag_mesh):
            fvalue=0.0
            steady=True
            if (flag_mesh == l1) : 
                fvalue=max(0,(r1-np.sqrt((x-c1[0])**2+(y-c1[1])**2))/r1)
            return fvalue,steady;
        sources.append(source)
        
        #sink 
        # sw
        c2=[0.35,0.15]
        r2=0.1
        l2=-1
        p2, v2 = mt.CircleSegments(c2,r2,edge_length=length)
        def sink0(t,x,y,z,flag_mesh):
            fvalue=0.0
            steady=True
            if (flag_mesh == l2):
                 fvalue=max(0,(r2-np.sqrt((x-c2[0])**2+(y-c2[1])**2))/r2)
            return fvalue,steady;
        sinks.append(sink0)

        # se
        c3=[0.8,0.2]
        r3=0.12
        l3=-2
        p3, v3 = mt.CircleSegments(c3,r3,edge_length=length)
        def sink1(t,x,y,z,flag_mesh):
            fvalue=0.0
            steady=True
            if (flag_mesh == l3):
                 fvalue=max(0,(r3-np.sqrt((x-c3[0])**2+(y-c3[1])**2))/r3)
            return fvalue,steady;
        sinks.append(sink1)


        # ne
        c4=[0.85,0.85]
        r4=0.14
        l4=-3
        p4, v4 = mt.CircleSegments(c4,r4,edge_length=length)
        def sink2(t,x,y,z,flag_mesh):
            fvalue=0.0
            steady=True
            if (flag_mesh == l4):
                fvalue=max(0,(r4-np.sqrt((x-c4[0])**2+(y-c4[1])**2))/r4)
            return fvalue,steady;
        sinks.append(sink2)

        # nw
        c5=[0.15,0.75]
        r5=0.1
        l5=-4
        p5, v5 = mt.CircleSegments(c5,r5,edge_length=length)
        def sink3(t,x,y,z,flag_mesh):
            fvalue=0.0
            steady=True
            if (flag_mesh == l5):
                fvalue=max(0,(r5-np.sqrt((x-c5[0])**2+(y-c5[1])**2))/r5)
            return fvalue,steady;
        sinks.append(sink3)


        #union
        points,vertices=mt.AddCurves(points,vertices,p1,v1)
        points,vertices=mt.AddCurves(points,vertices,p2,v2)
        points,vertices=mt.AddCurves(points,vertices,p3,v3)
        points,vertices=mt.AddCurves(points,vertices,p4,v4)
        points,vertices=mt.AddCurves(points,vertices,p5,v5)



        #set info mesh
        info = make_triangle.MeshInfo()
        info.set_points(points)
        info.set_facets(vertices)
        info.regions.resize(5)
        info.regions[0] = (c1 + [l1, 1])
        info.regions[1] = (c2 + [l2, 1]) 
        info.regions[2] = (c3 + [l3, 1]) 
        info.regions[3] = (c4 + [l4, 1]) 
        info.regions[4] = (c5 + [l5, 1]) 
    

        # mesh
        mesh = make_triangle.build(info, 
                              attributes=True, 
                              max_volume=(length**2/2.0),
                              min_angle=minimum_angle)

    if (flag_grid == 'one2foursquares'):
        # domain
        points,vertices=mt.MyRectangle([0,0],[1,1])
        
        #source
        c1=[0.5,0.5]
        r1=0.25
        l1=1
        base1_sw=[c1[0]-r1/2.0,c1[1]-r1/2.0]
        base1_ne=[c1[0]+r1/2.0,c1[1]+r1/2.0]
        p1, v1 = mt.MyRectangle(base1_sw,base1_ne)
        def source(t,x,y,z,flag_mesh):
            fvalue=0.0
            steady=True
            if (flag_mesh == l1) : 
                fvalue=max(0,(r1/2.0-max(abs(x-c1[0]),abs(y-c1[1])))/(r1/2.0))
            return fvalue,steady;
        sources.append(source)
        
        #sink 
        # sw
        c2=[0.35,0.15]
        r2=0.1
        l2=-1
        base2_sw=[c2[0]-r2/2.0,c2[1]-r2/2.0]
        base2_ne=[c2[0]+r2/2.0,c2[1]+r2/2.0]
        p2, v2 = mt.MyRectangle(base2_sw,base2_ne)
        def sink0(t,x,y,z,flag_mesh):
            fvalue=0.0
            steady=True
            if (flag_mesh == l2):
                fvalue=max(0,(r2/2.0-max(abs(x-c2[0]),abs(y-c2[1])))/(r2/2.0))
            return fvalue,steady;
        sinks.append(sink0)

        # se
        c3=[0.8,0.2]
        r3=0.12
        l3=-2
        base3_sw=[c3[0]-r3/2.0,c3[1]-r3/2.0]
        base3_ne=[c3[0]+r3/2.0,c3[1]+r3/2.0]
        p3, v3 = mt.MyRectangle(base3_sw,base3_ne)
        def sink1(t,x,y,z,flag_mesh):
            fvalue=0.0
            steady=True
            if (flag_mesh == l3):
                fvalue=max(0,(r3/2.0-max(abs(x-c3[0]),abs(y-c3[1])))/(r3/2.0))
            return fvalue,steady;
        sinks.append(sink1)


        # ne
        c4=[0.85,0.85]
        r4=0.14
        l4=-3
        base4_sw=[c4[0]-r4/2.0,c4[1]-r4/2.0]
        base4_ne=[c4[0]+r4/2.0,c4[1]+r4/2.0]
        p4, v4 = mt.MyRectangle(base4_sw,base4_ne)
        def sink2(t,x,y,z,flag_mesh):
            fvalue=0.0
            steady=True
            if (flag_mesh == l4):
                fvalue=max(0,(r4-np.sqrt((x-c4[0])**2+(y-c4[1])**2))/r4)
            return fvalue,steady;
        sinks.append(sink2)

        # nw
        c5=[0.2,0.7]
        r5=0.1
        l5=-4
        base5_sw=[c5[0]-r5/2.0,c5[1]-r5/2.0]
        base5_ne=[c5[0]+r5/2.0,c5[1]+r5/2.0]
        p5, v5 = mt.MyRectangle(base5_sw,base5_ne)
        def sink3(t,x,y,z,flag_mesh):
            fvalue=0.0
            steady=True
            if (flag_mesh == l5):
                fvalue=max(0,(r5/2.0-max(abs(x-c5[0]),abs(y-c5[1])))/(r5/2.0))
            return fvalue,steady;
        sinks.append(sink3)


        #union
        points,vertices=mt.AddCurves(points,vertices,p1,v1)
        points,vertices=mt.AddCurves(points,vertices,p2,v2)
        points,vertices=mt.AddCurves(points,vertices,p3,v3)
        points,vertices=mt.AddCurves(points,vertices,p4,v4)
        points,vertices=mt.AddCurves(points,vertices,p5,v5)



        #set info mesh
        info = make_triangle.MeshInfo()
        info.set_points(points)
        info.set_facets(vertices)
        info.regions.resize(5)
        info.regions[0] = (c1 + [l1, 1])
        info.regions[1] = (c2 + [l2, 1]) 
        info.regions[2] = (c3 + [l3, 1]) 
        info.regions[3] = (c4 + [l4, 1]) 
        info.regions[4] = (c5 + [l5, 1]) 
    

        # mesh
        mesh = make_triangle.build(info, 
                              attributes=True, 
                              max_volume=(length**2/2.0),
                              min_angle=minimum_angle)

    if (  flag_grid == 'ybranch'):
        # domain
        points,vertices=mt.MyRectangle([0,0],[1,1])

        # read points
        info_path=common.remove_comments(extra_info,'!')

        input_file=open(info_path,'r')
        input_lines = input_file.readlines()

        branch_points=[]
        branch_vertices=[]
        branch_values=[]
        for line in input_lines[0:3]:
            coord_dirac = [float(w) for w in line.split()[0:2]]
            branch_value = float(line.split()[2])
            branch_points.append(coord_dirac)
            branch_values.append(branch_value)            
            
        # mesh aligned or not with exact solution
        aligned=common.remove_comments(input_lines[3],'!')
        alpha=float(common.remove_comments(input_lines[4],'!'))
        if ( aligned == '1' ):
            # Build the optimal solution
            # add branching point (if exists)
            branch_points,branch_vertices=branch.y_branch(branch_points,branch_values,alpha)

        input_file.close()
            

        #Union
        points,vertices=mt.AddCurves(
            points,vertices,branch_points,branch_vertices)

        #set info mesh    
        info = make_triangle.MeshInfo()
        info.set_points(points)
        info.set_facets(vertices)
        info.regions.resize(1)
        info.regions[0] = ( [0.0, 0.0] + [0, 1])

        
        #build mesh
        mesh = make_triangle.build(
            info, attributes=True, max_volume=(length**2/2.0))
                

    if (  flag_grid == 'dirac'):
        # domain
        p,v=mt.MyRectangle([0,0],[1,1])
        # read dirac points
        input_file=open(dirac_path,'r')
        input_lines = input_file.readlines()
        ndirac       = int(common.remove_comments(input_lines[0],'!'))
        for line in input_lines[1:]:
            coord_dirac = [float(w) for w in line.split()[0:2]]
            p.append(coord_dirac)
            input_file.close()


    return points, vertices, mesh, sources, sinks, kappas,tdens0_functions,optdens_functions;





    # def optdens(x,y,z,base,width,height,shift,sink,source):
    #     value=0.0
    #     if (x >= base[0]) and (x <= base[0]+width) and \
    #        (y >= base[1]) and (y <= base[1]+height):
    #         value=(x-base[0])*source(x,y)
    #     if (x >= base[0]+shift) and (x <= base[0]+width+shift) and \
    #        (y >= base[1]) and (y <= base[1]+height):
    #         value=(base[0]+width+shift-x)*sink(x,y)
    #     if (x >= base[0]+width) and (x <= base[0]+shift) and \
    #        (y >= base[1]) and (y <= base[1]+height):
    #         value=(width)*source(x,y)
    #     return value;
