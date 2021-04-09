#!/usr/bin/env python
from __future__ import absolute_import
from six.moves import range

import re
import numpy as np
import meshtools as mt
import sys
import meshpy.triangle as triangle
from meshpy.tet import MeshInfo
from meshpy.tet import Options
from meshpy.tet import build
from meshpy.geometry import \
    generate_surface_of_revolution, EXT_OPEN, \
    GeometryBuilder


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

    if (flag_grid == 'rect_cnst'):
        # domain
        coord=[];
        topol=[];


        # supportrs info
        base_optdens=np.array([0.125,0.25,0.25])
        xwidth=0.25
        ywidth=0.5
        zwidth=0.5
        shift=0.5


        # domain
        base=np.array([0,0,0])
        xl=1
        yl=1
        zl=1
        coord,topol=mt.MyParallelepiped(base,xl,yl,zl)

        # source
        base=base_optdens
        xl=xwidth
        yl=ywidth
        zl=zwidth
        psource = base + np.array([xl,yl,zl])/2.0
        label_source = 1
        coord1,topol1=mt.MyParallelepiped(base,xl,yl,zl)

        # sink
        base=base_optdens+np.array([shift,0,0])
        xl=xwidth
        yl=ywidth
        zl=zwidth
        psink= base + np.array([xl,yl,zl])/2.0
        label_sink=-1
        coord2,topol2=mt.MyParallelepiped(base,xl,yl,zl)
        
        # join
        coord,topol=mt.Add3dSurface(coord,topol,coord1,topol1)
        coord,topol=mt.Add3dSurface(coord,topol,coord2,topol2)


        #set info mesh    
        info = MeshInfo()
        info.set_points(coord)
        info.set_facets(topol)
        info.regions.resize(3)
        info.regions[0] = ([0.5,0.5,0.5] + [0, 1])
        info.regions[1] = (list(psource) + [label_source, 1])
        info.regions[2] = (list(psink) + [label_sink,   1]) 

        #build mesh
        mesh = build(info, attributes=True, max_volume=(length**3/6.0))

    if (flag_grid == 'rect_cnst_aligned'):
        # domain
        coord=[];
        topol=[];


        # supportrs info
        base_optdens=np.array([0.125,0.25,0.25])
        xwidth=0.25
        ywidth=0.5
        zwidth=0.5
        shift=0.5

        # # supportrs info
        # base_optdens=np.array([0.1,0.1,0.1])
        # xwidth=0.25
        # ywidth=0.8
        # zwidth=0.8
        # shift=0.3


        # domain
        base=np.array([0,0,0])
        xl=1
        yl=1
        zl=1
        coord,topol=mt.MyParallelepiped(base,xl,yl,zl)

        # source
        base=base_optdens
        xl=xwidth
        yl=ywidth
        zl=zwidth
        psource = base + np.array([xl,yl,zl])/2.0
        label_source = 1
        coord1,topol1=mt.MyParallelepiped(base,xl,yl,zl)

        # sink
        base=base_optdens+np.array([shift,0,0])
        xl=xwidth
        yl=ywidth
        zl=zwidth
        psink= base + np.array([xl,yl,zl])/2.0
        label_sink=-1
        coord2,topol2=mt.MyParallelepiped(base,xl,yl,zl)

        # middle
        base=base_optdens+np.array([xwidth,0,0])
        xl=shift-xwidth
        yl=ywidth
        zl=zwidth
        pmiddle= base + np.array([xl,yl,zl])/2.0
        label_middle=0
        coord3,topol3=mt.MyParallelepiped(base,xl,yl,zl)

        
        # join
        coord,topol=mt.Add3dSurface(coord,topol,coord1,topol1)
        coord,topol=mt.Add3dSurface(coord,topol,coord3,topol3)
        coord,topol=mt.Add3dSurface(coord,topol,coord2,topol2)
        # print('coord')
        # for i in range (0,len(coord)):
        #     print i,coord[i]
        # print('topol')
        # for i in range (0,len(topol)):
        #     print i,topol[i]

        # topol=[
        #     [0, 1, 2, 3],
        #     [4, 5, 6, 7],
        #     [0, 1, 5, 4],
        #     [2, 6, 7, 3],
        #     [1, 2, 6, 5],
        #     [0, 4, 7, 3], # cubo
        #     [8, 9, 10, 11],
        #     [12, 13, 14, 15],
        #     [8, 9, 13, 12],
        #     [10, 14, 15, 11],
        #     [9, 10, 14, 13],
        #     [8, 12, 15, 11], # source
        #     [9, 16, 17, 10],
        #     [13, 18, 19, 14],
        #     [9, 16, 18, 13],
        #     [10, 14, 19, 17],
        #     [16, 18, 19, 17], # middle
        #     [16, 20, 21, 17],
        #     [18, 22, 23, 19],
        #     [20, 22, 23, 21],
        #     [16, 18, 22, 20],
        #     [17, 19, 23, 21], # sink
        # ]

            
        #set info mesh    
        info = MeshInfo()
        info.set_points(coord)
        info.set_facets(topol)
        info.regions.resize(4)
        info.regions[0] = ([0.005,0.005,0.005] + [0, 1])
        info.regions[1] = (list(psource) + [label_source, 1])
        info.regions[2] = (list(psink) + [label_sink,   1]) 
        info.regions[3] = (list(pmiddle) + [label_middle,   1])
        #build mesh
        mesh = build(info, attributes=True, max_volume=(length**3/6.0))


        cnst_value=2.0
        
        def source(time,coord,flag_mesh):
            fvalue=0.0
            x=coord[0]
            y=coord[1]
            z=coord[2]
            steady=True
            if (flag_mesh >0 ):
                fvalue=cnst_value
            return fvalue, steady;
        sources.append(source)

        def sink(time,coord,flag_mesh):
            fvalue=0.0
            x=coord[0]
            y=coord[1]
            z=coord[2]

            steady=True
            if (flag_mesh < 0 ):
                fvalue=cnst_value
            return fvalue, steady;
        sinks.append(sink)


    if (flag_grid == 'cube'):
        # domain
        coord=[];
        topol=[];

        # domain
        base=np.array([0,0,0])
        xl=1
        yl=1
        zl=1
        coord,topol=mt.MyParallelepiped(base,xl,yl,zl)

        #build mesh
        info = MeshInfo()
        info.set_points(coord)
        info.set_facets(topol)
        mesh = build(info, attributes=True, max_volume=(length**3/6.0))
        coord = np.array(mesh.points)


    if (flag_grid == 'balls'):
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
        
        points = int(1/length)
        #points = 100
        dphi = np.pi/points

        def truncate(radius):
            if abs(radius) < 1e-10:
                return 0
            else:
                return radius

        rz_l = [(truncate(radius_l*np.sin(i*dphi)), radius_l*np.cos(i*dphi)) for i in range(points+1)]
        rz_m = [(truncate(radius_m*np.sin(i*dphi)), radius_m*np.cos(i*dphi)) for i in range(points+1)]
        rz_s = [(truncate(radius_s*np.sin(i*dphi)), radius_s*np.cos(i*dphi)) for i in range(points+1)]

        
        geob = GeometryBuilder()
        geob.add_geometry(*generate_surface_of_revolution(rz_l,closure= EXT_OPEN, radial_subdiv=20))
        geob.add_geometry(*generate_surface_of_revolution(rz_s,closure= EXT_OPEN, radial_subdiv=20))
        geob.add_geometry(*generate_surface_of_revolution(rz_m,closure= EXT_OPEN, radial_subdiv=20))

        info = MeshInfo()
        geob.set(info)
        info.regions.resize(3)
        info.regions[0] = ([0.9,0.0,0.0] + [label_sink,   1]) 
        info.regions[1] = ([0.4,0.0,0.0] + [label_middle, 1])
        info.regions[2] = ([0.1,0.0,0.0] + [label_source, 1])

        mesh = build(info,attributes=True,max_volume=(length**3/6.0))
        
        coord=np.array(mesh.points)
        topol=np.array(mesh.elements)

        # cnst_value=1.0
        volume_source = (4.0/3.0)*np.pi*radius_s**3
        volume_sink   = (4.0/3.0)*np.pi*radius_l**3 - (4.0/3.0)*np.pi*radius_m**3
        cnst_value = 1.0 / volume_source
        cnst_sink_value = 1.0 / volume_sink
        cnst_sink_value = 1.2892
        # cnst_sink_value = (volume_source*cnst_value)/volume_sink
        c1 = cnst_value
        c2 = - cnst_sink_value
        print c1,c2

        
        def source(time,coord,flag_mesh):
            fvalue=0.0
            x=coord[0]
            y=coord[1]
            z=coord[2]

            steady=True
            if (flag_mesh>0):
                fvalue=cnst_value
            return fvalue, steady;
        sources.append(source)

        def sink(time,coord,flag_mesh):
            fvalue=0.0
            x=coord[0]
            y=coord[1]
            z=coord[2]

            steady=True
            if (flag_mesh<0):
                fvalue=cnst_sink_value
            return fvalue, steady;
        sinks.append(sink)
        
        # define optimal transport density
        # observe the center is in (0,0,0)
        def optdens(coord,flag_mesh):
            x=coord[0]
            y=coord[1]
            z=coord[2]
            value=0.0
            r_dist = np.sqrt(x**2+y**2+z**2)
            if (r_dist <= radius_s):
                value = - c1*0.5*r_dist
                value = np.abs(value)
                # if (r_dist > 0.30):
                #     print r_dist,value
            elif (r_dist <= radius_m):
                value =  - 1.0/r_dist*c1*0.5*(radius_s)**2
                value = np.abs(value)
                # if (r_dist > 0.64):
                #     print r_dist,value
            else:
                value = -1.0/r_dist*(0.5*c1*radius_s**2+0.5*c2*(r_dist**2-radius_m**2))
                value = np.abs(value)
                if (r_dist > 0.97):
                    print r_dist,value
            value = value**pflux
            return value;
        optdens_functions.append(optdens)
        

    if (flag_grid == 'sphere_plaplacian'):
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

        # to store middle points for each sphere
        middle_point_cache_l = {}
        middle_point_cache_m = {}
        middle_point_cache_s = {}

        # compute base coord and topol for each icosahedron
        coord_l,topol_l=mt.MyIcosahedron(center_domain,radius_l)
        coord_m,topol_m=mt.MyIcosahedron(center_domain,radius_m)
        coord_s,topol_s=mt.MyIcosahedron(center_domain,radius_s)

        # subdivide each sphere subdiv times
        subdiv = 4
        for i in range(subdiv):
            topol_l_subdiv = []
            for tri in topol_l:
                v1 = mt.middle_point(tri[0], tri[1], coord_l, radius_l, middle_point_cache_l)
                v2 = mt.middle_point(tri[1], tri[2], coord_l, radius_l, middle_point_cache_l)
                v3 = mt.middle_point(tri[2], tri[0], coord_l, radius_l, middle_point_cache_l)
                topol_l_subdiv.append([tri[0], v1, v3])
                topol_l_subdiv.append([tri[1], v2, v1])
                topol_l_subdiv.append([tri[2], v3, v2])
                topol_l_subdiv.append([v1, v2, v3])
            topol_l = topol_l_subdiv

            topol_m_subdiv = []
            for tri in topol_m:
                v1 = mt.middle_point(tri[0], tri[1], coord_m, radius_m, middle_point_cache_m)
                v2 = mt.middle_point(tri[1], tri[2], coord_m, radius_m, middle_point_cache_m)
                v3 = mt.middle_point(tri[2], tri[0], coord_m, radius_m, middle_point_cache_m)
                topol_m_subdiv.append([tri[0], v1, v3])
                topol_m_subdiv.append([tri[1], v2, v1])
                topol_m_subdiv.append([tri[2], v3, v2])
                topol_m_subdiv.append([v1, v2, v3])
            topol_m = topol_m_subdiv
            
            topol_s_subdiv = []
            for tri in topol_s:
                v1 = mt.middle_point(tri[0], tri[1], coord_s, radius_s, middle_point_cache_s)
                v2 = mt.middle_point(tri[1], tri[2], coord_s, radius_s, middle_point_cache_s)
                v3 = mt.middle_point(tri[2], tri[0], coord_s, radius_s, middle_point_cache_s)
                topol_s_subdiv.append([tri[0], v1, v3])
                topol_s_subdiv.append([tri[1], v2, v1])
                topol_s_subdiv.append([tri[2], v3, v2])
                topol_s_subdiv.append([v1, v2, v3])
            topol_s = topol_s_subdiv

        
        # join
        # I have no overlapping of volumes, surfaces and/or nodes,
        # so I just add each array within the others
        coord = coord_l + coord_m + coord_s
        topol = topol_l+(np.array(topol_m)+len(coord_l)).tolist()+(np.array(topol_s)+len(coord_l)+len(coord_m)).tolist()

        #set info mesh    
        info = MeshInfo()
        info.set_points(coord)
        info.set_facets(topol)
        info.regions.resize(3)
        info.regions[0] = ([0.8,0.0,0.0] + [label_sink,   1]) 
        info.regions[1] = ([0.5,0.0,0.0] + [label_middle, 1])
        info.regions[2] = ([0.1,0.0,0.0] + [label_source, 1])

        #build mesh
        mesh = build(info, attributes=True, max_volume=0.1*(length**3/6.0)/((1.2)**0))
        
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
        print 'Sink/Source - 1/19:',np.abs(c2/c1+1.0/19.0)
        
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



    return coord, topol, mesh, sources, sinks, kappas,tdens0_functions;
