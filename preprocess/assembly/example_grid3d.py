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
    if (flag_grid == 'rect_cnst'):
        # domain
        points=[];
        faces=[];


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
        points,faces=mt.MyParallelepiped(base,xl,yl,zl)

        # source
        base=base_optdens
        xl=xwidth
        yl=ywidth
        zl=zwidth
        psource = base + np.array([xl,yl,zl])/2.0
        label_source = 1
        points1,faces1=mt.MyParallelepiped(base,xl,yl,zl)

        # sink
        base=base_optdens+np.array([shift,0,0])
        xl=xwidth
        yl=ywidth
        zl=zwidth
        psink= base + np.array([xl,yl,zl])/2.0
        label_sink=-1
        points2,faces2=mt.MyParallelepiped(base,xl,yl,zl)
        
        # join
        points,faces=mt.Add3dSurface(points,faces,points1,faces1)
        points,faces=mt.Add3dSurface(points,faces,points2,faces2)


        #set info mesh    
        info = MeshInfo()
        info.set_points(points)
        info.set_facets(faces)
        info.regions.resize(3)
        info.regions[0] = ([0.5,0.5,0.5] + [0, 1])
        info.regions[1] = (list(psource) + [label_source, 1])
        info.regions[2] = (list(psink) + [label_sink,   1]) 

        #build mesh
        mesh = build(info, attributes=True, max_volume=(length**3/6.0))

    if (flag_grid == 'column'):
        # domain
        #points,faces=mt.MyParallelepiped([0,0,0],1,1,1)
        
        
        base_center=[0.5,0.5]
        radius=0.1
        height=3
        points, faces = mt.Cilinder(base_center, radius, height, length)
        
        info = MeshInfo()
        info.set_points(list(points))
        info.set_facets(faces)

        #build mesh
        mesh = build(info, attributes=True, max_volume=(length**3/6.0))
        
    if (flag_grid == 'rect_cnst_aligned'):
        # domain
        points=[];
        faces=[];


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
        points,faces=mt.MyParallelepiped(base,xl,yl,zl)

        # source
        base=base_optdens
        xl=xwidth
        yl=ywidth
        zl=zwidth
        psource = base + np.array([xl,yl,zl])/2.0
        label_source = 1
        points1,faces1=mt.MyParallelepiped(base,xl,yl,zl)

        # sink
        base=base_optdens+np.array([shift,0,0])
        xl=xwidth
        yl=ywidth
        zl=zwidth
        psink= base + np.array([xl,yl,zl])/2.0
        label_sink=-1
        points2,faces2=mt.MyParallelepiped(base,xl,yl,zl)

        # middle
        base=base_optdens+np.array([xwidth,0,0])
        xl=shift-xwidth
        yl=ywidth
        zl=zwidth
        pmiddle= base + np.array([xl,yl,zl])/2.0
        label_middle=0
        points3,faces3=mt.MyParallelepiped(base,xl,yl,zl)

        
        # join
        points,faces=mt.Add3dSurface(points,faces,points1,faces1)
        points,faces=mt.Add3dSurface(points,faces,points3,faces3)
        points,faces=mt.Add3dSurface(points,faces,points2,faces2)
       
            
        #set info mesh    
        info = MeshInfo()
        info.set_points(points)
        info.set_facets(faces)
        info.regions.resize(4)
        info.regions[0] = ([0.005,0.005,0.005] + [0, 1])
        info.regions[1] = (list(psource) + [label_source, 1])
        info.regions[2] = (list(psink) + [label_sink,   1]) 
        info.regions[3] = (list(pmiddle) + [label_middle,   1])
        #build mesh
        mesh = build(info, attributes=True, max_volume=(length**3/6.0))

    if (flag_grid == 'cube'):
        # domain
        points=[];
        faces=[];

        # domain
        base=np.array([0,0,0])
        xl=1
        yl=1
        zl=1
        points,faces=mt.MyParallelepiped(base,xl,yl,zl)

        #build mesh
        info = MeshInfo()
        info.set_points(points)
        info.set_facets(faces)
        mesh = build(info, attributes=True, max_volume=(length**3/6.0))
        points = np.array(mesh.points)


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
        
        points=np.array(mesh.points)
        faces=np.array(mesh.elements)
 

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

        # compute base points and faces for each icosahedron
        points_l,faces_l=mt.MyIcosahedron(center_domain,radius_l)
        points_m,faces_m=mt.MyIcosahedron(center_domain,radius_m)
        points_s,faces_s=mt.MyIcosahedron(center_domain,radius_s)

        # subdivide each sphere subdiv times
        subdiv = 4
        for i in range(subdiv):
            faces_l_subdiv = []
            for tri in faces_l:
                v1 = mt.middle_point(tri[0], tri[1], points_l, radius_l, middle_point_cache_l)
                v2 = mt.middle_point(tri[1], tri[2], points_l, radius_l, middle_point_cache_l)
                v3 = mt.middle_point(tri[2], tri[0], points_l, radius_l, middle_point_cache_l)
                faces_l_subdiv.append([tri[0], v1, v3])
                faces_l_subdiv.append([tri[1], v2, v1])
                faces_l_subdiv.append([tri[2], v3, v2])
                faces_l_subdiv.append([v1, v2, v3])
            faces_l = faces_l_subdiv

            faces_m_subdiv = []
            for tri in faces_m:
                v1 = mt.middle_point(tri[0], tri[1], points_m, radius_m, middle_point_cache_m)
                v2 = mt.middle_point(tri[1], tri[2], points_m, radius_m, middle_point_cache_m)
                v3 = mt.middle_point(tri[2], tri[0], points_m, radius_m, middle_point_cache_m)
                faces_m_subdiv.append([tri[0], v1, v3])
                faces_m_subdiv.append([tri[1], v2, v1])
                faces_m_subdiv.append([tri[2], v3, v2])
                faces_m_subdiv.append([v1, v2, v3])
            faces_m = faces_m_subdiv
            
            faces_s_subdiv = []
            for tri in faces_s:
                v1 = mt.middle_point(tri[0], tri[1], points_s, radius_s, middle_point_cache_s)
                v2 = mt.middle_point(tri[1], tri[2], points_s, radius_s, middle_point_cache_s)
                v3 = mt.middle_point(tri[2], tri[0], points_s, radius_s, middle_point_cache_s)
                faces_s_subdiv.append([tri[0], v1, v3])
                faces_s_subdiv.append([tri[1], v2, v1])
                faces_s_subdiv.append([tri[2], v3, v2])
                faces_s_subdiv.append([v1, v2, v3])
            faces_s = faces_s_subdiv

        
        # join
        # I have no overlapping of volumes, surfaces and/or nodes,
        # so I just add each array within the others
        points = points_l + points_m + points_s
        faces = faces_l+(np.array(faces_m)+ len(points_l)).tolist()+(np.array(faces_s)+len(points_l)+len(points_m)).tolist()

        #set info mesh    
        info = MeshInfo()
        info.set_points(points)
        info.set_facets(faces)
        info.regions.resize(3)
        info.regions[0] = ([0.8,0.0,0.0] + [label_sink,   1]) 
        info.regions[1] = ([0.5,0.0,0.0] + [label_middle, 1])
        info.regions[2] = ([0.1,0.0,0.0] + [label_source, 1])

        #build mesh
        mesh = build(info, attributes=True, max_volume=0.1*(length**3/6.0)/((1.2)**0))
        
    try:
        coordinates 
    except NameError:
        coordinates=mesh.points
    try:
        elements
    except NameError:
        elements=mesh.elements    
    try:
        elements_attributes
    except NameError:
        elements_attributes=mesh.element_attributes
    
    return points, faces, coordinates, elements, elements_attributes;#, nodes_attributes;

