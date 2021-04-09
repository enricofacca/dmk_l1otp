#!/usr/bin/env python
import re
import numpy as np
import meshtools as mt
import sys
import meshpy.triangle as make_triangle
import common 
import branch
import math
import os
import sys

# import timedata for IO operatorions
path_td= os.path.abspath(
    os.path.normpath(os.path.dirname(os.path.realpath(__file__))+
                     '/../../../globals/python_timedata'))
print (path_td)
print (os.getcwd())
print(os.path.dirname(os.path.realpath(__file__)))
sys.path.append(path_td)
import timedata as td


def example_grid(flag_grid,length,extra_info=None,minimum_angle=27.0):

    if (flag_grid == 'unitsquare'):
        # domain
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

    # Computational grid for example 1:
    # column with uniformly distributed load 
    # (traction and compression)
    if (flag_grid == 'example1_distributed_load'):
        # domain
        points,vertices=mt.MyRectangle([0,0],[0.3,2])

        #set info mesh    
        info = make_triangle.MeshInfo()
        info.set_points(points)
        info.set_facets(vertices)
        info.regions.resize(1)

        #build mesh
        mesh = make_triangle.build(info,
                                   attributes=True, 
                                   max_volume=(length**2/2.0))

    # Computational grid for example 1:
    # column with concentrated load 
    # (traction and compression)
    if (flag_grid == 'example1_concentrated_load'):
        # domain
        points,vertices=mt.MyRectangle([0,0],[0.3,2])

        coord = [[0.15,2.0]]
        coord=np.asarray(coord)
        list_coord=coord.tolist()
        points=points+list_coord

        #set info mesh    
        info = make_triangle.MeshInfo()
        info.set_points(points)
        info.set_facets(vertices)
        info.regions.resize(1)

        #build mesh
        mesh = make_triangle.build(info,
                                   attributes=True, 
                                   max_volume=(length**2/2.0))

    # Computational grid for example 1:
    # column with uniformly distributed load 
    # (traction and compression) in larger domain
    if (flag_grid == 'example1_distributed_load_inside'):
        # domain
        points,vertices=mt.MyRectangle([0,0],[1,3])

        p1,v1 = mt.MyRectangle([0.35,0.5],[0.65,2.5])

        points,vertices=mt.AddCurves(points,vertices,p1,v1)

        #set info mesh    
        info = make_triangle.MeshInfo()
        info.set_points(points)
        info.set_facets(vertices)
        info.regions.resize(1)

        #build mesh
        mesh = make_triangle.build(info,
                                   attributes=True, 
                                   max_volume=(length**2/2.0))

    # Computational grid for example 1:
    # column with concentrated load 
    # (traction and compression) in larger domain
    if (flag_grid == 'example1_concentrated_load_inside'):
        # domain
        points,vertices=mt.MyRectangle([0,0],[1,3])

        p1,v1 = mt.MyRectangle([0.35,0.5],[0.65,2.5])

        points,vertices=mt.AddCurves(points,vertices,p1,v1)

        coord = [[0.5,2.5]]
        coord=np.asarray(coord)
        list_coord=coord.tolist()
        points=points+list_coord

        #set info mesh    
        info = make_triangle.MeshInfo()
        info.set_points(points)
        info.set_facets(vertices)
        info.regions.resize(1)

        #build mesh
        mesh = make_triangle.build(info,
                                   attributes=True, 
                                   max_volume=(length**2/2.0))

    if (flag_grid == 'example1_distributed_load_semi_inside'):
        # domain
        points,vertices=mt.MyRectangle([0,0],[1,2])

        p1,v1 = mt.MyRectangle([0.35,0.0],[0.65,2.0])

        points,vertices=mt.AddCurves(points,vertices,p1,v1)

        #set info mesh    
        info = make_triangle.MeshInfo()
        info.set_points(points)
        info.set_facets(vertices)
        info.regions.resize(1)

        #build mesh
        mesh = make_triangle.build(info,
                                   attributes=True, 
                                   max_volume=(length**2/2.0))

    if (flag_grid == 'example2_distributed_load'):
        # domain
        points,vertices=mt.MyRectangle([0,0],[2,0.3])

        #set info mesh    
        info = make_triangle.MeshInfo()
        info.set_points(points)
        info.set_facets(vertices)
        info.regions.resize(1)

        #build mesh
        mesh = make_triangle.build(info,
                                   attributes=True, 
                                   max_volume=(length**2/2.0))

    if (flag_grid == 'example2_concentrated_load'):
        # domain
        points,vertices=mt.MyRectangle([0,0],[2,0.3])

        coord = [[2.0,0.15]]
        coord=np.asarray(coord)
        list_coord=coord.tolist()
        points=points+list_coord

        #set info mesh    
        info = make_triangle.MeshInfo()
        info.set_points(points)
        info.set_facets(vertices)
        info.regions.resize(1)

        #build mesh
        mesh = make_triangle.build(info,
                                   attributes=True, 
                                   max_volume=(length**2/2.0))

    if (flag_grid == 'example2_distributed_load_inside'):
        # domain
        points,vertices=mt.MyRectangle([0,0],[3,1])

        p1,v1 = mt.MyRectangle([0.5,0.35],[2.5,0.65])

        points,vertices=mt.AddCurves(points,vertices,p1,v1)

        #set info mesh    
        info = make_triangle.MeshInfo()
        info.set_points(points)
        info.set_facets(vertices)
        info.regions.resize(1)

        #build mesh
        mesh = make_triangle.build(info,
                                   attributes=True, 
                                   max_volume=(length**2/2.0))

    if (flag_grid == 'example2_concentrated_load_inside'):
        # domain
        points,vertices=mt.MyRectangle([0,0],[3,1])

        p1,v1 = mt.MyRectangle([0.5,0.35],[2.5,0.65])

        points,vertices=mt.AddCurves(points,vertices,p1,v1)

        coord = [[2.5,0.5]]
        coord=np.asarray(coord)
        list_coord=coord.tolist()
        points=points+list_coord

        #set info mesh    
        info = make_triangle.MeshInfo()
        info.set_points(points)
        info.set_facets(vertices)
        info.regions.resize(1)

        #build mesh
        mesh = make_triangle.build(info,
                                   attributes=True, 
                                   max_volume=(length**2/2.0))

        
    if (flag_grid == 'example2_distributed_load_semi_inside'):
        # domain
        points,vertices=mt.MyRectangle([0,0],[2,1.5])

        p1,v1 = mt.MyRectangle([0,0.6],[2,0.9])

        points,vertices=mt.AddCurves(points,vertices,p1,v1)

        #set info mesh    
        info = make_triangle.MeshInfo()
        info.set_points(points)
        info.set_facets(vertices)
        info.regions.resize(1)

        #build mesh
        mesh = make_triangle.build(info,
                                   attributes=True, 
                                   max_volume=(length**2/2.0))


    if (flag_grid == 'example3'):
        # domain
        points,vertices=mt.MyRectangle([0,0],[2.4,0.4])

        coord = [[1.2,0.4]]
        coord=np.asarray(coord)
        list_coord=coord.tolist()
        points=points+list_coord

        #set info mesh    
        info = make_triangle.MeshInfo()
        info.set_points(points)
        info.set_facets(vertices)
        info.regions.resize(1)

        #build mesh
        mesh = make_triangle.build(info,
                                   attributes=True, 
                                   max_volume=(length**2/2.0))

    if (flag_grid == 'example3_inside'):
        # domain
        points,vertices=mt.MyRectangle([0,0],[3,2])

        p1,v1 = mt.MyRectangle([0.3,0.8],[2.7,1.2])

        points,vertices=mt.AddCurves(points,vertices,p1,v1)

        coord = [[1.5,1.2]]
        coord=np.asarray(coord)
        list_coord=coord.tolist()
        points=points+list_coord

        #set info mesh    
        info = make_triangle.MeshInfo()
        info.set_points(points)
        info.set_facets(vertices)
        info.regions.resize(1)

        #build mesh
        mesh = make_triangle.build(info,
                                   attributes=True, 
                                   max_volume=(length**2/2.0))

    if (flag_grid == 'example3_semi_inside'):
        # domain
        points,vertices=mt.MyRectangle([0,0],[2.4,1.5])

        p1,v1 = mt.MyRectangle([0.0,0.0],[2.4,0.4])

        points,vertices=mt.AddCurves(points,vertices,p1,v1)

        coord = [[1.2,0.4]]
        coord=np.asarray(coord)
        list_coord=coord.tolist()
        points=points+list_coord

        #set info mesh    
        info = make_triangle.MeshInfo()
        info.set_points(points)
        info.set_facets(vertices)
        info.regions.resize(1)

        #build mesh
        mesh = make_triangle.build(info,
                                   attributes=True, 
                                   max_volume=(length**2/2.0))

    if (flag_grid == 'example_book'):
        # domain
        points,vertices=mt.MyRectangle([0,0],[1,1])

        coord = [[0.0,0.5],[0.5,0.5],[1.0,0.5]]
        coord=np.asarray(coord)
        list_coord=coord.tolist()
        points=points+list_coord

        #set info mesh    
        info = make_triangle.MeshInfo()
        info.set_points(points)
        info.set_facets(vertices)
        info.regions.resize(1)

        #build mesh
        mesh = make_triangle.build(info,
                                   attributes=True, 
                                   max_volume=(length**2/2.0))

    if (flag_grid == 'example_internet1'):
        # domain
        points,vertices=mt.MyRectangle([0,0],[1.5,0.75])

        #set info mesh    
        info = make_triangle.MeshInfo()
        info.set_points(points)
        info.set_facets(vertices)
        info.regions.resize(1)

        #build mesh
        mesh = make_triangle.build(info,
                                   attributes=True, 
                                   max_volume=(length**2/2.0))

    if (flag_grid == 'example_internet1_semi_inside'):
        # domain
        points,vertices=mt.MyRectangle([0,0],[1.75,1.0])

        coord = [[1.5,0.0]]
        coord=np.asarray(coord)
        list_coord=coord.tolist()
        points=points+list_coord

        #set info mesh    
        info = make_triangle.MeshInfo()
        info.set_points(points)
        info.set_facets(vertices)
        info.regions.resize(1)

        #build mesh
        mesh = make_triangle.build(info,
                                   attributes=True, 
                                   max_volume=(length**2/2.0))

    if (flag_grid == 'example_internet2'):
        # domain
        points,vertices=mt.MyRectangle([0,0],[1.5,0.75])

        coord = [[1.5,0.375]]
        coord=np.asarray(coord)
        list_coord=coord.tolist()
        points=points+list_coord

        #set info mesh    
        info = make_triangle.MeshInfo()
        info.set_points(points)
        info.set_facets(vertices)
        info.regions.resize(1)

        #build mesh
        mesh = make_triangle.build(info,
                                   attributes=True, 
                                   max_volume=(length**2/2.0))

    if (flag_grid == 'example_internet2_semi_inside'):
        # domain
        points,vertices=mt.MyRectangle([0,0],[1.75,1.25])

        p1,v1 = mt.MyRectangle([0.0,0.25],[1.5,1.0])

        points,vertices=mt.AddCurves(points,vertices,p1,v1)

        coord = [[1.5,0.625]]
        coord=np.asarray(coord)
        list_coord=coord.tolist()
        points=points+list_coord

        #set info mesh    
        info = make_triangle.MeshInfo()
        info.set_points(points)
        info.set_facets(vertices)
        info.regions.resize(1)

        #build mesh
        mesh = make_triangle.build(info,
                                   attributes=True, 
                                   max_volume=(length**2/2.0))

    if (flag_grid == 'regular'):
        # domain        
        input_file=open(extra_info,'r')
        input_lines = input_file.readlines()
        x_ndiv = int(input_lines[0].split()[0])
        y_ndiv = int(input_lines[0].split()[1])
        
        x_height = float(input_lines[1].split()[0])
        y_height = float(input_lines[1].split()[1])

        x_origin = float(input_lines[2].split()[0])
        y_origin = float(input_lines[2].split()[1])
        input_file.close()
        
        #print(x_ndiv,y_ndiv,x_height,y_height,x_origin,y_origin) 

        coordinates, elements=rectangule_grid(x_ndiv,y_ndiv,
                x_height,y_height,x_origin,y_origin)
        points=[
            [x_origin,y_origin,0],
            [x_origin+x_height,y_origin,0],
            [x_origin+x_height,y_origin+y_height,0],
            [x_origin,y_origin+y_height,0],
        ]
        vertices=[
            [0,1],
            [1,2],
            [2,3],
            [3,0]
        ]
        elements_attributes=range(len(elements))

    if (flag_grid == 'fix_points'):
        # domain
        points,vertices=mt.MyRectangle([0,0],[1,1])

        # import nodes coordinates
        coord = td.read_steady_timedata(extra_info)
        #print(coord)
        list_coord=coord[:,0:2].tolist()
        points=points+list_coord
        #points.append(coord[
        # add the point
        
        


        #set info mesh    
        info = make_triangle.MeshInfo()
        info.set_points(points)
        info.set_facets(vertices)
        info.regions.resize(1)

        #build mesh
        mesh = make_triangle.build(info,
                                   attributes=True, 
                                   max_volume=(length**2/2.0))

    if (flag_grid == 'three_forces_0'):
        # domain
        points,vertices=mt.MyRectangle([0,0],[1,1])

        print(np.pi/2+2*np.pi/3)
        # set point at vertices of equilateral triangles
        coord = [
            [np.cos(np.pi/2          ),np.sin(np.pi/2          )],
            [np.cos(np.pi/2+2*np.pi/3),np.sin(np.pi/2+2*np.pi/3)],
            [np.cos(np.pi/2+4*np.pi/3),np.sin(np.pi/2+4*np.pi/3)],
        ]
        coord=np.asarray(coord)
        coord[:,:]=0.3*coord[:,:]
        coord[:,:]=coord+[0.5,0.5]
        
        print(coord)
        list_coord=coord[:,0:2].tolist()
        print(list_coord)
        points=points+list_coord

        #set info mesh    
        info = make_triangle.MeshInfo()
        info.set_points(points)
        info.set_facets(vertices)
        info.regions.resize(1)

        #build mesh
        mesh = make_triangle.build(info,
                                   attributes=True, 
                                   max_volume=(length**2/2.0))

    if (flag_grid == 'three_forces_1to6'):
        # domain
        points,vertices=mt.MyRectangle([0,0],[1,1])

        # set point at vertices of equilateral triangles
        coord = [
            [np.cos(3*np.pi/4),np.sin(3*np.pi/4)],
            [np.cos(5*np.pi/4),np.sin(5*np.pi/4)],
            [np.cos(7*np.pi/4),np.sin(7*np.pi/4)],
        ]
        coord=np.asarray(coord)
        coord[:,:]=0.3*coord[:,:]
        coord[:,:]=coord+[0.5,0.5]

        print(coord)
        list_coord=coord[:,0:2].tolist()
        points=points+list_coord

        #set info mesh    
        info = make_triangle.MeshInfo()
        info.set_points(points)
        info.set_facets(vertices)
        info.regions.resize(1)

        #build mesh
        mesh = make_triangle.build(info,
                                   attributes=True, 
                                   max_volume=(length**2/2.0),
                                   min_angle=28)

        
    if (flag_grid == 'hole_square'):
        # domain
        
        


        
        source=[0,0]
      
        # input_file=open(extra_info,'r')
        # input_lines = input_file.readlines()
        # x_ndiv = int(input_lines[0].split()[0])
        # y_ndiv = int(input_lines[0].split()[1])
        
        # x_height = float(input_lines[1].split()[0])
        # y_height = float(input_lines[1].split()[1])

        # x_origin = float(input_lines[2].split()[0])
        # y_origin = float(input_lines[2].split()[1])

        # box = float(input_lines[3])
        # subdivision_circle = float(input_lines[3])

        # input_file.close()
        
        # print(x_ndiv,y_ndiv,x_height,y_height,x_origin,y_origin)

        box=0.1
        points_hole, vertices_hole = mt.CircleSegments(source,0.1,edge_length=length/30.0)


        #points_hole,vertices_hole=mt.MyRectangle([source[0]-box,source[1]-box],[source[0]+box,source[1]+box])
        points,vertices=mt.MyRectangle([-2,-2],[2,2])
        markers = [1,1,1,1]
        points,vertices=mt.AddCurves(points,vertices,points_hole,vertices_hole)
        for i in range(len(vertices_hole)):
            markers.extend([2])

        #print markers
        

        #print points
        #print vertices

        #circ_start = len(points)
        #points.extend(
        #    (3 * np.cos(angle), 3 * np.sin(angle))
        #for angle in np.linspace(0, 2*np.pi, 29, endpoint=False))

        #vertices.extend(round_trip_connect(circ_start, len(points)-1))

        def needs_refinement(vertices, area):
            bary = np.sum(np.array(vertices), axis=0)/3
            #print bary
            max_area = 0.1*length**2 + abs((np.linalg.norm(bary)))*length
            return bool(area > max_area)

        info = make_triangle.MeshInfo()
        info.set_points(points)
        info.set_holes([(0,0)])
        info.set_facets(vertices,markers)

        
        mesh = make_triangle.build(info, 
                                   attributes=True, 
                                   #max_volume=length**2)
                                   refinement_func=needs_refinement)

        
        all_edges,boundary_edges = mt.FindEdges(np.int_(mesh.elements))

        coord=np.array(mesh.points)

        flag_nodes=np.zeros(len(coord),dtype=int)
        for i in range(len(boundary_edges)):
            #print boundary_edges[i]
            for node in boundary_edges[:][i]:
                if ( abs(np.linalg.norm(coord[:][node]-source)-box) < 1e-10 ):
                    flag_nodes[node] = 1
                    #print (coord[:][node])
        

        #mesh = triangle.build(info) 
        
        #print len(mesh.element_attributes)
        #for i in range(len(mesh.element_attributes)):
        #    print i, mesh.elements_attributes[i]


    if (flag_grid == 'segala'):
        # domain
        #length=length*4
        points,vertices=mt.MyRectangle([-1,-1],[1,1])

        r=0.75
        points.append([0 ,  0])
        points.append([r*np.cos(np.pi/6),  r*np.sin(np.pi/6)])
        points.append([r*np.cos(np.pi*5/6),r*np.sin(np.pi*5/6)])
        points.append([r*np.cos(-np.pi/2), r*np.sin(-np.pi/2)])

        #set info mesh    
        info = make_triangle.MeshInfo()
        info.set_points(points)
        info.set_facets(vertices)
        info.regions.resize(1)

        #build mesh
        mesh = make_triangle.build(info,
                                   attributes=True, 
                                   max_volume=(length**2/2.0))

        


    if (flag_grid == 'plaplacian'):
        # grid use for testing convergence for
        # plaplacian equation in  
        # Branching Structures Emerging From A
        #Continuous Optimal Transport Model
        # Enrico Facca, Franco Cardin, Mario Putti
        
        # domain
        c=[0.0,0.0]
        r=1.0
        points, vertices = mt.CircleSegments(c,r,edge_length=length)

        c=[0.0,0.0]
        r=1.0/3.0
        p1, v1 = mt.CircleSegments(c,r,edge_length=length)

        c=[0.0,0.0]
        r=2.0/3.0
        p2, v2 = mt.CircleSegments(c,r,edge_length=length)

        points,vertices=mt.AddCurves(points,vertices,p1,v1)
        points,vertices=mt.AddCurves(points,vertices,p2,v2)

        #build mesh
        info = make_triangle.MeshInfo()
        info.set_points(points)
        info.set_facets(vertices)
        info.regions.resize(1)
        mesh = make_triangle.build(info,
                                   attributes=True, 
                                   max_volume=(length**2/2.0))

    
        
        
    # An accurate discontinuous Galerkin method for solving point-source
    # Eikonal equation in 2-D heterogeneous anisotropic media
    # Bouteiller,1 M. Benjemaa,2 L. Metivier1,3 and J. Virieux1

    if (flag_grid == 'bouteiller18_tc1'):
        # domain
        points,vertices=mt.MyRectangle([0,0],[4,4])
        
        points_source=([2.0,2.0])
        c1=[2.0,2.0]
        r1=length
        points_source, vertices_source = mt.CircleSegments(c1,r1,edge_length=length/30)
        
        points_source.append([2.0,2.0])
        nportion=len(vertices_source)
        for i in range(nportion):
            vertices_source.append((i,nportion))
        
        
        points,vertices=mt.AddCurves(points,vertices,points_source,vertices_source)
        

        
        #set info mesh    
        info = make_triangle.MeshInfo()
        info.set_points(points)
        info.set_facets(vertices)
        info.regions.resize(1)

        def needs_refinement(vertices, area ):
            vert_origin, vert_destination, vert_apex = vertices
            bary_x = (vert_origin.x + vert_destination.x + vert_apex.x) / 3
            bary_y = (vert_origin.y + vert_destination.y + vert_apex.y) / 3
            
            dist_center = math.sqrt( (bary_x-2)**2 + (bary_y-2)**2 )
            max_area = (dist_center) * length**2/2.0
            return area > max_area
        
        #build mesh
        mesh=make_triangle.build(info,
                                 attributes=True,
                                 max_volume=(length**2/2.0))#,refinement_func=needs_refinement)
    

        

            

    if (flag_grid == 'bouteiller18_tc2'):
        # domain
        points,vertices=mt.MyRectangle([0,0],[4,4])
        
        points.append([2.0,2.0])

        c1=[1.0,1.5]
        r1=0.5
        p1, v1 = mt.CircleSegments(c1,r1,edge_length=length)
        
        points,vertices=mt.AddCurves(points,vertices,p1,v1)



        
        #set info mesh    
        info = make_triangle.MeshInfo()
        info.set_points(points)
        info.set_facets(vertices)
        info.regions.resize(1)

        #build mesh
        mesh = make_triangle.build(info,
                                   attributes=True, 
                                   max_volume=(length**2/2.0))
        

    #  A parallel fast sweeping method for the Eikonal equation
    #  Miles Detrixhea, Frederic Giboua,b, Chohong Min
    if (flag_grid == 'example1fmt'):
        # domain
        p1,v1=mt.MyRectangle([-1,-1],[0,0])
        p2,v2=mt.MyRectangle([-1,1],[0,0])
        p3,v3=mt.MyRectangle([1,-1],[0,0])
        p4,v4=mt.MyRectangle([1,1],[0,0])

        points,vertices=mt.AddCurves(p1,v1,p2,v2)
        points,vertices=mt.AddCurves(points,vertices,p3,v3)
        points,vertices=mt.AddCurves(points,vertices,p4,v4)

        
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
    
    if (flag_grid == 'cnst_center'):
        #domain
        points,vertices=mt.MyRectangle([0,0],[1,1])
        
        # center square
        p1,v1=mt.MyRectangle([0.4,0.4],[0.6,0.6])
        

        # union   
        points,vertices=mt.AddCurves(points,vertices,p1,v1)
        
    

        #set info mesh    
        info = make_triangle.MeshInfo()
    
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
   
    if (flag_grid == 'prigozhin'):
        # domain
        c0=[0.1,0.1]
        l0=0
        points,vertices=mt.MyRectangle([0,0],[1,1])
        
        #source
        c1=[0.2,0.75]
        r1=0.15
        l1=1
        p1, v1 = mt.CircleSegments(c1,r1,edge_length=length)
        
        #sink
        c2=[0.85,0.35]
        r2=0.3
        l2=-1
        p2, v2 = mt.CircleSegments(c2,r2,edge_length=length)
        for point in p2:
            point[0] = c2[0] + ( point[0] - c2[0] ) / 3.0

        #kappas
        c3=[0.0,0.0]
        r3=0.4
        l3=2
        p3, v3 = mt.CircleSegments(c3,r3,edge_length=length)
        for point in p3:
            point[0] = point[0]/ 5.0
        theta=np.pi/4.0
        for point in p3:
            x=point[0]
            y=point[1]
            point[0] = np.cos(theta)*x + np.sin(theta)*y
            point[1] = -np.sin(theta)*x+ np.cos(theta)*y
            point[0] = point[0] + 0.5 
            point[1] = point[1] + 0.5
        c3=[0.5,0.5]


        #union
        points,vertices=mt.AddCurves(points,vertices,p1,v1)
        points,vertices=mt.AddCurves(points,vertices,p2,v2)
        points,vertices=mt.AddCurves(points,vertices,p3,v3)
        
        #set info mesh
        info = make_triangle.MeshInfo()
        info.set_points(points)
        info.set_facets(vertices)
        info.regions.resize(4)
        info.regions[0] = (c0 + [l0, 1])
        info.regions[1] = (c1 + [l1, 1]) 
        info.regions[2] = (c2 + [l2, 1]) 
        info.regions[3] = (c3 + [l3, 1]) 

    
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
        # def source(t,x,y,z,flag_mesh):
        #     fvalue=0.0
        #     steady=True
        #     if (flag_mesh == l1) : 
        #         fvalue=max(0,(r1-np.sqrt((x-c1[0])**2+(y-c1[1])**2))/r1)
        #     return fvalue,steady;
        # sources.append(source)
        
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
        #sinks.append(sink0)

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
        #sinks.append(sink1)


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
        #sinks.append(sink2)

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
        #sinks.append(sink3)


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
        # def source(t,x,y,z,flag_mesh):
        #     fvalue=0.0
        #     steady=True
        #     if (flag_mesh == l1) : 
        #         fvalue=max(0,(r1/2.0-max(abs(x-c1[0]),abs(y-c1[1])))/(r1/2.0))
        #     return fvalue,steady;
        # sources.append(source)
        
        #sink 
        # sw
        c2=[0.35,0.15]
        r2=0.1
        l2=-1
        base2_sw=[c2[0]-r2/2.0,c2[1]-r2/2.0]
        base2_ne=[c2[0]+r2/2.0,c2[1]+r2/2.0]
        p2, v2 = mt.MyRectangle(base2_sw,base2_ne)
        # def sink0(t,x,y,z,flag_mesh):
        #     fvalue=0.0
        #     steady=True
        #     if (flag_mesh == l2):
        #         fvalue=max(0,(r2/2.0-max(abs(x-c2[0]),abs(y-c2[1])))/(r2/2.0))
        #     return fvalue,steady;
        # sinks.append(sink0)

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


    # convert to numpy arrays 
    coordinates=np.array(coordinates)
    if ( coordinates.shape[1] == 2 ):
        zcoord      = np.zeros([coordinates.shape[0],1])
        coordinates = np.append(coordinates, zcoord, axis=1)
        elements    = np.array(elements)
    try:
        print (len(elements_attributes))
        elements_attributes=np.array(elements_attributes)
    except NameError:
        elements_attributes=np.int_(range(len(topol)))
        

    return points, vertices, coordinates, elements, elements_attributes;





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

def rectangule_grid(ndivx,ndivy,x_height,y_height,x_origin,y_origin):    
    # grid data

    ndivx=int(ndivx)
    ndivy=int(ndivy)
    
    ntria=2*ndivx*ndivy
    nnode=(ndivx+1)*(ndivy+1)
    coord=np.zeros([nnode,3])
    triang=np.zeros([ntria,3],dtype=int)
        
    # aux. data
    ntriax=2*ndivx
    nnodex=ndivx+1
    lenx = x_height / ndivx
    leny = y_height / ndivy

    inode=0
    itria=0
    for iy in range(ndivy): # 0,1, ndivy-1
        for ix in range(ndivx): # 0,1, ndivx-1
            coord[inode,:] = (ix*lenx,iy*leny,0)
            inode = inode+1
            # ------
            # |    |
            # |\   |
            # | \  |
            # |  \ |
            # |   \|
            # ------ 
            # sw triangle 
            n1= iy * nnodex + ix 
            n2= n1 + 1
            n3= n1 + nnodex 
            triang[itria,:]=(n1,n2,n3)
            
            #print 'itria=',itria, 'n=', n1,n2,n3
            

            itria=itria+1
            # copy
            n1old=n1
            n2old=n2
            n3old=n3

            # ne triangle 
            n1= n2old
            n2= n3old + 1
            n3= n3old

            #print 'itria=',itria, 'n=', n1,n2,n3

            triang[itria,:]=(n1,n2,n3)

            itria=itria+1

        # add last point in x direction
        coord[inode,:] = ((ix+1)*lenx,iy*leny,0.0)
        inode = inode+1
   

    # add line of top points 
    for ix in range(ndivx+1): # 0,1, ndivx
        coord[inode,:] = (ix*lenx,(iy+1)*leny,0)
        inode = inode+1

    coord[:,1]=-coord[:,1]+ndivy*leny

    coord[1,:] = coord[1,:] + x_origin
    coord[2,:] = coord[2,:] + y_origin
    
    return coord, triang;

