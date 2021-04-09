# -*- coding: utf-8 -*-
"""
Toolbox for generating a mesh

"""
import numpy as np
import scipy as sp
from scipy.spatial import cKDTree
import matplotlib.pyplot as plt
import meshpy.triangle as triangle
from pyvtk import *
from meshpy.tet import MeshInfo, build


# Extract the edges
# ouput, edges and boundary edges 
def FindEdges(t):
  #pdb.set_trace();  
  NE=t.shape[0]
  # generate an array of all edges
  tt=np.array([t[:,0],t[:,1],t[:,1],t[:,2],t[:,2],t[:,0]]).T.reshape(3*NE,2)
  ttt=np.sort(tt,1)
  
  # find all boundary edges
  all_edges=[ tuple(x) for x in ttt ]
  boundary_edges=[x for x in all_edges if all_edges.count(x)==1]
  
  # find all unique edges
  all_edges=list(set(all_edges))
  return all_edges,boundary_edges;





##################
#
#  Boundary Tools
#
##################

# given one segment 
# e.g.  (X,2) find segment (2,Y) and delete (2,Y) from list 
def FindNextSegment(all_segments,node):
  # find next connecting segment  
  help=[x for x in all_segments if x[0]==node]   
  
  new_bound=False
  if len(help)==0: #if connecting segment does not exist (=>new boundary) 
    ret=all_segments[0]
    new_bound=True    
  else:
    ret=help[0]
  
  del all_segments[all_segments.index(ret)]
  return ret,new_bound;  
  
  
# sort segments:  (3,6),(6,1),(1,12),(12,5),...
# on output: sorted segments and indices of the different boundaries
def SortSegments(all_segments):  
  count=len(all_segments)

  node=-1
  sorted_segments=[]
  boundaries=[]
  for j in range(len(all_segments)):
    seg,new_bound=FindNextSegment(all_segments,node)
    node=seg[1]
    sorted_segments.append(seg)    
    if new_bound==True:
      boundaries.append(j)
    
  if len(sorted_segments)!=count:
    print("Something is wrong, number of segments not the same")   
  return sorted_segments,boundaries;

# connect segments in a defined way
# (see SortSegments), but start sorting with a defined point p
# multiple p'2 for different closed boundaries are possible
def ConnectBoundary(boundary_segments,Pall,p=[]):
  
  # sort the boundary segments  
  allseg=boundary_segments[:]  
  allseg,boundaries=SortSegments(allseg)
  if p==[]:
    return allseg,boundaries;
    
  max_boundaries=len(boundaries)
   
  # find all nodes on the given boundary
  nodes=[x[0] for x in allseg]
  # find closest nodes to desired point list p  
  indices,distances=FindClosestNode(nodes,Pall,p)
  
  #change order within each closed boundary
  flag_sorted=[]
  for j in range(len(boundaries)):
   flag_sorted.append(False) 
   
  for j in range(len(indices)):
    # find position of node in the boundary list
    # indj gives the position of the segment in allseg
    indj=nodes.index(indices[j])
    # find the number of boundary the node belongs to
    this_boundary=(np.where((np.array(boundaries)<=indj))[0])[-1]
    
    if flag_sorted[this_boundary]==False:
      # define the indices for slicing      
      ind_1=boundaries[this_boundary]
      if this_boundary+1==max_boundaries:
        ind_2=len(allseg)
      else:
        ind_2=boundaries[this_boundary+1]  
      
      # rearange the segments in the corresponding boundary     
      allseg=allseg[:ind_1]+allseg[indj:ind_2]+allseg[ind_1:indj]+allseg[ind_2:]
      # resort only once      
      flag_sorted[this_boundary]=True
  
  return allseg,boundaries;


#
# find closest node to point p0 in a list of N nodes
# Pall coordinates of M nodes  M>=N
# constraint defines constraints on distance
def FindClosestNode(nodes,Pall,p0,constraint=-1,tree=None):
  # take those points of the node list
  
  if tree==None:
    p_nodes=np.array(Pall)
    p_nodes=p_nodes[nodes] 
    # look for minimum distance, define dist function
    mytree = cKDTree(p_nodes)
  else:
    mytree=tree
    
  dist, index = mytree.query(np.array(p0))
  
  node_closest=[nodes[j] for j in index]
   
  # check constraints
  num_p= len(p0)
  if constraint<0:
    return node_closest,dist;
  elif np.isscalar(constraint)==True:
    constraint=constraint*np.ones(num_p)
  elif len(p0)!=len(constraint):
    print('Error in constraint definition')
    return [],[]
  
  # check constraint for each node
  flags=[((dist[j]<=constraint[j]) | (constraint[j]<0)) for j in range(num_p)]
  for j in range(num_p):
    if flags[j]==False:
      node_closest[j]=-1
  return node_closest,dist;
  
   
# check relative position of two points   
def SamePoint(p1,p2,delta):
  dp=(np.array(p1)-np.array(p2))
  d=np.sqrt(dp[0]**2+dp[1]**2+dp[2]**2)
  ret=False  
  if d<delta:
    ret=True
  return ret;

# test whether tuples have the same content
def cmpT(t1,t2):
  if t1 == t2:
    return True
  else:
    return False



#####################
#
# Make simple curves
#
#####################
#
#
# 
# make a circle or part of it  
#
def CircleSegments(middle,radius,num_points=10,a_min=0.,a_max=2.*np.pi,edge_length=-1):  
  # check for closed loop
  number_points=num_points
  if edge_length>0:
    number_points=np.floor(abs(radius/edge_length*(a_max-a_min)))+1
    
  delta=(a_max-a_min)/number_points  
  closed=False;  
  if abs(a_max-a_min-2*np.pi)<0.1*delta:
    closed=True
    
  t=np.linspace(a_min,a_max,number_points,not closed)
  # define points
  points=[(middle[0]+radius*np.cos(angle),middle[1]+radius*np.sin(angle)) for angle in t]
  
  # define vertices
  vertices=[(j,j+1) for j in range(0,len(points)-1,1)]    
  if closed==True:
    vertices+=[(len(points)-1,0)]
  return points,vertices;



# Straight line
def LineSegments(P1,P2,num_points=10,edge_length=-1):
  number_points=num_points
  if edge_length>0:
    p1=np.array(P1)
    p2=np.array(P2)
    number_points=np.floor(np.sqrt(np.sum((p2-p1)**2))/edge_length)+1
  
  t=np.linspace(0,1,number_points)
  points=[(P1[0]+param*(P2[0]-P1[0]),P1[1]+param*(P2[1]-P1[1])) for param in t]
  vertices=[(j,j+1) for j in range(0,len(points)-1,1)]
  return points,vertices;

# Straight line
def MyLineSegments(P1,P2,edge_length=-1):
  p1=np.array(P1)
  p2=np.array(P2)
  if edge_length>0:
    number_points=np.floor(np.sqrt(np.sum((p2-p1)**2))/edge_length)+1
  else:
    number_points=2
  
  t=np.linspace(0,1,number_points)
  points=[[P1[0]+param*(P2[0]-P1[0]),P1[1]+param*(P2[1]-P1[1])] for param in t]
  vertices=[(j,j+1) for j in range(0,len(points)-1,1)]
  return points,vertices;



# Rectangle
def RectangleSegments(P1,P2,num_points=60,edge_length=-1):
  P11=[P2[0],P1[1]]
  P22=[P1[0],P2[1]]  
  npoints=np.floor(num_points/4)
  p_1,v_1=LineSegments(P1,P11,npoints,edge_length)
  p_2,v_2=LineSegments(P11,P2,npoints,edge_length)  
  p_3,v_3=LineSegments(P2,P22,npoints,edge_length)
  p_4,v_4=LineSegments(P22,P1,npoints,edge_length)
  p,v=AddSegments(p_1,p_2)
  p,v=AddSegments(p,p_3)
  p,v=AddSegments(p,p_4)
  return p,v

# Rectangle
def MyRectangle(P1,P2,edge_length=-1):
  P11=[P2[0],P1[1]]
  P22=[P1[0],P2[1]]
  
    
  p_1,v_1=MyLineSegments(P1,P11,edge_length)
  p_2,v_2=MyLineSegments(P11,P2,edge_length)  
  p_3,v_3=MyLineSegments(P2,P22,edge_length)
  p_4,v_4=MyLineSegments(P22,P1,edge_length)
  p,v=AddSegments(p_1,p_2)
  p,v=AddSegments(p,p_3)
  p,v=AddSegments(p,p_4)
  return p,v

# Rectangle
def MyParallelepiped(basen,xl,yl,zl):
  xshift=np.array([xl,0,0])
  yshift=np.array([0,yl,0])
  zshift=np.array([0,0,zl])

  coord=[];

  # down face 
  coord.append(list(basen))
  coord.append(list(basen+xshift))
  coord.append(list(basen+xshift+yshift))
  coord.append(list(basen+yshift))

  # up face   
  coord.append(list(basen + zshift))
  coord.append(list(basen+xshift + zshift))
  coord.append(list(basen+xshift+yshift + zshift))
  coord.append(list(basen+yshift + zshift))

  
  topol=[
        [0, 1, 2, 3],
        [4, 5, 6, 7],
        [0, 4, 5 ,1],
        [1, 2, 6, 5],
        [2, 6, 7, 3],
        [0, 4, 7, 3],
        ]

  
  # topol=[
  #   [0, 1],
  #   [1, 2],
  #   [2, 3],
  #   [3, 0],
  #   [0, 4],
  #   [1, 5],
  #   [2, 6],
  #   [3, 7],
  #   [4, 5],
  #   [5, 6],
  #   [6, 7],
  #   [7, 0],
  # ]
  
  return coord,topol

# Append Curves (algorithm for surfaces of four vertices)
def Add3dSurface(coord1,topol1,coord2,topol2):
  # make one list   
  coord=coord1+coord2
  topol=topol1+(np.array(topol2)+len(coord1)).tolist()
  nc1=len(coord1)
  nc2=len(coord2)
  nt1=len(topol1)
  nt2=len(topol2)
  
  todelete = range(nc2)
  replace  = range(nc1+nc2)
  m=0
  for i in range(nc2):
    j=0
    delete = False 
    while (j<nc1 and (not delete) ):
      delete=SamePoint(coord[i+nc1],coord[j],1e-9)
      if delete == True:
        todelete[m]= i+nc1
        replace[nc1+i] = j
        m=m+1
      j=j+1

  j=1
  for i in range(nc2):
    if (replace[nc1+i] >= nc1):
      replace[nc1+i]=nc1-1+j
      j=j+1

  for j in range(m-1,-1,-1):
    del coord[ todelete[j] ]

  for i in range(nt2):
    n_old = topol[i+nt1][0] 
    topol[i+nt1][0] = replace[n_old]
    n_old = topol[i+nt1][1] 
    topol[i+nt1][1] = replace[n_old]
    n_old = topol[i+nt1][2] 
    topol[i+nt1][2] = replace[n_old]
    n_old = topol[i+nt1][3] 
    topol[i+nt1][3] = replace[n_old]

  aux = [row[:] for row in topol]
  for i in range(len(aux)):
    aux[i].sort()

  for i in range(nt1+nt2-1,nt1-1,-1):
    for j in range(nt1-1,-1,-1):
      if (cmpT(aux[i],aux[j])):
        del (topol[i])
        del (aux[i])
        break

    
  #topol = [ tuple(topol[j]) for j in range (len(topol)) ]
  
  return coord,topol;


# return vertex coordinates fixed to the sphere of radius radius
def vertex(x, y, z, radius):
    
  dist = np.sqrt(x**2 + y**2 + z**2)

  return [(i * radius) / dist for i in (x,y,z)]
  
# make the base icosahedron given center coordinates and radius
def MyIcosahedron(center,radius):

  xshift=np.array([center[0],0,0])
  yshift=np.array([0,center[1],0])
  zshift=np.array([0,0,center[2]])
  
  # Golden ratio
  PHI = (1 + np.sqrt(5)) / 2
  
  verts = [
    vertex(-1,  PHI, 0, radius),
    vertex( 1,  PHI, 0, radius),
    vertex(-1, -PHI, 0, radius),
    vertex( 1, -PHI, 0, radius),
    
    vertex(0, -1,  PHI, radius),
    vertex(0,  1,  PHI, radius),
    vertex(0, -1, -PHI, radius),
    vertex(0,  1, -PHI, radius),
    
    vertex( PHI, 0, -1, radius),
    vertex( PHI, 0,  1, radius),
    vertex(-PHI, 0, -1, radius),
    vertex(-PHI, 0,  1, radius),
  ]

  for i in range (0,len(verts)):
      verts[i] = list(verts[i]+xshift)
      verts[i] = list(verts[i]+yshift)
      verts[i] = list(verts[i]+zshift)


  faces = [
      # 5 faces around point 0
      [0, 11, 5],
      [0, 5, 1],
      [0, 1, 7],
      [0, 7, 10],
      [0, 10, 11],
      # Adjacent faces
      [1, 5, 9],
      [5, 11, 4],
      [11, 10, 2],
      [10, 7, 6],
      [7, 1, 8],
      # 5 faces around 3
      [3, 9, 4],
      [3, 4, 2],
      [3, 2, 6],
      [3, 6, 8],
      [3, 8, 9],  
      # Adjacent faces
      [4, 9, 5],
      [2, 4, 11],
      [6, 2, 10],
      [8, 6, 7],
      [9, 8, 1],
  ]

  return verts,faces


# Find a middle point and project to the unit sphere
def middle_point(point_1, point_2, verts, radius, middle_point_cache):
  # We check if we have already cut this edge first
  # to avoid duplicated verts
  smaller_index = min(point_1, point_2)
  greater_index = max(point_1, point_2)

  key = '{0}-{1}'.format(smaller_index, greater_index)
  if key in middle_point_cache:
    return middle_point_cache[key]
  
  # If it's not in cache, then we can cut it
  vert_1 = verts[point_1]
  vert_2 = verts[point_2]
  middle = [sum(i)/2 for i in zip(vert_1, vert_2)]
  verts.append(vertex(middle[0], middle[1], middle[2], radius))
  
  index = len(verts) - 1
  middle_point_cache[key] = index

  return index

# Generate mesh  
def Do3DGrid(coord,topol,edge_length=-1,tri_refine=None):
  info = MeshInfo()
  info.set_points(coord)
  info.set_facets(topol)

  if tri_refine!=None:
    mesh = build(info,refinement_func=tri_refine)#,min_angle=27.0)
  elif edge_length<=0:
    mesh = build(info)   
  else:
    mesh = build(info,max_volume=edge_length**3/6.0)#,min_angle=30.0)
  
  
  return mesh


# List of points
def PointSegments(p):
  p1=np.array(p)
  delta=np.min(np.sqrt(np.sum((p1[1:]-p1[:-1])**2,axis=1)))
  Pall=[(x[0],x[1]) for x in p]  
  closed=False  
  if SamePoint(p1[0],p1[-1],delta)==True:
    Pall=Pall[:-1]  
    closed=True    
    
  vertices=[(j,j+1) for j in range(0,len(Pall)-1,1)]
  if closed==True:
    vertices+=[(len(Pall)-1,0)]  
  
  return Pall,vertices

#Connect two different polygons  
def AddSegments(P1,P2,closed=False):  
  p1=np.array(P1)
  p2=np.array(P2)
  # find smallest distance within points p1 and p2
  min1=np.min(np.sqrt(np.sum((p1[1:]-p1[:-1])**2,axis=1)))
  min2=np.min(np.sqrt(np.sum((p2[1:]-p2[:-1])**2,axis=1)))
  delta=np.min([min1,min2])
  
  # Add second curve to first curve 
  del_first=SamePoint(p1[-1],p2[0],delta)
  Pall=P1[:]  
  if del_first==True:
    Pall+=P2[1:]
  else:
    Pall+=P2
  
  # check if Pall is closed 
  del_last=SamePoint(Pall[-1],p1[0],delta)
  if del_last==True:
    Pall=Pall[:-1]
    
  vertices=[(j,j+1) for j in range(0,len(Pall)-1,1)]
  if (del_last==True) or (closed==True):
    vertices+=[(len(Pall)-1,0)]
  
  return Pall,vertices;  


# Append Curves
def AddCurves(p1,v1,p2,v2):
  # make one list   
  p=p1+p2
  v2n=[(v2[j][0]+len(p1),v2[j][1]+len(p1)) for j in range(len(v2))]
  v=v1+v2n
  v_out=list ( list(v[j]) for j in range(len(v)) ) 

  np1=len(p1)
  np2=len(p2)
  nv1=len(v1)
  nv2=len(v2)

  todelete = range(np2)
  replace  = range(np1+np2)

  m=0
  for i in range(np2):
    j=0
    delete = False 
    while (j<np1 and (not delete) ):
      delete=SamePoint(p[i+np1],p[j],1e-9)
      if delete == True:
        todelete[m]= i+np1
        replace[np1+i] = j
        m=m+1
      j=j+1


  j=1
  for i in range(np2):
    if (replace[np1+i] >= np1):
      replace[np1+i]=np1-1+j
      j=j+1


  for j in range(m-1,-1,-1):
    del p[ todelete[j] ]


  for i in range(nv2):
    n_old = v_out[i+nv1][0] 
    v_out[i+nv1][0] = replace[n_old]
    n_old = v_out[i+nv1][1] 
    v_out[i+nv1][1] = replace[n_old]

 
  for i in range(nv2-1,-1,-1):
    for j in range(nv1-1,-1,-1):
      if ( v_out[i+nv1][:] == v_out[j][:]) or (v_out[i+nv1][:] == v_out[j][::-1] ):
        del (v_out[ i+nv1 ] )
        break
      

  
  v_out = [ tuple( v_out[j] ) for j in range( len(v_out) ) ]
  
  return p,v_out;

def distance_segment_point(segment_begin,segment_end,point):
  a=np.array(segment_begin)
  b=np.array(segment_end)
  p=np.array(point)
  
  v=b-a
  nrm2 = np.sum(v)**2
  u=np.dot(p-a,v)/np.dot(v,v) 
  if u > 1 :
    u = 1
  if u < 0 :
    u = 0
  q = a + u * v   
  dq = q-p
  dist = np.sqrt(np.dot( dq,  dq ) )

  return dist;

def refine_on_strip(points, vertices, tri_points, strip_width):
  num_p=len(points)
  center_tri = np.sum(np.array(tri_points), axis=0)/3.
  x=center_tri[0]
  y=center_tri[1]
  r =np.sqrt(x**2+y**2)  

  refine_flag=False
  for j in range(num_p):
    p=np.array(point[j])
    if ( np.sqrt(np.dot(center_tri-p)) <strip_width):
      refine_flag=True

  return flag_refine; 





  
# Generate mesh  
def DoTriMesh(points,vertices,edge_length=-1,holes=[],tri_refine=None):
  info = triangle.MeshInfo()
  info.set_points(points)
  if len(holes)>0:
    info.set_holes(holes)
  info.set_facets(vertices)


  if tri_refine!=None:
    mesh = triangle.build(info,refinement_func=tri_refine,min_angle=27.0)
  elif edge_length<=0:
    mesh = triangle.build(info)   
  else:
    mesh = triangle.build(info,max_volume=0.5*edge_length**2,min_angle=30.0)
  
  mesh_points = np.array(mesh.points)
  mesh_elements = np.array(mesh.elements)
  
  
#  plt.triplot(mesh_points[:, 0], mesh_points[:, 1], mesh_elements,) 
#  plt.show()
#  return mesh_points,mesh_elements;  
  return mesh

def write_poly(points,vertices,filename,out_format):
  #mesh_points = np.array(mesh.points)
  coord=[[float(p[0]),float(p[1]),0.0] for j, p in enumerate(points)]
  lines=[[v[0],v[1]] for j, v in enumerate(vertices)]
  if ( out_format=='dat'):
    out_filename=filename+".dat"
    
    f_poly = open(out_filename, 'w')
    f_poly.write(str(len(coord))+'\n')
    f_poly.write(str(len(lines))+"\n")
    for i, p in enumerate(coord):
      f_poly.write(str(p[0]) + " " + str(p[1])+" " + str(p[2])+"\n")
    for i, t in enumerate(lines):
      f_poly.write(str(t[0]+1) + " " + str(t[1]+1) +"\n")
            
  if ( out_format=='vtk'):
    # file name
    out_filename=filename+".vtk"
      
    # save all Vtk data for unstructured grid into a vtk object
    vtk_object = VtkData( UnstructuredGrid(coord,line=lines))
        
    # write vtk object to file (both ascii and binary versions)
    vtk_object.tofile(out_filename,'ascii')



def write_mesh(mesh,filename,out_format):
    #mesh_points = np.array(mesh.points)
    mesh_points=[[p[0],p[1],0] for j, p in enumerate(mesh.points)]
    triang=[[p[0],p[1],p[2]] for j, p in enumerate(mesh.elements)]

    
    if ( out_format=='dat'):
        # file name
        out_filename=filename+".dat"
        
        f_grid = open(out_filename, 'w')
        f_grid.write(str(len(mesh.points))+'! nodes number \n')
        f_grid.write(str(len(mesh.elements))+' 4 3d  ! cell number, node in cell, mesh type \n')
        for i, p in enumerate(mesh.points):
            f_grid.write(str(p[0]) + " " + str(p[1])+"\n")
        for itria in range(len(triang)):
          f_grid.write(str(triang[itria][0]+1) + " " + str(triang[itria][1]+1) + " " + str(triang[itria][2]+1)+ " "+str(itria+1)+"\n")
            
    if ( out_format=='vtk'):
        # file name
        out_filename=filename+".vtk"
        
        # save all Vtk data for unstructured grid into a vtk object
        vtk_object = VtkData( UnstructuredGrid(mesh_points,triangle=triang))
        
        # write vtk object to file (both ascii and binary versions)
        vtk_object.tofile(out_filename,'ascii')

def read_grid(filename):
  f_grid = open(filename, 'r')
  # writing data
  input_lines = f_grid.readlines()
  nnode= int(input_lines[0])
  ncell= int(input_lines[1])
  coord=np.zeros([nnode,3])
  topol=np.zeros([ncell,4],dtype=int)
  flags=np.zeros(ncell,dtype=int)
  inode=0
  for line in input_lines[2:2+nnode]:
    coord[inode][:]=[float(w) for w in line.split()[0:3]]
    inode=inode+1
  icell=0
  for line in input_lines[2+nnode:]:
    topol[icell][:]=[int(w)-1 for w in line.split()[0:4]]
    flags[icell]=line.split()[-1]
    icell=icell+1


  return coord,topol,flags,


def write_grid(coord,topol,filename,out_format,flags=None):
    if ( out_format=='dat'):
        f_grid = open(filename, 'w')
        # writing data
        f_grid.write(str(len(coord))+'\n')
        f_grid.write(str(len(topol))+"\n")
        for inode  in range(len(coord)):
            f_grid.write(
              str(coord[inode][0]) + " " +
              str(coord[inode][1]) + " " +
              str(coord[inode][2]) +"\n")
        if (flags is None):
          for itria in range(len(topol)):
            f_grid.write(
              str(topol[itria][0]+1) + " " + 
              str(topol[itria][1]+1) + " " + 
              str(topol[itria][2]+1) + " " + 
              str(topol[itria][3]+1) + " " + 
              str(itria+1)+"\n")
        else:
          for itria in range(len(topol)):
            f_grid.write(
              str(topol[itria][0]+1) + " " + 
              str(topol[itria][1]+1) + " " + 
              str(topol[itria][2]+1) + " " + 
              str(topol[itria][3]+1) + " " + 
              str(flags[itria])+"\n")
        f_grid.close()
          
            
    if ( out_format=='vtk'):
      coord_list=coord.tolist()
      topol_list=topol.tolist()      

      # file name
      
      celldata = CellData(Scalars(flags,name='materials'))
      
      # save all Vtk data for unstructured grid into a vtk object
      vtk_object = VtkData( 
        UnstructuredGrid(coord_list,tetra=topol_list),
        celldata)
      
      # write vtk object to file (both ascii and binary versions)
      vtk_object.tofile(filename,'ascii')


def make_bar(coord,topol):
  bar_tetra=np.empty([len(topol),3])
  for icell in range(len(topol)):
    nodes = topol[icell]
    bar_tetra[icell][:] = np.sum(np.array(coord[nodes]), axis=0)/4.0
  return bar_tetra 

def make_size(coord,topol):
    volume=np.empty(len(topol))
    for icell in range(len(topol)):
      nodes=topol[icell]
      n0=nodes[0]
      n1=nodes[1]
      n2=nodes[2]
      n3=nodes[3]
      s1n0=coord[n0][0]
      s2n0=coord[n0][1]
      s3n0=coord[n0][2]
      s1n1=coord[n1][0]
      s2n1=coord[n1][1]
      s3n1=coord[n1][2]
      s1n2=coord[n2][0]
      s2n2=coord[n2][1]
      s3n2=coord[n2][2]
      s1n3=coord[n3][0]
      s2n3=coord[n3][1]
      s3n3=coord[n3][2]
      
      volume[icell]=abs(
        (s1n3-s1n0)*((s2n1-s2n0)*(s3n2-s3n0)-(s3n1-s3n0)*(s2n2-s2n0)) +
        (s2n3-s2n0)*((s3n1-s3n0)*(s1n2-s1n0)-(s1n1-s1n0)*(s3n2-s3n0)) +
        (s3n3-s3n0)*((s1n1-s1n0)*(s2n2-s2n0)-(s2n1-s2n0)*(s1n2-s1n0)))/6.0

        
    return volume;


