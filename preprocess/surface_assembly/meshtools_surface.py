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
    
  node_closest=[nodes[index]]
   
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
  d=np.sqrt(dp[0]**2+dp[1]**2)
  ret=False  
  if d<delta:
    ret=True
  return ret;
 



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
  coord=[[float(p[0]),float(p[1]),float(p[2])] for j, p in enumerate(points)]
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
        f_grid.write(str(len(mesh.points))+'\n')
        f_grid.write(str(len(mesh.elements))+' 3 surface  ! cell number, node in cell, mesh type \n')
        for i, p in enumerate(mesh.points):
            f_grid.write(
              str(p[0]) + " " + 
              str(p[1]) + " " + 
              str(p[2]) + "\n")
        for itria in range(len(triang)):
          f_grid.write(str(triang[itria][0]+1) + " " + 
                       str(triang[itria][1]+1) + " " + 
                       str(triang[itria][2]+1) + " " + 
                       str(itria+1)+"\n")
            
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
  nnode= int(input_lines[0].split()[0])
  ntria= int(input_lines[1].split()[0])
  try:
    nnodeincell= int(input_lines[1].split()[1])
  except:
    print('No number of nodes defined')
  coord=np.zeros([nnode,3])
  triang=np.zeros([ntria,3],dtype=int)
  flags=np.zeros(ntria,dtype=int)
  inode=0
  for line in input_lines[2:2+nnode]:
    coord[inode][:]=[float(w) for w in line.split()[0:3]]
    inode=inode+1
  itria=0
  for line in input_lines[2+nnode:]:
    triang[itria][:]=[int(w)-1 for w in line.split()[0:3]]
    flags[itria]=line.split()[3]
    itria=itria+1
  


  return coord,triang,flags,


def write_grid(coord,triang,filename,out_format,flags=None):
    if ( out_format=='dat'):
        f_grid = open(filename, 'w')
        # writing data
        f_grid.write(str(len(coord))+'\n')
        f_grid.write(str(len(triang))+"\n")
        for inode  in range(len(coord)):
            f_grid.write(
              str(coord[inode][0]) + " " + 
              str(coord[inode][1]) + " " + 
              str(coord[inode][2]) +"\n")
        if (flags is None):
          for itria in range(len(triang)):
            f_grid.write(
              str(triang[itria][0]+1) + " " + 
              str(triang[itria][1]+1) + " " + 
              str(triang[itria][2]+1) + " " +
              str(itria+1)+"\n")
        else:
          for itria in range(len(triang)):
            f_grid.write(
              str(triang[itria][0]+1) + " " +
              str(triang[itria][1]+1) + " " +
              str(triang[itria][2]+1) + " " +
              str(flags[itria])+"\n")
        f_grid.close()
          
            
    if ( out_format=='vtk'):
      coord_list=[]
      triang_list=[]
      for i in range(len(coord)):
        coord_list.append([coord[i][0],coord[i][1],coord[i][2]])
      for i in range(len(triang)):
        triang_list.append([triang[i][0],triang[i][1],triang[i][2]])
      

      # file name
      
      celldata = CellData(Scalars(flags,name='materials'))
      
      # save all Vtk data for unstructured grid into a vtk object
      vtk_object = VtkData( UnstructuredGrid(coord_list,triangle=triang_list))
      
      # write vtk object to file (both ascii and binary versions)
      vtk_object.tofile(filename,'ascii')


def make_bar(coord,triang):
    bar=[]
    for itria in range(len(triang)):
        nodes = triang[itria]
        center_tri = np.sum(np.array(coord[nodes]), axis=0)/3.0
        x = center_tri[0]
        y = center_tri[1]
        z = center_tri[2]
        bar.append([x,y,z])
    nbar=np.array(bar)
    return nbar

def make_size(coord,triang):
    area=np.empty(len(triang))
    for itria in range(len(triang)):
        nodes = triang[itria]
        vec12=(coord[:][nodes[0]]-coord[:][nodes[1]])
        vec13=(coord[:][nodes[0]]-coord[:][nodes[2]])
        cross_prod = cross(vec12,vec13)
        area[itria]=0.5*np.sqrt(np.dot(cross_prod,cross_prod))
    return area

def make_node(coord,triang,area_tria):
    area_node=np.zeros(len(triang))
    for itria in range(len(triang)):
      for node in triang[itria][:]:
        area_node[node]=area_node[node]+area_tria[itria]/3.0
    return area_node;


def cross(vecA,vecB):
  cross=np.zeros(3)
  cross[0] = vecA[1]*vecB[2]-vecA[2]*vecB[1]
  cross[1] = vecA[2]*vecB[0]-vecA[0]*vecB[2]
  cross[2] = vecA[0]*vecB[1]-vecA[1]*vecB[0]
  return cross;
  
