#!python
#
# program to draw the cut locus of point on the external equator of 
# a  torus with radius R and r with R>r.
#

import sys
sys.path.append('../')
import meshtools_surface as mt
import numpy as np


R=2.0
r=1.0

cut_point=(-(R+r), 0, 0)
conjugate_angle=np.pi/np.sqrt(r**2*(R+r))
arc_angle=np.pi-conjugate_angle


hpar=0.01


#
# opposite meridian ( all points)
#

point_opposite_meridian=[]
edge_opposite_meridian=[]

npoint=int(2*np.pi*r/float(hpar))+1
t=np.linspace(0,2*np.pi,npoint)

for i in range(npoint):
    x=R+r*np.cos(t[i])
    y=0.0
    z=r*np.sin(t[i])
    point_opposite_meridian.append([x,y,z])
for i in range(npoint-1):
    edge_opposite_meridian.append([i,i+1])
edge_opposite_meridian.append([npoint-1,0])

points=point_opposite_meridian
vertices=edge_opposite_meridian

#
# internal equator (exclude point 1 0 0 ) 
#
point_internal_equator=[]
edge_internal_equator=[]

npoint=int(2*np.pi*(R-r)/float(hpar))+1
print npoint
if ( npoint % 2 != 1 ) :
    npoint=npoint+1
t=np.linspace(0,2*np.pi,npoint)

for i in range(npoint):
    x=(R-r)*np.cos(t[i])
    y=(R-r)*np.sin(t[i])
    z=0
    point_internal_equator.append([x,y,z])

for i in range(npoint-1):
    edge_internal_equator.append([i,i+1])



points,vertices=mt.AddCurves(points,vertices,point_internal_equator,edge_internal_equator)


#
# arc along external equator with negative y (exclude point 3 0 0 )
# 
point_arc_minus=[]
edge_arc_minus=[]

npoint=int(arc_angle*(R+r)/float(hpar))+1
t=np.linspace(0,arc_angle,npoint)

for i in range(npoint):
    x=(R+r)*np.cos(-t[i])
    y=(R+r)*np.sin(-t[i])
    z=0
    point_arc_minus.append([x,y,z])

for i in range(npoint-1):
    edge_arc_minus.append([i,i+1])


points,vertices=mt.AddCurves(points,vertices,point_arc_minus,edge_arc_minus)


#
# arc along external equator with positive y (exclude point 3 0 0 )
# 
point_arc_plus=[]
edge_arc_plus=[]

npoint=int(arc_angle*(R+r)/float(hpar))+1
t=np.linspace(0,arc_angle,npoint)

for i in range(npoint):
    x=(R+r)*np.cos(t[i])
    y=(R+r)*np.sin(t[i])
    z=0
    point_arc_plus.append([x,y,z])


for i in range(npoint-1):
    edge_arc_plus.append([i,i+1])

points,vertices=mt.AddCurves(points,vertices,point_arc_plus,edge_arc_plus)

for i in range(len(points)):
    print points[i]
    points[i][0]=-points[i][0]
    


cut_name='cut_locus_torus_R'+str(R)+'_r'+str(r)


mt.write_poly(points,vertices,str(cut_name),'vtk')


eps=0.1
Rout=(R+r+eps)
Rin = (R+r-eps)
points=[]
points.append([Rout*np.cos(conjugate_angle),Rout*np.sin(conjugate_angle),0])
points.append([Rin*np.cos(conjugate_angle),Rin*np.sin(conjugate_angle),0])
ang2=2*np.pi-conjugate_angle
points.append([Rout*np.cos(ang2),Rout*np.sin(ang2),0])
points.append([Rin*np.cos(ang2),Rin*np.sin(ang2),0])
# useless point to make visit draw 3d 
points.append([0,0,1])

vertices=[[0,1],[2,3]]

cut_name='bound_cut_locus_torus_R'+str(R)+'_r'+str(r)


mt.write_poly(points,vertices,str(cut_name),'vtk')
        
