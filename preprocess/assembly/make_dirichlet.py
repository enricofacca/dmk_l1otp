#!/usr/bin/env python
import numpy as np
import common
import meshtools as mt
import sys
import timecell

def define_dirichlet_node(coord,topol,flag_dirichlet,extra_dirichlet,flag_equation):
    if (flag_dirichlet =='zerobc'):
        dirichlet_nodes=mt.boundary_nodes(coord,topol)
        dirichlet_values=np.zeros(len(dirichlet_nodes))

    elif (flag_dirichlet =='leftPright0'):
        print(min(coord[:,0]), max(coord[:,0]))
        XL=min(coord[:,0])
        XR=max(coord[:,0])

        
        P=float(extra_dirichlet)
        dirichlet_nodes=[]
        dirichlet_values=[]
        is_dirichlet_dof=[]
        
        for inode,node in enumerate(coord):
            if ( node[0] == XL ):
                dirichlet_nodes.append(inode)
                dirichlet_values.append(P)
            if ( node[0] == XR ):
                dirichlet_nodes.append(inode)
                dirichlet_values.append(0.0)

        dirichlet_nodes=np.asarray( dirichlet_nodes,dtype='int')
        dirichlet_values=np.asarray(dirichlet_values)
        is_dirichlet_dof=np.asarray(is_dirichlet_dof)

    elif (flag_dirichlet =='example1'):
        dirichlet_nodes=[]
        dirichlet_values=[]
        is_dirichlet_dof=[]
        for inode,node in enumerate(coord):
            x=node[0]
            y=node[1]
            if ( abs(y)< 1e-12 ):
                if ( ( x>0-1e12 ) or (x < 0.3+1e-12) ) :
                     dirichlet_nodes.append([inode,inode])
                     dirichlet_values.append([0.0,0.0])
                     is_dirichlet_dof.append([1.0,1.0])
        dirichlet_nodes=np.asarray( dirichlet_nodes,dtype='int')
        dirichlet_values=np.asarray(dirichlet_values)
        is_dirichlet_dof=np.asarray(is_dirichlet_dof)
                     
    elif (flag_dirichlet =='example1_inside'):
        dirichlet_nodes=[]
        dirichlet_values=[]
        is_dirichlet_dof=[]
        for inode,node in enumerate(coord):
            x=node[0]
            y=node[1]
            if ( abs(y-0.5)< 1e-12 ):
                if ( ( x>0.35-1e-12 ) and (x < 0.65+1e-12) ) :
                     print(inode)
                     dirichlet_nodes.append([inode,inode])
                     dirichlet_values.append([0.0,0.0])
                     is_dirichlet_dof.append([1.0,1.0])
        dirichlet_nodes=np.asarray( dirichlet_nodes,dtype='int')
        dirichlet_values=np.asarray(dirichlet_values)
        is_dirichlet_dof=np.asarray(is_dirichlet_dof)
                     
    elif (flag_dirichlet =='example1_semi_inside'):
        dirichlet_nodes=[]
        dirichlet_values=[]
        is_dirichlet_dof=[]
        for inode,node in enumerate(coord):
            x=node[0]
            y=node[1]
            if ( abs(y-0.0)< 1e-12 ):
                if ( ( x>0.35-1e-12 ) and (x < 0.65+1e-12) ) :
                     print(inode)
                     dirichlet_nodes.append([inode,inode])
                     dirichlet_values.append([0.0,0.0])
                     is_dirichlet_dof.append([1.0,1.0])
        dirichlet_nodes=np.asarray( dirichlet_nodes,dtype='int')
        dirichlet_values=np.asarray(dirichlet_values)
        is_dirichlet_dof=np.asarray(is_dirichlet_dof)
                     
    elif (flag_dirichlet =='example2'):
        dirichlet_nodes=[]
        dirichlet_values=[]
        is_dirichlet_dof=[]
        for inode,node in enumerate(coord):
            x=node[0]
            y=node[1]
            if ( abs(x)< 1e-12 ):
                if ( ( y>-1e12 ) or (y < 0.3+1e-12) ) :
                     dirichlet_nodes.append([inode,inode])
                     dirichlet_values.append([0.0,0.0])
                     is_dirichlet_dof.append([1.0,1.0])
        dirichlet_nodes=np.asarray( dirichlet_nodes,dtype='int')
        dirichlet_values=np.asarray(dirichlet_values)
        is_dirichlet_dof=np.asarray(is_dirichlet_dof)

    elif (flag_dirichlet =='example2_inside'):
        dirichlet_nodes=[]
        dirichlet_values=[]
        is_dirichlet_dof=[]
        for inode,node in enumerate(coord):
            x=node[0]
            y=node[1]
            if ( abs(x-0.5)< 1e-12 ):
                if ( ( y>0.35-1e12 ) or (y < 0.65+1e-12) ) :
                     dirichlet_nodes.append([inode,inode])
                     dirichlet_values.append([0.0,0.0])
                     is_dirichlet_dof.append([1.0,1.0])
        dirichlet_nodes=np.asarray( dirichlet_nodes,dtype='int')
        dirichlet_values=np.asarray(dirichlet_values)
        is_dirichlet_dof=np.asarray(is_dirichlet_dof)

    elif (flag_dirichlet =='example2_semi_inside'):
        dirichlet_nodes=[]
        dirichlet_values=[]
        is_dirichlet_dof=[]
        for inode,node in enumerate(coord):
            x=node[0]
            y=node[1]
            if ( abs(x-0.0)< 1e-12 ):
                if ( ( y>0.6-1e-12 ) and (y < 0.9+1e-12) ) :
                     dirichlet_nodes.append([inode,inode])
                     dirichlet_values.append([0.0,0.0])
                     is_dirichlet_dof.append([1.0,1.0])
        dirichlet_nodes=np.asarray( dirichlet_nodes,dtype='int')
        dirichlet_values=np.asarray(dirichlet_values)
        is_dirichlet_dof=np.asarray(is_dirichlet_dof)

    elif (flag_dirichlet =='example3'):
        dirichlet_nodes=[]
        dirichlet_values=[]
        is_dirichlet_dof=[]
        for inode,node in enumerate(coord):
            x=node[0]
            y=node[1]
            if ( abs(y-0.0) < 1e-12 ):
                if ( abs(x-0.0) < 1e-12 ) : 
                     dirichlet_nodes.append([inode,inode])
                     dirichlet_values.append([0.0,0.0])
                     is_dirichlet_dof.append([1.0,1.0])
                if ( abs(x-2.4) < 1e-12 ) : 
                     dirichlet_nodes.append([inode,inode])
                     dirichlet_values.append([0.0,0.0])
                     is_dirichlet_dof.append([0.0,1.0])
        dirichlet_nodes=np.asarray( dirichlet_nodes,dtype='int')
        dirichlet_values=np.asarray(dirichlet_values)
        is_dirichlet_dof=np.asarray(is_dirichlet_dof)
        
    elif (flag_dirichlet =='example3_inside'):
        dirichlet_nodes=[]
        dirichlet_values=[]
        is_dirichlet_dof=[]
        for inode,node in enumerate(coord):
            x=node[0]
            y=node[1]
            if ( abs(y-0.8) < 1e-12 ):
                if ( abs(x-0.3) < 1e-12 ) : 
                     dirichlet_nodes.append([inode,inode])
                     dirichlet_values.append([0.0,0.0])
                     is_dirichlet_dof.append([1.0,1.0])
                if ( abs(x-2.7) < 1e-12 ) : 
                     dirichlet_nodes.append([inode,inode])
                     dirichlet_values.append([0.0,0.0])
                     is_dirichlet_dof.append([0.0,1.0])
        dirichlet_nodes=np.asarray( dirichlet_nodes,dtype='int')
        dirichlet_values=np.asarray(dirichlet_values)
        is_dirichlet_dof=np.asarray(is_dirichlet_dof)
        
    elif (flag_dirichlet =='example3_semi_inside'):
        dirichlet_nodes=[]
        dirichlet_values=[]
        is_dirichlet_dof=[]
        for inode,node in enumerate(coord):
            x=node[0]
            y=node[1]
            if ( abs(y-0.0) < 1e-12 ):
                if ( abs(x-0.0) < 1e-12 ) : 
                     dirichlet_nodes.append([inode,inode])
                     dirichlet_values.append([0.0,0.0])
                     is_dirichlet_dof.append([1.0,1.0])
                if ( abs(x-2.4) < 1e-12 ) : 
                     dirichlet_nodes.append([inode,inode])
                     dirichlet_values.append([0.0,0.0])
                     is_dirichlet_dof.append([0.0,1.0])
        dirichlet_nodes=np.asarray( dirichlet_nodes,dtype='int')
        dirichlet_values=np.asarray(dirichlet_values)
        is_dirichlet_dof=np.asarray(is_dirichlet_dof)
        
    elif (flag_dirichlet =='example_book'):
        dirichlet_nodes=[]
        dirichlet_values=[]
        is_dirichlet_dof=[]
        for inode,node in enumerate(coord):
            x=node[0]
            y=node[1]
            if ( abs(y-0.5) < 1e-12 ):
                if ( abs(x-0.0) < 1e-12 ) : 
                     dirichlet_nodes.append([inode,inode])
                     dirichlet_values.append([0.0,0.0])
                     is_dirichlet_dof.append([1.0,1.0])
                if ( abs(x-1.0) < 1e-12 ) : 
                     dirichlet_nodes.append([inode,inode])
                     dirichlet_values.append([0.0,0.0])
                     is_dirichlet_dof.append([1.0,1.0])
        dirichlet_nodes=np.asarray( dirichlet_nodes,dtype='int')
        dirichlet_values=np.asarray(dirichlet_values)
        is_dirichlet_dof=np.asarray(is_dirichlet_dof)
        
    elif (flag_dirichlet =='example_internet1'):
        dirichlet_nodes=[]
        dirichlet_values=[]
        is_dirichlet_dof=[]
        for inode,node in enumerate(coord):
            x=node[0]
            y=node[1]
            if ( abs(x-0.0)< 1e-12 ):
                if ( ( y>0.0-1e-12 ) and (y < 0.75+1e-12) ) :
                     dirichlet_nodes.append([inode,inode])
                     dirichlet_values.append([0.0,0.0])
                     is_dirichlet_dof.append([1.0,1.0])
        dirichlet_nodes=np.asarray( dirichlet_nodes,dtype='int')
        dirichlet_values=np.asarray(dirichlet_values)
        is_dirichlet_dof=np.asarray(is_dirichlet_dof)
        
    elif (flag_dirichlet =='example_internet1_semi_inside'):
        dirichlet_nodes=[]
        dirichlet_values=[]
        is_dirichlet_dof=[]
        for inode,node in enumerate(coord):
            x=node[0]
            y=node[1]
            if ( abs(x-0.0)< 1e-12 ):
                if ( ( y>0.0-1e-12 ) and (y < 0.75+1e-12) ) :
                     dirichlet_nodes.append([inode,inode])
                     dirichlet_values.append([0.0,0.0])
                     is_dirichlet_dof.append([1.0,1.0])
        dirichlet_nodes=np.asarray( dirichlet_nodes,dtype='int')
        dirichlet_values=np.asarray(dirichlet_values)
        is_dirichlet_dof=np.asarray(is_dirichlet_dof)
        
    elif (flag_dirichlet =='example_internet2'):
        dirichlet_nodes=[]
        dirichlet_values=[]
        is_dirichlet_dof=[]
        for inode,node in enumerate(coord):
            x=node[0]
            y=node[1]
            if ( abs(x-0.0)< 1e-12 ):
                if ( ( y>0.0-1e-12 ) and (y < 0.75+1e-12) ) :
                     dirichlet_nodes.append([inode,inode])
                     dirichlet_values.append([0.0,0.0])
                     is_dirichlet_dof.append([1.0,1.0])
        dirichlet_nodes=np.asarray( dirichlet_nodes,dtype='int')
        dirichlet_values=np.asarray(dirichlet_values)
        is_dirichlet_dof=np.asarray(is_dirichlet_dof)
        
    elif (flag_dirichlet =='example_internet2_semi_inside'):
        dirichlet_nodes=[]
        dirichlet_values=[]
        is_dirichlet_dof=[]
        for inode,node in enumerate(coord):
            x=node[0]
            y=node[1]
            if ( abs(x-0.0)< 1e-12 ):
                if ( ( y>0.25-1e-12 ) and (y < 1.0+1e-12) ) :
                     dirichlet_nodes.append([inode,inode])
                     dirichlet_values.append([0.0,0.0])
                     is_dirichlet_dof.append([1.0,1.0])
        dirichlet_nodes=np.asarray( dirichlet_nodes,dtype='int')
        dirichlet_values=np.asarray(dirichlet_values)
        is_dirichlet_dof=np.asarray(is_dirichlet_dof)
        
    else:
        print( 'Flag ', flag_dirichlet, ' not supported')


    
    return is_dirichlet_dof, dirichlet_nodes, dirichlet_values;
        


def make_dirichlet(time, node_coord, dirichlet_function, dir_nodes, dir_values):
    f_cell=0.0
    steady=True
    dir_nodes[:]=0
    dir_values[:]=0
    
    list_nodes, list_values, steady =dirichlet(t,infopath)
    for i in range(len(list_nodes)):
        inode,dist=mt.FindClosestNode(
            range(len(coord)),coord,list_nodes[i])
        dir_nodes[inode]=float(inode)
        dir_values[inode]=list_values[i]
    return dir_nodes,dir_value,steady;

# An accurate discontinuous Galerkin method for solving point-source
# Eikonal equation in 2-D heterogeneous anisotropic media
# Bouteiller,1 M. Benjemaa,2 L. Metivier1,3 and J. Virieux1
# if (flag_dirichlet == 'bouteiller18_tc1'):
#     def dirichlet(t,infopath):
#         coord_points=[[2.0,2.0,0.0]]
#         values=[[0.0]]
#         steady=True
#         return coord_points, values,steady;
#         dirichlet_functions.append(dirichlet)

def example_dirichlet(flag_dirichlet,coord,extra,dir_values):
    # An accurate discontinuous Galerkin method for solving point-source
    # Eikonal equation in 2-D heterogeneous anisotropic media
    # Bouteiller,1 M. Benjemaa,2 L. Metivier1,3 and J. Virieux1
    print (flag_dirichlet)
    if (flag_dirichlet == 'bouteiller18_tc1'):
        box=0.21
        nnode=len(coord)
        # file_out=open(filename, 'w')
        # file_out.write(
        #     str(2) + " " + 
        #     str((nnode)) + " !  dim ndata" 
        #     +"\n")
        # file_out.write("time    0.0  \n")
        # file_out.write(str(nnode)+" \n")
        
        for i in range(len(coord)):
            x=coord[i][0]
            y=coord[i][1]
            if ( ( abs(x-2.0) <= box ) and ( abs(y-2.0) <= box ) ) :
                value=np.sqrt((x-2.0)**2+(y-2.0)**2)
                dir_values[i][:] = [1.0, value] 
  

    if (flag_dirichlet == 'hole_square'):
        box=0.1
        nnode=len(coord)
        j=0        
        for i in range(len(coord)):
            x=coord[i][0]
            y=coord[i][1]
           
            if ( abs(np.sqrt(x**2+y**2) - box) < 1e-8 ):
                j=j+1
                value=np.sqrt((x)**2+(y)**2)
                print( i, x,y, value )
                dir_values[i][:] = [1.0, value] 
        print ('dir nodes',j)
    return dir_values;


    

if __name__ == "__main__":
    if len(sys.argv) > 1:
        fin_ctrl=sys.argv[1]  
        fin_grid=sys.argv[2]
        fout_dirichlet=sys.argv[3]
        
        ############################################
        # reads flags and other values
        ############################################
        ctrl = open(fin_ctrl, "r")
        ctrl_lines = ctrl.readlines()
        ctrl.close()
        
        # flags grid
        flag_grid=str(common.read_column('flag_grid',0,ctrl_lines))
        extra_grid=str(common.read_column('flag_grid',1,ctrl_lines))

        # flags grid
        flag_dirichlet=str(common.read_column('flag_dirichlet',0,ctrl_lines))
        extra_dirichlet=str(common.read_column('flag_dirichlet',1,ctrl_lines))
        
        # reads coord topol and flags 
        coord, topol, flags = mt.read_grid(fin_grid)
        if ( coord.shape[1] == 2 ):
            zcoord = np.zeros([coord.shape[0],1])
            coord=np.append(coord, zcoord, axis=1)            
        dir_out=np.zeros([len(coord),2]);
        nodedir,divalue=define_dirichlet_node(coord,topol,flag_dirichlet,extra_dirichlet)
        dir_out[nodedir[:],0]=1.0
        dir_out[nodedir[:],1]=divalue[:]
                              
        file_out = open(fout_dirichlet, "w")
        timecell.write2file(file_out,0.0,True,dir_out,steady=False)
        timecell.write2file(file_out,0.0,False,dir_out,steady=True)
                
        file_out.close()

    else:
        raise SystemExit("usage:  python make_dirichlet.py 'input ctrl file' 'grid file' 'output file'"  )
