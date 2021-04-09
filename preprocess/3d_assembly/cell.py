#!/usr/bin/env python
import numpy as np

#
# prgrams to build functions
# taking real value on grid elements
#


#
# evaluate function at certain time 
#
def build(cell_values,functions,bar_cell,flags):
    for icell in range(len(cell_values)):
        for ifun in range(len(functions)):
            value = functions[ifun](
                bar_cell[icell][:],
                flags[icell])
            cell_values[icell]+=value
    return cell_values;


# build examples of functions from
#
# function_scale [0,+\infty[ \mapsto \REAL^{\Dim}
# 
# where \Dim is the Cell number of the grid
# 
# flag_scalar (string) : define the example
# if flag_scalar is float number 
# function(t)=float(flag_scalar)
# other case aare defined by the flag used.
#
# 
# The string "extra_info" it may contain the path to 
# a file contiang additional information 
# defing the function ( constants, range of time, etc. )
#


def example(flag,extra_info,length):
    # if is const contained in flag_kappa
    try :
        cnst=float(flag)
        def funct(coord,flag_mesh):
            value=cnst
            steady=True
            return cnst;
    #
    # user defined example
    # 
    except ValueError:
        #############################################
        # TDENS0 examples 
        # They can be used also for other data
        #############################################
        
        #  parabola centered at origin
        if (flag == 'par'):
            def funct(coord,flag_mesh):
                x=coord[0]
                y=coord[1]
                value=1.0+x**2+y**2
                return value;

        # parabola centered at (0.5, 0.5)
        if (flag == 'par0505'):
            def funct(coord,flag_mesh):
                x=coord[0]
                y=coord[1]
                value=0.1+(x-0.5)**2+(y-0.5)**2
                return value;
        
        # ossilating data  
        if (flag == 'sin'):
            def funct(coord,flag_mesh):
                x=coord[0]
                y=coord[1]
                value=2.0+sin(8*np.pi*x)+sin(8*np.pi*y)
                return value;
      
        # double parabola
        if (flag == 'par_eastwest'):
            def funct(coord,flag_mesh):
                x=coord[0]
                y=coord[1]
                # west
                if ( x <= 0.5 ):
                    value=0.01+2*(x-0.5)**2+y**2;
                # east
                if ( x > 0.5 ):
                    value=0.01+2*(x-0.5)**2+(y-0.5)**2;
                return value;
    
        # y_branch
        if (flag == 'ybranch'):
            # read source-sink point
            info_path=common.remove_comments(extra_info,'!')
            coord_points=[]
            values=[]
            input_file=open(info_path,'r')
            input_lines = input_file.readlines()
            for line in input_lines[0:3]:
                coord_dirac = [float(w) for w in line.split()[0:2]]
                branch_value = float(line.split()[2])
                values.append(branch_value)
                coord_points.append(coord_dirac)                    
            
            
                alpha=float(common.remove_comments(input_lines[4],'!'))
            
            # Build the optimal solution
            # add branching point (if exists)
        
            branch_points,branch_vertices = branch.y_branch(
                coord_points,values,alpha)
            input_file.close()
        
            def funct(coord,flag_mesh):
                x=coord[0]
                y=coord[1]
                value = 1.0e-10
                for segments in branch_vertices :
                    distance=mt.distance_segment_point(
                        branch_points[segments[0]],
                        branch_points[segments[1]],
                        point)
                    if ( distance < length / 1.5 ):
                        value = 1.0e10
                    
                return value;

#############################################
# OPTDENS 
# They can be used also for other data
#############################################
        if (flag == 'optdens_rectangle'):
            def funct(coord,flag_mesh):
                x=coord[0]
                y=coord[1]
                value=0.0
                if ((x >= 1.0/8.0) and (x<=3.0/8.0) and 
                    (y >= 1.0/4.0) and (y<=3.0/4.0) ): 
                    value=(x-1.0/8.0)*2.0
                if ((x > 3.0/8.0) and (x<5.0/8.0) and 
                    (y >= 1.0/4.0) and (y<=3.0/4.0) ):
                    value=0.5
                if ((x >= 5.0/8.0) and (x<=7.0/8.0) and 
                    (y >= 1.0/4.0) and (y<=3.0/4.0) ):
                    value=(7.0/8.0-x)*2.0
                    return value;
        else:
            print 'flag =' + flag + 'not found'
    
    return funct;

#  function to write data to file
def write2file(filename,data):
    ndim=data.shape[1]
    ndata=data.shape[0]
    file_out = open(filename, 'w')
    file_out.write(str(ndim)+" "+str(ndata)+" !dim number data \n")
    for icell in range(ndata):
        file_out.write(" ".join(map(str,data[icell][:]))  +"\n")
    file_out.close()

#  function to read data from file
def readfromfile(filepath,data):
    ndata=len(data)
    filein=open(str(filepath), 'r')
    input_lines = filein.readlines()
    dim=int(input_lines[0].split()[0])
    
    ndata=int(input_lines[0].split()[1])
    if len(input_lines) == (ndata+1):
        i=0
        for line in input_lines[1:]:
            data[i]=float(line)
            i=i+1
    else:
        print 'Dimension mismatch-Ndata=',len(input_lines),'len(array)=',ndata    
    filein.close()   
    return data
