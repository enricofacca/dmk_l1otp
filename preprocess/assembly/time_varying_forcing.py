#!/usr/bin/env python
import numpy as np
import meshtools as mt
import sys

def build_dirac(locations, fluxes, grid, rhs ):
    #
    # read stations coordinates
    #
    input_file=open(locations,'r')
    input_lines = input_file.readlines()
    nstations=len(input_lines)
    coord_location=[];
    for line in input_lines:
        x = float(line.split()[0])
        y = float(line.split()[1])
        coord_location.append([x,y])
    input_file.close()

    #
    # read grid use for discrtization
    #
    
    coord, topol, flags = mt.read_grid(grid)

    index_stations=[];
    for point in coord_location:
        inode,dist=mt.FindClosestNode(
                range(len(coord)),coord,point)
        index_stations.append(int(inode[0]))
        
        



    #
    # read fluxess coordinates and write data into rhs
    #
    
    

    #
    # read all time sequence
    #
    input_file=open(fluxes,'r')
    input_lines = input_file.readlines()
    input_file.close()
    
    #
    # convert it into rhs
    #
    
    # write head timedata file
    file_out = open(rhs, "w")
    file_out.write(
        str(1) + " " + 
        str(len(coord)) + " !  dim ndata" 
        +"\n")
    #
    # write rhs sequentally
    #
    for line in input_lines:        
        time = float(line.split()[0])
        file_out.write(
            'time '+str(time)+"\n")
        file_out.write(
            str(nstations) +"\n")
                
        current_fluxes = np.array([float(w) for w in line.split()[1:nstations+1]])
        j=0
        for i in index_stations:
            file_out.write(str(i+1)+" " + str(current_fluxes[j])+" \n")
            j=j+1
    #
    # if only one time is passed, set data into steady state
    #       
    if (len(input:lines) == 1 ):
        file_out.write( "time 1.0e30 \n")
        
    file_out.close()


if __name__ == "__main__":
    if len(sys.argv) > 1:
        build_dirac(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4] )

    
    
    
    
    


         
         


         


