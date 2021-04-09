#!/usr/bin/env python
import numpy as np
import common
import sys

def pixel2dat(filein,fileout):
    
    f_in=open(str(filein), 'r')
    input_lines = f_in.readlines()
    data=[]
    for line in input_lines:
        data.append([float(w) for w in line.split(",")[:]])
    f_in.close()    
    temp=np.array(data).flatten()
    #square_data = common.readpuredata(filein)
    
    dim=np.size(temp)
    ntria = 2*dim
    tria_data = np.zeros(ntria)
    i=0
    for pixel_value in temp:
        tria_data[i]=pixel_value
        i=i+1
        tria_data[i]=pixel_value
        i=i+1
    common.writedata(fileout,tria_data) 

if __name__ == "__main__":
    print len(sys.argv)
    if len(sys.argv) > 1:
        pixel2dat(*sys.argv[1:])
    else:
        raise SystemExit("usage:  python pixel2dat.py ndivx ndivy <fileout>")
    


