#!/usr/bin/env python
import numpy as np
import sys

def dat2pixel(filein,fileout):
    f_in=open(str(filein), 'r')
    input_lines = f_in.readlines()
    size=int(input_lines[0].split()[0])
    ndata=int(input_lines[0].split()[1])
    ndiv=int(np.sqrt(ndata/2))
    data=[]
    print ndiv
    for i in range(ndata/2):
        left  =float(input_lines[1+2*i])
        right=float(input_lines[1+2*i+1])
        data.append((left+right)/2.0)
    f_in.close() 
                     
    f_out=open(str(fileout), 'w')
    for i in range(ndiv):
        for j in range(ndiv-1):
            f_out.write(str(data[i*ndiv+j])+",")
        f_out.write(str(data[i*ndiv+ndiv-1])+"\n")
    f_out.close()


if __name__ == "__main__":
    print len(sys.argv)
    if len(sys.argv) > 1:
        dat2pixel(*sys.argv[1:])
    else:
        raise SystemExit("usage:  python dat2pixel.py <filein fileout>")
    


