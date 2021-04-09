#!/usr/bin/env python
import re
import sys
import numpy as np
import fileinput



# function to erase comments after sep character from a string
def remove_comments(line, sep):
    line=str(line)
    for s in sep:
        i = line.find(s)
        if i >= 0:
            line = line[:i]
    return line.strip()


#  function to read data from file
def readdata(filepath,data):
    ndata=len(data)
    filein=open(str(filepath), 'r')
    input_lines = filein.readlines()
    if len(input_lines) == ndata:
        i=0
        for line in input_lines:
            data[i]=float(line)
            i=i+1
    else:
        print ('Dimension mismatch-Ndata=',len(input_lines),'len(array)=',ndata)
        filein.close()   
    return data



#  function to read data from file
def readmydata(filepath):
    filein=open(str(filepath), 'r')
    input_lines = filein.readlines()
    size_data=int(input_lines[0].split()[0])
    len_data=int(input_lines[0].split()[1])
    data=np.zeros([len_data,size_data])
    for i in range(len_data):
        data[i][:]=[float(w) for w in input_lines[i+1].split()[:]]
    filein.close()    
    return data;


#  function to read data from file
def writedata(filepath,data):
    fileout=open(str(filepath), 'w')
    for values in data:
        fileout.write(str(values)+"\n")
    fileout.close()

#  function to read data from file
def writemydata(filepath,data):
    fileout=open(str(filepath), 'w')
    fileout.write(str(data.shape[1])+" "+str(data.shape[0])+"\n")
    for i in range(data.shape[0]):
        fileout.write(" ".join(map(str,data[i][:]))  +"\n")
    fileout.close()




# funnction to read the ncol^{th}-string from the line containg
# the flag string
def read_column(flag,ncol,lines):
    stringa=''
    for line in lines:
        if str(flag) in line:
            clean=remove_comments(line,'!')
            stringa=clean.split()[ncol]
    if (stringa == ''):
        print ('Wrong Flag '+ flag +' not found')
    return stringa;

def replace(flag,filein, index_column, value):
    #fileinput = open(filein, "rw")
    for line in fileinput.input(str(filein), inplace=True): 
        line_str=line.rstrip().split()
        if ( flag in line_str ):
            line_str[index_column]=value
        print (" ".join([str(i) for i in line_str]))
    #fileinput.close()


#define time range and time stepping
deltat=0.1
def next_time(time):
    time=time+deltat
    return time;


