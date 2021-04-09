#
import numpy as np
import os



nrotate = 40

for i in range(0,nrotate):
    ang=2*np.pi*i/nrotate
    x=-0.6*np.cos(ang)
    y=0.6*np.sin(ang)
    # Read in the file
    with open('base00.session', 'r') as file :
        filedata = file.read()

    # Replace the target string
    filedata = filedata.replace('xnorm', str(x))
    filedata = filedata.replace('ynorm', str(y))
    filedata = filedata.replace('angolo_0', 'angolo_'+str(i))

    # Write the file out again
    fname='rotation'+str(i)+'.session'
    with open(fname, 'w') as file:
        file.write(filedata)
        
    #command=('visit -cli -nowin -s restore_print.py '+ str(fname))
    #os.system(command)

