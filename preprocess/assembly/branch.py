# -*- coding: utf-8 -*-
"""
Toolbox for branch transport computation

"""

#!/usr/bin/env python
import re
import numpy as np

def y_branch(p,m,alpha):
    points=np.array(p)
    masses=np.array(m)
    

    coord_O=points[0,:]
    coord_Q=points[1,:]
    coord_P=points[2,:]


    m_O=abs(masses[0])
    m_Q=abs(masses[1])
    m_P=abs(masses[2])

    OP=coord_P-coord_O
    OQ=coord_Q-coord_O
    QP=coord_P-coord_Q

    OQP=np.arccos( np.dot(-OQ, QP)/ ( np.sqrt(np.dot(OQ,OQ)) * np.sqrt(np.dot(QP,QP) )))
    QPO=np.arccos( np.dot(-OP,-QP)/ ( np.sqrt(np.dot(OP,OP)) * np.sqrt(np.dot(QP,QP) )))
    POQ=np.arccos( np.dot( OQ, OP)/ ( np.sqrt(np.dot(OQ,OQ)) * np.sqrt(np.dot(OP,OP) )))
    
    k_1=(m_P/m_O)**(2*alpha)
    k_2=(m_Q/m_O)**(2*alpha)

    theta_1=np.arccos( (k_2-k_1-1)/(2*np.sqrt(k_1))     )
    theta_2=np.arccos( (k_1-k_2-1)/(2*np.sqrt(k_2))     )
    theta_3=np.arccos( (1-k_1-k_2)/(2*np.sqrt(k_1*k_2)) )
    
    v=[];
    if (POQ>=theta_3):
        B_opt=coord_O
        v.append([0,1])
        v.append([0,2])
    elif ( (OQP>=theta_1) & (POQ<theta_3)):
        B_opt=coord_Q
        v[:,0]=[0,1]
        v[:,1]=[1,2]
    elif ( (QPO>=theta_2) & (POQ<theta_3)):
        B_opt=coord_P
        v[:,0]=[0,2]
        v[:,1]=[1,2]
    else:
        QM=np.dot(OP,OQ)/np.dot(OP,OP) * OP - OQ
        PH=np.dot(OP,OQ)/np.dot(OQ,OQ) * OQ - OP
        
        R=(coord_O+coord_P)/2.0 - (np.cos(theta_1)/np.sin(theta_1))/2.0 * np.sqrt(np.dot(OP,OP)/ np.dot(QM,QM) )* QM
        S=(coord_O+coord_Q)/2.0 - (np.cos(theta_2)/np.sin(theta_2))/2.0 * np.sqrt(np.dot(OQ,OQ)/ np.dot(PH,PH) ) * PH
        RO=coord_O-R
        RS=S-R
        
        B_opt=2*( (1-np.dot(RO,RS) / np.dot(RS,RS)) *R + np.dot(RO,RS)/np.dot(RS,RS)*S)-coord_O
        
        #p_B=tuple([tuple(i) for i in B_opt])
        #p_B=[i for i in enumerate(B_opt)]
        p_B=B_opt.tolist()

        p.append(p_B)
        v.append([0,3])
        v.append([1,3])
        v.append([2,3])
        
    return p,v;

p=[[0.5,0.1],[0.4,0.9],[0.6,0.9]]
m=[1.0,0.5,0.5]
alpha=0.0

p, v = y_branch(p,m,alpha)
