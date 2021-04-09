# -*- coding: utf-8 -*-
#!/usr/bin/env python
import numpy as np
import sys
import os

#
# geometry
#
import meshtools as mt
import example_grid as ex_grid

#
# inputs
#
import example_forcing as ex_forcing
import common as common
import timecell as timecell

def define_optimal_potential(flag_grid, flag_source, flag_sink,flag_kappa,flag_dirichlet,
                             extra_source=None,extra_sink=None):
    #####################################################################
    # optimal transport density definition=
    #####################################################################
    if ( (flag_grid=='bouteiller18_tc1') and 
         (flag_source=='bouteiller18_tc1') and 
         (flag_sink=='bouteiller18_tc1')  ) :
        def optpot(time,coord,flag):
            x=coord[0]
            y=coord[1]
            z=coord[2]
            #
            steady=True
            Szero=2
            gzero=(0.0,0.5)
            gzeronorm=np.sqrt(gzero[0]**2+gzero[1]**2)
            xzero=(2.0,2.0)
            S=1/(1/Szero + np.scal(gzero,coord[0:1]-xzero))
            
            f=1/Gzeronorm * np.arccosh(
                1+0.5*S*Szero*Gzeronrm**2*((x-xzero[0])**2+(y-xzer0[1])**2))
            return f,steady;

    if ( (flag_grid=='bouteiller18_tc1') and
         (flag_source=='no') and 
         (flag_sink=='bouteiller18_tc1') and
         (flag_kappa=='bouteiller18_tc1') ) :
        def optpot(time,coord,flag):
            x=coord[0]
            y=coord[1]
            z=coord[2]
            #
            steady=True    
            Szero=2
            gzero=(0.0,0.5)
            gzeronorm=np.sqrt(gzero[0]**2+gzero[1]**2)
            xzero=(2.0,2.0)
            S=1/(1/Szero + np.scal(gzero,coord[0:1]-xzero))
            
            f=1/Gzeronorm * np.arccosh(
                1+0.5*S*Szero*Gzeronrm**2*((x-xzero[0])**2+(y-xzer0[1])**2))
            return f,steady;

        return optpot;
    #
    # Equation 5.9
    # Finite element approximation of the $p$-Laplacian
    # Barrett, John W.Liu, W. B.
    if ( (flag_grid  =='plaplacain') and
         (flag_source=='plaplacian_dirichlet') and 
         (flag_sink  =='plaplacian_dirichlet') and
         (flag_dirichlet=='zero') ):
        plapl=float(extra_source)
        def optpot(time,coord,flag):
            x=coord[0]
            y=coord[1]
            z=coord[2]
            #
            r=np.srt(x**2+y**2)
            steady=True
            sigma=0
            p=plapl
            f=(p-1)/[1/(sigma+2)]**(1/(p-2))*[1-r**((sigma+p)/(p-1))]/(sigma+p)

            return f,steady;

        return optpot;

    
    
    if ( (flag_source=='zero') and 
         (flag_dirichlet=='hole_square') and 
         (flag_sink=='uniform') and 
         (float(flag_kappa)==1.0 ) ) :
        print ('ok')
        def optpot(time,coord,flag):
            x=coord[0]
            y=coord[1]
            z=coord[2]
            #
            steady=True
            f=np.sqrt((x)**2+(y)**2)
            return f,steady;
    

        return optpot;        

if __name__ == "__main__":
    if len(sys.argv) > 1:
        fin_ctrl=sys.argv[1]  
        fin_grid=sys.argv[2]
        fout_optpot=sys.argv[3]
        
        ############################################
        # reads flags and other values
        ############################################
        ctrl = open(fin_ctrl, "r")
        ctrl_lines = ctrl.readlines()
        ctrl.close()
        
        # flags grid
        flag_grid=str(common.read_column('flag_grid',0,ctrl_lines))
        extra_grid=str(common.read_column('flag_grid',1,ctrl_lines))
        
        # flag_source
        flag_source=str(common.read_column('flag_source',0,ctrl_lines))
        extra_source=str(common.read_column('flag_source',1,ctrl_lines))
        
        #flag sink
        flag_sink=str(common.read_column('flag_sink',0,ctrl_lines))
        extra_sink=str(common.read_column('flag_sink',1,ctrl_lines))

        #flag sink
        flag_kappa=str(common.read_column('flag_kappa',0,ctrl_lines))
        extra_kappa=str(common.read_column('flag_kappa',1,ctrl_lines))

        #flag sink
        flag_dirichlet=str(common.read_column('flag_dirichlet',0,ctrl_lines))
        extra_dirichlet=str(common.read_column('flag_dirichlet',1,ctrl_lines))


   
        opt_potential_functions=[]
        opt_potential_functions.append(
            define_optimal_potential(flag_grid, flag_source, flag_sink,flag_kappa,flag_dirichlet))
        if( len(opt_potential_functions) == 1):
            # reads coord topol and flags 
            coord, topol, flags = mt.read_grid(fin_grid)
            if ( coord.shape[1] == 2 ):
                zcoord = np.zeros([coord.shape[0],1])
                coord=np.append(coord, zcoord, axis=1)
            nnode=len(coord)
            ncell=len(topol)
            opt_potential_node = np.zeros([nnode,1])
            for inode in range(nnode):
                f,steady=opt_potential_functions[0](0.0,coord[inode],flags[inode])
                opt_potential_node[inode]=f

            
            #
            # imposed zero average 
            #
            size_cell=mt.make_size(coord,topol)
            size_node=np.zeros(nnode)
            for icell in range(ncell):
                nodes=topol[icell,:]
                for i in range(3):
                    size_node[nodes[i]]=size_node[nodes[i]]+size_cell[icell]/float(3)


            #scal=np.dot(opt_potential_node[:,0],size_node[:])
            #opt_potential_node[:]=opt_potential_node[:]-scal/np.sum(size_node)
            #print(np.dot(opt_potential_node[:,0],size_node[:]))
            
            # write to file
            file_out = open(fout_optpot, "w")
            timecell.write2file(file_out,0.0,True,opt_potential_node,steady=False)
            timecell.write2file(file_out,0.0,False,opt_potential_node,steady=True)
            file_out.close()

        else:
            print('Optimal transport density not defined')
    else:
        raise SystemExit("usage:  python make_optimal_potential.py 'input ctrl file' 'grid file' 'output file'"  )
