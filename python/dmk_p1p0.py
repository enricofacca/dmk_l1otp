#setup 
import numpy as np
import sys
import os

# Import I/O for timedata
try:
    sys.path.append('../../../../globals/python/timedata/')
    import timedata as td
except:
    print("Global repo non found")


#
# import fortran_pyhton_interface library 
#
current_source_dir=os.path.dirname(os.path.realpath(__file__))
relative_libpath='../../build/python/fortran_python_interface/'
dmk_lib_path=os.path.abspath(os.path.normpath(current_source_dir+'/'+relative_libpath))
sys.path.append(dmk_lib_path)
print(dmk_lib_path)
try:
    from dmk import (
    Dmkcontrols,          # controls for dmk simulations
    Abstractgeometry,     # grids and geometry
    data2grids,           # two grids initializer
    Dmkinputsdata,        # inputs data ( rhs, kappas,beta, etc)
    build_subgrid_rhs,    # integrate forcing term subroutine 
    Tdenspotentialsystem, # input/output result tdens, pot
    dmkp1p0_steady_data,   # main interface subroutine
    Timefunctionals # information of time/algorithm evolution)
    )
except:
    print("Compile with:")
    print("mkdir build; cd build; cmake f2py_interface_dmk../; make ")

#
# Compute Optimal Transport Problem for cost equal to Eucledain distance
# INPUTS:
#  topol, coord : mesh topology  and coordinates
#  forcing      : function f^+,f^- valued on mesh cells
#  tolerance    : tolerances to achieved (measured as tdens variation)
#  [ctrl]       : dmk controls
#  [tdens0]     : initial guess for optimal transport density

def solve_MongeKantorovichEquations(topol,coord,forcing,tolerance=1e-5,ctrl=None,tdens0=None):
    print(forcing.shape)
    return solve_MinFluxProblem(topol,coord,forcing,pflux=1.0,tolerance=1e-5,ctrl=ctrl,tdens0=tdens0);

#
# Run DMK dynamics for different exponent pflux
#
def solve_MinFluxProblem(topol,coord,forcing,pflux=1.0,tolerance=1e-5,
                         ctrl=None,
                         tdens0=None):
    #
    # set controls
    #
    if ( ctrl is None):
        ctrl = Dmkcontrols.DmkCtrl()
        Dmkcontrols.get_from_file(ctrl,'default_dmk.ctrl')
        ctrl.tolerance_system_variation=tolerance
        ctrl.fn_tdens='tdens.dat'
        ctrl.lun_tdens=0
        ctrl.fn_pot='pot.dat'
        ctrl.lun_pot=0
        ctrl.fn_statistics='dmk.log'
        ctrl.lun_statistics=13
    else:
        ctrl.tolerance_system_variation=tolerance


    # Init. meshes
    grid, subgrid = init_geometry(topol, coord, ctrl.id_subgrid)
    ntdens=grid.ncell
    npot=subgrid.nnode

    # Init tdens-pot variable
    tdpot=Tdenspotentialsystem.tdpotsys()
    Tdenspotentialsystem.tdpotsys_constructor(tdpot,0,ntdens, npot,1)
    if ( tdens0 is None ):
        tdpot.tdens[:]= 1.0
    else:
        tdpot.tdens[:]=tdens0[:]
    tdpot.pot[:]=0.0

    # Set inputs
    dmkin=Dmkinputsdata.DmkInputs()
    Dmkinputsdata.dmkinputs_constructor(dmkin,0,ntdens,npot,True)
    # integrate forcing term w.r.t. p1 base function
    build_subgrid_rhs(subgrid, dmkin.rhs, forcing, np.zeros(grid.nnode))
    dmkin.pflux = pflux


    #
    # init type for storing evolution/algorithm info
    #
    timefun=Timefunctionals.evolfun()
    Timefunctionals.evolfun_constructor(timefun, 0,
                                        ctrl.max_time_iterations,
                                        ctrl.max_nonlinear_iterations)
    
    # solve with dmk
    # solve with dmk
    info=0
    dmkp1p0_steady_data(grid, subgrid, tdpot, dmkin, ctrl, info,timefun=timefun)

    # Copy before freeing memory.
    # Use np.array() function because you will get errors otherwise
    tdens = np.array(tdpot.tdens)
    pot   = np.array(tdpot.pot)
    topol_subgrid=np.array(subgrid.topol)
    topol_coord=np.array(subgrid.topol)
    
    # free memory
    Tdenspotentialsystem.tdpotsys_destructor(tdpot, 0)
    Dmkinputsdata.dmkinputs_destructor(dmkin,0)
    Abstractgeometry.mesh_destructor(grid, 0)
    Abstractgeometry.mesh_destructor(subgrid, 0)
    
    return info, tdens,pot, topol_subgrid, topol_coord,timefun;

#
# Function initializing grids from raw data topol and coord
#
def init_geometry(topol, coord, id_subgrid):
    # Handle different shape and pyhton ordering 
    if ( topol.shape[0] > topol.shape[1]  ):
        topolfortran=topol.transpose()
    else: 
        topolfortran=topol
    if ( np.amin(topolfortran) == 0 ):
        topolfortran=topolfortran+1
        
    #
    # init grids
    # 
    if (id_subgrid == 1):
        grid=Abstractgeometry.abs_simplex_mesh()
        subgrid=Abstractgeometry.abs_simplex_mesh()
        coordT=coord.transpose()
        data2grids(6, 3, len(coord), len(topol), coordT, topolfortran, grid,subgrid)
    else:
        print('Use 2 level grids')

    return grid,subgrid;


#
# Run DMK dynamics for different exponent pflux
#
def solve_dmk(grid, subgrid, tdpot, dmkin, ctrl, info,timefun):
    # solve with dmk
    info=0
    dmkp1p0_steady_data(grid, subgrid, tdpot, dmkin, ctrl, info,timefun=timefun)

    return grid, subgrid, tdpot, dmkin, ctrl, info,timefun;
             


