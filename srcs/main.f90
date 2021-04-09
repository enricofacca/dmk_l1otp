program simple
use KindDeclaration
use Globals
use AbstractGeometry
use TimeInputs  
use DmkControls
use TdensPotentialSystem
use DmkP1P0
use DmkInputsData
implicit none
type(abs_simplex_mesh) :: grid,grid0,subgrid
type(file) :: fgrid,fforcing
type(TimeData) :: forcing
type(DmkCtrl) :: ctrl
integer :: lun_err=6,i,info
integer, allocatable :: topol(:,:)

real(kind=double), allocatable :: forcing_cell(:),tdens(:),pot(:)
logical :: endfile
real(kind=double) :: pflux

type(file) :: fsubgrid
type(tdpotsys), target :: tdpot
!type(p1p0_space_discretization), target :: p1p0
type(DmkInputs), target :: inputs_data
character (len=512) :: pathgrid,pathforcing,input,fname
real(kind=double), allocatable :: forcing_subgrid(:)
real(kind=double), allocatable :: dirac(:)
logical :: rc
integer :: res
type(file) :: fmat
integer :: ntdens,npot

!
! get grid path
!
CALL getarg(1, input)
read(input,*,iostat=res) pathgrid
write(*,*) 'Grid from', etb(pathgrid)

! 
! get forcing path
!
CALL getarg(2, input)
read(input,*,iostat=res) pathforcing
write(*,*) 'Forcing from', etb(pathforcing)

!
! get pflux
!
CALL getarg(3, input)
read(input,*,iostat=res) pflux
write(*,*) 'Pflux =', pflux

!
! read grid
! 
call fgrid%init(lun_err,'grid.dat',10,'in')
call grid%read_mesh(lun_err,fgrid)
call fgrid%kill(lun_err)

!
! read forcing
!
call fforcing%init(lun_err,'forcing.dat',11,'in')
call forcing%init(lun_err, fforcing, 1,grid%ncell)
call forcing%set(lun_err, fforcing, 0.0d0,endfile)
allocate(forcing_cell(grid%ncell),tdens(grid%ncell),pot(grid%nnode),topol(3,grid%ncell))
topol=grid%topol(1:3,1:grid%ncell)
do i=1,grid%ncell
   forcing_cell(i) =forcing%TDactual(1,i)
end do


!
! using optdmk
!
!tdens=one
!pot=zero
!call otpdmk(grid%nnodeincell,grid%nnode,grid%ncell,&
!     topol,grid%coord,&
!     pflux,forcing_cell,&
!     tdens,pot,&
!     ctrl,info)

!
! using steay
!
call data2grids(&
     lun_err,&
     grid%nnodeincell,grid%nnode,grid%ncell, &
     grid%coord,topol,&
     grid0,subgrid)

ntdens=grid%ncell
npot=subgrid%nnode

!
! inputs data
!
call inputs_data%init(lun_err, ntdens,npot,set2default=.True.)
inputs_data%pflux=pflux

!
! init set controls
!
! globals controls
ctrl%selection=0
ctrl%threshold_tdens=1.d-10
ctrl%debug=1
ctrl%min_tdens = 1.0d-13


! linear solver ctrl
ctrl%outer_solver_approach='ITERATIVE'
ctrl%outer_solver_approach='MG'
ctrl%outer_krylov_scheme='PCG'
ctrl%outer_prec_type='MG'
ctrl%outer_prec_n_fillin=30
ctrl%outer_prec_tol_fillin=1d-3
ctrl%relax4prec=0d-09
ctrl%relax_direct=1d-09
ctrl%iprt=1

ctrl%outer_imax=500
ctrl%outer_iprt=0
ctrl%tolerance_linear_solver=1d-12
ctrl%tolerance_linear_newton=1d-05
ctrl%inner_imax=5

! time  ctrl
! evolution controls
ctrl%max_time_iterations=1
ctrl%tolerance_system_variation=1d-04
! time stepping scheme
ctrl%time_discretization_scheme=7
! time stepping scheme controls
ctrl%newton_method=100
ctrl%max_nonlinear_iterations=20
ctrl%max_restart_update=10
ctrl%epsilon_W22=1d-8

! time step size controls
ctrl%deltat = 1
ctrl%deltat_control=2
ctrl%deltat_expansion_rate=2.0d0
ctrl%deltat_lower_bound=0.01
ctrl%deltat_upper_bound=100


! info , saving ctrl
ctrl%info_update=3
ctrl%info_state=2
ctrl%id_save_dat=3
ctrl%lun_statistics=10
ctrl%lun_out=6
ctrl%lun_tdens = 1234
ctrl%fn_tdens='tdens.dat'
ctrl%lun_pot=1235
ctrl%fn_pot='pot.dat'


ctrl%max_nonlinear_iterations=20
ctrl%tolerance_nonlinear=1d-09
ctrl%relax4prec=0.0d-10
ctrl%tolerance_system_variation=1d-10


allocate(&
     dirac(grid%nnode),&
     stat=res)
if (res.ne.0) rc = IOerr(lun_err, err_alloc, 'otpdmk', &
     ' work arrays forcing_subgrid dirac_subgrid ', res)


dirac=zero
call build_subgrid_rhs(subgrid,inputs_data%rhs,forcing_cell, dirac)

deallocate(&
     dirac,&
     stat=res)
if (res.ne.0) rc = IOerr(lun_err, err_dealloc, 'otpdmk', &
     ' work arrays forcing_subgrid dirac_subgrid ', res)
fname='rhs.dat'
call fmat%init(0,fname,10000,'out')
call write_steady(0, 10000, npot,inputs_data%rhs )
call fmat%kill(0)


!
! init tdpot
!
call tdpot%init( lun_err, ntdens, npot,1)
tdpot%tdens=one
tdpot%pot=zero



call dmkp1p0_steady_data(&
     grid0,subgrid,tdpot,inputs_data,&
     ctrl,info)


end program simple
