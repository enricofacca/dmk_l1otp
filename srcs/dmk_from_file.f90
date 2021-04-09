!>---------------------------------------------------------------------
!> 
!> MUFFE for the solution of Optimal Transport problems
!>
!> \author{Enrico Facca and Mario Putti}
!>
!> DESCRIPTION: 
!> The main program 
!> 
!>
!> REVISION HISTORY:
!> 22 Lug 2016 - Initial Version
!>
!> TODO_dd_mmm_yyyy - TODO_describe_appropriate_changes - TODO_name
!<---------------------------------------------------------------------

PROGRAM dmkfromfile 
  use Globals
  use Timing
  use TimeInputs
  use IOdefs
  use DmkOdeData
  use DataSequence
  use ControlParameters
  use MuffeTiming
  use Geometry2d
  use AbstractGeometry
  use P1Galerkin
  use TdensPotentialSystem
  use TimeFunctionals
  use SparseMatrix
  use StdSparsePrec
  use Matrix
  use LinearSolver
  use StdSparsePrec
  use DmkDiscretization
  use DmkInputsData
  use DmkP1P0
  use DmkControls


  implicit none
  
  type(IOfdescr) :: IOfiles
  type(DmkCtrl)  :: ctrl_read,ctrl
  type(codeTim)  :: CPU
  type(abs_simplex_mesh)     :: grid,subgrid
  type(tdpotsys) :: tdpot
  type(evolfun)  :: timefun 
  type(output_solver)  :: info_solver
  type(Tim) :: wasted_temp
  type(file) :: file_out,fout
  integer  :: info_prec, passed_reduced_jacobian
  type(OdeInp) :: inputs_from_file
  type(DmkInputs) :: inputs_data
  type(p1p0_space_discretization) :: p1p0
  
  
  ! varaible for projection
  type(p1gal) :: p1_grid
  type(spmat) :: stiff_grid
  real(kind=double), allocatable :: pot_grid(:), rhs_grid(:)
  type(TimeData) :: rhs_integrated
  type(stdprec) :: prec_stiff

  logical :: rc
  logical :: test_var
  logical :: test_maxit
  logical :: in_cycle
  logical :: endfile
  
  integer :: lun_err,lun_out,lun,res,lun_stat
  integer :: i, inode_sub, i1,i2,itemp,k,iloc,jloc
  integer :: nrestart, nrestart_max
  integer :: int_before_dat,  int_before_matrix
  integer :: max_time_iterations
  integer :: nupdate

  integer :: itemp_prec_calc
  integer :: total_kryvol=0, wasted_krylov=0
  integer :: total_update=0

  
  integer :: info,info_update

  integer :: npot, ntdens

  real(kind=double)              :: deltat,plapl
  real(kind=double), allocatable :: err_tdens(:)


  integer, parameter :: nsequence=2

  real(kind=double) :: ddot,dnrm2

  real(kind=double) :: tol_zero,pode,current_time,w1
  integer  :: ndir=1,flag
  integer,allocatable:: noddir(:)
  
    

  character(len=256) :: fname,folder,msg,msg2
  character(len=256) :: out_format,str
  
  !
  ! Initialization of Timing module
  !
  call CPU%init()
  call CPU%TOT%set('start')
  call CPU%OVH%set('start')

   
  !
  ! Initialization of INPUT_OUTPUT module
  !
  ctrl_read%debug=1
  ctrl_read%tolerance_linear_newton=1.0d-4
  if (ctrl_read%debug .eq. 1 ) then
     write(*,*) 'Init I/O from muffa.fnames'
  end if
  call IOfiles%init()
  lun_err  = IOfiles%stderr%lun
  lun_out  = IOfiles%stdout%lun
  lun_stat = IOfiles%statistic%lun

  !
  ! 1 -Init data structure with data stored in file
  !
  !
  ! 1.1 - Read controls
  !
  write(lun_stat,'(a)') ctrl_read%separator('controls simulation')
  call ctrl_read%readfromfilenew(IOfiles%controls%lun)
  ! set output file
  ctrl_read%lun_tdens=IOfiles%tdens_out%lun
  ctrl_read%fn_tdens=IOfiles%tdens_out%fn
  ctrl_read%lun_pot=IOfiles%pot_out%lun
  ctrl_read%fn_pot =IOfiles%pot_out%fn
  ctrl_read%lun_statistics=IOfiles%statistic%lun
  ctrl_read%fn_statistics =IOfiles%statistic%fn
  ! create a copy 
  ctrl = ctrl_read
  
  
  !
  ! 1.2 - Read grid
  !
  write(lun_stat,'(a)') ctrl_read%separator('Setup data begin')
  write(lun_stat,'(a)') ctrl_read%separator('Read grid')
  call grid%read_mesh(lun_err,IOfiles%grid)
  call grid%build_size_cell(lun_err)
  call grid%build_normal_cell(lun_err)
  
  !
  ! 1.3 - Read subgrid
  !
  write(lun_stat,'(a)') ctrl_read%separator('Read subgrid')
  call subgrid%read_mesh(lun_err,IOfiles%subgrid)
  call subgrid%read_parent(lun_err,IOfiles%parents)
  call subgrid%build_size_cell(lun_err)
  call subgrid%build_normal_cell(lun_err)

  
  write(lun_stat,'(a)') 'init inputs data'
  if (ctrl_read%debug .eq. 1) write(IOfiles%stdout%lun,'(a)') 'Init inputs data'

  if (  ctrl_read%id_subgrid .eq. 1) then
     ntdens = grid%ncell
     npot   = subgrid%nnode
  else
     ntdens = grid%ncell
     npot   = grid%nnode
  end if
  
  !
  ! 1.4 - Initialization inputs data from file
  !
  call inputs_from_file%init(&
       Iofiles,lun_out,&
       ntdens,npot,ctrl%id_subgrid,&
       ctrl%tzero)


  ctrl%lun_err=IOfiles%stderr%lun
  ctrl%lun_out=IOfiles%stdout%lun
  ctrl%lun_statistics=IOfiles%statistic%lun

  !
  ! 2 - Initializtion of data structure or the soultion of DMK
  !

  !
  ! 2.1 - Spatial discretization
  !
  if ( ctrl_read%id_subgrid == 1 ) then       
     call p1p0%init(&         
          ctrl_read,&
          ctrl_read%id_subgrid,grid, subgrid)
  else
     call p1p0%init(&
          ctrl_read,&
          ctrl_read%id_subgrid,grid, grid)
  end if
  if (ctrl_read%id_save_matrix > 0) call save_assembler(p1p0,IOfiles)

  
  !
  ! 2.2 - Data struture to pass inputs to dmk solver 
  !
  call inputs_data%init(lun_err, ntdens,npot)  
  allocate(&
       err_tdens(ntdens),&
       stat=res)
  if(res .ne. 0) rc = IOerr(lun_err, err_alloc , 'dmkfromfile', &
       ' temp array error_tdens ',res)
  if (IOfiles%opttdens%exist) then
     do i=1,ntdens
        err_tdens(i)=inputs_from_file%opttdens%TdActual(1,i)
     end do   
     call setopttdens(inputs_data,lun_err,err_tdens)
  end if
  


  !
  ! 2.3 -Initialized tdens-potential type
  !
  call tdpot%init( lun_err,ntdens, npot)
  do i=1,ntdens
     tdpot%tdens(i)=inputs_from_file%tdens0%TDactual(1,i)
  end do
  
  !
  ! 2.4 - Allocate array time, var_tdens and time functional
  !
  call timefun%init(lun_err,&
       ctrl_read%max_time_iterations, ctrl_read%max_nonlinear_iterations)
  
  msg=ctrl_read%separator('Setup data ended')
  if ( ( ctrl%info_state .gt. 1).and. (lun_out>0) )  write(lun_out,'(a)') etb(msg)
  if ( (lun_stat>0) ) write(lun_stat,'(a)') etb(msg)

  !------------------------------------------------------------------------------------ 
  !
  ! 3 Time evolution starts
  !
  msg=ctrl_read%separator('Evolution begins')
  if ( ( ctrl%info_state .gt. 1).and. (lun_out>0) )  write(lun_out,'(a)') etb(msg)
  if ( (lun_stat>0) ) write(lun_stat,'(a)') etb(msg)

  !
  ! run dmk time evolution
  !
  !
  ! start time cycle
  !
  info=0
  flag=1
  current_time=inputs_data%time
  ctrl%debug=1
  do while ( flag.ne.0)
     call p1p0%dmk_cycle_reverse_communication(&
          flag,info,current_time,& 
          inputs_data,tdpot,ctrl,ctrl_read)
     write(*,*) 'info', info
     select case (flag)
     case(2)
        !
        ! fill ode_inputs with data at current time
        ! reading data from files and copy them into inputs_data
        call inputs_from_file%set(IOfiles, current_time,endfile)
        call inputs_from_file%set_DmkInputs(inputs_data)        
     case(3)
        !
        ! flag==3 R
        ! Right before new update
        ! Here tdpot and inputs_data are syncronized 
        ! The user can compute print any information regarding
        ! the state of the system
        !
        if (inputs_data%pmass.ne. one) then
           plapl=2*inputs_data%pmass/(inputs_data%pmass-one)
           write(*,*) 'pmass ', inputs_data%pmass, ' Plapl =', plapl
           
           call p1p0%build_norm_grad_dyn(tdpot%pot,2.0d0,p1p0%norm_grad_dyn)
           write(*,*) minval(p1p0%norm_grad_dyn),'<|\Grad \Pot|^2 < ', maxval(p1p0%norm_grad_dyn)
           
           
           !
           ! tdens= beta |\grad \Pot|^{p-2} = one/ kappa^{1/(p-2)})|\grad \Pot|^{p-2} 
           !
           p1p0%scr_ntdens=one/(inputs_data%kappa**(plapl-two))
           write(*,*) minval( p1p0%scr_ntdens), ' beta  ', maxval( p1p0%scr_ntdens)
           err_tdens=p1p0%scr_ntdens*p1p0%norm_grad_dyn**((plapl-2.0d0)/2.0d0)-tdpot%tdens
         

           write(*,*) 'Err Plapl Strange Laplacian L1 norm= ', &
                grid%normp_cell(one,err_tdens)
           write(*,*) 'Err Plapl Strange Laplacian Eu norm = ',&
                dnrm2(ntdens,err_tdens,1)
           write(*,*) minval(err_tdens/inputs_data%beta), ' Err tdens/beta  ', maxval(err_tdens/inputs_data%beta)
           write(*,*) minval(err_tdens),                  ' Err tdens       ', maxval(err_tdens)

           p1p0%scr_ntdens = &
                tdpot%tdens*p1p0%norm_grad_dyn -  inputs_data%kappa**2*tdpot%tdens**inputs_data%pmass
           write(*,*) minval(p1p0%scr_ntdens), ' Mirrow grad   ', maxval(p1p0%scr_ntdens)

           p1p0%scr_ntdens = &
                p1p0%norm_grad_dyn -  inputs_data%kappa**2*tdpot%tdens**(inputs_data%pmass-one)
           write(*,*) minval(p1p0%scr_ntdens), ' Lyap   grad   ', maxval(p1p0%scr_ntdens)
           
           

           ! 1.1 - assembly matrix
           ! use beta |\grad pot|^p-2 + lambda
           p1p0%scr_ntdens=one/(inputs_data%kappa**(plapl-two)) ! this is beta
           p1p0%scr_ntdens = p1p0%scr_ntdens * p1p0%norm_grad_dyn**((plapl-two)/two)+ inputs_data%lambda 
           call p1p0%assembly_stiffness_matrix(lun_err,p1p0%scr_ntdens,stiff_grid)

           ! 1.2 set rhs  
           p1p0%rhs=inputs_data%rhs

           if ( inputs_data%ndir > 0 ) then
              !
              ! act on rhs
              !
              call p1p0%p1%dirichlet_bc(lun_err,&
                   stiff_grid,p1p0%rhs, p1p0%scr_npot,&
                   inputs_data%ndir,&
                   inputs_data%dirichlet_nodes,&
                   inputs_data%dirichlet_values)
           end if

           call stiff_grid%Mxv(tdpot%pot,p1p0%scr_npot)
           p1p0%scr_npot=p1p0%scr_npot-p1p0%rhs





           write(*,*) 'error A(lambda + tdens) pot = rhs', dnrm2(npot,p1p0%scr_npot,1)
           write(*,*) 'error A(lambda + tdens) pot = rhs/|rhs|', &
                dnrm2(npot,p1p0%scr_npot,1)/dnrm2(npot,p1p0%rhs,1)

           write(*,*) minval(p1p0%scr_npot),' error pot', maxval(p1p0%scr_npot)

        end if

        if (inputs_data%opttdens_exists) &
             write(*,*) 'error ', p1p0%grid_tdens%normp_cell(one,inputs_data%opttdens-tdpot%tdens)

        
        call p1p0%compute_functionals(tdpot,inputs_data,ctrl)
        write(msg,'(a)') ctrl%separator('INFO STATE') 
        if (ctrl%info_state .ge. 2) then
           if (lun_out>0) then
              write(lun_out,*)  ' '
              write(lun_out,*) etb(msg)
              call tdpot%info_functional(lun_out)
           end if
        end if
        if (lun_stat>0) then
           write(lun_stat,*)  ' '
           write(lun_stat,*) etb(msg)
           call tdpot%info_functional(lun_stat)
        end if

        call p1p0%set_evolfun(timefun,&
             tdpot%time_iteration,tdpot,inputs_data,ctrl)
        timefun%var_tdens(tdpot%time_iteration)= tdpot%system_variation
        

        if ( inputs_data%opttdens_exists) then
           plapl=2*inputs_data%pmass/(inputs_data%pmass-one)
           p1p0%norm_grad_dyn=p1p0%norm_grad_dyn**((plapl-two)/two)
           err_tdens=abs(p1p0%norm_grad_dyn-tdpot%tdens)
           write(*,*) 'Err Plapl = ', &
                grid%normp_cell(one,err_tdens),'pmass',inputs_data%pmass
           
           err_tdens=tdpot%tdens - inputs_data%opttdens
           write(msg2,*) 'Err tdens = ', &
                grid%normp_cell(one,err_tdens)/grid%normp_cell(one,inputs_data%opttdens)
           call p1p0%build_norm_grad_dyn(tdpot%pot,inputs_data%pode,p1p0%norm_grad_dyn)

           if (ctrl%info_state .ge. 2) then
              if (lun_out>0) then
                 write(lun_out,*)  ' '
                 write(lun_out,*) etb(msg2)

              end if
           end if
           if (lun_stat>0) then
              write(lun_stat,*)  ' '
              write(lun_stat,*) etb(msg2)
           end if
        end if
     end select
     !
     ! In case of negative info, break cycle, free memory 
     ! 
     if ( info.ne.0 )  flag=0
  end do


  !>------------------------------------------------------
  !> Save time, var_tdens, time functional
  !>------------------------------------------------------
  call CPU%SAVE%set('start')
  ! Saving other functional
  call timefun%write2dat(&
       lun_err,IOfiles%timefun%lun,IOfiles%timefun%fn,&
       tdpot%time_iteration,&
       ctrl%iformat,ctrl%rformat)
  !
  ! free memory
  !
  call CPU%OVH%set('start')
  call timefun%kill(lun_err)
  call tdpot%kill(lun_err)
  call inputs_data%kill(lun_err)
  deallocate(&
       err_tdens,&
       stat=res)
  if(res .ne. 0) rc = IOerr(lun_err, err_dealloc , 'dmkfromfile', &
         ' temp array error_tdens ',res)
  
  call grid%kill(lun_err)
  call subgrid%kill(lun_err)  
  call CPU%OVH%set('stop')

  

  !
  ! Print timing statistics
  !
  call CPU%TOT%set('stop')
  if (ctrl%info_state .gt. 1) then
     call CPU%info(lun_out, 'Global simulation times:')
  end if
  call CPU%info(lun_stat, 'Global simulation times:')
  call CPU%kill()
  call IOfiles%kill(lun_err)
  
  if ( info.ne.0 )  then
     rc = IOerr(lun_err, err_val, ' dmkfromfile ', &
          ' SIMULATION STOPPED - IN REVERSE COMMUNICATION CYCLE RECIVED NEGATIVE INFO',info)
  end if
   
end PROGRAM dmkfromfile
  




     

!>--------------------------------------------
!> Initialization of all varibles that
!> contains information stored in file
!>--------------------------------------------
subroutine init_p1p0_muffe_from_files(IOfiles, ctrl_read, odein, grid, subgrid)
  use Globals
  use IOdefs
  use ControlParameters
  use AbstractGeometry
  implicit none
  type(IOfdescr), intent(in   ) :: IOfiles
  type(CtrlPrm),  intent(inout) :: ctrl_read
  type(OdeInp),   intent(inout) :: odein
  type(abs_simplex_mesh),  intent(inout) :: grid
  type(abs_simplex_mesh),  intent(inout) :: subgrid
  ! local
  integer :: lun_err,lun_stat
  integer :: ntdens, npot
  
  lun_err  = IOfiles%stderr%lun
  lun_stat = IOfiles%statistic%lun
  
  !
  ! 1 - Read controls
  !
  call ctrl_read%readfromfile(IOfiles%controls%lun)
  write(lun_stat,'(a)') ctrl_read%separator('controls simulation')
  !call ctrl_read%info(IOfiles%statistic%lun)
  
  !
  ! 2 - Read geometry info
  !
  write(lun_stat,'(a)') ctrl_read%separator('setup data begin')
  ! 2.1 - Read grid
  write(lun_stat,'(a)') ctrl_read%separator('read grid')
  call grid%read_mesh(lun_err,IOfiles%grid)
  call grid%build_size_cell(lun_err)
  call grid%build_normal_cell(lun_err)
  
  
  ! 2.2 - Read subgrid
  write(lun_stat,'(a)') ctrl_read%separator('read subgrid')
  call subgrid%read_mesh(lun_err,IOfiles%subgrid)
  call subgrid%read_parent(lun_err,IOfiles%parents)
  call subgrid%build_size_cell(lun_err)
  call subgrid%build_normal_cell(lun_err)

  
  write(lun_stat,'(a)') 'init inputs data'
  if (ctrl_read%debug .eq. 1) write(IOfiles%stdout%lun,'(a)') 'init inputs data'

  if (  ctrl_read%id_subgrid .eq. 1) then
     ntdens = grid%ncell
     npot   = subgrid%nnode
  else
     ntdens = grid%ncell
     npot   = grid%nnode
  end if

  call odein%init(&
       IOfiles,ctrl_read%lun_statistics,&
       ntdens,npot,ctrl_read%id_subgrid,&
       ctrl_read%tzero)

end subroutine init_p1p0_muffe_from_files






subroutine kill_muffe( lun_err,&
     tdpot, grid, subgrid,  inputs_data, ctrl)
  use Globals
  use ControlParameters
  use AbstractGeometry
  use TdensPotentialSystem
  use DmkOdeData
  implicit none
  integer,                intent(in   ) :: lun_err
  type(tdpotsys),         intent(inout) :: tdpot
  type(abs_simplex_mesh), intent(inout) :: grid, subgrid
  type(OdeData),          intent(inout) :: inputs_data
  type(CtrlPrm),          intent(inout) :: ctrl
  
  ! Killing procedure
  call inputs_data%kill(lun_err)
  call tdpot%kill(lun_err)
  call grid%kill(lun_err)
  call subgrid%kill(lun_err)
  call ctrl%kill()
  

end subroutine kill_muffe

!!$subroutine main_subroutine(&
!!$     nnode,ncell,nnodeincell,&
!!$     coord, topol,&
!!$     ctrl_read,&
!!$     ntdens, tdens,&
!!$     pflux,pode,forcing,&
!!$     pot,&
!!$     nnode_pot,ncell_pot,coord_pot, topol_pot)
!!$  use Globals
!!$  use Timing
!!$  use DmkInputsData
!!$  use ControlParameters
!!$  use AbstractGeometry
!!$  use MuffeTiming
!!$  use TdensPotentialSystem
!!$  use TimeFunctionals
!!$  use DmkP1P0
!!$  implicit none
!!$  integer, intent(in) :: nnode
!!$  integer, intent(in) :: ncell
!!$  integer, intent(in) :: nnodeincell
!!$  real*8, intent(in)  :: coord(3,nnode)
!!$  integer, intent(in) :: topol(nnodeincell,ncell)
!!$  type(CtrlPrm),intent(in) :: ctrl_read
!!$  integer, intent(in) :: ntdens
!!$  real*8, intent(inout) :: tdens(ntdens)
!!$  real*8, intent(inout) :: pflux
!!$  real*8, intent(inout) :: pode
!!$  real*8, intent(inout) :: forcing(ntdens)
!!$  real*8, intent(inout) :: pot(4*nnode)
!!$  integer, intent(inout) :: nnode_pot
!!$  integer, intent(inout) :: ncell_pot
!!$  real*8, intent(inout)  :: coord_pot(3,4*nnode)
!!$  integer, intent(inout) :: topol_pot(nnodeincell,4*ncell)
!!$  ! local
!!$  type(codeTim) :: CPU
!!$  type(CtrlPrm) :: ctrl
!!$  type(abs_simplex_mesh)     :: grid,subgrid
!!$  type(tdpotsys) :: tdpot
!!$  type(evolfun)  :: timefun 
!!$  type(Tim) :: wasted_temp
!!$  type(p1p0_space_discretization) :: p1p0
!!$  type(DmkInputs)  :: inputs_data
!!$  logical :: rc
!!$  logical :: test_var
!!$  logical :: test_maxit
!!$  logical :: in_cycle
!!$  logical :: endfile
!!$  
!!$  integer :: lun_err=0,lun_out=6,lun,res,lun_stat=7
!!$  integer :: itemp,itemp_prec_calc
!!$  integer :: nrestart, nrestart_max,max_time_iterations,nupdate
!!$  integer :: info, info_update,npot
!!$  real(kind=double)              :: deltat
!!$  real(kind=double), allocatable :: time(:)
!!$  real(kind=double), allocatable :: var_tdens(:)
!!$  real(kind=double), allocatable :: var_tdens_linfty(:)
!!$  character(len=256) :: msg,str
!!$
!!$
!!$  !
!!$  ! Initialization of Timing module
!!$  !
!!$  call CPU%init()
!!$  call CPU%TOT%set('start')
!!$  call CPU%OVH%set('start')
!!$
!!$
!!$  !
!!$  ! 2 - Read geometry info
!!$  !
!!$  write(lun_stat,'(a)') ctrl_read%separator('setup data begin')
!!$  ! 2.1 - Read grid
!!$  call grid%init_from_data(lun_err,nnode, ncell, nnodeincell,'triangle', topol,coord)
!!$  call grid%build_size_cell(lun_err)
!!$  call grid%build_normal_cell(lun_err)
!!$  write(lun_stat,'(a)') ctrl_read%separator('Init. grid')
!!$
!!$  
!!$  ! 2.2 - Read subgrid
!!$  call subgrid%refine(lun_err,grid)
!!$  call subgrid%build_size_cell(lun_err)
!!$  call subgrid%build_normal_cell(lun_err)
!!$  write(lun_stat,'(a)') ctrl_read%separator('Subgrid Done')
!!$
!!$  
!!$  !
!!$  ! BUILD FEM SPACES
!!$  !
!!$  if (  ctrl_read%id_subgrid .eq. 1) then
!!$     npot   = subgrid%nnode
!!$     call p1p0%init(&
!!$          ctrl%lun_err,ctrl%lun_out,ctrl%lun_statistics,&
!!$          ctrl_read,&
!!$          ctrl_read%id_subgrid,grid, subgrid)
!!$  else
!!$     npot   = grid%nnode
!!$     call p1p0%init(&
!!$          ctrl%lun_err,ctrl%lun_out,ctrl%lun_statistics,&
!!$          ctrl_read,&
!!$          ctrl_read%id_subgrid,grid, grid)
!!$  end if
!!$  write(lun_stat,'(a)') 'SPATIAL DISCRETIZATION DONE'
!!$  if (ctrl_read%debug .eq. 1) write(lun_out,'(a)') 'SPATIAL DISCRETIZATION DONE'
!!$
!!$  !
!!$  ! build inputs data 
!!$  !
!!$  call inputs_data%init(lun_err,  p1p0%ntdens, p1p0%npot)
!!$  !
!!$  ! USER DEFINED
!!$  !
!!$  inputs_data%pflux=pflux
!!$  inputs_data%pode=pode
!!$  tdpot%tdens=tdens
!!$  if (  ctrl_read%id_subgrid .eq. 1) then
!!$     !
!!$     ! buils rhs integrated
!!$     !
!!$     call subgrid%proj_subgrid(forcing, p1p0%scr_nfull(1:subgrid%ncell))
!!$     call p1p0%p1%build_rhs_forcing(p1p0%scr_nfull(1:subgrid%ncell),inputs_data%rhs(:))
!!$     inputs_data%rhs=inputs_data%rhs(:)
!!$  else
!!$     !
!!$     ! buils rhs integrated
!!$     !
!!$     call inputs_data%init(lun_err,  p1p0%ntdens, p1p0%npot)
!!$     call p1p0%p1%build_rhs_forcing(forcing,inputs_data%rhs)
!!$     inputs_data%rhs=inputs_data%rhs
!!$  end if
!!$
!!$  !
!!$  ! DEFAULT OPTION
!!$  !
!!$  inputs_data%kappa=one
!!$  inputs_data%decay=one
!!$  inputs_data%pmass=one
!!$
!!$  !
!!$  ! initialized tdens-potential type
!!$  !
!!$  call tdpot%init(&
!!$       lun_err,lun_out,lun_stat,&
!!$       p1p0,ctrl_read)
!!$
!!$
!!$   !
!!$  ! 4 - Allocate array time, var_tdens and time functional
!!$  !
!!$  allocate(&
!!$       time(0:ctrl_read%max_time_iterations),&
!!$       var_tdens(0:ctrl_read%max_time_iterations),&
!!$       var_tdens_linfty(0:ctrl_read%max_time_iterations),&
!!$       stat=res)
!!$  if(res .ne. 0) rc = IOerr(lun_err, err_alloc , 'main', &
!!$         ' error alloc. time var_tdens',res)
!!$  time      = zero
!!$  var_tdens = zero
!!$  var_tdens(0) = huge
!!$  var_tdens_linfty = zero
!!$  var_tdens_linfty(0) = huge
!!$
!!$  call timefun%init(lun_err,&
!!$       ctrl_read%max_time_iterations, ctrl_read%max_nonlinear_iterations)
!!$  !
!!$  write(lun_stat,'(a)') ctrl_read%separator('setup data ended')
!!$  write(lun_out,'(a)') ctrl_read%separator('setup data ended')
!!$
!!$
!!$  !
!!$  ! create a work copy use to change orignal parameters
!!$  !
!!$  ctrl = ctrl_read
!!$  
!!$  !
!!$  ! Time evolution starts
!!$  !
!!$  write(lun_stat,'(a)') ctrl_read%separator('time evolution begins')
!!$  write(lun_out,'(a)') ctrl_read%separator('time evolution begins')
!!$  
!!$
!!$  !
!!$  ! Time=tzero
!!$  !
!!$  itemp=0
!!$  time(itemp)=ctrl%tzero
!!$  max_time_iterations=ctrl%max_time_iterations
!!$  in_cycle=.true.
!!$
!!$  !
!!$  ! Start evolution syncronazing tdens system and ode's inputs
!!$  !
!!$  !  
!!$  ! 1 - Eval ODE input varible at time $tzero$
!!$  ! Pflux, Pmass, Kappa, Decay and  Rhs_integrated kappa
!!$  !
!!$  msg =  ctrl%separator('Start Tdens/Pot syncronization at tzero}')
!!$  write(lun_stat,*) etb(msg)
!!$  write(lun_out,* ) etb(msg)
!!$  
!!$  !
!!$  ! Syncronize at tzero
!!$  !
!!$  tdpot%pot=zero
!!$  ctrl%build_prec=1
!!$  ctrl%ctrl_solver%scheme='PCG'
!!$  ctrl%ctrl_prec%prec_type='IC'
!!$  call tdpot%syncronize_at_tzero(info,&
!!$       lun_err,lun_out, lun_stat,&
!!$       p1p0,&
!!$       inputs_data,&
!!$       ctrl, &
!!$       CPU)
!!$  ctrl%ctrl_solver= ctrl_read%ctrl_solver
!!$  ctrl%ctrl_prec  = ctrl_read%ctrl_prec
!!$  call p1p0%evaluate_functionals(tdpot,inputs_data)
!!$  call timefun%set_evolfun(itemp,tdpot,p1p0,inputs_data,ctrl,CPU)
!!$  msg =  ctrl%separator('Tdens/Pot syncronized at tzero}')
!!$  write(lun_stat,*) etb(msg)
!!$  write(lun_out,* ) etb(msg)
!!$ 
!!$  !
!!$  ! Info: Time=tzero,  time_fucntionals (basic) 
!!$  !
!!$  write(lun_stat,'(a)')' '
!!$  call info_muffe(lun_stat,ctrl,itemp, time(itemp),zero, timefun,huge,huge,0,CPU)
!!$  call info_muffe(lun_out,ctrl,itemp, time(itemp),zero, timefun,huge,huge,0,CPU)
!!$  
!!$  !
!!$  ! Start time cycle
!!$  !
!!$  in_cycle    = ( max_time_iterations > 0) 
!!$  deltat      = ctrl%deltat
!!$  nupdate     = 0
!!$  
!!$  do while (in_cycle )
!!$     if ( ctrl%build_prec .eq. 1 ) itemp_prec_calc = itemp
!!$     !
!!$     ! Update spatial variables.
!!$     ! First save $\Tdens^\tstep$ and $\Pot^\tstep$, then 
!!$     ! Update sytem tdens-pot, and all variable inside
!!$     ! (gradients, ode input variable ) at next time iteration
!!$     ! After this sobroutine ALL spatial variables need to valued
!!$     ! at $t^{\tstepp}=t^{\tstep}+\Deltat$
!!$     !--------------------------------------------------------
!!$     itemp = itemp + 1
!!$     write(lun_out,'(a)') ' '
!!$     write(lun_stat,'(a)') ' '
!!$     write(lun_out,'(a)') ' '
!!$     write(lun_stat,'(a)') ' ' 
!!$     write(str, '(a,I5)') 'UPDATE begin: itemp= ', itemp
!!$     write(lun_stat,'(a)') ctrl%separator(etb(str))
!!$     write(lun_out,'(a)') ctrl%separator(etb(str))
!!$     
!!$     info_update = -1
!!$     nrestart = 0     
!!$     
!!$     if ( ctrl%debug .eq. 1) write(lun_out,*) 'Set controls before update'
!!$     !
!!$     ! set all controls for update system
!!$     !
!!$     call  set_controls_update (&
!!$          itemp,&
!!$          deltat,&
!!$          tdpot, &
!!$          p1p0,&
!!$          inputs_data,&
!!$          timefun,&
!!$          ctrl_read,&
!!$          ctrl)
!!$
!!$     do while ( ( info_update .ne. 0 ) .and. (nrestart < ctrl%nrestart_max) ) 
!!$        !
!!$        ! update system
!!$        !
!!$        wasted_temp = CPU%WASTED
!!$        call CPU%ALGORITHM%set('start')
!!$        call CPU%WASTED%set('start')
!!$        call tdpot%update(&
!!$             lun_err,lun_out,lun_stat,&
!!$             ctrl,&
!!$             deltat,&
!!$             itemp,&
!!$             time(itemp-1),&
!!$             CPU,info_update,&
!!$             p1p0,&
!!$             inputs_data)
!!$        call CPU%WASTED%set('stop') 
!!$        call CPU%ALGORITHM%set('stop')
!!$
!!$        if ( info_update .ne. 0 ) then
!!$           if (ctrl%id_time_discr .eq. 1) exit
!!$           
!!$           !
!!$           ! in case of failure 
!!$           ! 1- recover data
!!$           ! 2- reset controls
!!$           !
!!$           nrestart = nrestart + 1 
!!$
!!$           !
!!$           ! restores previous data
!!$           !
!!$           tdpot%tdens = tdpot%tdens_old 
!!$           tdpot%gfvar = tdpot%gfvar_old
!!$           tdpot%pot   = tdpot%pot_old
!!$
!!$           !
!!$           ! reset controls for next attemp of update
!!$           !
!!$           write(*,*) 'info_update=', info_update
!!$           call handle_failure_update ( lun_err,&
!!$                info_update,&
!!$                nrestart,&
!!$                deltat,&
!!$                p1p0,&
!!$                tdpot, &
!!$                ctrl_read,&
!!$                ctrl) 
!!$           write(*,*) 'build_prec=', ctrl%build_prec
!!$        end if
!!$     end do
!!$     !
!!$     ! reset original controls
!!$     !
!!$     ctrl = ctrl_read
!!$     
!!$     tdpot%nrestart_newton = nrestart
!!$     nupdate = nupdate + 1 + nrestart
!!$     
!!$     !
!!$     ! no time was wasted
!!$     !
!!$     if ( nrestart .eq. 0 ) then
!!$        CPU%WASTED=wasted_temp
!!$     end if
!!$
!!$     !
!!$     ! if update procedure failed break 
!!$     !
!!$     if ( info_update .ne. 0) then
!!$        rc = IOerr(lun_err, wrn_inp, 'main', &
!!$             ' UPDATE procedure failed')
!!$        exit
!!$     else
!!$        !
!!$        ! shift data from second to first slot, thus
!!$        ! first slot and tdpot are syncronized
!!$        !
!!$        !call inputs_data%shift()
!!$     end if
!!$     !
!!$     write(str, '(a,I5,a,I2,a)') &
!!$          ' UPDATE end itemp= ', itemp, &
!!$          ' nrestart=', nrestart
!!$     write(lun_stat,'(a)') ctrl%separator(etb(str))
!!$     write(lun_out,'(a)') ctrl%separator(etb(str))
!!$     
!!$     !
!!$     ! Update time varying functionals
!!$     ! Evaluate time, var_tdens and time funtionals
!!$     !
!!$     time(itemp) = time(itemp-1) + deltat  
!!$     var_tdens(itemp) = p1p0%eval_var_tdens(tdpot,deltat)
!!$     var_tdens_linfty(itemp) = p1p0%eval_var_tdens(tdpot,deltat,0.0d0)
!!$     tdpot%loc_var_tdens = var_tdens(itemp)
!!$     call p1p0%evaluate_functionals(tdpot,inputs_data)
!!$     call timefun%set_evolfun(itemp,tdpot,p1p0,inputs_data,ctrl,CPU)
!!$
!!$     !
!!$     ! Info: Time, var_tdens,  time_fucntionals (basic) 
!!$     !
!!$     write(lun_stat,'(a)') ' ' 
!!$     write(lun_out,'(a)') ' '
!!$     call info_muffe(lun_stat,&
!!$          ctrl,itemp, time(itemp),deltat, timefun, var_tdens(itemp),var_tdens_linfty(itemp),nupdate,CPU)
!!$     call info_muffe(lun_out,&
!!$          ctrl,itemp, time(itemp),deltat, timefun, var_tdens(itemp),var_tdens_linfty(itemp),nupdate,CPU)
!!$     
!!$
!!$     !
!!$     ! if continuing cyclying     
!!$     ! Test max time iterations or convergence achieved
!!$     test_var   = ( var_tdens_linfty(itemp) .gt. ctrl%tol_var_tdens )
!!$     test_maxit = ( itemp            .lt. ctrl%max_time_iterations    )
!!$     in_cycle   = ( test_var .and. test_maxit)
!!$          
!!$  end do
!!$  
!!$  call p1p0%kill(lun_err)
!!$
!!$  tdens=tdpot%tdens
!!$  pot(1:npot)=tdpot%pot
!!$  
!!$
!!$end subroutine main_subroutine


  
  


!!$program main2
!!$   use Globals
!!$  use Timing
!!$  use TimeInputs
!!$  use IOdefs
!!$  use DataSequence
!!$  use vtkloc
!!$  use ControlParameters
!!$  use MuffeTiming
!!$  use Geometry2d
!!$  use AbstractGeometry
!!$  use P1Galerkin
!!$  use TdensPotentialSystem
!!$  use TimeFunctionals
!!$  use SparseMatrix
!!$  use StdSparsePrec
!!$  use OdeInputs
!!$  use Matrix
!!$  use Preconditioner
!!$  use LinearSolver
!!$  use StdSparsePrec
!!$    
!!$  implicit none
!!$  
!!$  type(IOfdescr) :: IOfiles
!!$  type(OdeInp)   :: odein
!!$  type(OdeData)  :: inputs_data
!!$  type(CtrlPrm)  :: ctrl_read,ctrl
!!$  type(codeTim)  :: CPU
!!$  type(abs_simplex_mesh)     :: grid,subgrid
!!$  type(tdpotsys) :: tdpot
!!$  
!!$
!!$  
!!$  type(p1p0_space_discretization) :: p1p0
!!$  type(evolfun)  :: timefun
!!$  ! work arrray
!!$  type(Tim) :: wasted_temp
!!$  
!!$
!!$  !
!!$  ! Initialization of Timing module
!!$  !
!!$  call CPU%init()
!!$  call CPU%TOT%set('start')
!!$  call CPU%OVH%set('start')
!!$
!!$   
!!$  !
!!$  ! Initialization of INPUT_OUTPUT module
!!$  !
!!$  write(6,*) 'Init INPUT/OUTPUT files'
!!$  call IOfiles%init()
!!$  lun_err   = IOfiles%stderr%lun
!!$  lun_out   = IOfiles%stdout%lun
!!$  lun_stat = IOfiles%statistic%lun
!!$  write(6,*) 'All files are open'
!!$
!!$  
!!$  !
!!$  ! Initialization of  ctrl, odein, grid/subgrid from file
!!$  !
!!$  call init_p1p0_muffe_from_files(IOfiles, ctrl_read, odein, grid, subgrid )
!!$
!!$
!!$  !
!!$  ! init inputs_data containg all inputs data
!!$  !
!!$  call inputs_data%init(lun_err,  p1p0%ntdens, p1p0%npot)
!!$
!!$
!!$end program main2






!!$subroutine reverse_comunication_cycle(&
!!$     flag,ntime,info,time_iteration,current_time,&
!!$     inputs_data,ctrl, p1p0,tdpot,)
!!$  use Globals
!!$  use ControlParameters
!!$  use AbstractGeometry
!!$  use TdensPotentialSystem
!!$  use DmkOdeData
!!$  implicit none
!!$  integer,                         intent(inout) :: flag
!!$  integer,                         intent(inout) :: ntime
!!$  integer,                         intent(inout) :: info
!!$  integer,                         intent(inout) :: time_iter
!!$  real(kind=double),               intent(inout) :: current_time
!!$  type(OdeData),                   intent(in   ) :: inputs_data
!!$  type(CtrlPrm),                   intent(inout) :: ctrl
!!$  type(p1p0_space_discretization), intent(in   )  :: p1p0
!!$  type(tdpotsys),                  intent(inout) :: tdpot
!!$  
!!$  
!!$  
!!$  if ( flag==1) then
!!$     !
!!$     !  assign inputs
!!$     !
!!$     flag=2
!!$     ntime=1
!!$     current_time=cltr%t_zero
!!$     return
!!$  end if
!!$
!!$  !
!!$  ! inputs_data assigned
!!$  !
!!$  if ( flag==2) then
!!$     if ( iter .eq. 0 ) then 
!!$        call tdpot%syncronize_tzero(p1p0,inputs_data)
!!$        info=0
!!$        flag=4
!!$     else
!!$        if ( info_update == 0) then
!!$           call set_controls_before_update()
!!$        else
!!$           call ctrl%set_controls_after_failure()
!!$        end if
!!$        
!!$        call tdpot%update(&
!!$             lun_err,lun_out,lun_stat,&
!!$             ctrl,&
!!$             deltat,&
!!$             itemp,&
!!$             current_time,&
!!$             CPU,info_update,inputs_data)   
!!$        if ( info_update .ne. 0 ) then
!!$           !
!!$           nrestart = nrestart + 1 
!!$           !
!!$           ! restores previous data
!!$           !
!!$           tdpot%tdens  = tdpot%tdens_old 
!!$           tdpot%pot    = tdpot%pot_old
!!$
!!$           call tpot%build_grad_vars 
!!$
!!$           current_time = current_time
!!$
!!$           flag = 2
!!$           info_reverse = -1
!!$
!!$           if (nrestart .ge. ctrl%nrestart_max) then
!!$              flag=-1
!!$              info=-1
!!$           end if
!!$        else
!!$           ! syncronize inputs_data eval variation
!!$           flag=3
!!$        end if
!!$
!!$        if ( flag .eq. 3) then  
!!$           ! eval functional 
!!$           flag=4
!!$        end if
!!$
!!$        if ( flag .eq. 4) then  
!!$           ! eval functional 
!!$           flag=2
!!$           ntime=2
!!$           info_update=0
!!$        end if
!!$
!!$
!!$     end if
!!$  end if
!!$end subroutine reverse_comunication_cycle
!!$  
!!$
!!$
!!$end subroutine reserse_comunication_cycle
!!$
!!$! init iofiles
!!$! 
!!$
!!$
!!$
!!$info = 1
!!$do while( info > 0 ) !
!!$   select case ( info )
!!$   case (1)
!!$      call inputs_data%set(&
!!$           1,  time(itemp-1)+deltat,&
!!$           IOfiles,&
!!$           tdpot%odein)
!!$      
!!$      call tdpot%syncrozined_tzero()
!!$
!!$   case (2) 
!!$      if ( info_update .eq. 0) then
!!$         itemp = itemp + 1
!!$         write(str,'(a)') ' '
!!$         write(lun_out,'(a)') ctrl%separator(etb(str))
!!$         write(lun_stat,'(a)') ctrl%separator(etb(str))
!!$         write(str, '(a,I5)') 'UPDATE begin: itemp= ', itemp
!!$         write(lun_stat,'(a)') ctrl%separator(etb(str))
!!$         write(lun_out,'(a)') ctrl%separator(etb(str))
!!$
!!$         !
!!$         ! set all controls for update system
!!$         !
!!$         call  set_controls_update (&
!!$              itemp,&
!!$              deltat,&
!!$              tdpot, &
!!$              inputs_data,&
!!$              timefun,&
!!$              ctrl_read,&
!!$              ctrl)
!!$      end if
!!$
!!$      !
!!$      ! act only in this section and data inputs at time t+deltat
!!$      ! in the second slot
!!$      call inputs_data%set(&
!!$           2,  time(itemp-1)+deltat,&
!!$           IOfiles,&
!!$           tdpot%odein)
!!$
!!$      !
!!$      ! update whole system
!!$      !
!!$      call tdpot%update(info_update)
!!$
!!$      if ( info_update .eq. 0 ) then
!!$         !
!!$         ! in case of succesfull update, set time to 
!!$         !
!!$         write(str, '(a,I5,a,I2,a)') &
!!$              ' UPDATE end itemp= ', itemp, ' nrestart=', nrestart
!!$         write(lun_stat,'(a)') ctrl%separator(etb(str))
!!$         write(lun_out,'(a)') ctrl%separator(etb(str))
!!$         
!!$
!!$         !
!!$         ! Update time varying functionals
!!$         ! Evaluate time, var_tdens and time funtionals
!!$         !
!!$         time(itemp) = time(itemp-1) + deltat  
!!$         var_tdens(itemp) = tdpot%eval_var_tdens(deltat)
!!$         tdpot%loc_var_tdens = var_tdens(itemp)
!!$
!!$         
!!$         !
!!$         ! if continuing cyclying     
!!$         ! Test max time iterations or convergence achieved
!!$         test_var   = ( var_tdens(itemp) .gt. ctrl%tol_var_tdens )
!!$         test_maxit = ( itemp            .lt. ctrl%max_time_iterations    )
!!$         in_cycle   = ( test_var .and. test_maxit)
!!$      end if
!!$   end select
!!$   !
!!$   ! this section is skipped if update failed
!!$   !
!!$   if (info_update .eq.0) then
!!$      !
!!$      ! Print iteration number, current_time and var_tdens
!!$      ! 
!!$      write(lun,'(a)') ctrl%separator('info evolution ')
!!$      if ( current_var_tdens .lt. huge ) then 
!!$         out_format=ctrl%formatting('aiaaaraaaraaar')
!!$         write(lun,out_format) &
!!$              'iter=', itemp,&
!!$              sep,'time ',' = ', current_time,&
!!$              sep,'deltat',' = ', deltat,& 
!!$              sep,'var  ',' = ', current_var_tdens
!!$      else
!!$         out_format=ctrl%formatting('aiaaar')
!!$         write(lun,out_format) &
!!$              'iter=', itemp,&
!!$              sep,'time ',' = ', current_time
!!$      end if
!!$
!!$      !
!!$      ! Print info functional at current_time
!!$      !
!!$      if ( print_info_system ) then
!!$         write(lun,'(a)') ctrl%separator('info system')
!!$         call timefun%info(lun,1,ctrl,itemp)
!!$      end if
!!$      write(lun,'(a)') ctrl%separator(' info algoritm')
!!$      call CPU%ALGORITHM%info(lun,'Algorithm Time : ')
!!$      write(lun,*) 'time step =',itemp,' nupdate =', nupdate
!!$      write(lun,'(a)') ctrl%separator()
!!$
!!$      !
!!$      ! Save data 
!!$      !
!!$      call CPU%SAVE%set('start')
!!$      call evol2dat(tdpot,&
!!$           ctrl,&
!!$           in_cycle,&           
!!$           itemp, &
!!$           int_before_dat, int_before_matrix, &
!!$           time(itemp),&
!!$           var_tdens(itemp),&
!!$           IOfiles)
!!$      call CPU%SAVE%set('stop')
!!$   end if
!!$
!!$   
!!$end do
!!$   
!!$   
!!$
!!$   
!!$
!!$
!!$
!!$


  !
  ! flag = 1 : complete system at 
  !
!!$  
!!$  flag = 1
!!$  do while ( flag .ne. 0 ) 
!!$     
!!$
!!$     call reverse_comunication_cycle(&
!!$          flag,current_time,&
!!$          tdpot,odedata,ctrl, p1p0)
!!$     !
!!$     ! mandatory action to pass inputs tothe algorithm
!!$     ! In case of ode's inputs that are all fix in time,
!!$     ! initialiazed the values before reverse communication cycle.
!!$     ! 
!!$     !
!!$     select case ( flag )
!!$     case (2) 
!!$        !
!!$        ! set ode inputs in odedata
!!$        ! 
!!$        call odein%set(IOfiles,current_time+deltat,endfile)
!!$        call inputs_data%set(&
!!$             2,  current_time+deltat,&
!!$             IOfiles,&
!!$             odein,ctrl%id_ode)
!!$     case (3)
!!$        !
!!$        ! Syncronized data tdpot data with ode's inputs
!!$        ! in the first slot
!!$        !
!!$        call inputs_data%shift()
!!$
!!$
!!$     case (4)
!!$        !
!!$        ! print information on state of the system 
!!$        ! or save data
!!$        !
!!$         
!!$     
!!$        !
!!$        ! store data 
!!$        !
!!$        call CPU%SAVE%set('start')
!!$        call evol2dat(tdpot,&
!!$             ctrl,&
!!$             in_cycle,&           
!!$             itemp, &
!!$             int_before_dat, int_before_matrix, &
!!$             time(itemp),&
!!$             var_tdens(itemp),&
!!$             IOfiles)
!!$        call CPU%SAVE%set('stop')
!!$
!!$
!!$     case(5)
!!$        !
!!$        ! Evalute stop criteria
!!$        !
!!$        
!!$        !
!!$        ! Evaluate time, var_tdens and time funtionals
!!$        !
!!$        time(itemp) = time(itemp-1) + deltat  
!!$        var_tdens(itemp) = tdpot%eval_var_tdens(p1p0,deltat)
!!$        
!!$        !
!!$        ! if continuing cyclying     
!!$        ! Test max time iterations or convergence achieved
!!$        test_var   = ( var_tdens(itemp) .gt. ctrl%tol_var_tdens )
!!$        test_maxit = ( itemp            .lt. ctrl%max_time_iterations    )
!!$        if ( test_var .and. test_maxit) then
!!$           flag=0
!!$        else
!!$           flag=2
!!$        end if
!!$     end select
!!$  end  do

subroutine save_assembler(p1p0,IOfiles)
  use Globals
  use TimeInputs
  use DmkP1P0
  use Iodefs
  use SparseMatrix
  type(p1p0_space_discretization), intent(in   ) :: p1p0
  type(IOfdescr), intent(in   ) :: IOfiles
  !local
  integer :: lun,i,j
  character(len=1024) :: fname
  type(spmat) :: gradx,grady


  fname=etb(etb(IOfiles%output_folder%fn)//'linsys'//'/stiff.dat')
  lun=IOfiles%linear_sys%lun
  open(lun,file=fname)
  call p1p0%stiff%write(lun,'matlab')
  close(lun)

  fname=etb(etb(IOfiles%output_folder%fn)//'linsys'//'/assembler_stiff.dat')   
  open(lun,file=fname)
  do i=1,p1p0%grid_pot%ncell
     do iloc=1,p1p0%grid_pot%nnodeincell
        do jloc=1,p1p0%grid_pot%nnodeincell
           write(lun,*) p1p0%p1%assembler_csr(jloc,iloc,i)
        end do
     end do
  end do
     close(lun)

     fname=etb(etb(IOfiles%output_folder%fn)//'linsys'//'/matrix_B.dat')
     open(lun,file=fname)
     call p1p0%B_matrix%write(lun,'matlab')
     close(lun)
  
     fname=etb(etb(IOfiles%output_folder%fn)//'linsys'//'/assembler_B.dat')
      open(lun,file=fname)
     do i=1,p1p0%grid_pot%ncell
        write(lun,*) p1p0%assembler_Bmatrix_subgrid(:,i)
     end do
     close(lun)

     fname=etb(etb(IOfiles%output_folder%fn)//'linsys'//'/sizesubcell.dat')
     open(lun,file=fname)
     call write_steady(6,lun,p1p0%grid_pot%ncell,p1p0%grid_pot%size_cell)
     close(lun)

     fname=etb(etb(IOfiles%output_folder%fn)//'linsys'//'/sizecell.dat')
     open(lun,file=fname)
     call write_steady(6,lun,p1p0%grid_tdens%ncell,p1p0%grid_tdens%size_cell)
     close(lun)

     call p1p0%p1%gradbase2spmat(lun_err,gradx,grady)
     fname=etb(etb(IOfiles%output_folder%fn)//'linsys'//'/gradx.dat')
     open(lun,file=fname)
     call gradx%write(lun,'matlab')
     close(lun)

     fname=etb(etb(IOfiles%output_folder%fn)//'linsys'//'/grady.dat')
     open(lun,file=fname)
     call grady%write(lun,'matlab')
     close(lun)

     
     call gradx%kill(lun_err)
     call grady%kill(lun_err)

   end subroutine save_assembler
   
