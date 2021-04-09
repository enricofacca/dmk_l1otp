!>---------------------------------------------------------------------
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

PROGRAM subgrid_var
  use Globals
  use IOdefs
  use GeometrySurface
  use SparseMatrix
  use TimeInputs
  use TimeOutput
  use P1GalerkinSurface
  
  implicit none
  
  type(mesh)     :: grid,subgrid
  type(IOfdescr) :: IOfiles
  type(TimeData) :: source,sink,dirac_source,dirac_sink
  type(TDOut)    :: rhs_subgrid,rhs_grid

  integer  :: nref
  logical  :: rc
  integer  :: res
  character(len=256) :: file_fnames

  integer :: stderr,stdout,debug !,res
  integer :: i, icell
  integer :: ngrids

  real(kind=double), allocatable  :: boundary_flux(:,:)
  real(kind=double), allocatable  :: tdens(:)
  real(kind=double) :: time,imbalance

  logical :: steady
  logical :: end_reached
  logical :: end_files
  logical :: end_reached_source
  logical :: end_reached_sink  
  logical :: end_reached_dirac_source 
  logical :: end_reached_dirac_sink   

  type(file) :: stiff
  type(spmat) :: stiffmat
  type(p1gal) :: p1

  stderr=6
  stdout=6


  !>----------------------------------------------------------------------------
  
  !>----------------------------------------------------------------------------
  !> Geometry info
  !> Read original grid
  ! no renumbering
  ! no connection of second level
  call getarg(1,file_fnames)

  call IOfiles%init()

  ! init grid
  call grid%init(stderr,0,Iofiles%grid)
  ! init subgrid
  call subgrid%init(stderr,1,input_mesh=grid,flag_reorder=1)
  ! write geom. info
  call subgrid%write(stderr,IOfiles%subgrid) 
  call subgrid%write_parent(stderr, IOfiles%parent)


!!$  call p1%init(6,subgrid)
!!$  allocate(tdens(subgrid%ncell))
!!$  tdens=one
!!$  call p1%build_stiff(6,'csr',tdens,stiffmat)
!!$  call p1%kill(6)
!!$  stiff%lun=111
!!$  stiff%fn='stiff.eps'
!!$  open(stiff%lun,file=stiff%fn)
!!$  call stiffmat%plot2ps(1,stiff%lun)
!!$  !call stiffmat%write(stiff%lun)
!!$  close(stiff%lun)
  

  

  ! piecewise source and sink file
  call source%init(stderr, IOfiles%source, 1, grid%ncell)
  call sink%init(stderr, IOfiles%sink, 1, grid%ncell)

  ! dirac source and sink
  call dirac_source%init(stderr, IOfiles%dirac_source, 1, grid%nnode)
  call dirac_sink%init(stderr, IOfiles%dirac_sink, 1, grid%nnode)

  ! boundary flux data
  !call boundary_flux%init(stderr, fforcing, 2, grid%nedge_bc)
  allocate(boundary_flux(2,grid%nedge_bc),stat=res)
  if(res .ne. 0) rc = IOerr(stderr, err_alloc, 'main',& 
       'work arrays boundary flux',res)
  boundary_flux=zero
  
  ! init output rhs_subgrid
  call rhs_subgrid%init(stderr, 1,subgrid%nnode)
  call rhs_grid%init(stderr, 1,grid%nnode)



  ! build
  end_files=.false.
  time=source%TDtime(1)
  steady=.false.
  write(IOfiles%rhs_subgrid_integrated%lun,*) 1,subgrid%nnode,' ! dim data' 
  write(IOfiles%rhs_grid_integrated%lun,*) 1,grid%nnode,' ! dim data' 


  do while ( (time .eq. source%TDtime(1) ) .or. &
       ( (.not. steady ) .and. ( .not. end_files  )  ) ) 
     ! read inputs
     call source%set(stderr, IOfiles%source,time,end_reached_source)
     call sink%set(stderr, IOfiles%sink,time,end_reached_sink)
     call dirac_source%set(stderr, &
          IOfiles%dirac_source,time,end_reached_dirac_source)
     call dirac_sink%set(stderr, &
          IOfiles%dirac_sink,time,end_reached_dirac_sink)
     
     end_files = ( &
          end_reached_source       .and. &
          end_reached_sink         .and. &
          end_reached_dirac_source .and. &
          end_reached_dirac_sink   )

     steady = ( &
          source%steadyTD       .and. &
          sink%steadyTD         .and. &
          dirac_source%steadyTD .and. &
          dirac_sink%steadyTD   )


     ! assembly rhs subgrid 
     rhs_subgrid%time=time
     call assembly_rhs_subgrid_integrated(grid,subgrid,&
          source%TDactual,&
          sink%TDactual,&
          dirac_source%TDactual,&
          dirac_sink%TDactual,&
          boundary_flux,&
          rhs_subgrid%TDactual)
     rhs_subgrid%steadyTD=steady
     
     ! check imbalance
     imbalance = sum(rhs_subgrid%TDactual)
     write(*,*) 'inbalance =', imbalance
     if (abs(imbalance) .gt. 1.0d-12) then
        write(*,*) 'ibalance =', imbalance
        rc = IOerr(stderr, err_inp, 'main', 'inbalanced inputs')
     end if

     !> write rhs_subgrid_integrated
     call rhs_subgrid%write2dat(IOfiles%rhs_subgrid_integrated%lun)


     !
     ! assembly rhs grid
     !
     rhs_subgrid%time=time
     call assembly_rhs_grid_integrated(grid,&
          source%TDactual,&
          sink%TDactual,&
          dirac_source%TDactual,&
          dirac_sink%TDactual,&
          boundary_flux,&
          rhs_grid%TDactual)
     rhs_grid%steadyTD=steady
     
     ! check imbalance
     imbalance = sum(rhs_grid%TDactual)
     if (imbalance .gt. 1.0d-12) then
        rc = IOerr(stderr, err_val, 'main', 'inbalanced inputs')
     end if
     

     ! write rhs_subgrid_integrated
     call rhs_grid%write2dat(IOfiles%rhs_grid_integrated%lun)



     ! next time
     time=source%TDtime(2)
  end do

  call rhs_subgrid%kill(stderr)
    
  call source%kill(stderr)
  call sink%kill(stderr)
  call dirac_source%kill(stderr)
  call dirac_sink%kill(stderr)

  call subgrid%kill(stderr)
  call grid%kill(stderr)

end PROGRAM subgrid_var

!>---------------------------------------------------------
!> Assembly rhs_subgrid of elliptic equation given forcing and Neumann terms
!>----------------------------------------------------------
subroutine assembly_rhs_subgrid_integrated(grid,subgrid,&
     source,sink,&
     dirac_source,dirac_sink,&
     boundary_flux,rhs_subgrid_forcing)
    use Globals
    use GeometrySurface
    implicit none
    type(mesh),        intent(in ) :: grid, subgrid
    real(kind=double), intent(in ) :: source(grid%ncell)
    real(kind=double), intent(in ) :: sink(grid%ncell)
    real(kind=double), intent(in ) :: dirac_source(grid%nnode)
    real(kind=double), intent(in ) :: dirac_sink(grid%nnode)
    real(kind=double), intent(in ) :: boundary_flux(2,grid%nedge_bc)
    real(kind=double), intent(out) :: rhs_subgrid_forcing(subgrid%nnode)

    !local 
    integer :: inode, icell, iedge, iloc,ifather,inode_parent
    integer :: n_sub(3)
    real(kind=double) :: neum_contribution,ddot
    real(kind=double) :: plus_mass,minus_mass

    rhs_subgrid_forcing  = zero
    do icell = 1, subgrid%ncell
       ifather=subgrid%cell_parent(icell)
       do iloc = 1,3
          inode = subgrid%topol(iloc,icell)
          rhs_subgrid_forcing(inode) =  rhs_subgrid_forcing(inode) + &
               onethird * ( source(ifather) - sink(ifather) ) * &
               subgrid%size_cell(icell)
       end do
    end do

    plus_mass=zero
    minus_mass=zero
    do inode = 1, subgrid%nnode
       if ( subgrid%node_parent(1,inode) .eq. &
            subgrid%node_parent(2,inode) ) then
          inode_parent = subgrid%node_parent(1,inode)
          if ( abs(dirac_source(inode_parent) ) .gt. small ) then
             rhs_subgrid_forcing(inode) =  dirac_source(inode_parent)
          end if
          if ( abs(dirac_sink(inode_parent)) .gt. small ) then
             rhs_subgrid_forcing(inode) =  -dirac_sink(inode_parent)
          end if
       end if
       if ( rhs_subgrid_forcing(inode) > zero ) then
          plus_mass=plus_mass+rhs_subgrid_forcing(inode)
       end if
       if ( rhs_subgrid_forcing(inode) < zero ) then
          minus_mass=minus_mass+rhs_subgrid_forcing(inode)
       end if
    end do
    
    if ( abs(sum(rhs_subgrid_forcing) )> 1.0d0-12) then
       do inode = 1, subgrid%nnode
          if ( rhs_subgrid_forcing(inode) < zero ) then
             rhs_subgrid_forcing(inode) = rhs_subgrid_forcing(inode) * &
                  abs(plus_mass/minus_mass) 
          end if
       end do
    end if

       


    ! add neumman contribution
!!$    do iedge = 1,grid%nedge_bc
!!$       neum_contribution = ddot(2,&
!!$            boundary_flux(:,iedge),1,&
!!$            grid%normal(:,iedge),1) * &
!!$            grid%leng_edge(iedge) * onehalf
!!$
!!$       n_sub(1) = subgrid%node_sons( grid%iside(1,iedge) )
!!$       n_sub(2) = subgrid%node_sons( grid%iside(2,iedge) )
!!$       n_sub(3) = subgrid%node_sons( subgrid%nnode_parent + iedge )
!!$
!!$       ! node on grid
!!$       do iloc=1,2
!!$          inode = n_sub(iloc)
!!$          rhs_subgrid_forcing(inode) = rhs_subgrid_forcing(inode) - neum_contribution / 2.0d0 
!!$       end do
!!$       inode = n_sub(3)
!!$       rhs_subgrid_forcing(inode) = rhs_subgrid_forcing(inode) - neum_contribution
!!$    end do

  end subroutine assembly_rhs_subgrid_integrated

!>---------------------------------------------------------
!> Assembly rhs_subgrid of elliptic equation given forcing and Neumann terms
!>----------------------------------------------------------
subroutine assembly_rhs_grid_integrated(grid,&
     source,sink,&
     dirac_source,dirac_sink,&
     boundary_flux,rhs_grid_forcing)
    use Globals
    use GeometrySurface
    implicit none
    type(mesh),        intent(in ) :: grid
    real(kind=double), intent(in ) :: source(grid%ncell)
    real(kind=double), intent(in ) :: sink(grid%ncell)
    real(kind=double), intent(in ) :: dirac_source(grid%nnode)
    real(kind=double), intent(in ) :: dirac_sink(grid%nnode)
    real(kind=double), intent(in ) :: boundary_flux(2,grid%nedge_bc)
    real(kind=double), intent(out) :: rhs_grid_forcing(grid%nnode)

    !local 
    integer :: inode, icell, iedge, iloc,ifather,inode_parent
    integer :: n_sub(3)
    real(kind=double) :: neum_contribution,ddot
    real(kind=double) :: plus_mass, minus_mass

    rhs_grid_forcing  = zero
    do icell = 1, grid%ncell
       do iloc = 1,3
          inode = grid%topol(iloc,icell)
          rhs_grid_forcing(inode) =  rhs_grid_forcing(inode) + &
               onethird * ( source(icell) - sink(icell) ) * &
               grid%size_cell(icell)
       end do
    end do

    do inode = 1, grid%nnode
       if ( abs(dirac_source(inode) ) .gt. small ) then
          rhs_grid_forcing(inode) =  dirac_source(inode)
       end if
       if ( abs(dirac_sink(inode)) .gt. small ) then
          rhs_grid_forcing(inode) =  -dirac_sink(inode)
       end if
       
       if ( rhs_grid_forcing(inode) > zero ) then
          plus_mass=plus_mass+rhs_grid_forcing(inode)
       end if
       if ( rhs_grid_forcing(inode) < zero ) then
          minus_mass=minus_mass+rhs_grid_forcing(inode)
       end if
    end do
    
    if ( abs(sum(rhs_grid_forcing) )> 1.0d0-12) then
       do inode = 1, grid%nnode
          if ( rhs_grid_forcing(inode) < zero ) then
             rhs_grid_forcing(inode) = rhs_grid_forcing(inode) * &
                  abs(plus_mass/minus_mass) 
          end if
       end do
    end if


    ! add neumman contribution
!!$    do iedge = 1,grid%nedge_bc
!!$       neum_contribution = ddot(2,&
!!$            boundary_flux(:,iedge),1,&
!!$            grid%normal(:,iedge),1) * &
!!$            grid%leng_edge(iedge) * onehalf
!!$
!!$       n_sub(1) = subgrid%node_sons( grid%iside(1,iedge) )
!!$       n_sub(2) = subgrid%node_sons( grid%iside(2,iedge) )
!!$       n_sub(3) = subgrid%node_sons( subgrid%nnode_parent + iedge )
!!$
!!$       ! node on grid
!!$       do iloc=1,2
!!$          inode = n_sub(iloc)
!!$          rhs_grid_forcing(inode) = rhs_grid_forcing(inode) - neum_contribution / 2.0d0 
!!$       end do
!!$       inode = n_sub(3)
!!$       rhs_grid_forcing(inode) = rhs_grid_forcing(inode) - neum_contribution
!!$    end do

  end subroutine assembly_rhs_grid_integrated
