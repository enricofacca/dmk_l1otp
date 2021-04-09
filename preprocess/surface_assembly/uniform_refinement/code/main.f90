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

PROGRAM uniform_refinement
  use Globals
  use Geometry
  

  implicit none
  
  type(mesh), allocatable :: grids(:)
  type(file) :: fgrid
  type(file) :: fsubgrid
  integer  :: nref
  logical  :: rc
  integer  :: res
  character(len=256) :: input,path_grid,path_subgrid

  integer :: stderr,stdout,debug !,res
  integer :: i, itria
  integer :: ngrids

  !>----------------------------------------------------------------------------
  
  !>----------------------------------------------------------------------------
  !> Geometry info
  !> Read original grid
  ! no renumbering
  ! no connection of second level
  call getarg(1,path_grid)
  call getarg(2,path_subgrid)
  call getarg(3,input)
  read(input,*) nref


  fgrid%fn=etb(path_grid)
  fgrid%lun=10
  fgrid%exist=.true.

  fsubgrid%fn=etb(path_subgrid)
  fsubgrid%lun=11
  fsubgrid%exist=.true.

  
  stderr=6
  
  allocate(grids(0:nref),stat=res)
  if(res .ne. 0) rc = IOerr(6, err_alloc, 'main', &
         ' subgrids',res)

  open(fgrid%lun,file=fgrid%fn)
  call grids(0)%init(6,0,fgrid)
  close(fgrid%lun)
    
  do i=1,nref 
     write(*,'(a,I2)') ' refinement nmb = ', i
     call grids(i)%init(6,1,input_mesh=grids(i-1),flag_reorder=0)
  end do  

  
  open(fsubgrid%lun,file=fsubgrid%fn)
  call grids(nref)%write(6,fsubgrid) 
  close(fsubgrid%lun)

end PROGRAM uniform_refinement
