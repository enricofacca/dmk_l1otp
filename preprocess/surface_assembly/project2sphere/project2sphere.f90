program project2sphere
  implicit none
  integer nnode,ntria
  integer i,j,k,itria,inode,ninput
  integer, allocatable :: triang(:,:),zone(:)
  real*8, allocatable  :: coord(:,:)
  real*8 :: radius,factor

  character(len=256) :: arg
  character(len=256) :: gridin
  character(len=256) :: gridout

  
  CALL getarg(1, arg)
  read(arg,'(a256)') gridin
  CALL getarg(2, arg)
  read(arg,*) radius
  CALL getarg(3, arg)
  read(arg,'(a256)') gridout


  open (10,file=etb(gridin))
  read (10,*) nnode
  read (10,*) ntria
  allocate (triang(4,ntria),coord(3,nnode))
  do inode=1,nnode
     read (10,*) coord(1,inode),coord(2,inode),coord(3,inode)
  end do
  do itria=1,ntria
     read (10,*) (triang(j,itria), j=1,4)
  end do
  close (10)


  open (10,file=etb(gridout))
  write (10,*) nnode
  write (10,*) ntria
  do inode=1,nnode
     if ( abs(coord(3,inode))> 0.8) then 
        factor=radius/sqrt(coord(1,inode)**2+coord(2,inode)**2+coord(3,inode)**2)
        write (10,*) factor*coord(1,inode),factor*coord(2,inode),factor*coord(3,inode)
     else
        factor=sqrt((radius**2-coord(3,inode)**2)/(coord(1,inode)**2+coord(2,inode)**2))
        write (10,*) factor*coord(1,inode),factor*coord(2,inode),coord(3,inode)
     end if
  end do
  do itria=1,ntria
     write (10,*) (triang(j,itria), j=1,4)
  end do  
  close(10)

  deallocate (triang,coord)
  

contains
  function etb(strIn) result(strOut)
    implicit none
    ! vars
    character(len=*), intent(in) :: strIn
    character(len=len_trim(adjustl(strIn))) :: strOut

    strOut=trim(adjustl(strIn))
  end function etb



end program project2sphere
