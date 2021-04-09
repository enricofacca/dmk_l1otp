!>-------------------------------------------------------------
!> Define global vars
!>
!> err numbers:
!> errno = [0:1][01-10] I/O errors (file existence etc...)
!>       = [0:1][11-20] input errors (reads)
!>       = [0:1][21-30] output errors (writes)
!>       = [0:1][31-40] alloc errors 
!>       = [0:1][41-50] dealloc errors
!>       = [0:1][51-60] vtk libray errors
!<------------------------------------------------------------- 
module Globals
  implicit none
  !> single precision parameter for real vars
  integer, parameter ::  single = kind(1.0e0)
  !> double precision parameter for real vars
  integer, parameter ::  double = kind(1.0d0)
  !> error codes 
  integer, parameter :: wrn_IO=1,wrn_read=11,wrn_write=21
  integer, parameter :: wrn_val=3,wrn_inp=13,wrn_out=23
  integer, parameter :: err_IO=101,err_inp=111,err_out=121
  integer, parameter :: err_read=161,err_write=171
  integer, parameter :: err_alloc=131,err_dealloc=141
  integer, parameter :: err_vtk=151
  integer, parameter :: err_val=201
  !> double parameters for useful constants
  real(kind=double), parameter  :: zero=0.0d0
  real(kind=double), parameter  :: one=1.0d0
  real(kind=double), parameter  :: onehalf=0.5d0
  real(kind=double), parameter  :: onethird=1.0d0/3.0d0
  real(kind=double), parameter  :: onefourth=1.0d0/4.0d0
  real(kind=double), parameter  :: onesixth=1.0d0/6.0d0
  real(kind=double), parameter  :: verysmall=1.0d-40
  real(kind=double), parameter  :: small=1.0d-15  
  real(kind=double), parameter  :: large=1.0d10
  real(kind=double), parameter  :: huge=1.0d30
  real(kind=double), parameter  :: pigreco=4.0d0*atan(one)

  type, public :: file
     logical :: exist             !< Logical of existence
     integer :: lun               !< I/O unit number
     character (len=256) :: fn    !< I/O file name or directory name
     character (len=256) :: fnloc !< I/O file name or directory name, local
   !  integer(hid_t)  :: hdf5_id=0 !< ID numbe for corresponding dataset in h5df
   contains
     !> output the content of the variable
     !> (public for type file)
     procedure, public, pass :: info => fn_print
  end type file
   

!!$     !> turn string into lowercase
!!$     !> (private for type IOfdescr)
!!$     procedure, private, nopass :: to_lower
!!$     !> turn string into uppercase
!!$     !> (private for type IOfdescr)
!!$     procedure, private, nopass :: to_upper
!!$     !> return string with no preceding and trailing spaces
!!$     !> (private for type IOfdescr)
!!$     procedure, public, nopass :: etb
!!$     !> erase comments from string (everything
!!$     !> after character "!" or "%") 
!!$     !> (private for type IOfdescr)
!!$     procedure, private, nopass :: erase_comment
contains

  !>-------------------------------------------------------------
  !> Handle and write alert I/O errors.
  !> (public global procedure)
  !>
  !> usage:
  !>    rc= IOerr(lun, errno, call_proc [, add_msg] [, add_int])
  !> where:
  !> \param[in] lun -> integer. I/O unit for output of error messages
  !> \param[in] errno -> integer. error number
  !> \param[in] call_proc -> character. name of the procedure where
  !>                          the error occured
  !> \param[in] (optional) add_msg -> character.
  !>                                   additional message to be printed
  !> \param[in] (optional) add_int -> integer.
  !>                                   integer to be output as string
  !>                                   in the message
  !>
  !> \return rc -> logical. true if execution continues
  !>                (obviously no return otherwise as execution stops)
  !>
  !> errno convention:
  !> 0<errno<=100 -> warning: execution continues
  !> errno > 100  -> severe: execution terminates 
  !>
  !>
  !> contains:
  !> private function errmsg storing the message values
  !> 
  !<-------------------------------------------------------------
  function IOerr(lun, errno, call_proc, add_msg, add_int) result(rc)
    implicit none
    integer, intent(in) :: lun,errno
    integer, optional, intent(in) :: add_int
    logical :: rc
    character(len=*), intent(in) :: call_proc
    character(len=*), optional, intent(in) :: add_msg
    ! local vars
    character(len=256) :: msg
    integer :: int

    if (present(add_msg)) then
       msg=trim(add_msg)
       if(present(add_int)) then
          int=add_int
       else
          int = 0
       end if
    else
       msg=''
       int = 0
    end if
    write(lun,*) ''
    write(lun,*) ' *************'
    if (errno .gt. 100) then
       write(lun,fmt=100) '  SEVERE ERROR:'
       write(lun,fmt=101) ' in procedure: ',trim(call_proc)
       write(lun,fmt=101) trim(errmsg(errno, trim(msg), int))
       stop 'EXECUTION TERMINATED IN PROC IOerr'
    else if (errno .gt. 0) then
       write(lun,fmt=100) ' WARNING ERROR:'
       write(lun,fmt=101) 'in procedure: ',trim(call_proc)
       write(lun,fmt=101) trim(errmsg(errno, trim(msg), int))
       rc=.true.
    else
       write(lun,fmt=100) ' IO ERROR NUMBER not found'
       stop 'EXECUTION TERMINATED IN PROC IOerr'
    end if
    
100 format(1x,a,1x,i5)
101 format(5(1x,a))
  contains

    !>-------------------------------------------------------------
    !> return msg value indexed by a number
    !> (private procedure of IOerr function)
    !>
    !> usage:
    !>     msg = errmsg(errno, addmsg, int)
    !>
    !> where:
    !> \param[in] errno -> integer, error number
    !> \param[in] addmsg -> character, additional message to be printed
    !> \param[in] int -> integer, to be included in msg if nonzero
    !> \return msg -> character, error message
    !>
    !<-------------------------------------------------------------
    function errmsg(errno, addmsg, addint) result(msg)
      implicit none
      integer, intent(in) :: errno
      character(len=*), intent(in) :: addmsg
      integer, intent(in) :: addint
      character(len=256) :: msg
      ! local vars
      character(len=5) :: rdwr

      if(addint.eq.0) then
         rdwr=''
      else
         write(rdwr,'(i5)') addint
      end if
      
      select case (errno)
      case (wrn_IO) ! file/dir not found but continues execution
         msg='File or Directory '//trim(addmsg)//' (unit '//rdwr//') does not exist.'
      case (wrn_read) ! error in reading file
         msg='Error Reading File'//trim(addmsg)//' (unit '//rdwr//')'
      case (wrn_write) ! error in writing file
         msg='Error Writing File'//trim(addmsg)//' (unit '//rdwr//')'
      case (wrn_val) ! error in inputs data
         msg='Error In Input/Output Parameter'//trim(addmsg)
      case (wrn_inp) ! error in inputs data
         msg='Error In Input Parameter'//trim(addmsg)
      case (wrn_out) ! error in output
         msg='Error In Input Parameter'//trim(addmsg)
      case (err_IO) ! file/dir not found but stops execution
         msg='File or Directory '//trim(addmsg)//' (unit '//rdwr//') does not exist.'
      case (err_read) ! error reading a number in input
         msg = '    ...error reading file '//trim(addmsg)//' (iostat='//rdwr//')'
      case (err_write) ! error writing a number in input
         msg = '    ...error writing file '//trim(addmsg)//' (iostat='//rdwr//')'
      case (err_inp) ! error reading a number in input
         msg = 'Error In Input Parameter'//trim(addmsg)
      case (err_out) ! error writing a number in input
         msg = 'Error In Output Parameter'//trim(addmsg)
      case (err_alloc) ! error allocating an array
         msg = '    ...alloc fail var '//trim(addmsg)//' (stat='//rdwr//')'
      case (err_dealloc) ! error deallocating an array
         msg = '    ...dealloc fail var '//trim(addmsg)//' (stat='//rdwr//')'
      case (err_vtk) ! error writing vtk files
         msg = '    ...VTK library: '//trim(addmsg)//rdwr
      case (err_val) ! error in parameter values
         msg = '    ...wrong parameter value: '//trim(addmsg)//rdwr
      case default
         msg = '  no known error number'
      end select
    end function errmsg
  end function IOerr

  !>-------------------------------------------------------------
  !> Transform string strIn to upper case
  !> (private procedure for type IOfdescr)
  !>
  !> usage:
  !>     str = to_upper(strIn)
  !> 
  !> where:
  !> \param[in] strIn -> input string
  !> \return strOut -> output string
  !>
  !>
  !> Adapted from http://rosettacode.org/wiki/String_case#Fortran
  !> Original author: Clive Page
  !<-------------------------------------------------------------
  function to_upper(strIn) result(strOut)
    implicit none
    ! vars
    character(len=*), intent(in) :: strIn
    character(len=len(strIn)) :: strOut
    ! local vars
    integer :: i
    
    do i = 1, len(strIn)
       select case(strIn(i:i))
       case("a":"z")
          strOut(i:i) = achar(iachar(strIn(i:i))-32)
        end select
    end do
    
  end function to_upper
  
  !>-------------------------------------------------------------
  !> Transform string strIn to lowercase case
  !> (private procedure for type IOfdescr)
  !>
  !> usage:
  !>   str = to_lower(strIn)
  !>
  !> where:
  !> \param[in] strIn -> input string
  !> \return strOut - > output string
  !>
  !> (we'll see if it is needed somewhere else.
  !> In this case it will become public, but possibly must
  !> be moved into the Globals module (MP)
  !>
  !> Adapted from http://rosettacode.org/wiki/String_case#Fortran
  !> Original author: Clive Page
  !<-------------------------------------------------------------
  function to_lower(strIn) result(strOut)
    implicit none
    ! vars
    character(len=*), intent(in) :: strIn
    character(len=len(strIn)) :: strOut
    ! local vars
    integer :: i

    do i = 1, len(strIn)
       select case(strIn(i:i))
       case("A":"Z")
          strOut(i:i) = achar(iachar(strIn(i:i))+32)
       case default
          strOut(i:i) = strIn(i:i)
       end select
    end do
  end function to_lower

  !>-------------------------------------------------------------
  !> return string strIn with no preceding and trailing spaces
  !>     (just a short name for trim(adjustl......))
  !> (private procedure fortype IOfdescr)
  !> 
  !> usage:
  !>   str = etb(strIn)
  !>
  !> where:
  !> \param[in] strIn -> input string
  !> \return strOut -> output string
  !>
  !> (we'll see if it is needed somewhere else.
  !> In this case it will become public, but possibly must
  !> be moved into the Globals module (MP)
  !> 
  !<-------------------------------------------------------------
  function etb(strIn) result(strOut)
    implicit none
    ! vars
    character(len=*), intent(in) :: strIn
    character(len=len_trim(adjustl(strIn))) :: strOut

    strOut=trim(adjustl(strIn))
  end function etb
  
  !>-------------------------------------------------------------
  !> erase comments from a string
  !> (private procedure for type IOfdescr)
  !>
  !> usage:
  !>   str = erase_comment(strIn)
  !>
  !> where:
  !> \param[in] strIn -> input string
  !> \return strOut -> output string
  !>
  !> a comment is a trailing string beginning
  !> with character ! (F90-style) or 
  !> with character % (TeX-style) 
  !> with character # (script-bash-style) 
  !<-------------------------------------------------------------
  function erase_comment(strIn) result(strOut)
    implicit none
    ! vars
    character(len=*), intent(in) :: strIn
    character(len=len_trim(adjustl(strIn))) :: strOut
    ! local vars
    character(len=1) :: strloc
    integer :: i,j

    do i=1,len_trim(adjustl(strIn))

       strloc=strIn(i:i)
       if (strloc .eq. '!' .or. strloc.eq.'%' .or. strloc.eq.'#') then
          do j=i,len_trim(adjustl(strIn))
             strOut(j:j)=' '
          end do
          exit
       else
          strOut(i:i)=strloc
       end if
    end do
  end function erase_comment

  !>-------------------------------------------------------------
  !> Function that calculates cross-products
  !> (public global procedure)
  !> 
  !> usage:
  !>     res_cross = cross(vec1,vec2)
  !>
  !> where:
  !> \param[in] vec1, vec2 -> real(double). dimension(3)
  !>                          vectors 
  !>
  !> \return res_cross -> real(double). dimension(3)
  !>                     vector
  !>
  !<-------------------------------------------------------------
  function  cross(vecA, vecB) result (res_cross)   
    implicit none
    ! vars
    real(kind=double), intent(in) :: vecA(3), vecB(3)
    real(kind=double) :: res_cross(3)
    ! local vars

    res_cross(1) = vecA(2)*vecB(3)-vecA(3)*vecB(2)
    res_cross(2) = vecA(3)*vecB(1)-vecA(1)*vecB(3)
    res_cross(3) = vecA(1)*vecB(2)-vecA(2)*vecB(1)

  end function cross



  !>--------------------------------------------------------------
  !> Simple sort algorithm to sort in increasing order integer 
  !> array. Use only for small narray.
  !> Use global_heapsort for big arrays
  !>-------------------------------------------------------------
  subroutine isort(narray,array)
    implicit none
    integer, intent(in   ) :: narray
    integer, intent(inout) :: array(narray)
    ! 
    integer :: itemp
    integer :: i,j ,indx,isgn

    if (narray.lt.2) return
    ! Initialize.
    i = 0
    indx = 0
    isgn = 0
    j = 0
    do 
       call global_heapsort(narray, indx, i,j,isgn)
       if (indx .gt. 0 ) then
          ! SWAP ELEMENT 
          itemp    =  array(i)
          array(i) = array(j)
          array(j) = itemp
       else if ( indx .lt. 0) then
          ! COMPARE (array(i) and array(j) )
          isgn = 0
          if ( array(i) .lt. array(j)  ) isgn = -1
          if ( array(i) .gt. array(j)  ) isgn = 1
       else if ( indx .eq. 0 ) then
          exit
       end if
    end do
  end subroutine isort


  ! says if the first array a is in lexicographic order
  ! with respect to the second array b
  function lexicographic_order(n,a,b) result(a_before_b)
    implicit none
    integer, intent(in) :: n
    integer, intent(in) :: a(n),b(n)
    logical :: a_before_b
    !local 
    integer             :: i
    a_before_b = .true.
    do i=1,n
       if  ( a(i) .gt. b(i) ) then
          a_before_b = .false.    
          exit
       else if ( a(i) .lt. b(i) ) then
          a_before_b = .true.    
          exit
       end if
    end do
  end function lexicographic_order

  !>---------------------------------------------------------------
  !> Subroutine for ordering items (integer real etc) with the 
  !> heapsort algorithm via reverse comunication style
  !>--------------------------------------------------------------
  subroutine global_heapsort( n, indx, i, j, isgn )
    !
    ! SORT_HEAP_EXTERNAL externally sorts a list 
    !  of items into ascending order.
    !
    !  Discussion:
    !
    !    The actual list of data is not passed to the routine.  
    !    Hence this routine may be used to sort integers, reals,
    !    numbers, names, dates, shoe sizes, and so on.  
    !    After each call, the routine asks
    !    the user to compare or interchange two items, until a special
    !    return value signals that the sorting is completed.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    05 February 2004
    !
    !  Author:
    !
    !    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
    !    FORTRAN90 version by John Burkardt.
    !
    !  Reference:
    !
    !    Albert Nijenhuis and Herbert Wilf,
    !    Combinatorial Algorithms,
    !    Academic Press, 1978, second edition,
    !    ISBN 0-12-519260-6.
    !
    !  Parameters:
    !
    !    Input, integer ( kind = 4 ) N,
    !    the number of items to be sorted.
    !
    !    Input/output, integer ( kind = 4 ) INDX, 
    !    the main communication signal.
    !
    !    The user must set INDX to 0 before the first call.
    !    Thereafter, the user should not change the value of INDX until
    !    the sorting is done.
    !
    !    On return, if INDX is
    !
    !      greater than 0,
    !      * interchange items I and J;
    !      * call again.
    !
    !      less than 0,
    !      * compare items I and J;
    !      * set ISGN = -1 if I < J, ISGN = +1 if J < I;
    !      * call again.
    !
    !      equal to 0, the sorting is done.
    !
    !   Example of use : sort array of thing of legnth nnz
    !   Code
    !   
    !   if (nnz.lt.2) return
    !   Initialize.
    !   i = 0
    !   indx = 0
    !   isgn = 0
    !   j = 0
    !   do 
    !      call global_heapsort(nnz, indx, i,j,isgn)
    !      if (indx .gt. 0 ) then
    !          ! SWAP ELEMENT 
    !          array(i) and array(j)
    !       else if ( indx .lt. 0) then
    !          ! COMPARE (array(i) and array(j) )
    !          isgn = 0
    !          if ( array(i) .lt. array(j)  ) isgn = -1
    !          if ( array(i) .gt. array(j)  ) isgn = 1
    !       else if ( indx .eq. 0 ) then
    !          exit
    !       end if
    !    end do
    !
    !  
    !
    !    Output, integer ( kind = 4 ) I, J, the indices of two items.
    !    On return with INDX positive, elements I and J should 
    !    be interchanged.
    !    On return with INDX negative, elements I and J should 
    !    be compared, and
    !    the result reported in ISGN on the next call.
    !
    !    Input, integer ( kind = 4 ) ISGN, results of comparison of elements 
    !    I and J.  (Used only when the previous call returned INDX less than 0).
    !    ISGN <= 0 means I is less than or equal to J;
    !    0 <= ISGN means I is greater than or equal to J.
    !
    implicit none

    integer ( kind = 4 ) i
    integer ( kind = 4 ), save :: i_save = 0
    integer ( kind = 4 ) indx
    integer ( kind = 4 ) isgn
    integer ( kind = 4 ) j
    integer ( kind = 4 ), save :: j_save = 0
    integer ( kind = 4 ), save :: k = 0
    integer ( kind = 4 ), save :: k1 = 0
    integer ( kind = 4 ) n
    integer ( kind = 4 ), save :: n1 = 0
    !
    !  INDX = 0: This is the first call.
    !
    if ( indx == 0 ) then

       i_save = 0
       j_save = 0
       k = n / 2
       k1 = k
       n1 = n
       !
       !  INDX < 0: The user is returning the results of a comparison.
       !
    else if ( indx < 0 ) then

       if ( indx == -2 ) then

          if ( isgn < 0 ) then
             i_save = i_save + 1
          end if

          j_save = k1
          k1 = i_save
          indx = -1
          i = i_save
          j = j_save
          return

       end if

       if ( 0 < isgn ) then
          indx = 2
          i = i_save
          j = j_save
          return
       end if

       if ( k <= 1 ) then

          if ( n1 == 1 ) then
             i_save = 0
             j_save = 0
             indx = 0
          else
             i_save = n1
             n1 = n1 - 1
             j_save = 1
             indx = 1
          end if

          i = i_save
          j = j_save
          return

       end if

       k = k - 1
       k1 = k
       !
       !  0 < INDX, the user was asked to make an interchange.
       !
    else if ( indx == 1 ) then

       k1 = k

    end if

    do

       i_save = 2 * k1

       if ( i_save == n1 ) then
          j_save = k1
          k1 = i_save
          indx = -1
          i = i_save
          j = j_save
          return
       else if ( i_save <= n1 ) then
          j_save = i_save + 1
          indx = -2
          i = i_save
          j = j_save
          return
       end if

       if ( k <= 1 ) then
          exit
       end if

       k = k - 1
       k1 = k

    end do

    if ( n1 == 1 ) then
       i_save = 0
       j_save = 0
       indx = 0
       i = i_save
       j = j_save
    else
       i_save = n1
       n1 = n1 - 1
       j_save = 1
       indx = 1
       i = i_save
       j = j_save
    end if

    return
  end subroutine global_heapsort

  function ilocate(narray,array,tofind) result (i)
    implicit none
    integer, intent(in) :: narray
    integer, intent(in) :: tofind
    integer, intent(in) :: array(narray)
    integer i

    i = 0
    do while ( (.not. (array(i) .eq. tofind) )  .and. &
         (i .le. narray)) 
       i=i+1
    end do
    if ( i .eq. narray+1) i=0
  end function ilocate

  subroutine double_col_permute ( m, n, p, a )
    !****************************************************************************
    !
    !  ! R8COL_PERMUTE permutes an R8COL in place.
    !
    !  Discussion:
    !
    !    An R8COL is an M by N array of double precision values, regarded
    !    as an array of N columns of length M.
    !
    !    The same logic can be used to permute an array of objects of any
    !    arithmetic type, or an array of objects of any complexity.  The only
    !    temporary storage required is enough to store a single object.  The number
    !    of data movements made is N + the number of cycles of order 2 or more,
    !    which is never more than N + N/2.
    !
    !  Example:
    !
    !    Input:
    !
    !      M = 2
    !      N = 5
    !      P = (   2,    4,    5,    1,    3 )
    !      A = ( 1.0,  2.0,  3.0,  4.0,  5.0 )
    !          (11.0, 22.0, 33.0, 44.0, 55.0 )
    !
    !    Output:
    !
    !      A    = (  2.0,  4.0,  5.0,  1.0,  3.0 )
    !             ( 22.0, 44.0, 55.0, 11.0, 33.0 ).
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    04 December 2017 (Enrico Facca)
    !    simple modification and adaptation to my code
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Input, integer ( kind = 4 ) M, the dimension of objects.
    !
    !    Input, integer ( kind = 4 ) N, the number of objects.
    !
    !    Input, integer ( kind = 4 ) P(N), the permutation.  P(I) = J means
    !    that the I-th element of the output array should be the J-th
    !    element of the input array.  P must be a legal permutation
    !    of the integers from 1 to N, otherwise the algorithm will
    !    fail catastrophically.
    !
    !    Input/output, real ( kind = 8 ) A(M,N), the array to be permuted.
    !
    implicit none

    integer ( kind = 4 ) m
    integer ( kind = 4 ) n

    real ( kind = double ) a(m,n)
    real ( kind = double ), allocatable ::  a_temp(:)

    integer ( kind = 4 ) res
    integer ( kind = 4 ) iget
    integer ( kind = 4 ) iput
    integer ( kind = 4 ) istart
    integer ( kind = 4 ) p(n)
    !
    !  Search for the next element of the permutation that has not been used.
    !
    allocate(a_temp(m),stat=res)
    if (res .ne. 0) write(*,*) ' Error allocation double double_col_permute'

    do istart = 1, n

       if ( p(istart) < 0 ) then

          cycle

       else if ( p(istart) == istart ) then

          p(istart) = -p(istart)
          cycle

       else

          a_temp(1:m) = a(1:m,istart)
          iget = istart
          !
          !  Copy the new value into the vacated entry.
          !
          do

             iput = iget
             iget = p(iget)

             p(iput) = -p(iput)

             if ( iget < 1 .or. n < iget ) then
                write ( *, '(a)' ) ' '
                write ( *, '(a)' ) 'R8COL_PERMUTE - Fatal error!'
                write ( *, '(a)' ) '  A permutation index is out of range.'
                write ( *, '(a,i8,a,i8)' ) '  P(', iput, ') = ', iget
                stop
             end if

             if ( iget == istart ) then
                a(1:m,iput) = a_temp(1:m)
                exit
             end if

             a(1:m,iput) = a(1:m,iget)

          end do

       end if

    end do
    !
    !  Restore the signs of the entries.
    !
    p(1:n) = -p(1:n)

    return

    deallocate(a_temp,stat=res)
    if (res .ne. 0) write(*,*) ' Error deallocation double double_col_permute'
  end subroutine double_col_permute


  subroutine integer_col_permute ( m, n, p, a )

    !**************************************************************************
    !
    !  ! R8COL_PERMUTE permutes an R8COL in place.
    !
    !  Discussion:
    !
    !    An R8COL is an M by N array of double precision values, regarded
    !    as an array of N columns of length M.
    !
    !    The same logic can be used to permute an array of objects of any
    !    arithmetic type, or an array of objects of any complexity.  The only
    !    temporary storage required is enough to store a single object.  The number
    !    of data movements made is N + the number of cycles of order 2 or more,
    !    which is never more than N + N/2.
    !
    !  Example:
    !
    !    Input:
    !
    !      M = 2
    !      N = 5
    !      P = (   2,    4,    5,    1,    3 )
    !      A = ( 1.0,  2.0,  3.0,  4.0,  5.0 )
    !          (11.0, 22.0, 33.0, 44.0, 55.0 )
    !
    !    Output:
    !
    !      A    = (  2.0,  4.0,  5.0,  1.0,  3.0 )
    !             ( 22.0, 44.0, 55.0, 11.0, 33.0 ).
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    04 December 2017 (Enrico Facca)
    !    simple modification and adaptation to my code
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Input, integer ( kind = 4 ) M, the dimension of objects.
    !
    !    Input, integer ( kind = 4 ) N, the number of objects.
    !
    !    Input, integer ( kind = 4 ) P(N), the permutation.  P(I) = J means
    !    that the I-th element of the output array should be the J-th
    !    element of the input array.  P must be a legal permutation
    !    of the integers from 1 to N, otherwise the algorithm will
    !    fail catastrophically.
    !
    !    Input/output, real ( kind = 8 ) A(M,N), the array to be permuted.
    !
    implicit none

    integer ( kind = 4 ) m
    integer ( kind = 4 ) n

    integer ( kind = 4 ) :: a(m,n)
    integer ( kind = 4 ) ,allocatable ::  a_temp(:)

    integer ( kind = 4 ) res
    integer ( kind = 4 ) iget
    integer ( kind = 4 ) iput
    integer ( kind = 4 ) istart
    integer ( kind = 4 ) p(n)
    !
    !  Search for the next element of the permutation that has not been used.
    !
    allocate(a_temp(m),stat=res)
    if (res .ne. 0) write(*,*) ' Error allocation double double_col_permute'

    do istart = 1, n

       if ( p(istart) < 0 ) then

          cycle

       else if ( p(istart) == istart ) then

          p(istart) = -p(istart)
          cycle

       else

          a_temp(1:m) = a(1:m,istart)
          iget = istart
          !
          !  Copy the new value into the vacated entry.
          !
          do

             iput = iget
             iget = p(iget)

             p(iput) = -p(iput)

             if ( iget < 1 .or. n < iget ) then
                write ( *, '(a)' ) ' '
                write ( *, '(a)' ) 'integer COL_PERMUTE - Fatal error!'
                write ( *, '(a)' ) '  A permutation index is out of range.'
                write ( *, '(a,i8,a,i8)' ) '  P(', iput, ') = ', iget
                stop
             end if

             if ( iget == istart ) then
                a(1:m,iput) = a_temp(1:m)
                exit
             end if

             a(1:m,iput) = a(1:m,iget)

          end do

       end if

    end do
    !
    !  Restore the signs of the entries.
    !
    p(1:n) = -p(1:n)

    return

    deallocate(a_temp,stat=res)
    if (res .ne. 0) write(*,*) &
         ' Error deallocation double double_col_permute'
  end subroutine integer_col_permute

  !>-----------------------------------------------------
  !> Procedure orthogonalizing x w.r.t. a series of vectors
  !>-----------------------------------------------------s
  subroutine ortogonalize(dim,nvectors,vectors,x)
    implicit none
    integer,           intent(in   ) :: dim
    integer,           intent(in   ) :: nvectors
    real(kind=double), intent(in   ) :: vectors(dim,nvectors)
    real(kind=double), intent(inout) :: x(dim)

    ! local
    integer :: i
    real(kind=double) :: ddot,alpha,beta

    do i=1,nvectors
       beta  = ddot(dim,vectors(1,i),1,x,1)
       alpha = -beta/ddot(dim,vectors(1,i),1,vectors(1,i),1)
       !x=x-alpha*vectors(:,i)
       call daxpy(dim,alpha,vectors(1,i),1,x,1)
    end do

  end subroutine ortogonalize


  !>-------------------------------------------------------------
  !> Info procedure.
  !> (public procedure for type file)
  !> Prints content of a variable of type file
  !> 
  !> usage:
  !>     call 'var'\%info(lun)
  !>
  !> where:
  !> \param[in] lun: output unit
  !>
  !> uses private function etb
  !<-------------------------------------------------------------
  subroutine fn_print(this, lun)
    implicit none
    class(file), intent(in) :: this
    integer,     intent(in) :: lun

    write(lun,'(a,a,a,i3)') ' filename ',etb(this%fn), &
         ' linked to lun ',this%lun

  end subroutine fn_print

end module Globals

