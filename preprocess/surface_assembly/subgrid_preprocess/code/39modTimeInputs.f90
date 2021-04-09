module TimeInputs
  use Globals
  implicit none
  
  private
  !> structure variable containing vars for variable in time input data
  !> (e.g., forcing functions, exponents, k(x,t)...)
  !> up to no, only forcing function input
  type, public :: TimeData
     !> true/false flag if is initialized
     logical :: built
     !> Dimension of real array
     integer :: dimdata
     !> Dimension of real array
     integer :: ndata
     !> Dimension (2,dimdata,ndata)
     !> Input values
     real(kind=double), allocatable :: TDval(:,:,:)
     !> Dimension (2)
     !> Input values
     real(kind=double), allocatable :: TDtime(:)
     !> Dimension (dimdata,Ndata)
     !> Actual time values     
     real(kind=double), allocatable :: TDactual(:,:)
     !> Dimension (dimdata)
     !> Scratch array for reading values     
     real(kind=double), allocatable :: val(:)
     !> Time of evaluated TDactual
     real(kind=double) :: time
     !> SteadyTD=.true.  : data do not change during time
     !> SteadyTD=.false. : data do    change during time
     logical :: steadyTD
     !> steadyTD_written=.true.  : closing time     written
     !> SteadyTD_written=.false. : closing time not written
     logical :: steadyTD_written=.false.
   contains
     !> static constructor
     !> (procedure public for type TimeData)
     procedure, public, pass :: init => init_TD
     !> static destructor
     !> (procedure public for type TimeData)
     procedure, public, pass :: kill => kill_TD
     !> Info procedure.
     !> (public procedure for type TimeData)
     procedure, public, pass :: info => info_TD
     !> Set (read if necessary) time data
     procedure, public, pass :: set => set_TD
     !> Set number of non zero element
     procedure, public, pass :: eval_ninput
  end type TimeData

contains
  !>-------------------------------------------------------------
  !> Static constructor.
  !> (procedure public for type TimeData)
  !> Instantiate (allocate)
  !> and initilize (by also reading from input file)
  !> variable of type TimeData
  !>
  !> usage:
  !>     call 'var'%init( InputFdescr, Ndata)
  !>
  !> where:
  !> \param[in] InputFdescr -> type(file) I/O file information for
  !>                           input variable 
  !> \param[in] Ndata -> number of data to be read in
  !>
  !<-------------------------------------------------------------------
  subroutine init_TD(this, stderr, InputFdescr, dimdata, Ndata)
    use Globals
    implicit none
    !vars
    class(TimeData),   intent(out) :: this
    integer,           intent(in ) :: stderr
    type(file),        intent(in ) :: InputFdescr
    integer, optional, intent(in ) :: dimdata
    integer, optional, intent(in ) :: Ndata
    ! local vars
    integer :: u_number
    integer :: NInput,ival
    integer :: res,i,k
    logical :: rc
    character(len=15) :: scratch
    character(len=256) :: str, fname

    this%built=.true.
    

    u_number = InputFdescr%lun
    fname    = InputFdescr%fn

    ! read file
    ! read dimdata, Ndata
    read(u_number,*,iostat=res) this%dimdata, this%Ndata
    if(res .ne. 0) &
         rc = IOerr(stderr, err_inp , 'read_TD', &
         trim(fname) // ' type TimeData member dimdata, Ndata ',res)

    ! check dimension
    if ( present (dimdata) ) then
       if ( dimdata .ne. this%dimdata) then
          rc = IOerr(stderr, err_inp , 'read_TD', &
               trim(fname) // ' mismatch between read and given dimdata')
       end if
    end if

    if ( present (Ndata) ) then
       if ( Ndata .ne. this%Ndata) then
          rc = IOerr(stderr, err_inp , 'read_TD', &
               trim(fname) // ' mismatch between read and given Ndata')
       end if
    end if

    ! allocate
    allocate(this%TDactual(this%dimdata,this%Ndata),stat=res)
    if(res .ne. 0) rc = IOerr(stderr, err_alloc, 'read_TD', &
         '  type TimeData member TDactual (array)',res)
    this%TDactual=zero
    
    allocate(this%TDval(this%dimdata,this%Ndata,2),stat=res)
    if(res .ne. 0) rc = IOerr(stderr, err_alloc, 'read_TD', &
         '  type TimeData member TDval (array)',res)
    this%TDval=zero
    
    allocate(this%TDtime(2),stat=res)
    if(res .ne. 0) rc = IOerr(stderr, err_alloc, 'read_TD', &
         '  type TimeData member TDtime (array)',res)
    this%TDtime=zero

    allocate(this%val(this%dimdata),stat=res)
    if(res .ne. 0) rc = IOerr(stderr, err_alloc, 'read_TD', &
         '  type TimeData member val (array)',res)
    
    this%TDtime=zero
    this%steadyTD=.false.


    ! read the first time for TD
    read(u_number,*,iostat=res) scratch, this%TDtime(1)
    if(res .ne. 0) &
         rc = IOerr(stderr, err_inp , 'read_TD', &
         trim(fname) // ' type TimeData member TDtime(1) ',res)

    ! read Ninputs
    read(u_number,*,iostat=res) NInput
    if(res .ne. 0) &
         rc = IOerr(stderr, err_inp , 'read_TD', &
         trim(fname) // ' type TimeData member NInput ',res)
    
    ! read Data
    do i = 1,NInput
       read(u_number,*,iostat=res) ival, (this%val(k), k=1,this%dimdata) 
       if(res .ne. 0) &
            rc = IOerr(stderr, err_inp , 'read_TD', &
            trim(fname) // ' type TimeDara member ival,val ',res)
       this%TDval(:,ival,1)=this%val(:)
    end do

    ! read the second time for TD
    read(u_number,*,iostat=res) scratch, this%TDtime(2)
    if(res .ne. 0) &
         rc = IOerr(stderr, err_inp , 'read_TD', &
         etb(trim(InputFdescr%fn)//' type TimeData member timein(2) '),res)
    

    if (this%TDtime(2).ge.huge)then
       ! if Time data is in steady state copy TD(:,:,1) into TD(:,:,2)

       this%TDactual(:,:)=this%TDval(:,:,1)
       this%TDval(:,:,2)=this%TDval(:,:,1)
       this%steadyTD=.true.
    else
       ! if not in stready state read TD(2,:,:)
       
       ! read Ninputs
       read(u_number,*,iostat=res) NInput
       if(res .ne. 0) &
            rc = IOerr(stderr, err_inp , 'read_TD', &
            trim(InputFdescr%fn) // ' type TimeData member NInput ',res)
       
       ! read Data
       do i = 1,NInput
          read(u_number,*,iostat=res) ival, (this%val(k), k=1,this%dimdata) 
          if(res .ne. 0) &
               rc = IOerr(stderr, err_inp , 'read_TD', &
               trim(InputFdescr%fn) // ' type TimeDara member ival,val ',res)
          this%TDval(:,ival,2)=this%val(:)
       end do
    end if

  end subroutine init_TD
  
  !>-------------------------------------------------------------
  !> Static destructor.
  !> (procedure public for type TimeData)
  !> deallocate all arrays for a var of type TimeData
  !>
  !> usage:
  !>     call 'var'%kill(lun)
  !>
  !> where:
  !> \param[in] lun -> integer. I/O unit for error message output.
  !<-----------------------------------------------------------
  subroutine kill_TD(this, lun)
    use Globals
    implicit none
    ! vars
    class(TimeData), intent(inout):: this
    integer, intent(in) :: lun
    ! local vars
    integer :: res
    logical :: rc
    
    this%built=.false.
    deallocate(this%TDval,this%TDtime,this%TDActual,this%val)
    if (res.ne.0) rc=IOerr(lun, err_dealloc, 'kill_TD', &
         ' inputs dealloc fail for TimeData var TDvar,TDtime,TDactual, val')
  end subroutine kill_TD
  
  !>-------------------------------------------------------------
  !> Info procedure.
  !> (public procedure for type TimeData)
  !> Prints content of a variable of type TimeData.
  !> Report if the time-range covered [t1,t2], actual time
  !> If (nsample .gt. 0) it prints the first nsample non-zero terms
  !> of TDval(1,:.:), TDactual(:,:) and TDval(2,:.:)   
  !> 
  !> usage:
  !>     call 'var'%info(lun,nsample)
  !>
  !> where:
  !> \param[in] lun -> integer. I/O unit for error message output.
  !> \param[in] lun -> integer. Number of first no zero samples to be printed
  !<-------------------------------------------------------------
  subroutine info_TD(this, lun, nsample)
    use Globals
    implicit none
    ! vars
    class(TimeData), intent(in) :: this
    integer,         intent(in) :: lun
    integer,         intent(in) :: nsample

    ! local vars
    integer :: i,j,k
    real(kind=double) :: dnrm2
        
    write(lun,*) ' '
    write(lun,*) ' Info: TimeData structure definition:'
    
    write(lun,*) 'ndata', this%ndata
    if( this%steadyTD       ) write(lun,*) ' Steady state'
    write(lun,*) ' t1          = ', this%TDtime(1)
    write(lun,*) ' actual time = ', this%time
    write(lun,*) ' t2          = ', this%TDtime(2)
    if ( nsample .gt. 0 ) then
       write(lun,*) ' First Data'
       i=0
       j=0
       do while ( (i .lt. this%ndata) .and. (j .lt. nsample) )
          i=i+1
          if( dnrm2(this%dimdata,this%TDval(:,i,1),1) .ne.zero ) then
             j=j+1
             write(lun,'(5(i5,e11.3))') i,(this%TDval(k,i,1),k=1,this%dimdata)
          end if
       end do
       write(lun,*)
       write(lun,*) ' Actual Data'
       i=0
       j=0
       do while ( (i .lt. this%ndata) .and. (j .lt. nsample) )
          i=i+1
          if( dnrm2(this%dimdata,this%TDactual(:,i),1).ne.zero ) then
             j=j+1
             write(lun,'(5(i5,e11.3))') i,(this%TDactual(k,i),k=1,this%dimdata)
          end if
       end do
       write(lun,*)
       write(lun,*) ' Second time  Data'
       i=0
       j=0
       do while ( (i .lt. this%ndata) .and. (j .lt. nsample) )
          i=i+1
          if( dnrm2(this%dimdata,this%TDval(:,i,2),1) .ne.zero ) then
             j=j+1
             write(lun,'(5(i5,e11.3))') i,(this%TDval(k,i,2),k=1,this%dimdata)
          end if
       end do
    end if
  end subroutine info_TD

  !>-------------------------------------------------------------
  !> Set (read if necessary) Time dependent variables
  !> (public procedure for type TimeData) at given 
  !> time 
  !> 
  !> usage:
  !>     call 'var'%set(stderr, InputFdescr, time)
  !>
  !> where:
  !> \param[in] stderr      -> integer. I/O err msg unit. 
  !> \param[in] InputFdescr -> type(file) I/O file information for
  !>                           input variable 
  !> \param[in] time        -> real(double). Require time
  !<-------------------------------------------------------------
  subroutine set_TD(this, stderr, InputFdescr, time,endfile)
    use Globals
    implicit none
    ! vars
    class(TimeData),   intent(inout) :: this
    integer,           intent(in   ) :: stderr
    type(file),        intent(in   ) :: InputFdescr
    real(kind=double), intent(in   ) :: time
    logical,           intent(inout) :: endfile
    ! local vars
    integer :: res, u_number
    integer :: NInput
    integer :: i,k,ival
    real(kind=double) ::  TDt1,TDt2,tperc,next_time
    logical :: rc, read_next
    character(len=15)  :: scratch
    character(len=256) :: fname

    endfile=.false.

    this%time = time
    u_number = InputFdescr%lun
    fname    = InputFdescr%fn
    
    if (.not. this%steadyTD) then
       if ( time .ge. this%TDtime(2) ) then
          read_next = .true.
          do while ( read_next )
             ! Read time2 
             read(u_number,*,iostat=res) scratch, next_time
             if(res .ne. 0) then
                if ( res .eq. -1 ) then
                   endfile=.true.
                   return
                else
                   rc = IOerr(stderr, err_inp , 'set_TD', &
                        trim(fname)&
                        //' type TimeData member TDtime(2)',res)
                end if
             else
                this%TDtime(1)=this%TDtime(2)
                this%TDtime(2)=next_time
             end if
             ! copy TD2 into TD1 before overwritten
             this%TDval(:,:,1)=this%TDval(:,:,2)

             ! If steady state (time2 .ge. huge) then frezee
             if ( this%TDtime(2) .ge. huge ) then
                this%steadyTD = .true.
                this%TDval(:,:,2)  = this%TDval(:,:,1)
                this%TDactual(:,:) = this%TDval(:,:,1)
                read_next = .false.
             else
                ! Otherwise read new data TD2
                
                ! read ninputs
                read(u_number,*,iostat=res) NInput
                if(res .ne. 0) &
                     rc = IOerr(stderr, err_inp , 'set_TD', &
                     trim(fname)//' type TimeData member NInput ',res)
                
                ! read data
                do i = 1,NInput
                   read(u_number,*,iostat=res) ival, (this%val(k),k=1,this%dimdata)
                   if(res .ne. 0) &
                        rc = IOerr(stderr, err_inp , 'set_TD', &
                        trim(fname)&
                        //' type TimeData aux varibles i,val',res) 
                   this%TDval(:,ival,2) = this%val(:)
                end do
                !> Test if continue reading file 
                read_next = ( time .ge. this%TDtime(2) ) 
             end if
          end do
       end if
       
       TDt1=this%TDtime(1)
       TDt2=this%TDtime(2)
       tperc=(time-TDt1)/(TDt2-TDt1)


       
       do i=1,this%ndata
          this%TDactual(:,i)=(one-tperc)*this%TDval(:,i,1)+&
               tperc*this%TDval(:,i,2)
       end do

       ! EF TODO replace with
!!$       call daxpy(this%dimdata*this%ndata, &
!!$            (one-tperc),this%TDval(1,1,1),1,&
!!$            this%TDactual(1,1),1)
!!$       call daxpy(this%dimdata*this%ndata, &
!!$            tperc,this%TDval(1,1,2),1,&
!!$            this%TDactual(1,1),1)
    end if
  end subroutine set_TD

  function eval_ninput(this) result(ninput)
    use Globals
    implicit none
    class(TimeData), intent(in) :: this
    integer :: ninput
    !local
    integer :: i
    real(kind=double) :: dnrm2
    ninput= 0
    do i = 1, this%NData
       if ( dnrm2(this%dimdata,this%TDactual(:,i),1) > small ) then 
          ninput=ninput+1
       end if
    end do
  end function eval_ninput

end module TimeInputs
       
  
