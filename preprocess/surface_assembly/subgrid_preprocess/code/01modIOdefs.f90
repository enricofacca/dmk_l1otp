!> Defines I/O units, directories, filenames
module IOdefs
  use Globals
  implicit none
  private
  type, public :: IOfdescr
     !> max unit number used in the code.
     integer :: maxUsedLun

     !> terminal (generally stdout)
     type(file) :: term
     !> obvious defs
     type(file) :: stdin,stdout,stderr


     !> Dynamic var used to avoid conflicts in opening new units.
     !>---------------------------------------------------------
     !> Inputs 
     !>---------------------------------------------------------
     !> General input directory
     type(file) :: input_folder
     !>---------------------------------------------------------
     !> Mandatory Input file
     !>---------------------------------------------------------
     !> input file for grid
     type(file) :: grid
     !> input file for piecewise source term
     type(file) :: source
     !> input file for piecewise sink term
     type(file) :: sink
     !> input file for dirac source term
     type(file) :: dirac_source
     !> input file for dirac sink term
     type(file) :: dirac_sink
     !> input file for output boundary flux
     !type(file) :: boundary_flux
     !>--------------------------------------------------------- 
     !> Outputs 
     !>--------------------------------------------------------- 
     !> output subgrid 
     type(file) :: subgrid
     !> output file of grid/subgrid relation
     type(file) :: parent
     !> output file of $\Forcing$ integrated w.r.t.
     !> subgrid p1-basis fucntions
     type(file) :: rhs_subgrid_integrated
     !> output file of $\Forcing$ integrated w.r.t.
     !> grid p1-basis fucntions
     type(file) :: rhs_grid_integrated
!!$     !> generic debug file
!!$     type(file) :: debug
   contains
     !> static constructor
     !> (public for type IOfdescr)
     procedure, public, pass :: init => IO_initialize
     !> output the content of the variable
     !> (public for type IOfdescr)
     procedure, public, pass :: info => IO_print_info 
     !> static destructor
     !> (public for type IOfdescr)
     procedure, public, pass :: kill => IO_close_all_units
     !> Saving real data from start to finsih
     !> (public for type IOfdescr)
     procedure, public, nopass :: save_rdatafromto
     !> Saving real data from start to finsih
     !> (public for type IOfdescr)
     procedure, public, nopass :: save_idatafromto
     !> Saving real data from start to finsih
     !> (public for type IOfdescr)
     procedure, public, nopass :: rwrite     
  end type IOfdescr

contains
  !>-------------------------------------------------------------
  !> Static constructor. 
  !> (public procedure for type IOfdescr)
  !> Instantiate and initialize a variable of type IOfdescr
  !>
  !> usage:
  !>     call 'var'\%init()
  !>
  !> - defines stdin,stderr,stdout
  !> - reads from "preprocess.fnames" working directory and filenames
  !> - open all file units checking for consistency
  !<-------------------------------------------------------------
  subroutine IO_initialize(this)
    use Globals
    use, intrinsic :: iso_fortran_env, only : stdin=>input_unit, &
         stdout=>output_unit, &
         stderr=>error_unit
    class(IOfdescr) :: this

    ! local vars
    integer :: opstat,nlun
    logical :: fv_exist=.false.,file_exist=.false.,dir_exist
    character (len=256) :: fname='subgrid.fnames', dir
    !character (len=256) :: add_msg,rdwr

    integer :: lun 
    
    this%stderr%exist = .true.
    this%stderr%lun=stderr
    this%stderr%fn='stderr'
    
    this%stdin%exist = .true.
    this%stdin%lun=stdin
    this%stdin%fn='stdin'
    
    this%stdout%exist = .true.
    this%stdout%lun=stdout
    this%stdout%fn='stdout'
    
    this%term%exist = .true.
    this%term%lun=this%stdout%lun
    this%term%fn='terminal'

    nlun=6

    ! Want to avoid default lun 0, 5 and 6
    ! as they usually are stderr,stdin,stdout
    
    nlun=nlun+1
    this%grid%exist = .false.
    this%grid%lun = nlun
    this%grid%fn = 'grid'


    nlun=nlun+1
    this%source%exist = .false.
    this%source%lun = nlun
    this%source%fn = 'source'


    nlun=nlun+1
    this%sink%exist = .false.
    this%sink%lun = nlun
    this%sink%fn = 'sink'

    nlun=nlun+1
    this%dirac_source%exist = .false.
    this%dirac_source%lun = nlun
    this%dirac_source%fn = 'dirac_source'


    nlun=nlun+1
    this%dirac_sink%exist = .false.
    this%dirac_sink%lun = nlun
    this%dirac_sink%fn = 'dirac_sink'

    
!!$    nlun=nlun+1
!!$    this%boundary_flux%exist = .false.
!!$    this%boundary_flux%lun = nlun
!!$    this%boundary_flux%fn = 'boundary_flux' 


    ! output 
    nlun=nlun+1
    this%subgrid%exist = .false.
    this%subgrid%lun = nlun
    this%subgrid%fn = 'grid'



    nlun=nlun+1
    this%parent%exist = .false.
    this%parent%lun = nlun
    this%parent%fn = 'parent'

    nlun=nlun+1
    this%rhs_subgrid_integrated%exist = .false.
    this%rhs_subgrid_integrated%lun = nlun
    this%rhs_subgrid_integrated%fn = 'rhs_subgrid_integrated'

    nlun=nlun+1
    this%rhs_grid_integrated%exist = .false.
    this%rhs_grid_integrated%lun = nlun
    this%rhs_grid_integrated%fn = 'rhs_grid_integrated'

    
!!$    nlun=nlun+1
!!$    this%debug%exist = .false.
!!$    this%debug%lun=nlun
!!$    this%debug%fn ='debug'

    ! remember to change this if we add a lun
    ! this is used to get free luns
    this%maxUsedLun=nlun


    ! read in from standard file (or input) the IO control file
    ! that overrides above default values
    do while (.not. fv_exist)
       inquire(file=fname,exist=fv_exist)
       if (fv_exist) then
          open(1,file=fname)
          read(1,'(a)') dir
          dir=etb(erase_comment(dir))
          ! inputs
          read(1,'(a)') this%grid%fn
          read(1,'(a)') this%source%fn
          read(1,'(a)') this%sink%fn
          read(1,'(a)') this%dirac_source%fn
          read(1,'(a)') this%dirac_sink%fn
          !read(1,'(a)') this%boundary_flux%fn
          

          ! output data
          read(1,'(a)') this%subgrid%fn
          read(1,'(a)') this%parent%fn
          read(1,'(a)') this%rhs_subgrid_integrated%fn
          read(1,'(a)') this%rhs_grid_integrated%fn
          !read(1,'(a)') this%debug%fn
          close(1)

          inquire(file=dir,exist=dir_exist)
          if(.not.dir_exist) then
             file_exist = IOerr(stderr, wrn_IO,'IOinitialize',etb(dir))
             cycle
          end if
         
          ! prepend character DIR to all filenames
          ! inputs
          call prepend_dir(this%grid,dir)
          call prepend_dir(this%source,dir)
          call prepend_dir(this%sink,dir)
          call prepend_dir(this%dirac_source,dir)
          call prepend_dir(this%dirac_sink,dir)
          !call prepend_dir(this%boundary_flux,dir)
          
          ! output
          call prepend_dir(this%subgrid,dir)
          call prepend_dir(this%parent,dir)
          call prepend_dir(this%rhs_subgrid_integrated,dir)
          call prepend_dir(this%rhs_grid_integrated,dir)
          !call prepend_dir(this%debug,dir)

          lun=this%stderr%lun
          !>--------------------------------------------------------------------
          ! Check files existence and open units
          ! inputs
          call check_open(this%grid,stderr,1,1,1)
          call check_open(this%source,stderr,1,1,1)
          call check_open(this%sink,stderr,1,1,1)
          call check_open(this%dirac_source,stderr,1,1,1)
          call check_open(this%dirac_sink,stderr,1,1,1)
          !call check_open(this%boundary_flux,stderr,1,2,1)         
          ! optional inputs
          call check_open(this%subgrid,stderr,1,1,2)
          call check_open(this%parent,stderr,1,1,2)
          call check_open(this%rhs_subgrid_integrated,stderr,1,1,2)
          call check_open(this%rhs_grid_integrated,stderr,1,1,2)
          !call check_open(this%debug,stderr,1,1,2)
  
       else
          write(lun,*) ' Control file ',etb(fname),&
               ' does not exist in the current directory!!'
          if (lun.eq.this%stdin%lun) then
             write(lun,*) ' input another file name or END/QUIT to stop'
             read(*,'(a)') fname
             fname=to_lower(fname)
             if ( &
                  fname.eq.'end' .or. &
                  fname.eq.'quit' .or. &
                  fname.eq.'q' .or. &
                  fname.eq.'e' &
                  ) then
                stop ' simulation ended'
             end if
             ! everything has failed, restart loop
          else 
             stop ' simulation ended'
          end if
       end if
    end do
  contains
    subroutine prepend_dir(file_name,dir)
      use Globals
      implicit none
      type(file),       intent(inout) :: file_name
      character(len=*), intent(in   ) ::dir
      
      file_name%fn=etb(dir)//'/'//etb(erase_comment(file_name%fn))
    end subroutine prepend_dir

    subroutine check_open(file2read,lun,&
         file_folder_flag,&
         mandatory_optional_flag,&
         input_output_flag)
      use Globals
      implicit none
      type(file), intent(inout) :: file2read
      integer,    intent(in   ) :: lun
      integer,    intent(in   ) :: file_folder_flag
      integer,    intent(in   ) :: mandatory_optional_flag
      integer,    intent(in   ) :: input_output_flag

      ! rhs_subgrid_integrated
      ! local
      logical :: file_exist
      
      if   ( file_folder_flag .eq. 1 ) then
         if ( input_output_flag .eq. 1 ) then
            inquire(file=file2read%fn,exist=file_exist)
            if ( file_exist .and. mandatory_optional_flag .eq. 1 ) then 
               file2read%exist = .true.
               open(file2read%lun,file=file2read%fn,iostat=opstat)
               if(opstat.ne.0) then
                  write(lun,*) ' failed to open unit ',file2read%lun,&
                       ' to be linked to file ',etb(file2read%fn)
                  stop ' simulation ended'
               end if
            else 
               write(lun,*) 'File ', etb(file2read%fn), ' not found'
               write(lun,*) 'Optional for the simulation'
               write(lun,*) 'Continue'
            end if
         else if ( input_output_flag .eq. 2 )  then
            if ( mandatory_optional_flag .eq. 1 ) then 
               open(file2read%lun,file=file2read%fn,iostat=opstat)
               if(opstat.ne.0) then
                  write(lun,*) ' failed to open unit ',file2read%lun,&
                       ' to be linked to file ',etb(file2read%fn)
                  stop ' simulation ended'
               else
                  file2read%exist = .true.
               end if
            end if
         end if
      else if ( file_folder_flag .eq. 2 ) then        
         inquire(file=etb(file2read%fn),exist=file_exist)
         if(.not.dir_exist) then
            call execute_command_line ('mkdir -p ' // etb(file2read%fn))
         end if
         file2read%exist= .true.
         open(file2read%lun,file=etb(file2read%fn)//'/deleteme',&
              iostat=opstat)
         if (opstat.ne.0) then
            write(lun,*) ' problems creating or accessing input_folder directory ', &
                 etb(file2read%fn)
            stop ' simulation ended'
         end if
         close(file2read%lun)
         call execute_command_line ('rm -f ' // etb(file2read%fn)//'/deleteme')  
      end if

    end subroutine check_open

  end subroutine IO_initialize

  !-------------------------------------------------------------
  !> Info procedure.
  !> (public procedure for type IOfdescr)
  !> Prints content of a variable of type IOfdescr
  !>
  !> usage:
  !>    call 'var'%info(lun)
  !>
  !> where:
  !> \param[in] lun: output unit
  !-------------------------------------------------------------
  subroutine IO_print_info(this, lun)
    class(IOfdescr) :: this
    integer :: lun

    write(lun,*) ' IO Units used in the code'
    
    write(lun,*) ' '
    write(lun,*) ' Standard units'
    call this%term%info(lun)
    call this%stderr%info(lun)
    call this%stdin%info(lun)
    call this%stdout%info(lun)

    write(lun,*) ' '
    write(lun,*) ' Inputs units'
    call this%input_folder%info(lun)
    !
    call this%grid%info(lun)
    call this%source%info(lun)
    call this%sink%info(lun)
    call this%dirac_source%info(lun)
    call this%dirac_sink%info(lun)
    !call this%boundary_flux%info(lun)
 
    write(lun,*) ' Output units'
    call this%subgrid%info(lun)
    call this%parent%info(lun)
    call this%rhs_subgrid_integrated%info(lun)
    call this%rhs_grid_integrated%info(lun)
    !call this%debug%info(lun)

  end subroutine IO_print_info
  
  !-------------------------------------------------------------
  !> Static destructor.
  !> (public procedure for type IOfdescr)
  !> Deallocate a variable of type IOfdescr and close I/O units
  !>
  !> usage:
  !>    call 'var'%kill(lun)
  !>
  !> where:
  !> \param[in] lun: output unit for error messages
  !-------------------------------------------------------------
  subroutine IO_close_all_units(this)
    class(IOfdescr) :: this
    ! input units
    close(1)
    close (this%term%lun)
    close (this%stderr%lun)
    close (this%stdin%lun)
    close (this%stdout%lun)

    !
    close (this%grid%lun)
    close (this%source%lun)
    close (this%sink%lun)
    close (this%dirac_source%lun)
    close (this%dirac_sink%lun)
    !close (this%boundary_flux%lun)
    !
    close (this%subgrid%lun)
    close (this%parent%lun)
    close (this%rhs_subgrid_integrated%lun)
    close (this%rhs_grid_integrated%lun)
    !close (this%debug%lun)

    
  end subroutine IO_close_all_units  
  
  !>---------------------------------------------------------
  !> Saving procedure for saving real data array from start 
  !> to finish into given directory with name 'datanme'.dat
  !> (procedure public for type IOfdescr)
  !>
  !> usage: call savedatafromto(lun,&
  !>                       size,start,finish,&
  !>                       dir,dataname,data)
  !> where
  !> \param[in] lun      -> integer. Logic unit 
  !> \param[in] size     -> integer. Size of the array of data
  !> \param[in] start    -> integer. Starting save position
  !> \param[in] finish   -> integer. Final save position
  !> \param[in] dir      -> character. Directory name
  !> \param[in] dataname -> character. Data label
  !> \param[in] data     -> real. Data array
  !<----------------------------------------------------------
  subroutine save_rdatafromto(lun,&
       size,start,finish,&
       dir,dataname,data)
    use Globals
    implicit none
    integer, intent(in) :: lun
    integer, intent(in) :: size, start, finish
    character(len=*), intent(in) :: dir,dataname
    real(kind=double), intent(in) :: data(size)
    !local
    integer :: isav
    character(256) :: fname

    fname=etb(etb(dir)//'/'//etb(dataname)//'.dat')
    open(lun,file=fname,access='append')
    do isav=start,finish
       write(lun,'(1e15.6)')  data(isav)
    end do
    close(lun)
  end subroutine save_rdatafromto

  !>---------------------------------------------------------
  !> Saving procedure for saving integer data array from start 
  !> to finish into given directory with name 'datanme'.dat
  !> (procedure public for type IOfdescr)
  !>
  !> usage: call save_idatafromto(lun,&
  !>                       size,start,finish,&
  !>                       dir,dataname,data)
  !> where
  !> \param[in] lun      -> integer. Logic unit 
  !> \param[in] size     -> integer. Size of the array of data
  !> \param[in] start    -> integer. Starting save position
  !> \param[in] finish   -> integer. Final save position
  !> \param[in] dir      -> character. Directory name
  !> \param[in] dataname -> character. Data label
  !> \param[in] data     -> integer. Data array
  !<----------------------------------------------------------
  subroutine save_idatafromto(lun,&
       size,start,finish,&
       dir,dataname,data)
    use Globals
    implicit none
    integer, intent(in) :: lun
    integer, intent(in) :: size, start, finish
    character(len=*), intent(in) :: dir,dataname
    integer, intent(in) :: data(size)
    !local
    integer :: isav
    character(len=256) :: fname

    fname=etb(etb(dir)//'/'//etb(dataname)//'.dat')
    open(lun,file=fname,access='append')
    do isav=start,finish
       write(lun,'(I8)')  data(isav)
    end do
    close(lun)
  end subroutine save_idatafromto



  subroutine rwrite(&
       lun,&
       lun_err,&
       ndata,&
       current_iteration, &
       fname,&
       label_grid,&
       label_data,&
       dataname,&
       current_time,&
       rdata)
    use Globals
    implicit none
    integer,            intent(in) :: lun
    integer,            intent(in) :: lun_err
    integer,            intent(in) :: ndata
    integer,            intent(in) :: current_iteration
    character(len=*),   intent(in) :: fname
    character(len=*),   intent(in) :: label_grid
    character(len=*),   intent(in) :: label_data
    character(len=*),   intent(in) :: dataname
    real(kind=double),  intent(in) :: current_time
    real(kind=double),  intent(in) :: rdata(ndata)

    ! local
    integer :: i, opstat
    
    open(lun,file=etb(fname),iostat=opstat)
    if(opstat.ne.0) then
       write(lun_err,*) ' failed to open unit ',lun,&
            ' to be linked to file ',etb(fname)
       stop ' simulation ended'
    end if
    write(*,*) etb(label_grid),&
         current_iteration, &
         current_time
    
    write(lun,*) etb(label_data), ndata 
    write(lun,*) 'SCALARS', etb(dataname),' double 1'
    write(lun,*) 'LOOKUP_TABLE default' 
    do i=1,ndata
       write(lun,*) rdata(i)
    end do
    close(lun)


  end subroutine rwrite

!!$  subroutine write_head(&
!!$       lun,&
!!$       ndata,&
!!$       current_iteration, &
!!$       current_time,&
!!$       current_var,&
!!$       label_grid,&
!!$       label_data,&
!!$       dataname)
!!$
!!$    use Globals
!!$    implicit none
!!$    integer,            intent(in   ) :: lun
!!$    integer,            intent(in   ) :: lun_err
!!$    integer,            intent(in   ) :: ndata
!!$    integer,            intent(in   ) :: current_iteration
!!$    character(len=*),   intent(in   ) :: fname,label_grid,label_data,dataname
!!$    real(kind=double),  intent(in   ) :: current_time
!!$    real(kind=double), optional,  intent(in   ) :: current_var
!!$    ! local
!!$    integer :: i, opstat
!!$    
!!$    if ( present(current_var) ) then
!!$       write(lun,*) etb(label_grid),&
!!$            current_iteration, &
!!$            current_time,&
!!$            curernt_var
!!$    else
!!$       write(lun,*) etb(label_grid),&
!!$            current_iteration, &
!!$            current_time
!!$    end if
!!$          
!!$    write(lun,*) etb(label_data), ndata 
!!$    write(lun,*) 'SCALARS', etb(dataname),' double 1'
!!$    write(lun,*) 'LOOKUP_TABLE default' 
    
!!$  end subroutine write_head
  
  
    
end module IOdefs
