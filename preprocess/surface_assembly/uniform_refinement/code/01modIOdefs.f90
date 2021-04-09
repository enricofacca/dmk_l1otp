!> Defines I/O units, directories, filenames
module IOdefs
  implicit none 
 
  private
  type, public :: file
     logical :: exist          !< Logical of existence
     integer :: lun            !< I/O unit number
     character (len=256) :: fn !< I/O file name or directory name
   contains
     !> output the content of the variable
     !> (public for type file)
     procedure, public, pass :: info => fn_print
  end type file
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
     !> input file for control
     type(file) :: control
     !> input file for grid
     type(file) :: grid
     !> input file for building initial data
     type(file) :: initial
     !>---------------------------------------------------------
     !> Optional Input files
     !>---------------------------------------------------------
     !> input for the muffa ode
     !> input file for power on flux
     type(file) :: pflux
     !> input file for power on mass
     type(file) :: pmass
     !> input file for time decay
     type(file) :: decay
     !> input file for spatial decay
     type(file) :: kappa
     !> input file for initial tdens
     type(file) :: tdens0
     !> input file for initial tdens
     type(file) :: optdens
     !> input file of $\Forcing$ integrated w.r.t.
     !> subgrid p1-basis fucntions
     type(file) :: rhs_integrated
     !> input file for forcing term ode
     type(file) :: forcing
     !> input file for output boundary flux
     type(file) :: boundary_flux     

     !>--------------------------------------------------------- 
     !> Outputs 
     !>--------------------------------------------------------- 
     !> General output directory
     type(file) :: output_folder
     !>--------------------------------------------------------- 
     !> Mandatory Output file (Necessary for muffa)
     !>--------------------------------------------------------- 
     !> Geomtery output
     !> Output file for subgrid
     type(file) :: grid_out
     !> Output file for subgrid
     type(file) :: subgrid
     !> Output file for grid/subgrid conenctions
     type(file) :: sons
     !>--------------------------------------------------------- 
     !> Output for the muffa ode
     !> output file for power on flux
     type(file) :: pflux_out
     !> output file for power on mass
     type(file) :: pmass_out
     !> output file for time decay
     type(file) :: decay_out
     !> output file for spatial decay
     type(file) :: kappa_out
     !> output file for initial tdens
     type(file) :: tdens0_out
     !> output file of $\Forcing$ integrated w.r.t.
     !> subgrid p1-basis fucntions
     type(file) :: rhs_integrated_out
     !>---------------------------------------------------------
     !> Optional output files
     !> output file for optimal tdens
     type(file) :: optdens_out
     !> output file for forcing term ode
     type(file) :: forcing_out
     !> output file for output boundary flux
     type(file) :: boundary_flux_out
     !>---------------------------------------------------------
     !> generic debug file
     type(file) :: debug
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
    use Globals
    class(file) :: this
    integer, intent(in) :: lun

    write(lun,'(a,a,a,i3)') ' filename ',etb(this%fn), &
         ' linked to lun ',this%lun

  end subroutine fn_print

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
    character (len=256) :: fname='preprocess.fnames', dir
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
    
    ! General input folder
    nlun=nlun+1
    this%input_folder%exist = .false.
    this%input_folder%lun = nlun
    this%input_folder%fn = 'input_folder'

    ! controls 
    nlun=nlun+1
    this%control%exist = .false.
    this%control%lun = nlun
    this%control%fn = 'control'

    nlun=nlun+1
    this%grid%exist = .false.
    this%grid%lun = nlun
    this%grid%fn = 'grid'

    nlun=nlun+1
    this%sons%exist = .false.
    this%sons%lun = nlun
    this%sons%fn = 'sons'

    nlun=nlun+1
    this%initial%exist = .false.
    this%initial%lun = nlun
    this%initial%fn = 'initial'
    
    ! optional input
    nlun=nlun+1
    this%pflux%exist = .false.
    this%pflux%lun = nlun
    this%pflux%fn = 'pflux'

    nlun=nlun+1
    this%pmass%exist = .false.
    this%pmass%lun = nlun
    this%pmass%fn = 'pmass'

    nlun=nlun+1
    this%decay%exist = .false.
    this%decay%lun = nlun
    this%decay%fn = 'decay'

    nlun=nlun+1
    this%tdens0%exist = .false.
    this%tdens0%lun = nlun
    this%tdens0%fn = 'tdens0'

    nlun=nlun+1
    this%kappa%exist = .false.
    this%kappa%lun = nlun
    this%kappa%fn = 'kappa'    

    nlun=nlun+1
    this%rhs_integrated%exist = .false.
    this%rhs_integrated%lun = nlun
    this%rhs_integrated%fn = 'rhs_integrated'


    nlun=nlun+1
    this%optdens%exist = .false.
    this%optdens%lun = nlun
    this%optdens%fn = 'optdens'

    nlun=nlun+1
    this%forcing%exist = .false.
    this%forcing%lun = nlun
    this%forcing%fn = 'forcing'

    nlun=nlun+1
    this%boundary_flux%exist = .false.
    this%boundary_flux%lun = nlun
    this%boundary_flux%fn = 'boundary_flux' 
    
    !>--------------------------------------------------------------------------
    !> Output files  
    nlun=nlun+1
    this%output_folder%exist = .false.
    this%output_folder%lun = nlun
    this%output_folder%fn = 'output_folder'

    nlun=nlun+1
    this%grid_out%exist = .false.
    this%grid_out%lun = nlun
    this%grid_out%fn = 'grid_out'
    
    nlun=nlun+1
    this%subgrid%exist = .false.
    this%subgrid%lun = nlun
    this%subgrid%fn = 'subgrid'

    nlun=nlun+1
    this%sons%exist = .false.
    this%sons%lun = nlun
    this%sons%fn = 'sons'

   
    nlun=nlun+1
    this%pflux_out%exist = .false.
    this%pflux_out%lun = nlun
    this%pflux_out%fn = 'pflux_out'

    nlun=nlun+1
    this%pmass_out%exist = .false.
    this%pmass_out%lun = nlun
    this%pmass_out%fn = 'pmass_out'

    nlun=nlun+1
    this%decay_out%exist = .false.
    this%decay_out%lun = nlun
    this%decay_out%fn = 'decay_out'

    nlun=nlun+1
    this%tdens0_out%exist = .false.
    this%tdens0_out%lun = nlun
    this%tdens0_out%fn = 'tdens0_out'

    nlun=nlun+1
    this%kappa_out%exist = .false.
    this%kappa_out%lun = nlun
    this%kappa_out%fn = 'kappa_out'    

    nlun=nlun+1
    this%rhs_integrated_out%exist = .false.
    this%rhs_integrated_out%lun = nlun
    this%rhs_integrated_out%fn = 'rhs_integrated_out'


    nlun=nlun+1
    this%optdens_out%exist = .false.
    this%optdens_out%lun = nlun
    this%optdens_out%fn = 'optdens_out'

    nlun=nlun+1
    this%forcing_out%exist = .false.
    this%forcing_out%lun = nlun
    this%forcing_out%fn = 'forcing_out'

    nlun=nlun+1
    this%boundary_flux_out%exist = .false.
    this%boundary_flux_out%lun = nlun
    this%boundary_flux_out%fn = 'boundary_flux_out' 
    
    nlun=nlun+1
    this%debug%exist = .false.
    this%debug%lun=nlun
    this%debug%fn ='debug'

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
          read(1,'(a)') this%input_folder%fn
          ! mandatory
          read(1,'(a)') this%control%fn
          read(1,'(a)') this%grid%fn
          read(1,'(a)') this%initial%fn
          ! optionals
          read(1,'(a)') this%pflux%fn
          read(1,'(a)') this%pmass%fn
          read(1,'(a)') this%decay%fn
          read(1,'(a)') this%kappa%fn
          read(1,'(a)') this%tdens0%fn
          read(1,'(a)') this%optdens%fn          
          read(1,'(a)') this%rhs_integrated%fn
          read(1,'(a)') this%forcing%fn
          read(1,'(a)') this%boundary_flux%fn
          

          ! output data
          read(1,'(a)') this%output_folder%fn
          ! mandatory
          read(1,'(a)') this%grid_out%fn
          read(1,'(a)') this%subgrid%fn
          read(1,'(a)') this%sons%fn

          read(1,'(a)') this%pflux_out%fn
          read(1,'(a)') this%pmass_out%fn
          read(1,'(a)') this%decay_out%fn
          read(1,'(a)') this%kappa_out%fn
          read(1,'(a)') this%tdens0_out%fn
          read(1,'(a)') this%rhs_integrated_out%fn
          ! optionals
          read(1,'(a)') this%optdens_out%fn
          read(1,'(a)') this%forcing_out%fn
          read(1,'(a)') this%boundary_flux_out%fn
          read(1,'(a)') this%debug%fn
          close(1)

          inquire(file=dir,exist=dir_exist)
          if(.not.dir_exist) then
             file_exist = IOerr(stderr, wrn_IO,'IOinitialize',etb(dir))
             cycle
          end if
          ! prepend character DIR to all filenames
          ! inputs
!!$          this%input_folder%fn=etb(dir)//'/'//etb(erase_comment(this%input_folder%fn))
!!$          ! mandatoryinputs
!!$          this%control%fn=etb(dir)//'/'//etb(erase_comment(this%control%fn))
!!$          this%grid%fn=etb(dir)//'/'//etb(erase_comment(this%grid%fn))
!!$          this%initial%fn=etb(dir)//'/'//etb(erase_comment(this%initial%fn))
!!$          ! optional inputs
!!$          this%pflux%fn=etb(dir)//'/'//etb(erase_comment(this%pflux%fn))
!!$          this%pmass%fn=etb(dir)//'/'//etb(erase_comment(this%pmass%fn))
!!$          this%decay%fn=etb(dir)//'/'//etb(erase_comment(this%decay%fn))
!!$          this%kappa%fn=etb(dir)//'/'//etb(erase_comment(this%kappa%fn))
!!$          this%tdens0%fn=etb(dir)//'/'//etb(erase_comment(this%tdens0%fn))
!!$          this%optdens%fn=etb(dir)//'/'//etb(erase_comment(this%optdens%fn))
!!$          this%forcing%fn=etb(dir)//'/'//etb(erase_comment(this%forcing%fn))
!!$          this%boundary_flux%fn= etb(dir)//'/'//etb(erase_comment(this%boundary_flux%fn))
!!$          this%rhs_integrated%fn= etb(dir)//'/'//etb(erase_comment(this%rhs_integrated%fn))
!!$          
!!$          ! output
!!$          this%output_folder%fn=etb(dir)//'/'//etb(erase_comment(this%output_folder%fn))
!!$          ! output mandatory 
!!$          this%grid_out%fn=etb(dir)//'/'//etb(erase_comment(this%grid_out%fn))
!!$          this%subgrid%fn=etb(dir)//'/'//etb(erase_comment(this%subgrid%fn))
!!$          this%sons%fn=etb(dir)//'/'//etb(erase_comment(this%sons%fn))
!!$          !
!!$          this%pflux_out%fn=etb(dir)//'/'//etb(erase_comment(this%pflux_out%fn))
!!$          this%pmass_out%fn=etb(dir)//'/'//etb(erase_comment(this%pmass_out%fn))
!!$          this%decay_out%fn=etb(dir)//'/'//etb(erase_comment(this%decay_out%fn))
!!$          this%kappa_out%fn=etb(dir)//'/'//etb(erase_comment(this%kappa_out%fn))
!!$          this%tdens0_out%fn=etb(dir)//'/'//etb(erase_comment(this%tdens0_out%fn))
!!$          
!!$          ! optional
!!$          this%forcing_out%fn=etb(dir)//'/'//etb(erase_comment(this%forcing_out%fn))
!!$          this%boundary_flux_out%fn= etb(dir)//'/'//etb(erase_comment(this%boundary_flux_out%fn))
!!$          this%rhs_integrated_out%fn= etb(dir)//'/'//etb(erase_comment(this%rhs_integrated_out%fn))
!!$
!!$          
!!$
!!$          ! optional output
!!$          this%optdens_out%fn=etb(dir)//'/'//etb(erase_comment(this%optdens_out%fn))
!!$          this%forcing_out%fn=etb(dir)//'/'//etb(erase_comment(this%forcing_out%fn))
!!$          this%boundary_flux_out%fn=etb(dir)//'/'//etb(erase_comment(this%boundary_flux_out%fn))
!!$          this%debug%fn=etb(dir)//'/'//etb(erase_comment(this%debug%fn))
          
          ! prepend character DIR to all filenames
          ! inputs
          call prepend_dir(this%input_folder,dir)
          ! mandatory inputs
          call prepend_dir(this%control,dir)
          call prepend_dir(this%grid,dir)
          call prepend_dir(this%initial,dir)
          ! optional inputs
          call prepend_dir(this%pflux,dir)
          call prepend_dir(this%pmass,dir)
          call prepend_dir(this%decay,dir)
          call prepend_dir(this%kappa,dir)
          call prepend_dir(this%tdens0,dir)
          call prepend_dir(this%optdens,dir)
          call prepend_dir(this%forcing,dir)
          call prepend_dir(this%boundary_flux,dir)
          call prepend_dir(this%rhs_integrated,dir)
          
          ! output
          call prepend_dir(this%output_folder,dir)
          ! output mandatory 
          call prepend_dir(this%grid_out,dir)
          call prepend_dir(this%subgrid,dir)
          call prepend_dir(this%sons,dir)
          !
          call prepend_dir(this%pflux_out,dir)
          call prepend_dir(this%pmass_out,dir)
          call prepend_dir(this%decay_out,dir)
          call prepend_dir(this%kappa_out,dir)
          call prepend_dir(this%tdens0_out,dir)
          call prepend_dir(this%rhs_integrated_out,dir)

          ! optional output
          call prepend_dir(this%optdens_out,dir)
          call prepend_dir(this%forcing_out,dir)
          call prepend_dir(this%boundary_flux_out,dir)
          call prepend_dir(this%debug,dir)

          lun=this%stderr%lun
          !>--------------------------------------------------------------------
          ! Check files existence and open units
          ! inputs
          call check_open(this%input_folder,stderr,2,1,1)
          ! mandatory inputs
          call check_open(this%control,stderr,1,1,1)
          call check_open(this%grid,stderr,1,1,1)
          call check_open(this%initial,stderr,1,1,1)
          ! optional inputs
          call check_open(this%pflux,stderr,1,2,1)
          call check_open(this%pmass,stderr,1,2,1)
          call check_open(this%decay,stderr,1,2,1)
          call check_open(this%kappa,stderr,1,2,1)
          call check_open(this%tdens0,stderr,1,2,1)
          call check_open(this%optdens,stderr,1,2,1)
          call check_open(this%forcing,stderr,1,2,1)
          call check_open(this%boundary_flux,stderr,1,2,1)
          call check_open(this%rhs_integrated,stderr,1,2,1)
          
          ! output
          call check_open(this%output_folder,stderr,2,1,2)
          ! output mandatory 
          call check_open(this%grid_out,stderr,1,1,2)
          call check_open(this%subgrid,stderr,1,1,2)
          call check_open(this%sons,stderr,1,1,2)
          !
          call check_open(this%pflux_out,stderr,1,1,2)
          call check_open(this%pmass_out,stderr,1,1,2)
          call check_open(this%decay_out,stderr,1,1,2)
          call check_open(this%kappa_out,stderr,1,1,2)
          call check_open(this%tdens0_out,stderr,1,1,2)
          call check_open(this%rhs_integrated_out,stderr,1,1,2)

          ! optional output
          call check_open(this%optdens_out,stderr,1,1,2)
          call check_open(this%forcing_out,stderr,1,1,2)
          call check_open(this%boundary_flux_out,stderr,1,1,2)
          call check_open(this%debug,stderr,1,1,2)
  

!!$          lun=this%stderr%lun
!!$          !>--------------------------------------------------------------------
!!$          !> inputs
!!$          ! input_folder
!!$          inquire(file=etb(this%input_folder%fn),exist=dir_exist)
!!$          if(.not.dir_exist) then
!!$             call execute_command_line ('mkdir -p ' // etb(this%input_folder%fn))
!!$             this%input_folder%exist= .true.
!!$          end if
!!$          open(this%input_folder%lun,file=etb(this%input_folder%fn)//'/deleteme',&
!!$               iostat=opstat)
!!$          if (opstat.ne.0) then
!!$             write(lun,*) ' problems creating or accessing input_folder directory ', &
!!$                  etb(this%input_folder%fn)
!!$             stop ' simulation ended'
!!$          end if
!!$          close(this%input_folder%lun)
!!$          call execute_command_line ('rm -f ' // etb(this%input_folder%fn)//'/deleteme')
!!$
!!$
!!$          ! control
!!$          inquire(file=this%control%fn,exist=file_exist)
!!$          if ( file_exist ) then 
!!$             this%control%exist = .true.
!!$             open(this%control%lun,file=this%control%fn,iostat=opstat)
!!$             if(opstat.ne.0) then
!!$                rc = IOerr(stderr,err_IO,'IOinitialize', &
!!$                     etb(this%control%fn), &
!!$                     this%control%lun)
!!$             end if
!!$          else
!!$             write(lun,*) 'File', etb(this%control%fn), 'not found'
!!$             write(lun,*) 'Necessary for the simulation'
!!$             stop 'Simulations ended'
!!$          end if
!!$
!!$
!!$          !grid
!!$          inquire(file=this%grid%fn,exist=file_exist)
!!$          if ( file_exist ) then 
!!$             this%grid%exist = .true.
!!$             open(this%grid%lun,file=this%grid%fn,iostat=opstat)
!!$             if(opstat.ne.0) then
!!$                write(lun,*) ' failed to open unit ',this%grid%lun,&
!!$                     ' to be linked to file ',etb(this%grid%fn)
!!$                stop ' simulation ended'
!!$             end if
!!$          else
!!$             write(lun,*) 'File', etb(this%grid%fn), 'not found'
!!$             write(lun,*) 'Necessary for the simulation'
!!$             stop 'Simulations ended'
!!$          end if
!!$
!!$          ! initial
!!$          inquire(file=this%initial%fn,exist=file_exist)
!!$          if ( file_exist ) then 
!!$             this%initial%exist = .true.
!!$             open(this%initial%lun,file=this%initial%fn,iostat=opstat)
!!$             if(opstat.ne.0) then
!!$                rc = IOerr(stderr,err_IO,'IOinitialize', &
!!$                     etb(this%initial%fn), &
!!$                     this%initial%lun)
!!$             end if
!!$          else
!!$             write(lun,*) 'File ', etb(this%initial%fn), ' not found'
!!$             write(lun,*) 'Necessary for the simulation'
!!$             stop 'Simulations ended'
!!$          end if
!!$         
!!$          ! pflux          
!!$          open(this%pflux%lun,file=this%pflux%fn,iostat=opstat)
!!$          if(opstat.ne.0) then
!!$             write(lun,*) ' failed to open unit ',this%pflux%lun,&
!!$                  ' to be linked to file ',etb(this%pflux%fn)
!!$             stop ' simulation ended'
!!$          else
!!$             write(lun,*) 'File ', etb(this%pflux%fn), ' not found'
!!$             write(lun,*) 'Optional for the simulation'
!!$             write(lun,*) 'Simulation conintues'
!!$             write(lun,*) 'Data will be built by the preprocess'
!!$          end if
!!$
!!$          ! pmass          
!!$          open(this%pmass%lun,file=this%pmass%fn,iostat=opstat)
!!$          if(opstat.ne.0) then
!!$             write(lun,*) ' failed to open unit ',this%pmass%lun,&
!!$                  ' to be linked to file ',etb(this%pmass%fn)
!!$             stop ' simulation ended'
!!$          else
!!$             write(lun,*) 'File ', etb(this%pmass%fn), ' not found'
!!$             write(lun,*) 'Optional for the simulation'
!!$             write(lun,*) 'Simulation conintues'
!!$             write(lun,*) 'Data will be built by the preprocess'
!!$
!!$          end if
!!$
!!$          ! decay          
!!$          open(this%decay%lun,file=this%decay%fn,iostat=opstat)
!!$          if(opstat.ne.0) then
!!$             write(lun,*) ' failed to open unit ',this%decay%lun,&
!!$                  ' to be linked to file ',etb(this%decay%fn)
!!$             stop ' simulation ended'
!!$          else
!!$             write(lun,*) 'File ', etb(this%decay%fn), ' not found'
!!$             write(lun,*) 'Optional for the simulation'
!!$             write(lun,*) 'Simulation conintues'
!!$             write(lun,*) 'Data will be built by the preprocess'
!!$          end if
!!$
!!$          ! kappa          
!!$          open(this%kappa%lun,file=this%kappa%fn,iostat=opstat)
!!$          if(opstat.ne.0) then
!!$             write(lun,*) ' failed to open unit ',this%kappa%lun,&
!!$                  ' to be linked to file ',etb(this%kappa%fn)
!!$             stop ' simulation ended'
!!$          else
!!$             write(lun,*) 'File ', etb(this%kappa%fn), ' not found'
!!$             write(lun,*) 'Optional for the simulation'
!!$             write(lun,*) 'Simulation conintues'
!!$             write(lun,*) 'Data will be built by the preprocess'
!!$          end if
!!$
!!$          ! optdens          
!!$          open(this%optdens%lun,file=this%optdens%fn,iostat=opstat)
!!$          if(opstat.ne.0) then
!!$             write(lun,*) ' failed to open unit ',this%optdens%lun,&
!!$                  ' to be linked to file ',etb(this%optdens%fn)
!!$             stop ' simulation ended'
!!$          else
!!$             write(lun,*) 'File ', etb(this%optdens%fn), ' not found'
!!$             write(lun,*) 'Optional for the simulation'
!!$             write(lun,*) 'Simulation conintues'
!!$             write(lun,*) 'Data will be built by the preprocess'
!!$          end if
!!$
!!$
!!$          ! optdens
!!$          open(this%optdens%lun,file=this%optdens%fn,iostat=opstat)
!!$          if(opstat.ne.0) then
!!$             write(lun,*) ' failed to open unit ',this%optdens%lun,&
!!$                  ' to be linked to file ',etb(this%optdens%fn)
!!$             write(lun,*) 'Optimal density does not exist' 
!!$          else
!!$             write(lun,*) 'File ', etb(this%optdens%fn), ' not found'
!!$             write(lun,*) 'Optional for the simulation'
!!$             write(lun,*) 'Simulation conintues'
!!$             write(lun,*) 'Data will be built by the preprocess'
!!$          end if
!!$          
!!$
!!$          ! rhs_integrated
!!$          inquire(file=this%rhs_integrated%fn,exist=file_exist)
!!$          if ( file_exist ) then 
!!$             this%rhs_integrated%exist = .true.
!!$             open(this%rhs_integrated%lun,file=this%rhs_integrated%fn,iostat=opstat)
!!$             if(opstat.ne.0) then
!!$                write(lun,*) ' failed to open unit ',this%rhs_integrated%lun,&
!!$                     ' to be linked to file ',etb(this%rhs_integrated%fn)
!!$                stop ' simulation ended'
!!$             end if
!!$          else
!!$             write(lun,*) 'File ', etb(this%rhs_integrated%fn), ' not found'
!!$             write(lun,*) 'Optional for the simulation'
!!$             write(lun,*) 'Simulation conintues'
!!$             write(lun,*) 'Data will be built by the preprocess'
!!$          end if
!!$
!!$
!!$          !>----------------------------------------------------------------------
!!$          !> Output Files 
!!$          !> grid_out        
!!$          open(this%grid_out%lun,file=this%grid_out%fn,iostat=opstat)
!!$          if(opstat.ne.0) then
!!$             write(lun,*) ' failed to open unit ',this%grid_out%lun,&
!!$                  ' to be linked to file ',etb(this%grid_out%fn)
!!$             stop ' simulation ended'
!!$          else
!!$             this%grid_out%exist = .true.
!!$          end if
!!$
!!$          !> subgrid        
!!$          open(this%subgrid%lun,file=this%subgrid%fn,iostat=opstat)
!!$          if(opstat.ne.0) then
!!$             write(lun,*) ' failed to open unit ',this%subgrid%lun,&
!!$                  ' to be linked to file ',etb(this%subgrid%fn)
!!$             stop ' simulation ended'
!!$          else
!!$             this%subgrid%exist = .true.
!!$          end if
!!$          
!!$          !> sons        
!!$          open(this%sons%lun,file=this%sons%fn,iostat=opstat)
!!$          if(opstat.ne.0) then
!!$             write(lun,*) ' failed to open unit ',this%sons%lun,&
!!$                  ' to be linked to file ',etb(this%sons%fn)
!!$             stop ' simulation ended'
!!$          else
!!$             this%sons%exist = .true.
!!$          end if
!!$
!!$          ! pflux          
!!$          open(this%pflux_out%lun,file=this%pflux_out%fn,iostat=opstat)
!!$          if(opstat.ne.0) then
!!$             write(lun,*) ' failed to open unit ',this%pflux_out%lun,&
!!$                  ' to be linked to file ',etb(this%pflux_out%fn)
!!$             stop ' simulation ended'
!!$          else
!!$             this%pflux_out%exist = .true.
!!$          end if
!!$
!!$          ! pmass          
!!$          open(this%pmass_out%lun,file=this%pmass_out%fn,iostat=opstat)
!!$          if(opstat.ne.0) then
!!$             write(lun,*) ' failed to open unit ',this%pmass_out%lun,&
!!$                  ' to be linked to file ',etb(this%pmass_out%fn)
!!$             stop ' simulation ended'
!!$          else
!!$             this%pmass_out%exist = .true.
!!$          end if
!!$
!!$          ! decay          
!!$          open(this%decay_out%lun,file=this%decay_out%fn,iostat=opstat)
!!$          if(opstat.ne.0) then
!!$             write(lun,*) ' failed to open unit ',this%decay_out%lun,&
!!$                  ' to be linked to file ',etb(this%decay_out%fn)
!!$             stop ' simulation ended'
!!$          else
!!$             this%decay_out%exist = .true.
!!$          end if
!!$
!!$          ! kappa          
!!$          open(this%kappa_out%lun,file=this%kappa_out%fn,iostat=opstat)
!!$          if(opstat.ne.0) then
!!$             write(lun,*) ' failed to open unit ',this%kappa_out%lun,&
!!$                  ' to be linked to file ',etb(this%kappa_out%fn)
!!$             stop ' simulation ended'
!!$          else
!!$             this%kappa_out%exist = .true.
!!$          end if
!!$
!!$          ! tdens0          
!!$          open(this%tdens0_out%lun,file=this%tdens0_out%fn,iostat=opstat)
!!$          if(opstat.ne.0) then
!!$             write(lun,*) ' failed to open unit ',this%tdens0_out%lun,&
!!$                  ' to be linked to file ',etb(this%tdens0_out%fn)
!!$             stop ' simulation ended'
!!$          else
!!$             this%tdens0_out%exist = .true.
!!$          end if
!!$
!!$          ! rhs_integrated
!!$          inquire(file=this%rhs_integrated_out%fn,exist=file_exist)
!!$          if ( file_exist ) then 
!!$             this%rhs_integrated_out%exist = .true.
!!$             open(this%rhs_integrated_out%lun,file=this%rhs_integrated_out%fn,iostat=opstat)
!!$             if(opstat.ne.0) then
!!$                write(lun,*) ' failed to open unit ',this%rhs_integrated_out%lun,&
!!$                     ' to be linked to file ',etb(this%rhs_integrated_out%fn)
!!$                stop ' simulation ended'
!!$             end if
!!$          else
!!$             this%rhs_integrated_out%exist = .true.   
!!$          end if
!!$
!!$
!!$          ! optdens
!!$          open(this%optdens_out%lun,file=this%optdens_out%fn,iostat=opstat)
!!$          if(opstat.ne.0) then
!!$             write(lun,*) ' failed to open unit ',this%optdens_out%lun,&
!!$                  ' to be linked to file ',etb(this%optdens_out%fn)
!!$             write(lun,*) 'Optimal density does not exist' 
!!$          else
!!$             this%optdens_out%exist = .true.   
!!$          end if
!!$          
!!$          ! forcing          
!!$          open(this%forcing_out%lun,file=this%forcing_out%fn,iostat=opstat)
!!$          if(opstat.ne.0) then
!!$             write(lun,*) ' failed to open unit ',this%forcing_out%lun,&
!!$                  ' to be linked to file ',etb(this%forcing_out%fn)
!!$             stop ' simulation ended'
!!$          else
!!$             this%forcing_out%exist = .true.
!!$          end if
!!$
!!$          ! boundary_flux        
!!$          open(this%boundary_flux_out%lun,file=this%boundary_flux_out%fn,iostat=opstat)
!!$          if(opstat.ne.0) then
!!$             write(lun,*) ' failed to open unit ',this%boundary_flux_out%lun,&
!!$                  ' to be linked to file ',etb(this%boundary_flux_out%fn)
!!$             stop ' simulation ended'
!!$          else
!!$             this%boundary_flux_out%exist = .true.
!!$          end if
!!$
!!$          ! debug
!!$          open(this%debug%lun,file=this%debug%fn,iostat=opstat)
!!$          this%debug%exist = .true.
!!$          if(opstat.ne.0) then
!!$             write(lun,*) ' failed to open unit ',this%debug%lun,&
!!$                  ' to be linked to file ',etb(this%debug%fn)
!!$             stop ' simulation ended'
!!$          end if
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

      ! rhs_integrated
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
    call this%control%info(lun)
    call this%grid%info(lun)
    call this%initial%info(lun)
    !
    call this%pflux%info(lun)
    call this%pmass%info(lun)
    call this%decay%info(lun)
    call this%kappa%info(lun)
    call this%tdens0%info(lun)
    call this%optdens%info(lun)
    call this%rhs_integrated%info(lun)
    call this%forcing%info(lun)
    call this%boundary_flux%info(lun)
 
    write(lun,*) ' Output units'
    call this%output_folder%info(lun)
    call this%grid_out%info(lun)
    call this%subgrid%info(lun)
    call this%sons%info(lun)
    !
    call this%pflux_out%info(lun)
    call this%pmass_out%info(lun)
    call this%decay_out%info(lun)
    call this%kappa_out%info(lun)
    call this%tdens0_out%info(lun)
    call this%optdens_out%info(lun)
    call this%rhs_integrated_out%info(lun)
    call this%forcing_out%info(lun)
    call this%boundary_flux_out%info(lun)

    call this%debug%info(lun)

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

    close (this%input_folder%lun)
    !
    close (this%control%lun)
    close (this%grid%lun)
    close (this%initial%lun)
    !
    close (this%pflux%lun)
    close (this%pmass%lun)
    close (this%decay%lun)
    close (this%kappa%lun)
    close (this%tdens0%lun)
    close (this%optdens%lun)
    close (this%rhs_integrated%lun)
    close (this%forcing%lun)
    close (this%boundary_flux%lun)
 
    close (this%output_folder%lun)
    close (this%grid_out%lun)
    close (this%subgrid%lun)
    close (this%sons%lun)
    !
    close (this%pflux_out%lun)
    close (this%pmass_out%lun)
    close (this%decay_out%lun)
    close (this%kappa_out%lun)
    close (this%tdens0_out%lun)
    close (this%optdens_out%lun)
    close (this%rhs_integrated_out%lun)
    close (this%forcing_out%lun)
    close (this%boundary_flux_out%lun)

    close (this%debug%lun)

    
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
