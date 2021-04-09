module SparseMatrix
  use Globals
  use Matrix
  implicit none
  private
  public :: MtimesN
  !>----------------------------------------------------------------------
  !> Structure variable containg member storing sparse real matrices
  !> in csr ( compress sparse row ) and ssr (symmetric sparse row format )
  !>----------------------------------------------------------------------
  type, extends(abs_matrix) ,public :: spmat
     !> Flag if spmat has been initialized
     logical :: is_initialized =.false.
     !> Flag for ordered coefficient
     logical :: is_sorted=.true.
     !> Flag for storage style
     !> 'csr'= Compressed Storage Row
     !> 'ssr'= Symmetric  Storege Row (use only
     !> for Symmetric Matrices, it stores only the
     !> upper triangular part)
     character(len=3) :: storage_system
     !> Number of non zero terms
     integer :: nterm=0
     ! Dimension (nrow+1)
     ! Pointers to the first term of sysamt in row
     integer, allocatable :: ia(:)
     ! Dimension (nrow)
     ! Pointers to the diagonal term in coeff
     ! May be not always be allocated
     integer, allocatable :: ind_diag(:)
     ! Dimension (nterm)
     ! Pointers to the column index of coeff
     integer, allocatable :: ja(:)
     ! Dimension (nterm)
     ! Non zero elements of the sparse matrix
     real(kind=double), allocatable :: coeff(:)
   contains
     !> static constructor 
     !> (procedure public for type spmat)
     procedure, public, pass :: init => init_spmat
     !> static destructor
     !> (procedure public for type spmat)
     procedure, public, pass :: kill => kill_spmat
     !> Reading procedure.
     !> (public procedure for type spmat)
     procedure, public, pass :: read => read_spmat
     !> Writing procedure.
     !> (public procedure for type spmat)
     procedure, public, pass :: write => write_spmat
     !> Info procedure  sparse matrix
     !> (public procedure for type spmat)
     procedure, public, pass :: info => info_spmat
     !> Build the pointers to the diagonal terms,
     !> for a specific row, given ia and ja
     !> (public procedure for type spmat)
     procedure, public, pass :: idiag
     !> Build the pointers to the diagonal terms,
     !> given ia and ja
     !> (public procedure for type spmat)
     procedure, public, pass :: build_ind_diag
     !> Build the Diagonal Matrix array
     !> (public procedure for type spmat)
     procedure, public, pass :: get_diagonal
     !> Convert CSR mat. into SSR
     !> It can initilize a new matrix
     !> (public procedure for type spmat)
     procedure, public , pass :: csr2ssr
     !> Convert SSR mat. into CSR
     !> It can initilize a new matrix
     !> (public procedure for type spmat)
     procedure, public , pass :: ssr2csr
!!$     !> Permute a matrix
!!$     !> (public procedure for type spmat)
!!$     procedure, public , pass :: perm => perm_mat
     !> Procedure to fix the irow-term of the solution to 
     !> given value, cmodifieng the matrix and the rhs
     !> (public procedure for type spmat)
     procedure, public , pass :: fix_irow
     !> Procedure to operate the diagonal scaling
     !> of matrix M by the diagonal matrix D
     !> out_matrix = left_D M right_D
     !> (public procedure for type spmat)
     procedure, public , pass :: diagonal_scale
     !> Procedure to  multiply sparse matrix M 
     !> by the diagonal matrix D 
     !> out_matrix = D M D
     !> (public procedure for type spmat)
     procedure, public , pass :: MxD
     !> Procedure to multiply  diagonal matrix D 
     !> by the sparse matrix M 
     !> out_matrix = D M
     !> (public procedure for type spmat)
     procedure, public , pass :: DxM
     !> Procedure to compute 
     !>         y = M * x 
     !> (public procedure for type spmat)
     procedure, public, pass :: Mxv
     !> Procedure to compute 
     !>         y = M^T * x 
     !> with M^T the transposed of a matrix M
     !> (public procedure for type spmat)
     procedure, public, pass :: MTxv
     !> Logical output for dimensions matching
     !> (public procedure for type spmat)
     procedure, public, pass :: check => check_spmat
     !> Procedure to sort arrays ja and coeff
     !> so that the colums indeces of ja of each row
     !> (j = ja(ia(irow),ja(ia(irow+1)-1) 
     !> are in increasing order 
     !> (public procedure for type spmat)
     procedure, public, pass :: sort => sort_spmat
     !> Procedure to sort arrays the passed portion
     !> ja and coeff so that the colums indeces of ja
     !> are sorted in increasing order 
     !> (public procedure for type spmat)
     procedure, private, nopass :: sort_row
     !> Procedure for trasposition (operate in place)
     !> (public procedure for type spmat)
     procedure, public, pass :: transpose => transpose_spmat 
     !> Procedure to solve system
     !>           M x = b
     !> with M triagular
     !> (public procedure for type spmat)
     procedure, public, pass :: solve_triangular
     !> Procedure to solve system
     !>           M x = b
     !> with M triagular
     !> (public procedure for type spmat)
     procedure, public, pass :: solve_triangular_unitdiag
     !> Procedure Build the approximate Cholesky factorization 
     !> of a symmetric sparse matrix M in the form
     !>          M = U^T U
     !> (public procedure for type spmat)
     procedure, public, pass :: incomplete_cholesky
     !> Procedure Build the approximate factorization 
     !> of a sparse matrix M in the form
     !>          M = L U
     !> with L and U lower and upper triangular matrix
     !> (public procedure for type spmat)
     procedure, public, pass :: incomplete_lu
     !> Procedure Build the approximate factorization 
     !> of a sparse matrix M in the form
     !>          M = L U
     !> with L and U lower and upper triangular matrix
     !> (public procedure for type spmat)
     procedure, public, nopass :: MtimesN
     !> Procedure that adds a row to a matrix
     !> ( A | 0 )
     !> ( w | 1 )
     !> where A is the matrix, w is a vector (dim=nrow)
     procedure, public, pass :: add_row
     !> Procedure Build the reverse Cuthill-Mckee ordering 
     !> for the spmat contiang the connection of a
     !> a genreal graph.
     !> (public procedure for type spmat)
     procedure, public, pass :: genrcm
     !> Procedure to compute the maximum bandwidth 
     !> for a general matrix.
     !> (public procedure for type spmat)
     procedure, public, pass :: bandwidth
     !> Procedure to create ps file containg the
     !> sparsity pattern of a sparse matrix.
     !> (public procedure for type spmat)
     procedure, public, pass :: plot2ps
  end type spmat
contains
  !>-------------------------------------------------------------
  !> Static constructor.
  !> (procedure public for type spmat)
  !> Instantiate (allocate if necessary)
  !> and initilize (by also reading from input file)
  !> variable of type spmat
  !>
  !> usage:
  !>     call 'var'%init(lun_err, nrow, nterm )
  !>
  !> where:
  !> \param[in] lun_err               -> integer. Error logical unit
  !> \param[in] nrow                  -> integer. Number of rows 
  !> \param[in] nrow                  -> integer. Number of columns
  !> \param[in] nterm                 -> integer. Number of non-zero term
  !>                                     stored in the matrix
  !> \param[in] (optional) dim_kernel -> integer. Number of non-zero term
  !> \param[in] (optional) is_sym     -> Logical. T/F flag for symmetric matrix
  !> \param[in] (optional) triangular -> Character. 
  !>                                       'N' = not triangular ( the default)
  !>                                       'U' = uppertriangular
  !>                                       'L' = lower_triangular' 
  !<-------------------------------------------------------------
  subroutine init_spmat(this, lun_err, &
       nrow, ncol, nterm,&
       storage_system,&
                                ! optional arguments
       dim_kernel, &
       is_symmetric, &
       triangular,unitary_diag)
    use Globals
    implicit none
    !var
    class(spmat),                intent(inout) :: this
    integer,                     intent(in   ) :: lun_err
    integer,                     intent(in   ) :: nrow
    integer,                     intent(in   ) :: ncol
    integer,                     intent(in   ) :: nterm
    character(len=*),            intent(in   ) :: storage_system
    integer,           optional, intent(in   ) :: dim_kernel
    logical,           optional, intent(in   ) :: is_symmetric
    character (len=1), optional, intent(in   ) :: triangular
    logical,           optional, intent(in   ) :: unitary_diag

    ! local vars
    integer :: res
    logical :: rc


    this%is_initialized  = .true.
    this%nrow            = nrow
    this%ncol            = ncol
    this%nterm           = nterm

    if ( this%ncol .eq. this%nrow) this%is_squared=.true.


    allocate(this%ia(nrow+1),stat=res)
    if(res .ne. 0) rc = IOerr(lun_err, err_alloc, 'init_spmat', &
         '  type spmat member ia (array)',res)

    allocate(this%ja(nterm),stat=res)
    if(res .ne. 0) rc = IOerr(lun_err, err_alloc, 'init_spmat', &
         '  type spmat member ja (array)',res)

    allocate(this%coeff(nterm),stat=res)
    if(res .ne. 0) rc = IOerr(lun_err, err_alloc, 'init_spmat', &
         '  type spmat member coeff (array)',res)


    if (present(dim_kernel) ) then
       this%dim_kernel = dim_kernel
       allocate(this%kernel(ncol,dim_kernel),stat=res)
       if(res .ne. 0) rc = IOerr(lun_err, err_alloc, 'init_spmat', &
            '  type spmat member kernel',res)
    end if

    this%storage_system = storage_system

    if (this%storage_system .eq. 'ssr') then
       this%is_symmetric=.true.
       this%triangular  = 'N'
       if ( present(is_symmetric) ) then
          if ( (.not.is_symmetric) .and. this%storage_system.eq.'ssr') &
               rc = IOerr(lun_err, err_inp, 'init_spmat', &
               '  conflict between storage_system and is_symmetric)')
          this%is_symmetric   = is_symmetric
       end if

    else
       
       if ( present(is_symmetric) ) then
          this%is_symmetric   = is_symmetric
       end if

       if ( present(triangular) ) then
          if ( ( triangular .eq. 'N') .or. &
               ( triangular .eq. 'U') .or. &
               ( triangular .eq. 'L') ) then
             this%triangular = triangular
          else
             rc = IOerr(lun_err, err_inp, 'init_spmat', &
                  '  wrong triangular passed ='//etb(triangular))
          end if
       end if
    end if

    if ( present(unitary_diag) ) this%unitary_diag = unitary_diag

  end subroutine init_spmat

  !>-------------------------------------------------------------
  !> Static destructor.
  !> (procedure public for type spmat)
  !> deallocate all arrays for a var of type spmat
  !>
  !> usage:
  !>     call 'var'%kill(lun)
  !>
  !> where:
  !> \param[in] lun -> integer. I/O unit for error message output
  !<-----------------------------------------------------------
  subroutine kill_spmat(this, lun_err)
    implicit none
    ! vars
    class(spmat),intent(inout) :: this
    integer,     intent(in   ) :: lun_err
    ! local vars
    integer :: res
    logical :: rc


    deallocate(this%ia,stat=res)
    if (res.ne.0) rc=IOerr(lun_err, err_dealloc, 'kill_spmat', &
         'dealloc fail for spmat members ia',res)
    deallocate(this%ja,stat=res)
    if (res.ne.0) rc=IOerr(lun_err, err_dealloc, 'kill_spmat', &
         'dealloc fail for spmat members ja',res)

    deallocate(this%coeff,stat=res)
    if (res.ne.0) rc=IOerr(lun_err, err_dealloc, 'kill_spmat', &
         'dealloc fail for spmat member coeff',res)


    if ( allocated(this%ind_diag) ) then
       deallocate(this%ind_diag, stat=res)
       if(res .ne. 0) rc = IOerr(lun_err, err_dealloc, 'kill_spmat', &
            '  type spmat member ia (array)',res)

    end if


    if ( allocated(this%kernel) ) then
       deallocate(this%kernel,stat=res)
       if (res.ne.0) rc=IOerr(lun_err, err_dealloc, 'kill_spmat', &
            'dealloc fail for spmat member kernel')
    end if

    this%dim_kernel  = 0
    this%nrow        = 0
    this%nterm       = 0
    this%triangular     = 'N'
    this%is_symmetric   = .false.
    this%is_squared     = .false.
    this%is_initialized =.false.

  end subroutine kill_spmat

  !>-------------------------------------------------------------
  !> Reading procedure.
  !> (public procedure for type spmat)
  !> Read content of a variable of type spamat
  !> 
  !> usage:
  !>     call 'var'%read(lun)
  !>
  !> where:
  !> \param[in] lun -> integer. I/O unit for error message output
  !>
  !<-------------------------------------------------------------
  subroutine read_spmat(this,lun_err,input_file)
    use Globals
    implicit none
    class(spmat),      intent(inout) :: this
    integer,           intent(in   ) :: lun_err
    type(file),        intent(in   ) :: input_file
    ! loc. var
    logical :: rc
    integer :: res, lun
    integer:: iterm,icol,irow,irow_tmp

    lun = input_file%lun
    
    read(lun,*) this%nrow, this%ncol! '! number of rows and colums '
    read(lun,*) this%nterm !'! number of non-zero terms '
    allocate(this%ia(this%nrow+1),this%ja(this%nterm),&
         this%coeff(this%nterm), stat=res)
    if(res .ne. 0) rc = IOerr(lun_err, err_alloc, 'read_spmat', &
         '  type spmat member ia, ja, coeff (array)',res)
    irow = 1
    this%ia(1) = 1
    do iterm=1,this%nterm
       read(lun,*) irow_tmp,this%ja(iterm),this%coeff(iterm)
       if(irow_tmp.gt.irow) then
          irow = irow+1
          this%ia(irow) = iterm
       end if
    end do
    this%ia(this%nrow+1) = this%nterm+1

    this%is_initialized = .true.
    
  end subroutine read_spmat

  !>-------------------------------------------------------------
  !> Writing procedure.
  !> (public procedure for type spmat)
  !> Prints content of a variable of type spamat
  !> 
  !> usage:
  !>     call 'var'%write(lun)
  !>
  !> where:
  !> \param[in] lun -> integer. I/O unit for error message output
  !>
  !<-------------------------------------------------------------
  subroutine write_spmat(this, lun)
    use Globals
    implicit none
    class(spmat),      intent(in) :: this
    integer,           intent(in) :: lun
    !real(kind=double), optional, intent(in) :: rho
    ! loc. var
    integer i,m,n,ind

    write(lun,*) this%nrow, '! number of rows     '
    write(lun,*) this%nterm,'! number of non-zero terms '
    do i=1,this%nrow
       m=this%ia(i)
       n=this%ia(i+1)-1
       write(lun,1010) (i,this%ja(ind),this%coeff(ind),ind=m,n)
    end do
1010 format(2i15,1pe24.16)

  end subroutine write_spmat

  !>-------------------------------------------------------------
  !> Info procedure.
  !> (public procedure for type spmat)
  !> Prints content of a variable of type spamat
  !> 
  !> usage:
  !>     call 'var'%info(lun_out)
  !>
  !> where:
  !> \param[in] lun_out -> integer. I/O unit for output message
  !>
  !<-------------------------------------------------------------
  subroutine info_spmat(this, lun)
    use Globals
    implicit none
    class(spmat), intent(in) :: this
    integer :: lun
    ! loc. var
    character(len=256) :: state,mem,prop
    character(len=3) :: sep=' | '

    if ( this%is_initialized ) then
       write(state,'(a)') 'Spmat allocated'
       mem=this%storage_system
       if ( this%triangular .ne. 'N') then
          if ( this%triangular .eq. 'U' ) then
             write(prop,'(a)') 'upper-triangular'
          else 
             write(prop,'(a)') 'lower-triangular'
          end if
       else
          if ( this%is_symmetric ) then  
             write(prop,'(a)') 'symmetric'
          else
             write(prop,'(a)') ''
          end if
       end if

       write(lun,'(a)') etb(etb(state)//sep//etb(mem)//sep//etb(prop))

       write(lun,'(a,I8,a,a,I8,a,a,I8)') 'nrows= ', this%nrow,sep,&      
            'ncol= ', this%ncol,sep,&      
            'non-zero terms= ', this%nterm     
       if (allocated(this%kernel)) then
          write(lun,'(a,I8)') 'kernel-dim= ', this%dim_kernel
       end if
    else
       write(lun,*) 'Spmat not initialized'
    end if

  end subroutine info_spmat


  !>-------------------------------------------------------------
  !> check procedure.
  !> (public procedure for type spmat)
  !> verifies if a matrix structure sastisfy the array bounds
  !> for given, nrow, ncol, nterm
  !> 
  !> usage:
  !>     sym = 'var'%check(nrow,ncol,nterm)
  !>
  !> \param[in] (optional) nrow          -> integer. Nmb. of rows required
  !> \param[in] (optional) ncol          -> integer. Nmb. of columnss required
  !> \param[in] (optional) nterm         -> integer. Nmb. of non-zero 
  !>                                       terms required
  !> \param[in] (optional)            -> integer. Storage system
  !> \param[in] (optional) is_symmetric  -> Logical. pass ".true." to  
  !>                                        ask if spmat is symmetric
  !> \param[in] (optional) triangular         -> Character(len=1). pass "U" or "L"
  !>                                        to ask if spmat is upper or lower 
  !>                                        triangular
  !> \return    check_spmat -> logical. true/false flag if Spmat respects
  !>                                   dimension bound   e               
  !<------------------------------------------------------------
  function check_spmat(this, &
       nrow, ncol, nterm,&
       storage_system,dim_kernel,is_symmetric,triangular)
    use Globals
    implicit none
    class(spmat),                intent(in) :: this
    integer,           optional, intent(in) :: nrow
    integer,           optional, intent(in) :: ncol
    integer,           optional, intent(in) :: nterm
    character(len=3),  optional, intent(in) :: storage_system
    integer,           optional, intent(in) :: dim_kernel
    logical          , optional, intent(in) :: is_symmetric
    character (len=1), optional, intent(in) :: triangular


    logical                       :: check_spmat

    check_spmat = .true.

    if ( present(nrow) .and. (nrow .ne. this%nrow ) ) then
       check_spmat=.false. 
    end if

    if ( present(ncol) .and. (ncol .ne. this%ncol ) ) then
       check_spmat=.false. 
    end if

    if ( present(nterm) .and. (nterm .ne. this%nterm) ) then
       check_spmat = .false.
    end if


    if ( present(storage_system) .and. &
         (storage_system .ne. this%storage_system) ) then
       check_spmat = .false.
    end if

    if ( present(dim_kernel) .and. (dim_kernel .ne. this%dim_kernel) ) then
       check_spmat = .false.
    end if


    if ( ( present(is_symmetric)  ) .and. &
         ( .not. this%is_symmetric ) ) then
       check_spmat = .false.
    end if

    if ( ( present (triangular)  ) .and. &
         ( this%triangular .ne. triangular  ) ) then
       check_spmat = .false.
    end if




  end function check_spmat


  !>-------------------------------------------------------------
  !> Procedure for identify the pointer the diagonal 
  !> term(irow) of a matrix 
  !>
  !> usage : this%idiag(irow)
  !> 
  !> where:
  !> \param[in] irow  -> integer. Index of the row
  !> \return    idiag -> integer. Index of the diagonal term of the
  !>                              irow^th row
  !>                              idiag=0 if the diagonal term is null
  !>--------------------------------------------------------------
  function idiag( this, irow )
    use Globals
    implicit none
    class(spmat), intent(in) :: this
    integer,      intent(in) :: irow
    integer :: idiag
    !local 
    integer :: j

    if (allocated(this%ind_diag) ) then
       idiag = this%ind_diag(irow)
    else
       if ( this%storage_system .eq. 'ssr' ) then
          !
          ! use ssr structure
          !
          idiag = this%ia(irow)

          if ( this%ja(idiag) .ne. irow ) idiag = 0 
       else if ( this%storage_system == 'csr' ) then 
          if ( this%triangular .ne. 'N' ) then
             ! use triangular structure
             select case( this%triangular ) 
             case ('L') 
                idiag = this%ia(irow+1)-1
             case ('U') 
                idiag = this%ia(irow)
             end select
             if ( this%ja(idiag) .ne. irow ) idiag = 0
          else
             ! search in row
             j = this%ia(irow)
             idiag = 0
             do j = this%ia(irow),this%ia(irow+1) - 1
                if (this%ja(j) .eq. irow) then
                   idiag=j
                   exit
                end if
             end do
          end if
       end if
    end if
  end function idiag


  !>-------------------------------------------------------------
  !> Procedure to initilized and build the member ind_diag
  !> of structure spamt     
  !>
  !> usage : this%build_diag(lun_err)
  !>
  !> \param[in] lun_err -> integer. I/O unit for err. msg. 
  !>--------------------------------------------------------------
  subroutine build_ind_diag(this,lun_err)
    use Globals
    implicit none
    class(spmat),      intent(inout) :: this
    integer,           intent(in   ) :: lun_err
    !local 
    logical rc
    integer res
    integer i

    if ( .not. allocated(this%ind_diag)) then
       allocate (this%ind_diag(this%nrow),stat=res) 
       if(res .ne. 0) rc = IOerr(lun_err, err_alloc, &
            'build_ind_diag', &
            'type spmat member ind_diag',res)
    end if

    do i = 1,this%nrow
       this%ind_diag(i) = this%idiag(i)
    end do

  end subroutine build_ind_diag

  !>-------------------------------------------------------------
  !> Procedure to get the diagonal of a structure varible spmat
  !> of a matrix in ssr format
  !>
  !> usage : this%get_diagonal(diag)
  !> 
  !> \param[out] diag -> real. (dimension = this%nrow)
  !>                       Diagonal of spmat
  !>--------------------------------------------------------------
  subroutine get_diagonal(this,diagonal)
    use Globals
    implicit none
    class(spmat),      intent(in   ) :: this
    real(kind=double), intent(inout) :: diagonal(min(this%nrow,this%ncol))
    !local 
    integer i, ind
    
    diagonal=zero   
    do i=1,min(this%nrow,this%ncol)
       ind = this%idiag(i)
       if (ind .ne. 0 ) diagonal(i) = this%coeff(ind)
    end do

  end subroutine get_diagonal


  !>------------------------------------------------------------------
  !> Procedure to transform matrix structure from ssr into csr format.
  !> Reinitialized the matrix.
  !>
  !> usage: call var%ssr2csr(lun_err)
  !>
  !> where:
  !> \param[in ] lun_err -> integer. I\O unit for error message
  !>------------------------------------------------------------------   
  subroutine ssr2csr(this,lun_err)
    use Globals
    implicit none
    class(spmat), intent(inout) :: this
    integer,      intent(in   ) :: lun_err

    !local
    integer i, j ,m
    integer nrow, nterm, nterm_out
    type(spmat) :: upper, lower

    nrow =this%nrow
    nterm=this%nterm
    nterm_out=2*nterm-nrow

    ! create a copy and change to csr 
    select type (this)
    type is (spmat) 
       upper = this
    end select
    upper%storage_system = 'csr'
    upper%triangular     = 'U'

    ! create workin matrix
    ! We used the following trick: change the storage system
    ! to upper to build a lower triangular matrix scr_transpose
    lower =  upper
    call lower%transpose(lun_err)

    ! renitialized 
    call this%kill(lun_err)
    call this%init(lun_err, &
         nrow, nrow, nterm_out,&
         'csr',&
         upper%dim_kernel, is_symmetric=.true.)
    this%kernel=upper%kernel

    ! the matrix is build this=upper+lower-diag
    m=0
    this%ia(1)=1
    do i = 1,nrow
       do j = lower%ia(i),lower%ia(i+1)-2              
          m=m+1
          this%ja(m) = lower%ja(j)
          this%coeff(m) = lower%coeff(j)              
       end do
       do j = upper%ia(i),upper%ia(i+1)-1
          m=m+1
          this%ja(m) = upper%ja(j)
          this%coeff(m) = upper%coeff(j)
       end do
       this%ia(i+1) = m+1
    end do


    ! free memory 
    call upper%kill(lun_err)
    call lower%kill(lun_err)


  end subroutine ssr2csr

  !>---------------------------------------------------------------
  !> Procedure to transform a symmetric matrix from csr into ssr
  !> format. Initialize the varaible mat_out if it is not
  !> 
  !> usage: call var%ssr2csr(lun_err,mat_out)
  !>
  !!> where:
  !> \param[in ] lun_err            -> integer. I\O unit for error message
  !> \param[out] (optional) mat_out -> type(spmat). Matrix in csr format
  !<---------------------------------------------------------------------
  subroutine csr2ssr(this, lun_err)
    use Globals
    implicit none
    class(spmat), intent(inout) :: this
    integer,      intent(in   ) :: lun_err

    !local
    integer :: i, j , m
    integer :: nrow, ncol, nterm_out
    type(spmat) :: ssr        

    if ( .not. this%is_symmetric  ) then
       write(lun_err,*) ' Not symmetric matrix'
       write(lun_err,*) ' Convertion not possible'
       stop
    end if

    nrow      = this%nrow
    ncol      = this%ncol
    nterm_out = (this%nterm + nrow) / 2

    ! initialized and set properteis of working spmat
    call ssr%init(lun_err,nrow,ncol, nterm_out, &
         'ssr',this%dim_kernel)
    ssr%kernel = this%kernel

    m=0
    ssr%ia(1)=1
    do i = 1, nrow       
       do j = this%idiag(i), this%ia(i+1)-1
          m=m+1
          ssr%ja(m)    = this%ja(j)
          ssr%coeff(m) = this%coeff(j)
       end do
       ssr%ia(i+1) = m + 1
    end do

    ! copy the working spmat into this
    select type (this)
    type is (spmat) 
       this = ssr
    end select

    ! free memory
    call ssr%kill(lun_err)

  end subroutine csr2ssr

!!$    !>-----------------------------------------------------------
!!$    !> Procedure to build the matrix PAP' for given sparse matrix 
!!$    !> A and a permutation perm. It initializes spmat mat_out.
!!$    !>
!!$    !> usage : call var%perm(lun_err,nrow,perm,iperm)
!!$    !>
!!$    !> where:
!!$    !> \param[in ] lun_err -> integer. I\O unit for error message
!!$    !> \param[in ] nrow    -> integer. Nmb of rows and columns
!!$    !> \param[in ] perm    -> integer(nrow). Permutation
!!$    !> \param[in ] iperm   -> integer(nrow). Inverse Permutation
!!$    !>-----------------------------------------------------------
!!$    subroutine perm_mat(this, lun_err, nrow, perm , inv_perm )
!!$      use Globals
!!$      implicit none
!!$      class(spmat), intent(inout)  :: this
!!$      integer,      intent(in   )  :: lun_err, nrow
!!$      integer,      intent(in   )  :: perm(nrow)
!!$      integer,      intent(in   )  :: inv_perm(nrow)
!!$
!!$      !local
!!$      logical :: rc
!!$      integer :: res,info
!!$      integer :: i
!!$      type(spmat) :: full, sorted
!!$
!!$      if ( this%storage_system .eq. ssr ) then
!!$         ! local copy in ssr format
!!$         full = this
!!$         call full%ssr2csr(lun_err)
!!$
!!$         ! create a sorted the local copy
!!$         call sorted%init(lun_err,&
!!$              full%nrow,full%nterm,full%dim_kernel)
!!$         
!!$         call apply_perm(lun_err,&
!!$              full%nrow,full%nterm,&
!!$              perm,inv_perm,&
!!$              full%ia,full%ja,full%coeff,&
!!$              sorted%ia,sorted%ja,sorted%coeff,&
!!$              info)
!!$
!!$         do i=1,this%dim_kernel 
!!$            call double_permute(
!!$            this%kernel(:,i) = this%kernel  ( inv_perm, i)
!!$         end do
!!$         
!!$         call full%kill(lun_err)
!!$         ! convert sorted mat into ssr format and copy into this
!!$         call sorted%csr2ssr(lun_err) 
!!$         this = sorted
!!$         call sorted%kill(lun_err)
!!$
!!$         
!!$
!!$      else         
!!$         call sorted%init(&
!!$              lun_err,1,nrow,this%nterm,this%nterm,this%dim_kernel)
!!$         call apply_perm(lun_err,&
!!$              this%nrow,this%nterm,&
!!$              perm,iperm,&
!!$              this%ia,this%ja,this%coeff,&
!!$              sorted%ia,sorted%ja,sorted%coeff,&
!!$              info)                  
!!$         this = sorted
!!$         call sorted%kill(lun_err)
!!$         do i=1,this%dim_kernel 
!!$            this%kernel(:,i) = this%kernel  ( iperm, i)
!!$         end do
!!$      end if
!!$      if (info .ne. 0) then
!!$         write(lun_err,*) 'Error in procedure perm_mat'
!!$         stop
!!$      end if
!!$      
!!$    contains
!!$      subroutine apply_perm(lun_err,&
!!$           nn,nt,&
!!$           perm,iperm,&
!!$           iat_in,ja_in,coef_in,&
!!$           iat_out,ja_out,coef_out,&
!!$           info)
!!$        !--------------------------------------
!!$        !
!!$        ! Creates a permuted matrix
!!$        !
!!$        !--------------------------------------
!!$
!!$        use Globals
!!$
!!$        implicit none
!!$
!!$        ! Interfaces
!!$        integer                        :: lun_err
!!$        integer, intent(in)            :: nn,nt
!!$        integer, intent(in)            :: perm(nn),iperm(nn)
!!$        integer, intent(in)            :: iat_in(nn+1),ja_in(nt)
!!$        real(kind=double), intent(in)  :: coef_in(nt)
!!$
!!$        integer, intent(out)           :: info
!!$        integer, intent(out)           :: iat_out(nn+1),ja_out(nt)
!!$        real(kind=double), intent(out) :: coef_out(nt)
!!$
!!$        integer                        :: i,j,k,m1,m2
!!$        integer, allocatable           :: irow_tmp(:)
!!$        !
!!$        integer, allocatable           :: iat_tmp(:),ja_tmp(:)
!!$        real(kind=double), allocatable :: coef_tmp(:)
!!$        !
!!$        !local
!!$        logical :: rc
!!$        integer :: res
!!$        info = 0
!!$
!!$        allocate(irow_tmp(nt),iat_tmp(nn+1),ja_tmp(nt),coef_tmp(nt),stat=res)
!!$        if(res .ne. 0) rc = IOerr(lun_err, err_alloc, 'init_spmat', &
!!$             '  type spmat member ia (array)',res)
!!$
!!$        k = 0
!!$        do i = 1,nn
!!$           m1 = iat_in(iperm(i))
!!$           m2 = iat_in(iperm(i)+1)-1
!!$           do j = m1,m2
!!$              k = k + 1
!!$              irow_tmp(k) = i
!!$              ja_out(k) = perm(ja_in(j))
!!$              coef_out(k) = coef_in(j)
!!$           end do
!!$        end do
!!$        call irow2iat(nn,nt,irow_tmp,iat_out)
!!$        call TranspMat(nn,nt,&
!!$             iat_out,ja_out,coef_out,&
!!$             nn,k,iat_tmp,ja_tmp,coef_tmp,&
!!$             irow_tmp)
!!$        call TranspMat(nn,nt,&
!!$             iat_tmp,ja_tmp,coef_tmp,&
!!$             nn,k,iat_out,ja_out,coef_out,&
!!$             irow_tmp)
!!$
!!$        deallocate(irow_tmp,iat_tmp,ja_tmp,coef_tmp,stat=info)
!!$        if (info .ne. 0) return
!!$
!!$      end subroutine apply_perm
!!$      subroutine irow2iat(n,nterm,irow,iat)
!!$        !-----------------------------------------------------------------
!!$        !
!!$        !  Subroutine: irow2iat
!!$        !
!!$        !  Coded by Carlo Janna
!!$        !  April 2012
!!$        !
!!$        !  Purpose: Build the topology vector iat 
!!$        !  from a row indices list irow
!!$        !
!!$        !  Variables:
!!$        !
!!$        !  n     : # of rows of a given matrix stored in CSR format
!!$        !  nterm : # of non-zeroes of a given matrix stored in CSR format
!!$        !  irow  : row indices of non-zeroes of a matrix in CSR format
!!$        !  iat   : integer array of the pointers to the begin of each row
!!$        !
!!$        !----------------------------------------------------------------
!!$
!!$        implicit none
!!$
!!$        !Input variables
!!$        integer, intent(in) :: n,nterm
!!$        integer, intent(in) :: irow(nterm)
!!$
!!$        !Output variables
!!$        integer, intent(out) :: iat(n+1)
!!$
!!$        !Local variables
!!$        integer             :: j,k,irow_old,irow_new
!!$
!!$        !------------------------------------------------------------
!!$
!!$        irow_old = 0
!!$        do k=1,nterm
!!$           irow_new = irow(k)
!!$           if ( irow_new .gt. irow_old ) then
!!$              do j = irow_old+1,irow_new
!!$                 iat(j)=k
!!$              enddo
!!$              irow_old = irow_new
!!$           end if
!!$        end do
!!$        k = nterm+1
!!$        do j = irow_old+1,n+1
!!$           iat(j) = k
!!$        enddo
!!$
!!$      end subroutine irow2iat
!!$
!!$    end subroutine perm_mat
!!$    

  !>-------------------------------------------------------------
  !> Procedure for reset a row and optionally compute
  !> the correction of rhs in order to fix the solution
  !> at irow with value sol_irow
  !> (public procedure for type spmat, works in place)
  !> 
  !> usage call var%fix_irow(irow, nrow, sol_irow,rhs)
  !> 
  !> where:
  !> \param[in ] irow              -> integer. Index of 
  !>                                     row where fix the diag.
  !> \param[in ] nrow              -> integer. Nmb of eqs
  !> \param[in ] (optin.) sol_irow -> real. Sol in irow
  !> \param[out] (optin.) rhs      -> real(nrow). rhs to correct
  !>----------------------------------------------------------------
  subroutine fix_irow(this,irow,sol_irow,rhs)
    use Globals
    implicit none   

    class(spmat),                intent(inout) :: this
    integer,                     intent(in   ) :: irow
    real(kind=double), optional, intent(in   ) :: sol_irow
    real(kind=double), optional, intent(inout) :: rhs(this%nrow)
    !local 
    logical :: corr
    integer :: j

    corr = (present(rhs) .and. present(sol_irow))

    if (this%storage_system .eq. 'ssr' ) then
       ! elements before diagonal
       if (irow > 1) then
          do j = 1, this%ia(irow)-1
             if ( this%ja(j) .eq. irow ) then
                if ( corr ) then
                   rhs(this%ja(j)) = rhs(this%ja(j)) - &
                        this%coeff(j) * sol_irow
                end if
                this%coeff(j) = zero
             end if
          end do
       end if
       ! diagonal element
       this%coeff(this%ia(irow)) = one
       if ( corr ) rhs(irow) = sol_irow

       ! elements after diagonal 
       do j=this%ia(irow)+1, this%ia(irow+1)-1
          if ( corr ) then
             rhs(this%ja(j)) = rhs(this%ja(j)) - &
                  this%coeff(j) * sol_irow
          end if
          this%coeff(j) = zero
       end do

    else if ( this%storage_system .eq. 'csr' ) then
       this%coeff(this%idiag(irow)) = one
       do j=this%ia(irow), this%ia(irow+1)-1
          if ( this%ja(j) .ne. irow ) then
             if ( corr ) then
                rhs(this%ja(j)) = rhs(this%ja(j)) - &
                     this%coeff(j)*sol_irow
             end if
             this%coeff(j) = zero
          end if
       end do
    end if

  end subroutine fix_irow


  !>-------------------------------------------------------
  !> Procedure to operate diagonal scaling of spmat
  !>     M => left_D M right_D
  !> where D(M) is a diagonal matrix
  !> (pubblic procedure for type spmat, operate in place)
  !> 
  !> usage: var%diag_scale(diag_matrix)
  !>
  !> where: 
  !>\param[in] diag_matrix -> real(nrow). Diagonal mat. values
  !<-------------------------------------------------------
  subroutine diagonal_scale(this, lun_err, diagonal)
    use Globals
    implicit none
    class(spmat),      intent(inout) :: this
    integer,           intent(in   ) :: lun_err
    real(kind=double), intent(in   ) :: diagonal(this%ncol)
    !local
    logical rc
    integer res,i,j
    real(kind=double) ::  di

    if (.not. this%is_squared) then 
       rc = IOerr(lun_err, err_val, 'diagonal_scale', &
            'spmat is not squared',res)
    end if


    do i = 1, this%nrow
       di = diagonal(i)
       do j = this%ia(i), this%ia(i+1)-1
          this%coeff(j) = &
               di * this%coeff(j) * diagonal( this%ja(j) )
       end do
    end do

    do i=1,this%dim_kernel
       this%kernel(:,i)=this%kernel(:,i)/diagonal
    end do

  end subroutine diagonal_scale


  !>-------------------------------------------------------
  !> Procedure to operate matrix / matrix product.
  !>   B = D * A
  !> where D is a diagonal matrix and A is a sparse matrix
  !> (pubblic procedure for type spmat, operate in place)
  !> It converts the storage scheme in case of 
  !> A=symmetric matrix in ssr
  !> 
  !> 
  !> usage: var%DxM(lun_err,diag)
  !>
  !> where: 
  !>\param[in] lun_err -> integer. I/O unit for err. msg. 
  !>\param[in] diag    -> real. dimension(nequ). 
  !>                         Diagonal mat. values
  !<-------------------------------------------------------
  subroutine DxM( this, lun_err, diag)
    use Globals
    implicit none
    class(spmat),       intent(inout) :: this
    integer,            intent(in   ) :: lun_err
    real(kind= double), intent(in   ) :: diag(this%nrow)
    !local 
    integer i
    real(kind=double) :: diag_scal

    if (this%storage_system .eq. 'ssr' ) then
       call this%ssr2csr(lun_err)
    end if
    do i=1, this%nrow
       diag_scal = diag(i)
       this%coeff( this%ia(i): this%ia(i+1)-1) =  &
            diag_scal * &
            this%coeff( this%ia(i):this%ia(i+1)-1 )
    end do

    !
    ! in general this is not true
    ! 
    this%is_symmetric = .false.
    this%unitary_diag = .false.
  end subroutine DxM

  !>-------------------------------------------------------
  !> Procedure to operate matrix / matrix product.
  !>   B = A * D
  !> where D is a diagonal matrix and A is a sparse matrix
  !> (pubblic procedure for type spmat, operate in place)
  !> It converts the storage scheme in case of 
  !> A=symmetric matrix in ssr
  !> 
  !> 
  !> usage: var%DxM(lun_err,diag)
  !>
  !> where: 
  !>\param[in] lun_err -> integer. I/O unit for err. msg. 
  !>\param[in] diag    -> real. dimension(nequ). 
  !>                         Diagonal mat. values
  !<-------------------------------------------------------
  subroutine MxD( this, lun_err, diag)
    use Globals
    implicit none
    class(spmat),       intent(inout) :: this
    integer,            intent(in   ) :: lun_err
    real(kind= double), intent(in   ) :: diag(this%ncol)
    !local 
    integer i, j

    if (this%storage_system .eq. 'ssr' ) then
       call this%ssr2csr(lun_err)
    end if

    do i=1, this%nrow
       do j = this%ia(i), this%ia(i+1)-1
          this%coeff(j) =  &
               diag(this%ja(j)) * this%coeff( j )
       end do
    end do
    
    !
    ! right multiplication affect kernel
    !
    do i=1,this%dim_kernel
       this%kernel(:,i)=this%kernel(:,i)/diag
    end do


    !
    ! in general this is not true
    ! 
    this%is_symmetric = .false.
    this%unitary_diag = .false.
  end subroutine MxD



  !>-------------------------------------------------------------
  !> Procedure to compute Matrix vector product
  !>         vec_out = (M) times (vec_in)
  !> (public procedure for type spmat)
  !> 
  !> usage:
  !>     call 'var'%Mxv(vec_in,vec_out,[info])
  !>
  !> where 
  !> \param[in   ] vec_in          -> real. dimension('var'%ncol)
  !>                                  vector to be multiplied
  !> \param[inout] vec_out         -> real. dimension('var'%nrow)
  !>                                  vector (M) times (vec_in) 
  !> \param[in   ] (optional) info -> integer. Info number
  !>                                  in case of error   
  !<-------------------------------------------------------------
  subroutine Mxv(this,vec_in,vec_out, info)
    use Globals
    implicit none
    class(spmat),      intent(inout) :: this
    real(kind=double), intent(in   ) :: vec_in(this%ncol)
    real(kind=double), intent(inout) :: vec_out(this%nrow)
    integer, optional, intent(inout) :: info

    if (present(info)) info=0

    select case(this%storage_system) 
    case('ssr') 
       call axbsym(this%nrow,this%ncol,this%nterm,this%ia,this%ja,&
            this%coeff,vec_in,vec_out)
    case('csr') 
       call axbnsym(this%nrow,this%ncol,this%nterm,this%ia,this%ja,&
            this%coeff,vec_in,vec_out)
    end select
  contains

    ! vec_out = M * vec_in 
    subroutine axbsym(nrow,ncol,nterm,ia,ja,coef1,xvec,bvec)
      use Globals
      implicit none
      integer :: nrow,ncol,nterm
      integer :: ja(nterm),ia(nrow+1)
      real(kind=double) ::coef1(nterm),xvec(ncol),bvec(nrow)
      integer :: i,k,m,mm

      bvec=zero
      do k=1,nrow
         m=ia(k)
         mm=ia(k+1)-1
         bvec(k)=bvec(k)+coef1(m)*xvec(ja(m))
         do i=m+1,mm
            bvec(k)=bvec(k)+coef1(i)*xvec(ja(i))
            bvec(ja(i))=bvec(ja(i))+coef1(i)*xvec(k)
         end do
      end do

    end subroutine axbsym

    ! vec_out = M * vec_in 
    subroutine axbnsym(nrow,ncol,nterm,ia,ja,coef1,xvec,bvec)
      use GLobals
      implicit none
      integer nrow,ncol,nterm
      integer i,k
      real(kind=double)  :: coef1(nterm),xvec(ncol),bvec(nrow)
      integer ja(nterm),ia(nrow+1)

      bvec = zero
      do k=1,nrow
         do i=ia(k),ia(k+1)-1
            bvec(k)=bvec(k)+coef1(i)*xvec(ja(i))
         end do
      end do
    end subroutine axbnsym

  end subroutine Mxv

  !>-------------------------------------------------------------
  !> Procedure to compute Matrix vector product
  !>         vec_out = (M) times (vec_in)
  !> (public procedure for type spmat)
  !> 
  !> usage:
  !>     call 'var'%MTxv(vec_in,vec_out,[info])
  !>
  !> where 
  !> \param[in   ] vec_in          -> real. dimension('var'%ncol)
  !>                                  vector to be multiplied
  !> \param[inout] vec_out         -> real. dimension('var'%nrow)
  !>                                  vector (M) times (vec_in) 
  !> \param[in   ] (optional) info -> integer. Info number
  !>                                  in case of error   
  !<-------------------------------------------------------------
  subroutine MTxv(this,vec_in,vec_out, info)
    use Globals
    implicit none
    class(spmat),      intent(inout) :: this
    real(kind=double), intent(in   ) :: vec_in(this%nrow)
    real(kind=double), intent(inout) :: vec_out(this%ncol)
    integer, optional, intent(inout) :: info

    if (present(info)) info=0

    select case ( this%storage_system ) 
    case ('ssr')
       call axbsym(this%nrow,this%ncol,this%nterm,this%ia,this%ja,&
            this%coeff,vec_in,vec_out)
    case('csr')
       call axbnsym_transpose(this%nrow,this%ncol,&
            this%nterm,this%ia,this%ja,&
            this%coeff,vec_in,vec_out)
    end select
  contains

    ! vec_out = M * vec_in 
    subroutine axbsym(nrow,ncol,nterm,ia,ja,coef1,xvec,bvec)
      use Globals
      implicit none
      integer :: nrow,ncol,nterm
      integer :: ja(nterm),ia(nrow+1)
      real(kind=double) ::coef1(nterm),xvec(ncol),bvec(nrow)
      integer :: i,k,m,mm

      bvec=zero
      do k=1,nrow
         m=ia(k)
         mm=ia(k+1)-1
         bvec(k)=bvec(k)+coef1(m)*xvec(ja(m))
         do i=m+1,mm
            bvec(k)=bvec(k)+coef1(i)*xvec(ja(i))
            bvec(ja(i))=bvec(ja(i))+coef1(i)*xvec(k)
         end do
      end do

    end subroutine axbsym

    ! vec_out = M^t * vec_in 
    subroutine axbnsym_transpose(nrow,ncol,nterm,ia,ja,coef1,xvec,bvec)
      use GLobals
      implicit none
      integer nrow,ncol,nterm
      integer i,k,m,mm
      real(kind=double)  coef1(nterm),xvec(nrow),bvec(ncol)
      integer ja(nterm),ia(nrow+1)

      bvec = zero
      do k=1,nrow
         m=ia(k)
         mm=ia(k+1)-1
         do i=m,mm
            bvec(ja(i))=bvec(ja(i))+coef1(i)*xvec(k)
         end do
      end do
    end subroutine axbnsym_transpose

  end subroutine MTxv


  subroutine transpose_spmat(this, lun_err,perm)
    use Globals
    implicit none
    ! Input variables
    class(spmat),     intent(inout) :: this
    integer,          intent(in   ) :: lun_err
    integer,optional, intent(inout) :: perm(this%nterm)

    ! Local variables
    logical :: rc
    integer :: res
    integer :: i,j,ind
    character(len=1) :: new_triangular
    integer, allocatable :: IW(:),permutation(:)
    type(spmat) :: transpose
    
    ! check if the matrix id symmetric
    if ( this%is_symmetric ) return
    ! check if the matrix uses csr storage
    if ( this%storage_system  .ne. 'csr' ) then
       write(lun_err,*) 'Matrix not in csr format'
       stop
    end if
    ! wok array
    allocate (&
         IW(this%ncol+1),&
         permutation(this%nterm),&
         stat=res) 
    if(res .ne. 0) rc = IOerr(lun_err, err_alloc, 'transpose matrix', &
         'work array IW',res)


    ! set properties of transpose matrix
    select case (this%triangular)
    case ('N')
       new_triangular = 'N'
    case ('U')
       new_triangular = 'L'
    case ('L')
       new_triangular = 'U'
    end select

    ! initialized working matrix
    call transpose%init(lun_err,&
         this%ncol,this%nrow,& 
         this%nterm,&
         'csr',&
         triangular=new_triangular)  

    ! Initialize pointers
    transpose%ia = 0

    ! Count non-zeroes for each column of the input matrix
    do i = 1,this%nrow
       do j = this%ia(i),this%ia(i+1)-1
          transpose%ia(this%ja(j)) = transpose%ia(this%ja(j)) + 1
       enddo
    enddo

    ! Reset pointers
    IW(1) = 1
    do i = 2,transpose%nrow+1
       IW(i) = IW(i-1) + transpose%ia(i-1)
    enddo
    call SCOPY(transpose%nrow+1,IW,1,transpose%ia,1)


    ! Transpose coeffiecients
    do i = 1,this%nrow
       do j = this%ia(i),this%ia(i+1)-1
          ind = IW(this%ja(j))
          transpose%ja(ind) = i
          transpose%coeff(ind) = this%coeff(j)
          permutation(j) = ind 
          IW(this%ja(j)) = ind+1
       end do
    end do

    ! copy the working spmat into this
    select type(this)
    type is (spmat)
       this = transpose
    end select
    ! copy if present
    if (present(perm) ) perm=permutation


    ! free memory
    call transpose%kill(lun_err)
    deallocate (IW,permutation,stat=res) 
    if(res .ne. 0) rc = IOerr(lun_err, err_dealloc, 'transpose matrix', &
         'work array IW perm',res)

  end subroutine transpose_spmat


  !>-----------------------------------------------------------------
  !> Subroutine to solve linear system in the form
  !>      M x = b    or     M^T x = b
  !> with M triangular
  !>-------------------------------------------------------
  subroutine solve_triangular(this,lun_err,rhs,sol,&
                                ! optional arguments
       transpose, &
       inverse_diagonal)
    use Globals
    !use Timing
    implicit none
    class(spmat),                intent(in   ) :: this
    integer,                     intent(in   ) :: lun_err
    real(kind=double),           intent(in   ) :: rhs(this%nrow)
    real(kind=double),           intent(inout) :: sol(this%ncol)
    character(len=1),  optional, intent(in   ) :: transpose    
    real(kind=double), optional, intent(in   ) :: inverse_diagonal(this%nrow)
    ! local
    logical :: rc
    integer :: i,j,k,m,nequ
    integer :: idiag,ifirst,ilast,irow
    real(kind=double) :: scr
    !type(Tim) :: trisol1,trisol2



    nequ = this%nrow
    !
    ! checks
    !
    if ( (this%triangular .eq. 'N') ) then
       rc = IOerr(lun_err, err_inp, 'solve_triangular', &
            'non-triangualr matrix passed')
    end if

    if ( .not. this%is_squared ) then
       rc = IOerr(lun_err, err_inp, 'solve_triangular', &
            'non-squared matrix passed')
    end if


    !
    ! case without inverse of the diagonal
    !
    if ( .not. present (inverse_diagonal) ) then 
       if ( present(transpose) .and. (transpose .eq. 'T') ) then
          select case (this%triangular)
          case('U')
             !
             ! solve system U^t x = rhs
             !
             !call trisol1%init()
             !call trisol1%set('start')
             sol = zero
             do k = 1, nequ
                idiag = this%ia(k)
                ilast = this%ia(k+1) - 1
                sol(k) = ( rhs(k) - sol(k) ) / this%coeff(this%ia(k))
                do m = idiag+1, ilast
                   sol(this%ja(m)) = sol(this%ja(m)) + &
                        this%coeff(m) * sol(k)
                end do
             end do
             !call  trisol1%set('stop')
             !call  trisol1%info(lun_err,'      no  diag inv U^T')
             !call  trisol1%kill()
          case('L')
             !
             rc = IOerr(lun_err, err_inp, 'solve_triangular', &
                  ' not implemated for L^t x = b')
          end select
       else
          !
          ! solve with given matrix 
          !
          select case (this%triangular)
          case('L')
             !
             ! solve system L x = rhs
             !
             sol=zero
             do k = 1, nequ
                ! select range
                ifirst = this%ia(k)
                idiag  = this%ia(k+1) - 1  
                sol(k) = rhs(k)
                do m = ifirst,idiag-1
                   sol(k) = sol(k) - this%coeff(m) * sol(this%ja(m))
                end do
                sol(k) = sol(k) / this%coeff(idiag)
             end do
          case('U')
             !
             ! solve system U x = rhs
             !
             !call trisol1%init()
             !call trisol1%set('start')
             do k = nequ,1,-1
                idiag  = this%ia(k  )
                ilast  = this%ia(k+1) - 1
                ! sol(k) =  a_{k,k}^{-1} (
                !           rhs(k) - 
                !           \sum_{j} a_{irow,j} sol_{j} )
                scr = zero
                !
                ! EF change the order of the loop ??
                ! use do m = ilast, idiag+1,-1
                !
                do m = ilast, idiag+1,-1
                   scr = scr + this%coeff(m)*sol(this%ja(m))
                end do
                sol(k) = ( rhs(k) - scr ) / this%coeff(idiag)
             end do
             !call  trisol1%set('stop')
             !call  trisol1%info(lun_err,'     no  diag inv U')
             !call  trisol1%kill()
          end select
       end if
    else
       !
       ! use given inverse_diagonal 
       ! 
       if ( present(transpose) .and. (transpose.eq.'T') ) then
          select case (this%triangular)
          case('U')
             !
             ! solve system U^t x = rhs
             !
             sol = zero
             !
             do k = 1, nequ
                ! line range
                idiag = this%ia(k)
                ilast = this%ia(k+1) - 1
                sol(k) = ( rhs(k) - sol(k) ) * inverse_diagonal(k)
                do m = idiag+1, ilast
                   sol(this%ja(m)) = sol(this%ja(m)) + &
                        this%coeff(m) * sol(k) 
                end do
             end do
!!$             do k = 1, nequ
!!$                sol(k) = sol(k) * inverse_diagonal(k)
!!$             end do
          case('L')
             !
             rc = IOerr(lun_err, err_inp, 'solve_triangular', &
                  ' not implemated for L^t x = b')
          end select
       else
          !
          ! solve with given matrix 
          !
          select case (this%triangular)
          case('L')
             !
             ! solve system L x = rhs
             !
             sol=zero

             do irow = 1, nequ
                ! select line range
                ifirst = this%ia(irow)
                idiag  = this%ia(irow+1) - 1  

                ! commpute  
                !   sol(k) =  a_{k,k}^{-1} (
                !               rhs(k) - 
                !               \sum_{j \neq irow } a_{irow,j} sol_{j} )
                sol(irow)=rhs(irow)
                do m = ifirst,idiag-1
                   sol(irow) = sol(irow) - this%coeff(m) * sol(this%ja(m))
                end do
                sol(irow) = sol(irow) * inverse_diagonal(irow)
             end do
          case('U')
             !
             ! solve system U x = rhs
             !
             do irow = nequ,1,-1
                idiag  = this%ia(irow  )
                ilast  = this%ia(irow+1) - 1
                ! sol(k) =  a_{k,k}^{-1} (
                !           rhs(k) - 
                !           \sum_{j} a_{irow,j} sol_{j} )
                sol(irow) = rhs(irow)
                do m = ilast, idiag+1, -1
                   sol(irow) = sol(irow) - this%coeff(m) * sol(this%ja(m))
                end do
                sol(irow) = sol(irow) * inverse_diagonal(irow) 
             end do
          end select
       end if
    end if

  end subroutine solve_triangular


  !>-----------------------------------------------------------------
  !> Subroutine to solve triangular linear system in the form
  !>      M x = b    or     M^T x = b
  !> with M triangular with diagonal equal one
  !>-------------------------------------------------------
  subroutine solve_triangular_unitdiag(this,lun_err,rhs,sol,transpose)
    use Globals
    !use Timing
    implicit none
    class(spmat),                intent(in   ) :: this
    integer,                     intent(in   ) :: lun_err
    real(kind=double),           intent(in   ) :: rhs(this%nrow)
    real(kind=double),           intent(inout) :: sol(this%ncol)
    character(len=1),  optional, intent(in   ) :: transpose    
    ! local
    logical :: rc
    integer :: i,j,k,m,nequ
    integer :: idiag,ifirst,ilast,irow
    real(kind=double) :: scr
    !type(Tim) :: trisol1,trisol2
        

    nequ = this%nrow
    !
    ! checks
    !
    if ( (this%triangular .eq. 'N') ) then
       rc = IOerr(lun_err, err_inp, 'solve_triangular', &
            'non-triangualr matrix passed')
    end if

    if ( .not. this%is_squared ) then
       rc = IOerr(lun_err, err_inp, 'solve_triangular', &
            'non-squared matrix passed')
    end if

    if ( .not. this%unitary_diag ) then
       rc = IOerr(lun_err, err_inp, 'solve_triangular', &
            'not unitary diagonal matrix passed')
    end if

    if ( present(transpose) .and. (transpose .eq. 'T') ) then
       select case (this%triangular)
       case('U')
          !call trisol1%init()
          !call trisol1%set('start')
          !
          ! solve system U^t x = rhs
          !
          sol = zero
          do k = 1, nequ
             idiag = this%ia(k)
             ilast = this%ia(k+1) - 1
             ! removed inversion
             !sol(k) = ( rhs(k) - sol(k) ) / this%coeff(idiag)
             sol(k) = rhs(k) - sol(k) 
             do m = idiag+1, ilast
                sol(this%ja(m)) = sol(this%ja(m)) + &
                     this%coeff(m) * sol(k)
             end do
          end do
          !call  trisol1%set('stop')
          !call  trisol1%info(lun_err,'    yes diag inv U^T')
          !call  trisol1%kill()
          return
       case('L')
          !
          rc = IOerr(lun_err, err_inp, 'solve_triangular', &
               ' not implemated for L^t x = b')
       end select
    else
       !
       ! solve with given matrix 
       !
       select case (this%triangular)
       case('L')
          !
          ! solve system L x = rhs
          !
          sol=zero
          do k = 1, nequ
             ! select range
             ifirst = this%ia(k)
             idiag  = this%ia(k+1) - 1  
             sol(k)=rhs(k)
             do m = ifirst,idiag-1
                sol(k) = sol(k) - this%coeff(m) * sol(this%ja(m))
             end do
             ! removed inversion{
             ! sol(k) = sol(k) / this%coeff(idiag)
             ! }
          end do
       case('U')
          !
          ! solve system U x = rhs
          !
          do k = nequ,1,-1
             ! sol(k) =  a_{k,k}^{-1} (
             !           rhs(k) - 
             !           \sum_{j} a_{k,j} sol_{j} )
             !
             idiag  = this%ia(k  )
             ilast  = this%ia(k+1) - 1

             sol(k) = rhs(k)
             do m = ilast, idiag+1, -1
                sol(k) = sol(k) - this%coeff(m)*sol(this%ja(m))
             end do
             ! removed{
             !sol(k) = sol(k) / this%coeff(idiag)
             !} 
          end do
       end select
    end if
    write(103,*) sol

  end subroutine solve_triangular_unitdiag

  !>-------------------------------------------------------
  !> Subroutine for incomplete cholesky factorization
  !> of a symmetric matrix A. 
  !> It computes the approximated factor U giving
  !>                  A ~ U^t U
  !> with U a upper triangualr matrix in csr format.
  !> If the optional array "diagonal" is passed the
  !> the output are
  !> 1- an upper triangular matrix with one on the diagonal
  !>            upper    = diag^{-1}(U) U
  !> 2 -the diagonal array containg the diagonal of U 
  !             diagonal = diag(U) 
  !>-------------------------------------------------------
  subroutine incomplete_cholesky (this,&
       lun_err,n_fillin,tol_fillin,info,&
       upper,&
       ! optional arguements
       diagonal)
    use Globals
    implicit none
    class(spmat), target,        intent(in   ) :: this
    integer,                     intent(in   ) :: lun_err
    integer,                     intent(in   ) :: n_fillin
    real(kind=double),           intent(in   ) :: tol_fillin
    integer,                     intent(inout) :: info
    type(spmat),                 intent(inout) :: upper
    real(kind=double), optional, intent(inout) :: diagonal(this%nrow)
    !local
    logical :: rc
    integer :: res
    integer :: nterm_max, nterm ,nterm_out
    integer :: irow,nrow
    integer,           allocatable :: jlu(:),jw(:),iscr(:),iend(:)
    integer,           allocatable :: iw1(:),jw1(:),jw2(:),jw3(:)
    integer, target,   allocatable :: iwork(:)
    real(kind=double), allocatable :: alu(:),w(:)
    integer, pointer :: idiag(:)
    type(spmat) :: ssr_temp


    !
    ! checks
    !
    if ( .not. ( this%is_squared ) ) &
         rc = IOerr(lun_err, err_inp, 'incomplete_cholesky ', &
         'not square matrix passed')

    if ( .not. (this%is_symmetric) ) &
         rc = IOerr(lun_err, err_inp, 'incomplete_cholesky ', &
         'not symmmetric matrix passed')
    
    !
    ! dimensions
    !
    nrow      = this%nrow
    nterm_max = this%nterm + (n_fillin+1) * nrow
    !
    ! work arrays for both storage system
    !
    allocate(&
         alu(nterm_max), &
         jlu(nterm_max),&
         jw1(nrow),&
         jw2(nrow),&
         jw3(nrow),&
         w(nrow+1),&
         jw(2*nrow),&
         iscr(nrow+1),&
         stat=res)
    if ( res .ne. 0) &
         rc = IOerr(lun_err, err_alloc, 'incomplete_cholesky ', &
         'work arrays alu, jlu w jw',res)

    
    if ( n_fillin .gt. 0) then      
       if ( this%storage_system .eq. 'ssr' ) then
          !
          ! 1 : assign integer pointer idiag, ia, ja, iw1
          ! 

          !
          ! 1.1 : allocate extra work arrays 
          !
          allocate(&
               iw1(this%nterm),&     ! ja_transpose
               iwork(nrow+1),&  ! ia_work to cycle ja_transpose
               stat=res)
          if ( res .ne. 0) &
               rc = IOerr(lun_err, err_alloc, 'incomplete_cholesky ', &
               'work arrays iw1, iwork',res)

          !
          ! 1.2 : build iw1, iwork
          ! iw1 shuold act like ja for the lower part of matrix
          ! iwork contains the incides of the startin indices for each row
          !
          call my_unpackja(nrow,this%nterm, this%ia, this%ja, iscr(1:nrow),&
               iw1, iwork)

          !
          ! 2 : compute factorization
          !
          
          info = 0
          call my_ssr_incomplete_cheolesky(&
               this%nrow,&
               this%nterm,&
               this%coeff,&
               this%ja,&
               this%ia,&
               iwork,&
               n_fillin,tol_fillin,&
               alu,&
               jlu,&
               nterm_max,&
               w,&
               iw1,& ! already computed 
               jw1,&
               jw2,&
               jw3,&
               info,lun_err)

          !
          ! free local memory
          !
          !write(0,*) 'problema: deallocate(iw1,iwork) - 22modSparseMatrix'
          deallocate(iw1, iwork, stat=res)
          if ( res .ne. 0) &
               rc = IOerr(lun_err, err_dealloc, 'incomplete_cholesky ', &
               'work arrays iw1, iwork',res)

       else if ( this%storage_system .eq. 'csr') then
          !
          ! 1 : assign integer pointer idiag, ia, ja, iw1
          ! 
          
          !
          ! diagonal pointer
          !
          if ( .not. allocated(this%ind_diag) ) then
             !
             ! 1.1 : allocate extra work arrays 
             !
             allocate(&
                  iwork(nrow),&  ! idiag
                  stat=res)
             if ( res .ne. 0) &
                  rc = IOerr(lun_err, err_alloc, 'incomplete_cholesky ', &
                  'work arrays ',res)
             !
             ! build diagonal
             !
             do irow = 1, nrow
                iwork(irow) = this%idiag(irow)
                if ( iwork(irow) .eq. 0) &
                     rc = IOerr(lun_err, err_inp, 'incomplete_cholesky ', &
                     'not zero diagonal term on line', irow)
             end do
             idiag => iwork
          else
             idiag => this%ind_diag
          end if
          !write(0,*) idiag

          !
          ! compute factorization
          !
          info = 0
          call my_csr_incomplete_cheolesky(&
               this%nrow,&
               this%nterm,&
               this%coeff,&
               this%ja,&
               this%ia,&
               idiag,&
               n_fillin,tol_fillin,&
               alu,&
               jlu,&
               nterm_max,&
               w,&
               jw1,&
               jw2,&
               jw3,&
               info,lun_err)

          !write(0,*) info
                              
          !
          ! free local memory
          !
          idiag => null()
          if ( allocated(iwork) ) then
             deallocate(iwork,stat=res)
              if ( res .ne. 0) &
                  rc = IOerr(lun_err, err_dealloc, 'incomplete_cholesky ', &
                  'work array iwork',res)
          end if
       end if
                        
       if (info.eq.0) then
          !
          ! convert mssr to csr
          !
          call mss2scaled_upper(lun_err,&
               nrow,&
               nterm_max,&            
               jlu,&
               alu,&
               upper)

          !
          ! aluU(1:nrow) contains the square of the diagonal term
          ! of the upper factor of the choelsky factorization
          !

          if ( present (diagonal) ) then
             !
             ! case for D=diag(U),  \tilde{U}=diag(U)^{-1} U 
             !
             diagonal(1:nrow) = sqrt( alu(1:nrow) )
          else
             !
             ! case for U factorization
             !
             ! use alu as scratch array to compute roots
             alu(1:nrow) = sqrt( alu(1:nrow) )
             call upper%DxM( lun_err,alu(1:nrow) )
             upper%unitary_diag = .false.
          end if

       else if ( info .lt. 0) then
          !
          ! info in case of lu factorization error 
          !
          write(lun_err,*) ' Error in IC construction' 
          select case (info) 
          case (-1)
             rc = IOerr(lun_err, err_val, 'build_ic_fillin', &
                  ' The elimination process has generated'//&
                  ' a row in U whose length is .gt.  n)')
          case (-2)  
             rc = IOerr(lun_err, err_val, 'build_ic_fillin', &
                  'The matrix U overflows the array U')
          case (-3)  
             rc = IOerr(lun_err, err_val, 'build_ic_fillin', &
                  ' Illegal value for n_fillin')
          case(-4)   
             rc = IOerr(lun_err, err_val, 'build_ic_fillin', &
                  'zero row encountered')
          case(-5) 
             rc = IOerr(lun_err, err_val, 'build_ic_fillin', &
                  'Non Positive Diagonal Element')
          end select
       end if

    else
       !
       ! No-fillin algorithm 
       !
       if ( this%storage_system .eq. 'ssr' ) then
          !
          ! init upper factor
          !
          if ( upper%is_initialized ) then !!!!!
             call upper%kill(lun_err)
          end if
          call upper%init(lun_err, &
               nrow, nrow, this%nterm,&
               'csr',&
               dim_kernel=0,&
               is_symmetric=.false.,&
               triangular='U')

          !
          ! build U factor 
          !
          upper%ia = this%ia
          upper%ja = this%ja
          call kersh_loc(lun_err,nrow,this%nterm,&
               this%ia,&
               this%ja,&
               this%coeff,&
               upper%coeff)
       else
          !
          ! assign pointer to diagonal terms of this
          !
          select type( this )
          type is (spmat)
             ssr_temp = this
          end select
          call ssr_temp%csr2ssr(lun_err)

          !
          ! init upper factor
          !
          if ( upper%is_initialized ) then!!!!!
             call upper%kill(lun_err)
          end if
          call upper%init(lun_err, &
               nrow, nrow, ssr_temp%nterm,&
               'csr',&
               dim_kernel=0,&
               is_symmetric=.false.,&
               triangular='U')
          

          !
          ! build U factor 
          !
          call kersh_loc(lun_err,nrow,ssr_temp%nterm,&
               ssr_temp%ia,&
               ssr_temp%ja,&
               ssr_temp%coeff,&
               upper%coeff)
          upper%ia = ssr_temp%ia
          upper%ja = ssr_temp%ja


          !
          ! free memory
          ! 
          call ssr_temp%kill(lun_err)

       end if

       !
       ! U => diag(U) , diag(U)^{-1} U 
       !
       if ( present (diagonal) ) then
          !
          ! case for D=diag(U),  \tilde{U}=diag(U)^{-1} U 
          !
          call upper%get_diagonal(diagonal)
          diagonal = one / diagonal
          call upper%DxM(lun_err, diagonal)
          upper%unitary_diag = .true. 
          diagonal = one / diagonal
       end if
    end if

    !
    ! free memory
    !
    deallocate(&
         alu, &
         jlu,&
         jw1,&
         jw2,&
         jw3,&
         w,&
         jw,&
         iscr,&
         stat=res)
    if ( res .ne. 0) rc = IOerr(lun_err, err_dealloc, &
         'incomplete_cholesky', &
         'work arrays alu, jlu, w, jw etc.',res)
 
  contains
    !--------------------------------------------------------------
    ! Subroutine computing incomplete upper triangular factor
    ! of the Cholesky decomposition 
    !            A ~ U^T U
    ! Output matrix is stored in MSS( modified storage stystem) 
    ! containg 
    !         diag(U)**2     [ in array U(1     :nequ) ]
    !         diag(U)^{-1} U [ in array U(nequ+2:iwk ) ]
    !---------------------------------------------------------------
    subroutine my_csr_incomplete_cheolesky(n,nterm,&
         a,ja,ia,idiag,&
         lfil,droptol,&
         U,ju,iwk,&
         w,&
         jw1,jw2,jw3,ierr,iout)
      !------------------------------------------------------------
      !      
      implicit none 
      !      
      integer n, nterm 
      integer lfil, iwk, ierr, iout
      integer ja(nterm), ia(n+1),idiag(n),iend(n)
      integer ju(iwk), iw1(nterm), jw1(n), jw2(n), jw3(n+1)   
      !      
      real*8  a(nterm), U(iwk), w(n+1), droptol
      !      
      !------------------------------------------------------------
      !                    *** SYMILUT preconditioner ***         
      !      incomplete U^tU factorization with dual truncation mechanism 
      !-------------------------------------------------------------------
      !     Author: Carlo Janna *December, 12, 2005                       
      !-------------------------------------------------------------------
      ! PARAMETERS                                                        
      !
      ! on entry:
      !========== 
      ! n       = integer. The row dimension of the matrix A.  
      !
      ! nterm   = integer. The number of non-zero of the matrix A. 
      !
      ! a       = ceofficient of matrix 
      ! ja      = column index 
      ! first_row = indecx of the first nonzero element of a for each row
      ! idaig     = 
      !
      ! lfil    = integer. Fill-in parameter. Each row of L and each row
      !           of U will have a maximum of lfil elements (excluding the 
      !           diagonal element). lfil must be .ge. 0.
      !           ** WARNING: THE MEANING OF LFIL HAS CHANGED 
      !           WITH RESPECT T0 EARLIER VERSIONS. 
      !
      ! droptol = real*8. Sets the threshold for dropping small terms in the
      !           factorization. See below for details on dropping strategy.
      !
      !  
      ! iwk     = integer. The lengths of arrays U and ju. If the arrays
      !           are not big enough to store the U^tU factorizations, symilut
      !           will stop with an error message. Last term in ju is used 
      !           to deactivate jw3 work array. The space avaliable for the
      !           factorization is (iwk-1).
      !
      ! On return:
      !===========
      !
      ! U       = matrix stored in Modified Symmetric Sparse Row (MSSR) format
      !           containing the U factor. 
      !           U(1:n) contains the squared of the diagonal of the factor U.
      !           Extra-Diagonal elements are stored in
      !           U(n+2:nterm) and they scaled by the inverse of the diagonal of U
      !           M = D+T with D=diag(U)^2 T=D(U)^{-1} U
      !           U(n+1) is unused.
      !
      ! ju      = integer array. The first n+1 position contain the pointers to
      !           the beginning of each row of U after the diagonal. ju(n+2:*)
      !           are the column indeces of the matrix factor U (excluding the
      !           diagonal entry). [ju(n+1)-1] is the dimension of vectors U and
      !           ju.
      !
      ! ierr    = integer. Error message with the following meaning.
      !           ierr  = 0    --> successful return.
      !           ierr  = -1   --> Error. input matrix may be wrong.
      !                            (The elimination process has generated a
      !                            row in U whose length is .gt.  n).
      !           ierr  = -2   --> The matrix U overflows the array U.
      !           ierr  = -3   --> Illegal value for lfil.
      !           ierr  = -4   --> zero row encountered.
      !           ierr  = -5   --> Non Positive Diagonal Element
      !
      ! work arrays:
      !=============
      ! iw1     = integer work array of lenght nterm+1
      ! jw1     = integer work array of length n (max number of non-zero).
      ! jw2     = integer work array of length n.
      ! jw3     = integer work array of length n.
      ! w       = real work array of length n+1 (max number of non-zero).
      !  
      !----------------------------------------------------------------------
      ! w, jw1      store the working array [1:ii-1 = L-part, ii:n = u] 
      ! jw2         stores nonzero indicators
      ! jw3         stores pointers to the "interesting" part of row of U 
      !             factor
      ! iw1         stores lines to be eliminated for each line.
      !             first n positions are lenght of each line of iw1
      !                                 "POINTERS???"
      ! 
      ! Notes:
      ! -----
      ! The diagonal elements of the input matrix must be  nonzero (at least
      ! 'structurally'). 
      !
      !----------------------------------------------------------------------* 
      !---- Dual drop strategy works as follows.                             *
      !                                                                      *
      !     1) Theresholding in L and U as set by droptol. Any element whose *
      !        magnitude is less than some tolerance (relative to the abs    *
      !        value of diagonal element in u) is dropped.                   *
      !                                                                      *
      !     2) Keeping only the largest lfil elements in the i-th row of U   * 
      !        (excluding diagonal elements).                                *
      !                                                                      *
      ! Flexibility: one  can use  droptol=0  to get  a strategy  based on   *
      ! keeping  the largest  elements in  each row  of U. Taking            *
      ! droptol .ne.  0 but lfil=n will give  the usual threshold strategy   *
      ! (however, fill-in is then umpredictible).                            *
      !----------------------------------------------------------------------*
      !
      !  Locals
      integer  ju0, rowHead, rowMid , rowEnd, jrow , jcol, jpos
      integer  lenu , lenu0 , lenl , lenght, ncut
      integer  i,  ii, j,  jj, k, kk, j1, j2 , j3 , j4
      !      
      real*8   tnorm, abstol, s, fact, zero 
      !      
      parameter(zero = 0.0)
      !
      !-----------------------------------------------------------------------
      !  Initialize ju0 (points to next element to be added to U,ju)
      !  and pointer array.
      !-----------------------------------------------------------------------
      ju0 = n+2 
      ju(1) = ju0
      jw3(1) = ju0 
      !
      !  Initialize nonzero indicator array. 
      !
      do j=1,n 
         jw2(j)  = 0
      enddo
      !
      !  Set ju(iwk) to deactivate jw3 pointers
      !
      ju(iwk) = -1      
      !--------------------------------------------------------------
      !  Beginning of main loop.
      !-------------------------------------------------------------
      jw1=0

      do ii = 1, n
         j1 = idiag(ii)
         j2 = ia(ii+1) - 1
         !j2 = iend(ii)
         j3 = ia(ii)
         j4 = idiag(ii)-1
         
         !  Calculate the norm of the ii-th row         
         tnorm = 0.0d0
         do  k=j1,j2
            tnorm = tnorm+abs(a(k))
         end do
         !         
         if (tnorm .eq. 0.0) goto 999
         tnorm = tnorm/real(j2-j1+1)
         !     
         !  Unpack Upper and Lower part of row of A in array w 
         !   
         lenu  = j2 - j1 + 1
         lenu0 = lenu
         lenl  = j4 - j3 + 1         
         k = ii-1          
         do j = j1,j2         
            w(j-j1+ii)   = a(j)
            jw1(j-j1+ii) = ja(j)
            k = k + 1         
            jw2(ja(j))  = k
         end do
          

         k = 0         
         do j = j3,j4   
            ! in case of csr iw1=ja
            jw1(j-j3+1)  = ja(j)
            k = k + 1         
            jw2(ja(j))  = k
         end do
!!$         write(*,'(10(1x,e8.1))') (w(k),k=1,n)  
!!$         write(*,'(10I6)') (jw1(k),k=1,n) 
!!$         write(*,'(10I6)') (jw2(k),k=1,n) 
!!$         write(*,*) ' ' 
!
         !----BLAS-Method-----------------------------------------------------
         !  Unpack Upper and Lower part of row of A in array w 
         !         lenu  = j2 - j1 + 1
         !         lenu0 = lenu
         !         lenl  = j4 - j3 + 1
         !         call DCOPY(lenu,a(j1),1,w(ii),1)
         !         call SCOPY(lenu,ja(j1),1,jw1(ii),1)
         !         call SCOPY(lenl,iw1(j3),1,jw1(1),1)
         !   
         !  Set non-zero indicators
         !         
         !         k = ii-1
         !         do j = j1,j2
         !            k = k + 1
         !            jw2(ja(j)) = k
         !         enddo
         !         k = 0
         !         do j = j3,j4
         !            k = k + 1
         !            jw2(iw1(j)) = k
         !         enddo
         !----End of BLAS-Method----------------------------------------------
         !
         !  Eliminate previous rows
         !
         jj = 0 
150      jj = jj+1
         if (jj .gt. lenl) goto 160
         !
         !  Determine the smallest column index
         !         
         jrow = jw1(jj)
         k = jj
         do j = jj+1,lenl
            if (jw1(j) .lt. jrow) then
               jrow = jw1(j)
               k = j
            endif
         enddo

         if (k .ne. jj) then
            ! Exchange in jw1
            j = jw1(jj)
            jw1(jj) = jw1(k)
            jw1(k) = j
            ! Exchange in jw1
            jw2(jrow) = jj
            jw2(j) = k
         endif
         !write(*,*) 'ii, jrow',ii, jrow

         !
         !  Zero out element in row by setting jw2(jrow) to zero
         !
         jw2(jrow) = 0
         !         
         !  Get the leading term of the row if it exists a(jj,ii)
         !  
         rowHead = ju(jrow)
         rowMid  = jw3(jrow)
         if (ju(rowMid) .eq. ii) then  
            !------------------------------------------------------------------
            !  Mettere la possibilità di saltare 
            !  questo ciclo se il termine pivotale
            !  è troppo piccolo: fact = U(rowMid)
            !-----------------------------------------------------------------
            !
            ! Exists, mark fill-in terms of the current row before the
            ! diagonal            
            !
            do k = rowHead,rowMid-1
               jcol = ju(k)
               jpos = jw2(jcol)
               if (jpos .eq. 0) then
                  ! This is a fill-in term
                  lenl = lenl+1
                  if (lenl .gt. n) goto 996
                  jw1(lenl) = jcol
                  jw2(jcol) = lenl
               endif
            enddo
            !                    
            ! Combine current row ii with row jrow
            !

            fact = U(rowMid)*U(jrow)
            ! fact = U(rowMid) * U(jrow) * sqrt(U(jrow))
            rowEnd = ju(jrow+1) - 1
            do k = rowMid,rowEnd
               !
               !write(*,*) ii, rowMid,ju(rowMid), jrow,  ju(k), k
               s = fact * U(k)
               jcol = ju(k)
               jpos = jw2(jcol)
               if (jpos .eq. 0) then
                  !                          
                  ! This is a fill-in element
                  !
                  lenu = lenu + 1
                  if (lenu .gt. n) goto 996
                  i = ii+lenu-1
                  jw1(i) = jcol
                  jw2(jcol) = i
                  w(i) = -s
                  !                     
               else
                  !
                  ! This is not a fill-in element
                  !
                  w(jpos) = w(jpos) - s
                  !                     
               endif
               !                  
            enddo
            !               
            ! Update pointer jw3
            !
            jw3(jrow) = jw3(jrow) + 1
            if (jw3(jrow) .gt. rowEnd) then
               ! Deactivate pointer
               jw3(jrow) = iwk
            endif
            !               
         endif
         !           
         goto 150
160      continue
         !
         !  Reset non-zero indicators
         !
         do k = 1,lenu
            jw2(jw1(ii+k-1)) = 0
         enddo
         !
         !  Store the diagonal term
         !
         s = w(ii)         
         if (s .le. zero) then
            write(iout,100) ii,s
            write(6,100) ii,s
            write(6,*) 'ilu error'
            !stop
            ierr = -5  
            write(iout,101) tnorm
            return
         endif
         U(ii) = s
         s = 1.d0/s
         !
         !  Apply dropping strategy
         !
         abstol = s * tnorm * droptol  ! Set value of true tolerance
         lenght = 0 
         do k = ii+1,ii+lenu-1
            if (abs(w(k)) .gt. abstol) then
               lenght = lenght+1
               w(ii+lenght)   = w(k)
               jw1(ii+lenght) = jw1(k)
            endif
         enddo
         !
         !  Set lenu = number of elements bigger than tolerance        
         !
         lenu = lenght
         ncut = min0(lenu,lenu0+lfil-1)
         !----------------------------------------------------------------
         !  May be interesting to save a fixed number of non-zero 
         !  for each line
         ! ncut = min0(lenu,lfil)         
         !----------------------------------------------------------------
         !
         !  Select the ncut biggest elements only if lenu gt zero
         !
         if (lenu .gt. 0) then
            jpos = ii+1
            call qsplit(w(jpos),jw1(jpos),lenu,ncut)
            !          
            !  Order them in increasing jw1
            !
            !call HPSORT(jw1(jpos),w(jpos),ncut)
            call sort_row(ncut,&
                 jw1(jpos:jpos+ncut-1),&
                 w(jpos:jpos+ncut-1))
            !
            !  Scale and store in the factor
            !         
            if ((ju0+ncut-1) .gt. iwk-1) goto 997
            do k = 0,ncut-1
               ! EF
               U(ju0+k)  = s*w(jpos+k)
               ju(ju0+k) = jw1(jpos+k) 
            enddo
            !            
            !----BLAS-Method-------------------------------------------------
            !  Scale by the pivotal term and store
            !            if ((ju0+ncut-1) .gt. iwk-1) goto 997
            !            call DCOPY(ncut,w(jpos),1,U(ju0),1)
            !            call DSCAL(ncut,s,w(jpos),1)
            !            call SCOPY(ncut,jw1(jpos),1,ju(ju0),1)
            !----End of BLAS-Method------------------------------------------
            !            
         endif
         !
         !  Update pointer to beginning of next row of U    
         !
         ju0 = ju0 + ncut
         ju(ii+1) = ju0 
         jw3(ii+1) = ju0 
         !
      enddo


      !-----------------------------------------------------------------------
      !     end main loop
      !-----------------------------------------------------------------------
      !
      !     Successful run
      !
      ierr = 0
      !     
      return
      !
      !     incomprehensible error. Matrix must be wrong.
      !     
996   ierr = -1
      return
      !    
      !     insufficient storage in U.
      !     
997   ierr = -2
      return
      !     
      !     illegal lfil entered.
      !     
998   ierr = -3
      return
      !     
      !     zero row encountered
      !     
999   ierr = -4
      return
      !----------------end-of-ilut--------------------------------------------
100   format(' NULL OR NEGATIVE DIAGONAL ELEMENT: I,J =',I8,2X,E16.5)
101   format(' NORM OF ROW =  ',E16.8)
      !-----------------------------------------------------------------------
    end subroutine my_csr_incomplete_cheolesky

    

    !--------------------------------------------------------------
    ! Subroutine computing incomplete upper triangular factor
    ! of the Cholesky decomposition 
    !            A ~ U^T U
    ! Output matrix is stored in MSS( modified storage stystem) 
    ! containg 
    !         diag(U)**2     [ in array U(1     :nequ) ]
    !         diag(U)^{-1} U [ in array U(nequ+2:iwk ) ]
    !---------------------------------------------------------------
    subroutine my_ssr_incomplete_cheolesky(n,nterm,&
         a,ja,ia,ia_work,&
         lfil,droptol,&
         U,ju,iwk,&
         w,&
         iw1,jw1,jw2,jw3,ierr,iout)
      !------------------------------------------------------------
      !      
      implicit none 
      !      
      integer n, nterm 
      integer lfil, iwk, ierr, iout
      integer ja(nterm), ia(n+1), ia_work(n+1)
      integer ju(iwk), iw1(nterm), jw1(n), jw2(n), jw3(n+1)   
      !      
      real*8  a(nterm), U(iwk), w(n+1), droptol
      !      
      !------------------------------------------------------------
      !                    *** SYMILUT preconditioner ***         
      !      incomplete U^tU factorization with dual truncation mechanism 
      !-------------------------------------------------------------------
      !     Author: Carlo Janna *December, 12, 2005                       
      !-------------------------------------------------------------------
      ! PARAMETERS                                                        
      !
      ! on entry:
      !========== 
      ! n       = integer. The row dimension of the matrix A.  
      !
      ! nterm   = integer. The number of non-zero of the matrix A. 
      !
      ! a       = ceofficient of matrix 
      ! ja      = column index 
      ! first_row = indecx of the first nonzero element of a for each row
      ! idaig     = 
      !
      ! lfil    = integer. Fill-in parameter. Each row of L and each row
      !           of U will have a maximum of lfil elements (excluding the 
      !           diagonal element). lfil must be .ge. 0.
      !           ** WARNING: THE MEANING OF LFIL HAS CHANGED 
      !           WITH RESPECT T0 EARLIER VERSIONS. 
      !
      ! droptol = real*8. Sets the threshold for dropping small terms in the
      !           factorization. See below for details on dropping strategy.
      !
      !  
      ! iwk     = integer. The lengths of arrays U and ju. If the arrays
      !           are not big enough to store the U^tU factorizations, symilut
      !           will stop with an error message. Last term in ju is used 
      !           to deactivate jw3 work array. The space avaliable for the
      !           factorization is (iwk-1).
      !
      ! On return:
      !===========
      !
      ! U       = matrix stored in Modified Symmetric Sparse Row (MSSR) format
      !           containing the U factor. 
      !           U(1:n) contains the squared of the diagonal of the factor U.
      !           Extra-Diagonal elements are stored in
      !           U(n+2:nterm) and they scaled by the inverse of the diagonal of U
      !           M = D+T with D=diag(U)^2 T=D(U)^{-1} U
      !           U(n+1) is unused.
      !
      ! ju      = integer array. The first n+1 position contain the pointers to
      !           the beginning of each row of U after the diagonal. ju(n+2:*)
      !           are the column indeces of the matrix factor U (excluding the
      !           diagonal entry). [ju(n+1)-1] is the dimension of vectors U and
      !           ju.
      !
      ! ierr    = integer. Error message with the following meaning.
      !           ierr  = 0    --> successful return.
      !           ierr  = -1   --> Error. input matrix may be wrong.
      !                            (The elimination process has generated a
      !                            row in U whose length is .gt.  n).
      !           ierr  = -2   --> The matrix U overflows the array U.
      !           ierr  = -3   --> Illegal value for lfil.
      !           ierr  = -4   --> zero row encountered.
      !           ierr  = -5   --> Non Positive Diagonal Element
      !
      ! work arrays:
      !=============
      ! iw1     = integer work array of lenght nterm+1
      ! jw1     = integer work array of length n (max number of non-zero).
      ! jw2     = integer work array of length n.
      ! jw3     = integer work array of length n.
      ! w       = real work array of length n+1 (max number of non-zero).
      !  
      !----------------------------------------------------------------------
      ! w, jw1      store the working array [1:ii-1 = L-part, ii:n = u] 
      ! jw2         stores nonzero indicators
      ! jw3         stores pointers to the "interesting" part of row of U 
      !             factor
      ! iw1         stores lines to be eliminated for each line.
      !             first n positions are lenght of each line of iw1
      !                                 "POINTERS???"
      ! 
      ! Notes:
      ! -----
      ! The diagonal elements of the input matrix must be  nonzero (at least
      ! 'structurally'). 
      !
      !----------------------------------------------------------------------* 
      !---- Dual drop strategy works as follows.                             *
      !                                                                      *
      !     1) Theresholding in L and U as set by droptol. Any element whose *
      !        magnitude is less than some tolerance (relative to the abs    *
      !        value of diagonal element in u) is dropped.                   *
      !                                                                      *
      !     2) Keeping only the largest lfil elements in the i-th row of U   * 
      !        (excluding diagonal elements).                                *
      !                                                                      *
      ! Flexibility: one  can use  droptol=0  to get  a strategy  based on   *
      ! keeping  the largest  elements in  each row  of U. Taking            *
      ! droptol .ne.  0 but lfil=n will give  the usual threshold strategy   *
      ! (however, fill-in is then umpredictible).                            *
      !----------------------------------------------------------------------*
      !
      !  Locals
      integer  ju0, rowHead, rowMid , rowEnd, jrow , jcol, jpos
      integer  lenu , lenu0 , lenl , lenght, ncut
      integer  i,  ii, j,  jj, k, kk, j1, j2 , j3 , j4
      !      
      real*8   tnorm, abstol, s, fact, zero 
      !      
      parameter(zero = 0.0)
      !
      !-----------------------------------------------------------------------
      !  Initialize ju0 (points to next element to be added to U,ju)
      !  and pointer array.
      !-----------------------------------------------------------------------
      ju0 = n+2 
      ju(1) = ju0
      jw3(1) = ju0 
      !
      !  Initialize nonzero indicator array. 
      !
      do j=1,n 
         jw2(j)  = 0
      enddo
      !
      !  Set ju(iwk) to deactivate jw3 pointers
      !
      ju(iwk) = -1      
      !--------------------------------------------------------------
      !  Beginning of main loop.
      !-------------------------------------------------------------
      jw1=0

      do ii = 1, n
         j1 = ia(ii)
         j2 = ia(ii+1) - 1
         j3 = ia_work(ii)
         j4 = ia_work(ii+1)-1

         !  Calculate the norm of the ii-th row         
         tnorm = 0.0d0
         do  k=j1,j2
            tnorm = tnorm+abs(a(k))
         end do
         !         
         if (tnorm .eq. 0.0) goto 999
         tnorm = tnorm/real(j2-j1+1)
         !     
         !  Unpack Upper and Lower part of row of A in array w 
         !   
         lenu  = j2 - j1 + 1
         lenu0 = lenu
         lenl  = j4 - j3 + 1         
         k = ii-1          
         do j = j1,j2         
            w(j-j1+ii)   = a(j)
            jw1(j-j1+ii) = ja(j)
            k = k + 1         
            jw2(ja(j))  = k
         end do
          

         k = 0         
         do j = j3,j4   
            ! in case of csr iw1=ja
            jw1(j-j3+1)  = iw1(j)
            k = k + 1         
            jw2(iw1(j))  = k
         end do
!!$         write(*,*) ii
!!$         write(*,*) 'range',j1,j2, j3,j4
!!$         write(*,'(10(1x,e8.1))') (w(k),k=1,n)  
!!$         write(*,'(10I6)') (jw1(k),k=1,n) 
!!$         write(*,'(10I6)') (jw2(k),k=1,n) 
!!$         write(*,*) ' ' 
!
         !----BLAS-Method-----------------------------------------------------
         !  Unpack Upper and Lower part of row of A in array w 
         !         lenu  = j2 - j1 + 1
         !         lenu0 = lenu
         !         lenl  = j4 - j3 + 1
         !         call DCOPY(lenu,a(j1),1,w(ii),1)
         !         call SCOPY(lenu,ja(j1),1,jw1(ii),1)
         !         call SCOPY(lenl,iw1(j3),1,jw1(1),1)
         !   
         !  Set non-zero indicators
         !         
         !         k = ii-1
         !         do j = j1,j2
         !            k = k + 1
         !            jw2(ja(j)) = k
         !         enddo
         !         k = 0
         !         do j = j3,j4
         !            k = k + 1
         !            jw2(iw1(j)) = k
         !         enddo
         !----End of BLAS-Method----------------------------------------------
         !
         !  Eliminate previous rows
         !
         jj = 0 
150      jj = jj+1
         if (jj .gt. lenl) goto 160
         !
         !  Determine the smallest column index
         !         
         jrow = jw1(jj)
         k = jj
         do j = jj+1,lenl
            if (jw1(j) .lt. jrow) then
               jrow = jw1(j)
               k = j
            endif
         enddo

         if (k .ne. jj) then
            ! Exchange in jw1
            j = jw1(jj)
            jw1(jj) = jw1(k)
            jw1(k) = j
            ! Exchange in jw1
            jw2(jrow) = jj
            jw2(j) = k
         endif
         !write(*,*) 'ii, jrow',ii, jrow

         !
         !  Zero out element in row by setting jw2(jrow) to zero
         !
         jw2(jrow) = 0
         !         
         !  Get the leading term of the row if it exists a(jj,ii)
         !  
         rowHead = ju(jrow)
         rowMid  = jw3(jrow)
         if (ju(rowMid) .eq. ii) then  
            !------------------------------------------------------------------
            !  Mettere la possibilità di saltare 
            !  questo ciclo se il termine pivotale
            !  è troppo piccolo: fact = U(rowMid)
            !-----------------------------------------------------------------
            !
            ! Exists, mark fill-in terms of the current row before the
            ! diagonal            
            !
            do k = rowHead,rowMid-1
               jcol = ju(k)
               jpos = jw2(jcol)
               if (jpos .eq. 0) then
                  ! This is a fill-in term
                  lenl = lenl+1
                  if (lenl .gt. n) goto 996
                  jw1(lenl) = jcol
                  jw2(jcol) = lenl
               endif
            enddo
            !                    
            ! Combine current row ii with row jrow
            !

            fact = U(rowMid)*U(jrow)
            ! fact = U(rowMid) * U(jrow) * sqrt(U(jrow))
            rowEnd = ju(jrow+1) - 1
            do k = rowMid,rowEnd
               !
               !write(*,*) ii, rowMid,ju(rowMid), jrow,  ju(k), k
               s = fact * U(k)
               jcol = ju(k)
               jpos = jw2(jcol)
               if (jpos .eq. 0) then
                  !                          
                  ! This is a fill-in element
                  !
                  lenu = lenu + 1
                  if (lenu .gt. n) goto 996
                  i = ii+lenu-1
                  jw1(i) = jcol
                  jw2(jcol) = i
                  w(i) = -s
                  !                     
               else
                  !
                  ! This is not a fill-in element
                  !
                  w(jpos) = w(jpos) - s
                  !                     
               endif
               !                  
            enddo
            !               
            ! Update pointer jw3
            !
            jw3(jrow) = jw3(jrow) + 1
            if (jw3(jrow) .gt. rowEnd) then
               ! Deactivate pointer
               jw3(jrow) = iwk
            endif
            !               
         endif
         !           
         goto 150
160      continue
         !
         !  Reset non-zero indicators
         !
         do k = 1,lenu
            jw2(jw1(ii+k-1)) = 0
         enddo
         !
         !  Store the diagonal term
         !
         s = w(ii)         
         if (s .le. zero) then
            write(iout,100) ii,s
            write(6,100) ii,s
            write(6,*) 'ilu error'
            !stop
            ierr = -5  
            write(iout,101) tnorm
            return
         endif
         U(ii) = s
         s = 1.d0/s
         !
         !  Apply dropping strategy
         !
         abstol = s * tnorm * droptol  ! Set value of true tolerance
         lenght = 0 
         do k = ii+1,ii+lenu-1
            if (abs(w(k)) .gt. abstol) then
               lenght = lenght+1
               w(ii+lenght)   = w(k)
               jw1(ii+lenght) = jw1(k)
            endif
         enddo
         !
         !  Set lenu = number of elements bigger than tolerance        
         !
         lenu = lenght
         ncut = min0(lenu,lenu0+lfil-1)
         !----------------------------------------------------------------
         !  May be interesting to save a fixed number of non-zero 
         !  for each line
         ! ncut = min0(lenu,lfil)         
         !----------------------------------------------------------------
         !
         !  Select the ncut biggest elements only if lenu gt zero
         !
         if (lenu .gt. 0) then
            jpos = ii+1
            call qsplit(w(jpos),jw1(jpos),lenu,ncut)
            !          
            !  Order them in increasing jw1
            !
            !call HPSORT(jw1(jpos),w(jpos),ncut)
            call sort_row(ncut,&
                 jw1(jpos:jpos+ncut-1),&
                 w(jpos:jpos+ncut-1))
            !
            !  Scale and store in the factor
            !         
            if ((ju0+ncut-1) .gt. iwk-1) goto 997
            do k = 0,ncut-1
               ! EF
               U(ju0+k)  = s*w(jpos+k)
               ju(ju0+k) = jw1(jpos+k) 
            enddo
            !            
            !----BLAS-Method-------------------------------------------------
            !  Scale by the pivotal term and store
            !            if ((ju0+ncut-1) .gt. iwk-1) goto 997
            !            call DCOPY(ncut,w(jpos),1,U(ju0),1)
            !            call DSCAL(ncut,s,w(jpos),1)
            !            call SCOPY(ncut,jw1(jpos),1,ju(ju0),1)
            !----End of BLAS-Method------------------------------------------
            !            
         endif
         !
         !  Update pointer to beginning of next row of U    
         !
         ju0 = ju0 + ncut
         ju(ii+1) = ju0 
         jw3(ii+1) = ju0 
         !
      enddo


      !-----------------------------------------------------------------------
      !     end main loop
      !-----------------------------------------------------------------------
      !
      !     Successful run
      !
      ierr = 0
      !     
      return
      !
      !     incomprehensible error. Matrix must be wrong.
      !     
996   ierr = -1
      return
      !    
      !     insufficient storage in U.
      !     
997   ierr = -2
      return
      !     
      !     illegal lfil entered.
      !     
998   ierr = -3
      return
      !     
      !     zero row encountered
      !     
999   ierr = -4
      return
      !----------------end-of-ilut--------------------------------------------
100   format(' NULL OR NEGATIVE DIAGONAL ELEMENT: I,J =',I8,2X,E16.5)
101   format(' NORM OF ROW =  ',E16.8)
      !-----------------------------------------------------------------------
    end subroutine my_ssr_incomplete_cheolesky

    subroutine my_unpackja(n,nterm,ia,ja,point,iw,first_row)
      !----------------------------------------------------------------------  
      !  Subroutine to unpack array ja and create array iw = "transposed"
      !  of array ja      
      !-----------------------------------------------------------------------
      implicit none
      !Input variables
      integer n,nterm
      integer ia(n+1),ja(nterm),point(n)
      !Output variables
      integer iw(nterm+1)      
      integer first_row(n+1) 
      !Local variables
      integer i,j,k,mm,nn   ,m   
      !
      !
      !Initialise pointers
      !
      do i = 1,n
         point(i) = 0
      enddo
      !
      ! Evaluate lenght of pre-diagonal terms for
      ! each line of array iw
      !
      do i = 1,n
         mm = ia(i)+1
         nn = ia(i+1)-1
         do j = mm,nn
            point(ja(j)) = point(ja(j)) + 1
         enddo
      enddo
      
      !
      ! Set pointers first_row
      !    
      first_row(1) = 1
      do i = 2,n+1
         !write(*,*) i-1, point(i-1)
         first_row(i) = first_row(i-1) + point(i-1)
         point(i-1)   = first_row(i-1)-1  
         !write(*,*) i, first_row(i)         
      enddo

      
      !
      !  Set iw1
      !
      m=0
      iw=0
      do i = 1,n
         mm = ia(i) + 1
         nn = ia(i+1)-1
         do j = mm,nn
            k = ja(j)
            point(k) = point(k) + 1            
            iw( point(k) ) = i
         enddo
      enddo
      
!!$      iw(1)=iw(2)-1
!!$      do i = 1,n
!!$         mm = first_row(i)
!!$         nn = first_row(i+1)-1
!!$         !write(*,*) 'range,',   mm, nn 
!!$         do j = mm,nn
!!$            write(*,*) i,iw(j)                    
!!$         enddo
!!$      enddo

      !
      return
      !
    end subroutine my_unpackja

    !----------------------------------------------------------------------- 
    subroutine qsplit(a,ind,n,ncut)
      use Globals
      integer,           intent(in   ):: n
      integer,           intent(inout):: ind(n)
      integer,           intent(in   ):: ncut
      real(kind=double), intent(inout) :: a(n)
      !-----------------------------------------------------------------------
      !     does a quick-sort split of a real array.
      !     on input a(1:n). is a real array
      !     on output a(1:n) is permuted such that its elements satisfy:
      !
      !     abs(a(i)) .ge. abs(a(ncut)) for i .lt. ncut and
      !     abs(a(i)) .le. abs(a(ncut)) for i .gt. ncut
      !
      !     ind(1:n) is an integer array which permuted in the same way as a(*).
      !-----------------------------------------------------------------------
      !local
      real(kind=double) :: tmp, abskey
      integer itmp, first, last
      integer i,j,mid
      logical keep_cycling
      !-----
      first = 1
      last = n
      if (ncut .lt. first .or. ncut .gt. last) return
      !
      !     outer loop -- while mid .ne. ncut do
      !

      keep_cycling = .true.
      do while ( keep_cycling ) 
         mid = first
         abskey = abs(a(mid))
         do j=first+1, last
            if (abs(a(j)) .gt. abskey) then
               mid = mid+1
               !     interchange
               tmp = a(mid)
               itmp = ind(mid)
               a(mid) = a(j)
               ind(mid) = ind(j)
               a(j)  = tmp
               ind(j) = itmp
            endif
         end do
         !
         !     interchange
         !
         tmp = a(mid)
         a(mid) = a(first)
         a(first)  = tmp
         !
         itmp = ind(mid)
         ind(mid) = ind(first)
         ind(first) = itmp
         !
         !     test for while loop
         !
         if (mid .eq. ncut) keep_cycling = .false.
         if (mid .gt. ncut) then
            last = mid-1
         else
            first = mid+1
         end if
      end do
    end subroutine qsplit

    subroutine  mss2scaled_upper(lun_err,&
         nrow,&
         iwk,&
         ju,&
         U,&
         upper)
      use Globals
      !-------------------------------------------------------------------
      ! Change the format of a given matrix from Modified Symmetric Sparse
      ! Row format to CSR. The upper matrix is 
      !             upper = diag(U)^{-1} U
      ! thus its has one on the diagonal
      !-------------------------------------------------------------------
      implicit none
      ! Input variables      
      integer,           intent(in   ) :: lun_err
      integer,           intent(in   ) :: nrow
      integer,           intent(in   ) :: iwk
      integer,           intent(in   ) :: ju(iwk)
      real(kind=double), intent(in   ) :: U(iwk)
      ! Output variables
      type(spmat),       intent(inout) :: upper
      !  Local variables
      integer i,j,k
      integer ind,lenght,nterm
      real(kind=double) :: d_term


      !
      ! init upper factor
      !
      if ( upper%is_initialized ) then
         call upper%kill(lun_err)
      end if
      nterm = ju(nrow+1)-1
      call upper%init(lun_err, &
           nrow, nrow, nterm,&
           'csr',&
           dim_kernel=0,&
           is_symmetric=.false.,&
           unitary_diag=.true.,&
           triangular='U')

      !
      ! copy into the csr the matrix diag(U)^{-1} U
      ! 
      upper%ia(1) = 1
      do i = 1,nrow
         ! diagonal
         k      = ju(i)
         lenght = ju(i+1)-k
         ind           = upper%ia(i)
         upper%ia(i+1) = ind+lenght+1
         upper%coeff(ind) = one
         upper%ja(ind)    = i
         ! extra diagonal
         ind = ind+1
         do j = 0,lenght-1
            upper%coeff(ind+j)  = U(k+j)
            upper%ja(ind+j)     = ju(k+j)
         end do
      end do

    end subroutine mss2scaled_upper


    !
    ! incomplete Choleski decompostion with no-fillin
    !
    subroutine kersh_loc(iout,nequ,nterm,ia,ja,sysmat,prec)
      use Globals
      implicit none
      integer  iout,nequ,nterm
      integer  ia(nequ+1),ja(nterm)
      integer  i,j,k,kk,k1,i1,j1,k2
      integer :: nterm_prec
      real(kind=double) :: prec(nterm),sysmat(nterm)
      real(kind=double) :: a

      info=0
      !
      do k=1,nterm
         prec(k) = zero
      end do
      !
      do kk=1,nequ-1
         k = ia(kk)
         a = sysmat(k) - prec(k)
         if (a.le.zero) then
            info = info+1
            write(iout,100) kk,a
            write(iout,101) prec(ia(kk-1))
            a = (prec(ia(kk-1)))**2
         end if
         prec(k) = sqrt(a)
         !
         i = ia(kk) + 1
         j = ia(kk+1) - 1
         !
         do k1 = i,j
            prec(k1) = (sysmat(k1) - prec(k1)) / prec(k)
         end do
         !
         do k2 = i,j-1
            j1 = ia(ja(k2))
            prec(j1) = prec(j1) + prec(k2)**2
            i1 = k2 + 1
            j1 = j1 + 1
            do while (j1.lt.ia(ja(k2)+1).and.i1.le.j)
               if (ja(j1).eq.ja(i1)) then
                  prec(j1) = prec(j1) + prec(k2) * prec(i1)
                  i1 = i1 + 1
                  j1 = j1 + 1
               else if (ja(j1).lt.ja(i1)) then
                  j1 = j1 + 1
               else if (ja(j1).gt.ja(i1)) then
                  i1 = i1 + 1
               end if
            end do
         end do
         !
         if (j.ge.i) prec(ia(ja(j))) = prec(ia(ja(j))) + prec(k2)**2
      end do
      !
      k = ia(nequ) 
      a = sysmat(k) - prec(k)
      if (a.le.zero) then
         write(iout,100) nequ ,a
         write(iout,101) prec(ia(nequ-1))
         a = (prec(ia(nequ-1)))**2
      end if
      prec(k) = sqrt(a)
      !
      return
100   format(' ELEMENTO DIAGONALE DI L NULLO,I,J =',I5,2X,E16.5)
101   format(' ELEMENTO DIAGONALE PRECEDENTE =  ',E16.8)
    end subroutine kersh_loc


  end subroutine incomplete_cholesky



  !-------------------------------------------------------
  ! Subroutine for incomplete LU factorization
  !-------------------------------------------------------
  subroutine incomplete_lu(this,&
       lun_err,n_fillin,tol_fillin,info,&
       lower, upper)
    use Globals
    implicit none
    class(spmat),               intent(in   ) :: this
    integer,                    intent(in   ) :: lun_err
    integer,                    intent(in   ) :: n_fillin
    real(kind=double),          intent(in   ) :: tol_fillin
    integer,                    intent(inout) :: info
    type(spmat),                intent(inout) :: lower
    type(spmat),                intent(inout) :: upper

    !local
    logical :: rc
    integer :: res, i, j
    integer :: nterm_max, nterm 
    integer :: nrow
    integer,           allocatable :: jlu(:),ju(:),jw(:),iscr(:)
    integer,           allocatable :: iw1(:),jw1(:),jw2(:),jw3(:)
    integer,           allocatable :: idiag(:),iwork(:)
    real(kind=double), allocatable :: alu(:),w(:),rwork(:)
    type(spmat)  :: lower_loc


    nrow = this%nrow
    nterm = this%nterm
    
    if ( this%storage_system .eq. 'ssr') then
       rc = IOerr(lun_err, err_inp, 'incomplete_lu', &
            'storage system ssr not supported yet')
    end if
    

    if ( n_fillin .eq. 0 ) then
       !
       ! use ILU0 from Saad
       !

       !
       ! prepare work arrays and diagonal pointer
       !
       allocate (iwork(nrow), idiag(nrow),rwork(nterm),stat=res)
       if( res .ne. 0) rc = IOerr(lun_err, err_alloc, 'incomplete_lu', &
            'work arrays idiag iwork rwork' )

       !
       ! build LU decomposition in Modified Compressed Storage
       !
       call loc_ilu0(lun_err,&
            nrow,nterm,info,this%ia,this%ja,idiag,iwork,&
            this%coeff,rwork)
!!$       do i=1,nrow
!!$          write(60,*) i,this%ja(idiag(i))
!!$          do j=this%ia(i),this%ia(i+1)-1
!!$             write(60,*) i,this%ja(j),rwork(j)
!!$          end do
!!$       end do
       
       !
       ! convert to scr format
       ! init. lower and upper matrix
       !
       call ilu02lu(lun_err,nrow,nterm,&
            rwork,this%ia,this%ja,idiag,&
            lower, upper)

       !
       ! prepare work arrays and diagonal pointer
       !
       deallocate (iwork, idiag,rwork,stat=res)
       if( res .ne. 0)  rc = IOerr(lun_err, err_dealloc, 'ilu0', &
            'work arrays idiag iwork rwork' )
    else
       !
       ! use ILUT from Saad
       !
       ! prepare work arrays 
       nrow      = this%nrow
       nterm_max = 2 * (n_fillin+1) * nrow

              
       allocate(&
            alu(nterm_max), &
            jlu(nterm_max),&
            ju(nrow),&
            w(nrow+1),&
            jw(2*nrow),&
            iscr(nrow+1),&
            stat=res)
       if (res.ne.0) rc = IOerr(lun_err, err_alloc, 'incomplete_lu', &
            'work arrays alu, jlu ju w jw',res)

       !
       ! non-symmetric case 
       !
       info =0
       call loc_ilut(nrow,&
            this%coeff,&
            this%ja,&
            this%ia,&
            n_fillin,&
            tol_fillin,&
            alu,jlu,&
            ju,&
            nterm_max,&
            w,jw,info)

       if ( info .eq. 0 ) then
          !
          ! convert msr into csr format
          ! init. lower and upper matrix
          !    
          call ilut2lu(lun_err,nrow,nterm_max,&
               alu,jlu,ju,&! jw, &
               lower, upper)
       else
          !
          ! info in case of lu factorization error 
          !
          write(lun_err,*) ' Error in ILUT construction' 
          select case (info) 
          case (-1)
             rc = IOerr(lun_err, err_val, 'build_ic_fillin', &
                  ' The elimination process has generated'//&
                  ' a row in U whose length is .gt.  n)')
          case (-2)  
             rc = IOerr(lun_err, err_val, 'build_ic_fillin', &
                  'The matrix U overflows the array U')
          case (-3)  
             rc = IOerr(lun_err, err_val, 'build_ic_fillin', &
                  ' Illegal value for n_fillin')
          case(-4)   
             rc = IOerr(lun_err, err_val, 'build_ic_fillin', &
                  'zero row encountered')
          case(-5) 
             rc = IOerr(lun_err, err_val, 'build_ic_fillin', &
                  'Non Positive Diagonal Element')
          end select
          if (info .gt. 0) then
             rc = IOerr(lun_err, err_val, 'build_ic_fillin', &
                  'Zero pivot encountered at step number=',info)
          end if
          return
       end if


       !
       ! free memory
       !
       if ( lower_loc%is_initialized ) call lower_loc%kill(lun_err)

       deallocate(&
            alu, &
            jlu,&
            ju,&
            w,&
            jw,&
            stat=res)
       if (res.ne.0) rc = IOerr(lun_err, err_dealloc, 'incomplete_lu', &
            'work arrays alu, jlu ju w jw',res)  
    end if

  contains
    !----------------------------------------------------------------------- 
    subroutine qsplit(a,ind,n,ncut)
      use Globals
      integer,           intent(in   ):: n
      integer,           intent(inout):: ind(n)
      integer,           intent(in   ):: ncut
      real(kind=double), intent(inout) :: a(n)
      !-----------------------------------------------------------------------
      !     does a quick-sort split of a real array.
      !     on input a(1:n). is a real array
      !     on output a(1:n) is permuted such that its elements satisfy:
      !
      !     abs(a(i)) .ge. abs(a(ncut)) for i .lt. ncut and
      !     abs(a(i)) .le. abs(a(ncut)) for i .gt. ncut
      !
      !     ind(1:n) is an integer array which permuted in the same way as a(*).
      !-----------------------------------------------------------------------
      !local
      real(kind=double) :: tmp, abskey
      integer itmp, first, last
      integer i,j,mid
      logical keep_cycling
      !-----
      first = 1
      last = n
      if (ncut .lt. first .or. ncut .gt. last) return
      !
      !     outer loop -- while mid .ne. ncut do
      !

      keep_cycling = .true.
      do while ( keep_cycling ) 
         mid = first
         abskey = abs(a(mid))
         do j=first+1, last
            if (abs(a(j)) .gt. abskey) then
               mid = mid+1
               !     interchange
               tmp = a(mid)
               itmp = ind(mid)
               a(mid) = a(j)
               ind(mid) = ind(j)
               a(j)  = tmp
               ind(j) = itmp
            endif
         end do
         !
         !     interchange
         !
         tmp = a(mid)
         a(mid) = a(first)
         a(first)  = tmp
         !
         itmp = ind(mid)
         ind(mid) = ind(first)
         ind(first) = itmp
         !
         !     test for while loop
         !
         if (mid .eq. ncut) keep_cycling = .false.
         if (mid .gt. ncut) then
            last = mid-1
         else
            first = mid+1
         end if
      end do
    end subroutine qsplit

    !
    !************************* ILU0 *************************************
    !
    !>--------------------------------------------------------------
    !> Subroutine for the computation of the lu decompostion
    !> in the mss format
    !> DIAGONAL IS NOT INVERTED
    !>--------------------------------------------------------------
    subroutine loc_ilu0(lun_err,&
         nequ,nterm,info,ia,ja,idiag,iwork,&
         sysmat,prec)
      !
      !  calculates the ILU(0) factorization of matrix SYSMAT
      !  (algorithm 10.4 page 277 in Saad 1996)
      !
      use Globals
      implicit none
      integer :: lun_err,nequ,nterm,info   
      integer :: irow,jrow,ind1,ind2,j,jnd,jw
      integer :: ia(nequ+1),ja(nterm),idiag(nequ),iwork(nequ)
      real(kind=double) :: sysmat(nterm),prec(nterm)
      real(kind=double) :: mult


      info = 0 
      !
      !  initializes prec to sysmat
      !
      call dcopy(nterm,sysmat,1,prec,1)

      do j=1,nequ
         iwork(j)=0
      end do

      do irow=1,nequ

         ind1 = ia(irow)
         ind2 = ia(irow+1)-1

         do j=ind1,ind2
            iwork(ja(j)) = j
         end do

         jnd=ind1
         jrow = ja(jnd)
         do while (jnd.le.ind2 .and. jrow.lt.irow)
            
            mult = prec(jnd)*prec(idiag(jrow))
            prec(jnd) = mult

            do j=idiag(jrow)+1,ia(jrow+1)-1
               jw=iwork(ja(j))
               if(jw.ne.0) prec(jw)=prec(jw)-mult*prec(j)
            end do

            jnd=jnd+1
            jrow = ja(jnd)

         end do

         !
         !  set pointer to diagonal element
         !
         idiag(irow)=jnd
         !
         !   check if diagonal element exists and MULT is nonzero 
         !
         if(jrow.ne.irow .or. prec(jnd).eq.zero) then
            write(lun_err,*)' ILU(0) error at row',irow
            info=irow
            return
         end if

         prec(jnd)=one/prec(jnd)

         do j=ind1,ind2
            iwork(ja(j))=0
         end do

      end do

      !
      ! invert diagonal back
      !
      do j=1,nequ
         prec(idiag(j)) = one / prec(idiag(j))
      end do

    end subroutine loc_ilu0

    !>---------------------------------------------------
    !> Converts the lud factors from the mss format to 
    !> two matrix (upper, lower)  in csr format 
    !>       A ~ L U
    !> L having one on the diagonal
    !>---------------------------------------------------
    subroutine ilu02lu(lun_err,nrow,nterm,&
         a_lu,ia_lu,ja_lu,idiag,&
         lower, upper)
      use Globals
      implicit none

      integer,           intent(in   ) :: lun_err
      integer,           intent(in   ) :: nrow
      integer,           intent(in   ) :: nterm
      real(kind=double), intent(in   ) :: a_lu(nterm)
      integer,           intent(in   ) :: ia_lu(nrow+1)
      integer,           intent(in   ) :: ja_lu(nterm)
      integer,           intent(in   ) :: idiag(nrow)
      type(spmat),       intent(inout) :: lower, upper
      !local
      logical :: rc
      integer :: res 
      integer :: irow, j,start,i
      integer :: low_ind,up_ind
      integer, allocatable :: lower_ia(:),upper_ia(:)
      integer, allocatable :: lower_ja(:),upper_ja(:)
      real(kind=double), allocatable :: lower_coeff(:),upper_coeff(:)


      !
      ! allocation of work arrays 
      !
      allocate(&
           lower_ia(nrow+1),&
           upper_ia(nrow+1),&
           lower_ja(nterm),&
           upper_ja(nterm),&
           lower_coeff(nterm),&
           upper_coeff(nterm),&
           stat=res)
      if (res.ne.0) rc = IOerr(lun_err, err_alloc, 'msr2ldu', &
           'work arrays lower_ia, lower_ja, lower_coeff'//&
           ' upper_ia, upper_ja, upper_coeff',res) 

      !
      ! init counter
      !
      low_ind= 0
      up_ind = 0
      lower_ia(1) = 1
      upper_ia(1) = 1

      do irow=1,nrow     
         !
         ! copy lower matrix
         !
         do j = ia_lu(irow), idiag(irow)-1 
            low_ind              = low_ind+1
            lower_ja(low_ind)    = ja_lu(j)
            lower_coeff(low_ind) = a_lu(j)
         end do
         ! add one on the diagonal
         low_ind              = low_ind+1
         lower_ja(low_ind)    = irow
         lower_coeff(low_ind) = one
         ! set next index of first element 
         lower_ia(irow+1) = low_ind+1
         !write(*,*) 'irow, lower_ia(irow+1)',irow, lower_ia(irow+1)
         !write(*,*) lower_coeff(lower_ia(irow):lower_ia(irow+1)-1)

         !
         ! copy upper matrix
         !
         ! the diagonal term
         up_ind              = up_ind+1
         upper_ja(up_ind)    = irow
         upper_coeff(up_ind) = a_lu(idiag(irow))
         ! extra diagonal terms
         do j = idiag(irow)+1, ia_lu(irow+1)-1
            up_ind              = up_ind+1
            upper_ja(up_ind)    = ja_lu(j)
            upper_coeff(up_ind) = a_lu(j)
         end do
         ! set next index of first element
         upper_ia(irow+1) = up_ind+1
      end  do

      !
      ! init and copy lower sparse matrix
      !
      if (lower%is_initialized) call lower%kill(lun_err)
      call lower%init(lun_err, &
           nrow, nrow, low_ind,&
           'csr',&
           dim_kernel=0,&
           is_symmetric=.false.,&
           unitary_diag = .true.,&
           triangular='L')

      lower%ia               = lower_ia
      lower%ja   (1:low_ind) = lower_ja   (1:low_ind)
      lower%coeff(1:low_ind) = lower_coeff(1:low_ind)
      lower%is_sorted = .true.
      call lower%sort()


      !
      ! init and copy upper sparse matrix
      !
      if (upper%is_initialized) call  upper%kill(lun_err)
      !
      call upper%init(lun_err, &
           nrow, nrow, up_ind,&
           'csr',&
           dim_kernel=0,&
           is_symmetric=.false.,&
           triangular='U')

      upper%ia              = upper_ia
      upper%ja(1:up_ind)    = upper_ja(1:up_ind)
      upper%coeff(1:up_ind) = upper_coeff(1:up_ind)
      upper%is_sorted = .true.
      call upper%sort()

      !
      ! free memory 
      !
      deallocate(&
           lower_ia,&
           upper_ia,&
           lower_ja,&
           upper_ja,&
           lower_coeff,&
           upper_coeff,&
           stat=res)
      if (res.ne.0) rc = IOerr(lun_err, err_dealloc, 'msr2ldu', &
           'work arrays lower_ia, lower_ja, lower_coeff'//&
           ' upper_ia, upper_ja, upper_coeff',res) 

    end subroutine ilu02lu


    
    !---------------------------------------------------------
    !                          S P A R S K I T  
    !---------------------------------------------------------
    !                   ITERATIVE SOLVERS MODULE 
    !---------------------------------------------------------
    ! This Version Dated: August 13, 1996. 
    !  Warning: meaning of some  arguments have changed 
    ! w.r.t. earlier versions. 
    ! Some Calling sequences may also have changed
    !--------------------------------------------------------
    subroutine loc_ilut(n,a,ja,ia,lfil,droptol,alu,jlu,ju,iwk,w,jw,ierr)
      !--------------------------------------------------------------------
      implicit none 
      integer n 
      real*8 a(*),alu(*),w(n+1),droptol
      integer ja(*),ia(n+1),jlu(*),ju(n),jw(2*n),lfil,iwk,ierr
      !--------------------------------------------------------------------*
      !                      *** ILUT preconditioner ***                   *
      !      incomplete LU factorization with dual truncation mechanism    *
      !--------------------------------------------------------------------*
      !     Author: Yousef Saad *May, 5, 1990, Latest revision, August 1996*
      !--------------------------------------------------------------------*
      ! PARAMETERS                                                           
      !-----------                                                           
      !
      ! on entry:
      !========== 
      ! n       = integer. The row dimension of the matrix A. The matrix 
      !
      ! a,ja,ia = matrix stored in Compressed Sparse Row format.           
      !
      ! lfil    = integer. The fill-in parameter. Each row of L and each row
      !           of U will have a maximum of lfil elements (excluding the 
      !           diagonal element). lfil must be .ge. 0.
      !           ** WARNING: THE MEANING OF LFIL HAS CHANGED WITH RESPECT TO
      !           EARLIER VERSIONS. 
      !
      ! droptol = real*8. Sets the threshold for dropping small terms in the
      !           factorization. See below for details on dropping strategy.
      !
      !  
      ! iwk     = integer. The lengths of arrays alu and jlu. If the arrays
      !           are not big enough to store the ILU factorizations,i lut
      !           will stop with an error message. 
      !
      ! On return:
      !===========
      !
      ! alu,jlu = matrix stored in Modified Sparse Row (MSR)
      !           format containing the L and U factors together.
      !           The diagonal (stored in alu(1:n) ) is not inverted.
      !           Each i-th row of the alu,jlu matrix
      !           contains the i-th row of L (excluding the diagonal entry=1)
      !           followed by the i-th row of U.
      !
      ! ju      = integer array of length n containing the pointers to
      !           the beginning of each row of U in the matrix alu,jlu.
      !
      ! ierr    = integer. Error message with the following meaning.
      !           ierr  = 0   --> successful return.
      !           ierr .gt. 0 --> zero pivot encountered at step number ierr.
      !           ierr  = -1  --> Error. input matrix may be wrong.
      !                            (The elimination process has generated a
      !                            row in L or U whose length is .gt.  n.)
      !           ierr  = -2  --> The matrix L overflows the array al.
      !           ierr  = -3  --> The matrix U overflows the array alu.
      !           ierr  = -4  --> Illegal value for lfil.
      !           ierr  = -5  --> zero row encountered.
      !
      ! work arrays:
      !=============
      ! jw      = integer work array of length 2*n.
      ! w       = real work array of length n 
      !  
      !----------------------------------------------------------------------
      ! w, ju (1:n) store the working array [1:ii-1 = L-part, ii:n = u] 
      ! jw(n+1:2n)  stores nonzero indicators
      ! 
      ! Notes:
      ! ------
      ! The diagonal elements of the input matrix must be  nonzero (at least
      ! 'structurally'). 
      !
      !--------------------------------------------------------------------* 
      !---- Dual drop strategy works as follows.                           *
      !                                                                    *
      !     1) Theresholding in L and U as set by droptol. Any element     *
      !        whose magnitude is less than some tolerance (relative to    *
      !        the abs value of diagonal element in u) is dropped.         *
      !                                                                    *
      !     2) Keeping only the largest lfil elements in the i-th row of L * 
      !        and the largest lfil elements in the i-th row of U          *
      !        (excluding diagonal elements).                              *
      !                                                                    *
      ! Flexibility: one  can use  droptol=0  to get  a strategy  based on *
      ! keeping  the largest  elements in  each row  of L  and U.   Taking *
      ! droptol .ne.  0 but lfil=n will give  the usual threshold strategy *
      ! (however, fill-in is then unpredictible).                          *
      !--------------------------------------------------------------------*
      !     locals
      integer ju0,k,j1,j2,j,ii,i,lenl,lenu,jj,jrow,jpos,len 
      real*8 tnorm, t, abs, s, fact 
      if (lfil .lt. 0) goto 998
      !---------------------------------------------------------------------
      !     initialize ju0 (points to next element to be added to alu,jlu)
      !     and pointer array.
      !---------------------------------------------------------------------
      ju0 = n+2
      jlu(1) = ju0
      !
      !     initialize nonzero indicator array. 
      !
      do 1 j=1,n
         jw(n+j)  = 0
1     end do
      !---------------------------------------------------------------------
      !     beginning of main loop.
      !---------------------------------------------------------------------
      do 500 ii = 1, n
         j1 = ia(ii)
         j2 = ia(ii+1) - 1
         !write(0,*) j1,j2
         tnorm = 0.0d0
         do 501 k=j1,j2
            tnorm = tnorm+abs(a(k))
501      end do
         if (tnorm .eq. 0.0) goto 999
         tnorm = tnorm/real(j2-j1+1)
         !     
         ! unpack L-part and U-part of row of A in arrays w 
         !     
         lenu = 1
         lenl = 0
         jw(ii) = ii
         w(ii) = 0.0
         jw(n+ii) = ii
         !
         do 170  j = j1, j2
            k = ja(j)
            t = a(j)
            if (k .lt. ii) then
               lenl = lenl+1
               jw(lenl) = k
               w(lenl) = t
               jw(n+k) = lenl
            else if (k .eq. ii) then
               w(ii) = t
            else
               lenu = lenu+1
               jpos = ii+lenu-1 
               jw(jpos) = k
               w(jpos) = t
               jw(n+k) = jpos
            endif
170      end do
         jj = 0
         len = 0 
         !     
         !     eliminate previous rows
         !     
150      jj = jj+1
         if (jj .gt. lenl) goto 160
         !---------------------------------------------------
         ! in order to do the elimination in the correct 
         ! order we must select  the smallest 
         !  column index among jw(k), k=jj+1, ..., lenl.
         !---------------------------------------------------
         jrow = jw(jj)
         k = jj
         !     
         ! determine smallest column index
         !     
         do 151 j=jj+1,lenl
            if (jw(j) .lt. jrow) then
               jrow = jw(j)
               k = j
            end if
151      end do
         !
         if (k .ne. jj) then
            ! exchange in jw
            j = jw(jj)
            jw(jj) = jw(k)
            jw(k) = j
            ! exchange in jr
            jw(n+jrow) = jj
            jw(n+j) = k
            ! exchange in w
            s = w(jj)
            w(jj) = w(k)
            w(k) = s
         endif
         !
         ! zero out element in row by setting 
         ! jw(n+jrow) to zero.
         !     
         jw(n+jrow) = 0
         !
         ! get the multiplier for row to be eliminated (jrow).
         !     
         fact = w(jj)*alu(jrow)
         if (abs(fact) .le. droptol) goto 150
         !     
         ! combine current row and row jrow
         !
         do 203 k = ju(jrow), jlu(jrow+1)-1
            s = fact*alu(k)
            j = jlu(k)
            jpos = jw(n+j)
            if (j .ge. ii) then
               !     
               !     dealing with upper part.
               !     
               if (jpos .eq. 0) then
                  !
                  ! this is a fill-in element
                  !     
                  lenu = lenu+1
                  if (lenu .gt. n) goto 995
                  i = ii+lenu-1
                  jw(i) = j
                  jw(n+j) = i
                  w(i) = - s
               else
                  !
                  ! this is not a fill-in element 
                  !
                  w(jpos) = w(jpos) - s
               endif
            else
               !     
               ! dealing  with lower part.
               !     
               if (jpos .eq. 0) then
                  !
                  !     this is a fill-in element
                  !     
                  lenl = lenl+1
                  if (lenl .gt. n) goto 995
                  jw(lenl) = j
                  jw(n+j) = lenl
                  w(lenl) = - s
               else
                  !     
                  !     this is not a fill-in element 
                  !     
                  w(jpos) = w(jpos) - s
               endif
            endif
203      end do
         !     
         ! store this pivot element -- 
         ! (from left to right -- no danger of
         ! overlap with the working elements in L (pivots). 
         !     
         len = len+1 
         w(len) = fact
         jw(len)  = jrow
         goto 150

         !
         ! reset double-pointer to zero (U-part)
         !     
160      do k=1, lenu
            jw(n+jw(ii+k-1)) = 0
         end do
         !     
         ! update L-matrix
         !     
         lenl = len 
         len = min0(lenl,lfil)
         !     
         ! sort by quick-split
         !
         call qsplit (w,jw,lenl,len)
         !
         ! store L-part
         ! 
         do 204 k=1, len 
            if (ju0 .gt. iwk) goto 996
            alu(ju0) =  w(k)
            jlu(ju0) =  jw(k)
            ju0 = ju0+1
204      end do
         !     
         !save pointer to beginning of row ii of U
         !     
         ju(ii) = ju0
         !
         ! update U-matrix -- first apply dropping strategy 
         !
         len = 0
         do k=1, lenu-1
            if (abs(w(ii+k)) .gt. droptol*tnorm) then 
               len = len+1
               w(ii+len) = w(ii+k) 
               jw(ii+len) = jw(ii+k) 
            endif
         end do
         lenu = len+1
         len = min0(lenu,lfil)
         !
         call qsplit (w(ii+1), jw(ii+1), lenu-1,len)
         !
         !     copy
         ! 
         t = abs(w(ii))
         if (len + ju0 .gt. iwk) goto 997
         do 302 k=ii+1,ii+len-1 
            jlu(ju0) = jw(k)
            alu(ju0) = w(k)
            t = t + abs(w(k) )
            ju0 = ju0+1
302      end do
         !     
         !     store inverse of diagonal element of u
         !   

         if (w(ii) .eq. 0.0) w(ii) = (0.0001 + droptol)*tnorm
         alu(ii) = 1.0d0/w(ii) 
         !     
         !     update pointer to beginning of next row of U.
         !     
         jlu(ii+1) = ju0
         !--------------------------------------------------------------
         !     end main loop
         !--------------------------------------------------------------
500   end do
      ! EF
      ! get the not inverted diagonal
      do ii=1,n
         alu(ii)=one/alu(ii)
      end do

      ierr = 0
      return
      !
      !     incomprehensible error. Matrix must be wrong.
      !     
995   ierr = -1
      return
      !     
      !     insufficient storage in L.
      !     
996   ierr = -ju0
      return
      !     
      !     insufficient storage in U.
      !     
997   ierr = -len-ju0
      return
      !     
      !     illegal lfil entered.
      !     
998   ierr = -4
      return
      !     
      !     zero row encountered
      !     
999   ierr = -5
      return
      !----------------end-of-ilut----------------------
      !-------------------------------------------------

    end subroutine loc_ilut

    subroutine ilut2lu(lun_err,nrow,nterm,&
         a_lu,ja_lu,ia_u,&
         lower, upper)
      use Globals
      implicit none

      integer,           intent(in   ) :: lun_err
      integer,           intent(in   ) :: nrow
      integer,           intent(in   ) :: nterm
      real(kind=double), intent(in   ) :: a_lu(nterm)
      integer,           intent(in   ) :: ja_lu(nterm)
      integer,           intent(in   ) :: ia_u(nrow)
      type(spmat),  intent(inout) :: lower, upper
      !local
      logical :: rc
      integer :: res 
      integer :: irow, j,start,i
      integer :: low_ind,up_ind
      integer, allocatable :: lower_ia(:),upper_ia(:)
      integer, allocatable :: lower_ja(:),upper_ja(:)
      real(kind=double), allocatable :: lower_coeff(:),upper_coeff(:)


      !
      ! allocation of work arrays 
      !
      allocate(&
           lower_ia(nrow+1),&
           upper_ia(nrow+1),&
           lower_ja(nterm),&
           upper_ja(nterm),&
           lower_coeff(nterm),&
           upper_coeff(nterm),&
           stat=res)
      if (res.ne.0) rc = IOerr(lun_err, err_alloc, 'msr2ldu', &
           'work arrays lower_ia, lower_ja, lower_coeff'//&
           ' upper_ia, upper_ja, upper_coeff',res) 

      !
      ! init counter
      !
      low_ind=0
      up_ind=0
      lower_ia(1)=1
      upper_ia(1)=1

!!$    do i=1,nrow
!!$       write(*,*) i,a_lu(i),ja_lu(i)
!!$    end do

      do irow=1,nrow     
         !
         ! copy lower matrix
         !
         !write(*,*) 'lu',irow
         do j = ja_lu(irow), ia_u(irow)-1 
            low_ind              = low_ind+1
            lower_ja(low_ind)    = ja_lu(j)
            lower_coeff(low_ind) = a_lu(j)
            !write(*,*) irow, ja_lu(j), a_lu(j) 
            !write(*,*) irow, lower_ja(low_ind),lower_coeff(low_ind)
         end do
         ! add one on the diagonal
         low_ind              = low_ind+1
         lower_ja(low_ind)    = irow
         lower_coeff(low_ind) = one
         ! set next index of first elemen 
         lower_ia(irow+1) = low_ind+1
         !write(*,*) 'irow, lower_ia(irow+1)',irow, lower_ia(irow+1)
         !write(*,*) lower_coeff(lower_ia(irow):lower_ia(irow+1)-1)

         !
         ! copy upper matrix
         !
         ! add one on the diagonal
         up_ind              = up_ind+1
         upper_ja(up_ind)    = irow
         upper_coeff(up_ind) = a_lu(irow)
         do j = ia_u(irow), ja_lu(irow+1)-1
            up_ind              = up_ind+1
            upper_ja(up_ind)    = ja_lu(j)
            upper_coeff(up_ind) = a_lu(j)
         end do
         upper_ia(irow+1) = up_ind+1
      end  do

      !
      ! init and copy lower sparse matrix
      !
      if (lower%is_initialized) call lower%kill(lun_err)
      call lower%init(lun_err, &
           nrow, nrow, low_ind,&
           'csr',&
           dim_kernel=0,&
           is_symmetric=.false.,&
           unitary_diag=.true.,&
           triangular='L')
      lower%ia               = lower_ia
      lower%ja   (1:low_ind) = lower_ja   (1:low_ind)
      lower%coeff(1:low_ind) = lower_coeff(1:low_ind)
      lower%is_sorted = .false.

      call lower%sort()


      if (upper%is_initialized) call  upper%kill(lun_err)
      !
      ! init and copy upper sparse matrix
      !
      call upper%init(lun_err, &
           nrow, nrow, up_ind,&
           'csr',&
           dim_kernel=0,&
           is_symmetric=.false.,&
           triangular='U')

      upper%ia              = upper_ia
      upper%ja(1:up_ind)    = upper_ja(1:up_ind)
      upper%coeff(1:up_ind) = upper_coeff(1:up_ind)

      upper%is_sorted = .false.
      call upper%sort()

      !
      ! free memory 
      !
      deallocate(&
           lower_ia,&
           upper_ia,&
           lower_ja,&
           upper_ja,&
           lower_coeff,&
           upper_coeff,&
           stat=res)
      if (res.ne.0) rc = IOerr(lun_err, err_dealloc, 'msr2ldu', &
           'work arrays lower_ia, lower_ja, lower_coeff'//&
           ' upper_ia, upper_ja, upper_coeff',res) 

    end subroutine ilut2lu
  end subroutine incomplete_lu

  !>--------------------------------------------------------------
  !> Subroutine for sorting the arrays ja and coeff
  !> to be ordered according to column index
  !>--------------------------------------------------------------
  subroutine sort_spmat(this)
    use Globals
    implicit none
    class(spmat), intent(inout) :: this
    !local
    integer :: irow,start,finish,nnz

    if (this%is_sorted) return

    do irow = 1, this%nrow
       ! select matrix line
       ! works both for ssr an csr format
       start  = this%ia(irow)
       finish = this%ia(irow+1)-1
       nnz    = finish-start+1
       ! order 
       if (nnz .gt. 1 ) then 
          call sort_row(nnz,&
               this%ja(start:finish),&
               this%coeff(start:finish))
       end if
    end do
    this%is_sorted = .true.
  contains
    subroutine sort_matrix_line(nnz,ja,coeff)
      use Globals
      implicit none
      integer,           intent(in   ) :: nnz
      integer,           intent(inout) :: ja(nnz)
      real(kind=double), intent(inout) :: coeff(nnz)
      !local 
      integer :: i,j, indx,isgn, itemp
      real(kind=double) :: rtemp

      !  Initialize.
      i = 0
      indx = 0
      isgn = 0
      j = 0
      do 
         call global_heapsort(nnz, indx, i,j,isgn)
         if (indx .gt. 0 ) then
            ! SWAP ELEMENT 

            ! swap column indeces
            itemp = ja(i)
            ja(i) = ja(j)
            ja(j) = itemp

            ! swap real nnzero coeff
            rtemp    = coeff(i)
            coeff(i) = coeff(j)
            coeff(j) = rtemp
         else if ( indx .lt. 0) then
            ! COMPARE
            isgn = 0
            if ( ja(i) .lt.  ja(j) ) isgn = -1
            if ( ja(i) .gt.  ja(j) ) isgn = 1
         else if ( indx .eq. 0 ) then
            exit
         end if
      end do
    end subroutine sort_matrix_line
  end subroutine sort_spmat


  subroutine MtimesN(matrix_m,matrix_n, out_matrix)
    use Globals
    implicit none
    class(spmat), intent(in) :: matrix_m
    class(spmat), intent(in) :: matrix_n
    real(kind=double), intent(inout) :: out_matrix(matrix_m%nrow,matrix_n%ncol)
    !local
    integer irow,icol, k, ind,ind2
    real(kind=double) :: b,c

    out_matrix = zero
    do irow=1,matrix_m%nrow
       do ind=matrix_m%ia(irow),matrix_m%ia(irow+1)-1
          b = matrix_m%coeff(ind)
          k = matrix_m%ja(ind)
          do ind2 = matrix_n%ia(k),matrix_n%ia(k+1)-1
             c    = matrix_n%coeff(ind2)
             icol = matrix_n%ja(ind2)
             out_matrix(irow,icol) = out_matrix(irow,icol) + &
                  b * c
          end do
       end do
    end do

  end subroutine MtimesN
  
  subroutine add_row(this,lun_err,vec_w)
    use Globals
    implicit none
    ! vars
    class(spmat),      intent(inout) :: this
    integer,           intent(in   ) :: lun_err
    real(kind=double), intent(in   ) :: vec_w(this%ncol)
    ! local vars
    integer :: icol
    integer :: nrow, ncol, nterm
    type(spmat) :: matrix_out
    
    call matrix_out%init(lun_err,this%nrow+1,this%ncol+1,&
         this%nterm+this%ncol+1,this%storage_system,&
         this%dim_kernel,is_symmetric=.false.,&
         unitary_diag=.false.)

    if(this%triangular.eq.'L') matrix_out%triangular='L'

    nrow = this%nrow
    ncol = this%ncol
    nterm = this%nterm

    matrix_out%ia(1:nrow+1) = this%ia
    matrix_out%ia(nrow+2)   = nterm+ncol+2
    
    matrix_out%ja(1:nterm) = this%ja
    do icol = 1,ncol+1
       matrix_out%ja(icol+nterm) = icol
    end do

    matrix_out%coeff(1:nterm) = this%coeff
    do icol = 1,ncol
       matrix_out%coeff(icol+nterm) = vec_w(icol)
    end do
    matrix_out%coeff(nterm+ncol+1) = sum(vec_w)    

    select type (this)
    type is (spmat)
       this = matrix_out
    end select
    
    call matrix_out%kill(lun_err)   
    
  end subroutine add_row

    function bandwidth ( this , perm, iperm) 
    !************************************************************************
    !
    !! ADJ_BANDWIDTH computes the bandwidth of an adjacency matrix.
    !  Enrico Facca 2018-09-09
    !  Al inputs included the spamt format
    !
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    11 March 2005
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Reference:
    !
    !    Alan George and Joseph Liu,
    !    Computer Solution of Large Sparse Positive Definite Systems,
    !    Prentice Hall, 1981.
    !
    !  Parameters:
    !
    !    Input, integer ( kind = 4 ) NODE_NUM, the number of nodes.
    !
    !    Input, integer ( kind = 4 ) ADJ_NUM, the number of adjacency entries.
    !
    !    Input, integer ( kind = 4 ) ADJ_ROW(NODE_NUM+1).  Information about row I
    !    is stored in entries ADJ_ROW(I) through ADJ_ROW(I+1)-1 of ADJ.
    !
    !    Input, integer ( kind = 4 ) ADJ(ADJ_NUM), the adjacency structure.
    !    For each row, it contains the column indices of the nonzero entries.
    !
    !    Output, integer ( kind = 4 ) ADJ_BANDWIDTH, the bandwidth of the adjacency
    !    matrix.
    !
    use Globals
    implicit none
    class(spmat),      intent(in   ) :: this
    integer, optional, intent(in   ) :: perm(this%nrow),iperm(this%nrow)
    integer ( kind = 4 ) :: bandwidth

    ! local 
    integer ( kind = 4 ) band_hi
    integer ( kind = 4 ) band_lo
    integer ( kind = 4 ) col
    integer ( kind = 4 ) i
    integer ( kind = 4 ) j

    band_lo = 0
    band_hi = 0


    if ( present(perm) .and. present(iperm) ) then
       do i = 1, this%nrow
          do j = this%ia(perm(i)), this%ia(perm(i)+1)-1
             col = iperm(this%ja(j))
             band_lo = max ( band_lo, i - col )
             band_hi = max ( band_hi, col - i )
          end do

       end do
    else
       do i = 1, this%nrow
          do j = this%ia(i), this%ia(i+1)-1
             col = this%ja(j)
             band_lo = max ( band_lo, i - col )
             band_hi = max ( band_hi, col - i )
          end do
       end do
    end if
          

    bandwidth = band_lo + 1 + band_hi

    return
  end function bandwidth

  
  subroutine genrcm ( this, lun_err,perm, inv_perm)
    !*************************************************************************
    !
    !! GENRCM finds the reverse Cuthill-Mckee ordering for a general graph.
    !  Enrico Facca :
    !   2018-09-09 The inputs of the algorithm are now included in the 
    !   csr format. Each irow-line of the matrix contians the index of the 
    !   graph nodes that are connnected to irow. All the indeces are listed in 
    !   thsi%ja
    !  
    !
    !  Discussion:
    !
    !    For each connected component in the graph, the routine obtains
    !    an ordering by calling RCM.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    04 January 2003
    !
    !  Author:
    !
    !    Original FORTRAN77 version by Alan George, Joseph Liu.
    !    FORTRAN90 version by John Burkardt.
    !
    !  Reference:
    !
    !    Alan George and Joseph Liu,
    !    Computer Solution of Large Sparse Positive Definite Systems,
    !    Prentice Hall, 1981.
    !
    !  Parameters:
    !
    !    Input, integer ( kind = 4 ) NODE_NUM, the number of nodes.
    !
    !    Input, integer ( kind = 4 ) ADJ_NUM, the number of adjacency entries.
    !
    !    Input, integer ( kind = 4 ) ADJ_ROW(NODE_NUM+1).  Information about row
    !    I is stored in entries ADJ_ROW(I) through ADJ_ROW(I+1)-1 of ADJ.
    !
    !    Input, integer ( kind = 4 ) ADJ(ADJ_NUM), the adjacency structure.
    !    For each row, it contains the column indices of the nonzero entries.
    !
    !    Output, integer ( kind = 4 ) PERM(NODE_NUM), the RCM ordering.
    !
    !  Local Parameters:
    !
    !    Local, integer LEVEL_ROW(NODE_NUM+1), the index vector for a level
    !    structure.  The level structure is stored in the currently unused
    !    spaces in the permutation vector PERM.
    !
    !    Local, integer MASK(NODE_NUM), marks variables that have been numbered.
    !
    use GLobals
    implicit none

    class(spmat), intent(in   ) :: this
    integer,      intent(in   ) :: lun_err    
    integer,      intent(inout) :: perm(this%nrow)
    integer,      intent(inout) :: inv_perm(this%nrow)
        
    ! local
    logical :: rc
    integer :: res
    integer ( kind = 4 ) adj_num
    integer ( kind = 4 ) node_num
    integer ( kind = 4 ) i
    integer ( kind = 4 ) iccsze
    integer ( kind = 4 ), allocatable :: mask(:)
    integer ( kind = 4 ) level_num
    integer ( kind = 4 ), allocatable ::level_row(:)
    integer ( kind = 4 ) num
    integer ( kind = 4 ) root

    node_num = this%nrow
    adj_num  = this%nterm

    ! allocate work arrays
    allocate(&
         mask(node_num),&
         level_row(node_num+1),&
         stat=res)
    if (res .ne. 0) &
         rc = IOerr(lun_err, err_alloc, 'genrcm', &
         'work arrays mask, level_row',res) 
    
    !build permutation 
    
    mask(1:node_num) = 1

    num = 1
    
    do i = 1, node_num
       !
       !  For each masked connected component...
       !
       if ( mask(i) /= 0 ) then

          root = i
          !
          !  Find a pseudo-peripheral node ROOT.  The level structure found by
          !  ROOT_FIND is stored starting at PERM(NUM).
          !
          call root_find ( root, adj_num, this%ia, this%ja,&
               mask, level_num, &
               level_row, perm(num), node_num )
          !
          !  RCM orders the component using ROOT as the starting node.
          !
          call rcm ( root, adj_num, this%ia, this%ja,&

               mask, perm(num), iccsze, &
               node_num )

          num = num + iccsze
          !
          !  We can stop once every node is in one of the connected components.
          !
          if ( node_num < num ) then
             exit
          end if

       end if

    end do

    do i = 1, node_num
       inv_perm(perm(i)) = i
    end do

    
    ! clean memory
    deallocate(&
         mask,&
         level_row,&
         stat=res)

    if (res .ne. 0) &
         rc = IOerr(lun_err, err_dealloc, 'genrcm', &
         'work arrays mask, level_row',res) 

  contains
    subroutine root_find ( root, adj_num, adj_row, adj, mask, level_num, &
         level_row, level, node_num )

      !************************************************************************0
      !
      !! ROOT_FIND finds a pseudo-peripheral node.
      !
      !  Discussion:
      !
      !    The diameter of a graph is the maximum distance (number of edges)
      !    between any two nodes of the graph.
      !
      !    The eccentricity of a node is the maximum distance between that
      !    node and any other node of the graph.
      !
      !    A peripheral node is a node whose eccentricity equals the
      !    diameter of the graph.
      !
      !    A pseudo-peripheral node is an approximation to a peripheral node;
      !    it may be a peripheral node, but all we know is that we tried our
      !    best.
      !
      !    The routine is given a graph, and seeks pseudo-peripheral nodes,
      !    using a modified version of the scheme of Gibbs, Poole and
      !    Stockmeyer.  It determines such a node for the section subgraph
      !    specified by MASK and ROOT.
      !
      !    The routine also determines the level structure associated with
      !    the given pseudo-peripheral node; that is, how far each node
      !    is from the pseudo-peripheral node.  The level structure is
      !    returned as a list of nodes LS, and pointers to the beginning
      !    of the list of nodes that are at a distance of 0, 1, 2, ...,
      !    NODE_NUM-1 from the pseudo-peripheral node.
      !
      !  Licensing:
      !
      !    This code is distributed under the GNU LGPL license.
      !
      !  Modified:
      !
      !    28 October 2003
      !
      !  Author:
      !
      !    Original FORTRAN77 version by Alan George, Joseph Liu.
      !    FORTRAN90 version by John Burkardt.
      !
      !  Reference:
      !
      !    Alan George and Joseph Liu,
      !    Computer Solution of Large Sparse Positive Definite Systems,
      !    Prentice Hall, 1981.
      !
      !    Norman Gibbs, William Poole, Paul Stockmeyer,
      !    An Algorithm for Reducing the Bandwidth and Profile of a Sparse Matrix,
      !    SIAM Journal on Numerical Analysis,
      !    Volume 13, pages 236-250, 1976.
      !
      !    Norman Gibbs,
      !    Algorithm 509: A Hybrid Profile Reduction Algorithm,
      !    ACM Transactions on Mathematical Software,
      !    Volume 2, pages 378-387, 1976.
      !
      !  Parameters:
      !
      !    Input/output, integer ( kind = 4 ) ROOT.  On input, ROOT is a node in the
      !    the component of the graph for which a pseudo-peripheral node is
      !    sought.  On output, ROOT is the pseudo-peripheral node obtained.
      !
      !    Input, integer ( kind = 4 ) ADJ_NUM, the number of adjacency entries.
      !
      !    Input, integer ( kind = 4 ) ADJ_ROW(NODE_NUM+1).  Information about row I
      !    is stored in entries ADJ_ROW(I) through ADJ_ROW(I+1)-1 of ADJ.
      !
      !    Input, integer ( kind = 4 ) ADJ(ADJ_NUM), the adjacency structure.
      !    For each row, it contains the column indices of the nonzero entries.
      !
      !    Input, integer ( kind = 4 ) MASK(NODE_NUM), specifies a section subgraph.  
      !    Nodes for which MASK is zero are ignored by FNROOT.
      !
      !    Output, integer ( kind = 4 ) LEVEL_NUM, is the number of levels in the 
      !    level structure rooted at the node ROOT.
      !
      !    Output, integer ( kind = 4 ) LEVEL_ROW(NODE_NUM+1), 
      !    integer LEVEL(NODE_NUM), the level structure array pair containing the 
      !    level structure found.
      !
      !    Input, integer ( kind = 4 ) NODE_NUM, the number of nodes.
      !
      implicit none

      integer ( kind = 4 ) adj_num
      integer ( kind = 4 ) node_num

      integer ( kind = 4 ) adj(adj_num)
      integer ( kind = 4 ) adj_row(node_num+1)
      integer ( kind = 4 ) iccsze
      integer ( kind = 4 ) j
      integer ( kind = 4 ) jstrt
      integer ( kind = 4 ) k
      integer ( kind = 4 ) kstop
      integer ( kind = 4 ) kstrt
      integer ( kind = 4 ) level(node_num)
      integer ( kind = 4 ) level_num
      integer ( kind = 4 ) level_num2
      integer ( kind = 4 ) level_row(node_num+1)
      integer ( kind = 4 ) mask(node_num)
      integer ( kind = 4 ) mindeg
      integer ( kind = 4 ) nabor
      integer ( kind = 4 ) ndeg
      integer ( kind = 4 ) node
      integer ( kind = 4 ) root
      !
      !  Determine the level structure rooted at ROOT.
      !
      call level_set ( root, adj_num, adj_row, adj, mask, level_num, &
           level_row, level, node_num )
      !
      !  Count the number of nodes in this level structure.
      !
      iccsze = level_row(level_num+1) - 1
      !
      !  Extreme case:
      !    A complete graph has a level set of only a single level.
      !    Every node is equally good (or bad).
      !
      if ( level_num == 1 ) then
         return
      end if
      !
      !  Extreme case:
      !    A "line graph" 0--0--0--0--0 has every node in its only level.
      !    By chance, we've stumbled on the ideal root.
      !
      if ( level_num == iccsze ) then
         return
      end if
      !
      !  Pick any node from the last level that has minimum degree
      !  as the starting point to generate a new level set.
      !
      do

         mindeg = iccsze

         jstrt = level_row(level_num)
         root = level(jstrt)

         if ( jstrt < iccsze ) then

            do j = jstrt, iccsze

               node = level(j)
               ndeg = 0
               kstrt = adj_row(node)
               kstop = adj_row(node+1)-1

               do k = kstrt, kstop
                  nabor = adj(k)
                  if ( 0 < mask(nabor) ) then
                     ndeg = ndeg+1
                  end if
               end do

               if ( ndeg < mindeg ) then
                  root = node
                  mindeg = ndeg
               end if

            end do

         end if
         !
         !  Generate the rooted level structure associated with this node.
         !
         call level_set ( root, adj_num, adj_row, adj, mask, level_num2, &
              level_row, level, node_num )
         !
         !  If the number of levels did not increase, accept the new ROOT.
         !
         if ( level_num2 <= level_num ) then
            exit
         end if

         level_num = level_num2
         !
         !  In the unlikely case that ROOT is one endpoint of a line graph,
         !  we can exit now.
         !
         if ( iccsze <= level_num ) then
            exit
         end if

      end do

      return
    end subroutine root_find

    subroutine level_set ( root, adj_num, adj_row, adj, mask, level_num, &
         level_row, level, node_num )

      !*****************************************************************************80
      !
      !! LEVEL_SET generates the connected level structure rooted at a given node.
      !
      !  Discussion:
      !
      !    Only nodes for which MASK is nonzero will be considered.
      !
      !    The root node chosen by the user is assigned level 1, and masked.
      !    All (unmasked) nodes reachable from a node in level 1 are
      !    assigned level 2 and masked.  The process continues until there
      !    are no unmasked nodes adjacent to any node in the current level.
      !    The number of levels may vary between 2 and NODE_NUM.
      !
      !  Licensing:
      !
      !    This code is distributed under the GNU LGPL license.
      !
      !  Modified:
      !
      !    28 October 2003
      !
      !  Author:
      !
      !    Original FORTRAN77 version by Alan George, Joseph Liu.
      !    FORTRAN90 version by John Burkardt.
      !
      !  Reference:
      !
      !    Alan George and Joseph Liu,
      !    Computer Solution of Large Sparse Positive Definite Systems,
      !    Prentice Hall, 1981.
      !
      !  Parameters:
      !
      !    Input, integer ( kind = 4 ) ROOT, the node at which the level structure
      !    is to be rooted.
      !
      !    Input, integer ( kind = 4 ) ADJ_NUM, the number of adjacency entries.
      !
      !    Input, integer ( kind = 4 ) ADJ_ROW(NODE_NUM+1).  Information about row I
      !    is stored in entries ADJ_ROW(I) through ADJ_ROW(I+1)-1 of ADJ.
      !
      !    Input, integer ( kind = 4 ) ADJ(ADJ_NUM), the adjacency structure.
      !    For each row, it contains the column indices of the nonzero entries.
      !
      !    Input/output, integer ( kind = 4 ) MASK(NODE_NUM).  On input, only nodes 
      !    with nonzero MASK are to be processed.  On output, those nodes which were
      !    included in the level set have MASK set to 1.
      !
      !    Output, integer ( kind = 4 ) LEVEL_NUM, the number of levels in the level
      !    structure.  ROOT is in level 1.  The neighbors of ROOT
      !    are in level 2, and so on.
      !
      !    Output, integer ( kind = 4 ) LEVEL_ROW(NODE_NUM+1), LEVEL(NODE_NUM), the 
      !    rooted level structure.
      !
      !    Input, integer ( kind = 4 ) NODE_NUM, the number of nodes.
      !
      implicit none

      integer ( kind = 4 ) adj_num
      integer ( kind = 4 ) node_num

      integer ( kind = 4 ) adj(adj_num)
      integer ( kind = 4 ) adj_row(node_num+1)
      integer ( kind = 4 ) i
      integer ( kind = 4 ) iccsze
      integer ( kind = 4 ) j
      integer ( kind = 4 ) jstop
      integer ( kind = 4 ) jstrt
      integer ( kind = 4 ) lbegin
      integer ( kind = 4 ) level_num
      integer ( kind = 4 ) level_row(node_num+1)
      integer ( kind = 4 ) level(node_num)
      integer ( kind = 4 ) lvlend
      integer ( kind = 4 ) lvsize
      integer ( kind = 4 ) mask(node_num)
      integer ( kind = 4 ) nbr
      integer ( kind = 4 ) node
      integer ( kind = 4 ) root

      mask(root) = 0
      level(1) = root
      level_num = 0
      lvlend = 0
      iccsze = 1
      !
      !  LBEGIN is the pointer to the beginning of the current level, and
      !  LVLEND points to the end of this level.
      !
      do

         lbegin = lvlend + 1
         lvlend = iccsze
         level_num = level_num + 1
         level_row(level_num) = lbegin
         !
         !  Generate the next level by finding all the masked neighbors of nodes
         !  in the current level.
         !
         do i = lbegin, lvlend

            node = level(i)
            jstrt = adj_row(node)
            jstop = adj_row(node+1)-1

            do j = jstrt, jstop

               nbr = adj(j)

               if ( mask(nbr) /= 0 ) then
                  iccsze = iccsze + 1
                  level(iccsze) = nbr
                  mask(nbr) = 0
               end if

            end do

         end do
         !
         !  Compute the current level width (the number of nodes encountered.)
         !  If it is positive, generate the next level.
         !
         lvsize = iccsze - lvlend

         if ( lvsize <= 0 ) then
            exit
         end if

      end do

      level_row(level_num+1) = lvlend + 1
      !
      !  Reset MASK to 1 for the nodes in the level structure.
      !
      mask(level(1:iccsze)) = 1

      return
    end subroutine level_set

    subroutine rcm ( root, adj_num, adj_row, adj, mask, perm, iccsze, node_num )

      !*****************************************************************************80
      !
      !! RCM renumbers a connected component by the reverse Cuthill McKee algorithm.
      !
      !  Discussion:
      !
      !    The connected component is specified by a node ROOT and a mask.
      !    The numbering starts at the root node.
      !
      !    An outline of the algorithm is as follows:
      !
      !    X(1) = ROOT.
      !
      !    for ( I = 1 to N-1)
      !      Find all unlabeled neighbors of X(I),
      !      assign them the next available labels, in order of increasing degree.
      !
      !    When done, reverse the ordering.
      !
      !  Licensing:
      !
      !    This code is distributed under the GNU LGPL license.
      !
      !  Modified:
      !
      !    05 December 2008
      !
      !  Author:
      !
      !    Original FORTRAN77 version by Alan George, Joseph Liu.
      !    FORTRAN90 version by John Burkardt.
      !
      !  Reference:
      !
      !    Alan George and Joseph Liu,
      !    Computer Solution of Large Sparse Positive Definite Systems,
      !    Prentice Hall, 1981.
      !
      !  Parameters:
      !
      !    Input, integer ( kind = 4 ) ROOT, the node that defines the connected 
      !    component.  It is used as the starting point for the RCM ordering.
      !
      !    Input, integer ( kind = 4 ) ADJ_NUM, the number of adjacency entries.
      !
      !    Input, integer ( kind = 4 ) ADJ_ROW(NODE_NUM+1).  Information about row I
      !    is stored in entries ADJ_ROW(I) through ADJ_ROW(I+1)-1 of ADJ.
      !
      !    Input, integer ( kind = 4 ) ADJ(ADJ_NUM), the adjacency structure.
      !    For each row, it contains the column indices of the nonzero entries.
      !
      !    Input/output, integer ( kind = 4 ) MASK(NODE_NUM), a mask for the nodes. 
      !    Only those nodes with nonzero input mask values are considered by the
      !    routine.  The nodes numbered by RCM will have their mask values
      !    set to zero.
      !
      !    Output, integer ( kind = 4 ) PERM(NODE_NUM), the RCM ordering.
      !
      !    Output, integer ( kind = 4 ) ICCSZE, the size of the connected component
      !    that has been numbered.
      !
      !    Input, integer ( kind = 4 ) NODE_NUM, the number of nodes.
      !
      !  Local parameters:
      !
      !    Workspace, integer DEG(NODE_NUM), a temporary vector used to hold
      !    the degree of the nodes in the section graph specified by mask and root.
      !
      implicit none

      integer ( kind = 4 ) adj_num
      integer ( kind = 4 ) node_num

      integer ( kind = 4 ) adj(adj_num)
      integer ( kind = 4 ) adj_row(node_num+1)
      integer ( kind = 4 ) deg(node_num)
      integer ( kind = 4 ) fnbr
      integer ( kind = 4 ) i
      integer ( kind = 4 ) iccsze
      integer ( kind = 4 ) j
      integer ( kind = 4 ) jstop
      integer ( kind = 4 ) jstrt
      integer ( kind = 4 ) k
      integer ( kind = 4 ) l
      integer ( kind = 4 ) lbegin
      integer ( kind = 4 ) lnbr
      integer ( kind = 4 ) lperm
      integer ( kind = 4 ) lvlend
      integer ( kind = 4 ) mask(node_num)
      integer ( kind = 4 ) nbr
      integer ( kind = 4 ) node
      integer ( kind = 4 ) perm(node_num)
      integer ( kind = 4 ) root
      !
      !  Find the degrees of the nodes in the component specified by MASK and ROOT.
      !
      call degree ( root, adj_num, adj_row, adj, mask, deg, iccsze, perm, node_num )

      mask(root) = 0

      if ( iccsze <= 1 ) then
         return
      end if

      lvlend = 0
      lnbr = 1
      !
      !  LBEGIN and LVLEND point to the beginning and
      !  the end of the current level respectively.
      !
      do while ( lvlend < lnbr )

         lbegin = lvlend + 1
         lvlend = lnbr

         do i = lbegin, lvlend
            !
            !  For each node in the current level...
            !
            node = perm(i)
            jstrt = adj_row(node)
            jstop = adj_row(node+1) - 1
            !
            !  Find the unnumbered neighbors of NODE.
            !
            !  FNBR and LNBR point to the first and last neighbors
            !  of the current node in PERM.
            !
            fnbr = lnbr + 1

            do j = jstrt, jstop

               nbr = adj(j)

               if ( mask(nbr) /= 0 ) then
                  lnbr = lnbr + 1
                  mask(nbr) = 0
                  perm(lnbr) = nbr
               end if

            end do
            !
            !  If no neighbors, skip to next node in this level.
            !
            if ( lnbr <= fnbr ) then
               cycle
            end if
            !
            !  Sort the neighbors of NODE in increasing order by degree.
            !  Linear insertion is used.
            !
            k = fnbr

            do while ( k < lnbr )

               l = k
               k = k + 1
               nbr = perm(k)

               do while ( fnbr < l )

                  lperm = perm(l)

                  if ( deg(lperm) <= deg(nbr) ) then
                     exit
                  end if

                  perm(l+1) = lperm
                  l = l-1

               end do

               perm(l+1) = nbr

            end do

         end do

      end do
      !
      !  We now have the Cuthill-McKee ordering.  Reverse it.
      !
      call i4vec_reverse ( iccsze, perm )

      return
    end subroutine rcm

    subroutine i4vec_reverse ( n, a )

      !*****************************************************************************80
      !
      !! I4VEC_REVERSE reverses the elements of an integer vector.
      !
      !  Example:
      !
      !    Input:
      !
      !      N = 5,
      !      A = ( 11, 12, 13, 14, 15 ).
      !
      !    Output:
      !
      !      A = ( 15, 14, 13, 12, 11 ).
      !
      !  Licensing:
      !
      !    This code is distributed under the GNU LGPL license.
      !
      !  Modified:
      !
      !    26 July 1999
      !
      !  Author:
      !
      !    John Burkardt
      !
      !  Parameters:
      !
      !    Input, integer ( kind = 4 ) N, the number of entries in the array.
      !
      !    Input/output, integer ( kind = 4 ) A(N), the array to be reversed.
      !
      implicit none

      integer ( kind = 4 ) n

      integer ( kind = 4 ) a(n)
      integer ( kind = 4 ) i

      do i = 1, n/2
         call i4_swap ( a(i), a(n+1-i) )
      end do

      return
    end subroutine i4vec_reverse
    subroutine i4_swap ( i, j )

      !*****************************************************************************80
      !
      !! I4_SWAP swaps two integer values.
      !
      !  Licensing:
      !
      !    This code is distributed under the GNU LGPL license.
      !
      !  Modified:
      !
      !    30 November 1998
      !
      !  Author:
      !
      !    John Burkardt
      !
      !  Parameters:
      !
      !    Input/output, integer ( kind = 4 ) I, J.  On output, the values of I and
      !    J have been interchanged.
      !
      implicit none

      integer ( kind = 4 ) i
      integer ( kind = 4 ) j
      integer ( kind = 4 ) k

      k = i
      i = j
      j = k

      return
    end subroutine i4_swap

    subroutine degree ( root, adj_num, adj_row, adj, mask, deg, iccsze, ls, &
         node_num )

      !*****************************************************************************80
      !
      !! DEGREE computes the degrees of the nodes in the connected component.
      !
      !  Discussion:
      !
      !    The connected component is specified by MASK and ROOT.
      !    Nodes for which MASK is zero are ignored.
      !
      !  Licensing:
      !
      !    This code is distributed under the GNU LGPL license.
      !
      !  Modified:
      !
      !   05 January 2003
      !
      !  Author:
      !
      !    Original FORTRAN77 version by Alan George, Joseph Liu.
      !    FORTRAN90 version by John Burkardt.
      !
      !  Reference:
      !
      !    Alan George and Joseph Liu,
      !    Computer Solution of Large Sparse Positive Definite Systems,
      !    Prentice Hall, 1981.
      !
      !  Parameters:
      !
      !    Input, integer ( kind = 4 ) ROOT, the node that defines the 
      !    connected component.
      !
      !    Input, integer ( kind = 4 ) ADJ_NUM, the number of adjacency entries.
      !
      !    Input, integer ( kind = 4 ) ADJ_ROW(NODE_NUM+1).  Information about row I 
      !    is stored in entries ADJ_ROW(I) through ADJ_ROW(I+1)-1 of ADJ.
      !
      !    Input, integer ( kind = 4 ) ADJ(ADJ_NUM), the adjacency structure.
      !    For each row, it contains the column indices of the nonzero entries.
      !
      !    Input, integer ( kind = 4 ) MASK(NODE_NUM), is nonzero for those nodes 
      !    which are to be considered.
      !
      !    Output, integer ( kind = 4 ) DEG(NODE_NUM), contains, for each  node in 
      !    the connected component, its degree.
      !
      !    Output, integer ( kind = 4 ) ICCSIZE, the number of nodes in the connected
      !    component.
      !
      !    Output, integer ( kind = 4 ) LS(NODE_NUM), stores in entries 1 through 
      !    ICCSIZE the nodes in the connected component, starting with ROOT, and 
      !    proceeding by levels.
      !
      !    Input, integer ( kind = 4 ) NODE_NUM, the number of nodes.
      !
      implicit none

      integer ( kind = 4 ) adj_num
      integer ( kind = 4 ) node_num

      integer ( kind = 4 ) adj(adj_num)
      integer ( kind = 4 ) adj_row(node_num+1)
      integer ( kind = 4 ) deg(node_num)
      integer ( kind = 4 ) i
      integer ( kind = 4 ) iccsze
      integer ( kind = 4 ) ideg
      integer ( kind = 4 ) j
      integer ( kind = 4 ) jstop
      integer ( kind = 4 ) jstrt
      integer ( kind = 4 ) lbegin
      integer ( kind = 4 ) ls(node_num)
      integer ( kind = 4 ) lvlend
      integer ( kind = 4 ) lvsize
      integer ( kind = 4 ) mask(node_num)
      integer ( kind = 4 ) nbr
      integer ( kind = 4 ) node
      integer ( kind = 4 ) root
      !
      !  The sign of ADJ_ROW(I) is used to indicate if node I has been considered.
      !
      ls(1) = root
      adj_row(root) = -adj_row(root)
      lvlend = 0
      iccsze = 1
      !
      !  LBEGIN is the pointer to the beginning of the current level, and
      !  LVLEND points to the end of this level.
      !
      do

         lbegin = lvlend + 1
         lvlend = iccsze
         !
         !  Find the degrees of nodes in the current level,
         !  and at the same time, generate the next level.
         !
         do i = lbegin, lvlend

            node = ls(i)
            jstrt = -adj_row(node)
            jstop = abs ( adj_row(node+1) ) - 1
            ideg = 0

            do j = jstrt, jstop

               nbr = adj(j)

               if ( mask(nbr) /= 0 ) then

                  ideg = ideg + 1

                  if ( 0 <= adj_row(nbr) ) then
                     adj_row(nbr) = -adj_row(nbr)
                     iccsze = iccsze + 1
                     ls(iccsze) = nbr
                  end if

               end if

            end do

            deg(node) = ideg

         end do
         !
         !  Compute the current level width.
         !
         lvsize = iccsze - lvlend
         !
         !  If the current level width is nonzero, generate another level.
         !
         if ( lvsize == 0 ) then
            exit
         end if

      end do
      !
      !  Reset ADJ_ROW to its correct sign and return.
      !
      do i = 1, iccsze
         node = ls(i)
         adj_row(node) = -adj_row(node)
      end do

      return
    end subroutine degree

  end subroutine genrcm
  
  
  subroutine sort_row(nnz,ja,coeff)
      use Globals
      implicit none
      integer,           intent(in   ) :: nnz
      integer,           intent(inout) :: ja(nnz)
      real(kind=double), intent(inout) :: coeff(nnz)
      !local 
      integer :: i,j, indx,isgn, itemp
      real(kind=double) :: rtemp

      if (nnz.lt.2) return


      !  Initialize.
      i = 0
      indx = 0
      isgn = 0
      j = 0
      do 
         call global_heapsort(nnz, indx, i,j,isgn)
         if (indx .gt. 0 ) then
            ! SWAP ELEMENT 

            ! swap column indeces
            itemp = ja(i)
            ja(i) = ja(j)
            ja(j) = itemp

            ! swap real nnzero coeff
            rtemp    = coeff(i)
            coeff(i) = coeff(j)
            coeff(j) = rtemp
         else if ( indx .lt. 0) then
            ! COMPARE
            isgn = 0
            if ( ja(i) .lt.  ja(j) ) isgn = -1
            if ( ja(i) .gt.  ja(j) ) isgn = 1
         else if ( indx .eq. 0 ) then
            exit
         end if
      end do
    end subroutine sort_row


    subroutine plot2ps ( this, job, iounit )

      !*****************************************************************************80
      !
      !! PLTMTPS creates a PostScript plot of a sparse matrix.
      !
      !  Discussion:
      !
      !    This routine creates a 'PS' file for plotting the pattern of
      !    a sparse matrix stored in general sparse format. It can be used
      !    for inserting matrix plots in a text. The size of the plot can be
      !    7in x 7in or 5 in x 5in ..
      !
      !    1) Plots square as well as rectangular matrices.
      !    2) Does not writer a caption yet.
      !    3) No bounding box put in yet
      !
      !  Modified:
      !
      !    07 January 2004
      !
      !  Author:
      !
      !    Paul Frederickson
      !
      !  Parameters:
      !
      !    Input, integer ( kind = 4 ) NROW, the row dimension of the matrix.
      !
      !    Input, integer ( kind = 4 ) NCOL, the column dimension of the matrix.
      !
      !    Input, integer ( kind = 4 ) MODE, indicates the matrix storage mode:
      !    0, by rows;
      !    1, by columns.
      !
      ! ja     = column indices of nonzero elements when matrix is
      !         stored rowise. Row indices if stores column-wise.
      ! ia     = integer ( kind = 4 ) array of containing the pointers to the
      !         beginning of the columns in arrays a, ja.
      !
      ! title  = character*72 = title of matrix test ( character a*72 ).
      ! key    = character*8  = key of matrix
      ! type   = character*3  = type of matrix.
      !
      ! job, integer ( kind = 4 ). tells pltmt whether or not to reduce
      ! the plot.
      !           if enabled then the standard size of 7in will be
      !           replaced by a 5in plot.
      !          job = 0 : do not reduce
      !          job = 1 : reduce plot to 5 inches.
      !
      ! iounit = logical unit number where to write the matrix into.
      !
      use Globals
      implicit none
      
      class(spmat),     intent(in ) :: this
      integer,          intent(in ) :: job
      ! local
      real ( kind = 8 ) delta
      integer ( kind = 4 ) ii
      integer ( kind = 4 ) ilast
      integer ( kind = 4 ) iounit
      integer ( kind = 4 ) istart


      integer ( kind = 4 ) k
      character ( len = 8 ) key
      integer ( kind = 4 ) m
      integer ( kind = 4 ) maxdim
      integer ( kind = 4 ) mode
      integer ( kind = 4 ) n
      integer ( kind = 4 ) ncol
      integer ( kind = 4 ) nrow
      integer ( kind = 4 ) nnz
      character ( len = 3 ) type

      nrow = this%nrow
      ncol = this%ncol
      
      mode = 1
      type = this%storage_system

      if ( mode == 0 ) then
         n = nrow
      else
         n = ncol
      end if

      nnz = this%ia(n+1) - this%ia(1)
      maxdim = max ( nrow, ncol )
      m = 1 + maxdim
      !
      !  Keep this test as in old pltmt (for future changes).
      !
      if ( mod ( job, 10 ) == 1 ) then
         delta = 72.0D+00 * 5.0D+00 / ( 2.0D+00 + maxdim )
      else
         delta = 72.0D+00 * 7.0D+00 / (2.0D+00 + maxdim )
      end if

      write(iounit,*)'%!PS'
      write(iounit,*)' gsave 50 50 translate'
      write(iounit,*) delta, delta, ' scale'
      write(iounit,*) ' 0.25 setlinewidth'

      if ( mod ( job, 10 ) == 1 ) then
         write (iounit,*) ' 23 55 translate'
      else
         write (iounit,*) ' 2 35 translate'
      end if

      write(iounit,*) ' newpath'
      write(iounit,*) 0,0,' moveto'
      write(iounit,*) m,0,' lineto'
      write(iounit,*) m,m,' lineto'
      write(iounit,*) 0,m,' lineto'
      write(iounit,*) ' closepath stroke'
      write(iounit,*) ' 1 1 translate'
      write(iounit,*) ' 0.5 setlinewidth'
      write(iounit,*) ' /p {moveto 0 -.25 rmoveto '
      write(iounit,*) '            0  .50 rlineto stroke} def'
      !
      !  Plotting loop
      !
      do ii = 1, n

         istart = this%ia(ii)
         ilast  = this%ia(ii+1)-1

         if ( mode /= 0 ) then

            do k = istart, ilast
               write(iounit,*) ii-1, nrow-this%ja(k), ' p'
            end do

         else
            !             y = xnrow - real ( ii, kind = 8 )
            do k = istart, ilast
               !               x = real ( this%ja(k) - 1, kind = 8 )
               write(iounit,*) this%ja(k)-1, nrow-ii, ' p'
            end do

         end if

      end do

      write(iounit,*)' showpage grestore'
130   format('Dimension: ',i4,' x ',i4',  Nonzero elements: ',i5)
      return
    end subroutine plot2ps




  end module SparseMatrix
