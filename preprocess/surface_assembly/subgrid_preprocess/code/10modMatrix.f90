module Matrix
  use Globals
  implicit none
  private
  !> Abstract type definition for matrix-type, which represents
  !> any linear operator for REAL^{NCOL} to REAL^{NROW}
  !> Additional label are included to described the matrix structure
  !> that will limit or active matrix operations. For example only
  !> symmetric matrices can be passed to PCG procedure, or can have 
  !> a Cholesky decomposition. While triangular matric can be easily solved.
  !> The Abstract inferface requires that any type extdending the 
  !> class abs_matrix must have a procedure defining the operation
  !> matrix-times-vector and matrix(transpose)-times-vector 
  ! TODO:
  ! Note
  ! The intent(inout) for abs_matrix in the muliptlication
  ! procedure is due to the fact that scratch arrays 
  ! included in the concrete may be used 
  ! (e.g. for the matrix-vector operation M times v with
  !  M = A^T A a scratch array is required )
  ! This can be removed passing an optional scratch type with the 
  ! proper size. 
  ! PRO:    intent(in) used, that avoid modification 
  ! CONTRO: not self-contained 
  public :: abs_matrix,multiplication,multiplication_transpose
  type, abstract :: abs_matrix
     !> Number of rows
     integer :: nrow=0
     !> Number of columns
     integer :: ncol=0
     !> Flag for squared matrix
     !> i.e. if nrow .eq. ncol
     logical :: is_squared=.false.
     !> Logical flag for sysmmetric matrix
     logical :: is_symmetric=.false.
     !> Character indicating if the matrix is 
     !> lower ('L')  or upper ('U') or not ('N')
     !> which is the default case
     character(len=1) :: triangular='N'
     !> Dimension of the Kernel space
     integer :: dim_kernel=0
     ! Dimension (ncol,dim_kernel)
     ! Null space of the matrix
     real(kind=double), allocatable :: kernel(:,:)
     !> Logical flag matrix with one on the diagonal
     ! TODO: Used only in SparseMatrix Module
     ! remove?
     logical :: unitary_diag=.false.
   contains
     !>-----------------------------------------------------------
     !> Procedure for computation of (matrix) times (vector)
     procedure(multiplication), deferred :: Mxv
     !>-----------------------------------------------------------
     !> Procedure for computation of (matrix transposed) times (vector)
     procedure(multiplication_transpose), deferred :: MTxv
  end type abs_matrix
  abstract interface
     !>-------------------------------------------------------------
     !> Abstract procedure defining the interface for a general
     !> matrix-vector multiplication
     !>         vec_out = (M) times (vec_in)
     !> (public procedure for class abs_matrix)
     !> 
     !> usage:
     !>     call 'var'%Mxv(vec_in,vec_out,[info])
     !>
     !> where 
     !> \param[in   ] vec_in          -> real, dimension('var'%ncol)
     !>                                   vector to be multiplied
     !> \param[inout] vec_out         -> real, dimension('var'%nrow)
     !>                                    vector (M) times (vec_in) 
     !> \param[in   ] (optional) info -> integer. Info number
     !>                                  in case of error   
     !<-------------------------------------------------------------
     subroutine multiplication(this,vec_in,vec_out,info)
       use Globals
       import abs_matrix
       implicit none
       class(abs_matrix), intent(inout) :: this
       real(kind=double), intent(in   ) :: vec_in(this%ncol)
       real(kind=double), intent(inout) :: vec_out(this%nrow)
       integer, optional, intent(inout) :: info
     end subroutine multiplication
      !>-------------------------------------------------------------
     !> Abstract procedure defining the interface for a general
     !> matrix(transpose)-vector multiplication
     !>         vec_out = (M)^T times (vec_in)
     !> (public procedure for class abs_matrix)
     !> 
     !> usage:
     !>     call 'var'%MTxv(vec_in,vec_out,[info])
     !>
     !> where 
     !> \param[in   ] vec_in          -> real, dimension('var'%nrow)
     !>                                   vector to be multiplied
     !> \param[inout] vec_out         -> real, dimension('var'%ncol)
     !>                                    vector (M) times (vec_in) 
     !> \param[in   ] (optional) info -> integer. Info number
     !>                                  in case of error   
     !<-------------------------------------------------------------
     subroutine multiplication_transpose(this,vec_in,vec_out,info)
       use Globals
       import abs_matrix
       implicit none
       class(abs_matrix), intent(inout) :: this
       real(kind=double), intent(in   ) :: vec_in(this%nrow)
       real(kind=double), intent(inout) :: vec_out(this%ncol)
       integer, optional, intent(inout) :: info
     end subroutine multiplication_transpose
  end interface
  !> Derived type used to create list of pointers to abstract matrices
  !> (Fortan does not allow to create arrays of pointers)
  !> E.g. in a Block-matrix
  !> M = (A D) with A and D of different type, we can create
  !> an array type that points toward A and D.
  !> Declarations: 
  !>  type(spmat)     :: A
  !>  type(diagmat  ) :: D
  !>  type(array_mat) :: list(2)
  !> Assignment
  !>  list(1)%mat => A
  !>  list(2)%mat => D
  !> See BlockMatrix module for a concrete example
  type, public :: array_mat
     !> Dimension (nmats)
     !> Array that contains the non-zero blocks
     class(abs_matrix), pointer :: mat
  end type array_mat  
end module Matrix




