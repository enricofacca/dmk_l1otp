module GeneralSparseMatrix
  use Globals
  use Matrix
  implicit none
  private 
  public :: gen_spmat, extract_diagonal,scale_by_diagonal 
  type, abstract, extends(abs_matrix) :: gen_spmat
     ! empty type
   contains
     !>-----------------------------------------------------------
     !> Procedure to extract the diagonal from a matrix
     procedure(extract_diagonal), deferred :: get_diagonal
     !>-----------------------------------------------------------
     !> Procedure to scale a matrix M with a diagonal matrix D obtaing
     !>                M(out) = D M D
     !> Procedure restriced to square matrix
     !>-----------------------------------------------------------
     procedure(scale_by_diagonal), deferred :: diagonal_scale
  end type gen_spmat
  abstract interface
     subroutine extract_diagonal(this,diagonal)
       use Globals
       import gen_spmat
       class(gen_spmat),  intent(in)  :: this
       real(kind=double), intent(out) :: diagonal(min(this%nrow,this%ncol))
     end subroutine extract_diagonal
     subroutine scale_by_diagonal(this,lun_err,diagonal)
       use Globals
       import gen_spmat
       class(gen_spmat),  intent(inout) :: this
       integer,           intent(in   ) :: lun_err
       real(kind=double), intent(in   ) :: diagonal(this%ncol)
     end subroutine scale_by_diagonal
  end interface
end module GeneralSparseMatrix
