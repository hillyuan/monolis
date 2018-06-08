module mod_monolis_converge
  use mod_monolis_prm
  use mod_monolis_com
  use mod_monolis_mat
  use mod_monolis_linalg
  use mod_monolis_linalg_util
  use mod_monolis_linalg_com
  implicit none

  real(kind=kdouble), save :: B2

contains

  subroutine monolis_set_converge(monoPRM, monoCOM, monoMAT, B, tcomm)
    implicit none
    type(monolis_prm) :: monoPRM
    type(monolis_com) :: monoCOM
    type(monolis_mat) :: monoMAT
    real(kind=kdouble) :: B(:)
    real(kind=kdouble) :: t1, t2
    real(kind=kdouble), optional :: tcomm

    call monolis_inner_product_R(monoCOM, monoMAT, monoMAT%NDOF, B, B, B2, tcomm)

  end subroutine monolis_set_converge

  subroutine monolis_check_converge(monoPRM, monoCOM, monoMAT, R, iter, is_converge, tcomm)
    implicit none
    type(monolis_prm) :: monoPRM
    type(monolis_com) :: monoCOM
    type(monolis_mat) :: monoMAT
    integer(kind=kint) :: iter
    real(kind=kdouble) :: R(:), R2, resid
    real(kind=kdouble) :: t1, t2
    real(kind=kdouble), optional :: tcomm
    logical :: is_converge

    is_converge = .false.

    call monolis_inner_product_R(monoCOM, monoMAT, monoMAT%NDOF, R, R, R2, tcomm)
    resid = dsqrt(R2/B2)

    if(monoCOM%myrank == 0) write (*,"(i7, 1pe16.6)") iter, resid
    if(resid < monoPRM%tol) is_converge = .true.

  end subroutine monolis_check_converge

end module mod_monolis_converge