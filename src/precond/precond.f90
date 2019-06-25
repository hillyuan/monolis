module mod_monolis_precond
  use mod_monolis_prm
  use mod_monolis_com
  use mod_monolis_mat
  use mod_monolis_precond_diag
  use mod_monolis_precond_ilu
  use mod_monolis_precond_Jacobi
  use mod_monolis_precond_SOR
  use mod_monolis_util

  implicit none

contains

  subroutine monolis_precond_setup(monoPRM, monoCOM, monoMAT)
    implicit none
    type(monolis_prm) :: monoPRM
    type(monolis_com) :: monoCOM
    type(monolis_mat) :: monoMAT

    if(monoPRM%is_debug) call monolis_debug_header("monolis_precond_setup")

    if(monoPRM%precond == monolis_prec_DIAG)then
      call monolis_precond_diag_setup(monoPRM, monoCOM, monoMAT)
    elseif(monoPRM%precond == monolis_prec_ILU)then
      call monolis_precond_ilu_setup(monoPRM, monoCOM, monoMAT)
    elseif(monoPRM%precond == monolis_prec_JACOBI)then
      call monolis_precond_Jacobi_setup(monoPRM, monoCOM, monoMAT)
    elseif(monoPRM%precond == monolis_prec_SOR)then
      call monolis_precond_SOR_setup(monoPRM, monoCOM, monoMAT)
    endif
  end subroutine monolis_precond_setup

  subroutine monolis_precond_apply(monoPRM, monoCOM, monoMAT, X, Y)
    implicit none
    type(monolis_prm) :: monoPRM
    type(monolis_com) :: monoCOM
    type(monolis_mat) :: monoMAT
    integer(kind=kint) :: i
    real(kind=kdouble) :: X(:), Y(:)

    if(monoPRM%precond == monolis_prec_DIAG)then
      call monolis_precond_diag_apply(monoPRM, monoCOM, monoMAT, X, Y)
    elseif(monoPRM%precond == monolis_prec_ILU)then
      call monolis_precond_ilu_apply(monoPRM, monoCOM, monoMAT, X, Y)
    elseif(monoPRM%precond == monolis_prec_JACOBI)then
      call monolis_precond_Jacobi_apply(monoPRM, monoCOM, monoMAT, X, Y)
    elseif(monoPRM%precond == monolis_prec_SOR)then
      call monolis_precond_SOR_apply(monoPRM, monoCOM, monoMAT, X, Y)
    else
      do i = 1, monoMAT%N*monoMAT%NDOF
        Y(i) = X(i)
      enddo
    endif
  end subroutine monolis_precond_apply

  subroutine monolis_precond_clear(monoPRM, monoCOM, monoMAT)
    implicit none
    type(monolis_prm) :: monoPRM
    type(monolis_com) :: monoCOM
    type(monolis_mat) :: monoMAT

    if(monoPRM%is_debug) call monolis_debug_header("monolis_precond_clear")

    if(monoPRM%precond == monolis_prec_DIAG)then
      call monolis_precond_diag_clear(monoPRM, monoCOM, monoMAT)
    elseif(monoPRM%precond == monolis_prec_ILU)then
      call monolis_precond_ilu_clear(monoPRM, monoCOM, monoMAT)
    elseif(monoPRM%precond == monolis_prec_JACOBI)then
      call monolis_precond_Jacobi_clear(monoPRM, monoCOM, monoMAT)
    elseif(monoPRM%precond == monolis_prec_SOR)then
      call monolis_precond_SOR_clear(monoPRM, monoCOM, monoMAT)
    endif
  end subroutine monolis_precond_clear

end module mod_monolis_precond