module mod_monolis_solver_IR
  use mod_monolis_prm
  use mod_monolis_com
  use mod_monolis_mat
  use mod_monolis_precond
  use mod_monolis_matvec
  use mod_monolis_linalg
  use mod_monolis_linalg_util
  use mod_monolis_converge

  implicit none

contains

  subroutine monolis_solver_IR(monoPRM, monoCOM, monoMAT)
    use mod_monolis_prm
    use mod_monolis_com
    use mod_monolis_mat
    implicit none
    type(monolis_prm) :: monoPRM
    type(monolis_com) :: monoCOM
    type(monolis_mat) :: monoMAT
    integer(kind=kint) :: N, NP, NDOF, NNDOF
    integer(kind=kint) :: i, iter
    real(kind=kdouble) :: tol, resid, R2, B2
    real(kind=kdouble) :: t1, t2, tsol, tcomm
    real(kind=kdouble), pointer :: B(:), X(:)
    real(kind=kdouble), allocatable :: R(:), D(:)

    t1 = monolis_wtime()

    N     = monoMAT%N
    NP    = monoMAT%NP
    NDOF  = monoMAT%NDOF
    NNDOF = N*NDOF
    X => monoMAT%X; X = 0.0d0
    B => monoMAT%B
    tol = monoPRM%tol

    allocate(R(NDOF*NP)); R = 0.0d0
    allocate(D(NDOF*NP)); D = 0.0d0

    call monolis_IR_setup(monoPRM, monoCOM, monoMAT)
    call monolis_inner_product_R(monoCOM, monoMAT, NDOF, B, B, B2, tcomm)

    do i = 1,NNDOF
      R(i) = B(i)
    enddo

    do iter = 1, monoPRM%maxiter
      call monolis_IR_apply(monoPRM, monoCOM, monoMAT, R, D)

      do i=1,NNDOF
        X(i) = X(i) + D(i)
      enddo

      call monolis_residual(monoCOM, monoMAT, X, B, R, tcomm)
      call monolis_inner_product_R(monoCOM, monoMAT, NDOF, R, R, R2, tcomm)

      resid = dsqrt(R2/B2)

      if(monoCOM%myrank == 0) write (*,"(i7, 1pe16.6)") iter, resid
      if(resid <= tol) exit
    enddo

    call monolis_IR_clear(monoPRM, monoCOM, monoMAT)
    call monolis_update_R(monoCOM, NDOF, X, tcomm)

    deallocate(R)
    deallocate(D)

    t2 = monolis_wtime()
    tsol = t2 - t1
  end subroutine monolis_solver_IR

  subroutine monolis_IR_setup(monoPRM, monoCOM, monoMAT)
    implicit none
    type(monolis_prm) :: monoPRM
    type(monolis_com) :: monoCOM
    type(monolis_mat) :: monoMAT

    if(monoPRM%precond == monolis_prec_DIAG)then
      call monolis_precond_diag_setup(monoPRM, monoCOM, monoMAT)
    endif
  end subroutine monolis_IR_setup

  subroutine monolis_IR_apply(monoPRM, monoCOM, monoMAT, X, Y)
    implicit none
    type(monolis_prm) :: monoPRM
    type(monolis_com) :: monoCOM
    type(monolis_mat) :: monoMAT
    integer(kind=kint) :: i
    real(kind=kdouble) :: X(:), Y(:)

    if(monoPRM%precond == monolis_prec_DIAG)then
      call monolis_precond_diag_apply(monoPRM, monoCOM, monoMAT, X, Y)
    else
      do i = 1, monoMAT%N*monoMAT%NDOF
        Y(i) = X(i)
      enddo
    endif
  end subroutine monolis_IR_apply

  subroutine monolis_IR_clear(monoPRM, monoCOM, monoMAT)
    implicit none
    type(monolis_prm) :: monoPRM
    type(monolis_com) :: monoCOM
    type(monolis_mat) :: monoMAT

    if(monoPRM%precond == monolis_prec_DIAG)then
      call monolis_precond_diag_clear(monoPRM, monoCOM, monoMAT)
    endif
  end subroutine monolis_IR_clear

end module mod_monolis_solver_IR
