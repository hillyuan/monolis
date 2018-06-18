module mod_monolis_solver_PipeBiCGSTAB_noprec
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

  subroutine monolis_solver_PipeBiCGSTAB_noprec(monoPRM, monoCOM, monoMAT)
    implicit none
    type(monolis_prm) :: monoPRM
    type(monolis_com) :: monoCOM
    type(monolis_mat) :: monoMAT
    integer(kind=kint) :: N, NP, NDOF, NNDOF
    integer(kind=kint) :: i, iter, iter_RR
    real(kind=kdouble) :: t1, t2, tsol, tcomm, CG(5), RR, RW, RR1, RS, RZ, R2, QY, YY
    real(kind=kdouble) :: alpha, beta, rho, rho1, omega, omega1
    real(kind=kdouble), allocatable :: W(:,:)
    real(kind=kdouble), pointer :: B(:), X(:)
    integer(kind=kint), parameter :: R  = 1
    integer(kind=kint), parameter :: R0 = 2
    integer(kind=kint), parameter :: W0 = 3
    integer(kind=kint), parameter :: T  = 4
    integer(kind=kint), parameter :: S  = 5
    integer(kind=kint), parameter :: P  = 6
    integer(kind=kint), parameter :: Z  = 7
    integer(kind=kint), parameter :: Q  = 8
    integer(kind=kint), parameter :: Y  = 9
    integer(kind=kint), parameter :: V  = 10
    logical :: is_converge

    t1 = monolis_wtime()

    N     = monoMAT%N
    NP    = monoMAT%NP
    NDOF  = monoMAT%NDOF
    NNDOF = N*NDOF
    X => monoMAT%X; X = 1.0d0
    B => monoMAT%B
    iter_RR = 50

    allocate(W(NDOF*NP,10))
    W = 0.0d0

    call monolis_set_converge(monoPRM, monoCOM, monoMAT, B, tcomm)
    call monolis_residual(monoCOM, monoMAT, X, B, W(:,R), tcomm)

    do i = 1, NNDOF
      W(i,R0) = W(i,R)
    enddo

    call monolis_matvec(monoCOM, monoMAT, W(:,R), W(:,W0), tcomm)
    call monolis_matvec(monoCOM, monoMAT, W(:,W0), W(:,T), tcomm)
    call monolis_inner_product_R(monoCOM, monoMAT, NDOF, W(:,R), W(:,R ), RR, tcomm)
    call monolis_inner_product_R(monoCOM, monoMAT, NDOF, W(:,R), W(:,W0), RW, tcomm)

    alpha = RR / RW
    beta  = 0.0d0
    omega = 0.0d0

    do iter = 1, monoPRM%maxiter
      do i = 1, NNDOF
        W(i,P) = W(i,R ) + beta *(W(i,P) - omega*W(i,S))
        W(i,S) = W(i,W0) + beta *(W(i,S) - omega*W(i,Z))
        W(i,Z) = W(i,T ) + beta *(W(i,Z) - omega*W(i,V))
        W(i,Q) = W(i,R ) - alpha* W(i,S)
        W(i,Y) = W(i,W0) - alpha* W(i,Z)
      enddo

      call monolis_inner_product_R(monoCOM, monoMAT, NDOF, W(:,Q), W(:,Y), CG(1), tcomm)
      call monolis_inner_product_R(monoCOM, monoMAT, NDOF, W(:,Y), W(:,Y), CG(2), tcomm)

      call monolis_matvec(monoCOM, monoMAT, W(:,Z), W(:,V), tcomm)

      QY = CG(1)
      YY = CG(2)
      omega1 = QY / YY

      if(mod(iter, iter_RR) == 0)then
        do i = 1, NNDOF
          X(i) = X(i) + alpha*W(i,P) + omega1*W(i,T)
        enddo
        call monolis_residual(monoCOM, monoMAT, X, B, W(:,R), tcomm)
        call monolis_matvec(monoCOM, monoMAT, W(:,R), W(:,W0), tcomm)
      else
        do i = 1, NNDOF
          X(i)    = X(i)    + alpha* W(i,P) + omega1*W(i,Q)
          W(i,R ) = W(i,Q) - omega1* W(i,Y)
          W(i,W0) = W(i,Y) - omega1*(W(i,T) - alpha *W(i,V))
        enddo
      endif

      call monolis_inner_product_R(monoCOM, monoMAT, NDOF, W(:,R0), W(:,R), CG(1), tcomm)
      call monolis_inner_product_R(monoCOM, monoMAT, NDOF, W(:,R0), W(:,W0),CG(2), tcomm)
      call monolis_inner_product_R(monoCOM, monoMAT, NDOF, W(:,R0), W(:,S), CG(3), tcomm)
      call monolis_inner_product_R(monoCOM, monoMAT, NDOF, W(:,R0), W(:,Z), CG(4), tcomm)
      call monolis_inner_product_R(monoCOM, monoMAT, NDOF, W(:,R),  W(:,R), CG(5), tcomm)

      call monolis_matvec(monoCOM, monoMAT, W(:,W0), W(:,T), tcomm)

      RR1= CG(1)
      RW = CG(2)
      RS = CG(3)
      RZ = CG(4)
      R2 = CG(5)

      call monolis_check_converge_2(monoPRM, monoCOM, monoMAT, R2, iter, is_converge, tcomm)
      if(is_converge) exit

      beta  = (alpha*RR1) / (RR*omega1)
      alpha = RR1 / (RW + beta*RS - beta*omega1*RZ)
      omega = omega1
      RR    = RR1
    enddo

    call monolis_update_R(monoCOM, NDOF, X, tcomm)

    deallocate(W)

    t2 = monolis_wtime()
    tsol = t2 - t1
  end subroutine monolis_solver_PipeBiCGSTAB_noprec
end module mod_monolis_solver_PipeBiCGSTAB_noprec