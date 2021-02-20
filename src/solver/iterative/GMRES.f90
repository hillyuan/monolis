module mod_monolis_solver_GMRES
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

  subroutine monolis_solver_GMRES(monoPRM, monoCOM, monoMAT)
    implicit none
    type(monolis_prm) :: monoPRM
    type(monolis_com) :: monoCOM
    type(monolis_mat) :: monoMAT
    integer(kint), parameter :: R  = 1
    integer(kint), parameter :: ZP = R + 1
    integer(kint), parameter :: ZQ = R + 2
    integer(kint), parameter :: S  = R + 3
    integer(kint), parameter :: W  = S + 1
    integer(kint), parameter :: Y  = W
    integer(kint), parameter :: AV = Y  + 1
    integer(kint), parameter :: V  = AV + 1
    integer(kint) :: N, NP, NDOF, NNDOF, CS, SN, NRK
    integer(kint) :: i, iter, ik, irow, jj, k, kk, nrest
    real(kdouble), pointer :: B(:), X(:), SS(:)
    real(kdouble), pointer :: WW(:,:), H(:,:)
    real(kdouble) :: LDH, LDW, BNRM2, DNRM2, RNORM, RESID, tol
    real(kdouble) :: coef, val, VCS, VSN, DTEMP, AA, BB, R0, scale, RR
    logical :: is_converge

    N     = monoMAT%N
    NP    = monoMAT%NP
    NDOF  = monoMAT%NDOF
    NNDOF = N*NDOF
    X => monoMAT%X
    B => monoMAT%B

    NREST = 10
    if (NREST >= NDOF*NP-1) NREST = NDOF*NP-2

    NRK = NREST + 7
    allocate (H (NRK,NRK))
    allocate (WW(NDOF*NP,NRK))
    allocate (SS(NRK))

    LDH = NREST + 2
    LDW = N
    CS  = NREST + 1
    SN  = CS    + 1

    call monolis_residual(monoCOM, monoMAT, X, B, WW(:,R), monoPRM%tspmv, monoPRM%tcomm_spmv)
    call monolis_set_converge(monoPRM, monoCOM, monoMAT, B, BNRM2, is_converge, monoPRM%tdotp, monoPRM%tcomm_dotp)
    if(is_converge) return

    OUTER:do
      I = 0
      call monolis_inner_product_R(monoCOM, N, NDOF, WW(:,R), WW(:,R), DNRM2, monoPRM%tdotp, monoPRM%tcomm_dotp)
      if (DNRM2 == 0.d0) exit ! converged

      RNORM= dsqrt(DNRM2)
      coef= 1.0d0/RNORM
      do ik= 1, NNDOF
        WW(ik,V)= WW(ik,R) * coef
      enddo

      WW(1 ,S) = RNORM
      do k = 2, NNDOF
        WW(k,S) = 0.0d0
      enddo

      do I = 1, NREST
        ITER= ITER + 1

        call monolis_precond_apply(monoPRM, monoCOM, monoMAT, WW(:,V+I-1), WW(:,ZQ))

        call monolis_matvec(monoCOM, monoMAT, WW(:,ZQ), WW(:,W), monoPRM%tspmv, monoPRM%tcomm_spmv)

        do K = 1, I
          call monolis_inner_product_R(monoCOM, N, NDOF, WW(:,W), WW(:,V+K-1), val, monoPRM%tdotp, monoPRM%tcomm_dotp)

          do ik= 1, NNDOF
            WW(ik,W)= WW(ik,W) - val * WW(ik,V+K-1)
          enddo
          H(K,I)= val
        enddo

        call monolis_inner_product_R(monoCOM, N, NDOF, WW(:,W), WW(:,W), val, monoPRM%tdotp, monoPRM%tcomm_dotp)
        if (val == 0.d0) exit ! converged

        H(I+1,I)= dsqrt(val)
        coef= 1.0d0 / H(I+1,I)
        do ik= 1, NNDOF
          WW(ik,V+I+1-1)= WW(ik,W) * coef
        enddo

        do k = 1, I-1
          VCS= H(k,CS)
          VSN= H(k,SN)
          DTEMP   = VCS*H(k  ,I) + VSN*H(k+1,I)
          H(k+1,I)= VCS*H(k+1,I) - VSN*H(k  ,I)
          H(k  ,I)= DTEMP
        enddo

        AA = H(I  ,I)
        BB = H(I+1,I)
        R0= BB
        if (dabs(AA).gt.dabs(BB)) R0= AA
        scale= dabs(AA) + dabs(BB)

        if (scale.ne.0.d0) then
          RR= scale * dsqrt((AA/scale)**2+(BB/scale)**2)
          RR= dsign(1.d0,R0)*RR
          H(I,CS)= AA/RR
          H(I,SN)= BB/RR
        else
          H(I,CS)= 1.d0
          H(I,SN)= 0.d0
          RR     = 0.d0
        endif

        VCS= H(I,CS)
        VSN= H(I,SN)
        DTEMP    = VCS*H(I  ,I) + VSN*H(I+1,I)
        H (I+1,I)= VCS*H(I+1,I) - VSN*H(I  ,I)
        H (I  ,I)= DTEMP

        DTEMP    = VCS*WW(I  ,S) + VSN*WW(I+1,S)
        WW(I+1,S)= VCS*WW(I+1,S) - VSN*WW(I  ,S)
        WW(I  ,S)= DTEMP

        RESID = dabs ( WW(I+1,S))/dsqrt(BNRM2)

!        if (my_rank.eq.0 .and. ITERlog.eq.1)                              &
!          &    write (*, '(2i8, 1pe16.6)') iter,I+1, RESID

        if ( RESID.le.TOL ) then
          !C-- [H]{y}= {s_tld}
          do ik= 1, I
            SS(ik)= WW(ik,S)
          enddo
          IROW= I
          WW(IROW,Y)= SS(IROW) / H(IROW,IROW)

          do kk= IROW-1, 1, -1
            do jj= IROW, kk+1, -1
              SS(kk)= SS(kk) - H(kk,jj)*WW(jj,Y)
            enddo
            WW(kk,Y)= SS(kk) / H(kk,kk)
          enddo

          !C-- {x}= {x} + {y}{V}
          do kk= 1, NNDOF
            WW(kk, AV)= 0.d0
          enddo

          jj= IROW
          do jj= 1, IROW
            do kk= 1, NNDOF
              WW(kk,AV)= WW(kk,AV) + WW(jj,Y)*WW(kk,V+jj-1)
            enddo
          enddo

          call monolis_precond_apply(monoPRM, monoCOM, monoMAT, WW(:,AV), WW(:,ZQ))

          do kk= 1, NNDOF
            X(kk)= X(kk) + WW(kk,ZQ)
          enddo

          exit OUTER
        endif

!        if ( ITER.gt.MAXIT ) then
!          exit OUTER
!        end if
      end do
      !C===

      !C-- [H]{y}= {s_tld}
      do ik= 1, NREST
        SS(ik)= WW(ik,S)
      enddo
      IROW= NREST
      WW(IROW,Y)= SS(IROW) / H(IROW,IROW)

      do kk= IROW-1, 1, -1
        do jj= IROW, kk+1, -1
          SS(kk)= SS(kk) - H(kk,jj)*WW(jj,Y)
        enddo
        WW(kk,Y)= SS(kk) / H(kk,kk)
      enddo

      !C-- {x}= {x} + {y}{V}
      do kk= 1, NNDOF
        WW(kk, AV)= 0.d0
      enddo

      jj= IROW
      do jj= 1, IROW
        do kk= 1, NNDOF
          WW(kk,AV)= WW(kk,AV) + WW(jj,Y)*WW(kk,V+jj-1)
        enddo
      enddo

      call monolis_precond_apply(monoPRM, monoCOM, monoMAT, WW(:,AV), WW(:,ZQ))

      do kk= 1, NNDOF
        X(kk)= X(kk) + WW(kk,ZQ)
      enddo

      !C-- Compute residual vector R, find norm, then check for tolerance.
      call monolis_residual(monoCOM, monoMAT, X, B, WW(:,R), monoPRM%tspmv, monoPRM%tcomm_spmv)

      call monolis_inner_product_R(monoCOM, N, NDOF, WW(:,R), WW(:,R), DNRM2, monoPRM%tdotp, monoPRM%tcomm_dotp)

      WW(I+1,S)= dsqrt(DNRM2/BNRM2)
      RESID    = WW( I+1,S )

      if ( RESID.le.TOL )   exit OUTER
    end do OUTER

    call monolis_update_R(monoCOM, NDOF, X, monoPRM%tcomm_spmv)

  end subroutine monolis_solver_GMRES

end module mod_monolis_solver_GMRES
