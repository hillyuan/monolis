module mod_monolis_prm
  implicit none

  integer(4), parameter :: kint    = 4
  integer(4), parameter :: kdouble = 8

  integer(4), parameter :: monolis_iter_CG       = 1
  integer(4), parameter :: monolis_iter_GropCG   = 2
  integer(4), parameter :: monolis_iter_PipeCG   = 3
  integer(4), parameter :: monolis_iter_PipeCR   = 4
  integer(4), parameter :: monolis_iter_BiCGSTAB = 5
  integer(4), parameter :: monolis_iter_PipeBiCGSTAB = 6
  integer(4), parameter :: monolis_iter_BiCGSTAB_noprec = 7
  integer(4), parameter :: monolis_iter_CABiCGSTAB_noprec = 8
  integer(4), parameter :: monolis_iter_PipeBiCGSTAB_noprec = 9
  integer(4), parameter :: monolis_iter_SOR      = 10
  integer(4), parameter :: monolis_iter_IR       = 11

  integer(4), parameter :: monolis_prec_DIAG   = 1
  integer(4), parameter :: monolis_prec_ILU    = 2
  integer(4), parameter :: monolis_prec_JACOBI = 3
  integer(4), parameter :: monolis_prec_SOR    = 4
  integer(4), parameter :: monolis_prec_SAINV  = 5
  integer(4), parameter :: monolis_prec_RIF    = 6
  integer(4), parameter :: monolis_prec_SPIKE  = 7
  integer(4), parameter :: monolis_prec_DIRECT = 8

  character*24, dimension(11) :: monolis_str_iter = (/&
  & "CG                 ", &
  & "GropCG             ", &
  & "PipeCG             ", &
  & "PipeCR             ", &
  & "BiCGSTAB           ", &
  & "PipeBiCGSTAB       ", &
  & "BiCGSTAB_noprec    ", &
  & "CABiCGSTAB_noprec  ", &
  & "PipeBiCGSTAB_noprec", &
  & "SOR                ", &
  & "IR                 "/)
  character*24, dimension(8)  :: monolis_str_prec = (/&
  & "Diag  ", &
  & "ILU   ", &
  & "Jacobi", &
  & "SOR   ", &
  & "SAINV ", &
  & "RIF   ", &
  & "SPIKE ", &
  & "Direct"/)

  type monolis_prm
    integer(kind=kint) :: method = 1
    integer(kind=kint) :: precond = 1
    integer(kind=kint) :: maxiter = 1000
    real(kind=kdouble) :: tol = 1.0d-8
    logical :: is_scaling    = .true.
    logical :: is_reordering = .true.
    logical :: is_init_x = .true.
    logical :: show_iteration = .true.
  end type monolis_prm

contains

  subroutine monolis_prm_initialize(monoPRM)
    implicit none
    type(monolis_prm) :: monoPRM

    monoPRM%method = 1
    monoPRM%precond = 1
    monoPRM%maxiter = 1000
    monoPRM%tol = 1.0d-8
    monoPRM%is_scaling = .true.
    monoPRM%is_reordering = .true.
    monoPRM%is_init_x = .true.
    monoPRM%show_iteration = .true.
  end subroutine monolis_prm_initialize

  subroutine monolis_prm_finalize(monoPRM)
    implicit none
    type(monolis_prm) :: monoPRM

  end subroutine monolis_prm_finalize

end module mod_monolis_prm