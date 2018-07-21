module Optimizer
  use MatrixDouble, only: DMat
  use VectorDouble, only: DVec
  use LinAlgLib
  implicit none
  private :: InitFixedPointProblem, FinFixedPointProblem, GetNextInput, DirectIteration, &
      & BroydenMethod, ModifiedBroydenMethod
  public :: FixedPointProblem
  type :: FixedPointProblem
    integer :: n = 0
    integer :: nite = 0
    real(8) :: si = 1.d0
    real(8) :: alpha = 1.0d0! step size
    integer :: m = 20
    type(DMat) :: Bj
    type(DVec) :: xj, fj
    type(DMat) :: Bi
    type(DVec) :: xi, fi
    type(DVec), allocatable :: df(:), dx(:)
    character(256) :: Method = 'direct'
  contains
    procedure :: init => InitFixedPointProblem
    procedure :: fin => FinFixedPointProblem
    procedure :: GetVector => GetNextInput
    procedure :: DirectIteration
    procedure :: BroydenMethod
    procedure :: ModifiedBroydenMethod
  end type FixedPointProblem
contains
  subroutine InitFixedPointProblem(this, n, alpha, method, m)
    class(FixedPointProblem), intent(inout) :: this
    real(8), optional, intent(in) :: alpha
    character(*), optional, intent(in) :: method
    integer, optional, intent(in) :: m
    integer, intent(in) :: n
    integer :: i
    this%n = n
    if(present(alpha)) this%alpha = alpha
    if(present(method)) this%method = method
    if(present(m)) this%m = m

    call this%xi%Ini(n)
    call this%fi%Ini(n)
    call this%xj%Ini(n)
    call this%fj%Ini(n)
    if(   this%method == 'broyden' .or. &
        & this%method == 'BFGS') then
      call this%Bi%Ini(n,n)
      call this%Bj%Ini(n,n)
      do i = 1, n
        this%Bi%m(i,i) =  - this%alpha
        this%Bj%m(i,i) =  - this%alpha
      end do
    end if
    if(   this%method == 'mbroyden' .or. &
        & this%method == 'LBFGS') then
      allocate(this%df(this%m))
      allocate(this%dx(this%m))
      do i = 1, this%m
        call this%df(i)%ini(n)
        call this%dx(i)%ini(n)
      end do
    end if
  end subroutine InitFixedPointProblem

  subroutine FinFixedPointProblem(this)
    class(FixedPointProblem), intent(inout) :: this
    integer :: i
    call this%xi%Fin()
    call this%fi%Fin()
    call this%Bi%Fin()
    call this%xj%Fin()
    call this%fj%Fin()
    call this%Bj%Fin()
    if(   this%method == 'mbroyden' .or. &
        & this%method == 'LBFGS') then
      do i = 1, this%m
        call this%df(i)%Fin()
        call this%dx(i)%Fin()
      end do
    end if
  end subroutine FinFixedPointProblem

  subroutine GetNextInput(this)
    class(FixedPointProblem), intent(inout) :: this

    if(this%Method == 'direct') call this%DirectIteration()
    if(this%Method == 'broyden') call this%BroydenMethod()
    if(this%Method == 'mbroyden') call this%ModifiedBroydenMethod()
  end subroutine GetNextInput

  subroutine DirectIteration(this)
    class(FixedPointProblem), intent(inout) :: this
    this%si = maxval(abs(this%fj%v - this%xj%v))
    this%xj = this%fj
  end subroutine DirectIteration

  subroutine BroydenMethod(this)
    ! find x* so that x* = f(x*)
    ! x_{j+1} = x_{j} - alpha B_{j} f_{j}
    ! B_{j+1} = B_{j} + (|dx> - Bj|df>) <dx|Bj / <dx|Bj|df>
    ! xj : input of j-th iteration
    ! fj : f(xj)
    ! Bj : Jacobian of j-th iteration
    class(FixedPointProblem), intent(inout) :: this
    type(DVec) :: df, dx
    real(8) :: a
    this%fj = (this%fj - this%xj)
    this%si = maxval(abs(this%fj%v))
    if(this%nite == 1) then
      this%xj = this%xi - (this%Bi * this%fj) * this%alpha
      this%fi = this%fj
      return
    else
      df = this%fj - this%fi
      dx = this%xj - this%xi
      a = 1.d0 / (dx * (this%Bi * df))
      this%Bj = this%Bi + ((dx - (this%Bi * df)) .x. (dx * this%Bi)) * a
      this%xi = this%xj
      this%fi = this%fj
      this%xj = this%xi - (this%Bj * this%fj)
      this%Bi = this%Bj
      call df%Fin()
      call dx%Fin()
      return
    end if
  end subroutine BroydenMethod

  subroutine ModifiedBroydenMethod(this)
    ! find x* so that x* = f(x*)
    ! x_{j+1} = x_{j} - alpha B_{j} f_{j}
    ! xj : input of j-th iteration
    ! fj : f(xj)
    ! Bj : Jacobian of j-th iteration
    class(FixedPointProblem), intent(inout) :: this
    type(DMat) :: beta, mat, u
    type(DVec) :: vec, gamm, c
    real(8) :: a, curv, w0
    integer :: iused, inext, icurr, i, j

    this%fj = (this%fj - this%xj)
    this%si = maxval(abs(this%fj%v))

    iused = min(this%nite - 1, this%m)
    icurr = this%nite - 1 - ((this%nite - 2) / this%m) * this%m
    inext = this%nite - ((this%nite - 1) / this%m) * this%m
    if(this%m == 0 .or. iused < 1) then
      this%xj = this%xj + this%fj * this%alpha
      this%df(inext) = this%fj
      this%dx(inext) = this%xj
      return
    end if
    w0 = 1.d-2
    this%df(icurr) = this%fj - this%df(icurr)
    this%dx(icurr) = this%xj - this%dx(icurr)
    a = 1.d0 / dsqrt((this%df(icurr) * this%df(icurr)))
    this%df(icurr) = this%df(icurr) * a
    this%dx(icurr) = this%dx(icurr) * a

    call mat%ini(iused, iused)
    call beta%ini(iused, iused)
    do i = 1, iused
      do j = 1, iused
        mat%m(i,j) = (this%df(i) * this%df(j))
      end do
      mat%m(i,i) = 1.d0 + w0*w0
    end do
    if(iused > 0) beta = mat%inv()
    call c%ini(iused)
    call u%ini(this%n, iused)
    do i = 1, iused
      c%v(i) = this%df(i) * this%fj
      u%m(:,i) = (this%df(i)%v(:) * this%alpha) + this%dx(i)%v(:)
    end do
    gamm = c * beta
    vec = (this%fj * this%alpha) - (u * gamm)
    curv = this%fj * vec
    this%df(inext) = this%fj
    this%dx(inext) = this%xj
    this%xj = this%xj + vec
    call mat%Fin()
    call vec%Fin()
    call beta%Fin()
    call gamm%Fin()
  end subroutine ModifiedBroydenMethod
end module Optimizer
