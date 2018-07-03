module VectorDouble
  use LinAlgParameters
  implicit none
  type :: DVec
    real(8), allocatable :: V(:)
  contains
    procedure :: Ini => iniV
    procedure :: Fin => FinV
    procedure :: prt => VectorPrint
    procedure :: Random => GetRandomVector
    procedure :: Nrm
    procedure :: Nrm2
  end type DVec
contains
  subroutine IniV(a, n)
  class(DVec), intent(inout) :: a
    integer, intent(in) :: n
    if(allocated(a%V)) deallocate(a%V)
    allocate(a%V(n))
    a%V(:) = 0.d0
  end subroutine IniV

  subroutine FinV(a)
  class(DVec), intent(inout) :: a
    if(allocated(a%V)) deallocate(a%V)
  end subroutine FinV

  subroutine VectorCopyD(b, a)
    type(DVec), intent(inout) :: b
    type(DVec), intent(in) :: a
    integer :: n
    n = size(a%V)
    call b%Ini(n)
    call dcopy(n, a%v, 1, b%v, 1)
  end subroutine VectorCopyD

  type(DVec) function VectorSumD(a, b) result(c)
    type(DVec), intent(in) :: a, b
    integer :: n
    if(size(a%v) /= size(b%v)) then
      write(*,'(a)') 'Error in DVectorSum'
      stop
    end if
    n = size(a%v)
    call VectorCopyD(c, a)
    call daxpy(n, 1.d0, b%v, 1, c%v, 1)
  end function VectorSumD

  type(DVec) function VectorSubtractD(a, b) result(c)
    type(DVec), intent(in) :: a, b
    integer :: n
    if(size(a%v) /= size(b%v)) then
      write(*,'(a)') 'Error in DVectorSubtract'
      stop
    end if
    n = size(a%v)
    call VectorCopyD(c, a)
    call daxpy(n, -1.d0, b%v, 1, c%v, 1)
  end function VectorSubtractD

  type(DVec) function VectorScaleRD(a, b) result(c)
    type(DVec), intent(in) :: a
    real(8), intent(in) :: b
    integer :: n
    n = size(a%v)
    call VectorCopyD(c, a)
    call dscal(n, b, c%v, 1)
  end function VectorScaleRD

  type(DVec) function VectorScaleLD(b, a) result(c)
    type(DVec), intent(in) :: a
    real(8), intent(in) :: b
    integer :: n
    n = size(a%v)
    call VectorCopyD(c, a)
    call dscal(n, b, c%v, 1)
  end function VectorScaleLD

  type(DVec) function VectorDivideD(a, b) result(c)
    type(DVec), intent(in) :: a
    real(8), intent(in) :: b
    integer :: n
    n = size(a%v)
    call VectorCopyD(c, a)
    call dscal(n, 1.d0 / b, c%v, 1)
  end function VectorDivideD

  real(8) function InnerProductD(a, b) result(c)
    type(DVec), intent(in) :: a, b
    integer :: n
    real(8) :: ddot
    if(size(a%v) /= size(b%v)) then
      write(*,'(a)') 'Error in InnerProduct'
      stop
    end if
    n = size(a%v)
    c = ddot(n, a%v, 1, b%v, 1)
  end function InnerProductD

  real(8) function Nrm(a) result(b)
  class(DVec), intent(in) :: a
    integer :: n
    real(8) :: dnrm2
    n = size(a%v)
    b = dnrm2(n, a%v, 1)
  end function Nrm

  real(8) function Nrm2(a) result(b)
  class(DVec), intent(in) :: a
    integer :: n
    real(8) :: ddot
    n = size(a%v)
    b = ddot(n, a%v, 1, a%v, 1)
  end function Nrm2

  subroutine VectorPrint(this, string)
  class(DVec),intent(in)::this
    integer :: n
    character(*), intent(in), optional :: string
    n = size(this%v, 1)
    write(*,*)
    if(present(string)) write(*,*) string
    write(*,'(10f10.4)') this%v(:)
  end subroutine VectorPrint

  subroutine GetRandomVector(v, n, dist)
    ! idist = 1: uniform (0, 1)
    ! idist = 2: uniform (-1, 1)
    ! idist = 3: normal (-1, 1)
  class(DVec), intent(inout) :: v
    integer, intent(in) :: n
    integer, intent(in), optional :: dist
    integer :: idist = 3
    if(present(dist)) idist = dist
    call v%ini(n)
    call dlarnv(idist, iseed, n, v%v)
  end subroutine GetRandomVector
end module VectorDouble
