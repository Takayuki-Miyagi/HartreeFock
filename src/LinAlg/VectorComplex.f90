module VectorComplex
  use LinAlgParameters
  implicit none
  type :: CVec
    complex(8), allocatable :: V(:)
  contains
    procedure :: Ini => iniV
    procedure :: Fin => FinV
    procedure :: CC => ComplexConjugate
    procedure :: prt => VectorPrint
    procedure :: Random => GetRandomVector
    procedure :: Nrm
    procedure :: Nrm2
  end type CVec
contains
  subroutine IniV(a, n)
  class(CVec), intent(inout) :: a
    integer, intent(in) :: n
    if(allocated(a%V)) deallocate(a%V)
    allocate(a%V(n))
    a%V(:) = 0.d0
  end subroutine IniV

  subroutine FinV(a)
  class(CVec), intent(inout) :: a
    if(allocated(a%V)) deallocate(a%V)
  end subroutine FinV

  subroutine VectorCopyC(b, a)
    type(CVec), intent(inout) :: b
    type(CVec), intent(in) :: a
    integer :: n
    n = size(a%V)
    call b%Ini(n)
    call zcopy(n, a%v, 1, b%v, 1)
  end subroutine VectorCopyC

  type(CVec) function VectorSumC(a, b) result(c)
    type(CVec), intent(in) :: a, b
    integer :: n
    if(size(a%v) /= size(b%v)) then
      write(*,'(a)') 'Error in DVectorSum'
      stop
    end if
    n = size(a%v)
    call VectorCopyC(c,a)
    call zaxpy(n, 1.d0, b%v, 1, c%v, 1)
  end function VectorSumC

  type(CVec) function VectorSubtractC(a, b) result(c)
    type(CVec), intent(in) :: a, b
    integer :: n
    if(size(a%v) /= size(b%v)) then
      write(*,'(a)') 'Error in DVectorSubtract'
      stop
    end if
    n = size(a%v)
    call VectorCopyC(c,a)
    call daxpy(n, -1.d0, b%v, 1, c%v, 1)
  end function VectorSubtractC

  type(CVec) function VectorScaleRC(a, b) result(c)
    type(CVec), intent(in) :: a
    complex(8), intent(in) :: b
    integer :: n
    n = size(a%v)
    call VectorCopyC(c,a)
    call zscal(n, b, c%v, 1)
  end function VectorScaleRC


  type(CVec) function VectorScaleLC(b, a) result(c)
    type(CVec), intent(in) :: a
    complex(8), intent(in) :: b
    integer :: n
    n = size(a%v)
    call VectorCopyC(c,a)
    call dscal(n, b, c%v, 1)
  end function VectorScaleLC

  type(CVec) function VectorDivideC(a, b) result(c)
    type(CVec), intent(in) :: a
    real(8), intent(in) :: b
    integer :: n
    n = size(a%v)
    call VectorCopyC(c,a)
    call zdscal(n, 1.d0 / b, c%v, 1)
  end function VectorDivideC

  type(CVec) function ComplexConjugate(a) result(b)
  class(CVec), intent(in) :: a
    integer :: n
    n = size(a%v)
    call b%ini(n)
    b%v = conjg(a%v)
  end function ComplexConjugate

  complex(8) function InnerProductC(a, b) result(c)
    type(CVec), intent(in) :: a, b
    integer :: n
    complex(8) :: zdotu
    if(size(a%v) /= size(b%v)) then
      write(*,'(a)') 'Error in InnerProduct'
      stop
    end if
    n = size(a%v)
    c = zdotu(n, a%v, 1, b%v, 1)
  end function InnerProductC

  real(8) function Nrm(a) result(b)
  class(CVec), intent(in) :: a
    integer :: n
    real(8) :: dznrm2
    n = size(a%v)
    b = dznrm2(n, a%v, 1)
  end function Nrm

  real(8) function Nrm2(a) result(b)
  class(CVec), intent(in) :: a
    integer :: n
    real(8) :: dznrm2
    n = size(a%v)
    b = dznrm2(n, a%v, 1) ** 2
  end function Nrm2

  subroutine VectorPrint(this, string)
  class(CVec),intent(in)::this
    integer :: n
    character(*), intent(in), optional :: string
    n = size(this%v, 1)
    write(*,*)
    if(present(string)) write(*,*) string
    write(*,'(a)',advance='no') 'Real:'
    write(*,'(10f10.4)') dble(this%v(:))
    write(*,'(a)',advance='no') 'Imag:'
    write(*,'(10f10.4)') dimag(this%v(:))
    !write(*,*) this%v(:)
  end subroutine VectorPrint

  subroutine GetRandomVector(v, n, dist)
    ! idist = = 1:  real and imaginary parts each uniform (0,1)
    ! idist = = 2:  real and imaginary parts each uniform (-1,1)
    ! idist = = 3:  real and imaginary parts each normal (0,1)
    ! idist = = 4:  uniformly distributed on the disc abs(z) < 1
    ! idist = = 5:  uniformly distributed on the circle abs(z) = 1
  class(CVec), intent(inout) :: v
    integer, intent(in) :: n
    integer, intent(in), optional :: dist
    integer :: idist = 4
    if(present(dist)) idist = dist
    call v%ini(n)
    call zlarnv(idist, iseed, n, v%v)
  end subroutine GetRandomVector
end module VectorComplex
