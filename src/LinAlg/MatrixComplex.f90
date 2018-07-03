module MatrixComplex
  implicit none
  type :: CMat
    complex(8), allocatable :: M(:,:)
  contains
    procedure :: Ini => IniM
    procedure :: Fin => FinM
    procedure :: eye
    procedure :: T => Trans
    procedure :: C => ComplexConjugate
    procedure :: H => HermiteConjugate
    procedure :: Inv => Inverse
    procedure :: Det
    procedure :: Random => GetRandomMatrix
    procedure :: prt => MatrixPrint
    procedure :: DiagMat
  end type CMat
contains
  subroutine iniM(a, m, n)
  class(CMat), intent(inout) :: a
    integer, intent(in) :: m,n
    if(allocated(a%m)) deallocate(a%m)
    allocate(a%m(m,n))
    a%m = 0.d0
  end subroutine iniM

  subroutine eye(a, n)
  class(CMat), intent(inout) :: a
    integer, intent(in) :: n
    integer :: i
    if(allocated(a%m)) deallocate(a%m)
    allocate(a%m(n,n))
    a%m = 0.d0
    do i = 1, n
      a%m(i,i) = (1.d0, 0.d0)
    end do
  end subroutine eye

  subroutine FinM(a)
  class(CMat), intent(inout) :: a
    if(allocated(a%m)) deallocate(a%m)
  end subroutine FinM

  subroutine MatrixCopyC(b, a)
    type(CMat), intent(inout) :: b
    type(CMat), intent(in) :: a
    integer :: m, n, i
    m = size(a%m, 1)
    n = size(a%m, 2)
    call b%Ini(m,n)
    do i = 1, n
      call zcopy(m, a%m(:,i), 1, b%m(:,i), 1)
    end do
  end subroutine MatrixCopyC

  type(CMat) function MatrixProductC(a, b) result(c)
    type(CMat), intent(in) :: a, b
    integer :: m, k, n
    m = size(a%m, 1)
    k = size(a%m, 2)
    if(size(a%m, 2) /= size(b%m, 1)) then
      write(*, '(a)') 'Error in MatrixProduct'
      stop
    end if
    n = size(b%m, 2)
    call c%Ini(m,n)
    call zgemm('n','n',m,n,k,1.d0,a%m,m,b%m,k,0.d0,c%m,m)
  end function MatrixProductC

  type(CMat) function MatrixSumC(a, b) result(c)
    type(CMat), intent(in) :: a, b
    integer :: m, n, i
    if(size(a%m, 1) /= size(b%m, 1) .or. size(a%m, 2) /= size(b%m, 2)) then
      write(*, '(a)') 'Error in MatrixSum'
      stop
    end if
    m = size(a%m, 1)
    n = size(a%m, 2)
    call MatrixCopyC(c, a)
    do i = 1, n
      call zaxpy(m, 1.d0, b%m(:,i), 1, c%m(:,i), 1)
    end do
  end function MatrixSumC

  type(CMat) function MatrixSubtractC(a, b) result(c)
    type(CMat), intent(in) :: a, b
    integer :: m, n, i
    if(size(a%m, 1) /= size(b%m, 1) .or. size(a%m, 2) /= size(b%m, 2)) then
      write(*, '(a)') 'Error in MatrixSum'
      stop
    end if
    m = size(a%m, 1)
    n = size(a%m, 2)
    call MatrixCopyC(c, a)
    do i = 1, n
      call zaxpy(m, -1.d0, b%m(:,i), 1, c%m(:,i), 1)
    end do
  end function MatrixSubtractC

  type(CMat) function MatrixScaleLC(b, a) result(c)
    type(CMat), intent(in) :: b
    complex(8), intent(in) :: a
    integer :: m, n, i
    m = size(b%m, 1)
    n = size(b%m, 2)
    call MatrixCopyC(c, b)
    do i = 1, n
      call dscal(m, a, c%m(:,i), 1)
    end do
  end function MatrixScaleLC


  type(CMat) function MatrixScaleRC(a, b) result(c)
    type(CMat), intent(in) :: b
    complex(8), intent(in) :: a
    integer :: m, n, i
    m = size(b%m, 1)
    n = size(b%m, 2)
    call MatrixCopyC(c, b)
    do i = 1, n
      call dscal(m, a, c%m(:,i), 1)
    end do
  end function MatrixScaleRC

  type(CMat) function MatrixScaleDivideC(b, a) result(c)
    type(CMat), intent(in) :: b
    real(8), intent(in) :: a
    integer :: m, n, i
    m = size(b%m, 1)
    n = size(b%m, 2)
    call MatrixCopyC(c, b)
    do i = 1, n
      call dscal(m, 1.d0 / a, c%m(:,i), 1)
    end do
  end function MatrixScaleDivideC

  type(CMat) function Trans(a) result(b)
    class(CMat), intent(in) :: a
    integer :: n, m
    m = size(a%m, 1)
    n = size(a%m, 2)
    call b%Ini(n,m)
    b%M = transpose(a%M)
  end function Trans

  type(CMat) function ComplexConjugate(a) result(b)
  class(CMat), intent(in) :: a
    integer :: n, m
    m = size(a%m, 1)
    n = size(a%m, 2)
    call b%Ini(n,m)
    b%M = conjg(a%M)
  end function ComplexConjugate

  type(CMat) function HermiteConjugate(a) result(b)
  class(CMat), intent(in) :: a
    b = a%C()
    b = b%T()
  end function HermiteConjugate

  type(CMat) function inverse(r) result(s)
  class(CMat), intent(in) :: r
    complex(8), allocatable :: a(:,:)
    complex(8), allocatable :: work(:)
    integer, allocatable :: ipvt(:)
    integer :: info, n
    n = size(r%m, 1)
    call s%Ini(n,n)
    allocate(work(n*n),ipvt(n))
    allocate(a(n,n))
    a = r%m
    call zgetrf(n,n,a,n,ipvt,info)
    call zgetri(n,a,n,ipvt,work,n**2,info)
    s%m = a
    deallocate(a,work,ipvt)
  end function inverse

  complex(8) function Det(r) result(d)
  class(CMat), intent(in) :: r
    integer :: n, i, info
    complex(8), allocatable :: a(:,:)
    integer, allocatable :: ipiv(:)
    n = size(r%m, 1)
    allocate(ipiv(n), a(n,n))
    a = r%m
    call zgetrf(n, n, a, n, ipiv, info)
    if(info /= 0) then
      write(*,'(a, i3)') "error in det: info = ", info
      stop
    end if
    d = 1.d0
    do i = 1, n
      if(ipiv(i) .ne. i) then
        d = -d * a(i, i)
      else
        d = d * a(i, i)
      end if
    end do
    deallocate(ipiv, a)
  end function Det

  subroutine MatrixPrint(this, string)
  class(CMat), intent(in) :: this
    character(12) :: cfmt
    integer :: i, n, m
    character(*), intent(in), optional :: string
    cfmt = '( xf10.4)'
    n = size(this%m, 1)
    m = size(this%m, 2)
    write(cfmt(2:3), '(I2)') m
    write(*,*)
    if(present(string)) write(*,*) string
    write(*,'(a)') 'Real:'
    do i=1,n
      write(*,cfmt) dble(this%m(i,:))
    end do
    write(*,'(a)') 'Imag:'
    do i=1,n
      write(*,cfmt) dimag(this%m(i,:))
    end do
  end subroutine MatrixPrint

  subroutine GetRandomMatrix(mat, m, n)
    use VectorComplex, only: CVec
  class(CMat), intent(inout) :: mat
    integer, intent(in) :: m, n
    integer :: i
    type(CVec) :: v
    call mat%ini(m,n)
    do i = 1, n
      call v%Random(m)
      mat%m(:,i) = v%v(:)
      call v%fin()
    end do
  end subroutine GetRandomMatrix

  subroutine DiagMat(b, a)
    use VectorComplex, only: CVec
  class(CMat), intent(inout) :: b
    type(CVec), intent(in) :: a
    integer :: n, i
    n = size(a%V)
    call b%Ini(n,n)
    do i = 1, n
      b%M(i,i) = a%V(i)
    end do
  end subroutine DiagMat
end module MatrixComplex
