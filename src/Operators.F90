module Operators
  use omp_lib
  use ModelSpace
  use NOperators
  implicit none

  type :: Op
    character(32) :: optr
    real(8) :: zero
    type(NBodyPart) :: one, two
#ifdef single_precision
    type(NBodyPartSp) :: thr
#else
    type(NBodyPart) :: thr
#endif
    logical :: is_nomal_ordred = .false.
    logical :: Scalar
    logical :: is_three_body = .false.
    integer :: jr, pr, zr
  contains
    procedure :: fin => FinOp
    procedure :: InitOp
    procedure :: InitOpFromString
    generic :: init => InitOp, InitOpFromString
  end type Op

  interface assignment(=)
    procedure :: CopyOp
  end interface assignment(=)

  interface operator(+)
    procedure :: SumOp
  end interface operator(+)

  interface operator(-)
    procedure :: SubtractOp
  end interface operator(-)

  interface operator(*)
    procedure :: ScaleLeftOp, ScaleRightOp
  end interface operator(*)
contains
  subroutine FinOp(this)
    class(Op), intent(inout) :: this
    call this%one%fin()
    call this%two%fin()
    if(.not. this%is_three_body) return
    call this%thr%fin()
  end subroutine FinOp

  subroutine InitOpFromString(this, optr, ms, is_three)
    use DefineOperators, only: GetOperatorRank
    class(Op), intent(inout) :: this
    type(MSpace), intent(in) :: ms
    character(*), intent(in) :: optr
    logical, intent(in) :: is_three
    integer :: jr, pr, zr

    call GetOperatorRank(optr,jr,pr,zr)
    call this%init(jr,pr,zr,optr,ms,is_three)
  end subroutine InitOpFromString

  subroutine InitOp(this, jr, pr, zr, optr, ms, is_three)
    use Profiler, only: timer
    class(Op), intent(inout) :: this
    integer, intent(in) :: jr, pr, zr
    character(*), intent(in) :: optr
    type(MSpace), intent(in) :: ms
    logical, intent(in) :: is_three
    real(8) :: ti

    ti = omp_get_wtime()
    call timer%cmemory()

    this%jr = jr
    this%pr = pr
    this%zr = zr
    this%optr = optr
    if(this%jr == 0 .and. this%pr == 1 .and. this%zr == 0) this%Scalar = .true.
    if(this%jr /= 0 .or.  this%pr /= 1 .or.  this%zr /= 0) this%Scalar = .false.

    call this%one%init(ms%one, this%Scalar, optr, jr, pr, zr)
    call this%two%init(ms%two, this%Scalar, optr, jr, pr, zr)
    if(.not. is_three) then
      call timer%countup_memory(trim(optr))
      call timer%Add('Construct '//trim(optr), omp_get_wtime()-ti)
      return
    end if
    this%is_three_body = .true.
    call this%thr%init(ms%thr, this%Scalar, optr, jr, pr, zr)

    call timer%countup_memory(trim(optr))
    call timer%Add('Construct '//trim(optr), omp_get_wtime()-ti)
  end subroutine InitOp

  subroutine CopyOp(a, b)
    type(Op), intent(inout) :: a
    type(Op), intent(in) :: b

    a%optr = b%optr
    a%Scalar = b%Scalar
    a%is_three_body = b%is_three_body
    a%jr = b%jr
    a%pr = b%pr
    a%zr = b%zr
    a%zero = b%zero
    a%one = b%one
    a%two = b%two
    if(.not. a%is_three_body) return
    a%thr = b%thr
  end subroutine CopyOp

  function SumOp(a, b) result(c)
    type(Op), intent(in) :: a, b
    type(Op) :: c

    c = a
    c%zero = a%zero + b%zero
    c%one = a%one + b%one
    c%two = a%two + b%two
    if(.not. a%is_three_body) return
    c%thr = a%thr + b%thr
  end function SumOp

  function SubtractOp(a, b) result(c)
    type(Op), intent(in) :: a, b
    type(Op) :: c

    c = a
    c%zero = a%zero - b%zero
    c%one = a%one - b%one
    c%two = a%two - b%two
    if(.not. a%is_three_body) return
    c%thr = a%thr - b%thr
  end function SubtractOp

  function ScaleLeftOp(a, b) result(c)
    type(Op), intent(in) :: a
    real(8), intent(in) :: b
    type(Op) :: c

    c = a
    c%zero = a%zero * b
    c%one = a%one * b
    c%two = a%two * b
    if(.not. a%is_three_body) return
    c%thr = a%thr * b
  end function ScaleLeftOp

  function ScaleRightOp(a, b) result(c)
    type(Op), intent(in) :: b
    real(8), intent(in) :: a
    type(Op) :: c

    c = b
    c%zero = b%zero * a
    c%one = b%one * a
    c%two = b%two * a
    if(.not. b%is_three_body) return
    c%thr = b%thr * a
  end function ScaleRightOp

  subroutine SetOperatorFromFile(this, ms, file_nn, file_3n, &
        & bound_2b_file, bound_3b_file)
    class(Op), intent(inout) :: this
    type(MSpace), intent(in) :: ms
    character(*), intent(in) :: file_nn, file_3n
    integer, intent(in), optional :: bound_2b_file(3), bound_3b_file(4)
    type(NBodyPart) :: Tcm_one, Tcm_two
    type(NBodyPart) :: Hcm_one, Hcm_two
    type(ReadFiles) :: rd

    if(file_nn == 'none' .and. file_3n == 'none') return
    ! -- boundary for two-body file
    rd%file_nn = file_nn
    rd%emax2 = ms%emax
    rd%e2max2 = ms%e2max
    rd%lmax2 = ms%lmax
    if(present(bound_2b_file)) then
      rd%emax2 = bound_2b_file(1)
      rd%e2max2= bound_2b_file(2)
      rd%lmax2 = bound_2b_file(3)
    end if

    ! -- boundary for three-body file
    rd%file_3n = file_3n
    rd%emax3 = ms%emax
    rd%e2max3 = ms%e2max
    rd%e3max3 = ms%e3max
    rd%lmax3 = ms%lmax
    if(present(bound_2b_file)) then
      rd%emax3 = bound_3b_file(1)
      rd%e2max3= bound_3b_file(2)
      rd%e3max3= bound_3b_file(3)
      rd%lmax3 = bound_3b_file(4)
    end if

    call this%one%set(ms%one, ms%sps, ms%hw, ms%A, ms%Z, ms%N)
    call this%two%ReadFile(ms%sps, ms%two, rd)
    if(this%is_three_body) then
      call this%thr%ReadFile(ms%isps,ms%thr, rd)
    end if

    select case(this%optr)
    case('hamil', 'Hamil')
      call Tcm_one%init(ms%one, .true., 'Tcm', 0, 1, 0)
      call Tcm_two%init(ms%two, .true., 'Tcm', 0, 1, 0)

      call Tcm_one%set(ms%one, ms%sps, ms%hw, ms%A, ms%Z, ms%N)
      call Tcm_two%set(ms%two, ms%sps, ms%hw, ms%A, ms%Z, ms%N)

      this%one = this%one - Tcm_one
      this%two = this%two - Tcm_two

      if(ms%beta > 1.d-4) then
        call Hcm_one%init(ms%one, .true., 'CMhamil', 0, 1, 0)
        call Hcm_two%init(ms%two, .true., 'CMhamil', 0, 1, 0)

        call Hcm_one%set(ms%one, ms%sps, ms%hw, ms%A, ms%Z, ms%N)
        call Hcm_two%set(ms%two, ms%sps, ms%hw, ms%A, ms%Z, ms%N)

        this%zero = this%zero - (1.5d0 * ms%hw * ms%beta)
        this%one = this%one + (Hcm_one * ms%beta)
        this%two = this%two + (Hcm_two * ms%beta)
      end if
    case default
      return
    end select
  end subroutine SetOperatorFromFile

  subroutine SetOperator(this, ms)
    class(Op), intent(inout) :: this
    type(MSpace), intent(in) :: ms

    call this%one%set(ms%one, ms%sps, ms%hw, ms%A, ms%Z, ms%N)
    call this%two%set(ms%two, ms%sps, ms%hw, ms%A, ms%Z, ms%N)
  end subroutine SetOperator

end module Operators
