module Operators
  use omp_lib
  use LinAlgLib
  use ModelSpace
  implicit none

  type, extends(DMat) :: NBodyChannel
    integer :: ndims(2)
    logical :: is = .false.
  contains
  end type NBodyChannel

  type :: NBodyPart
    type(NBodyChannel), allocatable :: MatCh(:,:)
    character(32) :: optr
    integer :: NChan
    logical :: Scalar
    integer :: jr, pr, zr
  contains
    procedure :: fin => FinNbodyPart
    procedure :: InitOnebodyPart
    procedure :: InitTwobodyPart
    procedure :: InitThreebodyPart
    procedure :: InitNonOrthIsospinThreebodyPart
    generic :: init => InitOneBodyPart, InitTwoBodyPart, InitThreeBodyPart, &
        & InitNonOrthIsospinThreeBodyPart
    procedure :: CopyNBodyPart
  end type NBodyPart

  type :: Op
    character(32) :: optr
    real(8) :: zero
    type(NBodyPart) :: one
    type(NBodyPart) :: two
    type(NBodyPart) :: thr
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
    module procedure :: CopyNBodyPart, CopyOp
  end interface assignment(=)

  interface operator(+)
    module procedure :: SumNBodyPart, SumOp
  end interface operator(+)

  interface operator(-)
    module procedure :: SubtractNBodyPart, SubtractOp
  end interface operator(-)

  interface operator(*)
    module procedure :: ScaleLeftNBodyPart, ScaleRightNBodyPart
    module procedure :: ScaleLeftOp, ScaleRightOp
  end interface operator(*)

contains

  subroutine FinOp(this)
    class(Op), intent(inout) :: this
    call this%one%fin()
    call this%two%fin()
    if(.not. this%is_three_body) return
    call this%thr%fin()
  end subroutine FinOp

  subroutine InitOpFromString(this, optr, ms)
    class(Op), intent(inout) :: this
    type(MSpace), intent(in) :: ms
    character(*), intent(in) :: optr
    integer :: jr, pr, zr

    select case(optr)
    case('hamil', 'Hamil')
      jr = 0
      pr = 1
      zr = 0
    end select

    call this%init(jr,pr,zr,optr,ms)

  end subroutine InitOpFromString

  subroutine InitOp(this, jr, pr, zr, optr, ms)
    use Profiler, only: timer
    class(Op), intent(inout) :: this
    integer, intent(in) :: jr, pr, zr
    character(*), intent(in) :: optr
    type(MSpace), intent(in) :: ms
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
    if(.not. ms%is_three_body) then
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

  subroutine FinNBodyPart(this)
    class(NBodyPart), intent(inout) :: this
    integer :: chbra, chket

    do chbra = 1, this%NChan
      do chket = 1, this%NChan
        if(.not. this%MatCh(chbra,chket)%is) cycle
        call this%MatCh(chbra,chket)%fin()
      end do
    end do
  end subroutine FinNBodyPart

  subroutine InitOneBodyPart(this, one, Scalar, optr, jr, pr, zr)
    use CommonLibrary, only: triag
    class(NBodyPart), intent(inout) :: this
    type(OneBodySpace), intent(in) :: one
    logical, intent(in) :: Scalar
    character(*), intent(in) :: optr
    integer, intent(in) :: jr, pr, zr
    integer :: chbra, chket
    integer :: jbra, pbra, zbra, nbra
    integer :: jket, pket, zket, nket

    this%optr = optr
    this%Scalar = Scalar
    this%jr = jr
    this%pr = pr
    this%zr = zr
    this%NChan = one%NChan

    allocate(this%MatCh(this%NChan,this%NChan))

    do chbra = 1, one%NChan
      jbra = one%jpz(chbra)%j
      pbra = one%jpz(chbra)%p
      zbra = one%jpz(chbra)%z
      nbra = one%jpz(chbra)%nst
      do chket = 1, one%NChan
        jket = one%jpz(chket)%j
        pket = one%jpz(chket)%p
        zket = one%jpz(chket)%z
        nket = one%jpz(chket)%nst

        if(triag(jbra, jket, 2*jr)) cycle
        if(pbra * pket * pr /= 1) cycle
        if(zbra - 2*zr - zket /= 0) cycle
        this%MatCh(chbra,chket)%is = .true.
        this%MatCh(chbra,chket)%ndims = [nbra,nket]
        call this%MatCh(chbra,chket)%zeros(nbra,nket)
      end do
    end do

  end subroutine InitOneBodyPart

  subroutine InitTwoBodyPart(this, two, Scalar, optr, jr, pr, zr)
    use CommonLibrary, only: triag
    class(NBodyPart), intent(inout) :: this
    type(TwoBodySpace), intent(in) :: two
    logical, intent(in) :: Scalar
    character(*), intent(in) :: optr
    integer, intent(in) :: jr, pr, zr
    integer :: chbra, chket
    integer :: jbra, pbra, zbra, nbra
    integer :: jket, pket, zket, nket

    this%optr = optr
    this%Scalar = Scalar
    this%jr = jr
    this%pr = pr
    this%zr = zr
    this%NChan = two%NChan

    allocate(this%MatCh(this%NChan,this%NChan))

    do chbra = 1, this%NChan
      jbra = two%jpz(chbra)%j
      pbra = two%jpz(chbra)%p
      zbra = two%jpz(chbra)%z
      nbra = two%jpz(chbra)%nst
      do chket = 1, this%NChan
        jket = two%jpz(chket)%j
        pket = two%jpz(chket)%p
        zket = two%jpz(chket)%z
        nket = two%jpz(chket)%nst

        if(triag(jbra, jket, jr)) cycle
        if(pbra * pket * pr /= 1) cycle
        if(zbra - zr - zket /= 0) cycle
        this%MatCh(chbra,chket)%is = .true.
        this%MatCh(chbra,chket)%ndims = [nbra,nket]
        call this%MatCh(chbra,chket)%zeros(nbra,nket)
      end do
    end do

  end subroutine InitTwoBodyPart

  subroutine InitThreeBodyPart(this, thr, Scalar, optr, jr, pr, zr)
    use CommonLibrary, only: triag
    class(NBodyPart), intent(inout) :: this
    type(ThreeBodySpace), intent(in) :: thr
    logical, intent(in) :: Scalar
    character(*), intent(in) :: optr
    integer, intent(in) :: jr, pr, zr
    integer :: chbra, chket
    integer :: jbra, pbra, zbra, nbra
    integer :: jket, pket, zket, nket

    this%optr = optr
    this%Scalar = Scalar
    this%jr = jr
    this%pr = pr
    this%zr = zr
    this%NChan = thr%NChan

    allocate(this%MatCh(this%NChan,this%NChan))

    do chbra = 1, this%NChan
      jbra = thr%jpz(chbra)%j
      pbra = thr%jpz(chbra)%p
      zbra = thr%jpz(chbra)%z
      nbra = thr%jpz(chbra)%nst
      do chket = 1, this%NChan
        jket = thr%jpz(chket)%j
        pket = thr%jpz(chket)%p
        zket = thr%jpz(chket)%z
        nket = thr%jpz(chket)%nst

        if(triag(jbra, jket, 2*jr)) cycle
        if(pbra * pket * pr /= 1) cycle
        if(zbra - 2*zr - zket /= 0) cycle
        this%MatCh(chbra,chket)%is = .true.
        this%MatCh(chbra,chket)%ndims = [nbra,nket]
        call this%MatCh(chbra,chket)%zeros(nbra,nket)
      end do
    end do

  end subroutine InitThreeBodyPart

  subroutine InitNonOrthIsospinThreeBodyPart(this, thr, Scalar, optr, jr, pr, zr)
    use CommonLibrary, only: triag
    class(NBodyPart), intent(inout) :: this
    type(NonOrthIsospinThreeBodySpace), intent(in) :: thr
    logical, intent(in) :: Scalar
    character(*), intent(in) :: optr
    integer, intent(in) :: jr, pr, zr
    integer :: chbra, chket
    integer :: jbra, pbra, zbra, nbra
    integer :: jket, pket, zket, nket

    this%optr = optr
    this%Scalar = Scalar
    this%jr = jr
    this%pr = pr
    this%zr = zr
    this%NChan = thr%NChan

    allocate(this%MatCh(this%NChan,this%NChan))

    do chbra = 1, this%NChan
      jbra = thr%jpt(chbra)%j
      pbra = thr%jpt(chbra)%p
      zbra = thr%jpt(chbra)%t
      nbra = thr%jpt(chbra)%nst
      do chket = 1, this%NChan
        jket = thr%jpt(chket)%j
        pket = thr%jpt(chket)%p
        zket = thr%jpt(chket)%t
        nket = thr%jpt(chket)%nst

        if(triag(jbra, jket, 2*jr)) cycle
        if(pbra * pket * pr /= 1) cycle
        if(zbra - 2*zr - zket /= 0) cycle
        this%MatCh(chbra,chket)%is = .true.
        this%MatCh(chbra,chket)%ndims = [nbra,nket]
        call this%MatCh(chbra,chket)%zeros(nbra,nket)
      end do
    end do

  end subroutine InitNonOrthIsospinThreeBodyPart

  subroutine CopyNBodyPart(a, b)
    class(NBodyPart), intent(inout) :: a
    type(NBodyPart), intent(in) :: b
    integer :: chbra, chket

    if(allocated(a%MatCh)) call a%fin()
    a%optr = b%optr
    a%NChan = b%NChan
    a%Scalar = b%Scalar
    a%jr = b%jr
    a%pr = b%pr
    a%zr = b%zr

    allocate(a%MatCh(b%NChan,b%NChan))
    do chbra = 1, b%NChan
      do chket = 1, b%NChan
        if(.not. b%MatCh(chbra,chket)%is) cycle
        a%MatCh(chbra,chket)%is = b%MatCh(chbra,chket)%is
        a%MatCh(chbra,chket)%ndims(:) = b%MatCh(chbra,chket)%ndims(:)
        a%MatCh(chbra,chket)%DMat = b%MatCh(chbra,chket)%DMat
      end do
    end do

  end subroutine CopyNBodyPart

  function SumNBodyPart(a, b) result(c)
    type(NBodyPart), intent(in) :: a, b
    type(NBodyPart) :: c
    integer :: chbra, chket
    c = a
    do chbra = 1, b%NChan
      do chket = 1, b%NChan
        if(.not. b%MatCh(chbra,chket)%is) cycle
        c%MatCh(chbra,chket)%is = b%MatCh(chbra,chket)%is
        c%MatCh(chbra,chket)%ndims(:) = b%MatCh(chbra,chket)%ndims(:)
        c%MatCh(chbra,chket)%DMat = &
            & a%MatCh(chbra,chket)%DMat + b%MatCh(chbra,chket)%DMat
      end do
    end do
  end function SumNBodyPart

  function SubtractNBodyPart(a, b) result(c)
    type(NBodyPart), intent(in) :: a, b
    type(NBodyPart) :: c
    integer :: chbra, chket
    c = a
    do chbra = 1, b%NChan
      do chket = 1, b%NChan
        if(.not. b%MatCh(chbra,chket)%is) cycle
        c%MatCh(chbra,chket)%is = b%MatCh(chbra,chket)%is
        c%MatCh(chbra,chket)%ndims(:) = b%MatCh(chbra,chket)%ndims(:)
        c%MatCh(chbra,chket)%DMat = &
            & a%MatCh(chbra,chket)%DMat - b%MatCh(chbra,chket)%DMat
      end do
    end do
  end function SubtractNBodyPart

  function ScaleLeftNBodyPart(a, b) result(c)
    type(NBodyPart), intent(in) :: a
    real(8), intent(in) :: b
    type(NBodyPart) :: c
    integer :: chbra, chket
    c = a
    do chbra = 1, a%NChan
      do chket = 1, a%NChan
        if(.not. a%MatCh(chbra,chket)%is) cycle
        c%MatCh(chbra,chket)%is = a%MatCh(chbra,chket)%is
        c%MatCh(chbra,chket)%ndims(:) = a%MatCh(chbra,chket)%ndims(:)
        c%MatCh(chbra,chket)%DMat = &
            & a%MatCh(chbra,chket)%DMat * b
      end do
    end do
  end function ScaleLeftNBodyPart

  function ScaleRightNBodyPart(a, b) result(c)
    type(NBodyPart), intent(in) :: b
    real(8), intent(in) :: a
    type(NBodyPart) :: c
    integer :: chbra, chket
    c = b
    do chbra = 1, b%NChan
      do chket = 1, b%NChan
        if(.not. b%MatCh(chbra,chket)%is) cycle
        c%MatCh(chbra,chket)%is = b%MatCh(chbra,chket)%is
        c%MatCh(chbra,chket)%ndims(:) = b%MatCh(chbra,chket)%ndims(:)
        c%MatCh(chbra,chket)%DMat = &
            & b%MatCh(chbra,chket)%DMat * a
      end do
    end do
  end function ScaleRightNBodyPart

end module Operators
