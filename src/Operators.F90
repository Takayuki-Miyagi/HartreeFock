module Operators
  use omp_lib
  use LinAlgLib
  use ModelSpace
  implicit none

  type, extends(DMat) :: NBodyChannel
    integer :: ndims(2)
    logical :: is = .false.
  contains
    procedure :: SetOneBodyChannel
    procedure :: SetTwoBodyChannel
    procedure :: SetThreeBodyChannel
    generic :: set => SetOneBodyChannel, setTwoBodyChannel, SetThreeBodyChannel
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

    procedure :: GetTwBME_general
    procedure :: GetTwBME_scalar
    generic :: GetTwBME => GetTwBME_general, GetTwBME_scalar
    procedure :: SetTwBME_general
    procedure :: SetTwBME_scalar
    generic :: SetTwBME => SetTwBME_general, SetTwBME_scalar
    procedure :: AddToTwBME_general
    procedure :: AddToTwBME_scalar
    generic :: AddToTwBME => AddToTwBME_general, AddToTwBME_scalar

    procedure :: GetThBMEpn_scalar
    procedure :: GetThBMEIso_scalar
    generic :: GetThBME => GetThBMEpn_scalar, GetThBMEIso_scalar

    procedure :: ReadTwoBodyFile
    procedure :: ReadIsospinThreeBodyFile
    generic :: ReadFile => ReadTwoBodyFile, ReadIsospinThreeBodyFile

    procedure :: SetOneBodyPart
    procedure :: SetTwoBodyPart
    procedure :: SetThreeBodyPart
    generic :: set => SetOneBodyPart, SetTwoBodyPart, SetThreeBodyPart
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

    procedure :: SetHamiltonian
    procedure :: SetOperators
    !procedure :: SetOperatorsFromFiles
  end type Op

  type, private :: ReadFiles
    character(:), allocatable :: file_nn, file_3n
    integer :: emax2, e2max2, lmax2
    integer :: emax3, e2max3, e3max3, lmax3
  contains
    procedure :: ReadScalar2BFile
    procedure :: ReadTensor2BFile
    procedure :: read_scalar_me2j_txt
    procedure :: read_scalar_me2j_bin
    procedure :: read_scalar_myg_txt
    procedure :: read_scalar_myg_bin
    procedure :: read_scalar_snt_txt
    procedure :: read_scalar_nv_txt
  end type ReadFiles

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
    use DefineOperators, only: GetOperatorRank
    class(Op), intent(inout) :: this
    type(MSpace), intent(in) :: ms
    character(*), intent(in) :: optr
    integer :: jr, pr, zr

    call GetOperatorRank(optr,jr,pr,zr)
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

  subroutine SetHamiltonian(this, ms, file_nn, file_3n, &
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
    rd%lmax2 = ms%emax
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
    rd%lmax3 = ms%emax
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
  end subroutine SetHamiltonian

  subroutine SetOperators(this, ms)
    class(Op), intent(inout) :: this
    type(MSpace), intent(in) :: ms

    call this%one%set(ms%one, ms%sps, ms%hw, ms%A, ms%Z, ms%N)
    call this%two%set(ms%two, ms%sps, ms%hw, ms%A, ms%Z, ms%N)
    return
  end subroutine SetOperators

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

  subroutine SetOneBodyPart(this, ms, sps, hw, A, Z, N)
    class(NBodyPart), intent(inout) :: this
    type(OneBodySpace), intent(in) :: ms
    type(Orbits), intent(in) :: sps
    real(8), intent(in) :: hw
    integer, intent(in) :: A, Z, N
    integer :: chbra
    integer :: chket

    do chbra = 1, this%NChan
      do chket = 1, this%NChan
        if(.not. this%MatCh(chbra,chket)%is) cycle
        call this%MatCh(chbra,chket)%set(this%optr, sps, &
            & ms%jpz(chbra), ms%jpz(chket), hw, A, Z, N)
      end do
    end do
  end subroutine SetOneBodyPart

  function GetTwBME_general(this,sps,ms,i1,i2,i3,i4,J12,J34) result(r)
    use CommonLibrary, only: triag
    real(8) :: r
    class(NBodyPart), intent(in) :: this
    type(Orbits), intent(in) :: sps
    type(TwoBodySpace), intent(in) :: ms
    integer, intent(in) :: i1,i2,i3,i4,J12,J34
    integer :: P12, P34
    integer :: Z12, Z34
    integer :: chbra, chket, bra, ket, iphase

    r = 0.d0
    P12 = (-1) ** (sps%orb(i1)%l + sps%orb(i2)%l)
    P34 = (-1) ** (sps%orb(i3)%l + sps%orb(i4)%l)
    Z12 = (-1) ** (sps%orb(i1)%z + sps%orb(i2)%z)
    Z34 = (-1) ** (sps%orb(i3)%z + sps%orb(i4)%z)

    if(triag(J12, J34, this%jr)) stop "Error, in GetTwBME_general: J"
    if(P12 * P34 * this%pr /= 1) stop "Error, in GetTwBME_general: P"
    if(Z12 - Z34 - this%zr /= 0) stop "Error, in GetTwBME_general: Tz"

    chbra = ms%jpz2ch(J12,P12,Z12)
    chket = ms%jpz2ch(J34,P34,Z34)
    if(chbra * chket == 0) return

    bra = ms%jpz(chbra)%spis2n(i1,i2)
    ket = ms%jpz(chket)%spis2n(i3,i4)
    if(bra * ket == 0) return

    iphase = ms%jpz(chbra)%iphase(i1,i2) * &
        &    ms%jpz(chket)%iphase(i3,i4)

    r = dble(iphase) * this%MatCh(chbra,chket)%m(bra,ket)
  end function GetTwBME_general

  function GetTwBME_scalar(this,sps,ms,i1,i2,i3,i4,J) result(r)
    use CommonLibrary, only: triag
    real(8) :: r
    class(NBodyPart), intent(in) :: this
    type(Orbits), intent(in) :: sps
    type(TwoBodySpace), intent(in) :: ms
    integer, intent(in) :: i1,i2,i3,i4,J
    integer :: P12, P34
    integer :: Z12, Z34
    integer :: ch, bra, ket, iphase

    r = 0.d0
    P12 = (-1) ** (sps%orb(i1)%l + sps%orb(i2)%l)
    P34 = (-1) ** (sps%orb(i3)%l + sps%orb(i4)%l)
    Z12 = (-1) ** (sps%orb(i1)%z + sps%orb(i2)%z)
    Z34 = (-1) ** (sps%orb(i3)%z + sps%orb(i4)%z)

    if(P12 * P34 /= 1) stop "Error, in GetTwBME_general: P"
    if(Z12 - Z34 /= 0) stop "Error, in GetTwBME_general: Tz"

    ch = ms%jpz2ch(J,P12,Z12)
    if(ch == 0) return

    bra = ms%jpz(ch)%spis2n(i1,i2)
    ket = ms%jpz(ch)%spis2n(i3,i4)
    if(bra * ket == 0) return

    iphase = ms%jpz(ch)%iphase(i1,i2) * &
        &    ms%jpz(ch)%iphase(i3,i4)

    r = dble(iphase) * this%MatCh(ch,ch)%m(bra,ket)
  end function GetTwBME_scalar

  subroutine SetTwBME_general(this,sps,ms,i1,i2,i3,i4,J12,J34,me)
    use CommonLibrary, only: triag
    class(NBodyPart), intent(inout) :: this
    type(Orbits), intent(in) :: sps
    type(TwoBodySpace), intent(in) :: ms
    integer, intent(in) :: i1,i2,i3,i4,J12,J34
    real(8), intent(in) :: me
    integer :: P12, P34
    integer :: Z12, Z34
    integer :: chbra, chket, bra, ket, iphase

    P12 = (-1) ** (sps%orb(i1)%l + sps%orb(i2)%l)
    P34 = (-1) ** (sps%orb(i3)%l + sps%orb(i4)%l)
    Z12 = (-1) ** (sps%orb(i1)%z + sps%orb(i2)%z)
    Z34 = (-1) ** (sps%orb(i3)%z + sps%orb(i4)%z)

    if(triag(J12, J34, this%jr)) stop "Error, in GetTwBME_general: J"
    if(P12 * P34 * this%pr /= 1) stop "Error, in GetTwBME_general: P"
    if(Z12 - Z34 - this%zr /= 0) stop "Error, in GetTwBME_general: Tz"

    chbra = ms%jpz2ch(J12,P12,Z12)
    chket = ms%jpz2ch(J34,P34,Z34)
    if(chbra * chket == 0) return

    bra = ms%jpz(chbra)%spis2n(i1,i2)
    ket = ms%jpz(chket)%spis2n(i3,i4)
    if(bra * ket == 0) return

    iphase = ms%jpz(chbra)%iphase(i1,i2) * &
        &    ms%jpz(chket)%iphase(i3,i4)

    this%MatCh(chbra,chket)%m(bra,ket) = dble(iphase) * me
  end subroutine SetTwBME_general

  subroutine SetTwBME_scalar(this,sps,ms,i1,i2,i3,i4,J,me)
    use CommonLibrary, only: triag
    real(8) :: r
    class(NBodyPart), intent(inout) :: this
    type(Orbits), intent(in) :: sps
    type(TwoBodySpace), intent(in) :: ms
    integer, intent(in) :: i1,i2,i3,i4,J
    real(8), intent(in) :: me
    integer :: P12, P34
    integer :: Z12, Z34
    integer :: ch, bra, ket, iphase

    r = 0.d0
    P12 = (-1) ** (sps%orb(i1)%l + sps%orb(i2)%l)
    P34 = (-1) ** (sps%orb(i3)%l + sps%orb(i4)%l)
    Z12 = (-1) ** (sps%orb(i1)%z + sps%orb(i2)%z)
    Z34 = (-1) ** (sps%orb(i3)%z + sps%orb(i4)%z)

    if(P12 * P34 /= 1) stop "Error, in SetTwBME_general: P"
    if(Z12 - Z34 /= 0) stop "Error, in SetTwBME_general: Tz"

    ch = ms%jpz2ch(J,P12,Z12)
    if(ch == 0) return

    bra = ms%jpz(ch)%spis2n(i1,i2)
    ket = ms%jpz(ch)%spis2n(i3,i4)
    if(bra * ket == 0) return

    iphase = ms%jpz(ch)%iphase(i1,i2) * &
        &    ms%jpz(ch)%iphase(i3,i4)

    this%MatCh(ch,ch)%m(bra,ket) = r * dble(iphase)
  end subroutine SetTwBME_scalar

  subroutine AddToTwBME_general(this,sps,ms,i1,i2,i3,i4,J12,J34,me)
    use CommonLibrary, only: triag
    class(NBodyPart), intent(inout) :: this
    type(Orbits), intent(in) :: sps
    type(TwoBodySpace), intent(in) :: ms
    integer, intent(in) :: i1,i2,i3,i4,J12,J34
    real(8), intent(in) :: me
    integer :: P12, P34
    integer :: Z12, Z34
    integer :: chbra, chket, bra, ket, iphase

    P12 = (-1) ** (sps%orb(i1)%l + sps%orb(i2)%l)
    P34 = (-1) ** (sps%orb(i3)%l + sps%orb(i4)%l)
    Z12 = (-1) ** (sps%orb(i1)%z + sps%orb(i2)%z)
    Z34 = (-1) ** (sps%orb(i3)%z + sps%orb(i4)%z)

    if(triag(J12, J34, this%jr)) stop "Error, in GetTwBME_general: J"
    if(P12 * P34 * this%pr /= 1) stop "Error, in GetTwBME_general: P"
    if(Z12 - Z34 - this%zr /= 0) stop "Error, in GetTwBME_general: Tz"

    chbra = ms%jpz2ch(J12,P12,Z12)
    chket = ms%jpz2ch(J34,P34,Z34)
    if(chbra * chket == 0) return

    bra = ms%jpz(chbra)%spis2n(i1,i2)
    ket = ms%jpz(chket)%spis2n(i3,i4)
    if(bra * ket == 0) return

    iphase = ms%jpz(chbra)%iphase(i1,i2) * &
        &    ms%jpz(chket)%iphase(i3,i4)

    this%MatCh(chbra,chket)%m(bra,ket) = &
        & this%MatCh(chbra,chket)%m(bra,ket) + dble(iphase) * me
  end subroutine AddToTwBME_general

  subroutine AddToTwBME_scalar(this,sps,ms,i1,i2,i3,i4,J,me)
    use CommonLibrary, only: triag
    real(8) :: r
    class(NBodyPart), intent(inout) :: this
    type(Orbits), intent(in) :: sps
    type(TwoBodySpace), intent(in) :: ms
    integer, intent(in) :: i1,i2,i3,i4,J
    real(8), intent(in) :: me
    integer :: P12, P34
    integer :: Z12, Z34
    integer :: ch, bra, ket, iphase

    r = 0.d0
    P12 = (-1) ** (sps%orb(i1)%l + sps%orb(i2)%l)
    P34 = (-1) ** (sps%orb(i3)%l + sps%orb(i4)%l)
    Z12 = (-1) ** (sps%orb(i1)%z + sps%orb(i2)%z)
    Z34 = (-1) ** (sps%orb(i3)%z + sps%orb(i4)%z)

    if(P12 * P34 /= 1) stop "Error, in AddToTwBME_general: P"
    if(Z12 - Z34 /= 0) stop "Error, in AddToTwBME_general: Tz"

    ch = ms%jpz2ch(J,P12,Z12)
    if(ch == 0) return

    bra = ms%jpz(ch)%spis2n(i1,i2)
    ket = ms%jpz(ch)%spis2n(i3,i4)
    if(bra * ket == 0) return

    iphase = ms%jpz(ch)%iphase(i1,i2) * &
        &    ms%jpz(ch)%iphase(i3,i4)

    this%MatCh(ch,ch)%m(bra,ket) = this%MatCh(ch,ch)%m(bra,ket) + &
        & r * dble(iphase)
  end subroutine AddToTwBME_scalar

  subroutine SetTwoBodyPart(this, ms, sps, hw, A, Z, N)
    class(NBodyPart), intent(inout) :: this
    type(TwoBodySpace), intent(in) :: ms
    type(Orbits), intent(in) :: sps
    real(8), intent(in) :: hw
    integer, intent(in) :: A, Z, N
    integer :: chbra
    integer :: chket

    do chbra = 1, this%NChan
      do chket = 1, this%NChan
        if(.not. this%MatCh(chbra,chket)%is) cycle
        call this%MatCh(chbra,chket)%set(this%optr, sps, &
            & ms%jpz(chbra), ms%jpz(chket), hw, A, Z, N)
      end do
    end do
  end subroutine SetTwoBodyPart

  function GetThBMEIso_scalar(this,sps,ms,i1,i2,i3,J12,T12,&
        & i4,i5,i6,J45,T45,J,T) result(r)
    class(NBodyPart), intent(in) :: this
    type(OrbitsIsospin), intent(in) :: sps
    type(NonOrthIsospinThreeBodySpace), intent(in) :: ms
    integer, intent(in) :: i1,i2,i3,i4,i5,i6
    integer, intent(in) :: J12,T12,J45,T45,J,T
    integer :: ch, idxbra, idxket, bra, ket
    integer :: P123, P456
    integer :: num_bra, num_ket
    integer :: isorted_bra, isorted_ket
    integer :: ibra, iket
    real(8) :: r

    r = 0.d0
    P123 = (-1) ** (sps%orb(i1)%l+sps%orb(i2)%l+sps%orb(i3)%l)
    P456 = (-1) ** (sps%orb(i4)%l+sps%orb(i5)%l+sps%orb(i6)%l)
    if(P123 * P456 /= 1) stop 'Error in GetThBMEIso_scalar: P'
    ch = ms%jpt2ch(J,P123,T)
    if(ch == 0) return
    idxbra = ms%jpt(ch)%sorting(i1,i2,i3)
    idxket = ms%jpt(ch)%sorting(i4,i5,i6)
    if(idxbra * idxket == 0) return
    if(i1 == i2 .and. mod(J12+T12,2) == 0) return
    if(i4 == i5 .and. mod(J45+T45,2) == 0) return
    isorted_bra = ms%jpt(ch)%sort(idxbra)%idx_sorted
    isorted_ket = ms%jpt(ch)%sort(idxket)%idx_sorted
    if(isorted_bra * isorted_ket == 0) return

    do ibra = 1, ms%jpt(ch)%sort(idxbra)%JT(J12,T12)%n
      bra = ms%jpt(ch)%sort(idxbra)%JT(J12,T12)%idx2num(ibra)
      do iket = 1, ms%jpt(ch)%sort(idxket)%JT(J45,T45)%n
        ket = ms%jpt(ch)%sort(idxket)%JT(J45,T45)%idx2num(iket)
        r = r + dble(this%MatCh(ch,ch)%m(bra,ket) * &
            & ms%jpt(ch)%sort(idxbra)%JT(J12,T12)%TrnsCoef(ibra) * &
            & ms%jpt(ch)%sort(idxket)%JT(J45,T45)%TrnsCoef(iket))
      end do
    end do
  end function GetThBMEIso_scalar

  function GetThBMEpn_scalar(this,ms,i1,i2,i3,J12,&
        & i4,i5,i6,J45,J) result(r)
    use CommonLibrary, only: dcg
    class(NBodyPart), intent(in) :: this
    type(MSpace), intent(in) :: ms
    integer, intent(in) :: i1,i2,i3,i4,i5,i6
    integer, intent(in) :: J12,J45,J
    integer :: z1, z2, z3, z4, z5, z6, T12, T45, T
    integer :: a, b, c, d, e, f
    integer :: P, Z, ch
    real(8) :: r
    r = 0.d0
    if(ms%sps%orb(i1)%e + ms%sps%orb(i2)%e + ms%sps%orb(i3)%e > ms%e3max) return
    if(ms%sps%orb(i4)%e + ms%sps%orb(i5)%e + ms%sps%orb(i6)%e > ms%e3max) return
    if(i1 == i2 .and. mod(J12, 2) == 1) return
    if(i3 == i4 .and. mod(J45, 2) == 1) return
    z1 = ms%sps%orb(i1)%z
    z2 = ms%sps%orb(i2)%z
    z3 = ms%sps%orb(i3)%z
    z4 = ms%sps%orb(i4)%z
    z5 = ms%sps%orb(i5)%z
    z6 = ms%sps%orb(i6)%z

    a = ms%isps%nlj2idx( ms%sps%orb(i1)%n, ms%sps%orb(i1)%l, ms%sps%orb(i1)%j )
    b = ms%isps%nlj2idx( ms%sps%orb(i2)%n, ms%sps%orb(i2)%l, ms%sps%orb(i2)%j )
    c = ms%isps%nlj2idx( ms%sps%orb(i3)%n, ms%sps%orb(i3)%l, ms%sps%orb(i3)%j )
    d = ms%isps%nlj2idx( ms%sps%orb(i4)%n, ms%sps%orb(i4)%l, ms%sps%orb(i4)%j )
    e = ms%isps%nlj2idx( ms%sps%orb(i5)%n, ms%sps%orb(i5)%l, ms%sps%orb(i5)%j )
    f = ms%isps%nlj2idx( ms%sps%orb(i6)%n, ms%sps%orb(i6)%l, ms%sps%orb(i6)%j )

    P = (-1) ** (ms%sps%orb(i1)%l+ms%sps%orb(i2)%l+ms%sps%orb(i3)%l)
    Z = z1 + z2 + z3
    do T12 = 0, 1
      if(abs(z1+z2) > T12) cycle
      do T45 = 0, 1
        if(abs(z4+z5) > T45) cycle
        do T = max(abs(2*T12-1),abs(2*T45-1)), min(2*T12+1,2*T45+1), 2
          if(abs(Z) > T) cycle
          ch = ms%thr%jpt2ch(J,P,T)
          if(ch == 0) cycle
          r = r + &
              & this%GetThBME(ms%isps,ms%thr,a,b,c,J12,T12,&
              & d,e,f,J45,T45,J,T) * &
              & dcg(1,z1,2,z2,2*T12,z1+z2) * dcg(2*T12,z1+z2,1,z3,T,Z) * &
              & dcg(1,z4,2,z5,2*T45,z4+z5) * dcg(2*T45,z4+z5,1,z6,T,Z)
        end do
      end do
    end do

  end function GetThBMEpn_scalar

  subroutine SetThreeBodyPart(this)
    class(NBodyPart), intent(inout) :: this

    write(*,*) "No such a simple operator in three-body system"
    return
  end subroutine SetThreeBodyPart

  subroutine SetOneBodyChannel(this, optr, sps, chbra, chket, hw, A, Z, N)
    class(NBodyChannel), intent(inout) :: this
    character(*), intent(in) :: optr
    type(Orbits), intent(in) :: sps
    type(OneBodyChannel), intent(in) :: chbra, chket
    real(8), intent(in) :: hw
    integer, intent(in) :: A, Z, N
    integer :: bra, ket, ia, ib
    do bra = 1, chbra%nst
      ia = chbra%n2spi(bra)
      do ket = 1, chket%nst
        ib = chbra%n2spi(ket)

        this%m(bra,ket) = mat_elm()

      end do
    end do
  contains
    function mat_elm() result(r)
      use DefineOperators, only: one_body_element
      real(8) :: r
      integer :: na, la, ja, za
      integer :: nb, lb, jb, zb

      r = 0.d0
      na = sps%orb(ia)%n
      la = sps%orb(ia)%l
      ja = sps%orb(ia)%j
      za = sps%orb(ia)%z

      nb = sps%orb(ib)%n
      lb = sps%orb(ib)%l
      jb = sps%orb(ib)%j
      zb = sps%orb(ib)%z

      r = one_body_element(optr, [na,la,ja,za], [nb,lb,jb,zb], &
          & hw, A, Z, N)
    end function mat_elm
  end subroutine SetOneBodyChannel

  subroutine SetTwoBodyChannel(this, optr, sps, chbra, chket, hw, A, Z, N)
    class(NBodyChannel), intent(inout) :: this
    character(*), intent(in) :: optr
    type(Orbits), intent(in) :: sps
    type(TwoBodyChannel), intent(in) :: chbra, chket
    real(8), intent(in) :: hw
    integer, intent(in) :: A, Z, N
    integer :: bra, ket, ia, ib, ic, id

    do bra = 1, chbra%nst
      ia = chbra%n2spi1(bra)
      ib = chbra%n2spi2(bra)
      do ket = 1, chket%nst
        ic = chbra%n2spi1(ket)
        id = chbra%n2spi2(ket)

        this%m(bra,ket) = mat_elm()

      end do
    end do
  contains
    function mat_elm() result(r)
      use DefineOperators, only: two_body_element
      real(8) :: r
      integer :: na, la, ja, za
      integer :: nb, lb, jb, zb
      integer :: nc, lc, jc, zc
      integer :: nd, ld, jd, zd

      r = 0.d0
      na = sps%orb(ia)%n
      la = sps%orb(ia)%l
      ja = sps%orb(ia)%j
      za = sps%orb(ia)%z

      nb = sps%orb(ib)%n
      lb = sps%orb(ib)%l
      jb = sps%orb(ib)%j
      zb = sps%orb(ib)%z

      nc = sps%orb(ic)%n
      lc = sps%orb(ic)%l
      jc = sps%orb(ic)%j
      zc = sps%orb(ic)%z

      nd = sps%orb(id)%n
      ld = sps%orb(id)%l
      jd = sps%orb(id)%j
      zd = sps%orb(id)%z

      r = two_body_element(optr, [na,la,ja,za], [nb,lb,jb,zb], &
          & [nc,lc,jc,zc], [nd,ld,jd,zd], chbra%J, chket%J, hw, A, Z, N)

    end function mat_elm
  end subroutine SetTwoBodyChannel

  subroutine SetThreeBodyChannel(this)
    class(NBodyChannel), intent(inout) :: this

    write(*,*) "No such a simple operator"
    return
  end subroutine SetThreeBodyChannel

  !
  !
  !     File reading methods
  !
  !

  subroutine ReadTwoBodyFile(this, sps, ms, rd)
    class(NBodyPart), intent(inout) :: this
    type(Orbits), intent(in) :: sps
    type(TwoBodySpace), intent(in) :: ms
    type(ReadFiles), intent(in) :: rd

    select case(rd%file_nn)
    case('None', 'NONE', 'none')
      write(*,*) "Error in ReadTwoBodyFile"

    case default

      if(this%Scalar) call rd%ReadScalar2BFile(this,sps,ms)
      if(.not. this%Scalar) call rd%ReadTensor2BFile(this,sps,ms)
    end select

  end subroutine ReadTwoBodyFile

  subroutine ReadScalar2BFile(this,two,sps,ms)
    use ClassSys, only: sys
    class(ReadFiles), intent(in) :: this
    type(NBodyPart), intent(inout) :: two
    type(Orbits), intent(in) :: sps
    type(TwoBodySpace), intent(in) :: ms
    type(sys) :: s

    if(s%find(this%file_nn, '.txt') .and. s%find(this%file_nn, '.me2j')) then
      call this%read_scalar_me2j_txt(two,sps,ms)
      return
    end if

    if(s%find(this%file_nn, '.bin') .and. s%find(this%file_nn, '.me2j')) then
      call this%read_scalar_me2j_bin(two,sps,ms)
      return
    end if

    if(s%find(this%file_nn, '.txt') .and. s%find(this%file_nn, '.myg')) then
      call this%read_scalar_myg_txt(two,sps,ms)
      return
    end if

    if(s%find(this%file_nn, '.bin') .and. s%find(this%file_nn, '.myg')) then
      call this%read_scalar_myg_bin(two,sps,ms)
      return
    end if

    if(s%find(this%file_nn, '.txt') .and. s%find(this%file_nn, '.snt')) then
      call this%read_scalar_snt_txt(two,sps,ms)
      return
    end if

    if(s%find(this%file_nn, '.txt') .and. s%find(this%file_nn, '.nv')) then
      call this%read_scalar_nv_txt(two,sps,ms)
      return
    end if
  end subroutine ReadScalar2BFile

  subroutine read_scalar_me2j_txt(this,two,sps,ms)
    class(ReadFiles), intent(in) :: this
    type(NBodyPart), intent(inout) :: two
    type(Orbits), intent(in) :: sps
    type(TwoBodySpace), intent(in) :: ms
    type(OrbitsIsospin) :: sps_me2j
    integer :: nelm
    integer :: io, runit = 22
    real(8), allocatable :: v(:)

    call sps_me2j%init(this%emax2, this%lmax2)
    nelm = count_scalar_me2j(sps_me2j,this%e2max2)
    allocate(v(nelm))
    open(runit, file=this%file_nn, action='read',iostat=io)
    if(io /= 0) then
      write(*,'(2a)') 'File open error: ', trim(this%file_nn)
      return
    end if
    call get_vector_me2j_formatted(runit,sps_me2j,this%e2max2,v)
    close(runit)

    deallocate(v)
    call sps_me2j%fin()
    return
  end subroutine read_scalar_me2j_txt

  subroutine read_scalar_me2j_bin(this,two,sps,ms)
    class(ReadFiles), intent(in) :: this
    type(NBodyPart), intent(inout) :: two
    type(Orbits), intent(in) :: sps
    type(TwoBodySpace), intent(in) :: ms
    type(OrbitsIsospin) :: sps_me2j
    integer :: nelm, io, runit = 22
#ifdef single_precisionA
    real(4), allocatable :: v(:)
#else
    real(8), allocatable :: v(:)
#endif

    call sps_me2j%init(this%emax2, this%lmax2)
    nelm = count_scalar_me2j(sps_me2j,this%e2max2)
    allocate(v(nelm))
    open(runit, file=this%file_nn, action='read',iostat=io)
    if(io /= 0) then
      write(*,'(2a)') 'File open error: ', trim(this%file_nn)
      return
    end if
    read(runit) v
    close(runit)

    deallocate(v)
    call sps_me2j%fin()

    return
  end subroutine read_scalar_me2j_bin

  function count_scalar_me2j(sps,e2max) result(r)
    type(OrbitsIsospin), intent(in) :: sps
    integer, intent(in) :: e2max
    integer :: r
    integer :: a, b, c, d, dmax
    integer :: la, ja, ea
    integer :: lb, jb, eb
    integer :: lc, jc, ec
    integer :: ld, jd, ed
    integer :: J, Jmin, Jmax
    r = 0
    do a = 1, sps%norbs
      la = sps%orb(a)%l
      ja = sps%orb(a)%j
      ea = sps%orb(a)%e
      do b = 1, a
        lb = sps%orb(b)%l
        jb = sps%orb(b)%j
        eb = sps%orb(b)%e
        if(ea + eb > e2max) cycle
        do c = 1, a
          lc = sps%orb(c)%l
          jc = sps%orb(c)%j
          ec = sps%orb(c)%e
          dmax = c
          if(a == c) dmax = b
          do d = 1, dmax
            ld = sps%orb(d)%l
            jd = sps%orb(d)%j
            ed = sps%orb(d)%e
            if(ec + ed > e2max) cycle
            if( (-1)**(la+lb) /= (-1)**(lc+ld) ) cycle
            Jmin = max(abs(ja-jb), abs(jc-jd))/2
            Jmax = min(   (ja+jb),    (jc+jd))/2
            if(Jmin > Jmax) cycle
            do J = Jmin, Jmax
              r = r + 4
            end do
          end do
        end do
      end do
    end do
  end function count_scalar_me2j

  subroutine get_vector_me2j_formatted(ut,sps,e2max,v)
    integer, intent(in) :: ut, e2max
    type(OrbitsIsospin), intent(in) :: sps
    real(8), intent(inout) :: v(:)
    integer :: nelm, lines, i
    nelm = size(v)
    lines = nelm/10
    do i = 1, lines
      read(ut,*) v( (i-1)*10+1 : i*10)
    end do
    read(ut,*) v(lines*10:)
  end subroutine get_vector_me2j_formatted

  subroutine read_scalar_myg_txt(this,two,sps,ms)
    class(ReadFiles), intent(in) :: this
    type(NBodyPart), intent(inout) :: two
    type(Orbits), intent(in) :: sps
    type(TwoBodySpace), intent(in) :: ms

    return
  end subroutine read_scalar_myg_txt

  subroutine read_scalar_myg_bin(this,two,sps,ms)
    class(ReadFiles), intent(in) :: this
    type(NBodyPart), intent(inout) :: two
    type(Orbits), intent(in) :: sps
    type(TwoBodySpace), intent(in) :: ms

    return
  end subroutine read_scalar_myg_bin

  subroutine read_scalar_snt_txt(this,two,sps,ms)
    class(ReadFiles), intent(in) :: this
    type(NBodyPart), intent(inout) :: two
    type(Orbits), intent(in) :: sps
    type(TwoBodySpace), intent(in) :: ms

    return
  end subroutine read_scalar_snt_txt

  subroutine read_scalar_nv_txt(this,two,sps,ms)
    class(ReadFiles), intent(in) :: this
    type(NBodyPart), intent(inout) :: two
    type(Orbits), intent(in) :: sps
    type(TwoBodySpace), intent(in) :: ms

    return
  end subroutine read_scalar_nv_txt



  subroutine ReadTensor2BFile(this,two,sps,ms)
    use ClassSys, only: sys
    class(ReadFiles), intent(in) :: this
    type(NBodyPart), intent(inout) :: two
    type(Orbits), intent(in) :: sps
    type(TwoBodySpace), intent(in) :: ms
  end subroutine ReadTensor2BFile

  subroutine ReadIsospinThreeBodyFile(this, sps, ms, rd)
    class(NBodyPart), intent(inout) :: this
    type(OrbitsIsospin), intent(in) :: sps
    type(NonOrthIsospinThreeBodySpace), intent(in) :: ms
    type(ReadFiles), intent(in) :: rd
  end subroutine ReadIsospinThreeBodyFile

end module Operators

! test for operators
program test
  use Profiler, only: timer
  use CommonLibrary, only: &
      &init_dbinomial_triangle, fin_dbinomial_triangle
  use ModelSpace, only: MSpace
  use Operators

  type(MSpace) :: ms
  type(Op) :: Hin

  call timer%init()
  call init_dbinomial_triangle()

  call ms%init('O16', 20.d0, 8, 16, e3max=8)
  call Hin%init("hamil",ms)

  call Hin%fin()
  call ms%fin()

  call fin_dbinomial_triangle()
  call timer%fin()
end program test