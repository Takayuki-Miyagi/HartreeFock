module NOperators
  use omp_lib
  use LinAlgLib
  use ModelSpace
  implicit none

  public :: NBodyChannelSp
  public :: NBodyChannel
  public :: NBodyPartSp
  public :: NBodyPart

  private :: FinNBodyPart
  private :: FinNBodyPartSp
  private :: InitOneBodyPart
  private :: InitTwoBodyPart
  private :: InitThreeBodyPart
  private :: InitNonOrthIsospinThreeBodyPart
  private :: InitNonOrthIsospinThreeBodyPartSp

  private :: SetOneBOdyPart
  private :: GetTwBME_general
  private :: GetTwBME_scalar
  private :: SetTwBME_general
  private :: SetTwBME_scalar
  private :: AddToTwBME_general
  private :: AddToTwBME_scalar
  private :: SetTwoBodyPart

  private :: GetThBMEIso_scalar
  private :: GetThBMEpn_scalar
  private :: GetThBMEIso_scalar_sp
  private :: GetThBMEpn_scalar_sp

  private :: GetThBMEIso_general
  private :: GetThBMEpn_general
  private :: GetThBMEIso_general_sp
  private :: GetThBMEpn_general_sp

  private :: SetOneBodyChannel
  private :: SetTwoBodyChannel
  private :: PrintNBodyPartSp
  private :: PrintNBodyPart

  private :: ReadTwoBodyFile
  private :: ReadIsospinThreeBodyFileSp
  private :: ReadIsospinThreeBodyFile

  private :: CopyNBodyPartSp
  private :: CopyNBodyPart
  private :: SumNBodyPartSp
  private :: SumNBodyPart
  private :: SubtractNBodyPartSp
  private :: SubtractNBodyPart
  private :: ScaleNBodyPartSp
  private :: ScaleNBodyPart

  ! reading methdos
  ! two-body scalar
  private :: ReadScalar2BFile
  private :: read_scalar_me2j_txt
  private :: read_scalar_me2j_bin
  private :: count_scalar_me2j
  private :: get_vector_me2j_formatted
  private :: read_scalar_myg_txt
  private :: read_Scalar_myg_bin
  private :: read_scalar_snt_txt
  private :: read_scalar_nv_txt
  ! two-body tensor
  private :: ReadTensor2BFile
  ! three-body scalar
  private :: ReadScalar3BFile
  private :: read_scalar_3bme_txt
  private :: read_scalar_3bme_bin
  private :: ReadScalar3BFileSp
  private :: read_scalar_3bme_txt_sp
  private :: read_scalar_3bme_bin_sp
  private :: count_scalar_3bme
  private :: store_scalar_3bme
  private :: store_scalar_3bme_sp
  ! three-body tensor
  private :: ReadTensor3BFile
  private :: read_tensor_3bme_txt
  private :: read_tensor_3bme_bin
  private :: ReadTensor3BFileSp
  private :: read_tensor_3bme_txt_sp
  private :: read_tensor_3bme_bin_sp


  type, extends(SMat) :: NBodyChannelSp
    integer :: ndims(2)
    logical :: is = .false.
  end type NBodyChannelSp

  type :: NBodyPartSp
    type(NBodyChannelSp), allocatable :: MatCh(:,:)
    character(32) :: optr
    integer :: NChan
    logical :: Scalar
    integer :: jr, pr, tr, zr
  contains
    procedure :: fin => FinNbodyPartSp
    procedure :: init => InitNonOrthIsospinThreebodyPartSp

    procedure :: GetThBMEpn_scalar_Sp
    procedure :: GetThBMEIso_scalar_Sp
    procedure :: GetThBMEpn_general_Sp
    procedure :: GetThBMEIso_general_Sp
    generic :: GetThBME => GetThBMEpn_scalar_Sp, GetThBMEIso_scalar_Sp, &
        & GetThBMEIso_general_Sp, GetThBMEpn_general_Sp
    procedure :: ReadFile => ReadIsospinThreeBodyFileSp
    procedure :: prt => PrintNBodyPartSp
    procedure :: CopyNBodyPartSp
    procedure :: SumNBodyPartSp
    procedure :: SubtractNBodyPartSp
    procedure :: ScaleNBodyPartSp
    generic :: assignment(=) => CopyNBodyPartSp
    generic :: operator(+) => SumNBodyPartSp
    generic :: operator(-) => SubtractNBodyPartSp
    generic :: operator(*) => ScaleNBodyPartSp

    procedure :: NormalOrderingFromSp3To2
  end type NBodyPartSp

  type, extends(DMat) :: NBodyChannel
    integer :: ndims(2)
    logical :: is = .false.
  contains
    procedure :: SetOneBodyChannel
    procedure :: SetTwoBodyChannel
    generic :: set => SetOneBodyChannel, setTwoBodyChannel
  end type NBodyChannel

  type :: NBodyPart
    type(NBodyChannel), allocatable :: MatCh(:,:)
    character(32) :: optr = 'none'
    integer :: NChan
    logical :: Scalar
    integer :: jr, pr, tr, zr
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
    procedure :: GetThBMEpn_general
    procedure :: GetThBMEIso_general
    generic :: GetThBME => GetThBMEpn_scalar, GetThBMEIso_scalar, &
        & GetThBMEIso_general, GetThBMEpn_general

    procedure :: ReadTwoBodyFile
    procedure :: ReadIsospinThreeBodyFile
    generic :: ReadFile => ReadTwoBodyFile, ReadIsospinThreeBodyFile

    procedure :: SetOneBodyPart
    procedure :: SetTwoBodyPart
    generic :: set => SetOneBodyPart, SetTwoBodyPart
    procedure :: prt => PrintNBodyPart

    procedure :: NormalOrderingFrom3To2
    procedure :: NormalOrderingFrom2To1
    procedure :: NormalOrderingFrom1To0

    procedure :: CopyNBodyPart
    procedure :: SumNBodyPart
    procedure :: SubtractNBodyPart
    procedure :: ScaleNBodyPart
    generic :: assignment(=) => CopyNBodyPart
    generic :: operator(+) => SumNBodyPart
    generic :: operator(-) => SubtractNBodyPart
    generic :: operator(*) => ScaleNBodyPart
  end type NBodyPart

  type :: ReadFiles
    character(:), allocatable :: file_nn, file_3n
    integer :: emax2, e2max2, lmax2
    integer :: emax3, e2max3, e3max3, lmax3
  contains
    ! methods for two-body matrix element
    procedure :: ReadScalar2BFile
    procedure :: ReadTensor2BFile
    procedure :: read_scalar_me2j_txt
    procedure :: read_scalar_me2j_bin
    procedure :: read_scalar_myg_txt
    procedure :: read_scalar_myg_bin
    procedure :: read_scalar_snt_txt
    procedure :: read_scalar_nv_txt

    ! methods for three-body marix element
    procedure :: ReadScalar3BFile
    procedure :: ReadScalar3BFileSp
    procedure :: ReadTensor3BFile
    procedure :: ReadTensor3BFileSp
    procedure :: read_scalar_3bme_txt
    procedure :: read_scalar_3bme_bin
    procedure :: read_scalar_3bme_txt_sp
    procedure :: read_scalar_3bme_bin_sp
    procedure :: read_tensor_3bme_txt
    procedure :: read_tensor_3bme_bin
    procedure :: read_tensor_3bme_txt_sp
    procedure :: read_tensor_3bme_bin_sp

  end type ReadFiles

contains

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

  subroutine FinNBodyPartSp(this)
    class(NBodyPartSp), intent(inout) :: this
    integer :: chbra, chket

    do chbra = 1, this%NChan
      do chket = 1, this%NChan
        if(.not. this%MatCh(chbra,chket)%is) cycle
        call this%MatCh(chbra,chket)%fin()
      end do
    end do
  end subroutine FinNBodyPartSp

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
    this%tr = 0
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
    this%tr = 0
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
    this%tr = 0
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

  subroutine InitNonOrthIsospinThreeBodyPart(this, thr, Scalar, optr, jr, pr, tr)
    use CommonLibrary, only: triag
    class(NBodyPart), intent(inout) :: this
    type(NonOrthIsospinThreeBodySpace), intent(in) :: thr
    logical, intent(in) :: Scalar
    character(*), intent(in) :: optr
    integer, intent(in) :: jr, pr, tr
    integer :: chbra, chket
    integer :: jbra, pbra, tbra, nbra
    integer :: jket, pket, tket, nket

    this%optr = optr
    this%Scalar = Scalar
    this%jr = jr
    this%pr = pr
    this%tr = tr
    this%zr = 0
    this%NChan = thr%NChan

    allocate(this%MatCh(this%NChan,this%NChan))

    do chbra = 1, this%NChan
      jbra = thr%jpt(chbra)%j
      pbra = thr%jpt(chbra)%p
      tbra = thr%jpt(chbra)%t
      nbra = thr%jpt(chbra)%nst
      do chket = 1, this%NChan
        jket = thr%jpt(chket)%j
        pket = thr%jpt(chket)%p
        tket = thr%jpt(chket)%t
        nket = thr%jpt(chket)%nst

        if(triag(jbra, jket, 2*jr)) cycle
        if(pbra * pket * pr /= 1) cycle
        if(triag(tbra, tket, 2*tr)) cycle
        this%MatCh(chbra,chket)%is = .true.
        this%MatCh(chbra,chket)%ndims = [nbra,nket]
        call this%MatCh(chbra,chket)%zeros(nbra,nket)
      end do
    end do
  end subroutine InitNonOrthIsospinThreeBodyPart

  subroutine InitNonOrthIsospinThreeBodyPartSp(this, thr, Scalar, optr, jr, pr, tr)
    use CommonLibrary, only: triag
    class(NBodyPartSp), intent(inout) :: this
    type(NonOrthIsospinThreeBodySpace), intent(in) :: thr
    logical, intent(in) :: Scalar
    character(*), intent(in) :: optr
    integer, intent(in) :: jr, pr, tr
    integer :: chbra, chket
    integer :: jbra, pbra, tbra, nbra
    integer :: jket, pket, tket, nket

    this%optr = optr
    this%Scalar = Scalar
    this%jr = jr
    this%pr = pr
    this%tr = tr
    this%zr = 0
    this%NChan = thr%NChan

    allocate(this%MatCh(this%NChan,this%NChan))

    do chbra = 1, this%NChan
      jbra = thr%jpt(chbra)%j
      pbra = thr%jpt(chbra)%p
      tbra = thr%jpt(chbra)%t
      nbra = thr%jpt(chbra)%nst
      do chket = 1, this%NChan
        jket = thr%jpt(chket)%j
        pket = thr%jpt(chket)%p
        tket = thr%jpt(chket)%t
        nket = thr%jpt(chket)%nst

        if(triag(jbra, jket, 2*jr)) cycle
        if(pbra * pket * pr /= 1) cycle
        if(triag(tbra, tket, 2*tr)) cycle
        this%MatCh(chbra,chket)%is = .true.
        this%MatCh(chbra,chket)%ndims = [nbra,nket]
        call this%MatCh(chbra,chket)%zeros(nbra,nket)
      end do
    end do
  end subroutine InitNonOrthIsospinThreeBodyPartSp

  subroutine CopyNBodyPartSp(a, b)
    class(NBodyPartSp), intent(inout) :: a
    type(NBodyPartSp), intent(in) :: b
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
        a%MatCh(chbra,chket)%SMat = b%MatCh(chbra,chket)%SMat
      end do
    end do
  end subroutine CopyNBodyPartSp

  subroutine CopyNBodyPart(a, b)
    class(NBodyPart), intent(out) :: a
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
    class(NBodyPart), intent(in) :: a, b
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

  function SumNBodyPartSp(a, b) result(c)
    class(NBodyPartSp), intent(in) :: a, b
    type(NBodyPartSp) :: c
    integer :: chbra, chket
    c = a
    do chbra = 1, b%NChan
      do chket = 1, b%NChan
        if(.not. b%MatCh(chbra,chket)%is) cycle
        c%MatCh(chbra,chket)%is = b%MatCh(chbra,chket)%is
        c%MatCh(chbra,chket)%ndims(:) = b%MatCh(chbra,chket)%ndims(:)
        c%MatCh(chbra,chket)%SMat = &
            & a%MatCh(chbra,chket)%SMat + b%MatCh(chbra,chket)%SMat
      end do
    end do
  end function SumNBodyPartSp

  function SubtractNBodyPart(a, b) result(c)
    class(NBodyPart), intent(in) :: a, b
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

  function SubtractNBodyPartSp(a, b) result(c)
    class(NBodyPartSp), intent(in) :: a, b
    type(NBodyPartSp) :: c
    integer :: chbra, chket
    c = a
    do chbra = 1, b%NChan
      do chket = 1, b%NChan
        if(.not. b%MatCh(chbra,chket)%is) cycle
        c%MatCh(chbra,chket)%is = b%MatCh(chbra,chket)%is
        c%MatCh(chbra,chket)%ndims(:) = b%MatCh(chbra,chket)%ndims(:)
        c%MatCh(chbra,chket)%SMat = &
            & a%MatCh(chbra,chket)%SMat - b%MatCh(chbra,chket)%SMat
      end do
    end do
  end function SubtractNBodyPartSp

  function ScaleNBodyPart(a, b) result(c)
    class(NBodyPart), intent(in) :: a
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
  end function ScaleNBodyPart

  function ScaleNBodyPartSp(a, b) result(c)
    class(NBodyPartSp), intent(in) :: a
    real(8), intent(in) :: b
    type(NBodyPartSp) :: c
    integer :: chbra, chket
    c = a
    do chbra = 1, a%NChan
      do chket = 1, a%NChan
        if(.not. a%MatCh(chbra,chket)%is) cycle
        c%MatCh(chbra,chket)%is = a%MatCh(chbra,chket)%is
        c%MatCh(chbra,chket)%ndims(:) = a%MatCh(chbra,chket)%ndims(:)
        c%MatCh(chbra,chket)%SMat = &
            & a%MatCh(chbra,chket)%SMat * real(b)
      end do
    end do
  end function ScaleNBodyPartSp

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
    Z12 = (sps%orb(i1)%z + sps%orb(i2)%z)/2
    Z34 = (sps%orb(i3)%z + sps%orb(i4)%z)/2

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
    Z12 = (sps%orb(i1)%z + sps%orb(i2)%z)/2
    Z34 = (sps%orb(i3)%z + sps%orb(i4)%z)/2

    if(P12 * P34 /= 1) stop "Error, in GetTwBME_scalar: P"
    if(Z12 - Z34 /= 0) stop "Error, in GetTwBME_scalar: Tz"

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
    Z12 = (sps%orb(i1)%z + sps%orb(i2)%z) / 2
    Z34 = (sps%orb(i3)%z + sps%orb(i4)%z) / 2

    if(triag(J12, J34, this%jr)) stop "Error, in SetTwBME_general: J"
    if(P12 * P34 * this%pr /= 1) stop "Error, in SetTwBME_general: P"
    if(Z12 - Z34 - this%zr /= 0) stop "Error, in SetTwBME_general: Tz"

    chbra = ms%jpz2ch(J12,P12,Z12)
    chket = ms%jpz2ch(J34,P34,Z34)
    if(chbra * chket == 0) return

    bra = ms%jpz(chbra)%spis2n(i1,i2)
    ket = ms%jpz(chket)%spis2n(i3,i4)
    if(bra * ket == 0) return

    iphase = ms%jpz(chbra)%iphase(i1,i2) * &
        &    ms%jpz(chket)%iphase(i3,i4)

    this%MatCh(chbra,chket)%m(bra,ket) = dble(iphase) * me
    this%MatCh(chket,chbra)%m(ket,bra) = dble(iphase) * me
  end subroutine SetTwBME_general

  subroutine SetTwBME_scalar(this,sps,ms,i1,i2,i3,i4,J,me)
    class(NBodyPart), intent(inout) :: this
    type(Orbits), intent(in) :: sps
    type(TwoBodySpace), intent(in) :: ms
    integer, intent(in) :: i1,i2,i3,i4,J
    real(8), intent(in) :: me
    integer :: P12, P34
    integer :: Z12, Z34
    integer :: ch, bra, ket, iphase

    P12 = (-1) ** (sps%orb(i1)%l + sps%orb(i2)%l)
    P34 = (-1) ** (sps%orb(i3)%l + sps%orb(i4)%l)
    Z12 = (sps%orb(i1)%z + sps%orb(i2)%z) / 2
    Z34 = (sps%orb(i3)%z + sps%orb(i4)%z) / 2

    if(P12 * P34 /= 1) stop "Error, in SetTwBME_general: P"
    if(Z12 - Z34 /= 0) stop "Error, in SetTwBME_general: Tz"

    ch = ms%jpz2ch(J,P12,Z12)
    if(ch == 0) return

    bra = ms%jpz(ch)%spis2n(i1,i2)
    ket = ms%jpz(ch)%spis2n(i3,i4)
    if(bra * ket == 0) return

    iphase = ms%jpz(ch)%iphase(i1,i2) * &
        &    ms%jpz(ch)%iphase(i3,i4)

    this%MatCh(ch,ch)%m(bra,ket) = dble(iphase) * me
    this%MatCh(ch,ch)%m(ket,bra) = dble(iphase) * me
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
    Z12 = (sps%orb(i1)%z + sps%orb(i2)%z) / 2
    Z34 = (sps%orb(i3)%z + sps%orb(i4)%z) / 2

    if(triag(J12, J34, this%jr)) stop "Error, in AddToTwBME_general: J"
    if(P12 * P34 * this%pr /= 1) stop "Error, in AddToTwBME_general: P"
    if(Z12 - Z34 - this%zr /= 0) stop "Error, in AddToTwBME_general: Tz"

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
    this%MatCh(chket,chbra)%m(ket,bra) = this%MatCh(chbra,chket)%m(bra,ket)
  end subroutine AddToTwBME_general

  subroutine AddToTwBME_scalar(this,sps,ms,i1,i2,i3,i4,J,me)
    class(NBodyPart), intent(inout) :: this
    type(Orbits), intent(in) :: sps
    type(TwoBodySpace), intent(in) :: ms
    integer, intent(in) :: i1,i2,i3,i4,J
    real(8), intent(in) :: me
    integer :: P12, P34
    integer :: Z12, Z34
    integer :: ch, bra, ket, iphase

    P12 = (-1) ** (sps%orb(i1)%l + sps%orb(i2)%l)
    P34 = (-1) ** (sps%orb(i3)%l + sps%orb(i4)%l)
    Z12 = (sps%orb(i1)%z + sps%orb(i2)%z) / 2
    Z34 = (sps%orb(i3)%z + sps%orb(i4)%z) / 2

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
        & me * dble(iphase)
    this%MatCh(ch,ch)%m(ket,bra) = this%MatCh(ch,ch)%m(bra,ket)
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

  function GetThBMEIso_general(this,sps,ms,i1,i2,i3,J12,T12,&
        & i4,i5,i6,J45,T45,Jbra,Jket,Tbra,Tket) result(r)
    use CommonLibrary, only: triag
    class(NBodyPart), intent(in) :: this
    type(OrbitsIsospin), intent(in) :: sps
    type(NonOrthIsospinThreeBodySpace), intent(in) :: ms
    integer, intent(in) :: i1,i2,i3,i4,i5,i6
    integer, intent(in) :: J12,T12,J45,T45,Jbra,Jket,Tbra,Tket
    integer :: chbra, chket, idxbra, idxket, bra, ket
    integer :: P123, P456
    integer :: isorted_bra, isorted_ket
    integer :: ibra, iket
    real(8) :: r

    r = 0.d0
    P123 = (-1) ** (sps%orb(i1)%l+sps%orb(i2)%l+sps%orb(i3)%l)
    P456 = (-1) ** (sps%orb(i4)%l+sps%orb(i5)%l+sps%orb(i6)%l)
    if(P123 * P456 * this%pr /= 1) stop 'Error in GetThBMEIso_general: P'
    if(triag(Tbra,Tket,2*this%zr)) stop 'Error in GetThBMEIso_general: T'
    chbra = ms%jpt2ch(Jbra,P123,Tbra)
    chket = ms%jpt2ch(Jket,P456,Tket)
    if(chbra*chket == 0) return
    idxbra = ms%jpt(chbra)%sorting(i1,i2,i3)
    idxket = ms%jpt(chket)%sorting(i4,i5,i6)
    if(idxbra * idxket == 0) return
    if(i1 == i2 .and. mod(J12+T12,2) == 0) return
    if(i4 == i5 .and. mod(J45+T45,2) == 0) return
    isorted_bra = ms%jpt(chbra)%sort(idxbra)%idx_sorted
    isorted_ket = ms%jpt(chket)%sort(idxket)%idx_sorted
    if(isorted_bra * isorted_ket == 0) return

    do ibra = 1, ms%jpt(chbra)%sort(idxbra)%JT(J12,T12)%n
      bra = ms%jpt(chbra)%sort(idxbra)%JT(J12,T12)%idx2num(ibra)
      do iket = 1, ms%jpt(chket)%sort(idxket)%JT(J45,T45)%n
        ket = ms%jpt(chket)%sort(idxket)%JT(J45,T45)%idx2num(iket)
        r = r + dble(this%MatCh(chbra,chket)%m(bra,ket) * &
            & ms%jpt(chbra)%sort(idxbra)%JT(J12,T12)%TrnsCoef(ibra) * &
            & ms%jpt(chket)%sort(idxket)%JT(J45,T45)%TrnsCoef(iket))
      end do
    end do
  end function GetThBMEIso_general

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
      if(abs(z1+z2) > 2*T12) cycle
      do T45 = 0, 1
        if(abs(z4+z5) > 2*T45) cycle
        do T = max(abs(2*T12-1),abs(2*T45-1)), min(2*T12+1,2*T45+1), 2
          if(abs(Z) > T) cycle
          ch = ms%thr%jpt2ch(J,P,T)
          if(ch == 0) cycle
          r = r + &
              & this%GetThBME(ms%isps,ms%thr,a,b,c,J12,T12,&
              & d,e,f,J45,T45,J,T) * &
              & dcg(1,z1,1,z2,2*T12,z1+z2) * dcg(2*T12,z1+z2,1,z3,T,Z) * &
              & dcg(1,z4,1,z5,2*T45,z4+z5) * dcg(2*T45,z4+z5,1,z6,T,Z)
        end do
      end do
    end do
  end function GetThBMEpn_scalar

  function GetThBMEpn_general(this,ms,i1,i2,i3,J12,&
        & i4,i5,i6,J45,Jbra,Jket) result(r)
    use CommonLibrary, only: dcg
    class(NBodyPart), intent(in) :: this
    type(MSpace), intent(in) :: ms
    integer, intent(in) :: i1,i2,i3,i4,i5,i6
    integer, intent(in) :: J12,J45,Jbra,Jket
    integer :: z1, z2, z3, z4, z5, z6, T12, T45, Tbra, Tket
    integer :: a, b, c, d, e, f
    integer :: Pbra, Pket, Zbra, Zket, chbra, chket
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

    Pbra = (-1) ** (ms%sps%orb(i1)%l+ms%sps%orb(i2)%l+ms%sps%orb(i3)%l)
    Pket = (-1) ** (ms%sps%orb(i4)%l+ms%sps%orb(i5)%l+ms%sps%orb(i6)%l)
    Zbra = z1 + z2 + z3
    Zket = z4 + z5 + z6
    do T12 = 0, 1
      if(abs(z1+z2) > 2*T12) cycle
      do T45 = 0, 1
        if(abs(z4+z5) > 2*T45) cycle
        do Tbra = abs(2*T12-1), 2*T12+1, 2
            if(abs(Zbra) > Tbra) cycle
          do Tket = abs(2*T45-1), 2*T45+1, 2
            if(abs(Zket) > Tket) cycle
            chbra = ms%thr%jpt2ch(Jbra,Pbra,Tbra)
            chket = ms%thr%jpt2ch(Jket,Pket,Tket)
            if(chbra*chket == 0) cycle
            r = r + &
                & this%GetThBME(ms%isps,ms%thr,a,b,c,J12,T12,&
                & d,e,f,J45,T45,Jbra,Jket,Tbra,Tket) * &
                & dcg(1,z1,1,z2,2*T12,z1+z2) * dcg(2*T12,z1+z2,1,z3,Tbra,Zbra) * &
                & dcg(1,z4,1,z5,2*T45,z4+z5) * dcg(2*T45,z4+z5,1,z6,Tket,Zket)
          end do
        end do
      end do
    end do
  end function GetThBMEpn_general

  function GetThBMEIso_scalar_sp(this,sps,ms,i1,i2,i3,J12,T12,&
        & i4,i5,i6,J45,T45,J,T) result(r)
    class(NBodyPartSp), intent(in) :: this
    type(OrbitsIsospin), intent(in) :: sps
    type(NonOrthIsospinThreeBodySpace), intent(in) :: ms
    integer, intent(in) :: i1,i2,i3,i4,i5,i6
    integer, intent(in) :: J12,T12,J45,T45,J,T
    integer :: ch, idxbra, idxket, bra, ket
    integer :: P123, P456
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
  end function GetThBMEIso_scalar_sp

  function GetThBMEIso_general_sp(this,sps,ms,i1,i2,i3,J12,T12,&
        & i4,i5,i6,J45,T45,Jbra,Jket,Tbra,Tket) result(r)
    use CommonLibrary, only: triag
    class(NBodyPartSp), intent(in) :: this
    type(OrbitsIsospin), intent(in) :: sps
    type(NonOrthIsospinThreeBodySpace), intent(in) :: ms
    integer, intent(in) :: i1,i2,i3,i4,i5,i6
    integer, intent(in) :: J12,T12,J45,T45,Jbra,Jket,Tbra,Tket
    integer :: chbra, chket, idxbra, idxket, bra, ket
    integer :: P123, P456
    integer :: isorted_bra, isorted_ket
    integer :: ibra, iket
    real(8) :: r

    r = 0.d0
    P123 = (-1) ** (sps%orb(i1)%l+sps%orb(i2)%l+sps%orb(i3)%l)
    P456 = (-1) ** (sps%orb(i4)%l+sps%orb(i5)%l+sps%orb(i6)%l)
    if(P123 * P456 * this%pr /= 1) stop 'Error in GetThBMEIso_general: P'
    if(triag(Tbra,Tket,2*this%zr)) stop 'Error in GetThBMEIso_general: T'
    chbra = ms%jpt2ch(Jbra,P123,Tbra)
    chket = ms%jpt2ch(Jket,P456,Tket)
    if(chbra*chket == 0) return
    idxbra = ms%jpt(chbra)%sorting(i1,i2,i3)
    idxket = ms%jpt(chket)%sorting(i4,i5,i6)
    if(idxbra * idxket == 0) return
    if(i1 == i2 .and. mod(J12+T12,2) == 0) return
    if(i4 == i5 .and. mod(J45+T45,2) == 0) return
    isorted_bra = ms%jpt(chbra)%sort(idxbra)%idx_sorted
    isorted_ket = ms%jpt(chket)%sort(idxket)%idx_sorted
    if(isorted_bra * isorted_ket == 0) return

    do ibra = 1, ms%jpt(chbra)%sort(idxbra)%JT(J12,T12)%n
      bra = ms%jpt(chbra)%sort(idxbra)%JT(J12,T12)%idx2num(ibra)
      do iket = 1, ms%jpt(chket)%sort(idxket)%JT(J45,T45)%n
        ket = ms%jpt(chket)%sort(idxket)%JT(J45,T45)%idx2num(iket)
        r = r + dble(this%MatCh(chbra,chket)%m(bra,ket) * &
            & ms%jpt(chbra)%sort(idxbra)%JT(J12,T12)%TrnsCoef(ibra) * &
            & ms%jpt(chket)%sort(idxket)%JT(J45,T45)%TrnsCoef(iket))
      end do
    end do
  end function GetThBMEIso_general_sp

  function GetThBMEpn_scalar_sp(this,ms,i1,i2,i3,J12,&
        & i4,i5,i6,J45,J) result(r)
    use CommonLibrary, only: dcg
    class(NBodyPartSp), intent(in) :: this
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

    a = ms%sps%pn2iso( ms%isps, i1 )
    b = ms%sps%pn2iso( ms%isps, i2 )
    c = ms%sps%pn2iso( ms%isps, i3 )
    d = ms%sps%pn2iso( ms%isps, i4 )
    e = ms%sps%pn2iso( ms%isps, i5 )
    f = ms%sps%pn2iso( ms%isps, i6 )
    if(a*b*c*d*e*f == 0) return

    P = (-1) ** (ms%sps%orb(i1)%l+ms%sps%orb(i2)%l+ms%sps%orb(i3)%l)
    Z = z1 + z2 + z3
    do T12 = 0, 1
      if(abs(z1+z2) > 2*T12) cycle
      do T45 = 0, 1
        if(abs(z4+z5) > 2*T45) cycle
        do T = max(abs(2*T12-1),abs(2*T45-1)), min(2*T12+1,2*T45+1), 2
          if(abs(Z) > T) cycle
          ch = ms%thr%jpt2ch(J,P,T)
          if(ch == 0) cycle
          r = r + &
              & this%GetThBME(ms%isps,ms%thr,a,b,c,J12,T12,&
              & d,e,f,J45,T45,J,T) * &
              & dcg(1,z1,1,z2,2*T12,z1+z2) * dcg(2*T12,z1+z2,1,z3,T,Z) * &
              & dcg(1,z4,1,z5,2*T45,z4+z5) * dcg(2*T45,z4+z5,1,z6,T,Z)
        end do
      end do
    end do

  end function GetThBMEpn_scalar_sp

  function GetThBMEpn_general_sp(this,ms,i1,i2,i3,J12,&
        & i4,i5,i6,J45,Jbra,Jket) result(r)
    use CommonLibrary, only: dcg
    class(NBodyPartSp), intent(in) :: this
    type(MSpace), intent(in) :: ms
    integer, intent(in) :: i1,i2,i3,i4,i5,i6
    integer, intent(in) :: J12,J45,Jbra,Jket
    integer :: z1, z2, z3, z4, z5, z6, T12, T45, Tbra, Tket
    integer :: a, b, c, d, e, f
    integer :: Pbra, Pket, Zbra, Zket, chbra, chket
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

    Pbra = (-1) ** (ms%sps%orb(i1)%l+ms%sps%orb(i2)%l+ms%sps%orb(i3)%l)
    Pket = (-1) ** (ms%sps%orb(i4)%l+ms%sps%orb(i5)%l+ms%sps%orb(i6)%l)
    Zbra = z1 + z2 + z3
    Zket = z4 + z5 + z6
    do T12 = 0, 1
      if(abs(z1+z2) > 2*T12) cycle
      do T45 = 0, 1
        if(abs(z4+z5) > 2*T45) cycle
        do Tbra = abs(2*T12-1), 2*T12+1, 2
            if(abs(Zbra) > Tbra) cycle
          do Tket = abs(2*T45-1), 2*T45+1, 2
            if(abs(Zket) > Tket) cycle
            chbra = ms%thr%jpt2ch(Jbra,Pbra,Tbra)
            chket = ms%thr%jpt2ch(Jket,Pket,Tket)
            if(chbra*chket == 0) cycle
            r = r + &
                & this%GetThBME(ms%isps,ms%thr,a,b,c,J12,T12,&
                & d,e,f,J45,T45,Jbra,Jket,Tbra,Tket) * &
                & dcg(1,z1,1,z2,2*T12,z1+z2) * dcg(2*T12,z1+z2,1,z3,Tbra,Zbra) * &
                & dcg(1,z4,1,z5,2*T45,z4+z5) * dcg(2*T45,z4+z5,1,z6,Tket,Zket)
          end do
        end do
      end do
    end do
  end function GetThBMEpn_general_sp

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

  ! Tensor case should be tested
  function NormalOrderingFromSp3To2(this,ms) result(r)
    use CommonLibrary, only: triag, sjs
    class(NBodyPartSp), intent(in) :: this
    type(MSpace), intent(in) :: ms
    type(NBodyPart) :: r
    integer :: i1, i2, i3, i4, ih, J12, J34, jh, Jbra, Jket
    integer :: e1, e2, e3, e4, eh
    integer :: chbra, chket, bra, ket, bra_max, ket_max
    real(8) :: v, fact, vsum, tfact

    call r%init(ms%two,this%Scalar,this%optr,this%jr,this%pr,this%zr)

    do chbra = 1, ms%two%NChan
      do chket = 1, chbra
        if( .not. r%MatCh(chbra,chket)%is) cycle

        bra_max = ms%two%jpz(chbra)%nst
        ket_max = ms%two%jpz(chket)%nst

        J12 = ms%two%jpz(chbra)%j
        J34 = ms%two%jpz(chket)%j

        do bra = 1, bra_max
          i1 = ms%two%jpz(chbra)%n2spi1(bra)
          i2 = ms%two%jpz(chbra)%n2spi2(bra)
          e1 = ms%sps%orb(i1)%e
          e2 = ms%sps%orb(i2)%e
          if(chbra == chket) ket_max = bra
          do ket = 1, ket_max
            i3 = ms%two%jpz(chket)%n2spi1(ket)
            i4 = ms%two%jpz(chket)%n2spi2(ket)
            e3 = ms%sps%orb(i3)%e
            e4 = ms%sps%orb(i4)%e

            fact = 1.d0
            if(i1==i2) fact = fact / dsqrt(2.d0)
            if(i3==i4) fact = fact / dsqrt(2.d0)

            vsum = 0.d0
            do ih = 1, ms%sps%norbs
              if(abs(ms%NOCoef(ih)) < 1.d-6) cycle
              jh = ms%sps%orb(ih)%j
              eh = ms%sps%orb(ih)%e
              if(e1+e2+eh > ms%e3max) cycle
              if(e3+e4+eh > ms%e3max) cycle

              v = 0.d0
              do Jbra = abs(2*J12-jh), (2*J12+jh), 2
                do Jket = abs(2*J34-jh), (2*J34+jh), 2
                  if(triag(Jbra,Jket,2*this%jr)) cycle
                  tfact = sqrt(dble( (Jbra+1) * (Jket+1) ))
                  if(.not. this%Scalar) then
                    tfact = sqrt( dble( (Jbra+1) * (Jket+1) ) ) * &
                        & (-1.d0) ** ( J12 + (jh+Jket)/2 + this%jr) * &
                        & sjs(2*J12,Jbra,jh,Jket,2*J34,2*this%jr)
                  end if
                  v = v + tfact * ms%NOCoef(ih) * &
                      & this%GetThBME(ms,i1,i2,ih,J12,i3,i4,ih,J34,Jbra,Jket)
                end do
              end do
            end do
            r%MatCh(chbra,chket)%m(bra,ket) = vsum * fact / sqrt( dble( (2*J12+1) * (2*J34+1) ) )
            r%MatCh(chket,chbra)%m(ket,bra) = r%MatCh(chbra,chket)%m(bra,ket) * &
                & (-1.d0) ** (J12-J34)
          end do
        end do

      end do
    end do

  end function NormalOrderingFromSp3To2

  ! Tensor case should be tested
  function NormalOrderingFrom3To2(this,ms) result(r)
    use CommonLibrary, only: triag, sjs
    class(NBodyPart), intent(in) :: this
    type(MSpace), intent(in) :: ms
    type(NBodyPart) :: r
    integer :: i1, i2, i3, i4, ih, J12, J34, jh, Jbra, Jket
    integer :: e1, e2, e3, e4, eh
    integer :: chbra, chket, bra, ket, bra_max, ket_max
    real(8) :: v, fact, vsum, tfact

    call r%init(ms%two,this%Scalar,this%optr,this%jr,this%pr,this%zr)

    do chbra = 1, ms%two%NChan
      do chket = 1, chbra
        if( .not. r%MatCh(chbra,chket)%is) cycle

        bra_max = ms%two%jpz(chbra)%nst
        ket_max = ms%two%jpz(chket)%nst

        J12 = ms%two%jpz(chbra)%j
        J34 = ms%two%jpz(chket)%j

        do bra = 1, bra_max
          i1 = ms%two%jpz(chbra)%n2spi1(bra)
          i2 = ms%two%jpz(chbra)%n2spi2(bra)
          e1 = ms%sps%orb(i1)%e
          e2 = ms%sps%orb(i2)%e
          if(chbra == chket) ket_max = bra
          do ket = 1, ket_max
            i3 = ms%two%jpz(chket)%n2spi1(ket)
            i4 = ms%two%jpz(chket)%n2spi2(ket)
            e3 = ms%sps%orb(i3)%e
            e4 = ms%sps%orb(i4)%e

            fact = 1.d0
            if(i1==i2) fact = fact / dsqrt(2.d0)
            if(i3==i4) fact = fact / dsqrt(2.d0)

            vsum = 0.d0
            do ih = 1, ms%sps%norbs
              if(abs(ms%NOCoef(ih)) < 1.d-6) cycle
              jh = ms%sps%orb(ih)%j
              eh = ms%sps%orb(ih)%e
              if(e1+e2+eh > ms%e3max) cycle
              if(e3+e4+eh > ms%e3max) cycle

              v = 0.d0
              do Jbra = abs(2*J12-jh), (2*J12+jh), 2
                do Jket = abs(2*J34-jh), (2*J34+jh), 2
                  if(triag(Jbra,Jket,2*this%jr)) cycle
                  tfact = sqrt(dble( (Jbra+1) * (Jket+1) ))
                  if(.not. this%Scalar) then
                    tfact = sqrt( dble( (Jbra+1) * (Jket+1) ) ) * &
                        & (-1.d0) ** ( J12 + (jh+Jket)/2 + this%jr) * &
                        & sjs(2*J12,Jbra,jh,Jket,2*J34,2*this%jr)
                  end if
                  v = v + tfact * ms%NOCoef(ih) * &
                      & this%GetThBME(ms,i1,i2,ih,J12,i3,i4,ih,J34,Jbra,Jket)
                end do
              end do
            end do
            r%MatCh(chbra,chket)%m(bra,ket) = vsum * fact / sqrt( dble( (2*J12+1) * (2*J34+1) ) )
            r%MatCh(chket,chbra)%m(ket,bra) = r%MatCh(chbra,chket)%m(bra,ket) * &
                & (-1.d0) ** (J12-J34)
          end do
        end do

      end do
    end do

  end function NormalOrderingFrom3To2

  ! Tensor case should be tested
  function NormalOrderingFrom2To1(this,ms) result(r)
    use CommonLibrary, only: sjs, triag
    class(NBodyPart), intent(in) :: this
    type(MSpace), intent(in) :: ms
    type(NBodyPart) :: r
    integer :: i1, i2, ih, j1, j2, l1, l2, z1, z2, e1, e2
    integer :: jh, eh, J_bra, J_ket
    integer :: chbra, chket, bra, ket, bra_max, ket_max
    real(8) :: vsum, v, fact, tfact

    call r%init(ms%one,this%Scalar,this%optr,this%jr,this%pr,this%zr)

    do chbra = 1, ms%one%NChan
      do chket = 1, chbra
        if( .not. r%MatCh(chbra,chket)%is) cycle


        bra_max = ms%one%jpz(chbra)%nst
        ket_max = ms%one%jpz(chket)%nst

        j1 = ms%one%jpz(chbra)%j
        j2 = ms%one%jpz(chket)%j

        do bra = 1, bra_max
          i1 = ms%one%jpz(chbra)%n2spi(bra)
          l1 = ms%sps%orb(i1)%l
          z1 = ms%sps%orb(i1)%z
          e1 = ms%sps%orb(i1)%e
          if(chbra == chket) ket_max = bra
          do ket = 1, ket_max
            i2 = ms%one%jpz(chket)%n2spi(ket)
            l2 = ms%sps%orb(i2)%l
            z2 = ms%sps%orb(i2)%z
            e2 = ms%sps%orb(i2)%e

            vsum = 0.d0
            do ih = 1, ms%sps%norbs
              if(abs(ms%NOCoef(ih)) < 1.d-6) cycle
              jh = ms%sps%orb(ih)%j
              eh = ms%sps%orb(ih)%e
              if(e1 + eh > ms%e2max) cycle
              if(e2 + eh > ms%e2max) cycle
              fact = 1.d0
              if(i1==ih) fact = fact * dsqrt(2.d0)
              if(i2==ih) fact = fact * dsqrt(2.d0)

              v = 0.d0
              do J_bra = abs(j1-jh)/2, (j1+jh)/2
                if(i1 == ih .and. mod(J_bra,2)==1) cycle
                do J_ket = abs(j2-jh)/2, (j2+jh)/2
                  if(i2 == ih .and. mod(J_bra,2)==1) cycle

                  if(triag(J_bra,J_ket,this%jr)) cycle
                  tfact = sqrt(dble( (2*J_bra+1) * (2*J_ket+1) ) )
                  if(.not. this%Scalar) then
                    tfact = sqrt(dble( (2*J_bra+1) * (2*J_ket+1) ) ) * &
                        (-1.d0) ** ((j1+jh)/2+J_ket+this%jr) * &
                        & sjs(j1,2*J_bra,jh,2*J_ket,j2,2*this%jr)
                  end if
                  v = v + tfact * &
                      & this%GetTwBME(ms%sps,ms%two,i1,ih,i2,ih,J_bra,J_ket) * &
                      & ms%NOcoef(ih)
                end do
              end do
              vsum = vsum + v * fact
            end do
            r%MatCh(chbra,chket)%m(bra,ket) = vsum
            if(this%Scalar) r%MatCh(chbra,chket)%m(bra,ket) = &
                & vsum / dble(j1+1)
          end do
        end do

      end do
    end do
  end function NormalOrderingFrom2To1

  function NormalOrderingFrom1To0(this,ms) result(r)
    class(NBodyPart), intent(in) :: this
    type(MSpace), intent(in) :: ms
    real(8) :: r
    integer :: ih, lh, jh, zh, ch, nh

    r = 0.d0
    if(.not. this%Scalar) return
    do ih = 1, ms%sps%norbs
      if(abs(ms%NOCoef(ih)) < 1.d-6) cycle
      jh = ms%sps%orb(ih)%j
      lh = ms%sps%orb(ih)%l
      zh = ms%sps%orb(ih)%z
      ch = ms%one%jpz2ch(jh,(-1)**lh,zh)
      nh = ms%one%jpz(ch)%spi2n(ih)
      r = r + dble(jh+1) * ms%NOCoef(ih) * &
          & this%MatCh(ch,ch)%m(nh,nh)
    end do
  end function NormalOrderingFrom1To0

  subroutine PrintNBodyPartSp(this, wunit)
    use ClassSys, only: sys
    class(NBodyPartSp), intent(in) :: this
    integer, intent(in), optional :: wunit
    integer :: chbra, chket
    type(sys) :: s
    character(:), allocatable :: msg


    do chbra = 1, this%NChan
      do chket = 1, this%NChan
        if(.not. this%MatCh(chbra,chket)%is) cycle
        msg = trim(this%optr) // " " // trim(s%str(chbra)) &
            &  // " " // trim(s%str(chket))
        call this%MatCh(chbra,chket)%prt(msg=msg,iunit=wunit)
      end do
    end do
  end subroutine PrintNBodyPartSp

  subroutine PrintNBodyPart(this, wunit)
    use ClassSys, only: sys
    class(NBodyPart), intent(in) :: this
    integer, intent(in), optional :: wunit
    integer :: chbra, chket
    type(sys) :: s
    character(:), allocatable :: msg

    do chbra = 1, this%NChan
      do chket = 1, this%NChan
        if(.not. this%MatCh(chbra,chket)%is) cycle
        msg = trim(this%optr) // " " // trim(s%str(chbra)) &
            &  // " " // trim(s%str(chket))
        call this%MatCh(chbra,chket)%prt(msg=msg,iunit=wunit)
      end do
    end do
  end subroutine PrintNBodyPart

  !
  !
  !     File reading methods
  !
  !

  subroutine ReadTwoBodyFile(this, sps, ms, rd)
    use Profiler, only: timer
    class(NBodyPart), intent(inout) :: this
    type(Orbits), intent(in) :: sps
    type(TwoBodySpace), intent(in) :: ms
    type(ReadFiles), intent(in) :: rd
    real(8) :: ti

    ti = omp_get_wtime()

    select case(rd%file_nn)
    case('None', 'NONE', 'none')
      write(*,*) "Error in ReadTwoBodyFile"

    case default

      if(this%Scalar) call rd%ReadScalar2BFile(this,sps,ms)
      if(.not. this%Scalar) call rd%ReadTensor2BFile(this,sps,ms)
    end select

    call timer%Add("Read from file", omp_get_wtime()-ti)

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
    integer :: a, b, c, d, dmax
    integer :: la, ja, ea
    integer :: lb, jb, eb
    integer :: lc, jc, ec
    integer :: ld, jd, ed
    integer :: J, Jmin, Jmax, icnt
    integer :: ap, bp, cp, dp
    integer :: an, bn, cn, dn
    real(8) :: me_00, me_pp, me_10, me_nn, fact

    call sps_me2j%init(this%emax2, this%lmax2)
    nelm = count_scalar_me2j(sps_me2j,this%e2max2)
    allocate(v(nelm))
    open(runit, file=this%file_nn, action='read',iostat=io)
    if(io /= 0) then
      write(*,'(2a)') 'File open error: ', trim(this%file_nn)
      return
    end if
    call get_vector_me2j_formatted(runit,v)
    close(runit)

    icnt = 0
    do a = 1, sps_me2j%norbs
      la = sps_me2j%orb(a)%l
      ja = sps_me2j%orb(a)%j
      ea = sps_me2j%orb(a)%e
      ap = sps_me2j%iso2pn(sps,a,-1)
      an = sps_me2j%iso2pn(sps,a, 1)
      do b = 1, a
        lb = sps_me2j%orb(b)%l
        jb = sps_me2j%orb(b)%j
        eb = sps_me2j%orb(b)%e
        bp = sps_me2j%iso2pn(sps,b,-1)
        bn = sps_me2j%iso2pn(sps,b, 1)
        if(ea + eb > this%e2max2) cycle
        do c = 1, a
          lc = sps_me2j%orb(c)%l
          jc = sps_me2j%orb(c)%j
          ec = sps_me2j%orb(c)%e
          cp = sps_me2j%iso2pn(sps,c,-1)
          cn = sps_me2j%iso2pn(sps,c, 1)
          dmax = c
          if(a == c) dmax = b
          do d = 1, dmax
            ld = sps_me2j%orb(d)%l
            jd = sps_me2j%orb(d)%j
            ed = sps_me2j%orb(d)%e
            dp = sps_me2j%iso2pn(sps,d,-1)
            dn = sps_me2j%iso2pn(sps,d, 1)
            if(ec + ed > this%e2max2) cycle
            if( (-1)**(la+lb) /= (-1)**(lc+ld) ) cycle
            Jmin = max(abs(ja-jb), abs(jc-jd))/2
            Jmax = min(   (ja+jb),    (jc+jd))/2
            if(Jmin > Jmax) cycle
            do J = Jmin, Jmax
              me_00 = v(icnt+1)
              me_nn = v(icnt+2)
              me_10 = v(icnt+3)
              me_pp = v(icnt+4)
              icnt = icnt + 4

              if(ap == 0) cycle
              if(bp == 0) cycle
              if(cp == 0) cycle
              if(dp == 0) cycle

              if(sps%orb(ap)%e + sps%orb(bp)%e > ms%e2max) cycle
              if(sps%orb(cp)%e + sps%orb(dp)%e > ms%e2max) cycle
#ifdef NOperatorsDebug
              write(*,'(5i3,4f12.6)') a, b, c, d, J, &
                  &  me_00, me_nn, me_10, me_pp
#endif
              fact = 1.d0
              if(a == b) fact = fact / dsqrt(2.d0)
              if(c == d) fact = fact / dsqrt(2.d0)

              if(abs(1.d0-fact) < 1.d-4 .or. mod(J,2) == 0) then
                call two%SetTwBME(sps,ms,ap,bp,cp,dp,J,me_pp*fact)
                call two%SetTwBME(sps,ms,an,bn,cn,dn,J,me_nn*fact)
                call two%AddToTwBME(sps,ms,ap,bn,cp,dn,J,0.5d0*me_10) ! pnpn
                if(c/=d) call two%AddToTwBME(sps,ms,ap,bn,cn,dp,J,0.5d0*me_10) ! pnnp
                if(a/=b .and. c/=d) &
                    &    call two%AddToTwBME(sps,ms,an,bp,cn,dp,J,0.5d0*me_10) ! npnp
                if(a/=b .and. (a/=c .or. b/=d)) &
                    &    call two%AddToTwBME(sps,ms,an,bp,cp,dn,J,0.5d0*me_10) ! nppn
              end if

              if(abs(1.d0-fact) < 1.d-4 .or. mod(J,2) == 1) then
                call two%AddToTwBME(sps,ms,ap,bn,cp,dn,J,0.5d0*me_00) ! pnpn
                if(c/=d) call two%AddToTwBME(sps,ms,ap,bn,cn,dp,J,-0.5d0*me_00) ! pnnp
                if(a/=b .and. c/=d) &
                    &    call two%AddToTwBME(sps,ms,an,bp,cn,dp,J, 0.5d0*me_00) ! npnp
                if(a/=b .and. (a/=c .or. b/=d)) &
                    &    call two%AddToTwBME(sps,ms,an,bp,cp,dn,J,-0.5d0*me_00) ! nppn
              end if
            end do
          end do
        end do
      end do
    end do

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
#ifdef single_precision
    real(4), allocatable :: v(:)
#else
    real(8), allocatable :: v(:)
#endif
    integer :: a, b, c, d, dmax
    integer :: la, ja, ea
    integer :: lb, jb, eb
    integer :: lc, jc, ec
    integer :: ld, jd, ed
    integer :: J, Jmin, Jmax, icnt
    integer :: ap, bp, cp, dp
    integer :: an, bn, cn, dn
    real(8) :: me_00, me_pp, me_10, me_nn, fact

    call sps_me2j%init(this%emax2, this%lmax2)
    nelm = count_scalar_me2j(sps_me2j,this%e2max2)
    allocate(v(nelm))
    open(runit, file=this%file_nn, action='read',iostat=io, &
        & form='unformatted', access='stream')
    if(io /= 0) then
      write(*,'(2a)') 'File open error: ', trim(this%file_nn)
      return
    end if
    read(runit) v
    close(runit)

    icnt = 0
    do a = 1, sps_me2j%norbs
      la = sps_me2j%orb(a)%l
      ja = sps_me2j%orb(a)%j
      ea = sps_me2j%orb(a)%e
      ap = sps_me2j%iso2pn(sps,a,-1)
      an = sps_me2j%iso2pn(sps,a, 1)
      do b = 1, a
        lb = sps_me2j%orb(b)%l
        jb = sps_me2j%orb(b)%j
        eb = sps_me2j%orb(b)%e
        bp = sps_me2j%iso2pn(sps,b,-1)
        bn = sps_me2j%iso2pn(sps,b, 1)
        if(ea + eb > this%e2max2) cycle
        do c = 1, a
          lc = sps_me2j%orb(c)%l
          jc = sps_me2j%orb(c)%j
          ec = sps_me2j%orb(c)%e
          cp = sps_me2j%iso2pn(sps,c,-1)
          cn = sps_me2j%iso2pn(sps,c, 1)
          dmax = c
          if(a == c) dmax = b
          do d = 1, dmax
            ld = sps_me2j%orb(d)%l
            jd = sps_me2j%orb(d)%j
            ed = sps_me2j%orb(d)%e
            dp = sps_me2j%iso2pn(sps,d,-1)
            dn = sps_me2j%iso2pn(sps,d, 1)
            if(ec + ed > this%e2max2) cycle
            if( (-1)**(la+lb) /= (-1)**(lc+ld) ) cycle
            Jmin = max(abs(ja-jb), abs(jc-jd))/2
            Jmax = min(   (ja+jb),    (jc+jd))/2
            if(Jmin > Jmax) cycle
            do J = Jmin, Jmax
              me_00 = v(icnt+1)
              me_nn = v(icnt+2)
              me_10 = v(icnt+3)
              me_pp = v(icnt+4)
              icnt = icnt + 4

              if(ap == 0) cycle
              if(bp == 0) cycle
              if(cp == 0) cycle
              if(dp == 0) cycle

              if(sps%orb(ap)%e + sps%orb(bp)%e > ms%e2max) cycle
              if(sps%orb(cp)%e + sps%orb(dp)%e > ms%e2max) cycle

              fact = 1.d0
              if(a == b) fact = fact / dsqrt(2.d0)
              if(c == d) fact = fact / dsqrt(2.d0)

              if(abs(1.d0-fact) < 1.d-4 .or. mod(J,2) == 0) then
                call two%SetTwBME(sps,ms,ap,bp,cp,dp,J,me_pp*fact)
                call two%SetTwBME(sps,ms,an,bn,cn,dn,J,me_nn*fact)
                call two%AddToTwBME(sps,ms,ap,bn,cp,dn,J,0.5d0*me_10)
                if(c/=d) call two%AddToTwBME(sps,ms,ap,bn,cn,dp,J,0.5d0*me_10)
                if(a/=b .and. c/=d) &
                    &    call two%AddToTwBME(sps,ms,an,bp,cn,dp,J,0.5d0*me_10)
                if(a/=b .and. (a/=c .or. b/=d)) &
                    &    call two%AddToTwBME(sps,ms,an,bp,cp,dn,J,0.5d0*me_10)
              end if

              if(abs(1.d0-fact) < 1.d-4 .or. mod(J,2) == 1) then
                call two%AddToTwBME(sps,ms,ap,bn,cp,dn,J,0.5d0*me_00)
                if(c/=d) call two%AddToTwBME(sps,ms,ap,bn,cn,dp,J,-0.5d0*me_00)
                if(a/=b .and. c/=d) &
                    &    call two%AddToTwBME(sps,ms,an,bp,cn,dp,J, 0.5d0*me_00)
                if(a/=b .and. (a/=c .or. b/=d)) &
                    &    call two%AddToTwBME(sps,ms,an,bp,cp,dn,J,-0.5d0*me_00)
              end if
            end do
          end do
        end do
      end do
    end do
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

  subroutine get_vector_me2j_formatted(ut,v)
    integer, intent(in) :: ut
    real(8), intent(inout) :: v(:)
    integer :: nelm, lines, i
    nelm = size(v)
    lines = nelm/10
    read(ut,*) ! header
    do i = 1, lines
      read(ut,*) v( (i-1)*10+1 : i*10)
    end do
    read(ut,*) v(lines*10+1:nelm)
  end subroutine get_vector_me2j_formatted

  subroutine read_scalar_myg_txt(this,two,sps,ms)
    class(ReadFiles), intent(in) :: this
    type(NBodyPart), intent(inout) :: two
    type(Orbits), intent(in) :: sps
    type(TwoBodySpace), intent(in) :: ms
    integer :: runit = 22, io
    integer :: lines, n
    integer :: a, b, c, d, J
    real(8) :: me

    open(runit,file=this%file_nn, action='read', iostat=io)
    if(io /= 0) then
      write(*,'(2a)') 'File open error: ', trim(this%file_nn)
      return
    end if

    read(runit,*) lines ! header
    read(runit,*)       ! header
    do n = 1, lines
      read(runit,*) a, b, c, d, J, me
#ifdef NOperatorsDebug
      write(*,'(5i3,f12.6)') a, b, c, d, J, me
#endif
      if(a > sps%norbs) cycle
      if(b > sps%norbs) cycle
      if(c > sps%norbs) cycle
      if(d > sps%norbs) cycle
      call two%SetTwBME(sps,ms,a,b,c,d,J,me)
    end do
    close(runit)
  end subroutine read_scalar_myg_txt

  subroutine read_scalar_myg_bin(this,two,sps,ms)
    class(ReadFiles), intent(in) :: this
    type(NBodyPart), intent(inout) :: two
    type(Orbits), intent(in) :: sps
    type(TwoBodySpace), intent(in) :: ms
    integer :: runit = 22, io
    integer :: lines, n
    integer, allocatable :: a(:), b(:), c(:), d(:), J(:)
    real(8), allocatable :: me(:)

    open(runit,file=this%file_nn, action='read', iostat=io, form='unformatted',&
        & access='stream')
    if(io /= 0) then
      write(*,'(2a)') 'File open error: ', trim(this%file_nn)
      return
    end if

    read(runit) lines ! header
    allocate(a(lines))
    allocate(b(lines))
    allocate(c(lines))
    allocate(d(lines))
    allocate(J(lines))
    allocate(me(lines))
    read(runit) a, b, c, d, J, me
    close(runit)

    !$omp parallel
    !$omp do private(n)
    do n = 1, lines
      if(a(n) > sps%norbs) cycle
      if(b(n) > sps%norbs) cycle
      if(c(n) > sps%norbs) cycle
      if(d(n) > sps%norbs) cycle
      call two%SetTwBME(sps,ms,a(n),b(n),c(n),d(n),J(n),me(n))
    end do
    !$omp end do
    !$omp end parallel
  end subroutine read_scalar_myg_bin

  subroutine read_scalar_snt_txt(this,two,sps,ms)
    use CommonLibrary, only: skip_comment
    class(ReadFiles), intent(in) :: this
    type(NBodyPart), intent(inout) :: two
    type(Orbits), intent(in) :: sps
    type(TwoBodySpace), intent(in) :: ms
    type(SingleParticleOrbit), allocatable :: ssnt(:)
    integer :: a, b, c, d, J, lines, idx, n, l, z
    integer :: ia, ib, ic, id
    integer :: runit=22, io
    integer :: porbs, norbs, pc, nc, line
    real(8) :: me

    open(runit, file=this%file_nn, action='read',iostat=io)
    if(io /= 0) then
      write(*,'(2a)') 'File open error: ', trim(this%file_nn)
      return
    end if
    call skip_comment(runit,'#')
    read(runit,*) porbs, norbs, pc, nc
    lines = porbs+norbs
    allocate(ssnt(lines))
    do line = 1, lines
      read(runit,*) idx, n, l, j, z
      call ssnt(line)%set(n,l,j,z,idx)
    end do

    call skip_comment(runit,'#')
    read(runit,*)  me ! zero-body part
    call skip_comment(runit,'#')
    read(runit,*) lines
    call skip_comment(runit,'#')

    do line = 1, lines
      read(runit,*) a, b, me ! one-body part
    end do

    call skip_comment(runit,'#')
    read(runit,*) lines
    call skip_comment(runit,'#')
    do line = 1, lines
      read(runit,*) a, b, c, d, J, me
#ifdef NOperatorsDebug
      write(*,'(5i3,f12.6)') a, b, c, d, J, me
#endif
      ia = sps%nljz2idx(ssnt(a)%n, ssnt(a)%l, ssnt(a)%j, ssnt(a)%z)
      ib = sps%nljz2idx(ssnt(b)%n, ssnt(b)%l, ssnt(b)%j, ssnt(b)%z)
      ic = sps%nljz2idx(ssnt(c)%n, ssnt(c)%l, ssnt(c)%j, ssnt(c)%z)
      id = sps%nljz2idx(ssnt(d)%n, ssnt(d)%l, ssnt(d)%j, ssnt(d)%z)
      call two%SetTwBME(sps,ms,ia,ib,ic,id,J,me)
    end do
    close(runit)
  end subroutine read_scalar_snt_txt

  subroutine read_scalar_nv_txt(this,two,sps,ms)
    class(ReadFiles), intent(in) :: this
    type(NBodyPart), intent(inout) :: two
    type(Orbits), intent(in) :: sps
    type(TwoBodySpace), intent(in) :: ms

    write(*,*) "reading this format has not been implemented yet."
    return
  end subroutine read_scalar_nv_txt

  !
  !
  ! reading two-body tensor
  !
  !

  subroutine ReadTensor2BFile(this,two,sps,ms)
    use ClassSys, only: sys
    class(ReadFiles), intent(in) :: this
    type(NBodyPart), intent(inout) :: two
    type(Orbits), intent(in) :: sps
    type(TwoBodySpace), intent(in) :: ms
    return
  end subroutine ReadTensor2BFile


  !
  !
  ! reading three-body scalar
  !
  !

  subroutine ReadIsospinThreeBodyFile(this, sps, ms, rd)
    use Profiler, only: timer
    class(NBodyPart), intent(inout) :: this
    type(OrbitsIsospin), intent(in) :: sps
    type(NonOrthIsospinThreeBodySpace), intent(in) :: ms
    type(ReadFiles), intent(in) :: rd
    real(8) :: ti

    ti = omp_get_wtime()

    select case(rd%file_3n)
    case('None', 'NONE', 'none')
      write(*,*) "No three-body matrix element."
      return
    case default

      if(this%Scalar) call rd%ReadScalar3BFile(this,sps,ms)
      if(.not. this%Scalar) call rd%ReadTensor3BFile(this,sps,ms)
    end select

    call timer%Add('Read from file', omp_get_wtime()-ti)

  end subroutine ReadIsospinThreeBodyFile

  subroutine ReadIsospinThreeBodyFileSp(this, sps, ms, rd)
    class(NBodyPartSp), intent(inout) :: this
    type(OrbitsIsospin), intent(in) :: sps
    type(NonOrthIsospinThreeBodySpace), intent(in) :: ms
    type(ReadFiles), intent(in) :: rd

    select case(rd%file_3n)
    case('None', 'NONE', 'none')
      write(*,*) "No three-body matrix element."
      return
    case default

      if(this%Scalar) call rd%ReadScalar3BFileSp(this,sps,ms)
      if(.not. this%Scalar) call rd%ReadTensor3BFileSp(this,sps,ms)
    end select

  end subroutine ReadIsospinThreeBodyFileSp

  subroutine ReadScalar3BFile(this,thr,sps,ms)
    use ClassSys, only: sys
    class(ReadFiles), intent(in) :: this
    type(NBodyPart), intent(inout) :: thr
    type(OrbitsIsospin), intent(in) :: sps
    type(NonOrthIsospinThreeBodySpace), intent(in) :: ms
    type(sys) :: s

    write(*,*) trim(this%file_3n)
    if(s%find(this%file_3n,'.txt')) then
      call this%read_scalar_3bme_txt(thr,sps,ms)
      return
    end if

    if(s%find(this%file_3n,'.bin')) then
      call this%read_scalar_3bme_bin(thr,sps,ms)
      return
    end if
  end subroutine ReadScalar3BFile

  subroutine read_scalar_3bme_txt(this,thr,sps,ms)
    class(ReadFiles), intent(in) :: this
    type(NBodyPart), intent(inout) :: thr
    type(OrbitsIsospin), intent(in) :: sps
    type(NonOrthIsospinThreeBodySpace), intent(in) :: ms
    type(OrbitsIsospin) :: spsf
    real(8), allocatable :: v(:)
    integer(8) :: nelm, n
    integer :: runit = 22, io
    call spsf%init(this%emax3, this%lmax3)
    nelm = count_scalar_3bme(spsf, this%e2max3, this%e3max3)
    allocate(v(nelm))
    open(runit, file=this%file_3n, action='read', iostat=io)
    if(io /= 0) then
      write(*,'(2a)') "File opening error: ", trim(this%file_3n)
      return
    end if
    do n = 1, nelm
      read(runit,*) v(n)
    end do
    close(runit)

    call store_scalar_3bme(thr,sps,ms,v,spsf,this%e2max3,this%e3max3)

    deallocate(v)
    call spsf%fin()
  end subroutine read_scalar_3bme_txt

  subroutine read_scalar_3bme_bin(this,thr,sps,ms)
    class(ReadFiles), intent(in) :: this
    type(NBodyPart), intent(inout) :: thr
    type(OrbitsIsospin), intent(in) :: sps
    type(NonOrthIsospinThreeBodySpace), intent(in) :: ms
    type(OrbitsIsospin) :: spsf
    real(8), allocatable :: v(:)
    integer(8) :: nelm, n
    integer :: runit = 22, io
    call spsf%init(this%emax3, this%lmax3)
    nelm = count_scalar_3bme(spsf, this%e2max3, this%e3max3)
    allocate(v(nelm))
    open(runit, file=this%file_3n, action='read', iostat=io)
    if(io /= 0) then
      write(*,'(2a)') "File opening error: ", trim(this%file_3n)
      return
    end if
    do n = 1, nelm
      read(runit) v
    end do
    close(runit)

    call store_scalar_3bme(thr,sps,ms,v,spsf,this%e2max3,this%e3max3)

    deallocate(v)
    call spsf%fin()
  end subroutine read_scalar_3bme_bin

  subroutine ReadScalar3BFileSp(this,thr,sps,ms)
    use ClassSys, only: sys
    class(ReadFiles), intent(in) :: this
    type(NBodyPartSp), intent(inout) :: thr
    type(OrbitsIsospin), intent(in) :: sps
    type(NonOrthIsospinThreeBodySpace), intent(in) :: ms
    type(sys) :: s

    if(s%find(this%file_3n,'.txt')) then
      call this%read_scalar_3bme_txt_sp(thr,sps,ms)
      return
    end if

    if(s%find(this%file_3n,'.bin')) then
      call this%read_scalar_3bme_bin_sp(thr,sps,ms)
      return
    end if
  end subroutine ReadScalar3BFileSp

  subroutine read_scalar_3bme_txt_sp(this,thr,sps,ms)
    class(ReadFiles), intent(in) :: this
    type(NBodyPartSp), intent(inout) :: thr
    type(OrbitsIsospin), intent(in) :: sps
    type(NonOrthIsospinThreeBodySpace), intent(in) :: ms
    type(OrbitsIsospin) :: spsf
    real(4), allocatable :: v(:)
    integer(8) :: nelm, n
    integer :: runit = 22, io
    call spsf%init(this%emax3, this%lmax3)
    nelm = count_scalar_3bme(spsf, this%e2max3, this%e3max3)
    allocate(v(nelm))
    open(runit, file=this%file_3n, action='read', iostat=io)
    if(io /= 0) then
      write(*,'(2a)') "File opening error: ", trim(this%file_3n)
      return
    end if
    do n = 1, nelm
      read(runit,*) v(n)
    end do
    close(runit)

    call store_scalar_3bme_sp(thr,sps,ms,v,spsf,this%e2max3,this%e3max3)

    deallocate(v)
    call spsf%fin()
  end subroutine read_scalar_3bme_txt_sp

  subroutine read_scalar_3bme_bin_sp(this,thr,sps,ms)
    class(ReadFiles), intent(in) :: this
    type(NBodyPartSp), intent(inout) :: thr
    type(OrbitsIsospin), intent(in) :: sps
    type(NonOrthIsospinThreeBodySpace), intent(in) :: ms
    type(OrbitsIsospin) :: spsf
    real(4), allocatable :: v(:)
    integer(8) :: nelm
    integer :: runit = 22, io
    call spsf%init(this%emax3, this%lmax3)
    nelm = count_scalar_3bme(spsf, this%e2max3, this%e3max3)
    allocate(v(nelm))
    open(runit, file=this%file_3n, action='read', iostat=io, &
        & form='unformatted', access='stream')
    if(io /= 0) then
      write(*,'(2a)') "File opening error: ", trim(this%file_3n)
      return
    end if
    read(runit) v
    close(runit)

    call store_scalar_3bme_sp(thr,sps,ms,v,spsf,this%e2max3,this%e3max3)

    deallocate(v)
    call spsf%fin()
  end subroutine read_scalar_3bme_bin_sp

  function count_scalar_3bme(spsf, e2max, e3max) result(r)
    type(OrbitsIsospin), intent(in) :: spsf
    integer, intent(in) :: e2max, e3max
    integer(8) :: r
    integer :: i1, l1, j1, e1
    integer :: i2, l2, j2, e2
    integer :: i3, l3, j3, e3
    integer :: i4, l4, j4, e4
    integer :: i5, l5, j5, e5, i5max
    integer :: i6, l6, j6, e6, i6max
    integer :: J12, T12, J45, T45, J, T, P123, P456

    r = 0
    do i1 = 1, spsf%norbs
      l1 = spsf%orb(i1)%l
      j1 = spsf%orb(i1)%j
      e1 = spsf%orb(i1)%e
      do i2 = 1, i1
        l2 = spsf%orb(i2)%l
        j2 = spsf%orb(i2)%j
        e2 = spsf%orb(i2)%e
        if(e1 + e2 > e2max) cycle
        do i3 = 1, i2
          l3 = spsf%orb(i3)%l
          j3 = spsf%orb(i3)%j
          e3 = spsf%orb(i3)%e
          if(e1 + e3 > e2max) cycle
          if(e2 + e3 > e2max) cycle
          if(e1 + e2 + e3 > e3max) cycle

          P123 = (-1) ** (l1+l2+l3)

          do i4 = 1, i1
            l4 = spsf%orb(i4)%l
            j4 = spsf%orb(i4)%j
            e4 = spsf%orb(i4)%e

            i5max = i4
            if(i1 == i4) i5max = i2

            do i5 = 1, i5max
              l5 = spsf%orb(i5)%l
              j5 = spsf%orb(i5)%j
              e5 = spsf%orb(i5)%e
              if(e4 + e5 > e2max) cycle

              i6max = i5
              if(i1 == i4 .and. i2 == i5) i6max = i3

              do i6 = 1, i6max
                l6 = spsf%orb(i6)%l
                j6 = spsf%orb(i6)%j
                e6 = spsf%orb(i6)%e
                if(e4 + e6 > e2max) cycle
                if(e5 + e6 > e2max) cycle
                if(e4 + e5 + e6 > e3max) cycle

                P456 = (-1) ** (l4+l5+l6)

                if(P123 /= P456) cycle
                do J12 = abs(j1-j2)/2, (j1+j2)/2
                  do J45 = abs(j4-j5)/2, (j4+j5)/2
                    do J = max(abs(2*J12-j3),abs(2*J45-j6)),&
                          &min(   (2*J12+j3),   (2*J45+j6)), 2

                      do T12 = 0, 1
                        do T45 = 0, 1
                          do T = max(abs(2*T12-1),abs(2*T45-1)),&
                                &min(   (2*T12+1),   (2*T45+1)), 2

                            r = r + 1

                          end do

                        end do
                      end do
                    end do
                  end do
                end do


              end do
            end do
          end do


        end do
      end do
    end do
  end function count_scalar_3bme

  subroutine store_scalar_3bme(thr,sps,ms,v,spsf,e2max,e3max)
    type(NBodyPart), intent(inout) :: thr
    type(NonOrthIsospinThreeBodySpace), intent(in) :: ms
    type(OrbitsIsospin), intent(in) :: spsf,sps
    real(8), intent(in) :: v(:)
    integer, intent(in) :: e2max, e3max
    integer(8) :: cnt
    integer :: i1, l1, j1, e1
    integer :: i2, l2, j2, e2
    integer :: i3, l3, j3, e3
    integer :: i4, l4, j4, e4
    integer :: i5, l5, j5, e5, i5max
    integer :: i6, l6, j6, e6, i6max
    integer :: J12, T12, J45, T45, J, T, P123, P456
    integer :: ch, idxb, idxk, bra, ket

    cnt = 0
    do i1 = 1, spsf%norbs
      l1 = spsf%orb(i1)%l
      j1 = spsf%orb(i1)%j
      e1 = spsf%orb(i1)%e
      do i2 = 1, i1
        l2 = spsf%orb(i2)%l
        j2 = spsf%orb(i2)%j
        e2 = spsf%orb(i2)%e
        if(e1 + e2 > e2max) cycle
        do i3 = 1, i2
          l3 = spsf%orb(i3)%l
          j3 = spsf%orb(i3)%j
          e3 = spsf%orb(i3)%e
          if(e1 + e3 > e2max) cycle
          if(e2 + e3 > e2max) cycle
          if(e1 + e2 + e3 > e3max) cycle

          P123 = (-1) ** (l1+l2+l3)

          do i4 = 1, i1
            l4 = spsf%orb(i4)%l
            j4 = spsf%orb(i4)%j
            e4 = spsf%orb(i4)%e

            i5max = i4
            if(i1 == i4) i5max = i2

            do i5 = 1, i5max
              l5 = spsf%orb(i5)%l
              j5 = spsf%orb(i5)%j
              e5 = spsf%orb(i5)%e
              if(e4 + e5 > e2max) cycle

              i6max = i5
              if(i1 == i4 .and. i2 == i5) i6max = i3

              do i6 = 1, i6max
                l6 = spsf%orb(i6)%l
                j6 = spsf%orb(i6)%j
                e6 = spsf%orb(i6)%e
                if(e4 + e6 > e2max) cycle
                if(e5 + e6 > e2max) cycle
                if(e4 + e5 + e6 > e3max) cycle

                P456 = (-1) ** (l4+l5+l6)

                if(P123 /= P456) cycle
                do J12 = abs(j1-j2)/2, (j1+j2)/2
                  do J45 = abs(j4-j5)/2, (j4+j5)/2
                    do J = max(abs(2*J12-j3),abs(2*J45-j6)),&
                          &min(   (2*J12+j3),   (2*J45+j6)), 2

                      do T12 = 0, 1
                        do T45 = 0, 1
                          do T = max(abs(2*T12-1),abs(2*T45-1)),&
                                &min(   (2*T12+1),   (2*T45+1)), 2
                            cnt = cnt + 1

                            if(e1 > ms%emax) cycle
                            if(e2 > ms%emax) cycle
                            if(e3 > ms%emax) cycle

                            if(e4 > ms%emax) cycle
                            if(e5 > ms%emax) cycle
                            if(e6 > ms%emax) cycle

                            if(e1 + e2 > ms%e2max) cycle
                            if(e2 + e3 > ms%e2max) cycle
                            if(e3 + e1 > ms%e2max) cycle

                            if(e4 + e5 > ms%e2max) cycle
                            if(e5 + e6 > ms%e2max) cycle
                            if(e6 + e4 > ms%e2max) cycle

                            if(e1 + e2 + e3 > ms%e3max) cycle
                            if(e4 + e5 + e6 > ms%e3max) cycle

                            if(i1==i2 .and. mod(J12+T12,2)==0) then
                              if(abs(v(cnt)) > 1.d-6) then
                                write(*,*) "Warning: something wrong, this three-body matrix element has to be zero."
                              end if
                              cycle
                            end if

                            if(i4==i5 .and. mod(J45+T45,2)==0) then
                              if(abs(v(cnt)) > 1.d-6) then
                                write(*,*) "Warning: something wrong, this three-body matrix element has to be zero."
                              end if
                              cycle
                            end if

                            ch = ms%jpt2ch(J,P123,T)
                            idxb = ms%jpt(ch)%spis2idx(i1,i2,i3)
                            idxk = ms%jpt(ch)%spis2idx(i4,i5,i6)
                            if(idxb*idxk == 0) cycle
                            bra = ms%jpt(ch)%idxqn(idxb)%JT2n(J12,T12)
                            ket = ms%jpt(ch)%idxqn(idxk)%JT2n(J45,T45)
                            if(bra*ket == 0) cycle
                            thr%MatCh(ch,ch)%m(bra,ket) = v(cnt)
                          end do

                        end do
                      end do
                    end do
                  end do
                end do


              end do
            end do
          end do


        end do
      end do
    end do
  end subroutine store_scalar_3bme

  subroutine store_scalar_3bme_sp(thr,sps,ms,v,spsf,e2max,e3max)
    type(NBodyPartSp), intent(inout) :: thr
    type(NonOrthIsospinThreeBodySpace), intent(in) :: ms
    type(OrbitsIsospin), intent(in) :: spsf,sps
    real(4), intent(in) :: v(:)
    integer, intent(in) :: e2max, e3max
    integer(8) :: cnt
    integer :: i1, l1, j1, e1
    integer :: i2, l2, j2, e2
    integer :: i3, l3, j3, e3
    integer :: i4, l4, j4, e4
    integer :: i5, l5, j5, e5, i5max
    integer :: i6, l6, j6, e6, i6max
    integer :: J12, T12, J45, T45, J, T, P123, P456
    integer :: ch, idxb, idxk, bra, ket

    cnt = 0
    do i1 = 1, spsf%norbs
      l1 = spsf%orb(i1)%l
      j1 = spsf%orb(i1)%j
      e1 = spsf%orb(i1)%e
      do i2 = 1, i1
        l2 = spsf%orb(i2)%l
        j2 = spsf%orb(i2)%j
        e2 = spsf%orb(i2)%e
        if(e1 + e2 > e2max) cycle
        do i3 = 1, i2
          l3 = spsf%orb(i3)%l
          j3 = spsf%orb(i3)%j
          e3 = spsf%orb(i3)%e
          if(e1 + e3 > e2max) cycle
          if(e2 + e3 > e2max) cycle
          if(e1 + e2 + e3 > e3max) cycle

          P123 = (-1) ** (l1+l2+l3)

          do i4 = 1, i1
            l4 = spsf%orb(i4)%l
            j4 = spsf%orb(i4)%j
            e4 = spsf%orb(i4)%e

            i5max = i4
            if(i1 == i4) i5max = i2

            do i5 = 1, i5max
              l5 = spsf%orb(i5)%l
              j5 = spsf%orb(i5)%j
              e5 = spsf%orb(i5)%e
              if(e4 + e5 > e2max) cycle

              i6max = i5
              if(i1 == i4 .and. i2 == i5) i6max = i3

              do i6 = 1, i6max
                l6 = spsf%orb(i6)%l
                j6 = spsf%orb(i6)%j
                e6 = spsf%orb(i6)%e
                if(e4 + e6 > e2max) cycle
                if(e5 + e6 > e2max) cycle
                if(e4 + e5 + e6 > e3max) cycle

                P456 = (-1) ** (l4+l5+l6)

                if(P123 /= P456) cycle
                do J12 = abs(j1-j2)/2, (j1+j2)/2
                  do J45 = abs(j4-j5)/2, (j4+j5)/2
                    do J = max(abs(2*J12-j3),abs(2*J45-j6)),&
                          &min(   (2*J12+j3),   (2*J45+j6)), 2

                      do T12 = 0, 1
                        do T45 = 0, 1
                          do T = max(abs(2*T12-1),abs(2*T45-1)),&
                                &min(   (2*T12+1),   (2*T45+1)), 2
                            cnt = cnt + 1

                            if(e1 > ms%emax) cycle
                            if(e2 > ms%emax) cycle
                            if(e3 > ms%emax) cycle

                            if(e4 > ms%emax) cycle
                            if(e5 > ms%emax) cycle
                            if(e6 > ms%emax) cycle

                            if(e1 + e2 > ms%e2max) cycle
                            if(e2 + e3 > ms%e2max) cycle
                            if(e3 + e1 > ms%e2max) cycle

                            if(e4 + e5 > ms%e2max) cycle
                            if(e5 + e6 > ms%e2max) cycle
                            if(e6 + e4 > ms%e2max) cycle

                            if(e1 + e2 + e3 > ms%e3max) cycle
                            if(e4 + e5 + e6 > ms%e3max) cycle

                            if(i1==i2 .and. mod(J12+T12,2)==0) then
                              if(abs(v(cnt)) > 1.d-6) then
                                write(*,*) "Warning: something wrong, this three-body matrix element has to be zero."
                              end if
                              cycle
                            end if

                            if(i4==i5 .and. mod(J45+T45,2)==0) then
                              if(abs(v(cnt)) > 1.d-6) then
                                write(*,*) "Warning: something wrong, this three-body matrix element has to be zero."
                              end if
                              cycle
                            end if
                            ch = ms%jpt2ch(J,P123,T)
                            idxb = ms%jpt(ch)%spis2idx(i1,i2,i3)
                            idxk = ms%jpt(ch)%spis2idx(i4,i5,i6)
                            if(idxb*idxk == 0) cycle
                            bra = ms%jpt(ch)%idxqn(idxb)%JT2n(J12,T12)
                            ket = ms%jpt(ch)%idxqn(idxk)%JT2n(J45,T45)
                            if(bra*ket == 0) cycle
                            thr%MatCh(ch,ch)%m(bra,ket) = v(cnt)
                          end do

                        end do
                      end do
                    end do
                  end do
                end do


              end do
            end do
          end do


        end do
      end do
    end do
  end subroutine store_scalar_3bme_sp

  !
  !
  ! reading three-body tensor
  !
  !

  subroutine ReadTensor3BFile(this,thr,sps,ms)
    use ClassSys, only: sys
    class(ReadFiles), intent(in) :: this
    type(NBodyPart), intent(inout) :: thr
    type(OrbitsIsospin), intent(in) :: sps
    type(NonOrthIsospinThreeBodySpace), intent(in) :: ms
    type(sys) :: s

    if(s%find(this%file_3n,'.txt')) then
      call this%read_tensor_3bme_txt(thr,sps,ms)
      return
    end if

    if(s%find(this%file_3n,'.bin')) then
      call this%read_tensor_3bme_bin(thr,sps,ms)
      return
    end if
  end subroutine ReadTensor3BFile

  subroutine read_tensor_3bme_txt(this,thr,sps,ms)
    class(ReadFiles), intent(in) :: this
    type(NBodyPart), intent(inout) :: thr
    type(OrbitsIsospin), intent(in) :: sps
    type(NonOrthIsospinThreeBodySpace), intent(in) :: ms
    write(*,*) "reading tensor is not ready"
    return
  end subroutine read_tensor_3bme_txt

  subroutine read_tensor_3bme_bin(this,thr,sps,ms)
    class(ReadFiles), intent(in) :: this
    type(NBodyPart), intent(inout) :: thr
    type(OrbitsIsospin), intent(in) :: sps
    type(NonOrthIsospinThreeBodySpace), intent(in) :: ms
    write(*,*) "reading tensor is not ready"
    return
  end subroutine read_tensor_3bme_bin

  subroutine ReadTensor3BFileSp(this,thr,sps,ms)
    use ClassSys, only: sys
    class(ReadFiles), intent(in) :: this
    type(NBodyPartSp), intent(inout) :: thr
    type(OrbitsIsospin), intent(in) :: sps
    type(NonOrthIsospinThreeBodySpace), intent(in) :: ms
    type(sys) :: s
    if(s%find(this%file_3n,'.txt')) then
      call this%read_tensor_3bme_txt_sp(thr,sps,ms)
      return
    end if

    if(s%find(this%file_3n,'.bin')) then
      call this%read_tensor_3bme_bin_sp(thr,sps,ms)
      return
    end if
  end subroutine ReadTensor3BFileSp

  subroutine read_tensor_3bme_txt_sp(this,thr,sps,ms)
    class(ReadFiles), intent(in) :: this
    type(NBodyPartSp), intent(inout) :: thr
    type(OrbitsIsospin), intent(in) :: sps
    type(NonOrthIsospinThreeBodySpace), intent(in) :: ms
    write(*,*) "reading tensor is not ready"
    return
  end subroutine read_tensor_3bme_txt_sp

  subroutine read_tensor_3bme_bin_sp(this,thr,sps,ms)
    class(ReadFiles), intent(in) :: this
    type(NBodyPartSp), intent(inout) :: thr
    type(OrbitsIsospin), intent(in) :: sps
    type(NonOrthIsospinThreeBodySpace), intent(in) :: ms
    write(*,*) "reading tensor is not ready"
    return
  end subroutine read_tensor_3bme_bin_sp

end module NOperators

! test for operators
!program test
!  use Profiler, only: timer
!  use CommonLibrary, only: &
!      &init_dbinomial_triangle, fin_dbinomial_triangle
!  use ModelSpace, only: MSpace
!  use NOperators
!
!  implicit none
!
!  type(MSpace) :: ms
!  type(NBodyPart) :: vnn_me2j, vnn_myg, vnn_diff
!  type(NBodyPart) :: vnn_myg_bin
!  type(ReadFiles) :: rd
!
!  call timer%init()
!  call init_dbinomial_triangle()
!
!  call ms%init('O16', 25.d0, 2, 4, e3max=2)
!
!  !rd%file_nn = '/home/takayuki/TwBME-HO_NN-only_N3LO_EM500_bare_hw25_emax6_e2max12.txt.me2j'
!  rd%file_nn = '/home/takayuki/TwBME-HO_NN-only_N3LO_EM500_bare_hw25_emax2_e2max4.txt.me2j'
!  rd%emax2 = 2
!  rd%e2max2 =4
!  rd%lmax2 = 2
!  call vnn_me2j%init(ms%two,.true.,'hamil',0,1,0)
!  call vnn_me2j%ReadFile(ms%sps, ms%two, rd)
!
!  rd%file_nn = '/home/takayuki/TwBME-HO_NN-only_N3LO_EM500_bare_hw25_emax2_e2max4.txt.myg'
!  rd%emax2 = 2
!  rd%e2max2 =4
!  rd%lmax2 = 2
!  call vnn_myg%init(ms%two,.true.,'hamil',0,1,0)
!  call vnn_myg%ReadFile(ms%sps, ms%two, rd)
!
!  rd%file_nn = '/home/takayuki/TwBME-HO_NN-only_N3LO_EM500_bare_hw25_emax2_e2max4.txt.snt'
!  rd%emax2 = 2
!  rd%e2max2 =4
!  rd%lmax2 = 2
!  call vnn_myg_bin%init(ms%two,.true.,'hamil',0,1,0)
!  call vnn_myg_bin%ReadFile(ms%sps, ms%two, rd)
!
!  vnn_diff = vnn_myg - vnn_myg_bin
!
!  call vnn_diff%prt()
!
!  call vnn_me2j%fin()
!  call vnn_myg%fin()
!  call vnn_diff%fin()
!  call ms%fin()
!
!  call fin_dbinomial_triangle()
!  call timer%fin()
!end program test
