module TwoBodyOperator
  use omp_lib
  use LinAlgLib
  use TwoBodyModelSpace
  use OneBodyOperator
  implicit none

  public :: TwoBodyPartChannel
  public :: TwoBodyPart

  private :: InitTwoBodyPartChannel
  private :: FinTwoBodyPartChannel
  private :: SetTwoBodyPartChannel

  private :: InitTwoBodyPart
  private :: FinTwoBodyPart
  private :: SetTwoBodyPart
  private :: PrintTwoBodyPart
  private :: SetTwBME_scalar
  private :: SetTwBME_tensor
  private :: GetTwBME_scalar
  private :: GetTwBME_tensor
  private :: AddToTwBME_scalar
  private :: AddToTwBME_tensor
  private :: CopyTwoBodyPart
  private :: SumTwoBodyPart
  private :: SubtractTwoBodyPart
  private :: ScaleTwoBodyPart
  private :: NormalOrderingFrom2To1

  private :: ExtractTwoBodyChannel_pppp
  private :: ExtractTwoBodyChannel_hhhh
  private :: ExtractTwoBodyChannel_pphh
  private :: ExtractTwoBodyChannel_phhh
  private :: ExtractTwoBodyChannel_phpp
  private :: SetTwoBodyChannel_pppp
  private :: SetTwoBodyChannel_hhhh
  private :: SetTwoBodyChannel_pphh
  private :: SetTwoBodyChannel_phhh
  private :: SetTwoBodyChannel_phpp
  private :: ExtractCrossCoupledChannel_pphh2phph
  private :: ExtractCrossCoupledChannel_hphp2phph
  private :: ExtractCrossCoupledChannel_phph2phph
  !private :: ExtractCrossCoupledChannel_phhh2phhh
  !private :: ExtractCrossCoupledChannel_ppph2ppph
  private :: GetCrossCoupledME1423
  private :: GetCrossCoupledME1324

  ! methods for two-body matrix element
  private :: ReadTwoBodyFile
  private :: Set2BodyReadFile
  private :: Set2BodyFileBoundaries
  private :: ReadScalar2BFile
  private :: ReadTensor2BFile
  private :: read_scalar_me2j_ascii
  private :: read_scalar_me2j_bin
  private :: read_scalar_myg_ascii
  private :: read_scalar_myg_bin
  private :: read_scalar_snt_ascii
  private :: read_scalar_navratil_ascii
  private :: read_scalar_navratil_ascii_gz

  type, extends(DMat) :: TwoBodyPartChannel
    type(TwoBodyChannel), pointer :: ch_bra, ch_ket
    logical :: is = .false.
  contains
    procedure :: InitTwoBodyPartChannel
    procedure :: FinTwoBodyPartChannel
    procedure :: SetTwoBodyPartChannel
    procedure :: SetTwoBodyPartChannelFromOneBody
    procedure :: ExtractTwoBodyChannel_pppp
    procedure :: ExtractTwoBodyChannel_hhhh
    procedure :: ExtractTwoBodyChannel_pphh
    procedure :: ExtractTwoBodyChannel_phhh
    procedure :: ExtractTwoBodyChannel_phpp
    procedure :: SetTwoBodyChannel_pppp
    procedure :: SetTwoBodyChannel_hhhh
    procedure :: SetTwoBodyChannel_pphh
    procedure :: SetTwoBodyChannel_phhh
    procedure :: SetTwoBodyChannel_phpp

    generic :: init => InitTwoBodyPartChannel
    generic :: release => FinTwoBodyPartChannel
    generic :: set => SetTwoBodyPartChannel, SetTwoBodyPartChannelFromOneBody
    generic :: get_pppp => ExtractTwoBodyChannel_pppp
    generic :: get_hhhh => ExtractTwoBodyChannel_hhhh
    generic :: get_pphh => ExtractTwoBodyChannel_pphh
    generic :: get_phhh => ExtractTwoBodyChannel_phhh
    generic :: get_phpp => ExtractTwoBodyChannel_phpp
    generic :: set_pppp => SetTwoBodyChannel_pppp
    generic :: set_hhhh => SetTwoBodyChannel_hhhh
    generic :: set_pphh => SetTwoBodyChannel_pphh
    generic :: set_phhh => SetTwoBodyChannel_phhh
    generic :: set_phpp => SetTwoBodyChannel_phpp
  end type TwoBodyPartChannel

  type :: TwoBodyPart
    type(TwoBodyPartChannel), allocatable :: MatCh(:,:)
    type(TwoBodySpace), pointer :: two
    character(:), allocatable :: oprtr
    logical :: Scalar
    integer :: jr, pr, tr, zr
  contains
    procedure :: InitTwoBodyPart
    procedure :: FinTwoBodyPart
    procedure :: SetTwoBodyPart
    procedure :: SetTwoBodyPartFromOneBody
    procedure :: SetTwBME_scalar
    procedure :: SetTwBME_tensor
    procedure :: GetTwBME_scalar
    procedure :: GetTwBME_tensor
    procedure :: GetTwBMEMon
    procedure :: AddToTwBME_scalar
    procedure :: AddToTwBME_tensor
    procedure :: CopyTwoBodyPart
    procedure :: SumTwoBodyPart
    procedure :: SubtractTwoBodyPart
    procedure :: ScaleTwoBodyPart
    procedure :: PrintTwoBodyPart
    procedure :: NormalOrderingFrom2To1
    procedure :: ExtractCrossCoupledChannel_pphh2phph
    procedure :: ExtractCrossCoupledChannel_hphp2phph
    procedure :: ExtractCrossCoupledChannel_phph2phph
    !procedure :: ExtractCrossCoupledChannel_phhh2phhh
    !procedure :: ExtractCrossCoupledChannel_ppph2ppph
    procedure :: GetCrossCoupledME1423
    procedure :: GetCrossCoupledME1324

    generic :: init => InitTwoBodyPart
    generic :: fin => FinTwoBodyPart
    generic :: set => SetTwoBodyPart, SetTwoBodyPartFromOneBody
    generic :: prt => PrintTwoBodyPart
    generic :: SetTwBME => SetTwBME_scalar, SetTwBME_tensor
    generic :: GetTwBME => GetTwBME_scalar, GetTwBME_tensor
    generic :: AddToTwBME => AddToTwBME_scalar, AddToTwBME_tensor
    generic :: assignment(=) => CopyTwoBodyPart
    generic :: operator(+) => SumTwoBodyPart
    generic :: operator(-) => SubtractTwoBodyPart
    generic :: operator(*) => ScaleTwoBodyPart
    generic :: get_xc_pphh2phph => ExtractCrossCoupledChannel_pphh2phph
    generic :: get_xc_hphp2phph => ExtractCrossCoupledChannel_hphp2phph
    generic :: get_xc_phph2phph => ExtractCrossCoupledChannel_phph2phph
    !generic :: get_xc_phhh2phhh => ExtractCrossCoupledChannel_phhh2phhh
    !generic :: get_xc_ppph2ppph => ExtractCrossCoupledChannel_ppph2ppph
    generic :: get_xc1423 => GetCrossCoupledME1423
    generic :: get_xc1324 => GetCrossCoupledME1324
  end type TwoBodyPart

  type :: Read2BodyFiles
    character(:), allocatable :: file_nn
    integer :: emax2=-1, e2max2=-1, lmax2=-1
  contains
    procedure :: Set2BodyReadFile       ! setter
    procedure :: Set2BodyFileBoundaries ! setter
    generic :: set => Set2BodyReadFile, Set2BodyFileBoundaries
    procedure :: ReadTwoBodyFile
    ! methods for two-body matrix element
    procedure :: ReadScalar2BFile
    procedure :: ReadTensor2BFile
    procedure :: read_scalar_me2j_ascii
    procedure :: read_scalar_me2j_bin
    procedure :: read_scalar_myg_ascii
    procedure :: read_scalar_myg_bin
    procedure :: read_scalar_snt_ascii
    procedure :: read_scalar_navratil_ascii
    procedure :: read_scalar_navratil_ascii_gz
  end type Read2BodyFiles

contains

  subroutine FinTwoBodyPart(this)
    use MyLibrary, only: triag
    class(TwoBodyPart), intent(inout) :: this
    integer :: chbra, chket

    if(.not. allocated(this%MatCh)) return
    do chbra = 1, this%Two%NChan
      do chket = 1, this%Two%NChan
        call this%MatCh(chbra,chket)%release()
      end do
    end do
    deallocate(this%MatCh)
    this%two => null()

  end subroutine FinTwoBodyPart

  subroutine InitTwoBodyPart(this, two, Scalar, oprtr, jr, pr, zr)
    use MyLibrary, only: triag
    class(TwoBodyPart), intent(inout) :: this
    type(TwoBodySpace), target, intent(in) :: two
    logical, intent(in) :: Scalar
    character(*), intent(in) :: oprtr
    integer, intent(in) :: jr, pr, zr
    integer :: chbra, chket
    integer :: jbra, pbra, zbra, nbra
    integer :: jket, pket, zket, nket

    if(allocated(this%MatCh)) call this%fin()
    this%two => two
    this%oprtr = oprtr
    this%Scalar = Scalar
    this%jr = jr
    this%pr = pr
    this%zr = zr

    allocate(this%MatCh(Two%NChan, Two%NChan))
    do chbra = 1, Two%NChan
      jbra = two%jpz(chbra)%j
      pbra = two%jpz(chbra)%p
      zbra = two%jpz(chbra)%z
      nbra = two%jpz(chbra)%n_state
      do chket = 1, Two%NChan
        jket = two%jpz(chket)%j
        pket = two%jpz(chket)%p
        zket = two%jpz(chket)%z
        nket = two%jpz(chket)%n_state

        if(triag(jbra, jket, 2*jr)) cycle
        if(pbra * pket * pr /= 1) cycle
        if(zbra - 2*zr - zket /= 0) cycle
        call this%MatCh(chbra,chket)%init(two%jpz(chbra),two%jpz(chket))
      end do
    end do

  end subroutine InitTwoBodyPart

  subroutine SetTwoBodyPart(this, hw, A, Z, N)
    class(TwoBodyPart), intent(inout) :: this
    real(8), intent(in) :: hw
    integer, intent(in) :: A, Z, N
    type(TwoBodySpace), pointer :: tbs
    type(Orbits), pointer :: sps
    integer :: chbra, chket

    tbs => this%two
    sps => tbs%sps
    do chbra = 1, tbs%NChan
      do chket = 1, tbs%NChan
        if(.not. this%MatCh(chbra,chket)%is) cycle
        call this%MatCh(chbra,chket)%set(this%oprtr, sps, hw, A, Z, N)
      end do
    end do
  end subroutine SetTwoBodyPart

  subroutine CopyTwoBodyPart(a, b)
    class(TwoBodyPart), intent(out) :: a
    type(TwoBodyPart), intent(in) :: b
    integer :: chbra, chket

    if(allocated(a%MatCh)) call a%fin()
    a%two => b%two
    a%oprtr = b%oprtr
    a%Scalar = b%Scalar
    a%jr = b%jr
    a%pr = b%pr
    a%zr = b%zr

    allocate(a%MatCh(b%two%NChan,b%two%NChan))
    do chbra = 1, b%two%NChan
      do chket = 1, b%two%NChan
        if(.not. b%MatCh(chbra,chket)%is) cycle
        a%MatCh(chbra,chket)%is = b%MatCh(chbra,chket)%is
        a%MatCh(chbra,chket)%DMat = b%MatCh(chbra,chket)%DMat
        a%MatCh(chbra,chket)%ch_bra => b%MatCh(chbra,chket)%ch_bra
        a%MatCh(chbra,chket)%ch_ket => b%MatCh(chbra,chket)%ch_ket
      end do
    end do
  end subroutine CopyTwoBodyPart

  function SumTwoBodyPart(a, b) result(c)
    class(TwoBodyPart), intent(in) :: a, b
    type(TwoBodyPart) :: c
    integer :: chbra, chket
    c = a
    do chbra = 1, b%two%NChan
      do chket = 1, b%two%NChan
        if(.not. b%MatCh(chbra,chket)%is) cycle
        c%MatCh(chbra,chket)%DMat = &
            & a%MatCh(chbra,chket)%DMat + b%MatCh(chbra,chket)%DMat
      end do
    end do
  end function SumTwoBodyPart

  function SubtractTwoBodyPart(a, b) result(c)
    class(TwoBodyPart), intent(in) :: a, b
    type(TwoBodyPart) :: c
    integer :: chbra, chket
    c = a
    do chbra = 1, b%two%NChan
      do chket = 1, b%two%NChan
        if(.not. b%MatCh(chbra,chket)%is) cycle
        c%MatCh(chbra,chket)%DMat = &
            & a%MatCh(chbra,chket)%DMat - b%MatCh(chbra,chket)%DMat
      end do
    end do
  end function SubtractTwoBodyPart

  function ScaleTwoBodyPart(a, b) result(c)
    class(TwoBodyPart), intent(in) :: a
    real(8), intent(in) :: b
    type(TwoBodyPart) :: c
    integer :: chbra, chket
    c = a
    do chbra = 1, a%two%NChan
      do chket = 1, a%two%NChan
        if(.not. a%MatCh(chbra,chket)%is) cycle
        c%MatCh(chbra,chket)%DMat = &
            & a%MatCh(chbra,chket)%DMat * b
      end do
    end do
  end function ScaleTwoBodyPart

  subroutine FinTwoBodyPartChannel(this)
    class(TwoBodyPartChannel), intent(inout) :: this
    if(.not. this%is) return
    call this%fin()
    this%is = .false.
    this%ch_bra => null()
    this%ch_ket => null()
  end subroutine FinTwoBodyPartChannel

  subroutine InitTwoBodyPartChannel(this, ch_bra, ch_ket)
    class(TwoBodyPartChannel), intent(inout) :: this
    type(TwoBodyChannel), target, intent(in) :: ch_bra, ch_ket
    this%ch_bra => ch_bra
    this%ch_ket => ch_ket
    call this%zeros(ch_bra%n_state, ch_ket%n_state)
    this%is = .true.
  end subroutine InitTwoBodyPartChannel

  subroutine SetTwoBodyPartChannel(this, optr, sps, hw, A, Z, N)
    class(TwoBodyPartChannel), intent(inout) :: this
    character(*), intent(in) :: optr
    type(Orbits), pointer, intent(in) :: sps
    real(8), intent(in) :: hw
    integer, intent(in) :: A, Z, N
    type(TwoBodyChannel), pointer :: chbra, chket
    integer :: bra, ket, ia, ib, ic, id

    chbra => this%ch_bra
    chket => this%ch_ket

    !$omp parallel
    !$omp do private(bra,ia,ib,ket,ic,id)
    do bra = 1, chbra%n_state
      ia = chbra%n2spi1(bra)
      ib = chbra%n2spi2(bra)
      do ket = 1, chket%n_state
        ic = chket%n2spi1(ket)
        id = chket%n2spi2(ket)

        this%m(bra,ket) = mat_elm(sps,optr,hw,A,Z,N,&
            & ia,ib,ic,id,chbra%J,chket%J)
      end do
    end do
    !$omp end do
    !$omp end parallel
  end subroutine SetTwoBodyPartChannel

  function mat_elm(sps,optr,hw,A,Z,N,ia,ib,ic,id,Jbra,Jket) result(r)
    use DefineOperators, only: two_body_element
    type(Orbits), intent(in) :: sps
    character(*), intent(in) :: optr
    real(8), intent(in) :: hw
    integer, intent(in) :: A, Z, N
    integer, intent(in) :: ia, ib, ic, id,Jbra,Jket
    real(8) :: r
    type(SingleParticleOrbit), pointer :: oa, ob, oc, od
    integer :: aa(4), bb(4), cc(4), dd(4)

    r = 0.d0
    oa => sps%GetOrbit(ia)
    ob => sps%GetOrbit(ib)
    oc => sps%GetOrbit(ic)
    od => sps%GetOrbit(id)
    aa = [oa%n, oa%l, oa%j, oa%z]
    bb = [ob%n, ob%l, ob%j, ob%z]
    cc = [oc%n, oc%l, oc%j, oc%z]
    dd = [od%n, od%l, od%j, od%z]
    r = two_body_element(optr, aa, bb, cc, dd, Jbra, Jket, hw, A, Z, N)
  end function mat_elm

  subroutine PrintTwoBodyPart(this, wunit)
    use ClassSys, only: sys
    class(TwoBodyPart), intent(in) :: this
    integer, intent(in), optional :: wunit
    integer :: chbra, chket
    type(sys) :: s
    character(:), allocatable :: msg
    integer :: jbra, pbra, zbra, jket, pket, zket

    msg = ""
    do chbra = 1, this%two%NChan
      jbra = this%two%jpz(chbra)%j
      pbra = this%two%jpz(chbra)%p
      zbra = this%two%jpz(chbra)%z
      do chket = 1, this%two%NChan
        jket = this%two%jpz(chket)%j
        pket = this%two%jpz(chket)%p
        zket = this%two%jpz(chket)%z
        if(.not. this%MatCh(chbra,chket)%is) cycle
        msg = trim(this%oprtr) // " (" // trim(s%str(jbra)) // &
            & "," // trim(s%str(pbra)) // "," // trim(s%str(zbra)) // &
            & ")  (" // trim(s%str(jket)) // "," // &
            & trim(s%str(pket)) // "," // trim(s%str(zket)) // ")"
        call this%MatCh(chbra,chket)%prt(msg=msg,iunit=wunit)
      end do
    end do
  end subroutine PrintTwoBodyPart

  function GetTwBME_tensor(this,i1,i2,i3,i4,J12,J34) result(r)
    use MyLibrary, only: triag
    real(8) :: r
    class(TwoBodyPart), intent(in) :: this
    integer, intent(in) :: i1,i2,i3,i4,J12,J34
    type(TwoBodySpace), pointer :: tbs
    type(SingleParticleOrbit), pointer :: o1, o2, o3, o4
    integer :: P12, P34
    integer :: Z12, Z34
    integer :: chbra, chket, bra, ket, iphase

    r = 0.d0
    tbs => this%two
    o1 => tbs%sps%orb(i1)
    o2 => tbs%sps%orb(i2)
    o3 => tbs%sps%orb(i3)
    o4 => tbs%sps%orb(i4)
    P12 = (-1) ** (o1%l + o2%l)
    P34 = (-1) ** (o3%l + o4%l)
    Z12 = (o1%z + o2%z)/2
    Z34 = (o3%z + o4%z)/2

    if(triag(J12, J34, this%jr)) then
      write(*,*) "Warning: in GetTwBME_tensor: J"
      return
    end if

    if(P12 * P34 * this%pr /= 1) then
      write(*,*) "Warning: in GetTwBME_tensor: P"
      return
    end if

    if(Z12 - Z34 - this%zr /= 0) then
      write(*,*) "Warning: in GetTwBME_tensor: Tz"
      return
    end if

    chbra = tbs%jpz2ch(J12,P12,Z12)
    chket = tbs%jpz2ch(J34,P34,Z34)
    if(chbra * chket == 0) return

    bra = tbs%jpz(chbra)%spis2n(i1,i2)
    ket = tbs%jpz(chket)%spis2n(i3,i4)
    if(bra * ket == 0) return

    iphase = tbs%jpz(chbra)%iphase(i1,i2) * &
        &    tbs%jpz(chket)%iphase(i3,i4)

    r = dble(iphase) * this%MatCh(chbra,chket)%m(bra,ket)
  end function GetTwBME_tensor

  function GetTwBME_scalar(this,i1,i2,i3,i4,J) result(r)
    use MyLibrary, only: triag
    real(8) :: r
    class(TwoBodyPart), intent(in) :: this
    integer, intent(in) :: i1,i2,i3,i4,J
    type(TwoBodySpace), pointer :: tbs
    type(SingleParticleOrbit), pointer :: o1, o2, o3, o4
    integer :: P12, P34
    integer :: Z12, Z34
    integer :: ch, bra, ket, iphase

    r = 0.d0
    tbs => this%two
    o1 => tbs%sps%orb(i1)
    o2 => tbs%sps%orb(i2)
    o3 => tbs%sps%orb(i3)
    o4 => tbs%sps%orb(i4)
    P12 = (-1) ** (o1%l + o2%l)
    P34 = (-1) ** (o3%l + o4%l)
    Z12 = (o1%z + o2%z)/2
    Z34 = (o3%z + o4%z)/2

    if(P12 * P34 /= 1) then
      write(*,*) "Warning: in GetTwBME_scalar: P"
      return
    end if

    if(Z12 - Z34 /= 0) then
      write(*,*) "Warning: in GetTwBME_scalar: Tz"
      return
    end if

    ch = tbs%jpz2ch(J,P12,Z12)
    if(ch == 0) return

    bra = tbs%jpz(ch)%spis2n(i1,i2)
    ket = tbs%jpz(ch)%spis2n(i3,i4)
    if(bra * ket == 0) return

    iphase = tbs%jpz(ch)%iphase(i1,i2) * &
        &    tbs%jpz(ch)%iphase(i3,i4)

    r = dble(iphase) * this%MatCh(ch,ch)%m(bra,ket)
  end function GetTwBME_scalar

  function GetTwBMEMon(this, i1, i2, i3, i4) result(r)
    real(8) :: r
    class(TwoBodyPart), intent(in) :: this
    integer, intent(in) :: i1,i2,i3,i4
    integer :: J, Jmin, Jmax
    type(TwoBodySpace), pointer :: tbs
    type(SingleParticleOrbit), pointer :: o1, o2, o3, o4
    integer :: P12, P34
    integer :: Z12, Z34
    real(8) :: norm, sumJ
    r = 0.d0
    tbs => this%two
    o1 => tbs%sps%GetOrbit(i1)
    o2 => tbs%sps%GetOrbit(i2)
    o3 => tbs%sps%GetOrbit(i3)
    o4 => tbs%sps%GetOrbit(i4)
    P12 = (-1) ** (o1%l + o2%l)
    P34 = (-1) ** (o3%l + o4%l)
    Z12 = (o1%z + o2%z)/2
    Z34 = (o3%z + o4%z)/2

    if(P12 * P34 /= 1) then
      write(*,*) "Warning: in GetTwBMEMon: P"
      return
    end if

    if(Z12 - Z34 /= 0) then
      write(*,*) "Warning: in GetTwBMEMon: Tz"
      return
    end if
    norm = 1.d0
    if(i1==i2) norm = norm * sqrt(2.d0)
    if(i3==i4) norm = norm * sqrt(2.d0)
    Jmin = max(abs(o1%j-o2%j), abs(o3%j-o4%j))/2
    Jmax = min(   (o1%j+o2%j),    (o3%j+o4%j))/2
    sumJ = 0.d0
    do J = Jmin, Jmax
      r = r + dble(2*J+1) * this%GetTwBME(i1, i2, i3, i4, J)
      sumJ = sumJ + dble(2*J+1)
    end do
    !r = r * norm / sqrt(dble( (o1%j+1) * (o2%j+1) * (o3%j+1) * (o4%j+1)) )
    r = r * norm / sumJ
  end function GetTwBMEMon

  subroutine SetTwBME_tensor(this,i1,i2,i3,i4,J12,J34,me)
    use MyLibrary, only: triag
    class(TwoBodyPart), intent(inout) :: this
    integer, intent(in) :: i1,i2,i3,i4,J12,J34
    real(8), intent(in) :: me
    type(TwoBodySpace), pointer :: tbs
    type(SingleParticleOrbit), pointer :: o1, o2, o3, o4
    integer :: P12, P34
    integer :: Z12, Z34
    integer :: chbra, chket, bra, ket, iphase

    tbs => this%two
    o1 => tbs%sps%orb(i1)
    o2 => tbs%sps%orb(i2)
    o3 => tbs%sps%orb(i3)
    o4 => tbs%sps%orb(i4)
    P12 = (-1) ** (o1%l + o2%l)
    P34 = (-1) ** (o3%l + o4%l)
    Z12 = (o1%z + o2%z)/2
    Z34 = (o3%z + o4%z)/2

    if(triag(J12, J34, this%jr)) then
      write(*,*) "Warning: in SetTwBME_general: J"
      return
    end if

    if(P12 * P34 * this%pr /= 1) then
      write(*,*) "Warning: in SetTwBME_general: P"
      return
    end if

    if(Z12 - Z34 - this%zr /= 0) then
      write(*,*) "Warning: in SetTwBME_general: Tz"
      return
    end if

    chbra = tbs%jpz2ch(J12,P12,Z12)
    chket = tbs%jpz2ch(J34,P34,Z34)
    if(chbra * chket == 0) return

    bra = tbs%jpz(chbra)%spis2n(i1,i2)
    ket = tbs%jpz(chket)%spis2n(i3,i4)
    if(bra * ket == 0) return

    iphase = tbs%jpz(chbra)%iphase(i1,i2) * &
        &    tbs%jpz(chket)%iphase(i3,i4)

    this%MatCh(chbra,chket)%m(bra,ket) = dble(iphase) * me
    this%MatCh(chket,chbra)%m(ket,bra) = dble(iphase) * me
  end subroutine SetTwBME_tensor

  subroutine SetTwBME_scalar(this,i1,i2,i3,i4,J,me)
    class(TwoBodyPart), intent(inout) :: this
    integer, intent(in) :: i1,i2,i3,i4,J
    real(8), intent(in) :: me
    type(TwoBodySpace), pointer :: tbs
    type(SingleParticleOrbit), pointer :: o1, o2, o3, o4
    integer :: P12, P34
    integer :: Z12, Z34
    integer :: ch, bra, ket, iphase

    tbs => this%two
    o1 => tbs%sps%orb(i1)
    o2 => tbs%sps%orb(i2)
    o3 => tbs%sps%orb(i3)
    o4 => tbs%sps%orb(i4)
    P12 = (-1) ** (o1%l + o2%l)
    P34 = (-1) ** (o3%l + o4%l)
    Z12 = (o1%z + o2%z)/2
    Z34 = (o3%z + o4%z)/2

    if(P12 * P34 /= 1) then
      write(*,*) "Warning: in SetTwBME_general: P"
      return
    end if

    if(Z12 - Z34 /= 0) then
      write(*,*) "Warning: in SetTwBME_general: Tz"
      return
    end if

    ch = tbs%jpz2ch(J,P12,Z12)
    if(ch == 0) return

    bra = tbs%jpz(ch)%spis2n(i1,i2)
    ket = tbs%jpz(ch)%spis2n(i3,i4)
    if(bra * ket == 0) return

    iphase = tbs%jpz(ch)%iphase(i1,i2) * &
        &    tbs%jpz(ch)%iphase(i3,i4)

    this%MatCh(ch,ch)%m(bra,ket) = dble(iphase) * me
    this%MatCh(ch,ch)%m(ket,bra) = dble(iphase) * me
  end subroutine SetTwBME_scalar

  subroutine AddToTwBME_tensor(this,i1,i2,i3,i4,J12,J34,me)
    use MyLibrary, only: triag
    class(TwoBodyPart), intent(inout) :: this
    integer, intent(in) :: i1,i2,i3,i4,J12,J34
    real(8), intent(in) :: me
    type(TwoBodySpace), pointer :: tbs
    type(SingleParticleOrbit), pointer :: o1, o2, o3, o4
    integer :: P12, P34
    integer :: Z12, Z34
    integer :: chbra, chket, bra, ket, iphase

    tbs => this%two
    o1 => tbs%sps%orb(i1)
    o2 => tbs%sps%orb(i2)
    o3 => tbs%sps%orb(i3)
    o4 => tbs%sps%orb(i4)
    P12 = (-1) ** (o1%l + o2%l)
    P34 = (-1) ** (o3%l + o4%l)
    Z12 = (o1%z + o2%z)/2
    Z34 = (o3%z + o4%z)/2

    if(triag(J12, J34, this%jr)) then
      write(*,*)"Warning: in AddToTwBME_general: J"
      return
    end if

    if(P12 * P34 * this%pr /= 1) then
      write(*,*) "Warning: in AddToTwBME_general: P"
      return
    end if

    if(Z12 - Z34 - this%zr /= 0) then
      write(*,*) "Warning: in AddToTwBME_general: Tz"
      return
    end if

    chbra = tbs%jpz2ch(J12,P12,Z12)
    chket = tbs%jpz2ch(J34,P34,Z34)
    if(chbra * chket == 0) return

    bra = tbs%jpz(chbra)%spis2n(i1,i2)
    ket = tbs%jpz(chket)%spis2n(i3,i4)
    if(bra * ket == 0) return

    iphase = tbs%jpz(chbra)%iphase(i1,i2) * &
        &    tbs%jpz(chket)%iphase(i3,i4)

    this%MatCh(chbra,chket)%m(bra,ket) = &
        & this%MatCh(chbra,chket)%m(bra,ket) + dble(iphase) * me
    this%MatCh(chket,chbra)%m(ket,bra) = this%MatCh(chbra,chket)%m(bra,ket)
  end subroutine AddToTwBME_tensor

  subroutine AddToTwBME_scalar(this,i1,i2,i3,i4,J,me)
    class(TwoBodyPart), intent(inout) :: this
    integer, intent(in) :: i1,i2,i3,i4,J
    real(8), intent(in) :: me
    type(TwoBodySpace), pointer :: tbs
    type(SingleParticleOrbit), pointer :: o1, o2, o3, o4
    integer :: P12, P34
    integer :: Z12, Z34
    integer :: ch, bra, ket, iphase

    tbs => this%two
    o1 => tbs%sps%orb(i1)
    o2 => tbs%sps%orb(i2)
    o3 => tbs%sps%orb(i3)
    o4 => tbs%sps%orb(i4)
    P12 = (-1) ** (o1%l + o2%l)
    P34 = (-1) ** (o3%l + o4%l)
    Z12 = (o1%z + o2%z)/2
    Z34 = (o3%z + o4%z)/2

    if(P12 * P34 /= 1) then
      write(*,*) "Warning: in AddToTwBME_general: P"
      return
    end if

    if(Z12 - Z34 /= 0) then
      write(*,*) "Warning: in AddToTwBME_general: Tz"
      return
    end if

    ch = tbs%jpz2ch(J,P12,Z12)
    if(ch == 0) return

    bra = tbs%jpz(ch)%spis2n(i1,i2)
    ket = tbs%jpz(ch)%spis2n(i3,i4)
    if(bra * ket == 0) return
    iphase = tbs%jpz(ch)%iphase(i1,i2) * &
        &    tbs%jpz(ch)%iphase(i3,i4)
    this%MatCh(ch,ch)%m(bra,ket) = this%MatCh(ch,ch)%m(bra,ket) + &
        & me * dble(iphase)
    this%MatCh(ch,ch)%m(ket,bra) = this%MatCh(ch,ch)%m(bra,ket)
  end subroutine AddToTwBME_scalar

  subroutine SetTwoBodyPartFromOneBody(two, one, mode)
    class(TwoBodyPart), intent(inout) :: two
    type(OneBodyPart), intent(in) :: one
    class(TwoBodySpace), pointer :: ms
    character(*), intent(in), optional :: mode
    integer :: chbra, chket
    ms => two%two
    do chbra = 1, ms%NChan
      do chket = 1, ms%NChan
        if(.not. two%MatCh(chbra,chket)%is) cycle
        call two%MatCh(chbra,chket)%set(one, mode)
      end do
    end do
  end subroutine SetTwoBodyPartFromOneBody

  subroutine SetTwoBodyPartChannelFromOneBody(two, one, mode)
    class(TwoBodyPartChannel), intent(inout) :: two
    type(OneBodyPart), intent(in) :: one
    character(*), intent(in), optional :: mode
    character(:), allocatable :: op

    op = "+"
    if(present(mode)) op = mode
    if(op == "+") then
      if(one%Scalar) call SetTwoBodyPartChannelFromOneBodyScalar(two, one)
      if(.not. one%Scalar) call SetTwoBodyPartChannelFromOneBodyTensor(two, one)
      return
    end if

    if(op == "*") then
      call SetTwoBodyChannelFromOneBodyTrans(two, one)
      return
    end if
  end subroutine SetTwoBodyPartChannelFromOneBody

  subroutine SetTwoBodyPartChannelFromOneBodyScalar(two, one)
    use Profiler, only: timer
    class(TwoBodyPartChannel), intent(inout) :: two
    type(OneBodyPart), intent(in) :: one
    class(TwoBodyChannel), pointer :: ch_bra, ch_ket
    integer :: bra, ket, a, b, c, d, jbra, jket
    type(SingleParticleOrbit), pointer :: oa, ob, oc, od
    real(8) :: norm, me, ti

    ti = omp_get_wtime()
    ch_bra => two%ch_bra
    ch_ket => two%ch_ket

    !$omp parallel
    !$omp do private(bra, jbra, a, b, oa, ob, ket, jket, c, d, oc, od, norm, me)
    do bra = 1, ch_bra%n_state
      jbra = ch_bra%j
      a = ch_bra%n2spi1(bra)
      b = ch_bra%n2spi2(bra)
      oa => one%one%sps%GetOrbit(a)
      ob => one%one%sps%GetOrbit(b)
      do ket = 1, ch_ket%n_state
        jket = ch_ket%j
        c = ch_bra%n2spi1(ket)
        d = ch_bra%n2spi2(ket)
        oc => one%one%sps%GetOrbit(c)
        od => one%one%sps%GetOrbit(d)
        norm = 1.d0
        if(a == b) norm = norm / dsqrt(2.d0)
        if(c == d) norm = norm / dsqrt(2.d0)
        me = 0.d0
        if(a == c) me = me + one%GetOBME(b,d)
        if(a == d) me = me - one%GetOBME(b,c) * (-1.d0)**((oa%j+ob%j)/2-jbra)
        if(b == c) me = me - one%GetOBME(a,d) * (-1.d0)**((oa%j+ob%j)/2-jbra)
        if(b == d) me = me + one%GetOBME(a,c)
        me = me * norm
        two%m(bra,ket) = me
      end do
    end do
    !$omp end do
    !$omp end parallel
    call timer%Add("SetTwoBodyPartChannelFromOneBodyScalar", omp_get_wtime()-ti)
  end subroutine SetTwoBodyPartChannelFromOneBodyScalar

  subroutine SetTwoBodyChannelFromOneBodyTrans(two, one)
    use Profiler, only: timer
    class(TwoBodyPartChannel), intent(inout) :: two
    type(OneBodyPart), intent(in) :: one
    class(TwoBodyChannel), pointer :: ch_bra, ch_ket
    integer :: bra, ket, a, b, c, d, j
    type(SingleParticleOrbit), pointer :: oa, ob, oc, od
    real(8) :: norm, me, ti

    ti = omp_get_wtime()
    ch_bra => two%ch_bra
    ch_ket => two%ch_ket
    j = ch_ket%j

    !$omp parallel
    !$omp do private(bra, a, b, oa, ob, ket, c, d, oc, od, norm, me)
    do bra = 1, ch_bra%n_state
      a = ch_bra%n2spi1(bra)
      b = ch_bra%n2spi2(bra)
      oa => one%one%sps%GetOrbit(a)
      ob => one%one%sps%GetOrbit(b)
      do ket = 1, ch_ket%n_state
        c = ch_bra%n2spi1(ket)
        d = ch_bra%n2spi2(ket)
        oc => one%one%sps%GetOrbit(c)
        od => one%one%sps%GetOrbit(d)
        norm = 1.d0
        if(a == b) norm = norm / dsqrt(2.d0)
        if(c == d) norm = norm / dsqrt(2.d0)
        me = 0.d0
        me = one%GetOBME(a,c) * one%GetOBME(b,d)
        me = me - one%GetOBME(a,d) * one%GetOBME(b,c) * &
            & (-1.d0)**((oa%j+ob%j)/2 - j)
        me = me * norm
        two%m(bra,ket) = me
      end do
    end do
    !$omp end do
    !$omp end parallel
    call timer%Add("SetTwoBodyChannelFromOneBodyTrans", omp_get_wtime()-ti)
  end subroutine SetTwoBodyChannelFromOneBodyTrans

  subroutine SetTwoBodyPartChannelFromOneBodyTensor(two, one)
    use MyLibrary, only: sjs
    use Profiler, only: timer
    class(TwoBodyPartChannel), intent(inout) :: two
    type(OneBodyPart), intent(in) :: one
    class(TwoBodyChannel), pointer :: ch_bra, ch_ket
    integer :: bra, ket, a, b, c, d, jbra, jket
    type(SingleParticleOrbit), pointer :: oa, ob, oc, od
    real(8) :: norm, me, ti

    ti = omp_get_wtime()
    ch_bra => two%ch_bra
    ch_ket => two%ch_ket

    !$omp parallel
    !$omp do private(bra, jbra, a, b, oa, ob, ket, jket, c, d, oc, od, norm, me)
    do bra = 1, ch_bra%n_state
      jbra = ch_bra%j
      a = ch_bra%n2spi1(bra)
      b = ch_bra%n2spi2(bra)
      oa => one%one%sps%GetOrbit(a)
      ob => one%one%sps%GetOrbit(b)
      do ket = 1, ch_ket%n_state
        jket = ch_ket%j
        c = ch_bra%n2spi1(ket)
        d = ch_bra%n2spi2(ket)
        oc => one%one%sps%GetOrbit(c)
        od => one%one%sps%GetOrbit(d)
        norm = 1.d0
        if(a == b) norm = norm / dsqrt(2.d0)
        if(c == d) norm = norm / dsqrt(2.d0)
        me = 0.d0
        ! tensor case not tested
        if(a == c) me = me + one%GetOBME(b,d) * &
            & (-1.d0)**((oc%j+od%j)/2-jbra) * sjs(2*jbra,2*jket,2*one%jr,od%j,ob%j,oa%j)
        if(a == d) me = me - one%GetOBME(b,c) * &
            & (-1.d0)**(         jket-jbra) * sjs(2*jbra,2*jket,2*one%jr,od%j,ob%j,oa%j)
        if(b == c) me = me - one%GetOBME(a,d) * &
            & (-1.d0)**(oa%j+ob%j+oc%j+od%j)/2 * sjs(2*jbra,2*jket,2*one%jr,od%j,oa%j,ob%j)
        if(b == d) me = me + one%GetOBME(a,c) * &
            & (-1.d0)**((oa%j+ob%j)/2+jket) * sjs(2*jbra,2*jket,2*one%jr,oc%j,oa%j,ob%j)
        me = me * norm * dsqrt(dble(2*jbra+1) * dble(2*jket+1)) * (-1.d0)**(one%jr)
        two%m(bra,ket) = me
      end do
    end do
    !$omp end do
    !$omp end parallel
    call timer%Add("SetTwoBodyPartChannelFromOneBodyTensor", omp_get_wtime()-ti)
  end subroutine SetTwoBodyPartChannelFromOneBodyTensor

  function ExtractTwoBodyChannel_pppp(opch, sps) result(Mat)
    class(TwoBodyPartChannel), intent(in) :: opch
    type(Orbits), intent(in) :: sps
    type(TwoBodyChannel), pointer :: ch_bra, ch_ket
    type(DMat) :: Mat
    integer :: ibra, iket, bra, ket
    integer :: a, b, c, d
    type(SingleParticleOrbit), pointer :: oa, ob, oc, od

    ch_bra => opch%ch_bra
    ch_ket => opch%ch_ket
    ibra = 0
    do bra = 1, ch_bra%n_state
      a = ch_bra%n2spi1(bra)
      b = ch_bra%n2spi2(bra)
      oa => sps%GetOrbit(a)
      ob => sps%GetOrbit(b)
      if(abs(oa%occ)+abs(ob%occ) < 1.d-6) ibra = ibra+1
    end do

    iket = 0
    do bra = 1, ch_ket%n_state
      a = ch_ket%n2spi1(bra)
      b = ch_ket%n2spi2(bra)
      oa => sps%GetOrbit(a)
      ob => sps%GetOrbit(b)
      if(abs(oa%occ)+abs(ob%occ) < 1.d-6) iket = iket+1
    end do
    call Mat%ini(ibra,iket)
    ibra = 0
    do bra = 1, ch_bra%n_state
      a = ch_bra%n2spi1(bra)
      b = ch_bra%n2spi2(bra)
      oa => sps%GetOrbit(a)
      ob => sps%GetOrbit(b)
      if(abs(oa%occ+ob%occ) > 1.d-6) cycle
      ibra = ibra + 1
      iket = 0
      do ket = 1, ch_ket%n_state
        c = ch_ket%n2spi1(ket)
        d = ch_ket%n2spi2(ket)
        oc => sps%GetOrbit(c)
        od => sps%GetOrbit(d)
        if(abs(oc%occ+od%occ) > 1.d-6) cycle
        iket = iket + 1
        Mat%m(ibra,iket) = opch%m(bra,ket)
      end do
    end do
  end function ExtractTwoBodyChannel_pppp

  function ExtractTwoBodyChannel_hhhh(opch, sps) result(Mat)
    class(TwoBodyPartChannel), intent(in) :: opch
    type(Orbits), intent(in) :: sps
    type(DMat) :: Mat
    type(TwoBodyChannel), pointer :: ch_bra, ch_ket
    integer :: ibra, iket, bra, ket
    integer :: a, b, c, d
    type(SingleParticleOrbit), pointer :: oa, ob, oc, od

    ch_bra => opch%ch_bra
    ch_ket => opch%ch_ket
    ibra = 0
    do bra = 1, ch_bra%n_state
      a = ch_bra%n2spi1(bra)
      b = ch_bra%n2spi2(bra)
      oa => sps%GetOrbit(a)
      ob => sps%GetOrbit(b)
      if(abs(oa%occ)*abs(ob%occ) > 1.d-6) ibra = ibra+1
    end do

    iket = 0
    do bra = 1, ch_ket%n_state
      a = ch_ket%n2spi1(bra)
      b = ch_ket%n2spi2(bra)
      oa => sps%GetOrbit(a)
      ob => sps%GetOrbit(b)
      if(abs(oa%occ)*abs(ob%occ) > 1.d-6) iket = iket+1
    end do

    call Mat%ini(ibra,iket)
    ibra = 0
    do bra = 1, ch_bra%n_state
      a = ch_bra%n2spi1(bra)
      b = ch_bra%n2spi2(bra)
      oa => sps%GetOrbit(a)
      ob => sps%GetOrbit(b)
      if(abs(oa%occ)*abs(ob%occ) < 1.d-6) cycle
      ibra = ibra + 1
      iket = 0
      do ket = 1, ch_ket%n_state
        c = ch_ket%n2spi1(ket)
        d = ch_ket%n2spi2(ket)
        oc => sps%GetOrbit(c)
        od => sps%GetOrbit(d)
        if(abs(oc%occ)*abs(od%occ) < 1.d-6) cycle
        iket = iket + 1
        Mat%m(ibra,iket) = opch%m(bra,ket)
      end do
    end do
  end function ExtractTwoBodyChannel_hhhh

  function ExtractTwoBodyChannel_pphh(opch, sps) result(Mat)
    class(TwoBodyPartChannel), intent(in) :: opch
    type(Orbits), intent(in) :: sps
    type(TwoBodyChannel), pointer :: ch_bra, ch_ket
    type(DMat) :: Mat
    integer :: ibra, iket, bra, ket
    integer :: a, b, c, d
    type(SingleParticleOrbit), pointer :: oa, ob, oc, od

    ch_bra => opch%ch_bra
    ch_ket => opch%ch_ket
    ibra = 0
    do bra = 1, ch_bra%n_state
      a = ch_bra%n2spi1(bra)
      b = ch_bra%n2spi2(bra)
      oa => sps%GetOrbit(a)
      ob => sps%GetOrbit(b)
      if(abs(oa%occ)+abs(ob%occ) < 1.d-6) ibra = ibra+1
    end do

    iket = 0
    do bra = 1, ch_ket%n_state
      a = ch_ket%n2spi1(bra)
      b = ch_ket%n2spi2(bra)
      oa => sps%GetOrbit(a)
      ob => sps%GetOrbit(b)
      if(abs(oa%occ)*abs(ob%occ) > 1.d-6) iket = iket+1
    end do

    call Mat%ini(ibra,iket)
    ibra = 0
    do bra = 1, ch_bra%n_state
      a = ch_bra%n2spi1(bra)
      b = ch_bra%n2spi2(bra)
      oa => sps%GetOrbit(a)
      ob => sps%GetOrbit(b)
      if(abs(oa%occ)+abs(ob%occ) > 1.d-6) cycle
      ibra = ibra + 1
      iket = 0
      do ket = 1, ch_ket%n_state
        c = ch_ket%n2spi1(ket)
        d = ch_ket%n2spi2(ket)
        oc => sps%GetOrbit(c)
        od => sps%GetOrbit(d)
        if(abs(oc%occ)*abs(od%occ) < 1.d-6) cycle
        iket = iket + 1
        Mat%m(ibra,iket) = opch%m(bra,ket)
      end do
    end do
  end function ExtractTwoBodyChannel_pphh

  function ExtractTwoBodyChannel_phhh(opch, sps) result(Mat)
    class(TwoBodyPartChannel), intent(in) :: opch
    type(Orbits), intent(in) :: sps
    type(DMat) :: Mat
    type(TwoBodyChannel), pointer :: ch_bra, ch_ket
    integer :: ibra, iket, bra, ket
    integer :: a, b, c, d
    type(SingleParticleOrbit), pointer :: oa, ob, oc, od

    ch_bra => opch%ch_bra
    ch_ket => opch%ch_ket
    ibra = 0
    do bra = 1, ch_bra%n_state
      a = ch_bra%n2spi1(bra)
      b = ch_bra%n2spi2(bra)
      oa => sps%GetOrbit(a)
      ob => sps%GetOrbit(b)
      if(abs(oa%occ)+abs(ob%occ) > 1.d-6 .and. abs(oa%occ)*abs(ob%occ) < 1.d-6) ibra = ibra+1
    end do

    iket = 0
    do bra = 1, ch_ket%n_state
      a = ch_ket%n2spi1(bra)
      b = ch_ket%n2spi2(bra)
      oa => sps%GetOrbit(a)
      ob => sps%GetOrbit(b)
      if(abs(oa%occ)*abs(ob%occ) > 1.d-6) iket = iket+1
    end do

    call Mat%ini(ibra,iket)
    ibra = 0
    do bra = 1, ch_bra%n_state
      a = ch_bra%n2spi1(bra)
      b = ch_bra%n2spi2(bra)
      oa => sps%GetOrbit(a)
      ob => sps%GetOrbit(b)
      if(abs(oa%occ)+abs(ob%occ) < 1.d-6 .or. abs(oa%occ)*abs(oa%occ) > 1.d-6) cycle
      ibra = ibra + 1
      iket = 0
      do ket = 1, ch_ket%n_state
        c = ch_ket%n2spi1(ket)
        d = ch_ket%n2spi2(ket)
        oc => sps%GetOrbit(c)
        od => sps%GetOrbit(d)
        if(abs(oc%occ)*abs(od%occ) < 1.d-6) cycle
        iket = iket + 1
        Mat%m(ibra,iket) = opch%m(bra,ket)
      end do
    end do
  end function ExtractTwoBodyChannel_phhh

  function ExtractTwoBodyChannel_phpp(opch, sps) result(Mat)
    class(TwoBodyPartChannel), intent(in) :: opch
    type(Orbits), intent(in) :: sps
    type(DMat) :: Mat
    type(TwoBodyChannel), pointer :: ch_bra, ch_ket
    integer :: ibra, iket, bra, ket
    integer :: a, b, c, d
    type(SingleParticleOrbit), pointer :: oa, ob, oc, od

    ch_bra => opch%ch_bra
    ch_ket => opch%ch_ket
    ibra = 0
    do bra = 1, ch_bra%n_state
      a = ch_bra%n2spi1(bra)
      b = ch_bra%n2spi2(bra)
      oa => sps%GetOrbit(a)
      ob => sps%GetOrbit(b)
      if(abs(oa%occ)+abs(ob%occ) > 1.d-6 .and. abs(oa%occ)*abs(ob%occ) < 1.d-6) ibra = ibra+1
    end do

    iket = 0
    do bra = 1, ch_ket%n_state
      a = ch_ket%n2spi1(bra)
      b = ch_ket%n2spi2(bra)
      oa => sps%GetOrbit(a)
      ob => sps%GetOrbit(b)
      if(abs(oa%occ)+abs(ob%occ) < 1.d-6) iket=iket+1
    end do

    call Mat%ini(ibra,iket)
    ibra = 0
    do bra = 1, ch_bra%n_state
      a = ch_bra%n2spi1(bra)
      b = ch_bra%n2spi2(bra)
      oa => sps%GetOrbit(a)
      ob => sps%GetOrbit(b)
      if(abs(oa%occ)+abs(ob%occ) < 1.d-6 .or. abs(oa%occ)*abs(oa%occ) > 1.d-6) cycle
      ibra = ibra + 1
      iket = 0
      do ket = 1, ch_ket%n_state
        c = ch_ket%n2spi1(ket)
        d = ch_ket%n2spi2(ket)
        oc => sps%GetOrbit(c)
        od => sps%GetOrbit(d)
        if(abs(oc%occ)+abs(od%occ) > 1.d-6) cycle
        iket = iket + 1
        Mat%m(ibra,iket) = opch%m(bra,ket)
      end do
    end do
  end function ExtractTwoBodyChannel_phpp

  subroutine SetTwoBodyChannel_pppp(opch, sps, Mat)
    class(TwoBodyPartChannel), intent(inout) :: opch
    type(Orbits), intent(in) :: sps
    type(DMat), intent(in) :: Mat
    type(TwoBodyChannel), pointer :: ch_bra, ch_ket
    integer :: ibra, iket, bra, ket
    integer :: a, b, c, d
    type(SingleParticleOrbit), pointer :: oa, ob, oc, od

    ch_bra => opch%ch_bra
    ch_ket => opch%ch_ket
    ibra = 0
    do bra = 1, ch_bra%n_state
      a = ch_bra%n2spi1(bra)
      b = ch_bra%n2spi2(bra)
      oa => sps%GetOrbit(a)
      ob => sps%GetOrbit(b)
      if(abs(oa%occ)+abs(ob%occ) < 1.d-6) ibra = ibra+1
    end do

    iket = 0
    do bra = 1, ch_ket%n_state
      a = ch_ket%n2spi1(bra)
      b = ch_ket%n2spi2(bra)
      oa => sps%GetOrbit(a)
      ob => sps%GetOrbit(b)
      if(abs(oa%occ)+abs(ob%occ) < 1.d-6) iket = iket+1
    end do

    ibra = 0
    do bra = 1, ch_bra%n_state
      a = ch_bra%n2spi1(bra)
      b = ch_bra%n2spi2(bra)
      oa => sps%GetOrbit(a)
      ob => sps%GetOrbit(b)
      if(abs(oa%occ+ob%occ) > 1.d-6) cycle
      ibra = ibra + 1
      iket = 0
      do ket = 1, ch_ket%n_state
        c = ch_ket%n2spi1(ket)
        d = ch_ket%n2spi2(ket)
        oc => sps%GetOrbit(c)
        od => sps%GetOrbit(d)
        if(abs(oc%occ+od%occ) > 1.d-6) cycle
        iket = iket + 1
        opch%m(bra,ket) = Mat%m(ibra,iket)
        opch%m(ket,bra) = Mat%m(ibra,iket)
      end do
    end do
  end subroutine SetTwoBodyChannel_pppp

  subroutine SetTwoBodyChannel_hhhh(opch, sps, Mat)
    class(TwoBodyPartChannel), intent(inout) :: opch
    type(Orbits), intent(in) :: sps
    type(DMat), intent(in) :: Mat
    type(TwoBodyChannel), pointer :: ch_bra, ch_ket
    integer :: ibra, iket, bra, ket
    integer :: a, b, c, d
    type(SingleParticleOrbit), pointer :: oa, ob, oc, od

    ch_bra => opch%ch_bra
    ch_ket => opch%ch_ket
    ibra = 0
    do bra = 1, ch_bra%n_state
      a = ch_bra%n2spi1(bra)
      b = ch_bra%n2spi2(bra)
      oa => sps%GetOrbit(a)
      ob => sps%GetOrbit(b)
      if(abs(oa%occ)*abs(ob%occ) > 1.d-6) ibra = ibra+1
    end do

    iket = 0
    do bra = 1, ch_ket%n_state
      a = ch_ket%n2spi1(bra)
      b = ch_ket%n2spi2(bra)
      oa => sps%GetOrbit(a)
      ob => sps%GetOrbit(b)
      if(abs(oa%occ)*abs(ob%occ) > 1.d-6) iket = iket+1
    end do

    ibra = 0
    do bra = 1, ch_bra%n_state
      a = ch_bra%n2spi1(bra)
      b = ch_bra%n2spi2(bra)
      oa => sps%GetOrbit(a)
      ob => sps%GetOrbit(b)
      if(abs(oa%occ)*abs(ob%occ) < 1.d-6) cycle
      ibra = ibra + 1
      iket = 0
      do ket = 1, ch_ket%n_state
        c = ch_ket%n2spi1(ket)
        d = ch_ket%n2spi2(ket)
        oc => sps%GetOrbit(c)
        od => sps%GetOrbit(d)
        if(abs(oc%occ)*abs(od%occ) < 1.d-6) cycle
        iket = iket + 1
        opch%m(bra,ket) = Mat%m(ibra,iket)
        opch%m(ket,bra) = Mat%m(ibra,iket)
      end do
    end do
  end subroutine SetTwoBodyChannel_hhhh

  subroutine SetTwoBodyChannel_pphh(opch, sps, Mat)
    class(TwoBodyPartChannel), intent(inout) :: opch
    type(Orbits), intent(in) :: sps
    type(DMat), intent(in) :: Mat
    type(TwoBodyChannel), pointer :: ch_bra, ch_ket
    integer :: ibra, iket, bra, ket
    integer :: a, b, c, d
    type(SingleParticleOrbit), pointer :: oa, ob, oc, od

    ch_bra => opch%ch_bra
    ch_ket => opch%ch_ket
    ibra = 0
    do bra = 1, ch_bra%n_state
      a = ch_bra%n2spi1(bra)
      b = ch_bra%n2spi2(bra)
      oa => sps%GetOrbit(a)
      ob => sps%GetOrbit(b)
      if(abs(oa%occ)+abs(ob%occ) < 1.d-6) ibra = ibra+1
    end do

    iket = 0
    do bra = 1, ch_ket%n_state
      a = ch_ket%n2spi1(bra)
      b = ch_ket%n2spi2(bra)
      oa => sps%GetOrbit(a)
      ob => sps%GetOrbit(b)
      if(abs(oa%occ)*abs(ob%occ) > 1.d-6) iket = iket+1
    end do

    ibra = 0
    do bra = 1, ch_bra%n_state
      a = ch_bra%n2spi1(bra)
      b = ch_bra%n2spi2(bra)
      oa => sps%GetOrbit(a)
      ob => sps%GetOrbit(b)
      if(abs(oa%occ)+abs(ob%occ) > 1.d-6) cycle
      ibra = ibra + 1
      iket = 0
      do ket = 1, ch_ket%n_state
        c = ch_ket%n2spi1(ket)
        d = ch_ket%n2spi2(ket)
        oc => sps%GetOrbit(c)
        od => sps%GetOrbit(d)
        if(abs(oc%occ)*abs(od%occ) < 1.d-6) cycle
        iket = iket + 1
        opch%m(bra,ket) = Mat%m(ibra,iket)
        opch%m(ket,bra) = Mat%m(ibra,iket)
      end do
    end do
  end subroutine SetTwoBodyChannel_pphh

  subroutine SetTwoBodyChannel_phhh(opch, sps, Mat)
    class(TwoBodyPartChannel), intent(inout) :: opch
    type(DMat), intent(in) :: Mat
    type(Orbits), intent(in) :: sps
    type(TwoBodyChannel), pointer :: ch_bra, ch_ket
    integer :: ibra, iket, bra, ket
    integer :: a, b, c, d
    type(SingleParticleOrbit), pointer :: oa, ob, oc, od

    ch_bra => opch%ch_bra
    ch_ket => opch%ch_ket
    ibra = 0
    do bra = 1, ch_bra%n_state
      a = ch_bra%n2spi1(bra)
      b = ch_bra%n2spi2(bra)
      oa => sps%GetOrbit(a)
      ob => sps%GetOrbit(b)
      if(abs(oa%occ)+abs(ob%occ) > 1.d-6 .and. abs(oa%occ)*abs(ob%occ) < 1.d-6) ibra = ibra+1
    end do

    iket = 0
    do bra = 1, ch_ket%n_state
      a = ch_ket%n2spi1(bra)
      b = ch_ket%n2spi2(bra)
      oa => sps%GetOrbit(a)
      ob => sps%GetOrbit(b)
      if(abs(oa%occ)*abs(ob%occ) > 1.d-6) iket = iket+1
    end do

    ibra = 0
    do bra = 1, ch_bra%n_state
      a = ch_bra%n2spi1(bra)
      b = ch_bra%n2spi2(bra)
      oa => sps%GetOrbit(a)
      ob => sps%GetOrbit(b)
      if(abs(oa%occ)+abs(ob%occ) < 1.d-6 .or. abs(oa%occ)*abs(oa%occ) > 1.d-6) cycle
      ibra = ibra + 1
      iket = 0
      do ket = 1, ch_ket%n_state
        c = ch_ket%n2spi1(ket)
        d = ch_ket%n2spi2(ket)
        oc => sps%GetOrbit(c)
        od => sps%GetOrbit(d)
        if(abs(oc%occ)*abs(od%occ) < 1.d-6) cycle
        iket = iket + 1
        opch%m(bra,ket) = Mat%m(ibra,iket)
        opch%m(ket,bra) = Mat%m(ibra,iket)
      end do
    end do
  end subroutine SetTwoBodyChannel_phhh

  subroutine SetTwoBodyChannel_phpp(opch, sps, Mat)
    class(TwoBodyPartChannel), intent(inout) :: opch
    type(Orbits), intent(in) :: sps
    type(DMat), intent(in) :: Mat
    type(TwoBodyChannel), pointer :: ch_bra, ch_ket
    integer :: ibra, iket, bra, ket
    integer :: a, b, c, d
    type(SingleParticleOrbit), pointer :: oa, ob, oc, od

    ch_bra => opch%ch_bra
    ch_ket => opch%ch_ket
    ibra = 0
    do bra = 1, ch_bra%n_state
      a = ch_bra%n2spi1(bra)
      b = ch_bra%n2spi2(bra)
      oa => sps%GetOrbit(a)
      ob => sps%GetOrbit(b)
      if(abs(oa%occ)+abs(ob%occ) > 1.d-6 .and. abs(oa%occ)*abs(ob%occ) < 1.d-6) ibra = ibra+1
    end do

    iket = 0
    do bra = 1, ch_ket%n_state
      a = ch_ket%n2spi1(bra)
      b = ch_ket%n2spi2(bra)
      oa => sps%GetOrbit(a)
      ob => sps%GetOrbit(b)
      if(abs(oa%occ)+abs(ob%occ) < 1.d-6) iket=iket+1
    end do

    ibra = 0
    do bra = 1, ch_bra%n_state
      a = ch_bra%n2spi1(bra)
      b = ch_bra%n2spi2(bra)
      oa => sps%GetOrbit(a)
      ob => sps%GetOrbit(b)
      if(abs(oa%occ)+abs(ob%occ) < 1.d-6 .or. abs(oa%occ)*abs(oa%occ) > 1.d-6) cycle
      ibra = ibra + 1
      iket = 0
      do ket = 1, ch_ket%n_state
        c = ch_ket%n2spi1(ket)
        d = ch_ket%n2spi2(ket)
        oc => sps%GetOrbit(c)
        od => sps%GetOrbit(d)
        if(abs(oc%occ)+abs(od%occ) > 1.d-6) cycle
        iket = iket + 1
        opch%m(bra,ket) = Mat%m(ibra,iket)
        opch%m(ket,bra) = Mat%m(ibra,iket)
      end do
    end do
  end subroutine SetTwoBodyChannel_phpp

  function ExtractCrossCoupledChannel_hphp2phph(op, ch_cc) result(Mat)
    !  only for scalar
    !  _________________
    !  <ph:J| V |p'h':J> = \sum_{J'} [J'] {jh  jp' J'} <hp':J'|V|h'p:J'>
    !                                     {jh' jp  J }
    use MyLibrary, only: sjs
    class(TwoBodyPart), intent(in) :: op
    type(CrossCoupledTwoBodyChannel), intent(in) :: ch_cc
    type(Orbits), pointer :: sps
    type(DMat) :: Mat
    integer :: a, b, c, d, K
    integer :: ibra, iket, bra, ket
    type(SingleParticleOrbit), pointer :: oa, ob, oc, od

    sps => op%two%sps
    K = ch_cc%j
    ibra = 0
    do bra = 1, ch_cc%n_state
      a = ch_cc%n2spi1(bra)
      b = ch_cc%n2spi2(bra)
      oa => sps%GetOrbit(a)
      ob => sps%GetOrbit(b)
      if(abs(oa%occ)+abs(ob%occ) > 1.d-6 .and. abs(oa%occ)*abs(ob%occ) < 1.d-6) ibra = ibra+1
    end do

    call Mat%zeros(ibra,ibra)
    if(ibra < 1) return
    ibra = 0
    do bra = 1, ch_cc%n_state
      a = ch_cc%n2spi1(bra) ! p
      b = ch_cc%n2spi2(bra) ! h
      oa => sps%GetOrbit(a)
      ob => sps%GetOrbit(b)
      if(abs(oa%occ)+abs(ob%occ) < 1.d-6 .or. abs(oa%occ)*abs(ob%occ) > 1.d-6) cycle
      ibra = ibra+1
      iket = 0
      do ket = 1, ch_cc%n_state
        c = ch_cc%n2spi1(ket) ! p
        d = ch_cc%n2spi2(ket) ! h
        oc => sps%GetOrbit(c)
        od => sps%GetOrbit(d)
        if(abs(oc%occ)+abs(od%occ) < 1.d-6 .or. abs(oc%occ)*abs(od%occ) > 1.d-6) cycle
        iket = iket+1
        if(abs(oa%occ)+abs(od%occ) < 1.d-6 .or. abs(oa%occ)*abs(od%occ) > 1.d-6) cycle
        if(abs(ob%occ)+abs(oc%occ) < 1.d-6 .or. abs(ob%occ)*abs(oc%occ) > 1.d-6) cycle
        if(oa%z+od%z /= ob%z+oc%z) cycle
        Mat%m(ibra,iket) = op%get_xc1423(b,c,d,a,K)
      end do
    end do
  end function ExtractCrossCoupledChannel_hphp2phph

  function ExtractCrossCoupledChannel_phph2phph(op, ch_cc) result(Mat)
    !  only for scalar
    !  _________________
    !  <ph:J| V |p'h':J> = \sum_{J'} [J'] {jp  jh' J'} <ph':J'|V|p'h:J'>
    !                                     {jp' jh  J }
    use MyLibrary, only: sjs
    class(TwoBodyPart), intent(in) :: op
    type(CrossCoupledTwoBodyChannel), intent(in) :: ch_cc
    type(Orbits), pointer :: sps
    type(DMat) :: Mat
    integer :: a, b, c, d, K
    integer :: ibra, iket, bra, ket
    type(SingleParticleOrbit), pointer :: oa, ob, oc, od

    sps => op%two%sps
    K = ch_cc%j
    ibra = 0
    do bra = 1, ch_cc%n_state
      a = ch_cc%n2spi1(bra)
      b = ch_cc%n2spi2(bra)
      oa => sps%GetOrbit(a)
      ob => sps%GetOrbit(b)
      if(abs(oa%occ)+abs(ob%occ) > 1.d-6 .and. abs(oa%occ)*abs(ob%occ) < 1.d-6) ibra = ibra+1
    end do

    call Mat%zeros(ibra,ibra)
    if(ibra < 1) return
    ibra = 0
    do bra = 1, ch_cc%n_state
      a = ch_cc%n2spi1(bra) ! p
      b = ch_cc%n2spi2(bra) ! h
      oa => sps%GetOrbit(a)
      ob => sps%GetOrbit(b)
      if(abs(oa%occ)+abs(ob%occ) < 1.d-6 .or. abs(oa%occ)*abs(ob%occ) > 1.d-6) cycle
      ibra = ibra+1
      iket = 0
      do ket = 1, ch_cc%n_state
        c = ch_cc%n2spi1(ket) ! p
        d = ch_cc%n2spi2(ket) ! h
        oc => sps%GetOrbit(c)
        od => sps%GetOrbit(d)
        if(abs(oc%occ)+abs(od%occ) < 1.d-6 .or. abs(oc%occ)*abs(od%occ) > 1.d-6) cycle
        iket = iket+1
        if(abs(oa%occ)+abs(od%occ) < 1.d-6 .or. abs(oa%occ)*abs(od%occ) > 1.d-6) cycle
        if(abs(ob%occ)+abs(oc%occ) < 1.d-6 .or. abs(ob%occ)*abs(oc%occ) > 1.d-6) cycle
        if(oa%z+od%z /= ob%z+oc%z) cycle
        Mat%m(ibra,iket) = op%get_xc1423(a,d,c,b,K)
      end do
    end do
  end function ExtractCrossCoupledChannel_phph2phph

  function ExtractCrossCoupledChannel_pphh2phph(op, ch_cc) result(Mat)
    !  only for scalar
    !  _________________
    !  <ph:J| V |p'h':J> = \sum_{J'} [J'] {jp  jp' J'} <pp':J'|V|h'h:J'>
    !                                     {jh' jh  J }
    class(TwoBodyPart), intent(in) :: op
    type(CrossCoupledTwoBodyChannel), intent(in) :: ch_cc
    type(Orbits), pointer :: sps
    type(DMat) :: Mat
    integer :: a, b, c, d, K
    integer :: ibra, iket, bra, ket
    type(SingleParticleOrbit), pointer :: oa, ob, oc, od
    real(8) :: v, norm

    sps => op%two%sps
    K = ch_cc%j
    ibra = 0
    do bra = 1, ch_cc%n_state
      a = ch_cc%n2spi1(bra)
      b = ch_cc%n2spi2(bra)
      oa => sps%GetOrbit(a)
      ob => sps%GetOrbit(b)
      if(abs(oa%occ)+abs(ob%occ) > 1.d-6 .and. abs(oa%occ)*abs(ob%occ) < 1.d-6) ibra = ibra+1
    end do
    call Mat%zeros(ibra,ibra)
    if(ibra < 1) return
    ibra = 0
    do bra = 1, ch_cc%n_state
      a = ch_cc%n2spi1(bra) ! p
      b = ch_cc%n2spi2(bra) ! h
      oa => sps%GetOrbit(a)
      ob => sps%GetOrbit(b)
      if(abs(oa%occ)+abs(ob%occ) < 1.d-6 .or. abs(oa%occ)*abs(oa%occ) > 1.d-6) cycle
      ibra = ibra+1
      iket = 0
      do ket = 1, ch_cc%n_state
        c = ch_cc%n2spi1(ket) ! p
        d = ch_cc%n2spi2(ket) ! h
        oc => sps%GetOrbit(c)
        od => sps%GetOrbit(d)
        if(abs(oc%occ)+abs(od%occ) < 1.d-6 .or. abs(oc%occ)*abs(od%occ) > 1.d-6) cycle
        iket = iket+1
        if(oa%z+oc%z /= ob%z+od%z) cycle
        if(abs(oa%occ)+abs(oc%occ) > 1.d-6) cycle
        if(abs(ob%occ)*abs(od%occ) < 1.d-6) cycle
        norm = 1.d0
        if(a==c) norm = norm*sqrt(2.d0)
        if(b==d) norm = norm*sqrt(2.d0)
        v = op%get_xc1423(a,c,d,b,K)
        Mat%m(ibra,iket) = v * norm
      end do
    end do
  end function ExtractCrossCoupledChannel_pphh2phph

!  function ExtractCrossCoupledChannel_phhh2phhh(op, ch_cc, f) result(Mat)
!    !  only for scalar
!    !  __________________
!    !  <ph:J| V |h'h'':J> = \sum_{J'} [J'] {jp   jh'  J'} <ph':J'|V|hh'':J'> / f_h - f_p
!    !                                      {jh'' jh   J }
!    use MyLibrary, only: sjs
!    class(TwoBodyPart), intent(in) :: op
!    type(CrossCoupledTwoBodyChannel), intent(in) :: ch_cc
!    type(OneBodyPart), intent(in), optional :: f
!    type(Orbits), pointer :: sps
!    type(DMat) :: Mat
!    integer :: a, b, c, d, K
!    integer :: Jmin, Jmax, J
!    integer :: ibra, iket, bra, ket
!    logical :: denom = .false.
!    type(SingleParticleOrbit), pointer :: oa, ob, oc, od
!    real(8) :: v, eps, norm
!
!    if(present(f)) denom = .true.
!    sps => op%two%sps
!    K = ch_cc%j
!    ibra = 0
!    do bra = 1, ch_cc%n_state
!      a = ch_cc%n2spi1(bra)
!      b = ch_cc%n2spi2(bra)
!      oa => sps%GetOrbit(a)
!      ob => sps%GetOrbit(b)
!      if(abs(oa%occ)+abs(ob%occ) > 1.d-6 .and. abs(oa%occ)*abs(ob%occ) < 1.d-6) ibra = ibra+1
!    end do
!
!    iket = 0
!    do bra = 1, ch_cc%n_state
!      a = ch_cc%n2spi1(bra)
!      b = ch_cc%n2spi2(bra)
!      oa => sps%GetOrbit(a)
!      ob => sps%GetOrbit(b)
!      if(abs(oa%occ)*abs(ob%occ) > 1.d-6) iket = iket+1
!    end do
!
!    call Mat%zeros(ibra,iket)
!    if(ibra < 1 .or. iket < 1) return
!    ibra = 0
!    do bra = 1, ch_cc%n_state
!      a = ch_cc%n2spi1(bra) ! p
!      b = ch_cc%n2spi2(bra) ! h
!      oa => sps%GetOrbit(a)
!      ob => sps%GetOrbit(b)
!      if(abs(oa%occ)+abs(ob%occ) < 1.d-6 .or. abs(oa%occ)*abs(ob%occ) > 1.d-6) cycle
!      ibra = ibra+1
!      iket = 0
!      do ket = 1, ch_cc%n_state
!        c = ch_cc%n2spi1(ket) ! p
!        d = ch_cc%n2spi2(ket) ! h
!        oc => sps%GetOrbit(c)
!        od => sps%GetOrbit(d)
!        if(abs(oc%occ)*abs(od%occ) < 1.d-6) cycle
!        iket = iket + 1
!
!        if(abs(oa%occ)+abs(oc%occ) < 1.d-6 .or. abs(oa%occ)*abs(oc%occ) > 1.d-6) cycle
!        if(abs(ob%occ)*abs(od%occ) < 1.d-6) cycle
!        if(oa%z+oc%z /= ob%z+od%z) cycle
!        norm = 1.d0
!        if(b==d) norm = norm*sqrt(2.d0)
!        v = op%get_xc1324(a,c,b,d,K)
!        eps = 1.d0
!        if(denom) eps = 1.d0 / f%GetDenominator1(b,a)
!        Mat%m(ibra,iket) = v * eps * norm
!      end do
!    end do
!  end function ExtractCrossCoupledChannel_phhh2phhh
!
!  function ExtractCrossCoupledChannel_ppph2ppph(op, ch_cc, f) result(Mat)
!    !  only for scalar
!    !  __________________
!    !  <pp':J| V |p''h:J> = \sum_{J'} [J'] {jp   jp''  J'} <pp'':J'|V|p'h:J'> / f_h - f_p''
!    !                                      {jh   jp'   J }
!    use MyLibrary, only: sjs
!    class(TwoBodyPart), intent(in) :: op
!    type(CrossCoupledTwoBodyChannel), intent(in) :: ch_cc
!    type(OneBodyPart), intent(in), optional :: f
!    type(Orbits), pointer :: sps
!    type(DMat) :: Mat
!    integer :: a, b, c, d, K
!    integer :: Jmin, Jmax, J
!    integer :: ibra, iket, bra, ket
!    logical :: denom = .false.
!    type(SingleParticleOrbit), pointer :: oa, ob, oc, od
!    real(8) :: v, eps, norm
!
!    if(present(f)) denom = .true.
!    sps => op%two%sps
!    K = ch_cc%j
!    ibra = 0
!    do bra = 1, ch_cc%n_state
!      a = ch_cc%n2spi1(bra)
!      b = ch_cc%n2spi2(bra)
!      oa => sps%GetOrbit(a)
!      ob => sps%GetOrbit(b)
!      if(abs(oa%occ)+abs(ob%occ) < 1.d-6) ibra = ibra+1
!    end do
!
!    iket = 0
!    do bra = 1, ch_cc%n_state
!      a = ch_cc%n2spi1(bra)
!      b = ch_cc%n2spi2(bra)
!      oa => sps%GetOrbit(a)
!      ob => sps%GetOrbit(b)
!      if(abs(oa%occ)+abs(ob%occ) > 1.d-6 .and. abs(oa%occ)*abs(ob%occ) < 1.d-6) iket = iket+1
!    end do
!
!    call Mat%zeros(ibra,iket)
!    if(ibra < 1 .or. iket < 1) return
!    ibra = 0
!    do bra = 1, ch_cc%n_state
!      a = ch_cc%n2spi1(bra) ! p
!      b = ch_cc%n2spi2(bra) ! h
!      oa => sps%GetOrbit(a)
!      ob => sps%GetOrbit(b)
!      if(abs(oa%occ)+abs(ob%occ) > 1.d-6) cycle
!      ibra = ibra+1
!      iket = 0
!      do ket = 1, ch_cc%n_state
!        c = ch_cc%n2spi1(ket) ! p
!        d = ch_cc%n2spi2(ket) ! h
!        oc => sps%GetOrbit(c)
!        od => sps%GetOrbit(d)
!        if(abs(oc%occ)+abs(od%occ) < 1.d-6 .or. abs(oc%occ)*abs(od%occ) > 1.d-6) cycle
!        iket = iket + 1
!
!        if(abs(oa%occ)+abs(oc%occ) > 1.d-6) cycle
!        if(abs(ob%occ)+abs(od%occ) < 1.d-6 .or. abs(ob%occ)*abs(od%occ) > 1.d-6) cycle
!        if(oa%z+oc%z /= ob%z+od%z) cycle
!        norm = 1.d0
!        if(a==c) norm = norm*sqrt(2.d0)
!        v = op%get_xc1324(a,c,b,d,K)
!        eps = 1.d0
!        if(denom) eps = 1.d0 / f%GetDenominator1(d,b)
!        Mat%m(ibra,iket) = v * eps * norm
!      end do
!    end do
!  end function ExtractCrossCoupledChannel_ppph2ppph

  ! Only for Scalar
  function GetCrossCoupledME1423(this, a, b, c, d, J) result(r)
    use MyLibrary, only: sjs
    class(TwoBodyPart), intent(in) :: this
    integer, intent(in) :: a, b, c, d, J
    real(8) :: r
    type(Orbits), pointer :: sps
    type(SingleParticleOrbit), pointer :: oa, ob, oc, od
    integer :: K, Kmin, Kmax
    sps => this%two%sps
    oa => sps%GetOrbit(a)
    ob => sps%GetOrbit(b)
    oc => sps%GetOrbit(c)
    od => sps%GetOrbit(d)
    r = 0.d0
    if((-1)**(oa%l+ob%l+oc%l+od%l) == -1) return
    if(oa%z+ob%z /= oc%z+od%z) return
    Kmin = max(abs(oa%j-ob%j), abs(oc%j-od%j))/2
    Kmax = min(    oa%j+ob%j ,     oc%j+od%j )/2
    do K = Kmin, Kmax
      if(a==b .and. mod(K,2)==1) cycle
      if(c==d .and. mod(K,2)==1) cycle
      r = r + dble(2*K+1) * &
          & sjs(oa%j, ob%j, 2*K, oc%j, od%j, 2*J) * &
          & this%GetTwBME(a,b,c,d,K)
    end do
  end function GetCrossCoupledME1423

  ! Only for Scalar
  function GetCrossCoupledME1324(this, a, b, c, d, J) result(r)
    use MyLibrary, only: sjs
    class(TwoBodyPart), intent(in) :: this
    integer, intent(in) :: a, b, c, d, J
    real(8) :: r
    type(Orbits), pointer :: sps
    type(SingleParticleOrbit), pointer :: oa, ob, oc, od
    integer :: K, Kmin, Kmax
    sps => this%two%sps
    oa => sps%GetOrbit(a)
    ob => sps%GetOrbit(b)
    oc => sps%GetOrbit(c)
    od => sps%GetOrbit(d)
    r = 0.d0
    if((-1)**(oa%l+ob%l+oc%l+od%l) == -1) return
    if(oa%z+ob%z /= oc%z+od%z) return
    Kmin = max(abs(oa%j-ob%j), abs(oc%j-od%j))/2
    Kmax = min(    oa%j+ob%j ,     oc%j+od%j )/2
    do K = Kmin, Kmax
      if(a==b .and. mod(K,2)==1) cycle
      if(c==d .and. mod(K,2)==1) cycle
      r = r + dble(2*K+1) * &
          & sjs(oa%j, ob%j, 2*K, od%j, oc%j, 2*J) * &
          & this%GetTwBME(a,b,c,d,K)
    end do
  end function GetCrossCoupledME1324

  function NormalOrderingFrom2To1(this, one) result(r)
    use MyLibrary, only: sjs, triag
    class(TwoBodyPart), intent(in) :: this
    type(OneBodySpace), intent(in) :: one
    type(OneBodyPart) :: r

    if(this%Scalar) r = ScalarNormalOrderingFrom2To1(this, one)
    if(.not. this%Scalar) r = TensorNormalOrderingFrom2To1(this, one)
  end function NormalOrderingFrom2To1

  function ScalarNormalOrderingFrom2To1(this, one) result(r)
    use MyLibrary, only: sjs, triag
    use Profiler, only: timer
    type(TwoBodyPart), intent(in) :: this
    type(OneBodySpace), intent(in) :: one
    type(OneBodyPart) :: r
    integer :: bra, ket, ih, J, ch, ibra, iket
    type(SingleParticleOrbit), pointer :: o1, o2, oh
    real(8) :: vsum, fact, v, ti

    call r%init(one,this%Scalar,this%oprtr,this%jr,this%pr,this%zr)

    ti = omp_get_wtime()
    !$omp parallel
    !$omp do private(bra, o1, ket, o2, vsum, ih, oh, fact, v, &
    !$omp &          J, ch, ibra, iket)
    do bra = 1, one%sps%norbs
      o1 => one%sps%GetOrbit(bra)
      do ket = 1, bra
        o2 => one%sps%GetOrbit(ket)

        if(o1%z /= o2%z) cycle
        if(o1%l /= o2%l) cycle
        if(o1%j /= o2%j) cycle
        vsum = 0.d0
        do ih = 1, one%sps%norbs
          oh => one%sps%GetOrbit(ih)
          if(oh%occ < 1.d-8) cycle
          if(o1%e + oh%e > this%two%e2max) cycle
          if(o1%e + oh%e > this%two%e2max) cycle
          fact = 1.d0
          if(bra==ih) fact = fact * dsqrt(2.d0)
          if(ket==ih) fact = fact * dsqrt(2.d0)
          v = 0.d0
          do J = abs(o1%j-oh%j)/2, (o1%j+oh%j)/2
            if(bra == ih .and. mod(J,2)==1) cycle
            if(ket == ih .and. mod(J,2)==1) cycle
            v = v + dble(2*J+1) * &
                & this%GetTwBME(bra,ih,ket,ih,J) * &
                & oh%occ
          end do
          vsum = vsum + v * fact
        end do
        ch = one%jpz2ch(o1%j, (-1)**o1%l, o1%z)
        ibra = one%jpz(ch)%spi2n(bra)
        iket = one%jpz(ch)%spi2n(ket)
        r%MatCh(ch,ch)%m(ibra,iket) = vsum / dble(o1%j+1)
        r%MatCh(ch,ch)%m(iket,ibra) = vsum / dble(o1%j+1)
      end do
    end do
    !$omp end do
    !$omp end parallel
    call timer%Add("ScalarOrderingFrom2To1", omp_get_wtime()-ti)
  end function ScalarNormalOrderingFrom2To1

  function TensorNormalOrderingFrom2To1(this, one) result(r)
    use MyLibrary, only: sjs, triag
    use Profiler, only: timer
    type(TwoBodyPart), intent(in) :: this
    type(OneBodySpace), intent(in) :: one
    type(OneBodyPart) :: r
    !type(SingleParticleOrbit), pointer :: o1, o2, oh
    !integer :: i1, i2, ih, j1, j2, l1, l2, z1, z2, e1, e2
    !integer :: jh, eh, J_bra, J_ket
    !integer :: chbra, chket, bra, ket, bra_max, ket_max
    !real(8) :: vsum, v, fact, tfact
    real(8) :: ti

    call r%init(one,this%Scalar,this%oprtr,this%jr,this%pr,this%zr)
    return
    ti = omp_get_wtime()
    !do chbra = 1, one%NChan
    !  do chket = 1, chbra
    !    if( .not. r%MatCh(chbra,chket)%is) cycle


    !    bra_max = one%jpz(chbra)%n_state
    !    ket_max = one%jpz(chket)%n_state

    !    j1 = one%jpz(chbra)%j
    !    j2 = one%jpz(chket)%j
    !    !$omp parallel
    !    !$omp do private(bra, i1, o1, l1, z1, e1, ket_max, ket, i2, o2, l2, z2, e2, vsum, &
    !    !$omp &          ih, oh, jh, eh, fact, v, J_bra, J_ket, tfact)
    !    do bra = 1, bra_max
    !      i1 = one%jpz(chbra)%n2spi(bra)
    !      o1 => this%two%sps%GetOrbit(i1)
    !      l1 = o1%l
    !      z1 = o1%z
    !      e1 = o1%e
    !      if(chbra == chket) ket_max = bra
    !      do ket = 1, ket_max
    !        i2 = one%jpz(chket)%n2spi(ket)
    !        o2 => this%two%sps%GetOrbit(i2)
    !        l2 = o2%l
    !        z2 = o2%z
    !        e2 = o2%e

    !        vsum = 0.d0
    !        do ih = 1, this%two%sps%norbs
    !          oh => this%two%sps%GetOrbit(ih)
    !          if(oh%occ < 1.d-8) cycle
    !          jh = oh%j
    !          eh = oh%e
    !          if(e1 + eh > this%two%e2max) cycle
    !          if(e2 + eh > this%two%e2max) cycle
    !          fact = 1.d0
    !          if(i1==ih) fact = fact * dsqrt(2.d0)
    !          if(i2==ih) fact = fact * dsqrt(2.d0)

    !          v = 0.d0
    !          do J_bra = abs(j1-jh)/2, (j1+jh)/2
    !            if(i1 == ih .and. mod(J_bra,2)==1) cycle
    !            do J_ket = abs(j2-jh)/2, (j2+jh)/2
    !              if(i2 == ih .and. mod(J_ket,2)==1) cycle

    !              if(triag(J_bra,J_ket,this%jr)) cycle
    !              tfact = sqrt(dble( (2*J_bra+1) * (2*J_ket+1) ) )
    !              if(.not. this%Scalar) then
    !                tfact = sqrt(dble( (2*J_bra+1) * (2*J_ket+1) ) ) * &
    !                    (-1.d0) ** ((j1+jh)/2+J_ket+this%jr) * &
    !                    & sjs(j1,2*J_bra,jh,2*J_ket,j2,2*this%jr)
    !              end if
    !              v = v + tfact * &
    !                  & this%GetTwBME(i1,ih,i2,ih,J_bra,J_ket) * &
    !                  & oh%occ
    !            end do
    !          end do
    !          vsum = vsum + v * fact
    !        end do
    !        r%MatCh(chbra,chket)%m(bra,ket) = vsum
    !        if(this%Scalar) then
    !          r%MatCh(chbra,chket)%m(bra,ket) = vsum / dble(j1+1)
    !          r%MatCh(chbra,chket)%m(ket,bra) = vsum / dble(j1+1)
    !        end if
    !      end do
    !    end do
    !    !$omp end do
    !    !$omp end parallel
    !  end do
    !end do
    call timer%Add("TensorNormalOrderingFrom2To1", omp_get_wtime()-ti)
  end function TensorNormalOrderingFrom2To1


  !
  !
  !     File reading methods
  !
  !
  subroutine Set2BodyReadFile(this, file_nn)
    use ClassSys, only: sys
    class(Read2BodyFiles), intent(inout) :: this
    character(*), intent(in), optional :: file_nn
    type(sys) :: s
    logical :: ex

    ! -- two-body file
    this%file_nn = 'none'
    if(present(file_nn)) this%file_nn = file_nn
    select case(this%file_nn)
    case("NONE", "none", "None")
    case default
      ex = s%isfile(this%file_nn, "SetReadFiles: two-body file")
    end select
  end subroutine Set2BodyReadFile

  subroutine Set2BodyFileBoundaries(this, emax, e2max, lmax)
    class(Read2BodyFiles), intent(inout) :: this
    integer, intent(in) :: emax, e2max, lmax

    this%emax2 = emax
    this%e2max2= e2max
    this%lmax2 = lmax

  end subroutine Set2BodyFileBoundaries

  subroutine ReadTwoBodyFile(this, V)
    use Profiler, only: timer
    class(Read2BodyFiles), intent(in) :: this
    class(TwoBodyPart), intent(inout) :: v
    real(8) :: ti

    ti = omp_get_wtime()

    select case(this%file_nn)
    case('None', 'NONE', 'none')
      write(*,*) "Error in ReadTwoBodyFile"

    case default

      if(V%Scalar) call this%ReadScalar2BFile(V)
      if(.not. V%Scalar) call this%ReadTensor2BFile(V)
    end select

    call timer%Add("Read from file", omp_get_wtime()-ti)
  end subroutine ReadTwoBodyFile

  subroutine ReadScalar2BFile(this,two)
    use ClassSys, only: sys
    class(Read2BodyFiles), intent(in) :: this
    type(TwoBodyPart), intent(inout) :: two
    type(sys) :: s


    if(s%find(this%file_nn, '.me2j.bin')) then
      call this%read_scalar_me2j_bin(two)
      return
    end if

    if(s%find(this%file_nn, '.me2j')) then
      call this%read_scalar_me2j_ascii(two)
      return
    end if

    if(s%find(this%file_nn, '.myg.bin')) then
      call this%read_scalar_myg_bin(two)
      return
    end if

    if(s%find(this%file_nn, '.myg')) then
      call this%read_scalar_myg_ascii(two)
      return
    end if

    if(s%find(this%file_nn, '.snt')) then
      call this%read_scalar_snt_ascii(two)
      return
    end if

    if(s%find(this%file_nn, 'TBMEA2') .and. s%find(this%file_nn, '.gz')) then
      call this%read_scalar_navratil_ascii_gz(two)
      return
    end if

    if(s%find(this%file_nn, 'TBMEA2')) then
      call this%read_scalar_navratil_ascii(two)
      return
    end if
    write(*,*) "In ReadScalar2BFile, file format cannot be detected."
  end subroutine ReadScalar2BFile

  subroutine read_scalar_me2j_ascii(this,two)
    use ClassSys, only: sys
    class(Read2BodyFiles), intent(in) :: this
    type(TwoBodyPart), intent(inout) :: two
    type(TwoBodySpace), pointer :: ms
    type(Orbits), pointer :: sps
    type(OrbitsIsospin) :: sps_me2j
    integer :: io, runit = 22
    real(8), allocatable :: v(:)
    integer :: a, b, c, d, dmax
    integer :: la, ja, ea
    integer :: lb, jb, eb
    integer :: lc, jc, ec
    integer :: ld, jd, ed
    integer :: J, Jmin, Jmax
    integer(8) :: nelm, icnt
    integer :: ap, bp, cp, dp
    integer :: an, bn, cn, dn
    real(8) :: me_00, me_pp, me_10, me_nn, fact
    type(sys) :: s

    ms => two%two
    sps => ms%sps
    call sps_me2j%init(this%emax2, this%lmax2)
    nelm = count_scalar_me2j(sps_me2j,this%e2max2)
    allocate(v(nelm))
    if( s%find(this%file_nn, '.gz') ) then
      call get_vector_me2j_gz(this%file_nn,v)
    end if

    if( .not. s%find(this%file_nn, '.gz') ) then
      open(runit, file=this%file_nn, action='read',iostat=io)
      if(io /= 0) then
        write(*,'(2a)') 'File open error: ', trim(this%file_nn)
        return
      end if
      call get_vector_me2j_formatted(runit,v)
      close(runit)
    end if

    icnt = 0
    do a = 1, sps_me2j%norbs
      la = sps_me2j%orb(a)%l
      ja = sps_me2j%orb(a)%j
      ea = sps_me2j%orb(a)%e
      if(ea > ms%emax) cycle
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

              if(a == b .and. mod(J,2) == 1 .and. abs(me_nn) > 1.d-8) call me2j_read_warning(a,b,c,d,J,me_nn)
              if(c == d .and. mod(J,2) == 1 .and. abs(me_nn) > 1.d-8) call me2j_read_warning(a,b,c,d,J,me_nn)
              if(a == b .and. mod(J,2) == 1 .and. abs(me_pp) > 1.d-8) call me2j_read_warning(a,b,c,d,J,me_pp)
              if(c == d .and. mod(J,2) == 1 .and. abs(me_pp) > 1.d-8) call me2j_read_warning(a,b,c,d,J,me_pp)

              if(ap == 0) cycle
              if(bp == 0) cycle
              if(cp == 0) cycle
              if(dp == 0) cycle

              if(sps%orb(ap)%e + sps%orb(bp)%e > ms%e2max) cycle
              if(sps%orb(cp)%e + sps%orb(dp)%e > ms%e2max) cycle
#ifdef TwoBodyOperatorDebug
              write(*,'(5i3,4f12.6)') a, b, c, d, J, &
                  &  me_00, me_nn, me_10, me_pp
#endif
              fact = 1.d0
              if(a == b) fact = fact / dsqrt(2.d0)
              if(c == d) fact = fact / dsqrt(2.d0)

              if(abs(1.d0-fact) < 1.d-4 .or. mod(J,2) == 0) then
                call two%SetTwBME(ap,bp,cp,dp,J,me_pp*fact)
                call two%SetTwBME(an,bn,cn,dn,J,me_nn*fact)
                call two%AddToTwBME(ap,bn,cp,dn,J,0.5d0*me_10) ! pnpn
                if(c/=d) call two%AddToTwBME(ap,bn,cn,dp,J,0.5d0*me_10) ! pnnp
                if(a/=b .and. c/=d) &
                    &    call two%AddToTwBME(an,bp,cn,dp,J,0.5d0*me_10) ! npnp
                if(a/=b .and. (a/=c .or. b/=d)) &
                    &    call two%AddToTwBME(an,bp,cp,dn,J,0.5d0*me_10) ! nppn
              end if

              if(abs(1.d0-fact) < 1.d-4 .or. mod(J,2) == 1) then
                call two%AddToTwBME(ap,bn,cp,dn,J,0.5d0*me_00) ! pnpn
                if(c/=d) call two%AddToTwBME(ap,bn,cn,dp,J,-0.5d0*me_00) ! pnnp
                if(a/=b .and. c/=d) &
                    &    call two%AddToTwBME(an,bp,cn,dp,J, 0.5d0*me_00) ! npnp
                if(a/=b .and. (a/=c .or. b/=d)) &
                    &    call two%AddToTwBME(an,bp,cp,dn,J,-0.5d0*me_00) ! nppn
              end if
            end do
          end do
        end do
      end do
    end do

    deallocate(v)
    call sps_me2j%fin()
    return
  end subroutine read_scalar_me2j_ascii

  subroutine me2j_read_warning(a,b,c,d,J,me)
    integer, intent(in) :: a, b, c, d, J
    real(8), intent(in) :: me
    write(*,'(a,5i3,f12.6)') "Warning: this TBME should be zero: ", a, b, c, d, J, me
  end subroutine me2j_read_warning

  subroutine read_scalar_me2j_bin(this,two)
    class(Read2BodyFiles), intent(in) :: this
    type(TwoBodyPart), intent(inout) :: two
    type(TwoBodySpace), pointer :: ms
    type(Orbits), pointer :: sps
    type(OrbitsIsospin) :: sps_me2j
    integer :: io, runit = 22
#ifdef single_precision_two_body_file
    real(4), allocatable :: v(:)
#else
    real(8), allocatable :: v(:)
#endif
    integer :: a, b, c, d, dmax
    integer :: la, ja, ea
    integer :: lb, jb, eb
    integer :: lc, jc, ec
    integer :: ld, jd, ed
    integer :: J, Jmin, Jmax
    integer(8) :: nelm, icnt
    integer :: ap, bp, cp, dp
    integer :: an, bn, cn, dn
    real(8) :: me_00, me_pp, me_10, me_nn, fact

    ms => two%two
    sps => ms%sps
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
      if(ea > ms%emax) cycle
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

              if(a == b .and. mod(J,2) == 1 .and. abs(me_nn) > 1.d-8) call me2j_read_warning(a,b,c,d,J,me_nn)
              if(c == d .and. mod(J,2) == 1 .and. abs(me_nn) > 1.d-8) call me2j_read_warning(a,b,c,d,J,me_nn)
              if(a == b .and. mod(J,2) == 1 .and. abs(me_pp) > 1.d-8) call me2j_read_warning(a,b,c,d,J,me_nn)
              if(c == d .and. mod(J,2) == 1 .and. abs(me_pp) > 1.d-8) call me2j_read_warning(a,b,c,d,J,me_nn)

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
                call two%SetTwBME(ap,bp,cp,dp,J,me_pp*fact)
                call two%SetTwBME(an,bn,cn,dn,J,me_nn*fact)
                call two%AddToTwBME(ap,bn,cp,dn,J,0.5d0*me_10)
                if(c/=d) call two%AddToTwBME(ap,bn,cn,dp,J,0.5d0*me_10)
                if(a/=b .and. c/=d) &
                    &    call two%AddToTwBME(an,bp,cn,dp,J,0.5d0*me_10)
                if(a/=b .and. (a/=c .or. b/=d)) &
                    &    call two%AddToTwBME(an,bp,cp,dn,J,0.5d0*me_10)
              end if

              if(abs(1.d0-fact) < 1.d-4 .or. mod(J,2) == 1) then
                call two%AddToTwBME(ap,bn,cp,dn,J,0.5d0*me_00)
                if(c/=d) call two%AddToTwBME(ap,bn,cn,dp,J,-0.5d0*me_00)
                if(a/=b .and. c/=d) &
                    &    call two%AddToTwBME(an,bp,cn,dp,J, 0.5d0*me_00)
                if(a/=b .and. (a/=c .or. b/=d)) &
                    &    call two%AddToTwBME(an,bp,cp,dn,J,-0.5d0*me_00)
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
    integer(8) :: r
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
    use ClassSys, only: sys
    integer, intent(in) :: ut
    real(8), intent(inout) :: v(:)
    integer(8) :: nelm, lines, i
    integer :: nelm_tail
#ifdef TwoBodyOperatorDebug
    type(sys) :: s
    character(20) :: fm
#endif
    nelm = size(v)
    lines = nelm/10
    nelm_tail = nelm - 10 * lines
    read(ut,*) ! header
    do i = 1, lines
      read(ut,*) v( (i-1)*10+1 : i*10)
#ifdef TwoBodyOperatorDebug
      write(*,'(10f12.6)') v((i-1)*10+1:i*10)
#endif
    end do
    if(mod(nelm,10) == 0) return
    read(ut,*) v(lines*10+1:nelm)
#ifdef TwoBodyOperatorDebug
    fm = "("//trim(s%str(nelm_tail)) // 'f12.6'//")"
    write(*,fm) v(lines*10+1:)
#endif
  end subroutine get_vector_me2j_formatted

  subroutine get_vector_me2j_gz(f,v)
    use, intrinsic :: iso_c_binding
    use ClassSys, only: sys
    use MyLibrary, only: gzip_open, gzip_readline, gzip_close
    character(*), intent(in) :: f
    real(8), intent(inout) :: v(:)
    integer(8) :: nelm, line, lines, cnt
    integer :: nelm_tail
    character(512) :: buffer = ''
    type(c_ptr) :: p, buf
#ifdef TwoBodyOperatorDebug
    type(sys) :: s
    character(20) :: fm
#endif
    nelm = size(v)
    lines = nelm/10
    nelm_tail = nelm - 10 * lines
    p = gzip_open(f,'r')
    buf = gzip_readline(p, buffer, len(buffer))
    buffer=""
    cnt = 0
    do line = 1, lines
      buf = gzip_readline(p, buffer, len(buffer))
      read(buffer,*) v(cnt+1:cnt+10)
#ifdef TwoBodyOperatorDebug
      write(*,'(10f12.6)') v(cnt+1:cnt+10)
#endif
      cnt = cnt + 10
    end do
    if(nelm_tail == 0) then
      buf = gzip_close(p)
      return
    end if
    buffer=""
    buf = gzip_readline(p, buffer, len(buffer))
    read(buffer,*) v(lines*10+1:nelm)
#ifdef TwoBodyOperatorDebug
    fm = "("//trim(s%str(nelm_tail)) // 'f12.6'//")"
    write(*,fm) v(lines*10+1:)
#endif
    buf = gzip_close(p)
  end subroutine get_vector_me2j_gz

  subroutine read_scalar_myg_ascii(this,two)
    class(Read2BodyFiles), intent(in) :: this
    type(TwoBodyPart), intent(inout) :: two
    type(TwoBodySpace), pointer :: ms
    type(Orbits), pointer :: sps
    integer :: runit = 22, io
    integer(8) :: n, lines
    integer :: a, b, c, d, J
    real(8) :: me

    ms => two%two
    sps => ms%sps
    open(runit,file=this%file_nn, action='read', iostat=io)
    if(io /= 0) then
      write(*,'(2a)') 'File open error: ', trim(this%file_nn)
      return
    end if

    read(runit,*) lines ! header
    read(runit,*)       ! header
    do n = 1, lines
      read(runit,*) a, b, c, d, J, me
#ifdef TwoBodyOperatorDebug
      write(*,'(5i3,f12.6)') a, b, c, d, J, me
#endif
      if(a > sps%norbs) cycle
      if(b > sps%norbs) cycle
      if(c > sps%norbs) cycle
      if(d > sps%norbs) cycle
      call two%SetTwBME(a,b,c,d,J,me)
    end do
    close(runit)
  end subroutine read_scalar_myg_ascii

  subroutine read_scalar_myg_bin(this,two)
    class(Read2BodyFiles), intent(in) :: this
    type(TwoBodyPart), intent(inout) :: two
    type(TwoBodySpace), pointer :: ms
    type(Orbits), pointer :: sps
    integer :: runit = 22, io
    integer(8) :: lines, n
    integer, allocatable :: a(:), b(:), c(:), d(:), J(:)
    real(8), allocatable :: me(:)

    ms => two%two
    sps => ms%sps
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
      call two%SetTwBME(a(n),b(n),c(n),d(n),J(n),me(n))
    end do
    !$omp end do
    !$omp end parallel
  end subroutine read_scalar_myg_bin

  subroutine read_scalar_snt_ascii(this,two)
    use MyLibrary, only: skip_comment
    class(Read2BodyFiles), intent(in) :: this
    type(TwoBodyPart), intent(inout) :: two
    type(TwoBodySpace), pointer :: ms
    type(Orbits), pointer :: sps
    type(SingleParticleOrbit), allocatable :: ssnt(:)
    integer :: a, b, c, d, J, idx, n, l, z
    integer :: ia, ib, ic, id
    integer :: runit=22, io
    integer :: porbs, norbs, pc, nc
    integer(8) :: lines, line
    real(8) :: me

    ms => two%two
    sps => ms%sps
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
#ifdef TwoBodyOperatorDebug
      write(*,'(5i3,f12.6)') a, b, c, d, J, me
#endif
      ia = sps%nljz2idx(ssnt(a)%n, ssnt(a)%l, ssnt(a)%j, ssnt(a)%z)
      ib = sps%nljz2idx(ssnt(b)%n, ssnt(b)%l, ssnt(b)%j, ssnt(b)%z)
      ic = sps%nljz2idx(ssnt(c)%n, ssnt(c)%l, ssnt(c)%j, ssnt(c)%z)
      id = sps%nljz2idx(ssnt(d)%n, ssnt(d)%l, ssnt(d)%j, ssnt(d)%z)
      call two%SetTwBME(ia,ib,ic,id,J,me)
    end do
    close(runit)
  end subroutine read_scalar_snt_ascii

  subroutine read_scalar_navratil_ascii_gz(this,two)
    use, intrinsic :: iso_c_binding
    use MyLibrary, only: gzip_open, gzip_readline, gzip_close
    class(Read2BodyFiles), intent(in) :: this
    type(TwoBodyPart), intent(inout) :: two
    type(TwoBodySpace), pointer :: ms
    type(Orbits), pointer :: sps
    type(OrbitsIsospin) :: isps
    integer :: a, b, c, d, J, T
    integer :: ap, bp, cp, dp
    integer :: an, bn, cn, dn
    integer(8) :: iline, nlines
    real(8) :: trel, horel, vcoul, vpn, vpp, vnn, fact
    integer :: emax_in, e2max_in
    real(8) :: hw_in, lambda_in
    type(c_ptr) :: p, buf
    character(512) :: buffer = ''


    ms => two%two
    sps => ms%sps
    call isps%init(this%emax2, this%lmax2)

    p = gzip_open(this%file_nn, 'r')
    buf = gzip_readline(p, buffer, len(buffer))
    read(buffer,*) nlines, emax_in, e2max_in, hw_in, lambda_in
    if(emax_in /= this%emax2 ) write(*,'(a)') "Warning: file emax doesn't match"
    if(e2max_in/= this%e2max2) write(*,'(a)') "Warning: file e2max doesn't match"
    buffer=""
    do iline = 1, nlines
      buf = gzip_readline(p, buffer, len(buffer))
      read(buffer,*) a, b, c, d, J, T, trel, horel, vcoul, vpn, vpp, vnn
      if(J > min(2*sps%lmax+1, ms%e2max+1)) exit
      if(isps%orb(a)%e > ms%emax) cycle
      if(isps%orb(b)%e > ms%emax) cycle
      if(isps%orb(c)%e > ms%emax) cycle
      if(isps%orb(d)%e > ms%emax) cycle
      if(isps%orb(a)%e + isps%orb(b)%e > ms%e2max) cycle
      if(isps%orb(c)%e + isps%orb(d)%e > ms%e2max) cycle
      ap = isps%iso2pn(sps,a,-1)
      an = isps%iso2pn(sps,a, 1)
      bp = isps%iso2pn(sps,b,-1)
      bn = isps%iso2pn(sps,b, 1)
      cp = isps%iso2pn(sps,c,-1)
      cn = isps%iso2pn(sps,c, 1)
      dp = isps%iso2pn(sps,d,-1)
      dn = isps%iso2pn(sps,d, 1)
      fact = 1.d0
      if(a/=b) fact = fact / dsqrt(2.d0)
      if(c/=d) fact = fact / dsqrt(2.d0)

      if(T==1) then
        call two%SetTwBME(ap,bp,cp,dp,J,vpp)
        call two%SetTwBME(an,bn,cn,dn,J,vnn)
        call two%AddToTwBME(ap,bn,cp,dn,J,fact*vpn)
        if(c/=d) call two%AddToTwBME(ap,bn,cn,dp,J,fact*vpn)
        if(a/=b .and. c/=d) &
            &    call two%AddToTwBME(an,bp,cn,dp,J,fact*vpn)
        if(a/=b .and. (a/=c .or. b/=d)) &
            &    call two%AddToTwBME(an,bp,cp,dn,J,fact*vpn)
      end if

      if(T==0) then
        call two%AddToTwBME(ap,bn,cp,dn,J,fact*vpn)
        if(c/=d) call two%AddToTwBME(ap,bn,cn,dp,J,-fact*vpn)
        if(a/=b .and. c/=d) &
            &    call two%AddToTwBME(an,bp,cn,dp,J,fact*vpn)
        if(a/=b .and. (a/=c .or. b/=d)) &
            &    call two%AddToTwBME(an,bp,cp,dn,J,-fact*vpn)
      end if
    end do
    buf = gzip_close(p)
  end subroutine read_scalar_navratil_ascii_gz

  subroutine read_scalar_navratil_ascii(this,two)
    class(Read2BodyFiles), intent(in) :: this
    type(TwoBodyPart), intent(inout) :: two
    type(TwoBodySpace), pointer :: ms
    type(Orbits), pointer :: sps
    type(OrbitsIsospin) :: isps
    integer :: io, runit = 22
    integer :: a, b, c, d, J, T
    integer :: ap, bp, cp, dp
    integer :: an, bn, cn, dn
    integer(8) :: iline, nlines
    real(8) :: trel, horel, vcoul, vpn, vpp, vnn, fact
    integer :: emax_in, e2max_in
    real(8) :: hw_in, lambda_in

    ms => two%two
    sps => ms%sps
    call isps%init(this%emax2, this%lmax2)
    open(runit, file=this%file_nn, action='read',iostat=io)
    read(runit, *) nlines, emax_in, e2max_in, hw_in, lambda_in
    if(emax_in /= this%emax2 ) write(*,'(a)') "Warning: file emax doesn't match"
    if(e2max_in/= this%e2max2) write(*,'(a)') "Warning: file e2max doesn't match"
    do iline = 1, nlines
      read(runit,*) a, b, c, d, J, T, trel, horel, vcoul, vpn, vpp, vnn
      if(J > min(2*sps%lmax+1, ms%e2max+1)) exit
      if(isps%orb(a)%e > ms%emax) cycle
      if(isps%orb(b)%e > ms%emax) cycle
      if(isps%orb(c)%e > ms%emax) cycle
      if(isps%orb(d)%e > ms%emax) cycle
      if(isps%orb(a)%e + isps%orb(b)%e > ms%e2max) cycle
      if(isps%orb(c)%e + isps%orb(d)%e > ms%e2max) cycle
      ap = isps%iso2pn(sps,a,-1)
      an = isps%iso2pn(sps,a, 1)
      bp = isps%iso2pn(sps,b,-1)
      bn = isps%iso2pn(sps,b, 1)
      cp = isps%iso2pn(sps,c,-1)
      cn = isps%iso2pn(sps,c, 1)
      dp = isps%iso2pn(sps,d,-1)
      dn = isps%iso2pn(sps,d, 1)
      fact = 1.d0
      if(a/=b) fact = fact / dsqrt(2.d0)
      if(c/=d) fact = fact / dsqrt(2.d0)

      if(T==1) then
        call two%SetTwBME(ap,bp,cp,dp,J,vpp)
        call two%SetTwBME(an,bn,cn,dn,J,vnn)
        call two%AddToTwBME(ap,bn,cp,dn,J,fact*vpn)
        if(c/=d) call two%AddToTwBME(ap,bn,cn,dp,J,fact*vpn)
        if(a/=b .and. c/=d) &
            &    call two%AddToTwBME(an,bp,cn,dp,J,fact*vpn)
        if(a/=b .and. (a/=c .or. b/=d)) &
            &    call two%AddToTwBME(an,bp,cp,dn,J,fact*vpn)
      end if

      if(T==0) then
        call two%AddToTwBME(ap,bn,cp,dn,J,fact*vpn)
        if(c/=d) call two%AddToTwBME(ap,bn,cn,dp,J,-fact*vpn)
        if(a/=b .and. c/=d) &
            &    call two%AddToTwBME(an,bp,cn,dp,J,fact*vpn)
        if(a/=b .and. (a/=c .or. b/=d)) &
            &    call two%AddToTwBME(an,bp,cp,dn,J,-fact*vpn)
      end if
    end do
    close(runit)
  end subroutine read_scalar_navratil_ascii

  subroutine ReadTensor2BFile(this,two)
    use ClassSys, only: sys
    class(Read2BodyFiles), intent(in) :: this
    type(TwoBodyPart), intent(inout) :: two
    return
  end subroutine ReadTensor2BFile

end module TwoBodyOperator

!program test
!  use TwoBodyOperator
!  type(Orbits), target :: sps
!  type(TwoBodySpace) :: two
!  type(TwoBodyPart) :: op
!  integer :: emax = 2
!
!  call sps%init(emax)
!  do i = 1, sps%norbs
!    call sps%orb(i)%SetHoleParticleValence(1)
!  end do
!  call two%init(sps, 2*emax)
!  call op%init(two, .true., 'Hcm', 0, 1, 0)
!  call op%set(sps, 20.d0, 4, 2, 2)
!
!  call op%fin()
!  call two%fin()
!  call sps%fin()
!end program test
