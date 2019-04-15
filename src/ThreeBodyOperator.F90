module ThreeBodyOperator
  use omp_lib
  use LinAlgLib
  use ThreeBodyModelSpace
  use OneBodyOperator
  use TwoBodyOperator
  implicit none
  public :: ThreeBodyPartChannel
  public :: ThreeBodyPart

  private :: InitThreeBodyPartChannel
  private :: FinThreeBodyPartChannel
  private :: SetThreeBodyPartChannelFromOneBody
  private :: SetThreeBodyPartChannelFromTwoBody
  private :: SetThreeBodyPartChannel
  private :: InitThreeBodyPart
  private :: FinThreeBodyPart
  private :: SetThreeBodyPartFromOneBody
  private :: SetThreeBodyPartFromTwoBody
  private :: SetThreeBodyPart
  private :: CopyThreeBodyPart
  private :: SumThreeBodyPart
  private :: SubtractThreeBodyPart
  private :: ScaleThreeBodyPart
  private :: GetThBME_i_scalar
  private :: GetThBME_i_tensor
  private :: GetThBME_J_scalar
  private :: GetThBME_J_tensor
  private :: PrintThreeBodyPart
  private :: NormalOrderingFrom3To2

  type, extends(DMat) :: ThreeBodyPartChannel
    type(ThreeBodyChannel), pointer :: ch_bra, ch_ket
    logical :: is = .false.
  contains
    procedure :: InitThreeBodyPartChannel
    procedure :: FinThreeBodyPartChannel
    procedure :: SetThreeBodyPartChannelFromOneBody
    procedure :: SetThreeBodyPartChannelFromTwoBody
    procedure :: SetThreeBodyPartChannel

    generic :: init => InitThreeBodyPartChannel
    generic :: release => FinThreeBodyPartChannel
    generic :: set => SetThreeBodyPartChannelFromOneBody, &
        & SetThreeBodyPartChannelFromTwoBody, &
        & SetThreeBodyPartChannel
  end type ThreeBodyPartChannel

  type :: ThreeBodyPart
    type(ThreeBodyPartChannel), allocatable :: MatCh(:,:)
    type(ThreeBodySpace), pointer :: thr
    character(:), allocatable :: oprtr
    logical :: Scalar
    integer :: jr, pr, zr
  contains
    procedure :: InitThreeBodyPart
    procedure :: FinThreeBodyPart
    procedure :: SetThreeBodyPartFromOneBody
    procedure :: SetThreeBodyPartFromTwoBody
    procedure :: SetThreeBodyPart
    procedure :: CopyThreeBodyPart
    procedure :: SumThreeBodyPart
    procedure :: SubtractThreeBodyPart
    procedure :: ScaleThreeBodyPart
    procedure :: GetThBME_i_scalar
    procedure :: GetThBME_i_tensor
    procedure :: GetThBME_J_scalar
    procedure :: GetThBME_J_tensor
    procedure :: PrintThreeBodyPart
    procedure :: NormalOrderingFrom3To2

    generic :: init => InitThreeBodyPart
    generic :: fin => FinThreeBodyPart
    generic :: set => SetThreeBodyPartFromOneBody, &
        & SetThreeBodyPartFromTwoBody, &
        & SetThreeBodyPart
    generic :: prt => PrintThreeBodyPart
    generic :: GetThBMEi => GetThBME_i_scalar, GetThBME_i_tensor
    generic :: GetThBMEJ => GetThBME_J_scalar, GetThBME_J_tensor
    generic :: assignment(=) => CopyThreeBodyPart
    generic :: operator(+) => SumThreeBodyPart
    generic :: operator(-) => SubtractThreeBodyPart
    generic :: operator(*) => ScaleThreeBodyPart
  end type ThreeBodyPart
contains
  subroutine FinThreeBodyPart(this)
    use MyLibrary, only: triag
    class(ThreeBodyPart), intent(inout) :: this
    integer :: chbra, chket

    do chbra = 1, this%thr%NChan
      do chket = 1, this%thr%NChan
        call this%MatCh(chbra,chket)%release()
      end do
    end do
    deallocate(this%MatCh)
    this%thr => null()

  end subroutine FinThreeBodyPart

  subroutine InitThreeBodyPart(this, thr, Scalar, oprtr, jr, pr, zr)
    use MyLibrary, only: triag
    class(ThreeBodyPart), intent(inout) :: this
    type(ThreeBodySpace), target, intent(in) :: thr
    logical, intent(in) :: Scalar
    character(*), intent(in) :: oprtr
    integer, intent(in) :: jr, pr, zr
    integer :: chbra, chket
    integer :: jbra, pbra, zbra, nbra
    integer :: jket, pket, zket, nket
    this%thr => thr
    this%oprtr = oprtr
    this%Scalar = Scalar
    this%jr = jr
    this%pr = pr
    this%zr = zr

    allocate(this%MatCh(thr%NChan, thr%NChan))
    do chbra = 1, thr%NChan
      jbra = thr%jpz(chbra)%j
      pbra = thr%jpz(chbra)%p
      zbra = thr%jpz(chbra)%z
      nbra = thr%jpz(chbra)%n_state
      do chket = 1, thr%NChan
        jket = thr%jpz(chket)%j
        pket = thr%jpz(chket)%p
        zket = thr%jpz(chket)%z
        nket = thr%jpz(chket)%n_state

        if(triag(jbra, jket, 2*jr)) cycle
        if(pbra * pket * pr /= 1) cycle
        if(zbra - 2*zr - zket /= 0) cycle
        call this%MatCh(chbra,chket)%init(thr%jpz(chbra),thr%jpz(chket))
      end do
    end do

  end subroutine InitThreeBodyPart

  subroutine SetThreeBodyPartFromOneBody(this, op1)
    class(ThreeBodyPart), intent(inout) :: this
    type(OneBodyPart), intent(in) ::op1
    integer :: chbra, chket

    do chbra = 1, this%thr%NChan
      do chket = 1, this%thr%NChan
        if(.not. this%MatCh(chbra,chket)%is) cycle
        call this%MatCh(chbra,chket)%set(op1)
      end do
    end do
  end subroutine SetThreeBodyPartFromOneBody

  function GetThBME_i_scalar(this,a,b,c,ibra,d,e,f,iket,J) result(r)
    use MyLibrary, only: triag
    class(ThreeBodyPart), intent(in) :: this
    integer, intent(in) :: a, b, c, ibra, d, e, f, iket, J
    real(8) :: r
    type(SingleParticleOrbit), pointer :: oa, ob, oc, od, oe, of
    integer :: Pbra, Pket, Zbra, Zket
    integer :: ch, idx_bra, idx_ket, nbra, nket
    r = 0.d0
    oa => this%thr%sps%GetOrbit(a)
    ob => this%thr%sps%GetOrbit(b)
    oc => this%thr%sps%GetOrbit(c)
    od => this%thr%sps%GetOrbit(d)
    oe => this%thr%sps%GetOrbit(e)
    of => this%thr%sps%GetOrbit(f)
    Pbra = (-1)**(oa%l + ob%l + oc%l)
    Pket = (-1)**(od%l + oe%l + of%l)
    Zbra = oa%z + ob%z + oc%z
    Zket = od%z + oe%z + of%z
    if(Pbra * Pket /= 1) return
    if(Zket /= Zbra) return
    ch = this%thr%jpz2ch(J,Pbra,Zbra)
    idx_bra = this%thr%jpz(ch)%spis2idx(a,b,c)
    idx_ket = this%thr%jpz(ch)%spis2idx(d,e,f)
    nbra = this%thr%jpz(ch)%idx(idx_bra)%idx2n(ibra)
    nket = this%thr%jpz(ch)%idx(idx_ket)%idx2n(iket)
    r = this%MatCh(ch,ch)%m(nbra,nket)
  end function GetThBME_i_scalar

  function GetThBME_i_tensor(this,a,b,c,ibra,d,e,f,iket,Jbra,Jket) result(r)
    use MyLibrary, only: triag
    class(ThreeBodyPart), intent(in) :: this
    integer, intent(in) :: a, b, c, ibra, d, e, f, iket, Jbra, Jket
    real(8) :: r
    type(SingleParticleOrbit), pointer :: oa, ob, oc, od, oe, of
    integer :: Pbra, Pket, Zbra, Zket
    integer :: chbra, chket, idx_bra, idx_ket, nbra, nket
    r = 0.d0
    if(triag(Jbra, Jket, 2*this%jr)) return
    oa => this%thr%sps%GetOrbit(a)
    ob => this%thr%sps%GetOrbit(b)
    oc => this%thr%sps%GetOrbit(c)
    od => this%thr%sps%GetOrbit(d)
    oe => this%thr%sps%GetOrbit(e)
    of => this%thr%sps%GetOrbit(f)
    Pbra = (-1)**(oa%l + ob%l + oc%l)
    Pket = (-1)**(od%l + oe%l + of%l)
    Zbra = oa%z + ob%z + oc%z
    Zket = od%z + oe%z + of%z
    if(Pbra * Pket * this%pr /= 1) return
    if(Zket + this%zr /= Zbra) return
    chbra = this%thr%jpz2ch(Jbra,Pbra,Zbra)
    chket = this%thr%jpz2ch(Jket,Pket,Zket)
    idx_bra = this%thr%jpz(chbra)%spis2idx(a,b,c)
    idx_ket = this%thr%jpz(chket)%spis2idx(d,e,f)
    nbra = this%thr%jpz(chbra)%idx(idx_bra)%idx2n(ibra)
    nket = this%thr%jpz(chket)%idx(idx_ket)%idx2n(iket)
    r = this%MatCh(chbra,chket)%m(nbra,nket)
  end function GetThBME_i_tensor

  function GetThBME_J_scalar(this, a, b, c, Jab, d, e, f, Jde, J) result(r)
    class(ThreeBodyPart), intent(in) :: this
    integer, intent(in) :: a, b, c, Jab, d, e, f, Jde, J
    real(8) :: r
    type(Orbits), pointer :: sps
    type(SingleParticleOrbit), pointer :: oa, ob, oc, od, oe, of
    type(AdditionalQN), pointer :: aqn_bra, aqn_ket
    integer :: Pbra, Pket, Zbra, Zket
    integer :: ch, ibra, iket
    integer :: idx_bra, idx_ket
    integer :: n_bra, n_ket
    type(DVec) :: vbra, vket
    type(DMat) :: m
    r = 0.d0
    sps => this%thr%sps
    oa => sps%GetOrbit(a)
    ob => sps%GetOrbit(b)
    oc => sps%GetOrbit(c)
    od => sps%GetOrbit(d)
    oe => sps%GetOrbit(e)
    of => sps%GetOrbit(f)
    Pbra = (-1)**(oa%l + ob%l + oc%l)
    Pket = (-1)**(od%l + oe%l + of%l)
    Zbra = oa%z + ob%z + oc%z
    Zket = od%z + oe%z + of%z
    if(Pbra * Pket /= 1) return
    if(Zket /= Zbra) return
    ch = this%thr%jpz2ch(J,Pbra,Zbra)
    idx_bra = this%thr%jpz(ch)%spis2idx(a, b, c)
    idx_ket = this%thr%jpz(ch)%spis2idx(d, e, f)
    aqn_bra => this%thr%jpz(ch)%idx(idx_bra)
    aqn_ket => this%thr%jpz(ch)%idx(idx_ket)
    n_bra = aqn_bra%find(a,b,c,Jab)
    n_ket = aqn_ket%find(d,e,f,Jde)
    call vbra%ini(aqn_bra%north)
    call vket%ini(aqn_ket%north)
    call m%ini(aqn_bra%north, aqn_ket%north)
    vbra%v = aqn_bra%cfp(:,n_bra)
    vket%v = aqn_ket%cfp(:,n_ket)
    do ibra = 1, aqn_bra%north
      do iket = 1, aqn_ket%north
        m%m(ibra,iket) = this%GetThBME_i_scalar(a,b,c,ibra,d,e,f,iket,J)
      end do
    end do
    r = vbra * (m * vket)
    r = r * 6.d0
    call vbra%fin()
    call vket%fin()
    call m%fin()
  end function GetThBME_J_scalar

  function GetThBME_J_tensor(this, a, b, c, Jab, d, e, f, Jde, Jbra, Jket) result(r)
    use MyLibrary, only: triag
    class(ThreeBodyPart), intent(in) :: this
    integer, intent(in) :: a, b, c, Jab, d, e, f, Jde, Jbra, Jket
    real(8) :: r
    type(Orbits), pointer :: sps
    type(SingleParticleOrbit), pointer :: oa, ob, oc, od, oe, of
    type(AdditionalQN), pointer :: aqn_bra, aqn_ket
    integer :: Pbra, Pket, Zbra, Zket
    integer :: chbra, chket, ibra, iket
    integer :: idx_bra, idx_ket
    integer :: n_bra, n_ket
    type(DVec) :: vbra, vket
    type(DMat) :: m
    r = 0.d0
    if(triag(Jbra, Jket, 2*this%jr)) return
    sps => this%thr%sps
    oa => sps%GetOrbit(a)
    ob => sps%GetOrbit(b)
    oc => sps%GetOrbit(c)
    od => sps%GetOrbit(d)
    oe => sps%GetOrbit(e)
    of => sps%GetOrbit(f)
    Pbra = (-1)**(oa%l + ob%l + oc%l)
    Pket = (-1)**(od%l + oe%l + of%l)
    Zbra = oa%z + ob%z + oc%z
    Zket = od%z + oe%z + of%z
    if(Pbra * Pket * this%pr /= 1) return
    if(Zket + this%zr /= Zbra) return
    chbra = this%thr%jpz2ch(Jbra,Pket,Zket)
    chket = this%thr%jpz2ch(Jket,Pket,Zket)
    idx_bra = this%thr%jpz(chbra)%spis2idx(a, b, c)
    idx_ket = this%thr%jpz(chket)%spis2idx(d, e, f)
    aqn_bra => this%thr%jpz(chbra)%idx(idx_bra)
    aqn_ket => this%thr%jpz(chket)%idx(idx_ket)
    n_bra = aqn_bra%find(a,b,c,Jab)
    n_ket = aqn_ket%find(d,e,f,Jde)
    call vbra%ini(aqn_bra%north)
    call vket%ini(aqn_ket%north)
    call m%ini(aqn_bra%north, aqn_ket%north)
    vbra%v = aqn_bra%cfp(:,n_bra)
    vket%v = aqn_ket%cfp(:,n_ket)
    do ibra = 1, aqn_bra%north
      do iket = 1, aqn_ket%north
        m%m(ibra,iket) = this%GetThBME_i_tensor(a,b,c,ibra,d,e,f,iket,Jbra, Jket)
      end do
    end do
    r = vbra * (m * vket)
    r = r * 6.d0
    call vbra%fin()
    call vket%fin()
    call m%fin()
  end function GetThBME_J_tensor

  subroutine CopyThreeBodyPart(a, b)
    class(ThreeBodyPart), intent(out) :: a
    type(ThreeBodyPart), intent(in) :: b
    integer :: chbra, chket

    if(allocated(a%MatCh)) call a%fin()
    a%thr => b%thr
    a%oprtr = b%oprtr
    a%Scalar = b%Scalar
    a%jr = b%jr
    a%pr = b%pr
    a%zr = b%zr

    allocate(a%MatCh(b%thr%NChan,b%thr%NChan))
    do chbra = 1, b%thr%NChan
      do chket = 1, b%thr%NChan
        if(.not. b%MatCh(chbra,chket)%is) cycle
        a%MatCh(chbra,chket)%is = b%MatCh(chbra,chket)%is
        a%MatCh(chbra,chket)%DMat = b%MatCh(chbra,chket)%DMat
        a%MatCh(chbra,chket)%ch_bra => b%MatCh(chbra,chket)%ch_bra
        a%MatCh(chbra,chket)%ch_ket => b%MatCh(chbra,chket)%ch_ket
      end do
    end do
  end subroutine CopyThreeBodyPart

  function SumThreeBodyPart(a, b) result(c)
    class(ThreeBodyPart), intent(in) :: a, b
    type(ThreeBodyPart) :: c
    integer :: chbra, chket
    c = a
    do chbra = 1, b%thr%NChan
      do chket = 1, b%thr%NChan
        if(.not. b%MatCh(chbra,chket)%is) cycle
        c%MatCh(chbra,chket)%DMat = &
            & a%MatCh(chbra,chket)%DMat + b%MatCh(chbra,chket)%DMat
      end do
    end do
  end function SumThreeBodyPart

  function SubtractThreeBodyPart(a, b) result(c)
    class(ThreeBodyPart), intent(in) :: a, b
    type(ThreeBodyPart) :: c
    integer :: chbra, chket
    c = a
    do chbra = 1, b%thr%NChan
      do chket = 1, b%thr%NChan
        if(.not. b%MatCh(chbra,chket)%is) cycle
        c%MatCh(chbra,chket)%DMat = &
            & a%MatCh(chbra,chket)%DMat - b%MatCh(chbra,chket)%DMat
      end do
    end do
  end function SubtractThreeBodyPart

  function ScaleThreeBodyPart(a, b) result(c)
    class(ThreeBodyPart), intent(in) :: a
    real(8), intent(in) :: b
    type(ThreeBodyPart) :: c
    integer :: chbra, chket
    c = a
    do chbra = 1, a%thr%NChan
      do chket = 1, a%thr%NChan
        if(.not. a%MatCh(chbra,chket)%is) cycle
        c%MatCh(chbra,chket)%DMat = &
            & a%MatCh(chbra,chket)%DMat * b
      end do
    end do
  end function ScaleThreeBodyPart

  subroutine SetThreeBodyPartFromTwoBody(this, op2)
    class(ThreeBodyPart), intent(inout) :: this
    type(TwoBodyPart), intent(in) ::op2
    integer :: chbra, chket

    do chbra = 1, this%thr%NChan
      do chket = 1, this%thr%NChan
        if(.not. this%MatCh(chbra,chket)%is) cycle
        call this%MatCh(chbra,chket)%set(op2)
      end do
    end do
  end subroutine SetThreeBodyPartFromTwoBody

  subroutine SetThreeBodyPart(this, op3_in)
    use ThreeBodyInteraction
    class(ThreeBodyPart), intent(inout) :: this
    type(ThreeBodyForce), intent(in) ::op3_in
    integer :: chbra, chket

    do chbra = 1, this%thr%NChan
      do chket = 1, this%thr%NChan
        if(.not. this%MatCh(chbra,chket)%is) cycle
        call this%MatCh(chbra,chket)%set(op3_in)
      end do
    end do
  end subroutine SetThreeBodyPart

  subroutine FinThreeBodyPartChannel(this)
    class(ThreeBodyPartChannel), intent(inout) :: this
    if(.not. this%is) return
    call this%fin()
    this%is = .false.
    this%ch_bra => null()
    this%ch_ket => null()
  end subroutine FinThreeBodyPartChannel

  subroutine InitThreeBodyPartChannel(this, ch_bra, ch_ket)
    class(ThreeBodyPartChannel), intent(inout) :: this
    type(ThreeBodyChannel), target, intent(in) :: ch_bra, ch_ket
    this%ch_bra => ch_bra
    this%ch_ket => ch_ket
    call this%zeros(ch_bra%n_state, ch_ket%n_state)
    this%is = .true.
  end subroutine InitThreeBodyPartChannel

  subroutine SetThreeBodyPartChannelFromOneBody(this, op1)
    class(ThreeBodyPartChannel), intent(inout) :: this
    type(OneBodyPart), intent(in) ::op1
    integer :: bra, ket
    type(ThreeBodyChannel), pointer :: ch_bra, ch_ket
    type(AdditionalQN), pointer :: qbra, qket
    type(DMat) :: cfp_bra, vin, cfp_ket, vout
    integer :: a, b, c, d, e, f, Jab, Jde
    integer :: bra_sub, ket_sub
    integer :: bra_start, bra_end, ket_start, ket_end

    ch_bra => this%ch_bra
    ch_ket => this%ch_ket

    !$omp parallel
    !$omp do private(bra, qbra, cfp_bra, ket, qket, cfp_ket, &
    !$omp &  vin, bra_sub, a, b, c, Jab, ket_sub, d, e, f, Jde, &
    !$omp &  vout, bra_start, bra_end, ket_start, ket_end)
    do bra = 1, ch_bra%n_idx
      qbra => ch_bra%idx(bra)
      call cfp_bra%ini(qbra%north,qbra%nphys)
      cfp_bra%m = transpose(qbra%cfp)
      do ket = 1, ch_ket%n_idx
        qket => ch_ket%idx(ket)

        call cfp_ket%ini(qket%nphys,qket%north)
        cfp_ket%m = qket%cfp

        call vin%zeros(qbra%nphys,qket%nphys)

        do bra_sub = 1, qbra%nphys
          a = qbra%n2spi1(bra_sub)
          b = qbra%n2spi2(bra_sub)
          c = qbra%n2spi3(bra_sub)
          Jab= qbra%n2J12(bra_sub)
          do ket_sub = 1, qket%nphys
            d = qket%n2spi1(ket_sub)
            e = qket%n2spi2(ket_sub)
            f = qket%n2spi3(ket_sub)
            Jde= qket%n2J12(ket_sub)
            if(c /= f) cycle
            if(b /= e) cycle
            vin%m(bra_sub,ket_sub) =  op1%GetOBME(a,c)

          end do
        end do

        bra_start = qbra%idx2n(1)
        bra_end = qbra%idx2n(qbra%north)
        ket_start = qket%idx2n(1)
        ket_end = qket%idx2n(qket%north)
        vout = cfp_bra * vin * cfp_ket
        vout = 3.d0 * vout
        this%m(bra_start:bra_end, ket_start:ket_end) = vout%m

        call cfp_ket%fin()
        call vout%fin()
        call vin%fin()

      end do
      call cfp_bra%fin()
    end do
    !$omp end do
    !$omp end parallel
  end subroutine SetThreeBodyPartChannelFromOneBody

  subroutine SetThreeBodyPartChannelFromTwoBody(this, op2)
    class(ThreeBodyPartChannel), intent(inout) :: this
    type(TwoBodyPart), intent(in) ::op2
    integer :: bra, ket
    type(ThreeBodyChannel), pointer :: ch_bra, ch_ket
    type(AdditionalQN), pointer :: qbra, qket
    type(DMat) :: cfp_bra, vin, cfp_ket, vout
    integer :: a, b, c, d, e, f, Jab, Jde
    integer :: bra_sub, ket_sub
    integer :: bra_start, bra_end, ket_start, ket_end
    real(8) :: norm

    ch_bra => this%ch_bra
    ch_ket => this%ch_ket

    !$omp parallel
    !$omp do private(bra, qbra, cfp_bra, ket, qket, cfp_ket, &
    !$omp &  vin, bra_sub, a, b, c, Jab, ket_sub, d, e, f, Jde, &
    !$omp &  vout, bra_start, bra_end, ket_start, ket_end)
    do bra = 1, ch_bra%n_idx
      qbra => ch_bra%idx(bra)
      call cfp_bra%ini(qbra%north,qbra%nphys)
      cfp_bra%m = transpose(qbra%cfp)
      do ket = 1, ch_ket%n_idx
        qket => ch_ket%idx(ket)

        call cfp_ket%ini(qket%nphys,qket%north)
        cfp_ket%m = qket%cfp

        call vin%zeros(qbra%nphys,qket%nphys)

        do bra_sub = 1, qbra%nphys
          a = qbra%n2spi1(bra_sub)
          b = qbra%n2spi2(bra_sub)
          c = qbra%n2spi3(bra_sub)
          Jab= qbra%n2J12(bra_sub)
          do ket_sub = 1, qket%nphys
            d = qket%n2spi1(ket_sub)
            e = qket%n2spi2(ket_sub)
            f = qket%n2spi3(ket_sub)
            Jde= qket%n2J12(ket_sub)
            if(c /= f) cycle
            if(Jab /= Jde) cycle
            norm = 1.d0
            if(a == b) norm = norm * dsqrt(2.d0)
            if(d == e) norm = norm * dsqrt(2.d0)
            vin%m(bra_sub,ket_sub) =  op2%GetTwBME(a,b,d,e,Jab) * norm

          end do
        end do

        bra_start = qbra%idx2n(1)
        bra_end = qbra%idx2n(qbra%north)
        ket_start = qket%idx2n(1)
        ket_end = qket%idx2n(qket%north)
        vout = cfp_bra * vin * cfp_ket
        vout = 1.5d0 * vout
        this%m(bra_start:bra_end, ket_start:ket_end) = vout%m

        call cfp_ket%fin()
        call vout%fin()
        call vin%fin()

      end do
      call cfp_bra%fin()
    end do
    !$omp end do
    !$omp end parallel
  end subroutine SetThreeBodyPartChannelFromTwoBody

  subroutine SetThreeBodyPartChannel(this, op3_in)
    use ThreeBodyInteraction
    class(ThreeBodyPartChannel), intent(inout) :: this
    type(ThreeBodyForce), intent(in) ::op3_in
    integer :: bra, ket
    type(ThreeBodyChannel), pointer :: ch_bra, ch_ket
    type(AdditionalQN), pointer :: qbra, qket
    type(DMat) :: cfp_bra, vin, cfp_ket, vout
    integer :: a, b, c, d, e, f, Jab, Jde, J
    integer :: bra_sub, ket_sub
    integer :: bra_start, bra_end, ket_start, ket_end

    ch_bra => this%ch_bra
    ch_ket => this%ch_ket
    J = ch_bra%j ! only scalar

    !$omp parallel
    !$omp do private(bra, qbra, cfp_bra, ket, qket, cfp_ket, &
    !$omp &  vin, bra_sub, a, b, c, Jab, ket_sub, d, e, f, Jde, &
    !$omp &  vout, bra_start, bra_end, ket_start, ket_end)
    do bra = 1, ch_bra%n_idx
      qbra => ch_bra%idx(bra)
      call cfp_bra%ini(qbra%north,qbra%nphys)
      cfp_bra%m = transpose(qbra%cfp)
      do ket = 1, ch_ket%n_idx
        qket => ch_ket%idx(ket)

        call cfp_ket%ini(qket%nphys,qket%north)
        cfp_ket%m = qket%cfp

        call vin%zeros(qbra%nphys,qket%nphys)

        do bra_sub = 1, qbra%nphys
          a = qbra%n2spi1(bra_sub)
          b = qbra%n2spi2(bra_sub)
          c = qbra%n2spi3(bra_sub)
          Jab= qbra%n2J12(bra_sub)
          do ket_sub = 1, qket%nphys
          d = qket%n2spi1(ket_sub)
          e = qket%n2spi2(ket_sub)
          f = qket%n2spi3(ket_sub)
          Jde= qket%n2J12(ket_sub)

          vin%m(bra_sub,ket_sub) = op3_in%GetThBME_pn(&
              & a,b,c,Jab,d,e,f,Jde,J)

          end do
        end do

        bra_start = qbra%idx2n(1)
        bra_end = qbra%idx2n(qbra%north)
        ket_start = qket%idx2n(1)
        ket_end = qket%idx2n(qket%north)
        vout = cfp_bra * vin * cfp_ket
        this%m(bra_start:bra_end, ket_start:ket_end) = vout%m

        call cfp_ket%fin()
        call vout%fin()
        call vin%fin()

      end do
      call cfp_bra%fin()
    end do
    !$omp end do
    !$omp end parallel
  end subroutine SetThreeBodyPartChannel

  subroutine PrintThreeBodyPart(this, wunit)
    use ClassSys, only: sys
    class(ThreeBodyPart), intent(in) :: this
    integer, intent(in), optional :: wunit
    integer :: chbra, chket
    type(sys) :: s
    character(:), allocatable :: msg

    do chbra = 1, this%thr%NChan
      do chket = 1, this%thr%NChan
        if(.not. this%MatCh(chbra,chket)%is) cycle
        msg = trim(this%oprtr) // " " // trim(s%str(chbra)) &
            &  // " " // trim(s%str(chket))
        call this%MatCh(chbra,chket)%prt(msg=msg,iunit=wunit)
      end do
    end do
  end subroutine PrintThreeBodyPart

  function NormalOrderingFrom3To2(this, two) result(r)
    use MyLibrary, only: triag, sjs
    class(ThreeBodyPart), intent(in) :: this
    type(TwoBodySpace), intent(in) :: two
    type(TwoBodyPart) :: r
    type(SingleParticleOrbit), pointer :: oh
    integer :: i1, i2, i3, i4, ih, J12, J34, jh, Jbra, Jket
    integer :: e1, e2, e3, e4, eh
    integer :: chbra, chket, bra, ket, bra_max, ket_max
    real(8) :: v, fact, vsum, tfact

    call r%init(two,this%Scalar,this%oprtr,this%jr,this%pr,this%zr)
    do chbra = 1, two%NChan
      do chket = 1, chbra
        if( .not. r%MatCh(chbra,chket)%is) cycle

        bra_max = two%jpz(chbra)%n_state
        ket_max = two%jpz(chket)%n_state

        J12 = two%jpz(chbra)%j
        J34 = two%jpz(chket)%j

        do bra = 1, bra_max
          i1 = two%jpz(chbra)%n2spi1(bra)
          i2 = two%jpz(chbra)%n2spi2(bra)
          e1 = two%sps%orb(i1)%e
          e2 = two%sps%orb(i2)%e
          if(chbra == chket) ket_max = bra
          do ket = 1, ket_max
            i3 = two%jpz(chket)%n2spi1(ket)
            i4 = two%jpz(chket)%n2spi2(ket)
            e3 = two%sps%orb(i3)%e
            e4 = two%sps%orb(i4)%e

            fact = 1.d0
            if(i1==i2) fact = fact / dsqrt(2.d0)
            if(i3==i4) fact = fact / dsqrt(2.d0)

            vsum = 0.d0
            do ih = 1, two%sps%norbs
              oh => two%sps%GetOrbit(ih)
              if(oh%occ < 1.d-6) cycle
              jh = oh%j
              eh = oh%e
              if(e1+e2+eh > this%thr%e3max) cycle
              if(e3+e4+eh > this%thr%e3max) cycle

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
                  v = v + tfact * oh%occ * &
                      & this%GetThBMEJ(i1,i2,ih,J12,i3,i4,ih,J34,Jket)
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

end module ThreeBodyOperator
