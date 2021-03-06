module OneBodyOperator
  use omp_lib
  use myfort
  use OneBodyModelSpace
  implicit none

  public :: OneBodyPartChannel
  public :: OneBodyPart

  private :: InitOneBodyPartChannel
  private :: FinOneBodyPartChannel
  private :: SetOneBodyPartChannel

  private :: InitOneBodyPart
  private :: FinOneBodyPart
  private :: SetOneBodyPart
  private :: GetOBME
  private :: SetOBME
  private :: CopyOneBodyPart
  private :: SumOneBodyPart
  private :: SubtractOneBodyPart
  private :: ScaleOneBodyPart
  private :: PrintOneBodyPart
  private :: NormalOrderingFrom1To0

  private :: GetDenominator1
  private :: GetDenominator2
  private :: GetDenominator3
  private :: GetDenominator4

  private :: ReadOneBodyFile
  private :: read_one_body_gennari

  type, extends(DMat) :: OneBodyPartChannel
    type(OneBodyChannel), pointer :: ch_bra, ch_ket
    logical :: is = .false.
  contains
    procedure :: InitOneBodyPartChannel
    procedure :: FinOneBodyPartChannel
    procedure :: SetOneBodyPartChannel

    generic :: init => InitOneBodyPartChannel
    generic :: release => FinOneBodyPartChannel
    generic :: set => SetOneBodyPartChannel
  end type OneBodyPartChannel

  type :: OneBodyPart
    type(OneBodyPartChannel), allocatable :: MatCh(:,:)
    type(OneBodySpace), pointer :: one
    character(:), allocatable :: oprtr
    logical :: Scalar
    integer :: jr, pr, tr, zr
  contains
    procedure :: InitOneBodyPart
    procedure :: FinOneBodyPart
    procedure :: SetOneBodyPart
    procedure :: GetOBME
    procedure :: SetOBME
    procedure :: CopyOneBodyPart
    procedure :: SumOneBodyPart
    procedure :: SubtractOneBodyPart
    procedure :: ScaleOneBodyPart
    procedure :: PrintOneBodyPart
    procedure :: NormalOrderingFrom1To0
    procedure :: GetDenominator1
    procedure :: GetDenominator2
    procedure :: GetDenominator3
    procedure :: GetDenominator4
    procedure :: ReadOneBodyFile

    generic :: init => InitOneBodyPart
    generic :: fin => FinOneBodyPart
    generic :: set => SetOneBodyPart
    generic :: prt => PrintOneBodyPart
    generic :: assignment(=) => CopyOneBodyPart
    generic :: operator(+) => SumOneBodyPart
    generic :: operator(-) => SubtractOneBodyPart
    generic :: operator(*) => ScaleOneBodyPart
  end type OneBodyPart
contains

  subroutine FinOneBodyPart(this)
    class(OneBodyPart), intent(inout) :: this
    integer :: chbra, chket

    if(.not. allocated(this%MatCh)) return
    do chbra = 1, this%one%NChan
      do chket = 1, this%one%NChan
        call this%MatCh(chbra,chket)%release()
      end do
    end do
    deallocate(this%MatCh)
    this%one => null()

  end subroutine FinOneBodyPart

  subroutine InitOneBodyPart(this, one, Scalar, oprtr, jr, pr, zr)
    class(OneBodyPart), intent(inout) :: this
    type(OneBodySpace), target, intent(in) :: one
    logical, intent(in) :: Scalar
    character(*), intent(in) :: oprtr
    integer, intent(in) :: jr, pr, zr
    integer :: chbra, chket
    integer :: jbra, pbra, zbra, nbra
    integer :: jket, pket, zket, nket

    if(allocated(this%MatCh)) call this%fin()

    this%one => one
    this%oprtr = oprtr
    this%Scalar = Scalar
    this%jr = jr
    this%pr = pr
    this%zr = zr
    allocate(this%MatCh(one%NChan, one%NChan))
    do chbra = 1, one%NChan
      jbra = one%jpz(chbra)%j
      pbra = one%jpz(chbra)%p
      zbra = one%jpz(chbra)%z
      nbra = one%jpz(chbra)%n_state
      do chket = 1, one%NChan
        jket = one%jpz(chket)%j
        pket = one%jpz(chket)%p
        zket = one%jpz(chket)%z
        nket = one%jpz(chket)%n_state

        if(triag(jbra, jket, 2*jr)) cycle
        if(pbra * pket * pr /= 1) cycle
        if(zbra - 2*zr - zket /= 0) cycle
        call this%MatCh(chbra,chket)%init(one%jpz(chbra),one%jpz(chket))
      end do
    end do

  end subroutine InitOneBodyPart

  subroutine SetOneBodyPart(this, hw, A, Z, N)
    class(OneBodyPart), intent(inout) :: this
    type(OneBodySpace), pointer :: obs
    type(Orbits), pointer :: sps
    real(8), intent(in) :: hw
    integer, intent(in) :: A, Z, N
    integer :: chbra
    integer :: chket
    obs => this%one
    sps => obs%sps
    do chbra = 1, obs%NChan
      do chket = 1, obs%NChan
        if(.not. this%MatCh(chbra,chket)%is) cycle
        call this%MatCh(chbra,chket)%set(this%oprtr, sps, hw, A, Z, N)
      end do
    end do
  end subroutine SetOneBodyPart

  subroutine CopyOneBodyPart(a, b)
    class(OneBodyPart), intent(out) :: a
    type(OneBodyPart), intent(in) :: b
    integer :: chbra, chket

    if(allocated(a%MatCh)) call a%fin()
    a%one => b%one
    a%oprtr = b%oprtr
    a%Scalar = b%Scalar
    a%jr = b%jr
    a%pr = b%pr
    a%zr = b%zr

    allocate(a%MatCh(b%one%NChan,b%one%NChan))
    do chbra = 1, b%one%NChan
      do chket = 1, b%one%NChan
        if(.not. b%MatCh(chbra,chket)%is) cycle
        a%MatCh(chbra,chket)%is = b%MatCh(chbra,chket)%is
        a%MatCh(chbra,chket)%DMat = b%MatCh(chbra,chket)%DMat
        a%MatCh(chbra,chket)%ch_bra => b%MatCh(chbra,chket)%ch_bra
        a%MatCh(chbra,chket)%ch_ket => b%MatCh(chbra,chket)%ch_ket
      end do
    end do
  end subroutine CopyOneBodyPart

  function SumOneBodyPart(a, b) result(c)
    class(OneBodyPart), intent(in) :: a, b
    type(OneBodyPart) :: c
    integer :: chbra, chket
    c = a
    do chbra = 1, b%one%NChan
      do chket = 1, b%one%NChan
        if(.not. b%MatCh(chbra,chket)%is) cycle
        c%MatCh(chbra,chket)%DMat = &
            & a%MatCh(chbra,chket)%DMat + b%MatCh(chbra,chket)%DMat
      end do
    end do
  end function SumOneBodyPart

  function SubtractOneBodyPart(a, b) result(c)
    class(OneBodyPart), intent(in) :: a, b
    type(OneBodyPart) :: c
    integer :: chbra, chket
    c = a
    do chbra = 1, b%one%NChan
      do chket = 1, b%one%NChan
        if(.not. b%MatCh(chbra,chket)%is) cycle
        c%MatCh(chbra,chket)%DMat = &
            & a%MatCh(chbra,chket)%DMat - b%MatCh(chbra,chket)%DMat
      end do
    end do
  end function SubtractOneBodyPart

  function ScaleOneBodyPart(a, b) result(c)
    class(OneBodyPart), intent(in) :: a
    real(8), intent(in) :: b
    type(OneBodyPart) :: c
    integer :: chbra, chket
    c = a
    do chbra = 1, a%one%NChan
      do chket = 1, a%one%NChan
        if(.not. a%MatCh(chbra,chket)%is) cycle
        c%MatCh(chbra,chket)%DMat = &
            & a%MatCh(chbra,chket)%DMat * b
      end do
    end do
  end function ScaleOneBodyPart

  subroutine FinOneBodyPartChannel(this)
    class(OneBodyPartChannel), intent(inout) :: this
    if(.not. this%is) return
    call this%fin()
    this%is = .false.
    this%ch_bra => null()
    this%ch_ket => null()
  end subroutine FinOneBodyPartChannel

  subroutine InitOneBodyPartChannel(this, ch_bra, ch_ket)
    class(OneBodyPartChannel), intent(inout) :: this
    type(OneBodyChannel), target, intent(in) :: ch_bra, ch_ket
    this%ch_bra => ch_bra
    this%ch_ket => ch_ket
    call this%zeros(ch_bra%n_state, ch_ket%n_state)
    this%is = .true.
  end subroutine InitOneBodyPartChannel

  subroutine SetOneBodyPartChannel(this, optr, sps, hw, A, Z, N)
    class(OneBodyPartChannel), intent(inout) :: this
    character(*), intent(in) :: optr
    type(Orbits), pointer, intent(in) :: sps
    type(OneBodyChannel), pointer :: chbra, chket
    real(8), intent(in) :: hw
    integer, intent(in) :: A, Z, N
    integer :: bra, ket, ia, ib
    chbra => this%ch_bra
    chket => this%ch_ket
    do bra = 1, chbra%n_state
      ia = chbra%n2spi(bra)
      do ket = 1, chket%n_state
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
  end subroutine SetOneBodyPartChannel

  function GetOBME(this,i1,i2) result(r)
    class(OneBodyPart), intent(in) :: this
    integer, intent(in) :: i1, i2
    type(OneBodySpace), pointer :: one
    type(SingleParticleOrbit), pointer :: o1, o2
    integer :: chbra, chket, bra, ket
    real(8) :: r

    r = 0.d0
    one => this%one
    o1 => one%sps%GetOrbit(i1)
    o2 => one%sps%GetOrbit(i2)
    chbra = one%jpz2ch(o1%j, (-1)**o1%l, o1%z)
    chket = one%jpz2ch(o2%j, (-1)**o2%l, o2%z)
    if(.not. this%MatCh(chbra,chket)%is) return
    bra = one%jpz(chbra)%spi2n(i1)
    ket = one%jpz(chket)%spi2n(i2)
    r = this%MatCh(chbra,chket)%m(bra,ket)
  end function GetOBME

  subroutine SetOBME(this,i1,i2,me,set_symmetric)
    class(OneBodyPart), intent(inout) :: this
    integer, intent(in) :: i1, i2
    real(8), intent(in) :: me
    logical, intent(in), optional :: set_symmetric
    type(OneBodySpace), pointer :: one
    type(SingleParticleOrbit), pointer :: o1, o2
    integer :: chbra, chket, bra, ket
    real(8) :: r

    r = 0.d0
    one => this%one
    o1 => one%sps%GetOrbit(i1)
    o2 => one%sps%GetOrbit(i2)
    chbra = one%jpz2ch(o1%j, (-1)**o1%l, o1%z)
    chket = one%jpz2ch(o2%j, (-1)**o2%l, o2%z)
    if(.not. this%MatCh(chbra,chket)%is) return
    bra = one%jpz(chbra)%spi2n(i1)
    ket = one%jpz(chket)%spi2n(i2)
    this%MatCh(chbra,chket)%m(bra,ket) = me
    if(present(set_symmetric)) then
      if(.not. set_symmetric) return
    end if
    this%MatCh(chket,chbra)%m(ket,bra) = me*(-1.d0)**( (o2%j-o1%j)/2 )
  end subroutine SetOBME

  subroutine PrintOneBodyPart(this, wunit)
    class(OneBodyPart), intent(in) :: this
    integer, intent(in), optional :: wunit
    integer :: chbra, chket
    type(sys) :: s
    character(:), allocatable :: msg
    integer :: jbra, pbra, zbra, jket, pket, zket

    msg = ""
    do chbra = 1, this%one%NChan
      jbra = this%one%jpz(chbra)%j
      pbra = this%one%jpz(chbra)%p
      zbra = this%one%jpz(chbra)%z
      do chket = 1, this%one%NChan
        jket = this%one%jpz(chket)%j
        pket = this%one%jpz(chket)%p
        zket = this%one%jpz(chket)%z
        if(.not. this%MatCh(chbra,chket)%is) cycle
        msg = trim(this%oprtr) // " (" // trim(s%str(jbra)) // &
            & "," // trim(s%str(pbra)) // "," // trim(s%str(zbra)) // &
            & ")  (" // trim(s%str(jket)) // "," // &
            & trim(s%str(pket)) // "," // trim(s%str(zket)) // ")"
        call this%MatCh(chbra,chket)%prt(msg=msg,iunit=wunit)
      end do
    end do
  end subroutine PrintOneBodyPart

  function NormalOrderingFrom1To0(this) result(zero)
    class(OneBodyPart), intent(in) :: this
    type(Orbits), pointer :: sps
    type(SingleParticleOrbit), pointer :: io
    real(8) :: zero
    integer :: ch, n, i, j
    zero = 0.d0
    if(.not. this%Scalar) return
    sps => this%one%sps
    do ch = 1, this%one%NChan
      j = this%one%jpz(ch)%j
      do n = 1, this%one%jpz(ch)%n_state
        i = this%one%jpz(ch)%n2spi(n)
        io => sps%GetOrbit(i)
        if(io%GetOccupation() < 1.d-8) cycle
        zero = zero + dble(j+1) * this%MatCh(ch,ch)%m(n,n) * io%GetOccupation()
      end do
    end do
  end function NormalOrderingFrom1To0

  function GetDenominator1(f,h1,p1) result(r)
    class(OneBodyPart), intent(in) :: f
    integer, intent(in) :: h1, p1
    real(8) :: r

    r = f%GetOBME(h1,h1) - f%GetOBME(p1,p1)
  end function GetDenominator1

  function GetDenominator2(f,h1,h2,p1,p2) result(r)
    class(OneBodyPart), intent(in) :: f
    integer, intent(in) :: h1, h2, p1, p2
    real(8) :: r

    r = f%GetOBME(h1,h1) + f%GetOBME(h2,h2) - &
        & f%GetOBME(p1,p1) - f%GetOBME(p2,p2)
  end function GetDenominator2

  function GetDenominator3(f,h1,h2,h3,p1,p2,p3) result(r)
    class(OneBodyPart), intent(in) :: f
    integer, intent(in) :: h1, h2, h3, p1, p2, p3
    real(8) :: r

    r = f%GetOBME(h1,h1) + f%GetOBME(h2,h2) + f%GetOBME(h3,h3) - &
        & f%GetOBME(p1,p1) - f%GetOBME(p2,p2) - f%GetOBME(p3,p3)
  end function GetDenominator3

  function GetDenominator4(f,h1,h2,h3,h4,p1,p2,p3,p4) result(r)
    class(OneBodyPart), intent(in) :: f
    integer, intent(in) :: h1, h2, h3, h4, p1, p2, p3, p4
    real(8) :: r

    r = f%GetOBME(h1,h1) + f%GetOBME(h2,h2) + &
        & f%GetOBME(h3,h3) + f%GetOBME(h4,h4) - &
        & f%GetOBME(p1,p1) - f%GetOBME(p2,p2) - &
        & f%GetOBME(p3,p3) - f%GetOBME(p4,p4)
  end function GetDenominator4

  subroutine ReadOneBodyFile(this, filename, emax, lmax)
    class(OneBodyPart), intent(inout) :: this
    character(*), intent(in) :: filename
    integer, intent(in) :: emax, lmax
    type(sys) :: s

    if(s%find(filename, ".gen")) then
      call read_one_body_gennari(this, filename, emax, lmax)
      return
    end if

  end subroutine ReadOneBodyFile

  subroutine read_one_body_gennari(this, filename, emax, lmax)
    class(OneBodyPart), intent(inout) :: this
    character(*), intent(in) :: filename
    integer, intent(in) :: emax, lmax
    type(OrbitsIsospin) :: isps
    type(Orbits), pointer :: sps
    integer :: iunit = 20, io
    integer :: ia, ib, ap, bp, an, bn
    real(8) :: me_pp, me_nn

    sps => this%one%sps
    call isps%init(emax, lmax)
    open(iunit, file=filename, action='read',iostat=io)
    read(iunit,*)
    do
      read(iunit, *, end=999) ia, ib, me_pp, me_nn
      me_pp = me_pp / dble(isps%orb(ia)%j+1)
      me_nn = me_nn / dble(isps%orb(ia)%j+1)
      ap = isps%iso2pn(sps, ia, -1)
      an = isps%iso2pn(sps, ia,  1)
      bp = isps%iso2pn(sps, ib, -1)
      bn = isps%iso2pn(sps, ib,  1)
      call this%SetOBME(ap,bp,me_pp)
      call this%SetOBME(an,bn,me_nn)
    end do
999 close(iunit)
    call isps%fin()
    return
  end subroutine read_one_body_gennari

end module OneBodyOperator
