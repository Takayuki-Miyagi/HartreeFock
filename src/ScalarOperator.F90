module ScalarOperator
  use InputParameters, only: parameters
  use ModelSpace, only: spo_pn, MSpace, OneBodySpace, TwoBodySpace, ThreeBodySpace, &
      & OneBodyChannel, TwoBodyChannel, ThreeBodyChannel, SpinParityTz
  use MatrixDouble, only: DMat
  use VectorDouble, only: DVec
  use LinAlgLib
  implicit none
  private :: CopyScalarOperators, CopyNBodyScalars, &
      & SumScalarOperators, SumNBodyScalars, &
      & SubtractScalarOperators, SubtractNBodyScalars, &
      & ScaleLScalarOperators, ScaleLNBodyScalars, &
      & ScaleRScalarOperators, ScaleRNBodyScalars, &
      & DivideScalarOperators, DivideNBodyScalars, &
      & one_body_element, read_2bme_pn_txt, read_2bme_pn_bin, &
      & read_2bme_snt_txt, read_2bme_snt_bin, calc_bare_2bme, &
      & two_body_element, r_dot_r, red_r_j, red_r_l, &
      & p_dot_p, red_nab_j, red_nab_l, skip_comment, Del, &
      & OneBodyScalarEmbedded3, OneBodyScalarEmbedded2, &
      & TwoBodyScalarEmbedded3, TransOrthoNormalThree, &
      & get_element_twobody, get_Uelement_twobody
  public :: assignment(=), operator(+), operator(-), &
      & operator(*), operator(/), NBodyScalars, ScalarOperators, &
      & ScalarEmbedded, UEmbedded

  type :: NBodyScalars
    type(DMat), allocatable :: jptz(:)
    real(8) :: usedmem = 0.d0
  contains
    procedure :: init => InitNBodyScalars
    procedure :: fin => FinNBodyScalars
    procedure :: SetOneBodyScalars
    procedure :: SetTwoBodyScalars
    procedure :: SetThreeBodyScalars
    generic :: set => SetOneBodyScalars, SetTwoBodyScalars, SetThreeBodyScalars
  end type NBodyScalars

  type :: ScalarOperators
    type(NBodyScalars) :: one, two, thr
    real(8) :: zero = 0.d0
    real(8) :: usedmem = 0.d0
  contains
    procedure :: init => InitScalarOperators
    procedure :: fin => FinScalarOperators
    procedure :: set => SetScalarOperators
  end type ScalarOperators

  interface assignment(=)
    module procedure :: CopyScalarOperators, CopyNBodyScalars
  end interface assignment(=)

  interface operator(+)
    module procedure :: SumScalarOperators, SumNBodyScalars
  end interface operator(+)

  interface operator(-)
    module procedure :: SubtractScalarOperators, SubtractNBodyScalars
  end interface operator(-)

  interface operator(*)
    module procedure :: ScaleLScalarOperators, ScaleLNBodyScalars, &
          & ScaleRScalarOperators, ScaleRNBodyScalars
  end interface operator(*)

  interface operator(/)
    module procedure :: DivideScalarOperators, DivideNBodyScalars
  end interface operator(/)

  interface ScalarEmbedded
    module procedure :: OneBodyScalarEmbedded3, &
          & OneBodyScalarEmbedded2, TwoBodyScalarEmbedded3
  end interface ScalarEmbedded
contains
  ! Constructor and Destructor
  subroutine InitScalarOperators(this, ms)
    class(ScalarOperators), intent(inout) :: this
    type(MSpace), intent(in) :: ms
    call this%one%init(ms%one%SpinParityTz)
    call this%two%init(ms%two%SpinParityTz)
    if(allocated(ms%thr%jptz)) then
      call this%thr%init(ms%thr%SpinParityTz)
    end if
    this%usedmem = this%one%usedmem + this%two%usedmem + this%thr%usedmem
  end subroutine InitScalarOperators

  subroutine FinScalarOperators(this)
    class(ScalarOperators), intent(inout) :: this
    call this%one%fin()
    call this%two%fin()
    call this%thr%fin()
  end subroutine FinScalarOperators

  subroutine InitNBodyScalars(this, nbody)
    class(NBodyScalars), intent(out) :: this
    type(SpinParityTz), intent(in) :: nbody
    integer :: ich, n
    real(8) :: cnt
    if(allocated(this%jptz)) call this%fin()
    n = nbody%n
    allocate(this%jptz(n))
    cnt = 0.d0
    do ich = 1, nbody%n
      n = nbody%ndim(ich)
      call this%jptz(ich)%ini(n,n)
      cnt = cnt + n ** 2
    end do
    this%usedmem = 8.d0 * cnt / (1024.d0) ** 3
  end subroutine InitNBodyScalars

  subroutine FinNBodyScalars(this)
    class(NBodyScalars), intent(inout) :: this
    integer :: ich, n
    if(.not. allocated(this%jptz)) return
    n = size(this%jptz)
    do ich = 1, n
      call this%jptz(ich)%fin()
    end do
    deallocate(this%jptz)
  end subroutine FinNBodyScalars

  subroutine CopyScalarOperators(a, b)
    type(ScalarOperators), intent(inout) :: a
    type(ScalarOperators), intent(in) :: b
    call CopyNBodyScalars(a%one, b%one)
    call CopyNBodyScalars(a%two, b%two)
    if(allocated(b%thr%jptz)) then
      call CopyNBodyScalars(a%thr, b%thr)
    end if
  end subroutine CopyScalarOperators

  subroutine CopyNBodyScalars(a, b)
    type(NBodyScalars), intent(inout) :: a
    type(NBodyScalars), intent(in) :: b
    integer :: ich, n
    if( allocated(a%jptz)) call a%fin()
    n = size(b%jptz)
    allocate(a%jptz(n))
    do ich = 1, n
      a%jptz(ich) = b%jptz(ich)
    end do
  end subroutine CopyNBodyScalars

  function SumScalarOperators(a, b) result(c)
    type(ScalarOperators) :: c
    type(ScalarOperators), intent(in) :: a, b
    call CopyNBodyScalars(c%one, SumNBodyScalars(a%one, b%one))
    call CopyNBodyScalars(c%two, SumNBodyScalars(a%two, b%two))
    if(allocated(a%thr%jptz) .and. allocated(b%thr%jptz)) then
      call CopyNBodyScalars(c%thr, SumNBodyScalars(a%thr, b%thr))
    end if
  end function SumScalarOperators

  function SumNBodyScalars(a, b) result(c)
    type(NBodyScalars) :: c
    type(NBodyScalars), intent(in) :: a, b
    integer :: ich, n
    call CopyNBodyScalars(c, a)
    n = size(a%jptz)
    do ich = 1, n
      c%jptz(ich) = a%jptz(ich) + b%jptz(ich)
    end do
  end function SumNBodyScalars

  function SubtractScalarOperators(a, b) result(c)
    type(ScalarOperators) :: c
    type(ScalarOperators), intent(in) :: a, b
    call CopyNBodyScalars(c%one, SubtractNBodyScalars(a%one, b%one))
    call CopyNBodyScalars(c%two, SubtractNBodyScalars(a%two, b%two))
    if(allocated(a%thr%jptz) .and. allocated(b%thr%jptz)) then
      call CopyNBodyScalars(c%thr, SubtractNBodyScalars(a%thr, b%thr))
    end if
  end function SubtractScalarOperators

  function SubtractNBodyScalars(a, b) result(c)
    type(NBodyScalars) :: c
    type(NBodyScalars), intent(in) :: a, b
    integer :: ich, n
    n = size(a%jptz)
    call CopyNBodyScalars(c, a)
    do ich = 1, n
      c%jptz(ich) = a%jptz(ich) - b%jptz(ich)
    end do
  end function SubtractNBodyScalars

  function ScaleLScalarOperators(a, b) result(c)
    type(ScalarOperators) :: c
    real(8), intent(in) :: a
    type(ScalarOperators), intent(in) :: b
    call CopyNBodyScalars(c%one, ScaleLNBodyScalars(a, b%one))
    call CopyNBodyScalars(c%two, ScaleLNBodyScalars(a, b%two))
    if(allocated(b%thr%jptz)) then
      call CopyNBodyScalars(c%thr, ScaleLNBodyScalars(a, b%thr))
    end if
  end function ScaleLScalarOperators

  function ScaleLNBodyScalars(a, b) result(c)
    type(NBodyScalars) :: c
    real(8), intent(in) :: a
    type(NBodyScalars), intent(in) :: b
    integer :: ich, n
    n = size(b%jptz)
    call CopyNBodyScalars(c, b)
    do ich = 1, n
      c%jptz(ich) = a * b%jptz(ich)
    end do
  end function ScaleLNBodyScalars

  function ScaleRScalarOperators(b, a) result(c)
    type(ScalarOperators) :: c
    real(8), intent(in) :: a
    type(ScalarOperators), intent(in) :: b
    call CopyNBodyScalars(c%one, ScaleRNBodyScalars(b%one, a))
    call CopyNBodyScalars(c%two, ScaleRNBodyScalars(b%two, a))
    if(allocated(b%thr%jptz)) then
      call CopyNBodyScalars(c%thr, ScaleRNBodyScalars(b%thr, a))
    end if
  end function ScaleRScalarOperators

  function ScaleRNBodyScalars(b, a) result(c)
    type(NBodyScalars) :: c
    real(8), intent(in) :: a
    type(NBodyScalars), intent(in) :: b
    integer :: ich, n
    n = size(b%jptz)
    call CopyNBodyScalars(c, b)
    do ich = 1, n
      c%jptz(ich) = a * b%jptz(ich)
    end do
  end function ScaleRNBodyScalars

  function DivideScalarOperators(b, a) result(c)
    type(ScalarOperators) :: c
    real(8), intent(in) :: a
    type(ScalarOperators), intent(in) :: b
    call CopyNBodyScalars(c%one, DivideNBodyScalars(b%one, a))
    call CopyNBodyScalars(c%two, DivideNBodyScalars(b%two, a))
    if(allocated(b%thr%jptz)) then
      call CopyNBodyScalars(c%thr, DivideNBodyScalars(b%thr, a))
    end if
  end function DivideScalarOperators

  function DivideNBodyScalars(b, a) result(c)
    type(NBodyScalars) :: c
    real(8), intent(in) :: a
    type(NBodyScalars), intent(in) :: b
    integer :: ich, n
    n = size(b%jptz)
    call CopyNBodyScalars(c, b)
    do ich = 1, n
      c%jptz(ich) = b%jptz(ich) / a
    end do
  end function DivideNBodyScalars

  subroutine SetScalarOperators(this, params, sps, ms, oprtr, f2, thbme)
    use read_3BME, only: iThreeBodyScalar
    class(ScalarOperators), intent(inout) :: this
    type(parameters), intent(in) :: params
    type(spo_pn), intent(in) :: sps
    type(MSpace), intent(in) :: ms
    character(*), intent(in) :: oprtr
    character(*), optional, intent(in) :: f2
    type(iThreeBodyScalar), optional, intent(in) :: thbme
    type(NBodyScalars) :: two

    call this%one%SetOneBodyScalars(params, sps, ms%one, oprtr)
    call this%two%SetTwoBodyScalars(params, sps, ms%two, oprtr, f2)
    if( allocated (ms%thr%jptz)) then
      call this%thr%SetThreeBodyScalars(params, sps, ms%thr, thbme)
    end if

    if(oprtr == 'hamil' .or. oprtr == 'Hamil') then
      call two%init(ms%two%SpinParityTz)
      call calc_bare_2bme('H_p_dot_p', params, sps, ms%two, two)
      call CopyNBodyScalars(this%two, SubtractNBodyScalars(this%two, two))
      call two%fin()
    end if
  end subroutine SetScalarOperators

  subroutine SetOneBodyScalars(this, params, sps, one, oprtr)
    class(NBodyScalars), intent(inout) :: this
    type(parameters), intent(in) :: params
    type(spo_pn), intent(in) :: sps
    type(OneBodySpace), intent(in) :: one
    character(*), intent(in) :: oprtr
    integer :: ich, j, p, itz, ichmax
    integer :: bra, ket, n
    integer :: i1, i2, l1, n1, n2
    real(8) :: e
    ichmax = one%n
    do ich = 1, ichmax
      j   = one%j(ich)
      p   = one%p(ich)
      itz = one%tz(ich)
      n   = one%jptz(ich)%n
      do bra = 1, n
        i1 = one%jptz(ich)%n2label(bra)
        n1 = sps%nn(i1)
        l1 = sps%ll(i1)
        do ket = 1, n
          i2 = one%jptz(ich)%n2label(ket)
          n2 = sps%nn(i2)
          e = one_body_element(n1, n2, l1, j, itz, params, oprtr)
          this%jptz(ich)%m(bra, ket) = e
        end do
      end do
    end do
  end subroutine SetOneBodyScalars

  subroutine SetTwoBodyScalars(this, params, sps, two, oprtr, f2)
    use class_sys, only: sy
    class(NBodyScalars), intent(inout) :: this
    type(parameters), intent(in) :: params
    type(spo_pn), intent(in) :: sps
    type(TwoBodySpace), intent(in) :: two
    character(*), intent(in) :: oprtr
    character(*), optional, intent(in) :: f2
    logical :: ex
    type(sy) :: sys

    if(oprtr == 'Hamil' .or. oprtr == 'hamil') then

      if(.not. present(f2)) then
        write(*,'(a)') 'Error NN int file is not detected!'
      end if

      inquire(file = f2, exist = ex)
      if(.not. ex) then
        write(*, '(a, " does not exist.")') trim(f2)
        stop
      end if

    else

      if(.not. present(f2)) then
        call calc_bare_2bme(oprtr, params, sps, two, this)
        return
      end if

    end if

    if(sys%find(f2, '.txt') .and. sys%find(f2, '.myg')) then
      call read_2bme_pn_txt(f2, params, sps, two, this)
    elseif(sys%find(f2, '.bin') .and. sys%find(f2, '.myg')) then
      call read_2bme_pn_bin(f2, params, sps, two, this)
    elseif(sys%find(f2, '.txt') .and. sys%find(f2, '.snt')) then
      call read_2bme_snt_txt(f2, sps, two, this)
    elseif(sys%find(f2, '.bin') .and. sys%find(f2, '.snt')) then
      call read_2bme_snt_bin(f2, sps, two, this)
    end if
  end subroutine SetTwoBodyScalars

  subroutine SetThreeBodyScalars(this, params, sps, thr, thbme)
    use read_3BME, only: iThreeBodyScalar
    class(NBodyScalars), intent(inout) :: this
    type(parameters), intent(in) :: params
    type(spo_pn), intent(in) :: sps
    type(ThreeBodySpace), intent(in) :: thr
    type(iThreeBodyScalar), optional, intent(in) :: thbme
    integer :: ich, j, p, tz
    if(.not. present(thbme)) return
    do ich = 1, thr%n
      j   = thr%j(ich)
      p   = thr%p(ich)
      tz = thr%tz(ich)
      call TransOrthoNormalThree(thr%jptz(ich), this%jptz(ich), thbme, sps, params, j, p, tz)
    end do

  end subroutine SetThreeBodyScalars

  real(8) function Get3BMEpn(sps, thbme, &
        & j, p, itz, a, b, c, jab, d, e, f, jde) result(v)
    use read_3BME, only: iThreeBodyScalar
    use RotationGroup, only: dcg
    use read_3BME, only: Get3BME
    type(iThreeBodyScalar), intent(in) :: thbme
    type(spo_pn), intent(in) :: sps
    integer, intent(in) :: j, p, itz
    integer, intent(in) :: a, b, c, d, e, f, jab, jde
    integer :: za, zb, zc, zd, ze, zf, tab, tde
    integer :: t
    integer :: ich
    real(8) :: viso
    v = 0.d0
    if(a == b .and. mod(jab, 2) == 1) return
    if(d == e .and. mod(jde, 2) == 1) return
    za = sps%itz(a); zb = sps%itz(b); zc = sps%itz(c)
    zd = sps%itz(d); ze = sps%itz(e); zf = sps%itz(f)
    do tab = 0, 1
      if(iabs(za + zb)/2 > tab) cycle
      do tde = 0, 1
        if(iabs(zd + ze)/2 > tde) cycle
        do t = max(iabs(2*tab - 1), iabs(2*tde - 1)), &
              & min(2*tab + 1, 2*tde + 1), 2
          if(iabs(itz) > t) cycle
          ich = thbme%jpt2n(j,p,t)
          if(ich == 0) cycle
          viso =  dble(Get3BME(thbme%jpt(ich), (a+1)/2, (b+1)/2, (c+1)/2, &
              &  jab, tab, (d+1)/2, (e+1)/2, (f+1)/2, jde, tde))
          v = v + viso * &
              &  dcg(1, za, 1, zb, 2*tab, za+zb) * dcg(2*tab, za+zb, 1, zc, t, itz) * &
              &  dcg(1, zd, 1, ze, 2*tde, zd+ze) * dcg(2*tde, zd+ze, 1, zf, t, itz)
        end do
      end do
    end do
  end function Get3BMEpn

  function one_body_element(n1, n2, l, j, tz, params, oprtr) result(e)
    real(8) :: e
    type(parameters), intent(in) :: params
    integer, intent(in) :: n1, n2, l, j, tz
    integer :: A, Z, N
    real(8) :: hw, hc, am
    character(*), intent(in) :: oprtr
    A = params%mass
    Z = params%pmass
    N = params%nmass
    hw = params%hw
    hc = params%hc
    am = params%amnucl

    e = 0.d0
    select case(oprtr)
    case('Hamil','hamil')
      if(n1 == n2) e = dble(2 * n1 + l) + 1.5d0
      if(n1 == n2 + 1) e = dsqrt(dble(n1) * (dble(n1 + l) + 0.5d0))
      if(n1 == n2 - 1) e = dsqrt(dble(n2) * (dble(n2 + l) + 0.5d0))
      e = e * (1.d0 - 1.d0 / dble(A)) * 0.5d0 * hw
    case('Rm','rm')
      if(n1 == n2) e = (dble(2 * n1 + l) + 1.5d0)
      if(n1 == n2 + 1) e = -dsqrt(dble(n1) * (dble(n1 + l) + 0.5d0))
      if(n1 == n2 - 1) e = -dsqrt(dble(n2) * (dble(n2 + l) + 0.5d0))
      e = e * (1.d0 - 1.d0 / dble(A)) * hc ** 2 / (am * hw * A)
    end select
  end function one_body_element

  subroutine calc_bare_2bme(oprtr, params, sps, two, tbme)
    type(NBodyScalars), intent(inout) :: tbme
    type(parameters), intent(in) :: params
    type(spo_pn), intent(in) :: sps
    type(TwoBodySpace), intent(in) :: two
    character(*), intent(in) :: oprtr
    integer :: j, ipar, itz
    integer :: a, b, c, d
    integer :: ich, i, k, n

    do ich = 1, two%n
      j    = two%j(ich)
      ipar = two%p(ich)
      itz  = two%tz(ich)
      n    = two%jptz(ich)%n
      !$omp parallel
      !$omp do private(i, a, b, k, c, d)
      do i = 1, n
        a = two%jptz(ich)%n2label1(i)
        b = two%jptz(ich)%n2label2(i)
        do k = 1, i
          c = two%jptz(ich)%n2label1(k)
          d = two%jptz(ich)%n2label2(k)
          tbme%jptz(ich)%m(i,k) = two_body_element(oprtr, params, sps, a, b, c, d, j)
          tbme%jptz(ich)%m(k,i) = tbme%jptz(ich)%m(i,k)
        end do
      end do
      !$omp end do
      !$omp end parallel
    end do
  end subroutine calc_bare_2bme

  function two_body_element(oprtr, params, sps, i1, i2, i3, i4, j) result(ele)
    real(8) :: ele
    type(parameters), intent(in) :: params
    type(spo_pn), intent(in) :: sps
    integer, intent(in) :: i1, i2, i3, i4, j
    character(*), intent(in) :: oprtr
    integer :: A, N, Z
    real(8) :: hc, hw, am
    A = params%mass
    Z = params%pmass
    N = params%nmass
    hc = params%hc
    hw = params%hw
    am = params%amnucl
    ele = 0.d0
    select case(oprtr)
    case('H_r_dot_r') ! rirj / A
      ele = r_dot_r(sps, i1, i2, i3, i4, j) * hw / dble(A)
    case('H_p_dot_p') ! pipj / A
      ele = p_dot_p(sps, i1, i2, i3, i4, j) * hw / dble(A)
    case('R_r_dot_r') ! - 4 * rirj / A**2 mc**2 hw
      ele = - 4.d0 * r_dot_r(sps, i1, i2, i3, i4, j) * hc ** 2 / (dble(A) * am * hw)
    end select
  end function two_body_element

  real(8) function r_dot_r(this, a, b, c, d, j) result(r)
    use RotationGroup, only:sjs
    type(spo_pn) :: this
    integer, intent(in) :: a, b, c, d, j
    integer :: ja, jb, jc, jd
    real(8) :: delab, delcd

    r = 0.d0
    ja = this%jj(a)
    jb = this%jj(b)
    jc = this%jj(c)
    jd = this%jj(d)

    delab = 1.d0
    delcd = 1.d0
    if(a == b) delab = 1.d0 / dsqrt(2.d0)
    if(c == d) delcd = 1.d0 / dsqrt(2.d0)

    r = (-1.d0) ** ((jb + jc) / 2 + j) * &
        & sjs(ja, jb, 2 * j, jd, jc, 2) * &
        & red_r_j(this, a, c) * red_r_j(this, b, d) + &
        & (-1.d0) ** ((jb + jc)/2) * &
        & sjs(ja, jb, 2 * j, jc, jd, 2) * &
        & red_r_j(this, a, d) * red_r_j(this, b, c)

    r = r * delab * delcd
  end function r_dot_r

  real(8) function red_r_j(this, a, b) result(rj)
    use RotationGroup, only: sjs
    type(spo_pn) :: this
    integer, intent(in) :: a, b
    integer :: na, la, ja, iza
    integer :: nb, lb, jb, izb

    rj = 0.d0

    na = this%nn(a)
    la = this%ll(a)
    ja = this%jj(a)
    iza = this%itz(a)

    nb = this%nn(b)
    lb = this%ll(b)
    jb = this%jj(b)
    izb = this%itz(b)

    if(iza /= izb) return

    rj = (-1.d0) ** ((3+2*la+jb)/2) * dsqrt(dble(ja + 1) * dble(jb + 1)) * &
        &  sjs(ja, 2, jb, 2 * lb, 1, 2 * la) * red_r_l(na, la, nb, lb)
  end function red_r_j

  real(8) function red_r_l(n1, l1, n2, l2) result(rl)
    integer, intent(in) :: n1, l1, n2, l2
    if (n1 == n2 .and. l1 == l2-1) then
      rl = -dsqrt(dble(l2)*(dble(n2 + l2) + 0.5d0))
    elseif (n1 == n2-1 .and. l1 == l2+1) then
      rl = -dsqrt(dble(l2 + 1)*dble(n2))
    elseif (n1 == n2+1 .and. l1 == l2-1) then
      rl = dsqrt(dble(l2)*(dble(n2  + 1)))
    elseif (n1 == n2 .and. l1==l2+1) then
      rl = dsqrt(dble(l2+1)*(dble(n2 +l2)+1.5d0))
    else
      rl = 0.d0
    end if
  end function red_r_l

  real(8) function p_dot_p(this, a, b, c, d, j) result(t)
    use RotationGroup, only:sjs
    type(spo_pn) :: this
    integer, intent(in) :: a, b, c, d, j
    integer :: ja, jb, jc, jd
    real(8) :: delab, delcd

    t = 0.d0
    ja = this%jj(a)
    jb = this%jj(b)
    jc = this%jj(c)
    jd = this%jj(d)

    delab = 1.d0
    delcd = 1.d0
    if(a == b) delab = 1.d0 / dsqrt(2.d0)
    if(c == d) delcd = 1.d0 / dsqrt(2.d0)

    t = - (-1.d0) ** ((jb + jc) / 2 + j) * &
        & sjs(ja, jb, 2 * j, jd, jc, 2) * &
        & red_nab_j(this, a, c) * red_nab_j(this, b, d) - &
        & (-1.d0) ** ((jb + jc) / 2) * &
        & sjs(ja, jb, 2 * j, jc, jd, 2) * &
        & red_nab_j(this, a, d) * red_nab_j(this, b, c)
    t = t * delab * delcd
  end function p_dot_p

  real(8) function red_nab_j(this, a, b) result(nj)
    use RotationGroup, only: sjs
    type(spo_pn) :: this
    integer, intent(in) :: a, b
    integer :: na, la, ja, iza
    integer :: nb, lb, jb, izb

    nj = 0.d0

    na = this%nn(a)
    la = this%ll(a)
    ja = this%jj(a)
    iza = this%itz(a)

    nb = this%nn(b)
    lb = this%ll(b)
    jb = this%jj(b)
    izb = this%itz(b)

    if(iza /= izb) return

    nj = (-1.d0) ** ((3+2*la+jb)/2) * dsqrt(dble(ja + 1) * dble(jb + 1)) * &
        &  sjs(ja, 2, jb, 2 * lb, 1, 2 * la) * red_nab_l(na, la, nb, lb)

  end function red_nab_j

  real(8) function red_nab_l(n1, l1, n2, l2) result(nl)
    integer, intent(in) :: n1, l1, n2, l2
    if(n1 == n2 .and. l1 == l2+1) then
      nl = -dsqrt(dble(l2 + 1)*(dble(n2 + l2) + 1.5d0))
    elseif(n1 == n2-1 .and. l1 == l2+1) then
      nl = -dsqrt(dble(l2 + 1)*dble(n2))
    elseif(n1 == n2 .and. l1 == l2-1) then
      nl = -dsqrt(dble(l2)*(dble(n2 + l2) + 0.5d0))
    elseif(n1 == n2+1 .and. l1==l2-1) then
      nl = -dsqrt(dble(l2)*dble(n2 + 1))
    else
      nl = 0.d0
    end if
  end function red_nab_l

  subroutine UEmbedded(two, this, one, op1)
    type(DMat), intent(inout) :: this
    type(OneBodySpace), intent(in) :: one
    type(TwoBodyChannel), intent(in) :: two
    type(NBodyScalars), intent(in) :: op1
    integer :: bra, ket, phb, a, b, c, d, n
    real(8) :: o12
    n = two%n
    !$omp parallel
    !$omp do private(bra, a, b, phb, ket, c, d, o12)
    do bra = 1, n
      a = two%n2label1(bra)
      b = two%n2label2(bra)
      phb = two%iphase(b,a)
      do ket = 1, n
        c = two%n2label1(ket)
        d = two%n2label2(ket)
        o12 = get_Uelement_twobody(a, b, c, d, phb, one, op1)
        this%m(bra,ket) = o12
      end do
    end do
    !$omp end do
    !$omp end parallel
  end subroutine UEmbedded

  function get_Uelement_twobody(a, b, c, d, phb, one, op1) result(o12)
    real(8) :: o12
    integer, intent(in) :: a, b, c, d, phb
    type(OneBodySpace), intent(in) :: one
    type(NBodyScalars), intent(in) :: op1
    integer :: icha, ichb, ichc, ichd
    real(8) :: u1, u2
    icha = one%label2jptz(a)
    ichb = one%label2jptz(b)
    ichc = one%label2jptz(c)
    ichd = one%label2jptz(d)
    o12 = 0.d0
    if(icha == ichc .and. ichb == ichd) then
      u1 = op1%jptz(icha)%m(one%jptz(icha)%label2n(a), one%jptz(icha)%label2n(c))
      u2 = op1%jptz(ichb)%m(one%jptz(ichb)%label2n(b), one%jptz(ichb)%label2n(d))
      o12 = u1 * u2
    end if

    if(icha == ichd .and. ichb == ichc) then
      u1 = op1%jptz(icha)%m(one%jptz(icha)%label2n(a), one%jptz(icha)%label2n(d))
      u2 = op1%jptz(ichb)%m(one%jptz(ichb)%label2n(b), one%jptz(ichb)%label2n(c))
      o12 = o12 + u1 * u2 * dble(phb)
    end if
    o12 = o12 / (Del(a,b) * Del(c, d))
  end function get_Uelement_twobody

  subroutine OneBodyScalarEmbedded2(two, this, one, op1, isgn)
    type(DMat), intent(inout) :: this
    type(OneBodySpace), intent(in) :: one
    type(TwoBodyChannel), intent(in) :: two
    type(NBodyScalars), intent(in) :: op1
    integer, intent(in) :: isgn
    integer :: bra, ket, phb, phk, a, b, c, d, n
    real(8) :: o12
    n = two%n
    !$omp parallel
    !$omp do private(bra, a, b, phb, ket, c, d, phk, o12)
    do bra = 1, n
      a = two%n2label1(bra)
      b = two%n2label2(bra)
      phb = two%iphase(b,a)
      do ket = 1, bra
        c = two%n2label1(ket)
        d = two%n2label2(ket)
        phk = two%iphase(d,c)
        o12 = get_element_twobody(a, b, c, d, phb, phk, one, op1)
        this%m(bra,ket) = o12
        this%m(ket,bra) = o12 * dble(isgn)
      end do
    end do
    !$omp end do
    !$omp end parallel
  end subroutine OneBodyScalarEmbedded2

  function get_element_twobody(a, b, c, d, phb, phk, one, op1) result(o12)
    real(8) :: o12
    integer, intent(in) :: a, b, c, d, phb, phk
    type(OneBodySpace), intent(in) :: one
    type(NBodyScalars), intent(in) :: op1
    integer :: icha, ichb, ichc, ichd
    icha = one%label2jptz(a)
    ichb = one%label2jptz(b)
    ichc = one%label2jptz(c)
    ichd = one%label2jptz(d)
    o12 = 0.d0
    if(b == d .and. icha == ichc) then
      o12 = o12 + op1%jptz(icha)%m(one%jptz(icha)%label2n(a), one%jptz(icha)%label2n(c))
    end if

    if(b == c .and. icha == ichd) then
      o12 = o12 + op1%jptz(icha)%m(one%jptz(icha)%label2n(a), one%jptz(icha)%label2n(d)) * dble(phk)
    end if

    if(a == d .and. ichb == ichc) then
      o12 = o12 + op1%jptz(ichb)%m(one%jptz(ichb)%label2n(b), one%jptz(ichb)%label2n(c)) * dble(phb)
    end if

    if(a == c .and. ichb == ichd) then
      o12 = o12 + op1%jptz(ichb)%m(one%jptz(ichb)%label2n(b), one%jptz(ichb)%label2n(d)) * dble(phb) * dble(phk)
    end if
    o12 = o12 / (Del(a,b) * Del(c, d))
  end function get_element_twobody

  subroutine OneBodyScalarEmbedded3(thr, this, one, op1)
    type(DMat), intent(inout) :: this
    type(OneBodySpace), intent(in) :: one
    type(ThreeBodyChannel), intent(in) :: thr
    type(NBodyScalars), intent(in) :: op1
    integer :: idxb
    integer :: njb
    type(DMat) :: cfpb, cfpk, vint, vfin
    integer :: iblck1, iblck2, iblck3, iblck4, blcksize = 5000
    integer :: n1, n2, m1, m2
    integer :: n1s, n1e, n2s, n2e, m1s, m1e, m2s, m2e
    integer :: ntot, mtot
    ntot = thr%n
    mtot = 0
    do idxb = 1, thr%nsub
      njb = thr%idx(idxb)%n
      mtot = mtot + njb
    end do

    do iblck1 = 1, mtot / blcksize + 1
      call block_dim(iblck1, blcksize, mtot, m1, m1s, m1e)
      if(m1 == 0) cycle
      do iblck2 = 1, mtot / blcksize + 1
        call block_dim(iblck2, blcksize, mtot, m2, m2s, m2e)
        if(m2 == 0) cycle
        call vint%ini(m1,m2)
        call MatOneBodyScalar(thr, vint, one, op1, m1s, m1e, m2s, m2e)

        do iblck3 = 1, ntot / blcksize + 1
          call block_dim(iblck3, blcksize, ntot, n1, n1s, n1e)
          if(n1 == 0) cycle
          call cfpb%ini(n1,m1)
          call SetCFPmat(cfpb, thr, n1s, n1e, m1s, m1e, 1)
          do iblck4 = 1, ntot / blcksize + 1
            call block_dim(iblck4, blcksize, ntot, n2, n2s, n2e)
            if(n2 == 0) cycle
            call cfpk%ini(m2,n2)
            call SetCFPmat(cfpk, thr, m2s, m2e, n2s, n2e, 0)
            vfin = (cfpb * (vint * cfpk))

            this%m(n1s:n1e, n2s:n2e) = this%m(n1s:n1e, n2s:n2e) + vfin%m(:,:)
          end do
        end do
      end do
    end do
  contains
    subroutine MatOneBodyScalar(thr, mat, one, op1, m1s, m1e, m2s, m2e)
      type(DMat), intent(inout) :: mat
      type(ThreeBodyChannel), intent(in) :: thr
      type(OneBodySpace), intent(in) :: one
      type(NBodyScalars), intent(in) :: op1
      integer, intent(in) :: m1s, m1e, m2s, m2e
      integer :: m1, m2, n, cnt1, cnt2, ich1, num1, num2
      integer :: idxb, njb, jb
      integer :: idxk, njk, jk
      integer :: a, b, c, jab, d, e, f, jde
      real(8) :: v

      n = thr%nsub
      cnt1 = 0
      m1 = 0
      !$omp parallel
      !$omp do private(idxb, njb, jb, a, b, c, jab, &
      !$omp &  cnt2, m2, idxk, njk, jk, d, e, f, jde, v, ich1, num1, num2) reduction(+: m1, cnt1)
      do idxb = 1, n
        njb = thr%idx(idxb)%n
        do jb = 1, njb
          m1 = m1 + 1
          if(m1 < m1s) cycle
          if(m1 > m1e) cycle
          cnt1 = cnt1 + 1
          a   = thr%idx(idxb)%n2label1(jb)
          b   = thr%idx(idxb)%n2label2(jb)
          c   = thr%idx(idxb)%n2label3(jb)
          jab = thr%idx(idxb)%n2label4(jb)

          cnt2 = 0
          m2 = 0
          do idxk = 1, n
            njk = thr%idx(idxk)%n
            do jk = 1, njk
              m2 = m2 + 1
              if(m2 < m2s) cycle
              if(m2 > m2e) cycle
              cnt2 = cnt2 + 1
              d   = thr%idx(idxk)%n2label1(jk)
              e   = thr%idx(idxk)%n2label2(jk)
              f   = thr%idx(idxk)%n2label3(jk)
              jde = thr%idx(idxk)%n2label4(jk)
              if(b /= e) cycle
              if(c /= f) cycle
              if(jab /= jde) cycle
              ich1 = one%label2jptz(a)
              if(ich1 /= one%label2jptz(d)) cycle
              num1 = one%jptz(ich1)%label2n(a)
              num2 = one%jptz(ich1)%label2n(d)
              v = op1%jptz(ich1)%m(num1, num2) * 3.d0
              mat%m(cnt1, cnt2) = v
            end do
          end do
        end do
      end do
      !$omp end do
      !$omp end parallel
    end subroutine MatOneBodyScalar
  end subroutine OneBodyScalarEmbedded3

  subroutine TwoBodyScalarEmbedded3(thr, this, two, op2, sps)
    type(DMat), intent(inout) :: this
    type(TwoBodySpace), intent(in) :: two
    type(ThreeBodyChannel), intent(in) :: thr
    type(NBodyScalars), intent(in) :: op2
    type(spo_pn), intent(in) :: sps
    integer :: idxb
    integer :: njb
    type(DMat) :: cfpb, cfpk, vint, vfin
    integer :: iblck1, iblck2, iblck3, iblck4, blcksize = 5000
    integer :: n1, n2, m1, m2
    integer :: n1s, n1e, n2s, n2e, m1s, m1e, m2s, m2e
    integer :: ntot, mtot
    ntot = thr%n
    mtot = 0
    do idxb = 1, thr%nsub
      njb = thr%idx(idxb)%n
      mtot = mtot + njb
    end do

    do iblck1 = 1, mtot / blcksize + 1
      call block_dim(iblck1, blcksize, mtot, m1, m1s, m1e)
      if(m1 == 0) cycle
      do iblck2 = 1, mtot / blcksize + 1
        call block_dim(iblck2, blcksize, mtot, m2, m2s, m2e)
        if(m2 == 0) cycle
        call vint%ini(m1,m2)
        call MatTwoBodyScalar(thr, vint, two, op2, sps, m1s, m1e, m2s, m2e)

        do iblck3 = 1, ntot / blcksize + 1
          call block_dim(iblck3, blcksize, ntot, n1, n1s, n1e)
          if(n1 == 0) cycle
          call cfpb%ini(n1,m1)
          call SetCFPmat(cfpb, thr, n1s, n1e, m1s, m1e, 1)
          do iblck4 = 1, ntot / blcksize + 1
            call block_dim(iblck4, blcksize, ntot, n2, n2s, n2e)
            if(n2 == 0) cycle
            call cfpk%ini(m2,n2)
            call SetCFPmat(cfpk, thr, m2s, m2e, n2s, n2e, 0)
            vfin = (cfpb * (vint * cfpk))

            this%m(n1s:n1e, n2s:n2e) = this%m(n1s:n1e, n2s:n2e) + vfin%m(:,:)
          end do
        end do
      end do
    end do
  contains
    subroutine MatTwoBodyScalar(thr, mat, two, op2, sps, m1s, m1e, m2s, m2e)
      type(DMat), intent(inout) :: mat
      type(spo_pn), intent(in) :: sps
      type(ThreeBodyChannel), intent(in) :: thr
      type(TwoBodySpace), intent(in) :: two
      type(NBodyScalars), intent(in) :: op2
      integer, intent(in) :: m1s, m1e, m2s, m2e
      integer :: m1, m2, n, cnt1, cnt2, ich2
      integer :: num1, num2, phase
      integer :: idxb, njb, jb
      integer :: idxk, njk, jk
      integer :: a, b, c, jab, d, e, f, jde
      real(8) :: v

      n = thr%nsub
      cnt1 = 0
      m1 = 0
      !$omp parallel
      !$omp do private(idxb, njb, jb, a, b, c, jab, &
      !$omp &  cnt2, m2, idxk, njk, jk, d, e, f, jde, v, ich2, num1, num2) reduction(+: m1, cnt1)
      do idxb = 1, n
        njb = thr%idx(idxb)%n
        do jb = 1, njb
          m1 = m1 + 1
          if(m1 < m1s) cycle
          if(m1 > m1e) cycle
          cnt1 = cnt1 + 1
          a   = thr%idx(idxb)%n2label1(jb)
          b   = thr%idx(idxb)%n2label2(jb)
          c   = thr%idx(idxb)%n2label3(jb)
          jab = thr%idx(idxb)%n2label4(jb)

          cnt2 = 0
          m2 = 0
          do idxk = 1, n
            njk = thr%idx(idxk)%n
            do jk = 1, njk
              m2 = m2 + 1
              if(m2 < m2s) cycle
              if(m2 > m2e) cycle
              cnt2 = cnt2 + 1
              d   = thr%idx(idxk)%n2label1(jk)
              e   = thr%idx(idxk)%n2label2(jk)
              f   = thr%idx(idxk)%n2label3(jk)
              jde = thr%idx(idxk)%n2label4(jk)
              if(c /= f) cycle
              if(jab /= jde) cycle
              ich2 = two%jptz2n(jab, (-1)**(sps%ll(a) + sps%ll(b)), (sps%itz(a) + sps%itz(b))/2)
              num1 = two%jptz(ich2)%labels2n(a,b)
              num2 = two%jptz(ich2)%labels2n(d,e)
              phase = two%jptz(ich2)%iphase(a,b) * two%jptz(ich2)%iphase(d,e)
              v = op2%jptz(ich2)%m(num1, num2) * dble(phase) * Del(a,b) * Del(d,e) * 1.5d0
              mat%m(cnt1, cnt2) = v
            end do
          end do
        end do
      end do
      !$omp end do
      !$omp end parallel
    end subroutine MatTwoBodyScalar
  end subroutine TwoBodyScalarEmbedded3

  subroutine TransOrthoNormalThree(thr, sc3, thbme, sps, params, j, p, itz)
    use read_3BME, only: iThreeBodyScalar
    type(DMat), intent(inout) :: sc3
    type(ThreeBodyChannel), intent(in) :: thr
    type(parameters), intent(in) :: params
    type(iThreeBodyScalar), intent(in) :: thbme
    type(spo_pn), intent(in) :: sps
    integer, intent(in) :: j, p, itz
    integer :: idxb
    integer :: njb
    type(DMat) :: cfpb, cfpk, vint, vfin
    integer :: iblck1, iblck2, iblck3, iblck4, blcksize = 5000
    integer :: n1, n2, m1, m2
    integer :: n1s, n1e, n2s, n2e, m1s, m1e, m2s, m2e
    integer :: ntot, mtot
    ntot = thr%n
    mtot = 0
    do idxb = 1, thr%nsub
      njb = thr%idx(idxb)%n
      mtot = mtot + njb
    end do

    do iblck1 = 1, mtot / blcksize + 1
      call block_dim(iblck1, blcksize, mtot, m1, m1s, m1e)
      if(m1 == 0) cycle
      do iblck2 = 1, mtot / blcksize + 1
        call block_dim(iblck2, blcksize, mtot, m2, m2s, m2e)
        if(m2 == 0) cycle
        call vint%ini(m1,m2)
        call MatThreeBodyScalar(vint, params, sps, thr, thbme, j, p, itz, m1s, m1e, m2s, m2e)

        do iblck3 = 1, ntot / blcksize + 1
          call block_dim(iblck3, blcksize, ntot, n1, n1s, n1e)
          if(n1 == 0) cycle
          call cfpb%ini(n1,m1)
          call SetCFPmat(cfpb, thr, n1s, n1e, m1s, m1e, 1)
          do iblck4 = 1, ntot / blcksize + 1
            call block_dim(iblck4, blcksize, ntot, n2, n2s, n2e)
            if(n2 == 0) cycle
            call cfpk%ini(m2,n2)
            call SetCFPmat(cfpk, thr, m2s, m2e, n2s, n2e, 0)
            vfin = (cfpb * (vint * cfpk))

            sc3%m(n1s:n1e, n2s:n2e) = sc3%m(n1s:n1e, n2s:n2e) + vfin%m(:,:)
          end do
        end do
      end do
    end do
  contains

    subroutine MatThreeBodyScalar(mat, params, sps, thr, thbme, j, p, itz, m1s, m1e, m2s, m2e)
    use read_3BME, only: iThreeBodyScalar
      type(DMat), intent(inout) :: mat
      type(parameters), intent(in) :: params
      type(spo_pn), intent(in) :: sps
      type(ThreeBodyChannel), intent(in) :: thr
      type(iThreeBodyScalar), intent(in) :: thbme
      integer, intent(in) :: j, p, itz, m1s, m1e, m2s, m2e
      integer :: m1, m2, n, cnt1, cnt2
      integer :: idxb, njb, jb
      integer :: idxk, njk, jk
      integer :: a, b, c, jab, d, e, f, jde
      real(8) :: v

      n = thr%nsub
      cnt1 = 0
      m1 = 0
      !$omp parallel
      !$omp do private(idxb, njb, jb, a, b, c, jab, &
      !$omp &  cnt2, m2, idxk, njk, jk, d, e, f, jde, v) reduction(+: m1, cnt1)
      do idxb = 1, n
        njb = thr%idx(idxb)%n
        do jb = 1, njb
          m1 = m1 + 1
          if(m1 < m1s) cycle
          if(m1 > m1e) cycle
          cnt1 = cnt1 + 1
          a   = thr%idx(idxb)%n2label1(jb)
          b   = thr%idx(idxb)%n2label2(jb)
          c   = thr%idx(idxb)%n2label3(jb)
          jab = thr%idx(idxb)%n2label4(jb)
          if(sps%nshell(a) > params%emax_3nf) cycle
          if(sps%nshell(b) > params%emax_3nf) cycle
          if(sps%nshell(c) > params%emax_3nf) cycle
          if(sps%nshell(a) + sps%nshell(b) > params%e2max_3nf) cycle
          if(sps%nshell(b) + sps%nshell(c) > params%e2max_3nf) cycle
          if(sps%nshell(c) + sps%nshell(a) > params%e2max_3nf) cycle
          if(sps%nshell(a) + sps%nshell(b) + sps%nshell(c) > params%e3max_3nf) cycle

          cnt2 = 0
          m2 = 0
          do idxk = 1, n
            njk = thr%idx(idxk)%n
            do jk = 1, njk
              m2 = m2 + 1
              if(m2 < m2s) cycle
              if(m2 > m2e) cycle
              cnt2 = cnt2 + 1
              d   = thr%idx(idxk)%n2label1(jk)
              e   = thr%idx(idxk)%n2label2(jk)
              f   = thr%idx(idxk)%n2label3(jk)
              jde = thr%idx(idxk)%n2label4(jk)
              if(sps%nshell(d) > params%emax_3nf) cycle
              if(sps%nshell(e) > params%emax_3nf) cycle
              if(sps%nshell(f) > params%emax_3nf) cycle
              if(sps%nshell(d) + sps%nshell(e) > params%e2max_3nf) cycle
              if(sps%nshell(e) + sps%nshell(f) > params%e2max_3nf) cycle
              if(sps%nshell(f) + sps%nshell(d) > params%e2max_3nf) cycle
              if(sps%nshell(d) + sps%nshell(e) + sps%nshell(f) > params%e3max_3nf) cycle
              v = Get3BMEpn(sps, thbme, &
                & j, p, itz, a, b, c, jab, d, e, f, jde) / 6.d0
              mat%m(cnt1, cnt2) = v
            end do
          end do
        end do
      end do
      !$omp end do
      !$omp end parallel
    end subroutine MatThreeBodyScalar
  end subroutine TransOrthoNormalThree

  subroutine block_dim(i, sizeb, n, nn, nstart, nend)
    integer, intent(in) :: i, sizeb, n
    integer, intent(out) :: nn, nstart, nend
    nn = min(sizeb, n - sizeb * (i-1))
    nstart = (i - 1) * sizeb + 1
    nend   = min(i * sizeb, n)
  end subroutine block_dim

  subroutine SetCFPmat(cfp, thr, n1s, n1e, n2s, n2e, nt)
    type(DMat), intent(inout) :: cfp
    type(ThreeBodyChannel), intent(in) :: thr
    integer, intent(in) :: n1s, n1e, n2s, n2e, nt
    integer :: m1s, m1e, m2s, m2e
    integer :: ket, m1, n, cnt1, cnt2
    integer :: kt, idxb, njb, jb
    integer :: i1, i2, i3

    m1s = n1s; m1e = n1e
    m2s = n2s; m2e = n2e
    if(nt == 1) then
      m1s = n2s; m1e = n2e
      m2s = n1s; m2e = n1e
    end if
    n = thr%nsub

    cnt1 = 0
    m1 = 0
    !$omp parallel
    !$omp do private(idxb, njb, jb, cnt2, ket, i1, i2, i3, kt) reduction(+: cnt1, m1)
    do idxb = 1, n
      njb = thr%idx(idxb)%n
      do jb = 1, njb
        m1 = m1 + 1
        if(m1 < m1s) cycle
        if(m1 > m1e) cycle
        cnt1 = cnt1 + 1

        cnt2 = 0
        do ket = 1, thr%n
          if(ket < m2s) cycle
          if(ket > m2e) cycle
          i1 = thr%n2label1(ket)
          i2 = thr%n2label2(ket)
          i3 = thr%n2label3(ket)
          kt = thr%n2label4(ket)
          cnt2 = cnt2 + 1
          if(thr%labels2nsub(i1,i2,i3) /= idxb) cycle
          if(nt == 0) cfp%m(cnt1, cnt2) = thr%idx(idxb)%cfp(jb, kt)
          if(nt == 1) cfp%m(cnt2, cnt1) = thr%idx(idxb)%cfp(jb, kt)
        end do
      end do
    end do
    !$omp end do
    !$omp end parallel
  end subroutine SetCFPmat

  !--------------------------------------------------
  !
  ! Read Two-Body Matrix element
  !
  !--------------------------------------------------

  subroutine read_2bme_pn_txt(f2, params, sps, two, tbme)
    type(NBodyScalars), intent(inout) :: tbme
    type(parameters), intent(in) :: params
    type(spo_pn), intent(in) :: sps
    type(TwoBodySpace), intent(in) :: two
    character(*), optional, intent(in) :: f2

    integer :: iunit = 19
    integer :: a, b, c, d
    integer :: ia, ib, ic, id
    integer :: iza, izb, izc, izd
    integer :: ij, numpot, num, ich
    integer :: itz, ip, bra, ket, phaseab, phasecd
    real(8) :: v12
    open(iunit, file = f2, status = "old")
    call skip_comment(iunit)
    read(iunit, *) numpot
    call skip_comment(iunit)
    do num = 1, numpot
      read(iunit, *) iza, ia, izb, ib, izc, ic, izd, id,  &
          & ij, v12
      a = 2 * ia + (iza - 1) / 2
      b = 2 * ib + (izb - 1) / 2
      c = 2 * ic + (izc - 1) / 2
      d = 2 * id + (izd - 1) / 2
      if(a > sps%n) cycle
      if(b > sps%n) cycle
      if(c > sps%n) cycle
      if(d > sps%n) cycle
      if(sps%nshell(a) + sps%nshell(b) > params%e2max) cycle
      if(sps%nshell(c) + sps%nshell(d) > params%e2max) cycle
      ip = (-1) ** (sps%ll(a) + sps%ll(b))
      itz = (iza + izb) / 2
      ich = two%jptz2n(ij, ip, itz)
      bra = two%jptz(ich)%labels2n(a, b)
      phaseab = two%jptz(ich)%iphase(a, b)
      ket = two%jptz(ich)%labels2n(c, d)
      phasecd = two%jptz(ich)%iphase(c, d)
      if(bra * ket == 0) cycle
      tbme%jptz(ich)%m(bra, ket) = v12 * dble(phaseab * phasecd)
      tbme%jptz(ich)%m(ket, bra) = v12 * dble(phaseab * phasecd)
    end do
    close(iunit)
  end subroutine read_2bme_pn_txt

  subroutine read_2bme_pn_bin(f2, params, sps, two, tbme)
    use class_stopwatch, only: start_stopwatch, stop_stopwatch, &
        & time_io_read
    type(NBodyScalars), intent(inout) :: tbme
    type(parameters), intent(in) :: params
    type(spo_pn), intent(in) :: sps
    type(TwoBodySpace), intent(in) :: two
    character(*), optional, intent(in) :: f2
    integer :: iunit = 19
    integer :: a, b, c, d
    integer :: ia, ib, ic, id
    integer :: iza, izb, izc, izd
    integer :: ij, numpot, num, ich
    integer :: itz, ip, bra, ket, phaseab, phasecd
    real(8) :: v12
    integer, allocatable :: aa(:), bb(:), cc(:), dd(:), jj(:)
    integer, allocatable :: az(:), bz(:), cz(:), dz(:)
    real(8), allocatable :: v1save(:)
    call start_stopwatch(time_io_read)
    open(iunit, form = "unformatted", file = f2, status = "old", access='stream')
    read(iunit) numpot
    allocate(aa(numpot), bb(numpot), cc(numpot), dd(numpot), jj(numpot))
    allocate(az(numpot), bz(numpot), cz(numpot), dz(numpot))
    allocate(v1save(numpot))
    read(iunit) az
    read(iunit) aa
    read(iunit) bz
    read(iunit) bb
    read(iunit) cz
    read(iunit) cc
    read(iunit) dz
    read(iunit) dd
    read(iunit) jj
    read(iunit) v1save
    close(iunit)
    call stop_stopwatch(time_io_read)
    !$omp parallel
    !$omp do private(iza, ia, izb, ib, izc, ic, izd, id, ij, v12, &
    !$omp &  a, b, c, d, ip, itz, ich, bra, phaseab, ket, phasecd)
    do num = 1, numpot
      iza = az(num)
      ia = aa(num)
      izb = bz(num)
      ib = bb(num)
      izc = cz(num)
      ic = cc(num)
      izd = dz(num)
      id = dd(num)
      ij = jj(num)
      v12 = v1save(num)
      a = 2 * ia + (iza - 1) / 2
      b = 2 * ib + (izb - 1) / 2
      c = 2 * ic + (izc - 1) / 2
      d = 2 * id + (izd - 1) / 2
      if(a > sps%n) cycle
      if(b > sps%n) cycle
      if(c > sps%n) cycle
      if(d > sps%n) cycle
      if(sps%nshell(a) + sps%nshell(b) > params%e2max) cycle
      if(sps%nshell(c) + sps%nshell(d) > params%e2max) cycle
      ip = (-1) ** (sps%ll(a) + sps%ll(b))
      itz = (iza + izb) / 2
      ich = two%jptz2n(ij, ip, itz)
      bra = two%jptz(ich)%labels2n(a, b)
      phaseab = two%jptz(ich)%iphase(a, b)
      ket = two%jptz(ich)%labels2n(c, d)
      phasecd = two%jptz(ich)%iphase(c, d)
      if(bra * ket == 0) cycle
      tbme%jptz(ich)%m(bra, ket) = v12 * dble(phaseab * phasecd)
      tbme%jptz(ich)%m(ket, bra) = v12 * dble(phaseab * phasecd)
    end do
    !$omp end do
    !$omp end parallel
    deallocate(aa,bb,cc,dd,jj,az,bz,cz,dz)
    deallocate(v1save)
  end subroutine read_2bme_pn_bin

  subroutine read_2bme_snt_txt(f2, sps, two, tbme)
    type(NBodyScalars), intent(inout) :: tbme
    type(spo_pn), intent(in) :: sps
    type(TwoBodySpace), intent(in) :: two
    character(*), optional, intent(in) :: f2

    integer :: num_p_orb, num_n_orb, num_p_core, num_n_core
    integer :: i, num, num_orb, method
    integer :: n, l, j, itz
    real(8) :: hw_read, t1, v2, zerobdy
    integer :: num_tbme
    integer :: a, b, c, d, ipar, ich
    integer :: bra, ket
    integer :: phase
    integer :: iunit = 19
    open(iunit, file = f2, status = "old")
    call skip_comment(iunit)
    read(iunit,*) num_p_orb,num_n_orb,num_p_core,num_n_core
    num_orb=num_p_orb+num_n_orb
    do i = 1, num_orb
      read(iunit,*)num,n,l,j,itz
    end do
    close(iunit)
    call skip_comment(iunit)
    read(iunit,*) zerobdy
    call skip_comment(iunit)
    read(iunit,*)num, method, hw_read
    call skip_comment(iunit)
    do i = 1, num
      read(iunit,*) a,b,t1
    end do
    call skip_comment(iunit)
    read(iunit,*) num_tbme
    call skip_comment(iunit)
    do i =1, num_tbme
      read(iunit,*)a,b,c,d,j,v2
      if(a > sps%n) cycle
      if(b > sps%n) cycle
      if(c > sps%n) cycle
      if(d > sps%n) cycle
      ipar = (-1) ** (sps%ll(a) + sps%ll(b))
      itz  = (sps%itz(a) + sps%itz(b))/2
      ich = two%jptz2n(j, ipar, itz)
      bra = two%jptz(ich)%labels2n(a,b)
      ket = two%jptz(ich)%labels2n(c,d)
      if(bra*ket==0)cycle
      phase=two%jptz(ich)%iphase(a,b)*two%jptz(ich)%iphase(c,d)
      tbme%jptz(ich)%m(bra, ket) = v2 *dble(phase)
      tbme%jptz(ich)%m(ket, bra) = v2 *dble(phase)
    end do
  end subroutine read_2bme_snt_txt

  subroutine read_2bme_snt_bin(f2, sps, two, tbme)
    use class_stopwatch, only: start_stopwatch, stop_stopwatch, &
        & time_io_read
    type(NBodyScalars), intent(inout) :: tbme
    type(spo_pn), intent(in) :: sps
    type(TwoBodySpace), intent(in) :: two
    character(*), optional, intent(in) :: f2

    integer :: iunit = 19
    integer :: num_p_orb, num_n_orb, num_p_core, num_n_core
    integer :: i, num, num_orb, method
    integer :: j, itz
    real(8) :: hw_read, v2, zerobdy
    integer :: num_tbme, ich
    integer :: a, b, c, d, ipar
    integer :: bra,ket
    integer :: phase
    integer, allocatable :: aa(:), bb(:), cc(:), dd(:)
    real(8), allocatable :: v1save(:), v2save(:)
    integer, allocatable :: nnum(:), nn(:), ll(:), jj(:), zz(:)
    call start_stopwatch(time_io_read)

    open(iunit, form = "unformatted", file = f2, status = "old")
    read(iunit)num_p_orb, num_n_orb, num_p_core, num_n_core
    num_orb = num_p_orb + num_n_orb
    allocate(nnum(num_orb), nn(num_orb), ll(num_orb), jj(num_orb), zz(num_orb))
    read(iunit) nnum
    read(iunit) nn
    read(iunit) ll
    read(iunit) jj
    read(iunit) zz
    deallocate(nnum, nn, ll, jj, zz)

    read(iunit) zerobdy
    read(iunit) num, method, hw_read
    allocate(aa(num), bb(num), v1save(num), v2save(num))
    read(iunit) aa
    read(iunit) bb
    read(iunit) v1save
    read(iunit) v2save
    deallocate(aa, bb, v1save, v2save)

    read(iunit) num_tbme
    allocate(aa(num_tbme), bb(num_tbme), cc(num_tbme), dd(num_tbme), jj(num_tbme))
    allocate(v1save(num_tbme))
    read(iunit) aa
    read(iunit) bb
    read(iunit) cc
    read(iunit) dd
    read(iunit) jj
    read(iunit) v1save
    close(iunit)
    call stop_stopwatch(time_io_read)
    do i =1, num_tbme
      a = aa(i)
      b = bb(i)
      c = cc(i)
      d = dd(i)
      j = jj(i)
      v2 = v1save(i)
      if(a > sps%n) cycle
      if(b > sps%n) cycle
      if(c > sps%n) cycle
      if(d > sps%n) cycle
      ipar = (-1) ** (sps%ll(a)+sps%ll(b))
      itz  = (sps%itz(a) + sps%itz(b))/2
      ich = two%jptz2n(j, ipar, itz)
      bra = two%jptz(ich)%labels2n(a,b)
      ket = two%jptz(ich)%labels2n(c,d)
      if(bra*ket==0)cycle
      phase = two%jptz(ich)%iphase(a,b)*two%jptz(ich)%iphase(c,d)
      tbme%jptz(ich)%m(bra, ket) = v2 *dble(phase)
      tbme%jptz(ich)%m(ket, bra) = v2 *dble(phase)
    end do
    deallocate(aa, bb, cc, dd, jj)
    deallocate(v1save)
  end subroutine read_2bme_snt_bin

  subroutine skip_comment(nfile)
    implicit none
    integer,intent(in)::nfile
    character(1),parameter::com1='!', com2='#'
    character(1)::c1,c2
    read(nfile,'(2a1)') c1, c2
    do while  (c1 == com1 .or. c1 == com2 .or. c2 == com1 .or. c2 == com2)
      read(nfile,'(2a1)') c1, c2
    end do
    backspace(nfile)
  end subroutine skip_comment

  real(8) function Del(i1, i2)
    integer, intent(in) :: i1, i2
    Del = 1.d0
    if(i1 == i2) Del = dsqrt(2.d0)
  end function Del
end module ScalarOperator
