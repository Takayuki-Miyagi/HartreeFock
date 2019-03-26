module TwoBodyOperator
  use omp_lib
  use LinAlgLib
  use TwoBodyModelSpace

  type, extends(DMat) :: TwoBodyPartChannel
    type(TwoBodyChannel), pointer :: ch_bra, ch_ket
    logical :: is = .false.
  contains
    procedure :: InitTwoBodyPartChannel
    procedure :: FinTwoBodyPartChannel
    procedure :: SetTwoBodyPartChannel

    generic :: init => InitTwoBodyPartChannel
    generic :: release => FinTwoBodyPartChannel
    generic :: set => SetTwoBodyPartChannel
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

    generic :: init => InitTwoBodyPart
    generic :: fin => FinTwoBodyPart
    generic :: set => SetTwoBodyPart
  end type TwoBodyPart
contains

  subroutine FinTwoBodyPart(this)
    use MyLibrary, only: triag
    class(TwoBodyPart), intent(inout) :: this
    integer :: chbra, chket

    do chbra = 1, this%Two%NChan
      do chket = 1, this%Two%NChan
        call this%MatCh(chbra,chket)%release()
      end do
    end do
    deallocate(this%MatCh)
    this%Two => null()

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

  subroutine SetTwoBodyPart(this, sps, hw, A, Z, N)
    class(TwoBodyPart), intent(inout) :: this
    type(Orbits), pointer, intent(in) :: sps
    real(8), intent(in) :: hw
    integer, intent(in) :: A, Z, N
    type(TwoBodySpace), pointer :: tbs
    integer :: chbra, chket

    tbs => this%two
    do chbra = 1, tbs%NChan
      do chket = 1, tbs%NChan
        if(.not. this%MatCh(chbra,chket)%is) cycle
        call this%MatCh(chbra,chket)%set(this%oprtr, sps, hw, A, Z, N)
      end do
    end do
  end subroutine SetTwoBodyPart

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
    call this%ini(ch_bra%n_state, ch_ket%n_state)
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
    !$omp do private(bra,ia,ib,ket,ic,id) schedule(dynamic)
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
    type(Orbits), pointer, intent(in) :: sps
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
end module TwoBodyOperator

program test
  use TwoBodyOperator
  type(Orbits), target :: sps
  type(Orbits), pointer :: p_sps
  type(TwoBodySpace) :: two
  type(TwoBodyPart) :: op
  integer :: emax = 4

  p_sps => sps
  call sps%init(emax)
  call two%init(sps, 2*emax)
  call op%init(two, .true., 'Hcm', 0, 1, 0)
  call op%set(sps, 20.d0, 4, 2, 2)

  call op%fin()
  call two%fin()
  call sps%fin()
end program test
