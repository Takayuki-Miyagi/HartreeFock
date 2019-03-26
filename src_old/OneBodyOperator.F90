module OneBodyOperator
  use omp_lib
  use LinAlgLib
  use OneBodyModelSpace

  implicit none
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

    generic :: init => InitOneBodyPart
    generic :: fin => FinOneBodyPart
    generic :: set => SetOneBodyPart
  end type OneBodyPart
contains

  subroutine FinOneBodyPart(this)
    use MyLibrary, only: triag
    class(OneBodyPart), intent(inout) :: this
    integer :: chbra, chket

    do chbra = 1, this%one%NChan
      do chket = 1, this%one%NChan
        call this%MatCh(chbra,chket)%release()
      end do
    end do
    deallocate(this%MatCh)
    this%one => null()

  end subroutine FinOneBodyPart

  subroutine InitOneBodyPart(this, one, Scalar, oprtr, jr, pr, zr)
    use MyLibrary, only: triag
    class(OneBodyPart), intent(inout) :: this
    type(OneBodySpace), target, intent(in) :: one
    logical, intent(in) :: Scalar
    character(*), intent(in) :: oprtr
    integer, intent(in) :: jr, pr, zr
    integer :: chbra, chket
    integer :: jbra, pbra, zbra, nbra
    integer :: jket, pket, zket, nket
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

  subroutine SetOneBodyPart(this, sps, hw, A, Z, N)
    class(OneBodyPart), intent(inout) :: this
    type(Orbits), pointer :: sps
    type(OneBodySpace), pointer :: obs
    real(8), intent(in) :: hw
    integer, intent(in) :: A, Z, N
    integer :: chbra
    integer :: chket
    obs => this%one
    do chbra = 1, obs%NChan
      do chket = 1, obs%NChan
        if(.not. this%MatCh(chbra,chket)%is) cycle
        call this%MatCh(chbra,chket)%set(this%oprtr, sps, hw, A, Z, N)
      end do
    end do
  end subroutine SetOneBodyPart

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
    call this%ini(ch_bra%n_state, ch_ket%n_state)
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
end module OneBodyOperator
