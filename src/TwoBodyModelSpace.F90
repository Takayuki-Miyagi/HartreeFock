module TwoBodyModelSpace
  use omp_lib
  use myfort
  use SingleParticleState
  implicit none

  public :: TwoBodyChannel
  public :: TwoBodySpace

  private :: InitTwoBodySpace
  private :: FinTwoBodySpace
  private :: GetTwoBodyChannelFromJPZ
  private :: GetTwoBodyChannelFromCh
  private :: InitTwoBodyChannel
  private :: FinTwoBodyChannel

  private :: InitCrossCoupledTwoBodyChannel
  private :: FinCrossCoupledTwoBodyChannel
  private :: InitCrossCoupledTwoBodySpace
  private :: FinCrossCoupledTwoBodySpace
  private :: GetCrossCoupledTwoBodyChannelFromJPZ
  private :: GetCrossCoupledTwoBodyChannelFromCh

  type :: TwoBodyChannel
    integer :: j = -1
    integer :: p = 0
    integer :: z = 100
    integer :: n_state = 0
    integer :: n_co_state = 0 ! core-outsidee (not relevant)
    integer :: n_cc_state = 0 ! core-core
    integer :: n_cv_state = 0 ! core-valence
    integer :: n_vv_state = 0 ! valence-valence
    integer :: n_vo_state = 0 ! valence-outside
    integer :: n_oo_state = 0 ! outside-outside
    integer :: n_hh_state = 0 ! hole-hole
    integer :: n_ph_state = 0 ! particle-hole
    integer :: n_pp_state = 0 ! particle-particle
    integer, allocatable :: n2spi1(:)
    integer, allocatable :: n2spi2(:)
    integer, allocatable :: spis2n(:,:)
    integer, allocatable :: iphase(:,:)
    integer, allocatable :: co_s(:)
    integer, allocatable :: cc_s(:)
    integer, allocatable :: cv_s(:)
    integer, allocatable :: vv_s(:)
    integer, allocatable :: vo_s(:)
    integer, allocatable :: oo_s(:)
    integer, allocatable :: hh_s(:)
    integer, allocatable :: ph_s(:)
    integer, allocatable :: pp_s(:)
  contains
    procedure :: InitTwoBodyChannel
    procedure :: FinTwoBodyChannel
    generic :: init => InitTwoBodyChannel
    generic :: fin => FinTwoBodyChannel
  end type TwoBodyChannel

  type :: TwoBodySpace
    type(TwoBodyChannel), allocatable :: jpz(:)
    type(Orbits), pointer :: sps
    integer, allocatable :: jpz2ch(:,:,:)
    integer :: emax, e2max
    integer :: NChan
  contains
    procedure :: InitTwoBodySpace
    procedure :: FinTwoBodySpace
    procedure :: GetTwoBodyChannelFromJPZ
    procedure :: GetTwoBodyChannelFromCh

    generic :: init => InitTwoBodySpace
    generic :: fin => FinTwoBodySpace
    generic :: GetTwoBodyChannel => &
        & GetTwoBodyChannelFromJPZ, &
        & GetTwoBodyChannelFromCh
  end type TwoBodySpace

  type :: CrossCoupledTwoBodyChannel
    ! This is for
    ! < ab:J | V | cd:J > -> < a\bar{c}: J' | V | b\bar{d}: J'>
    !                     or < a\bar{d}: J' | V | b\bar{c}: J' >
    ! Storing |ac: J'>, not actual two-body state (no Pauli principle)
    ! < pp | V | pp > -> < pp | V | pp >
    ! < nn | V | nn > -> < nn | V | nn >
    ! < pn | V | pn > -> < pn | V | pn > or < pp | V | nn >
    ! So, Tz is not good number any more, but |Tz| is still good.
    !
    integer :: j = -1
    integer :: p = 0
    integer :: z_abs = 100
    integer :: n_state = 0
    integer, allocatable :: n2spi1(:)
    integer, allocatable :: n2spi2(:)
    integer, allocatable :: spis2n(:,:)
    integer, allocatable :: iphase(:,:)
  contains
    procedure :: InitCrossCoupledTwoBodyChannel
    procedure :: FinCrossCoupledTwoBodyChannel
    generic :: init => InitCrossCoupledTwoBodyChannel
    generic :: fin => FinCrossCoupledTwoBodyChannel
  end type CrossCoupledTwoBodyChannel

  type :: CrossCoupledTwoBodySpace
    type(CrossCoupledTwoBodyChannel), allocatable :: jpz(:)
    type(Orbits), pointer :: sps
    integer, allocatable :: jpz2ch(:,:,:)
    integer :: emax, e2max
    integer :: NChan
  contains
    procedure :: InitCrossCoupledTwoBodySpace
    procedure :: FinCrossCoupledTwoBodySpace
    procedure :: GetCrossCoupledTwoBodyChannelFromJPZ
    procedure :: GetCrossCoupledTwoBodyChannelFromCh

    generic :: init => InitCrossCoupledTwoBodySpace
    generic :: fin => FinCrossCoupledTwoBodySpace
    generic :: GetCrossCoupledTwoBodyChannel => &
        & GetCrossCoupledTwoBodyChannelFromJPZ, &
        & GetCrossCoupledTwoBodyChannelFromCh
  end type CrossCoupledTwoBodySpace
contains
  subroutine FinTwoBodySpace(this)
    class(TwoBodySpace), intent(inout) :: this
    integer :: ch
    do ch = 1, this%NChan
      call this%jpz(ch)%fin()
    end do
    deallocate(this%jpz)
    deallocate(this%jpz2ch)
  end subroutine FinTwoBodySpace

  subroutine InitTwoBodySpace(this, sps, e2max)
    class(TwoBodySpace), intent(inout) :: this
    type(Orbits), target, intent(in) :: sps
    integer, intent(in) :: e2max
    integer :: j, p, z, n, ich
    integer :: i1, i2
    type(SingleParticleOrbit), pointer :: o1, o2
    integer, allocatable :: jj(:), pp(:), zz(:), nn(:)

    if(allocated(this%jpz)) call this%fin()

    this%sps => sps
    allocate(this%jpz2ch(0:min(2*sps%lmax,e2max)+1,-1:1,-1:1))
    this%jpz2ch(:,:,:) = 0
    ich = 0
    do j = 0, min(2*sps%lmax,e2max)+1
      do p = 1, -1, -2
        do z = -1, 1
          n = 0

          do i1 = 1, sps%norbs
            o1 => sps%GetOrbit(i1)
            do i2 = 1, i1
              o2 => sps%GetOrbit(i2)

              if(o1%e + o2%e > e2max) cycle
              if(triag(o1%j, o2%j, 2*j)) cycle
              if((-1) ** (o1%l+o2%l) /= p) cycle
              if(o1%z + o2%z /= 2*z) cycle
              if(i1 == i2 .and. mod(j,2) == 1) cycle
              n = n + 1

            end do
          end do
          if(n /= 0) then
            ich = ich + 1
            this%jpz2ch(j,p,z) = ich
          end if
        end do
      end do
    end do

    this%NChan = ich
    this%emax = sps%emax
    this%e2max = e2max
    allocate(this%jpz(this%NChan))
    allocate(jj(this%NChan))
    allocate(pp(this%NChan))
    allocate(zz(this%NChan))
    allocate(nn(this%NChan))

    ich = 0
    do j = 0, min(2*sps%lmax,e2max)+1
      do p = 1, -1, -2
        do z = -1, 1
          n = 0

          do i1 = 1, sps%norbs
            o1 => sps%GetOrbit(i1)
            do i2 = 1, i1
              o2 => sps%GetOrbit(i2)

              if(o1%e + o2%e > e2max) cycle
              if(triag(o1%j, o2%j, 2*j)) cycle
              if((-1) ** (o1%l+o2%l) /= p) cycle
              if(o1%z + o2%z /= 2*z) cycle
              if(i1 == i2 .and. mod(j,2) == 1) cycle
              n = n + 1

            end do
          end do
          if(n /= 0) then
            ich = ich + 1
            jj(ich) = j
            pp(ich) = p
            zz(ich) = z
            nn(ich) = n
          end if
        end do
      end do
    end do

    do ich = 1, this%NChan
      call this%jpz(ich)%init(jj(ich),pp(ich),zz(ich),nn(ich),sps,e2max)
    end do
    deallocate(jj,pp,zz,nn)
  end subroutine InitTwoBodySpace

  function GetTwoBodyChannelFromJPZ(this, j, p, z) result(ptr)
    class(TwoBodySpace), intent(in) :: this
    type(TwoBodyChannel), pointer :: ptr
    integer, intent(in) :: j, p, z
    integer :: ch
    ptr => null()
    ch = this%jpz2ch(j,p,z)
    if(ch == 0) return
    ptr => this%GetTwoBodyChannel(ch)
  end function GetTwoBodyChannelFromJPZ

  function GetTwoBodyChannelFromCh(this, ch) result(ptr)
    class(TwoBodySpace), target, intent(in) :: this
    type(TwoBodyChannel), pointer :: ptr
    integer, intent(in) :: ch
    ptr => this%jpz(ch)
  end function GetTwoBodyChannelFromCh

  subroutine FinTwoBodyChannel(this)
    class(TwoBodyChannel), intent(inout) :: this
    deallocate(this%n2spi1)
    deallocate(this%n2spi2)
    deallocate(this%spis2n)
    deallocate(this%iphase)
    deallocate(this%co_s)
    deallocate(this%cc_s)
    deallocate(this%cv_s)
    deallocate(this%vv_s)
    deallocate(this%vo_s)
    deallocate(this%oo_s)
    deallocate(this%hh_s)
    deallocate(this%ph_s)
    deallocate(this%pp_s)
  end subroutine FinTwoBodyChannel

  subroutine InitTwoBodyChannel(this, j, p, z, n, sps, e2max)
    class(TwoBodyChannel), intent(inout) :: this
    integer, intent(in) :: j, p, z, n, e2max
    type(Orbits), intent(in), target :: sps
    integer :: i1, i2
    type(SingleParticleOrbit), pointer :: o1, o2
    integer :: a, b, cnt, loop, cnt_sub

    this%j = j
    this%p = p
    this%z = z
    this%n_state = n
    allocate(this%n2spi1(n))
    allocate(this%n2spi2(n))
    allocate(this%spis2n(sps%norbs,sps%norbs))
    allocate(this%iphase(sps%norbs,sps%norbs))
    this%spis2n(:,:) = 0
    this%iphase(:,:) = 0
    cnt = 0
    do loop = 1, 9
      cnt_sub = 0
      do i1 = 1, sps%norbs
        o1 => sps%GetOrbit(i1)
        do i2 = 1, i1
          o2 => sps%GetOrbit(i2)

          if(loop == 1 .and. .not. (o1%GetCoreValenceOutside() == 0 .and. o2%GetCoreValenceOutside() == 2)) cycle ! co states
          if(loop == 2 .and. .not. (o1%GetCoreValenceOutside() == 2 .and. o2%GetCoreValenceOutside() == 0)) cycle ! oc states
          if(loop == 3 .and. .not. (o1%GetCoreValenceOutside() == 0 .and. o2%GetCoreValenceOutside() == 0)) cycle ! cc states
          if(loop == 4 .and. .not. (o1%GetCoreValenceOutside() == 0 .and. o2%GetCoreValenceOutside() == 1)) cycle ! cv states
          if(loop == 5 .and. .not. (o1%GetCoreValenceOutside() == 1 .and. o2%GetCoreValenceOutside() == 0)) cycle ! vc states
          if(loop == 6 .and. .not. (o1%GetCoreValenceOutside() == 1 .and. o2%GetCoreValenceOutside() == 1)) cycle ! vv states
          if(loop == 7 .and. .not. (o1%GetCoreValenceOutside() == 2 .and. o2%GetCoreValenceOutside() == 1)) cycle ! vo states
          if(loop == 8 .and. .not. (o1%GetCoreValenceOutside() == 1 .and. o2%GetCoreValenceOutside() == 2)) cycle ! ov states
          if(loop == 9 .and. .not. (o1%GetCoreValenceOutside() == 2 .and. o2%GetCoreValenceOutside() == 2)) cycle ! oo states

          if(o1%e + o2%e > e2max) cycle
          if(triag(o1%j, o2%j, 2*j)) cycle
          if((-1) ** (o1%l+o2%l) /= p) cycle
          if(o1%z + o2%z /= 2*z) cycle
          if(i1 == i2 .and. mod(j,2) == 1) cycle
          a = i1; b = i2
          if( o1%GetOccupation() > o2%GetOccupation() ) then
            a = i2; b = i1 ! hp -> ph
          end if
          cnt = cnt + 1
          cnt_sub = cnt_sub + 1
          this%n2spi1(cnt) = a
          this%n2spi2(cnt) = b
          this%spis2n(a,b) = cnt
          this%spis2n(b,a) = cnt
          this%iphase(a,b) = 1
          this%iphase(b,a) = -(-1) ** ((o1%j+o2%j)/2 - j)
        end do
      end do
      if(1 <= loop .and. loop <= 2) this%n_co_state = this%n_co_state + cnt_sub
      if(loop == 3) this%n_cc_state = cnt_sub
      if(4 <= loop .and. loop <= 5) this%n_cv_state = this%n_cv_state + cnt_sub
      if(loop == 6) this%n_vv_state = cnt_sub
      if(7 <= loop .and. loop <= 8) this%n_vo_state = this%n_vo_state + cnt_sub
      if(loop == 9) this%n_oo_state = cnt_sub
    end do

    allocate(this%co_s(this%n_co_state))
    allocate(this%cc_s(this%n_cc_state))
    allocate(this%cv_s(this%n_cv_state))
    allocate(this%vv_s(this%n_vv_state))
    allocate(this%vo_s(this%n_vo_state))
    allocate(this%oo_s(this%n_oo_state))

    cnt = 0
    cnt_sub = 0
    do loop = 1, 9
      if(loop==3) cnt_sub = 0
      if(loop==4) cnt_sub = 0
      if(loop==6) cnt_sub = 0
      if(loop==7) cnt_sub = 0
      if(loop==9) cnt_sub = 0
      do i1 = 1, sps%norbs
        o1 => sps%GetOrbit(i1)
        do i2 = 1, i1
          o2 => sps%GetOrbit(i2)

          if(loop == 1 .and. .not. (o1%GetCoreValenceOutside() == 0 .and. o2%GetCoreValenceOutside() == 2)) cycle ! co states
          if(loop == 2 .and. .not. (o1%GetCoreValenceOutside() == 2 .and. o2%GetCoreValenceOutside() == 0)) cycle ! oc states
          if(loop == 3 .and. .not. (o1%GetCoreValenceOutside() == 0 .and. o2%GetCoreValenceOutside() == 0)) cycle ! cc states
          if(loop == 4 .and. .not. (o1%GetCoreValenceOutside() == 0 .and. o2%GetCoreValenceOutside() == 1)) cycle ! cv states
          if(loop == 5 .and. .not. (o1%GetCoreValenceOutside() == 1 .and. o2%GetCoreValenceOutside() == 0)) cycle ! vc states
          if(loop == 6 .and. .not. (o1%GetCoreValenceOutside() == 1 .and. o2%GetCoreValenceOutside() == 1)) cycle ! vv states
          if(loop == 7 .and. .not. (o1%GetCoreValenceOutside() == 2 .and. o2%GetCoreValenceOutside() == 1)) cycle ! vo states
          if(loop == 8 .and. .not. (o1%GetCoreValenceOutside() == 1 .and. o2%GetCoreValenceOutside() == 2)) cycle ! ov states
          if(loop == 9 .and. .not. (o1%GetCoreValenceOutside() == 2 .and. o2%GetCoreValenceOutside() == 2)) cycle ! oo states

          if(o1%e + o2%e > e2max) cycle
          if(triag(o1%j, o2%j, 2*j)) cycle
          if((-1) ** (o1%l+o2%l) /= p) cycle
          if(o1%z + o2%z /= 2*z) cycle
          if(i1 == i2 .and. mod(j,2) == 1) cycle

          cnt = cnt + 1
          cnt_sub = cnt_sub + 1
          if(1 <= loop .and. loop <= 2) this%co_s(cnt_sub) = cnt
          if(loop == 3) this%cc_s(cnt_sub) = cnt
          if(4 <= loop .and. loop <= 5) this%cv_s(cnt_sub) = cnt
          if(6 == loop) this%vv_s(cnt_sub) = cnt
          if(7 <= loop .and. loop <= 8) this%vo_s(cnt_sub) = cnt
          if(loop == 9) this%oo_s(cnt_sub) = cnt
        end do
      end do
    end do

    cnt = 0
    do loop = 1, 3
      cnt_sub = 0
      do i1 = 1, sps%norbs
        o1 => sps%GetOrbit(i1)
        do i2 = 1, i1
          o2 => sps%GetOrbit(i2)

          if(loop == 1 .and. o1%GetHoleParticle() + o2%GetHoleParticle() /= 0) cycle ! hh states
          if(loop == 2 .and. o1%GetHoleParticle() + o2%GetHoleParticle() /= 1) cycle ! ph states
          if(loop == 3 .and. o1%GetHoleParticle() + o2%GetHoleParticle() /= 2) cycle ! pp states

          if(o1%e + o2%e > e2max) cycle
          if(triag(o1%j, o2%j, 2*j)) cycle
          if((-1) ** (o1%l+o2%l) /= p) cycle
          if(o1%z + o2%z /= 2*z) cycle
          if(i1 == i2 .and. mod(j,2) == 1) cycle

          cnt_sub = cnt_sub + 1
        end do
      end do
      if(loop==1) this%n_hh_state = cnt_sub
      if(loop==2) this%n_ph_state = cnt_sub
      if(loop==3) this%n_pp_state = cnt_sub
    end do
    allocate(this%hh_s(this%n_hh_state))
    allocate(this%ph_s(this%n_ph_state))
    allocate(this%pp_s(this%n_pp_state))

    do loop = 1, 3
      cnt_sub = 0
      do i1 = 1, sps%norbs
        o1 => sps%GetOrbit(i1)
        do i2 = 1, i1
          o2 => sps%GetOrbit(i2)

          if(loop == 1 .and. o1%GetHoleParticle() + o2%GetHoleParticle() /= 0) cycle ! hh states
          if(loop == 2 .and. o1%GetHoleParticle() + o2%GetHoleParticle() /= 1) cycle ! ph states
          if(loop == 3 .and. o1%GetHoleParticle() + o2%GetHoleParticle() /= 2) cycle ! pp states

          if(o1%e + o2%e > e2max) cycle
          if(triag(o1%j, o2%j, 2*j)) cycle
          if((-1) ** (o1%l+o2%l) /= p) cycle
          if(o1%z + o2%z /= 2*z) cycle
          if(i1 == i2 .and. mod(j,2) == 1) cycle

          cnt_sub = cnt_sub + 1
          if(loop==1) this%hh_s(cnt_sub) = this%spis2n(i1,i2)
          if(loop==2) this%ph_s(cnt_sub) = this%spis2n(i1,i2)
          if(loop==3) this%pp_s(cnt_sub) = this%spis2n(i1,i2)
        end do
      end do
    end do

#ifdef ModelSpaceDebug
    write(*,'(a,i3,a,i3,a,i3,a,i5)') "Two-body channel: J=", j, ", P=", p, ", Tz=", z, ", # of states=", this%n_state
#endif
  end subroutine InitTwoBodyChannel

  subroutine FinCrossCoupledTwoBodySpace(this)
    class(CrossCoupledTwoBodySpace), intent(inout) :: this
    integer :: ch
    do ch = 1, this%NChan
      call this%jpz(ch)%fin()
    end do
    deallocate(this%jpz)
    deallocate(this%jpz2ch)
  end subroutine FinCrossCoupledTwoBodySpace

  subroutine InitCrossCoupledTwoBodySpace(this, sps, e2max)
    class(CrossCoupledTwoBodySpace), intent(inout) :: this
    type(Orbits), target, intent(in) :: sps
    integer, intent(in) :: e2max
    integer :: j, p, z, n, ich
    integer :: i1, i2
    type(SingleParticleOrbit), pointer :: o1, o2
    integer, allocatable :: jj(:), pp(:), zz(:), nn(:)

    if(allocated(this%jpz)) call this%fin()

    this%sps => sps
    allocate(this%jpz2ch(0:min(2*sps%lmax,e2max)+1,-1:1,0:1))
    this%jpz2ch(:,:,:) = 0
    ich = 0
    do j = 0, min(2*sps%lmax,e2max)+1
      do p = 1, -1, -2
        do z = 0, 1
          n = 0

          do i1 = 1, sps%norbs
            o1 => sps%GetOrbit(i1)
            do i2 = 1, i1
              o2 => sps%GetOrbit(i2)

              if(o1%e + o2%e > e2max) cycle
              if(triag(o1%j, o2%j, 2*j)) cycle
              if((-1) ** (o1%l+o2%l) /= p) cycle
              if(abs(o1%z + o2%z) /= 2*z) cycle
              n = n + 1

            end do
          end do
          if(n /= 0) then
            ich = ich + 1
            this%jpz2ch(j,p,z) = ich
          end if
        end do
      end do
    end do

    this%NChan = ich
    this%emax = sps%emax
    this%e2max = e2max
    allocate(this%jpz(this%NChan))
    allocate(jj(this%NChan))
    allocate(pp(this%NChan))
    allocate(zz(this%NChan))
    allocate(nn(this%NChan))

    ich = 0
    do j = 0, min(2*sps%lmax,e2max)+1
      do p = 1, -1, -2
        do z = 0, 1
          n = 0

          do i1 = 1, sps%norbs
            o1 => sps%GetOrbit(i1)
            do i2 = 1, i1
              o2 => sps%GetOrbit(i2)

              if(o1%e + o2%e > e2max) cycle
              if(triag(o1%j, o2%j, 2*j)) cycle
              if((-1) ** (o1%l+o2%l) /= p) cycle
              if(abs(o1%z + o2%z) /= 2*z) cycle
              n = n + 1

            end do
          end do
          if(n /= 0) then
            ich = ich + 1
            jj(ich) = j
            pp(ich) = p
            zz(ich) = z
            nn(ich) = n
          end if
        end do
      end do
    end do

    do ich = 1, this%NChan
      call this%jpz(ich)%init(jj(ich),pp(ich),zz(ich),nn(ich),sps,e2max)
    end do
    deallocate(jj,pp,zz,nn)
  end subroutine InitCrossCoupledTwoBodySpace

  function GetCrossCoupledTwoBodyChannelFromJPZ(this, j, p, z) result(ptr)
    class(CrossCoupledTwoBodySpace), intent(in) :: this
    type(CrossCoupledTwoBodyChannel), pointer :: ptr
    integer, intent(in) :: j, p, z
    integer :: ch
    ptr => null()
    ch = this%jpz2ch(j,p,z)
    if(ch == 0) return
    ptr => this%GetCrossCoupledTwoBodyChannel(ch)
  end function GetCrossCoupledTwoBodyChannelFromJPZ

  function GetCrossCoupledTwoBodyChannelFromCh(this, ch) result(ptr)
    class(CrossCoupledTwoBodySpace), target, intent(in) :: this
    type(CrossCoupledTwoBodyChannel), pointer :: ptr
    integer, intent(in) :: ch
    ptr => this%jpz(ch)
  end function GetCrossCoupledTwoBodyChannelFromCh

  subroutine FinCrossCoupledTwoBodyChannel(this)
    class(CrossCoupledTwoBodyChannel), intent(inout) :: this
    deallocate(this%n2spi1)
    deallocate(this%n2spi2)
    deallocate(this%spis2n)
    deallocate(this%iphase)
  end subroutine FinCrossCoupledTwoBodyChannel

  subroutine InitCrossCoupledTwoBodyChannel(this, j, p, z, n, sps, e2max)
    class(CrossCoupledTwoBodyChannel), intent(inout) :: this
    integer, intent(in) :: j, p, z, n, e2max
    type(Orbits), intent(in), target :: sps
    integer :: i1, i2
    type(SingleParticleOrbit), pointer :: o1, o2
    integer :: a, b, cnt

    this%j = j
    this%p = p
    this%z_abs = z
    this%n_state = n
    allocate(this%n2spi1(n))
    allocate(this%n2spi2(n))
    allocate(this%spis2n(sps%norbs,sps%norbs))
    allocate(this%iphase(sps%norbs,sps%norbs))
    cnt = 0
    do i1 = 1, sps%norbs
      o1 => sps%GetOrbit(i1)
      do i2 = 1, i1
        o2 => sps%GetOrbit(i2)

        if(o1%e + o2%e > e2max) cycle
        if(triag(o1%j, o2%j, 2*j)) cycle
        if((-1) ** (o1%l+o2%l) /= p) cycle
        if(abs(o1%z + o2%z) /= 2*z) cycle

        a = i1; b = i2
        if( o1%GetOccupation() > o2%GetOccupation() ) then
          a = i2; b = i1 ! hp -> ph
        end if
        cnt = cnt + 1
        this%n2spi1(cnt) = a
        this%n2spi2(cnt) = b
        this%spis2n(a,b) = cnt
        this%spis2n(b,a) = cnt
        this%iphase(a,b) = 1
        this%iphase(b,a) = -(-1) ** ((o1%j+o2%j)/2 - j)
      end do
    end do
#ifdef ModelSpaceDebug
    write(*,'(a,i3,a,i3,a,i3,a,i5)') "Cross coupled two-body channel: J=", &
        & j, ", P=", p, ", |Tz|=", z, ", # of states=", this%n_state
#endif
  end subroutine InitCrossCoupledTwoBodyChannel

end module TwoBodyModelSpace
