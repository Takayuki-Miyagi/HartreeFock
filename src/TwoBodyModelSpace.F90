module TwoBodyModelSpace
  use omp_lib
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

  type :: TwoBodyChannel
    integer :: j = -1
    integer :: p = 0
    integer :: z = 100
    integer :: n_state = 0
    integer :: n_hp_state = 0 ! hole-particle (not relevant)
    integer :: n_hh_state = 0 ! hole-hole
    integer :: n_hv_state = 0 ! hole-valence
    integer :: n_vv_state = 0 ! valence-valence
    integer :: n_vp_state = 0 ! valence-particle
    integer :: n_pp_state = 0 ! particle-particle
    integer, allocatable :: n2spi1(:)
    integer, allocatable :: n2spi2(:)
    integer, allocatable :: spis2n(:,:)
    integer, allocatable :: iphase(:,:)
    integer, allocatable :: holes(:)
    integer, allocatable :: particles(:)
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
    use MyLibrary, only: triag
    class(TwoBodySpace), intent(inout) :: this
    type(Orbits), target, intent(in) :: sps
    integer, intent(in) :: e2max
    integer :: j, p, z, n, ich
    integer :: i1, i2
    type(SingleParticleOrbit), pointer :: o1, o2
    integer, allocatable :: jj(:), pp(:), zz(:), nn(:)

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
  end subroutine FinTwoBodyChannel

  subroutine InitTwoBodyChannel(this, j, p, z, n, sps, e2max)
    use MyLibrary, only: triag
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
    cnt = 0
    do loop = 1, 9
      cnt_sub = 0
      do i1 = 1, sps%norbs
        o1 => sps%GetOrbit(i1)
        do i2 = 1, i1
          o2 => sps%GetOrbit(i2)

          if(loop == 1 .and. .not. (o1%ph == 0 .and. o2%ph == 1)) cycle ! hp states
          if(loop == 2 .and. .not. (o1%ph == 1 .and. o2%ph == 0)) cycle ! ph states
          if(loop == 3 .and. .not. (o1%ph == 0 .and. o2%ph == 0)) cycle ! hh states
          if(loop == 4 .and. .not. (o1%ph == 0 .and. o2%ph == 2)) cycle ! hv states
          if(loop == 5 .and. .not. (o1%ph == 2 .and. o2%ph == 0)) cycle ! vh states
          if(loop == 6 .and. .not. (o1%ph == 2 .and. o2%ph == 2)) cycle ! vv states
          if(loop == 7 .and. .not. (o1%ph == 2 .and. o2%ph == 1)) cycle ! vp states
          if(loop == 8 .and. .not. (o1%ph == 1 .and. o2%ph == 2)) cycle ! pv states
          if(loop == 9 .and. .not. (o1%ph == 1 .and. o2%ph == 1)) cycle ! pp states

          if(o1%e + o2%e > e2max) cycle
          if(triag(o1%j, o2%j, 2*j)) cycle
          if((-1) ** (o1%l+o2%l) /= p) cycle
          if(o1%z + o2%z /= 2*z) cycle
          if(i1 == i2 .and. mod(j,2) == 1) cycle

          a = i1; b = i2
          if( o1%occ > o2%occ ) then
            ! hp -> ph
            a = i2; b = i1
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
      if(1 <= loop .and. loop <= 2) this%n_hp_state = this%n_hp_state + cnt_sub
      if(loop == 3) this%n_hh_state = cnt_sub
      if(4 <= loop .and. loop <= 5) this%n_hv_state = this%n_hv_state + cnt_sub
      if(loop == 6) this%n_vv_state = cnt_sub
      if(7 <= loop .and. loop <= 8) this%n_vp_state = this%n_vp_state + cnt_sub
      if(loop == 9) this%n_pp_state = cnt_sub
    end do
#ifdef ModelSpaceDebug
    write(*,'(a,i3,a,i3,a,i3,a,i5)') "Two-body channel: J=", j, ", P=", p, ", Tz=", z, ", # of states=", this%n_state
#endif
  end subroutine InitTwoBodyChannel
end module TwoBodyModelSpace
