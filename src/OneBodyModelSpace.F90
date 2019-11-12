module OneBodyModelSpace
  use omp_lib
  use SingleParticleState
  implicit none

  public :: OneBodySpace
  public :: OneBodyChannel

  private :: InitOneBodyChannel
  private :: FinOneBodyChannel
  private :: InitOneBodySpace
  private :: FinOneBodySpace
  private :: GetOneBodyChannelFromJPZ
  private :: GetOneBodyChannelFromCh

  type :: OneBodyChannel
    integer :: j = -1
    integer :: p = 0
    integer :: z = 100
    integer :: n_state = 0
    integer :: n_h_state = 0
    integer :: n_p_state = 0
    integer :: n_v_state = 0
    integer, allocatable :: n2spi(:)
    integer, allocatable :: spi2n(:)
  contains
    procedure :: InitOneBodyChannel
    procedure :: FinOneBodyChannel
    generic :: init => InitOneBodyChannel
    generic :: fin => FinOneBodyChannel
  end type OneBodyChannel

  type :: OneBodySpace
    type(OneBodyChannel), allocatable :: jpz(:)
    type(Orbits), pointer :: sps
    integer, allocatable :: jpz2ch(:,:,:)
    integer :: NChan
    integer :: emax
  contains
    procedure :: InitOneBodySpace
    procedure :: FinOneBodySpace
    procedure :: GetOneBodyChannelFromJPZ
    procedure :: GetOneBodyChannelFromCh

    generic :: init => InitOneBodySpace
    generic :: fin => FinOneBodySpace
    generic :: GetOneBodyChannel => &
        & GetOneBodyChannelFromJPZ, &
        & GetOneBodyChannelFromCh
  end type OneBodySpace
contains
  subroutine FinOneBodySpace(this)
    class(OneBodySpace), intent(inout) :: this
    integer :: ich
    do ich = 1, this%NChan
      call this%jpz(ich)%fin()
    end do
    deallocate(this%jpz)
    deallocate(this%jpz2ch)
  end subroutine FinOneBodySpace

  subroutine InitOneBodySpace(this, sps)
    class(OneBodySpace), intent(inout) :: this
    type(Orbits), target, intent(in) :: sps
    integer :: ich
    integer :: j, p, z, n, i
    integer, allocatable :: jj(:), pp(:), zz(:), nn(:)

    if(allocated(this%jpz)) call this%fin()

    this%sps => sps
    ich = 0
    allocate(this%jpz2ch(1:2*sps%lmax+1,-1:1,-1:1))
    this%jpz2ch(:,:,:) = 0
    do j = 1, 2*sps%lmax+1, 2
      do p = 1, -1, -2
        do z = -1, 1, 2
          n = 0
          do i = 1, sps%norbs
            if(sps%orb(i)%j /= j) cycle
            if((-1)**sps%orb(i)%l /= p) cycle
            if(sps%orb(i)%z /= z) cycle
            n = n + 1
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
    allocate(this%jpz(this%NChan))
    allocate(jj(this%NChan))
    allocate(pp(this%NChan))
    allocate(zz(this%NChan))
    allocate(nn(this%NChan))

    ich = 0
    do j = 1, 2*sps%lmax+1, 2
      do p = 1, -1, -2
        do z = -1, 1, 2
          n = 0
          do i = 1, sps%norbs
            if(sps%orb(i)%j /= j) cycle
            if((-1)**sps%orb(i)%l /= p) cycle
            if(sps%orb(i)%z /= z) cycle
            n = n + 1
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
      call this%jpz(ich)%init(jj(ich), pp(ich), zz(ich), nn(ich), sps)
    end do
  end subroutine InitOneBodySpace

  function GetOneBodyChannelFromJPZ(this, j, p, z) result(ptr)
    class(OneBodySpace), intent(in) :: this
    type(OneBodyChannel), pointer :: ptr
    integer, intent(in) :: j, p, z
    integer :: ch
    ptr => null()
    ch = this%jpz2ch(j,p,z)
    if(ch == 0) return
    ptr => this%GetOneBodyChannel(ch)
  end function GetOneBodyChannelFromJPZ

  function GetOneBodyChannelFromCh(this, ch) result(ptr)
    class(OneBodySpace), target, intent(in) :: this
    type(OneBodyChannel), pointer :: ptr
    integer, intent(in) :: ch
    ptr => this%jpz(ch)
  end function GetOneBodyChannelFromCh

  subroutine FinOneBodyChannel(this)
    class(OneBodyChannel), intent(inout) :: this
    deallocate(this%n2spi)
    deallocate(this%spi2n)
  end subroutine FinOneBodyChannel

  subroutine InitOneBodyChannel(this, j, p, z, n, sps)
    class(OneBodyChannel), intent(inout) :: this
    type(Orbits), intent(in), target :: sps
    type(SingleParticleOrbit), pointer :: o
    integer, intent(in) :: j, p, z, n
    integer :: cnt, i, loop, cnt_sub
    integer :: partition(3) = [0,1,2]
    this%j = j
    this%p = p
    this%z = z
    this%n_state = n
    allocate(this%n2spi(n))
    allocate(this%spi2n(sps%norbs))
    this%spi2n(:) = 0
    cnt = 0
    do loop = 1, size(partition)
      cnt_sub = 0
      do i = 1, sps%norbs
        o => sps%GetOrbit(i)
        if(o%GetCoreValenceOutside() /= partition(loop)) cycle
        if(o%j /= j) cycle
        if((-1)**o%l /= p) cycle
        if(o%z /= z) cycle
        cnt = cnt + 1
        cnt_sub = cnt_sub + 1
        this%n2spi(cnt) = i
        this%spi2n(i) = cnt
      end do
      if(partition(loop) == 0) this%n_h_state = cnt_sub
      if(partition(loop) == 1) this%n_p_state = cnt_sub
      if(partition(loop) == 2) this%n_v_state = cnt_sub
    end do
#ifdef ModelSpaceDebug
    write(*,'(a,i3,a,i3,a,i3,a,i5)') "One-body channel: J=", j, ", P=", p, ", Tz=", z, ", # of states=", this%n_state
#endif
  end subroutine InitOneBodyChannel
end module OneBodyModelSpace
