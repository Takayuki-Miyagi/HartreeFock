module SingleParticleState
  implicit none

  public :: SingleParticleOrbitIsospin
  public :: SingleParticleOrbit
  public :: OrbitsIsospin
  public :: Orbits

  private :: SetSingleParticleOrbitIsospin
  private :: SetSingleParticleOrbit
  private :: InitOrbitsIsospin
  private :: FinOrbitsIsospin
  private :: InitOrbits
  private :: FinOrbits

  type :: SingleParticleOrbitIsospin
    integer :: n = -1
    integer :: l = -1
    integer :: j = -1
    integer :: e = -1
    integer :: idx = -1
  contains
    procedure :: set => SetSingleParticleOrbitIsospin
  end type SingleParticleOrbitIsospin

  type :: SingleParticleOrbit
    integer :: n = -1
    integer :: l = -1
    integer :: j = -1
    integer :: z =  0
    integer :: e = -1
    integer :: idx = -1
  contains
    procedure :: set => SetSingleParticleOrbit
  end type SingleParticleOrbit

  type :: OrbitsIsospin
    integer, allocatable :: nlj2idx(:,:,:)
    type(SingleParticleOrbitIsospin), allocatable :: orb(:)
    logical :: is_constructed=.false.
    integer :: emax, lmax, norbs
  contains
    procedure :: init => InitOrbitsIsospin
    procedure :: fin => FinOrbitsIsospin
  end type OrbitsIsospin

  type :: Orbits
    integer, allocatable :: nljz2idx(:,:,:,:)
    type(SingleParticleOrbit), allocatable :: orb(:)
    logical :: is_constructed=.false.
    integer :: emax, lmax, norbs
  contains
    procedure :: init => InitOrbits
    procedure :: fin => FinOrbits
  end type Orbits

contains

  subroutine FinOrbitsIsospin(this)
    class(OrbitsIsospin), intent(inout) :: this
    if(.not. this%is_constructed) return
#ifdef SingleParticleStateDebug
    write(*,'(a)') 'In FinOrbitsIsospin'
#endif
    deallocate(this%orb)
    deallocate(this%nlj2idx)
  end subroutine FinOrbitsIsospin

  subroutine InitOrbitsIsospin(this, emax, lmax)
    class(OrbitsIsospin), intent(inout) :: this
    integer, intent(in) :: emax
    integer, intent(in), optional :: lmax
    integer :: e, l, n, s, j, cnt

#ifdef SingleParticleStateDebug
    write(*,'(a)') 'In InitOrbitsIsospin:'
#endif
    this%emax = emax
    this%lmax = emax
    if(present(lmax)) this%lmax = lmax
    allocate(this%nlj2idx(0:this%emax/2, 0:this%lmax, 1:2*this%lmax+1))
    this%nlj2idx(:,:,:) = 0
    cnt = 0
    do e = 0, this%emax
      do l = 0, min(e,this%lmax)
        if(mod(e - l, 2) == 1) cycle
        n = (e - l) / 2
        do s = -1, 1, 2
          j = 2*l + s
          if(j < 1) cycle
          cnt = cnt + 1
        end do
      end do
    end do
    this%norbs = cnt
    allocate(this%orb(this%norbs))
    cnt = 0
    do e = 0, this%emax
      do l = 0, min(e,this%lmax)
        if(mod(e - l, 2) == 1) cycle
        n = (e - l) / 2
        do s = -1, 1, 2
          j = 2*l + s
          if(j < 1) cycle
          cnt = cnt + 1
          call this%orb(cnt)%set(n,l,j,cnt)
          this%nlj2idx(n,l,j) = cnt
        end do
      end do
    end do
    this%is_constructed = .true.
  end subroutine InitOrbitsIsospin

  subroutine FinOrbits(this)
    class(Orbits), intent(inout) :: this
    if(.not. this%is_constructed) return
#ifdef SingleParticleStateDebug
    write(*,'(a)') 'In FinOrbits'
#endif
    deallocate(this%orb)
    deallocate(this%nljz2idx)
  end subroutine FinOrbits

  subroutine InitOrbits(this, emax, lmax)
    class(Orbits), intent(inout) :: this
    integer, intent(in) :: emax
    integer, intent(in), optional :: lmax
    integer :: e, l, n, s, j, z, cnt

#ifdef SingleParticleStateDebug
    write(*,'(a)') 'In InitOrbits:'
#endif
    this%emax = emax
    this%lmax = emax
    if(present(lmax)) this%lmax = lmax
    allocate(this%nljz2idx(0:this%emax/2, 0:this%lmax, 1:2*this%lmax+1, -1:1))
    this%nljz2idx(:,:,:,:) = 0
    cnt = 0
    do e = 0, this%emax
      do l = 0, min(e,this%lmax)
        if(mod(e - l, 2) == 1) cycle
        n = (e - l) / 2
        do s = -1, 1, 2
          j = 2*l + s
          if(j < 1) cycle
          do z = -1, 1, 2
            cnt = cnt + 1
          end do
        end do
      end do
    end do
    this%norbs = cnt
    allocate(this%orb(this%norbs))
    cnt = 0
    do e = 0, this%emax
      do l = 0, min(e,this%lmax)
        if(mod(e - l, 2) == 1) cycle
        n = (e - l) / 2
        do s = -1, 1, 2
          j = 2*l + s
          if(j < 1) cycle
          do z = -1, 1, 2
            cnt = cnt + 1
            call this%orb(cnt)%set(n,l,j,z,cnt)
            this%nljz2idx(n,l,j,z) = cnt
          end do
        end do
      end do
    end do
    this%is_constructed = .true.
  end subroutine InitOrbits

  subroutine SetSingleParticleOrbitIsospin(this, n, l, j, idx)
    class(SingleParticleOrbitIsospin), intent(inout) :: this
    integer, intent(in) :: n, l, j, idx
    this%n = n
    this%l = l
    this%j = j
    this%e = 2*n + l
    this%idx = idx
#ifdef SingleParticleStateDebug
    write(*,'(a,i3,a,i3,a,i3,a,i3,a,i3)') &
        & 'index=', idx, ', n=', n, ', l=', l, ', j=', j, ', e=', 2*n+l
#endif
  end subroutine SetSingleParticleOrbitIsospin

  subroutine SetSingleParticleOrbit(this, n, l, j, z, idx)
    class(SingleParticleOrbit), intent(inout) :: this
    integer, intent(in) :: n, l, j, z, idx
    this%n = n
    this%l = l
    this%j = j
    this%z = z
    this%e = 2*n + l
    this%idx = idx
#ifdef SingleParticleStateDebug
    write(*,'(a,i3,a,i3,a,i3,a,i3,a,i3,a,i3)') &
        & 'indx=',idx,', n=', n, ', l=', l, ', j=', j, ', z=', z, ', e=', 2*n+l
#endif
  end subroutine SetSingleParticleOrbit
end module SingleParticleState

! main program for check
!program test
!  use SingleParticleState
!  type(Orbits) :: o
!  type(OrbitsIsospin) :: io
!  call o%init(4)
!  call o%fin()
!
!  call o%init(4,2)
!  call o%fin()
!
!  call io%init(4)
!  call io%fin()
!
!  call io%init(4,2)
!  call io%fin()
!end program test


