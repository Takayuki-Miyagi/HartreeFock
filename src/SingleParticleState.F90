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
  private :: iso2pn
  private :: pn2iso

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
    real(8) :: occ = -1.d0
    integer :: ph = -1 ! 0:hole,1:particle,2:partial occ
  contains
    procedure :: set => SetSingleParticleOrbit
    procedure :: SetOccupation
    procedure :: SetParticleHole
  end type SingleParticleOrbit

  type :: OrbitsIsospin
    integer, allocatable :: nlj2idx(:,:,:)
    type(SingleParticleOrbitIsospin), allocatable :: orb(:)
    logical :: is_constructed=.false.
    integer :: emax, lmax, norbs
  contains
    procedure :: init => InitOrbitsIsospin
    procedure :: fin => FinOrbitsIsospin
    procedure :: iso2pn
  end type OrbitsIsospin

  type :: Orbits
    integer, allocatable :: nljz2idx(:,:,:,:)
    type(SingleParticleOrbit), allocatable :: orb(:)
    logical :: is_constructed=.false.
    integer :: emax, lmax, norbs
  contains
    procedure :: init => InitOrbits
    procedure :: fin => FinOrbits
    procedure :: pn2iso
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

  subroutine SetOccupation(this,occ)
    class(SingleParticleOrbit), intent(inout) :: this
    real(8), intent(in) :: occ
    this%occ = occ
  end subroutine SetOccupation

  subroutine SetParticleHole(this,ph_label)
    class(SingleParticleOrbit), intent(inout) :: this
    integer, intent(in) :: ph_label
    this%ph = ph_label
  end subroutine SetParticleHole

  function iso2pn(this, sps, idx, z) result(r)
    class(OrbitsIsospin), intent(in) :: this
    type(Orbits), intent(in) :: sps
    integer, intent(in) :: idx, z
    integer :: r

    r = 0
    if(z /= -1 .and. z /= 1) then
      write(*,'(a,i3)') "Error in iso2pn, tz has to be -1 or 1. tz = ", z
      stop
    end if
    if(this%orb(idx)%e > sps%emax) return
    if(this%orb(idx)%l > sps%lmax) return
    r=sps%nljz2idx(this%orb(idx)%n,this%orb(idx)%l,this%orb(idx)%j,z)
  end function iso2pn

  function pn2iso(this, sps, idx) result(r)
    class(Orbits), intent(in) :: this
    type(OrbitsIsospin), intent(in) :: sps
    integer, intent(in) :: idx
    integer :: r

    r = 0
    if(this%orb(idx)%e > sps%emax) return
    if(this%orb(idx)%l > sps%lmax) return
    r=sps%nlj2idx(this%orb(idx)%n,this%orb(idx)%l,this%orb(idx)%j)
  end function pn2iso
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


