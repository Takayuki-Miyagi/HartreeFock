! Define one-, two-, and three-body model space
! Note that the three-body space is orthonormalized
! with coefficient of fractional parentage (cfp), which will not be used
! in Hartree-Fock calculation. However, this definition may be useful for
! beyond Hartree-Fock calculation methods.

module ModelSpace
  use omp_lib
  use SingleParticleState
  implicit none

  public :: MSpace
  public :: OneBodySpace
  public :: OneBodyChannel
  public :: TwoBodySpace
  public :: TwoBodyChannel
  public :: ThreeBodySpace
  public :: ThreeBodyChannel
  public :: AdditionalQN

  private :: FinMSpace
  private :: InitMSpaceFromAZN
  private :: InitMSPaceFromReference
  private :: GetNOCoef
  private :: GetParticleHoleOrbits
  private :: GetAZNFromReference
  private :: FinOneBodySpace
  private :: InitOneBodySpace
  private :: FinOneBodyChannel
  private :: InitOneBodyChannel
  private :: FinTwoBodySpace
  private :: InitTwoBodySpace
  private :: FinTwoBodyChannel
  private :: InitTwoBodyChannel
  private :: FinThreeBodySpace
  private :: InitThreeBodySpace
  private :: FinThreeBodyChannel
  private :: InitThreeBodyChannel
  private :: FinAdditionalQN
  private :: InitAdditionalQN
  private :: CntDim
  private :: GenIter
  private :: Aop3, A3drct, A3exc1, A3exc2

  type, private :: jpz
    integer :: j = -1
    integer :: p = 0
    integer :: z = 100
  end type jpz

  type, extends(jpz) :: OneBodyChannel
    integer :: nst
    integer, allocatable :: n2spi(:)
    integer, allocatable :: spi2n(:)
  contains
    procedure :: init => InitOneBodyChannel
    procedure :: fin => FinOneBodyChannel
  end type OneBodyChannel

  type :: OneBodySpace
    type(OneBodyChannel), allocatable :: jpz(:)
    integer, allocatable :: jpz2ch(:,:,:)
    integer :: NChan
  contains
    procedure :: init => InitOneBodySpace
    procedure :: fin => FinOneBodySpace
  end type OneBodySpace

  type, extends(jpz) :: TwoBodyChannel
    integer :: nst
    integer, allocatable :: n2spi1(:)
    integer, allocatable :: n2spi2(:)
    integer, allocatable :: spis2n(:,:)
    integer, allocatable :: iphase(:,:)
  contains
    procedure :: init => InitTwoBodyChannel
    procedure :: fin => FinTwoBodyChannel
  end type TwoBodyChannel

  type :: TwoBodySpace
    type(TwoBodyChannel), allocatable :: jpz(:)
    integer, allocatable :: jpz2ch(:,:,:)
    integer :: NChan
  contains
    procedure :: init => InitTwoBodySpace
    procedure :: fin => FinTwoBodySpace
  end type TwoBodySpace

  type :: AdditionalQN
    integer :: north, nphys
    integer, allocatable :: idx2n(:)
    integer, allocatable :: n2spi1(:) ! permutations
    integer, allocatable :: n2spi2(:) ! permutations
    integer, allocatable :: n2spi3(:) ! permutations
    integer, allocatable :: n2J12(:)  ! Jab
    real(8), allocatable :: cfp(:,:)
  contains
    procedure :: init => InitAdditionalQN
    procedure :: fin => FinAdditionalQN
  end type AdditionalQN

  type, extends(jpz) :: ThreeBodyChannel
    integer :: nst, n_idx
    type(AdditionalQN), allocatable :: idx(:)
    integer, allocatable :: n2spi1(:)
    integer, allocatable :: n2spi2(:)
    integer, allocatable :: n2spi3(:)
    integer, allocatable :: n2labl(:)
    integer, allocatable :: spis2idx(:,:,:)
  contains
    procedure :: init => InitThreeBodyChannel
    procedure :: fin => FinThreeBodyChannel
  end type ThreeBodyChannel

  type :: ThreeBodySpace
    type(ThreeBodyChannel), allocatable :: jpz(:)
    integer, allocatable :: jpz2ch(:,:,:)
    integer :: NChan
  contains
    procedure :: init => InitThreeBodySpace
    procedure :: fin => FinThreeBodySpace
  end type ThreeBodySpace

  type :: MSpace
    type(Orbits) :: sps
    type(OneBodySpace) :: one
    type(TwoBodySpace) :: two
    type(ThreeBodySpace) :: thr
    logical :: is_constructed=.false.
    logical :: is_three_body =.false.
    real(8), allocatable :: NOCoef(:)
    integer, allocatable :: holes(:)
    integer, allocatable :: particles(:)
    character(:), allocatable :: Nucl
    integer :: A, Z, N
    integer :: np, nh
    integer :: emax = 0
    integer :: e2max = 0
    integer :: e3max = 0
    integer :: lmax = 0
  contains
    procedure :: fin => FinMSpace
    procedure :: InitMSpaceFromReference
    procedure :: InitMSpaceFromAZN
    generic :: init => InitMSPaceFromReference, InitMSpaceFromAZN
    procedure :: GetNOCoef
    procedure :: GetParticleHoleOrbits
  end type MSpace

  character(2), private, parameter :: elements(119) = &
      [ 'n ', &
      & 'H ', 'He', 'Li', 'Be', 'B ', 'C ', 'N ', 'O ', 'F ', 'Ne', &
      & 'Na', 'Mg', 'Al', 'Si', 'P ', 'S ', 'Cl', 'Ar', 'K ', 'Ca', &
      & 'Sc', 'Ti', 'V ', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', &
      & 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y ', 'Zr', &
      & 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', &
      & 'Sb', 'Te', 'I ', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', &
      & 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', &
      & 'Lu', 'Hf', 'Ta', 'W ', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', &
      & 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th', &
      & 'Pa', 'U ', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', &
      & 'Md', 'No', 'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds', &
      & 'Rg', 'Cn', 'Nh', 'Fl', 'Mc', 'Lv', 'Ts', 'Og' ]

  type, private :: iter3
    integer, allocatable :: ii1(:), ii2(:), ii3(:)
  contains
    procedure :: GenIter
  end type iter3
contains

  subroutine FinMSpace(this)
    class(MSpace), intent(inout) :: this
    if(.not. this%is_constructed) return
    deallocate(this%NOCoef)
    deallocate(this%holes)
    deallocate(this%particles)
    call this%sps%fin()
    call this%one%fin()
    call this%two%fin()
    if(.not. this%is_three_body) return
    call this%thr%fin()
  end subroutine FinMSpace

  subroutine InitMSpaceFromAZN(this, A, Z, N, emax, e2max, e3max, lmax)
    use Profiler, only: timer
    use ClassSys, only: sys
    class(MSpace), intent(inout) :: this
    integer, intent(in) :: A, Z, N
    integer, intent(in) :: emax, e2max
    integer, intent(in), optional :: lmax, e3max
    type(sys) :: s
    integer :: i, l
    character(20) :: Nucl
    real(8) :: ti

    ti = omp_get_wtime()
    call timer%cmemory()
    write(*,*)
    write(*,'(a)') "### Model-space constructor ###"
    Nucl = trim(s%str(A)) // trim(elements(Z+1))
    this%is_constructed = .true.
    this%A = A
    this%Z = Z
    this%N = N
    this%Nucl = adjustl(trim(Nucl))
    this%emax = emax
    this%e2max = e2max
    this%lmax = emax
    if(present(e3max)) this%e3max = e3max
    if(present(lmax)) this%lmax = lmax
    call this%sps%init(this%emax,this%lmax)
    call this%GetNOCoef()
    call this%GetParticleHoleOrbits()

    write(*,'(3a)') " Target Nuclide is ", trim(this%Nucl), ", Orbits:"
    write(*,'(a)') "      p/h, idx,  n,  l,  j, tz,   occupation"
    do i = 1, this%nh
      l = this%holes(i)
      write(*,'(a10,5i4,f14.6)') '     hole:', l, this%sps%orb(l)%n, this%sps%orb(l)%l, &
          & this%sps%orb(l)%j, this%sps%orb(l)%z, this%NOcoef(l)
    end do

    do i = 1, this%np
      l = this%particles(i)
      write(*,'(a10,5i4,f14.6)') ' particle:', l, this%sps%orb(l)%n, this%sps%orb(l)%l, &
          & this%sps%orb(l)%j, this%sps%orb(l)%z, this%NOcoef(l)
    end do
    write(*,*)

    call this%one%init(this%sps)
    call this%two%init(this%sps, this%e2max)
    write(*,'(a)') "  # of J, parity, and tz Channels:"
    write(*,'(a,i3)') "    OneBody: ", this%one%NChan
    write(*,'(a,i3)') "    TwoBody: ", this%two%NChan
    if(.not. present(e3max)) then
      write(*,*)
      call timer%Add('Construct Model Space', omp_get_wtime()-ti)
      return
    end if
    this%is_three_body = .true.
    call this%thr%init(this%sps, this%e2max, this%e3max)
    write(*,'(a,i3)') "  ThreeBody: ", this%thr%NChan
    write(*,*)
    call timer%tmemory("Model Space")
    call timer%Add('Construct Model Space', omp_get_wtime()-ti)
  end subroutine InitMSpaceFromAZN

  subroutine InitMSpaceFromReference(this, Nucl, emax, e2max, e3max, lmax)
    class(MSpace), intent(inout) :: this
    character(*), intent(in) :: Nucl
    integer, intent(in) :: emax, e2max
    integer, intent(in), optional :: e3max, lmax
    integer :: A, Z, N

    call GetAZNFromReference(Nucl,A,Z,N)
    call this%init(A,Z,N,emax,e2max,e3max,lmax)
  end subroutine InitMSpaceFromReference

  subroutine GetNOCoef(this)
    ! This should be called after obtaining A, Z, N
    class(MSpace), intent(inout) :: this
    integer :: Z, N
    integer :: e, l, j, g, ns
    integer :: zz, nn, vz, vn
    integer :: idxp, idxn

    allocate(this%NOCoef(this%sps%norbs))
    this%NOCoef(:) = 0.d0
    Z = this%Z
    N = this%N
    l = 0
    nn = 0
    zz = 0
    do e = 0, this%sps%emax
      do g = 2*e+1, -2*e+1, -4
        j = abs(g)
        if(g > 0) l = (j-1)/2
        if(g < 0) l = (j+1)/2
        ns = (e-l)/2
        if(l > this%sps%lmax) cycle

        nn = nn + (j+1)
        zz = zz + (j+1)
        idxp = this%sps%nljz2idx(ns,l,j,-1)
        idxn = this%sps%nljz2idx(ns,l,j, 1)
        if(Z - zz >= 0) this%NOcoef(idxp) = 1.d0
        if(Z - zz < 0) then
          vz = Z - zz + j + 1
          if(vz > 0) this%NOcoef(idxp) = dble(vz) / dble(j+1)
          if(vz < 0) this%NOcoef(idxp) = 0.d0
        end if

        if(N - nn >= 0) this%NOcoef(idxn) = 1.d0
        if(N - nn < 0) then
          vn = N - nn + j + 1
          if(vn > 0) this%NOcoef(idxn) = dble(vn) / dble(j+1)
          if(vn < 0) this%NOcoef(idxn) = 0.d0
        end if

      end do
    end do

    if(Z - zz > 0 .or. N - nn > 0) then
      write(*,'(a)', advance='no') "Error: emax is too small: "
      write(*,'(a,i4,a,i3,a,i3)') "A = ", this%A, ", Z = ", this%Z, ", N = ", this%N
      stop
    end if

#ifdef ModelSpaceDebug
    write(*,'(a)') "In GetNOCoef:"
    write(*,'(a)') "   n,  l,  j, tz,   occupation"
    do l = 1, this%sps%norbs
      write(*,'(4i4,f14.6)') this%sps%orb(l)%n, this%sps%orb(l)%l, &
          & this%sps%orb(l)%j, this%sps%orb(l)%z, this%NOcoef(l)
    end do
#endif
  end subroutine GetNOCoef

  subroutine GetParticleHoleOrbits(this)
    ! This should be called after obtaining NOCoef array
    class(MSpace), intent(inout) :: this
    integer :: n_h, n_p, i

    n_h = 0
    n_p = 0
    do i = 1, this%sps%norbs
      if(abs(1.d0 - this%NOCoef(i)) < 1.d-6) n_h = n_h + 1
      if(abs(1.d0 - this%NOCoef(i)) > 1.d-6) n_p = n_p + 1
    end do
    this%nh = n_h
    this%np = n_p

    allocate(this%holes(n_h))
    allocate(this%particles(n_p))

    n_h = 0
    n_p = 0
    do i = 1, this%sps%norbs
      if(abs(1.d0 - this%NOCoef(i)) < 1.d-6) then
        n_h = n_h + 1
        this%holes(n_h) = i
      end if

      if(abs(1.d0 - this%NOCoef(i)) > 1.d-6) then
        n_p = n_p + 1
        this%particles(n_p) = i
      end if
    end do

  end subroutine GetParticleHoleOrbits

  subroutine GetAZNFromReference(Nucl,A,Z,N)
    ! Nucl has to be "str(mass number)+element" or "element+str(mass number)"
    ! It must not include space
    character(*), intent(in) :: Nucl
    integer, intent(out) :: A, Z, N
    integer :: i
    character(:), allocatable :: str
    character(20) :: str_A, element

    str = Nucl
    str_A = ''
    element= ''
    ! "str(mass number)+element" case
    do
      if(scan(str,'1234567890') /= 1) exit
      str_A = trim(str_A) // str(1:1)
      str = str(2:)
    end do
    element = str

    if(str_A /= '') then
      do i = 1, size(elements)
        if(trim(element) == trim(elements(i))) then
          Z = i - 1
          exit
        end if

        if(i == size(elements)) then
          write(*,'(2a)') 'Error in GetAZNFromRefernece: Nucl is ', trim(Nucl)
          stop
        end if
      end do
      read(str_A,*) A
      N = A - Z
      return
    end if

    str = Nucl
    str_A = ''
    element= ''
    ! "element+str(mass number)" case
    do
      if(scan(str,'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ') /= 1) exit
      element = trim(element) // str(1:1)
      str = str(2:)
    end do
    str_A = str
    if(element /= '') then
      do i = 1, size(elements)
        if(trim(element) == trim(elements(i))) then
          Z = i - 1
          exit
        end if

        if(i == size(elements)) then
          write(*,'(2a)') 'Error in GetAZNFromRefernece: Nucl is ', trim(Nucl)
          stop
        end if
      end do
      read(str_A,*) A
      N = A - Z
      return
    end if

    write(*,'(2a)') 'Error in GetAZNFromRefernece: Nucl is ', trim(Nucl)
    stop
  end subroutine GetAZNFromReference

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
    type(Orbits), intent(in) :: sps
    integer :: ich
    integer :: j, p, z, n, i
    integer, allocatable :: jj(:), pp(:), zz(:), nn(:)

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

  subroutine FinOneBodyChannel(this)
    class(OneBodyChannel), intent(inout) :: this
    deallocate(this%n2spi)
    deallocate(this%spi2n)
  end subroutine FinOneBodyChannel

  subroutine InitOneBodyChannel(this, j, p, z, n, sps)
    class(OneBodyChannel), intent(inout) :: this
    type(Orbits), intent(in) :: sps
    integer, intent(in) :: j, p, z, n
    integer :: cnt, i
    this%j = j
    this%p = p
    this%z = z
    this%nst = n
#ifdef ModelSpaceDebug
    write(*,'(a,i3,a,i3,a,i3,a)') "One-body channel: J=", j, "/2, P=", p, ", Tz=", z, "/2"
#endif
    allocate(this%n2spi(n))
    allocate(this%spi2n(sps%norbs))
    this%spi2n(:) = 0
    cnt = 0
    do i = 1, sps%norbs
      if(sps%orb(i)%j /= j) cycle
      if((-1)**sps%orb(i)%l /= p) cycle
      if(sps%orb(i)%z /= z) cycle
      cnt = cnt + 1
      this%n2spi(cnt) = i
      this%spi2n(i) = cnt
#ifdef ModelSpaceDebug
      write(*,'(a,i3,a,i3)') "i1=",i,", Num=",cnt
#endif
    end do
#ifdef ModelSpaceDebug
    write(*,*)
#endif
  end subroutine InitOneBodyChannel

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
    use CommonLibrary, only: triag
    class(TwoBodySpace), intent(inout) :: this
    type(Orbits), intent(in) :: sps
    integer, intent(in) :: e2max
    integer :: j, p, z, n, ich
    integer :: i1, j1, l1, z1, e1
    integer :: i2, j2, l2, z2, e2
    integer, allocatable :: jj(:), pp(:), zz(:), nn(:)

    allocate(this%jpz2ch(0:min(2*sps%lmax,e2max)+1,-1:1,-1:1))
    this%jpz2ch(:,:,:) = 0
    ich = 0
    do j = 0, 2*sps%lmax+1
      do p = 1, -1, -2
        do z = -1, 1
          n = 0

          do i1 = 1, sps%norbs
            j1 = sps%orb(i1)%j
            l1 = sps%orb(i1)%l
            z1 = sps%orb(i1)%z
            e1 = sps%orb(i1)%e

            do i2 = 1, i1
              j2 = sps%orb(i2)%j
              l2 = sps%orb(i2)%l
              z2 = sps%orb(i2)%z
              e2 = sps%orb(i1)%e

              if(e1 + e2 > e2max) cycle
              if(triag(j1, j2, 2*j)) cycle
              if((-1) ** (l1+l2) /= p) cycle
              if(z1 + z2 /= 2*z) cycle
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
    allocate(this%jpz(this%NChan))
    allocate(jj(this%NChan))
    allocate(pp(this%NChan))
    allocate(zz(this%NChan))
    allocate(nn(this%NChan))

    ich = 0
    do j = 0, 2*sps%lmax+1
      do p = 1, -1, -2
        do z = -1, 1
          n = 0

          do i1 = 1, sps%norbs
            j1 = sps%orb(i1)%j
            l1 = sps%orb(i1)%l
            z1 = sps%orb(i1)%z
            e1 = sps%orb(i1)%e

            do i2 = 1, i1
              j2 = sps%orb(i2)%j
              l2 = sps%orb(i2)%l
              z2 = sps%orb(i2)%z
              e2 = sps%orb(i1)%e

              if(e1 + e2 > e2max) cycle
              if(triag(j1, j2, 2*j)) cycle
              if((-1) ** (l1+l2) /= p) cycle
              if(z1 + z2 /= 2*z) cycle
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
  end subroutine InitTwoBodySpace

  subroutine FinTwoBodyChannel(this)
    class(TwoBodyChannel), intent(inout) :: this
    deallocate(this%n2spi1)
    deallocate(this%n2spi2)
    deallocate(this%spis2n)
    deallocate(this%iphase)
  end subroutine FinTwoBodyChannel

  subroutine InitTwoBodyChannel(this, j, p, z, n, sps, e2max)
    use CommonLibrary, only: triag
    class(TwoBodyChannel), intent(inout) :: this
    integer, intent(in) :: j, p, z, n, e2max
    type(Orbits), intent(in) :: sps
    integer :: i1, j1, l1, z1, e1
    integer :: i2, j2, l2, z2, e2
    integer :: cnt

    this%j = j
    this%p = p
    this%z = z
    this%nst = n
    allocate(this%n2spi1(n))
    allocate(this%n2spi2(n))
    allocate(this%spis2n(sps%norbs,sps%norbs))
    allocate(this%iphase(sps%norbs,sps%norbs))
    cnt = 0
#ifdef ModelSpaceDebug
    write(*,'(a,i3,a,i3,a,i3)') "Two-body channel: J=", j, ", P=", p, ", Tz=", z
#endif
    do i1 = 1, sps%norbs
      j1 = sps%orb(i1)%j
      l1 = sps%orb(i1)%l
      z1 = sps%orb(i1)%z
      e1 = sps%orb(i1)%e
      do i2 = 1, i1
        j2 = sps%orb(i2)%j
        l2 = sps%orb(i2)%l
        z2 = sps%orb(i2)%z
        e2 = sps%orb(i2)%e

        if(e1 + e2 > e2max) cycle
        if(triag(j1, j2, 2*j)) cycle
        if((-1) ** (l1+l2) /= p) cycle
        if(z1 + z2 /= 2*z) cycle
        if(i1 == i2 .and. mod(j,2) == 1) cycle

        cnt = cnt + 1
        this%n2spi1(cnt) = i1
        this%n2spi2(cnt) = i2
        this%spis2n(i1,i2) = cnt
        this%spis2n(i2,i1) = cnt
        this%iphase(i1,i2) = 1
        this%iphase(i1,i2) = -(-1) ** ((j1+j2)/2 - j)
#ifdef ModelSpaceDebug
        write(*,'(a,i3,a,i3,a,i6)') "i1=",i1,", i2=",i2,", Num=",cnt
#endif
      end do
    end do
#ifdef ModelSpaceDebug
    write(*,*)
#endif
  end subroutine InitTwoBodyChannel

  subroutine FinThreeBodySpace(this)
    class(ThreeBodySpace), intent(inout) :: this
    integer :: ich
    do ich = 1, this%NChan
      call this%jpz(ich)%fin()
    end do
    deallocate(this%jpz)
    deallocate(this%jpz2ch)
  end subroutine FinThreeBodySpace

  subroutine InitThreeBodySpace(this, sps, e2max, e3max)
    class(ThreeBodySpace), intent(inout) :: this
    type(Orbits), intent(in) :: sps
    integer, intent(in) :: e2max, e3max
    integer :: j, p, z, ich, nidx, n
    integer :: i1, j1, l1, z1, e1
    integer :: i2, j2, l2, z2, e2
    integer :: i3, j3, l3, z3, e3
    integer :: ni, nj
    integer, allocatable :: jj(:), pp(:), zz(:), nn(:), nnidx(:)
    type :: spis2n
      integer, allocatable :: spis2nidx(:,:,:)
    end type spis2n
    type(spis2n), allocatable :: ch(:)

    allocate(this%jpz2ch(1:2*min(e3max,3*sps%lmax)+3,-1:1,-3:3))
    this%jpz2ch(:,:,:) = 0
    ich = 0
    do j = 1, 2*min(e3max,3*sps%lmax)+3, 2
      do p = 1, -1, -2
        do z = -3, 3, 2
          n = 0
          nidx = 0
          do i1 = 1, sps%norbs
            j1 = sps%orb(i1)%j
            l1 = sps%orb(i1)%l
            z1 = sps%orb(i1)%z
            e1 = sps%orb(i1)%e
            do i2 = 1, i1
              j2 = sps%orb(i2)%j
              l2 = sps%orb(i2)%l
              z2 = sps%orb(i2)%z
              e2 = sps%orb(i2)%e
              if(e1 + e2 > e2max) cycle
              do i3 = 1, i2
                j3 = sps%orb(i3)%j
                l3 = sps%orb(i3)%l
                z3 = sps%orb(i3)%z
                e3 = sps%orb(i3)%e
                if(e1 + e3 > e2max) cycle
                if(e2 + e3 > e2max) cycle
                if(e1 + e2 + e3 > e3max) cycle
                if(z1 + z2 + z3 /= z) cycle
                if((-1) ** (l1+l2+l3) /= p) cycle

                call CntDim(sps,i1,i2,i3,j,ni,nj)
                n = n + ni
                if(ni /= 0) nidx = nidx + 1
              end do
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
    allocate(this%jpz(ich))
    allocate(jj(ich))
    allocate(pp(ich))
    allocate(zz(ich))
    allocate(nn(ich))
    allocate(nnidx(ich))
    allocate(ch(ich))
    do ich = 1, this%NChan
      allocate(ch(ich)%spis2nidx(sps%norbs,sps%norbs,sps%norbs))
      ch(ich)%spis2nidx(:,:,:) = 0
    end do

    ich = 0
    do j = 1, 2*min(e3max,3*sps%lmax)+3, 2
      do p = 1, -1, -2
        do z = -3, 3, 2
          n = 0
          nidx = 0
          do i1 = 1, sps%norbs
            j1 = sps%orb(i1)%j
            l1 = sps%orb(i1)%l
            z1 = sps%orb(i1)%z
            e1 = sps%orb(i1)%e
            do i2 = 1, i1
              j2 = sps%orb(i2)%j
              l2 = sps%orb(i2)%l
              z2 = sps%orb(i2)%z
              e2 = sps%orb(i2)%e
              if(e1 + e2 > e2max) cycle
              do i3 = 1, i2
                j3 = sps%orb(i3)%j
                l3 = sps%orb(i3)%l
                z3 = sps%orb(i3)%z
                e3 = sps%orb(i3)%e
                if(e1 + e3 > e2max) cycle
                if(e2 + e3 > e2max) cycle
                if(e1 + e2 + e3 > e3max) cycle
                if(z1 + z2 + z3 /= z) cycle
                if((-1) ** (l1+l2+l3) /= p) cycle

                call CntDim(sps,i1,i2,i3,j,ni,nj)
                n = n + ni
                if(ni /= 0) then
                  nidx = nidx + 1
                  ch(ich+1)%spis2nidx(i1,i2,i3) = nidx
                end if
              end do
            end do
          end do
          if(n /= 0) then
            ich = ich + 1
            jj(ich) = j
            pp(ich) = p
            zz(ich) = z
            nn(ich) = n
            nnidx(ich) = nidx
          end if

        end do
      end do
    end do

    do ich = 1, this%NChan
      call this%jpz(ich)%init(jj(ich),pp(ich),zz(ich),nn(ich),nnidx(ich),&
          & ch(ich)%spis2nidx, sps, e2max, e3max)
      deallocate(ch(ich)%spis2nidx)
    end do
    deallocate(ch)
  end subroutine InitThreeBodySpace

  subroutine FinThreeBodyChannel(this)
    class(ThreeBodyChannel), intent(inout) :: this
    integer :: idx_sub

    do idx_sub = 1, this%n_idx
      call this%idx(idx_sub)%fin()
    end do
    deallocate(this%idx)
    deallocate(this%n2spi1)
    deallocate(this%n2spi2)
    deallocate(this%n2spi3)
    deallocate(this%n2labl)
    deallocate(this%spis2idx)
  end subroutine FinThreeBodyChannel

  subroutine InitThreeBodyChannel(this, j, p, z, n, nidx, spis2idx, sps, e2max, e3max)
    class(ThreeBodyChannel), intent(inout) :: this
    type(Orbits), intent(in) :: sps
    integer, intent(in) :: j, p, z, n, nidx, spis2idx(:,:,:), e2max, e3max
    integer :: i1, j1, l1, z1, e1
    integer :: i2, j2, l2, z2, e2
    integer :: i3, j3, l3, z3, e3
    integer :: ni, nj, idx, i, cnt
    integer, allocatable :: nni(:), nnj(:)

    this%j = j
    this%p = p
    this%z = z
    this%nst = n
    this%n_idx = nidx
    this%spis2idx = spis2idx
#ifdef ModelSpaceDebug
    write(*,'(a,i3,a,i3,a,i3,a)') "Three-body channel: J=", j, "/2, P=", p, ", Tz=", z, "/2"
#endif

    allocate(this%idx(this%n_idx))
    allocate(nni(this%n_idx))
    allocate(nnj(this%n_idx))
    allocate(this%n2spi1(n))
    allocate(this%n2spi2(n))
    allocate(this%n2spi3(n))
    allocate(this%n2labl(n))

    cnt = 0
    do i1 = 1, sps%norbs
      j1 = sps%orb(i1)%j
      l1 = sps%orb(i1)%l
      z1 = sps%orb(i1)%z
      e1 = sps%orb(i1)%e
      do i2 = 1, i1
        j2 = sps%orb(i2)%j
        l2 = sps%orb(i2)%l
        z2 = sps%orb(i2)%z
        e2 = sps%orb(i2)%e
        if(e1 + e2 > e2max) cycle
        do i3 = 1, i2
          j3 = sps%orb(i3)%j
          l3 = sps%orb(i3)%l
          z3 = sps%orb(i3)%z
          e3 = sps%orb(i3)%e
          if(e1 + e3 > e2max) cycle
          if(e2 + e3 > e2max) cycle
          if(e1 + e2 + e3 > e3max) cycle
          if(z1 + z2 + z3 /= z) cycle
          if((-1) ** (l1+l2+l3) /= p) cycle
          idx = this%spis2idx(i1,i2,i3)
          if(idx == 0) cycle
          call CntDim(sps,i1,i2,i3,j,ni,nj)
          nni(idx) = ni
          nnj(idx) = nj
          allocate(this%idx(idx)%idx2n(ni))
          do i = 1, ni
            cnt = cnt + 1
            this%n2spi1(cnt) = i1
            this%n2spi2(cnt) = i2
            this%n2spi3(cnt) = i3
            this%n2labl(cnt) = i
          end do
        end do
      end do
    end do

    cnt = 0
    do i1 = 1, sps%norbs
      j1 = sps%orb(i1)%j
      l1 = sps%orb(i1)%l
      z1 = sps%orb(i1)%z
      e1 = sps%orb(i1)%e

      do i2 = 1, i1
        j2 = sps%orb(i2)%j
        l2 = sps%orb(i2)%l
        z2 = sps%orb(i2)%z
        e2 = sps%orb(i2)%e

        if(e1 + e2 > e2max) cycle
        do i3 = 1, i2
          j3 = sps%orb(i3)%j
          l3 = sps%orb(i3)%l
          z3 = sps%orb(i3)%z
          e3 = sps%orb(i3)%e

          if(e1 + e3 > e2max) cycle
          if(e2 + e3 > e2max) cycle
          if(e1 + e2 + e3 > e3max) cycle
          if(z1 + z2 + z3 /= z) cycle
          if((-1) ** (l1+l2+l3) /= p) cycle
          idx = this%spis2idx(i1,i2,i3)
          if(idx == 0) cycle
          ni = nni(idx)
          nj = nnj(idx)
          call this%idx(idx)%init(sps,i1,i2,i3,j,ni,nj)
          do i = 1, ni
            cnt = cnt + 1
            this%n2spi1(cnt) = i1
            this%n2spi2(cnt) = i2
            this%n2spi3(cnt) = i3
            this%n2labl(cnt) = i
            this%idx(idx)%idx2n(i) = cnt
#ifdef ModelSpaceDebug
            write(*,'(a,i3,a,i3,a,i3,a,i3,a,i8)') &
                & "i1=",i1,", i2=",i2,", i3=",i3, ", label=", i, ", Num=",cnt
#endif
          end do

        end do
      end do
    end do

#ifdef ModelSpaceDebug
    write(*,*)
#endif

  end subroutine InitThreeBodyChannel

  subroutine FinAdditionalQN(this)
    class(AdditionalQN), intent(inout) :: this
    deallocate(this%idx2n)
    deallocate(this%n2spi1)
    deallocate(this%n2spi2)
    deallocate(this%n2spi3)
    deallocate(this%n2J12)
    deallocate(this%cfp)
  end subroutine FinAdditionalQN

  subroutine InitAdditionalQN(this, sps, i1, i2, i3, j, nni, nnj)
    use LinAlgLib
    use CommonLibrary, only: triag
    class(AdditionalQN), intent(inout) :: this
    type(Orbits), intent(in) :: sps
    integer, intent(in) :: i1, i2, i3, j, nni, nnj
    type(iter3) :: ite
    integer :: nj, nch, i, k, kk
    integer :: bra, ket
    integer :: ii1, ii2, ii3, jj1, jj2, jj3, jj12
    integer :: ii4, ii5, ii6, jj45
    type(DMat) :: mat
    type(EigenSolSymD) :: sol
    call ite%GenIter(i1, i2, i3, nch)
    this%nphys = nnj
    this%north = nni
    allocate(this%n2spi1(nnj))
    allocate(this%n2spi2(nnj))
    allocate(this%n2spi3(nnj))
    allocate(this%n2J12( nnj))
    allocate(this%cfp(nnj,nni))
    nj = 0
    do i = 1, nch
      ii1 = ite%ii1(i)
      ii2 = ite%ii2(i)
      ii3 = ite%ii3(i)
      jj1 = sps%orb(ii1)%j
      jj2 = sps%orb(ii2)%j
      jj3 = sps%orb(ii3)%j
      do jj12 = iabs(jj1 - jj2) / 2, (jj1 + jj2) / 2
        if(ii1 == ii2 .and. mod(jj12, 2) == 1) cycle
        if(triag(2*jj12, jj3, j)) cycle
        nj = nj + 1
        this%n2spi1(nj) = ii1
        this%n2spi2(nj) = ii2
        this%n2spi3(nj) = ii3
        this%n2J12( nj) = jj12
      end do
    end do
    if(this%nphys < 1) return
    call mat%Ini(this%nphys, this%nphys)
    do bra = 1, this%nphys
      ii1  = this%n2spi1(bra)
      ii2  = this%n2spi2(bra)
      ii3  = this%n2spi3(bra)
      jj12 = this%n2J12( bra)
      do ket = 1, this%nphys
        ii4  = this%n2spi1(ket)
        ii5  = this%n2spi2(ket)
        ii6  = this%n2spi3(ket)
        jj45 = this%n2J12( ket)
        mat%m(bra, ket) = Aop3(sps, ii1, ii2, ii3, jj12, &
            & ii4, ii5, ii6, jj45, j)
      end do
    end do
    call sol%init(mat)
    call sol%DiagSym(mat)

    do i = 1, this%nphys
      kk = 0
      do k = 1, this%nphys
        if(abs(1.d0-sol%eig%v(k)).le.1.d-4) then
          kk = kk + 1
          this%cfp(i,kk) = sol%vec%m(i, k)
        end if

        if(abs(1.d0 - sol%eig%v(k)) > 1.d-4 .and. abs(sol%eig%v(k)) > 1.d-4) then
          write(*,'(10f12.6)') sol%eig%v(k)
        end if
      end do
    end do
    call mat%Fin()
    call sol%fin()
  end subroutine InitAdditionalQN

  subroutine CntDim(sps, i1, i2, i3, j, ni, nj)
    use CommonLibrary, only: triag
    type(Orbits), intent(in) :: sps
    integer, intent(in) :: i1, i2, i3, j
    integer, intent(out) :: ni, nj
    type(iter3) :: ite
    integer :: ii1, ii2, ii3, jj1, jj2, jj3, jj12
    integer :: i, nch, jmin, jmax
    real(8) :: tr
    nj = 0; ni = 0
    call ite%GenIter(i1, i2, i3, nch)
    tr = 0.d0
    jmin = 100; jmax = -100
    do i = 1, nch
      ii1 = ite%ii1(i)
      ii2 = ite%ii2(i)
      ii3 = ite%ii3(i)
      jj1 = sps%orb(ii1)%j
      jj2 = sps%orb(ii2)%j
      jj3 = sps%orb(ii3)%j
      do jj12 = iabs(jj1 - jj2) / 2, (jj1 + jj2) / 2
        if(ii1 == ii2 .and. mod(jj12, 2) == 1) cycle
        if(triag(2*jj12, jj3, j)) cycle
        if(jmin > jj12) jmin = jj12
        if(jmax < jj12) jmax = jj12
        nj = nj + 1
        tr = tr + Aop3(sps, ii1, ii2, ii3, jj12, ii1, ii2, ii3, jj12, j)
      end do
    end do
    ni = int(tr + 1.d-2)
    deallocate(ite%ii1, ite%ii2, ite%ii3)
  end subroutine CntDim

  subroutine GenIter(ite, i1, i2, i3, num)
    class(iter3), intent(inout) :: ite
    integer, intent(in) :: i1, i2, i3
    integer, intent(out) :: num
    integer :: icase
    if(allocated(ite%ii1)) deallocate(ite%ii1)
    if(allocated(ite%ii2)) deallocate(ite%ii2)
    if(allocated(ite%ii3)) deallocate(ite%ii3)

    icase = -100
    if(i2 == i3 .and. i1 == i2) icase = 1
    if(i2 == i3 .and. i1 /= i2) icase = 2
    if(i2 /= i3 .and. i1 == i2) icase = 3
    if(i2 /= i3 .and. i1 == i3) icase = 4
    if(i2 /= i3 .and. i1 /= i2 .and. i1 /= i3) icase = 5

    select case(icase)
    case(1)
      num = 1
      allocate(ite%ii1(num), ite%ii2(num), ite%ii3(num))
      ite%ii1(1) = i1; ite%ii2(1) = i1; ite%ii3(1) = i1
    case(2)
      num = 3
      allocate(ite%ii1(num), ite%ii2(num), ite%ii3(num))
      ite%ii1(1) = i1; ite%ii2(1) = i2; ite%ii3(1) = i2
      ite%ii1(2) = i2; ite%ii2(2) = i1; ite%ii3(2) = i2
      ite%ii1(3) = i2; ite%ii2(3) = i2; ite%ii3(3) = i1
    case(3)
      num = 3
      allocate(ite%ii1(num), ite%ii2(num), ite%ii3(num))
      ite%ii1(1) = i1; ite%ii2(1) = i1; ite%ii3(1) = i3
      ite%ii1(2) = i3; ite%ii2(2) = i1; ite%ii3(2) = i1
      ite%ii1(3) = i1; ite%ii2(3) = i3; ite%ii3(3) = i1
    case(4)
      num = 3
      allocate(ite%ii1(num), ite%ii2(num), ite%ii3(num))
      ite%ii1(1) = i1; ite%ii2(1) = i2; ite%ii3(1) = i1
      ite%ii1(2) = i1; ite%ii2(2) = i1; ite%ii3(2) = i2
      ite%ii1(3) = i2; ite%ii2(3) = i1; ite%ii3(3) = i1
    case(5)
      num = 6
      allocate(ite%ii1(num), ite%ii2(num), ite%ii3(num))
      ite%ii1(1) = i1; ite%ii2(1) = i2; ite%ii3(1) = i3
      ite%ii1(2) = i3; ite%ii2(2) = i1; ite%ii3(2) = i2
      ite%ii1(3) = i2; ite%ii2(3) = i3; ite%ii3(3) = i1
      ite%ii1(4) = i2; ite%ii2(4) = i1; ite%ii3(4) = i3
      ite%ii1(5) = i3; ite%ii2(5) = i2; ite%ii3(5) = i1
      ite%ii1(6) = i1; ite%ii2(6) = i3; ite%ii3(6) = i2
    case default
      write(*,*) 'Error in GenIter:'
      return
    end select
  end subroutine GenIter

  real(8) function Aop3(sps, i1, i2, i3, j12, &
        & i4, i5, i6, j45, j) result(a)
    type(Orbits), intent(in) :: sps
    integer, intent(in) :: i1, i2, i3, i4, i5, i6
    integer, intent(in) :: j12, j45, j
    real(8) :: phase12, phase23, phase31
    phase12 = (-1.d0) ** ((sps%orb(i1)%j + sps%orb(i2)%j)/2 - j45)
    phase23 = (-1.d0) ** ((sps%orb(i2)%j + sps%orb(i3)%j)/2 - j45)
    phase31 = (-1.d0) ** ((sps%orb(i3)%j + sps%orb(i1)%j)/2 - j45)
    a = 0.d0
    a = a + A3drct(     i1, i2, i3, j12, i4, i5, i6, j45   )
    a = a - A3drct(     i1, i2, i3, j12, i5, i4, i6, j45   ) * phase12
    a = a + A3exc1(sps, i1, i2, i3, j12, i4, i5, i6, j45, j)
    a = a - A3exc1(sps, i1, i2, i3, j12, i5, i4, i6, j45, j) * phase31
    a = a + A3exc2(sps, i1, i2, i3, j12, i4, i5, i6, j45, j)
    a = a - A3exc2(sps, i1, i2, i3, j12, i5, i4, i6, j45, j) * phase23
    a = a / 6.d0
  end function Aop3

  real(8) function A3drct(i1, i2, i3, j12, i4, i5, i6, j45) result(d)
    integer, intent(in) :: i1, i2, i3, i4, i5, i6
    integer, intent(in) :: j12, j45
    d = 0.d0
    if(i1 /= i4) return
    if(i2 /= i5) return
    if(i3 /= i6) return
    if(j12 /= j45) return
    d = 1.d0
  end function A3drct

  real(8) function A3exc1(sps, i1, i2, i3, j12, i4, i5, i6, j45, j) result(e)
    use CommonLibrary, only: sjs
    type(Orbits), intent(in) :: sps
    integer, intent(in) :: i1, i2, i3, i4, i5, i6
    integer, intent(in) :: j12, j45, j
    integer :: j1, j2, j3
    e = 0.d0
    if(i1 /= i5) return
    if(i2 /= i6) return
    if(i3 /= i4) return
    j1 = sps%orb(i1)%j
    j2 = sps%orb(i2)%j
    j3 = sps%orb(i3)%j
    e = - (-1.d0) ** ((j1 + j2) / 2 + j12) * &
        & dsqrt(dble(2 * j12 + 1)) * dsqrt(dble(2 * j45 + 1)) * &
        & sjs(j1, j2, 2*j12, j, j3, 2*j45)
  end function A3exc1

  real(8) function A3exc2(sps, i1, i2, i3, j12, i4, i5, i6, j45, j) result(e)
    use CommonLibrary, only: sjs
    type(Orbits), intent(in) :: sps
    integer, intent(in) :: i1, i2, i3, i4, i5, i6
    integer, intent(in) :: j12, j45, j
    integer :: j1, j2, j3
    e = 0.d0
    if(i1 /= i6) return
    if(i2 /= i4) return
    if(i3 /= i5) return
    j1 = sps%orb(i1)%j
    j2 = sps%orb(i2)%j
    j3 = sps%orb(i3)%j
    e = - (-1.d0) ** ((j2 + j3) / 2 + j45) * &
        & dsqrt(dble(2 * j12 + 1)) * dsqrt(dble(2 * j45 + 1)) * &
        & sjs(j1, j2, 2*j12, j3, j, 2*j45)
  end function A3exc2
end module ModelSpace

! main for test
program main
  use Profiler, only: timer
  use ModelSpace, only: MSpace
  use CommonLibrary, only: &
      &init_dbinomial_triangle, fin_dbinomial_triangle
  type(MSpace) :: ms

  call timer%init()
  call init_dbinomial_triangle()
  !call ms%init('O18', 4, 8, e3max=4, lmax=4)
  call ms%init('16O', 10, 20, e3max=10)
  call ms%fin()
  call fin_dbinomial_triangle()
  call timer%fin()
end program main
