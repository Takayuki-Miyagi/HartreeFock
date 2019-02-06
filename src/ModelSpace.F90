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

  private :: FinNonOrthIsospinThreeBodySpace
  private :: InitNonOrthIsospinThreeBodySpace
  private :: FinNonOrthIsospinThreeBodyChannel
  private :: InitNonOrthIsospinThreeBodyChannel
  private :: SetSortingIndices
  private :: Sort123
  private :: fin_sort_index
  private :: init_sort_index

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
    integer :: emax
  contains
    procedure :: init => InitOneBodySpace
    procedure :: fin => FinOneBodySpace
    procedure :: GetOneBodyChannel
    procedure :: GetOneBodyChannelJPZ
    generic :: get => GetOneBodyChannel, GetOneBodyChannelJPZ
  end type OneBodySpace

  type, extends(jpz) :: TwoBodyChannel
    integer :: nst
    integer, allocatable :: n2spi1(:)
    integer, allocatable :: n2spi2(:)
    integer, allocatable :: spis2n(:,:)
    integer, allocatable :: iphase(:,:)
    integer, allocatable :: holes(:)
    integer, allocatable :: particles(:)
  contains
    procedure :: init => InitTwoBodyChannel
    procedure :: fin => FinTwoBodyChannel
  end type TwoBodyChannel

  type :: TwoBodySpace
    type(TwoBodyChannel), allocatable :: jpz(:)
    integer, allocatable :: jpz2ch(:,:,:)
    integer :: emax, e2max
    integer :: NChan
  contains
    procedure :: init => InitTwoBodySpace
    procedure :: fin => FinTwoBodySpace
    procedure :: GetTwoBodyChannel
    procedure :: GetTwoBodyChannelJPZ
    generic :: get => GetTwoBodyChannel, GetTwoBodyChannelJPZ
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
    integer :: emax
    integer :: e2max
    integer :: e3max
  contains
    procedure :: init => InitThreeBodySpace
    procedure :: fin => FinThreeBodySpace
    procedure :: GetThreeBodyChannel
    procedure :: GetThreeBodyChannelJPZ
    generic :: get => GetThreeBodyChannel, GetThreeBodyChannelJPZ
  end type ThreeBodySpace

  type :: coef
    integer :: n = 0
    integer, allocatable :: idx2num(:)
#ifdef single_precision
    real(4), allocatable :: TrnsCoef(:)
#else
    real(8), allocatable :: TrnsCoef(:)
#endif
  end type coef

  type :: sort_index
    integer :: idx_sorted = 0
    integer :: j12min, j12max
    integer :: t12min, t12max
    type(coef), allocatable :: jt(:,:)
  contains
    procedure :: fin => fin_sort_index
    procedure :: init => init_sort_index
  end type sort_index

  type :: NonOrthIsospinAdditionalQN
    integer :: j12min, j12max
    integer :: t12min, t12max
    integer, allocatable :: JT2n(:,:)
  end type NonOrthIsospinAdditionalQN

  type :: NonOrthIsospinThreeBodyChannel
    integer :: nst, n_idx, n_sort
    integer :: j, p, t
    type(sort_index), allocatable :: sort(:)
    type(NonOrthIsospinAdditionalQN), allocatable :: idxqn(:)
    integer, allocatable :: n2spi1(:)
    integer, allocatable :: n2spi2(:)
    integer, allocatable :: n2spi3(:)
    integer, allocatable :: n2J12( :)
    integer, allocatable :: n2T12( :)
    integer, allocatable :: spis2idx(:,:,:) ! => NonOrthIsospinAdditionalQN
    integer, allocatable :: sorting(:,:,:)
  contains
    procedure :: init => InitNonOrthIsospinThreeBodyChannel
    procedure :: fin => FinNonOrthIsospinThreeBodyChannel
    procedure :: set => SetSortingIndices
  end type NonOrthIsospinThreeBodyChannel

  type :: NonOrthIsospinThreeBodySpace
    type(NonOrthIsospinThreeBodyChannel), allocatable :: jpt(:)
    integer, allocatable :: jpt2ch(:,:,:)
    integer :: NChan
    integer :: emax, e2max, e3max
  contains
    procedure :: init => InitNonOrthIsospinThreeBodySpace
    procedure :: fin => FinNonOrthIsospinThreeBodySpace
  end type NonOrthIsospinThreeBodySpace

  type :: MSpace
    type(Orbits) :: sps
    type(OrbitsIsospin) :: isps
    type(OneBodySpace) :: one
    type(TwoBodySpace) :: two
    !type(ThreeBodySpace) :: thr ! Orthogonal pn three-body
    type(NonOrthIsospinThreeBodySpace) :: thr
    logical :: is_constructed=.false.
    logical :: is_three_body =.false.
    real(8), allocatable :: NOCoef(:)
    real(8) :: hw, beta = 0.d0
    integer, allocatable :: holes(:)
    integer, allocatable :: particles(:)
    character(:), allocatable :: Nucl
    integer :: A, Z, N
    integer :: np, nh
    integer :: emax = -1
    integer :: e2max = -1
    integer :: e3max = -1
    integer :: lmax = -1
  contains
    procedure :: fin => FinMSpace
    procedure :: InitMSpace
    procedure :: InitMSpaceFromReference
    procedure :: InitMSpaceFromAZN
    procedure :: InitMSpaceFromFile
    generic :: init => InitMspace, InitMSPaceFromReference, InitMSpaceFromAZN, &
        & InitMSpaceFromFile
    procedure :: GetConfFromFile
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
    call this%isps%fin()
    call this%thr%fin()
  end subroutine FinMSpace

  subroutine InitMSpaceFromReference(this, Nucl, hw, emax, e2max, e3max, lmax, beta)
    class(MSpace), intent(inout) :: this
    character(*), intent(in) :: Nucl
    real(8), intent(in) :: hw
    integer, intent(in) :: emax, e2max
    integer, intent(in), optional :: e3max, lmax
    real(8), intent(in), optional :: beta
    integer :: A, Z, N

    call GetAZNFromReference(Nucl,A,Z,N)
    call this%InitMSpaceFromAZN(A,Z,N,hw,emax,e2max,e3max,lmax,beta)
  end subroutine InitMSpaceFromReference

  subroutine InitMSpaceFromAZN(this, A, Z, N, hw, emax, e2max, e3max, lmax, beta)
    use Profiler, only: timer
    use ClassSys, only: sys
    class(MSpace), intent(inout) :: this
    real(8), intent(in) :: hw
    integer, intent(in) :: A, Z, N
    integer, intent(in) :: emax, e2max
    integer, intent(in), optional :: lmax, e3max
    real(8), intent(in), optional :: beta
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
    this%hw = hw
    if(present(e3max)) this%e3max = e3max
    if(present(lmax)) this%lmax = lmax
    if(present(beta)) this%beta = beta
    write(*,'(a,i3,a,i3,a,i3)',advance='no') " emax=",this%emax,", e2max=",this%e2max
    if(present(e3max)) then
      write(*,'(a,i3)',advance='no') ", e3max=",this%e3max
    end if
    write(*,'(a,i3)') ", lmax=",this%lmax
    write(*,'(a,f6.2,a)') " hw = ",this%hw, " MeV"
    call this%sps%init(this%emax,this%lmax)
    call this%GetNOCoef()
    call this%GetParticleHoleOrbits()

    write(*,'(3a)') " Target Nuclide is ", trim(this%Nucl), ", Orbits:"
    write(*,'(a)') "      p/h, idx,  n,  l,  j, tz,   occupation"
    do i = 1, this%nh ! print hole states
      l = this%holes(i)
      write(*,'(a10,5i4,f14.6)') '     hole:', l, this%sps%orb(l)%n, this%sps%orb(l)%l, &
          & this%sps%orb(l)%j, this%sps%orb(l)%z, this%NOcoef(l)
    end do

    !do i = 1, this%np ! print particle states
    !  l = this%particles(i)
    !  write(*,'(a10,5i4,f14.6)') ' particle:', l, this%sps%orb(l)%n, this%sps%orb(l)%l, &
    !      & this%sps%orb(l)%j, this%sps%orb(l)%z, this%NOcoef(l)
    !end do
    write(*,*)
    call this%InitMSpace()
    write(*,*)
    call timer%countup_memory("Model Space")
    call timer%Add('Construct Model Space', omp_get_wtime()-ti)
  end subroutine InitMSpaceFromAZN

  subroutine InitMSpaceFromFile(this, Nucl, filename, hw, emax, e2max, e3max, lmax, beta)
    use Profiler, only: timer
    use ClassSys, only: sys
    class(MSpace), intent(inout) :: this
    character(*), intent(in) :: Nucl, filename
    real(8), intent(in) :: hw
    integer, intent(in) :: emax, e2max
    integer, intent(in), optional :: lmax, e3max
    real(8), intent(in), optional :: beta
    type(sys) :: s
    integer :: A, Z, N, i, l
    real(8) :: ti

    ti = omp_get_wtime()
    call timer%cmemory()
    write(*,*)
    write(*,'(a)') "### Model-space constructor ###"
    this%is_constructed = .true.
    this%emax = emax
    this%e2max = e2max
    this%lmax = emax
    this%hw = hw
    if(present(e3max)) this%e3max = e3max
    if(present(lmax)) this%lmax = lmax
    if(present(beta)) this%beta = beta
    write(*,'(a,i3,a,i3,a,i3)',advance='no') " emax=",this%emax,", e2max=",this%e2max
    if(present(e3max)) then
      write(*,'(a,i3)',advance='no') ", e3max=",this%e3max
    end if
    write(*,'(a,i3)') ", lmax=",this%lmax
    write(*,'(a,f6.2,a)') " hw = ",this%hw, " MeV"

    call this%sps%init(this%emax,this%lmax)
    call this%GetConfFromFile(filename, A, Z, N)
    call this%GetParticleHoleOrbits()
    this%A = A
    this%Z = Z
    this%N = N
    this%Nucl = trim(elements(this%Z+1)) // trim(s%str(this%A))
    if(Nucl /= this%Nucl) then
      write(*,'(a)') "Warning: make sure the target nuclide"
    end if
    write(*,'(3a)') " Target Nuclide is ", trim(this%Nucl), ", Orbits:"
    write(*,'(a)') "      p/h, idx,  n,  l,  j, tz,   occupation"
    do i = 1, this%nh ! print hole states
      l = this%holes(i)
      write(*,'(a10,5i4,f14.6)') '     hole:', l, this%sps%orb(l)%n, this%sps%orb(l)%l, &
          & this%sps%orb(l)%j, this%sps%orb(l)%z, this%NOcoef(l)
    end do
    write(*,*)
    call this%InitMSpace()
    write(*,*)

    call timer%countup_memory("Model Space")
    call timer%Add('Construct Model Space', omp_get_wtime()-ti)
  end subroutine InitMSpaceFromFile

  subroutine InitMSpace(this)
    class(MSpace), intent(inout) :: this
    call this%one%init(this%sps)
    call this%two%init(this%sps, this%e2max)
    write(*,'(a)') "  # of J, parity, and tz Channels:"
    write(*,'(a,i3)') "    OneBody: ", this%one%NChan
    write(*,'(a,i3)') "    TwoBody: ", this%two%NChan
    if(this%e3max == -1) return
    this%is_three_body = .true.
    call this%isps%init(this%emax)
    call this%thr%init(this%isps, this%e2max, this%e3max)
    write(*,'(a,i3)') "  ThreeBody: ", this%thr%NChan
  end subroutine InitMSpace

  subroutine GetConfFromFile(this, filename, A, Z, N)
    use CommonLibrary, only: skip_comment
    class(MSpace), intent(inout) :: this
    character(*), intent(in) :: filename
    integer, intent(out) :: A, Z, N
    integer :: nn, ll, jj, zz, idx, occ_num
    integer :: runit=20

    Z = 0
    N = 0
    allocate(this%NOCoef(this%sps%norbs))
    this%NOCoef(:) = 0.d0
    open(runit, file=filename, status='old', action='read')
    call skip_comment(runit, '!')
    do
      read(runit,*,end=999) nn, ll, jj, zz, occ_num
      idx = this%sps%nljz2idx(nn,ll,jj,zz)
      this%NOcoef(idx) = dble(occ_num)/dble(jj+1)
      if(zz == -1) Z = Z + occ_num
      if(zz ==  1) N = N + occ_num
    end do
999 close(runit)
    A = Z + N

    do idx = 1, this%sps%norbs
      call this%sps%orb(idx)%SetOccupation(this%NOcoef(idx))
      if(this%NOcoef(idx) < 1.d-6) then
        call this%sps%orb(idx)%SetParticleHole(1)
      elseif(abs(1.d0 - this%NOCoef(idx)) < 1.d-6) then
        call this%sps%orb(idx)%SetParticleHole(0)
      else
        call this%sps%orb(idx)%SetParticleHole(2)
      end if
    end do

  end subroutine GetConfFromFile

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

    do l = 1, this%sps%norbs
      call this%sps%orb(l)%SetOccupation(this%NOcoef(l))
      if(this%NOcoef(l) < 1.d-6) then
        call this%sps%orb(l)%SetParticleHole(1)
      elseif(abs(1.d0 - this%NOCoef(l)) < 1.d-6) then
        call this%sps%orb(l)%SetParticleHole(0)
      else
        call this%sps%orb(l)%SetParticleHole(2)
      end if
    end do

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

  function GetOneBodyChannel(this,ch) result(one)
    class(OneBodySpace), intent(in) :: this
    integer, intent(in) :: ch
    type(OneBodyChannel) :: one
    one = this%jpz(ch)
  end function GetOneBodyChannel

  function GetOneBodyChannelJPZ(this,J,P,Z) result(one)
    class(OneBodySpace), intent(in) :: this
    integer, intent(in) :: J,P,Z
    integer :: ch
    type(OneBodyChannel) :: one
    ch = this%jpz2ch(J,P,Z)
    one = this%jpz(ch)
  end function GetOneBodyChannelJPZ

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
    do j = 0, min(2*sps%lmax,e2max)+1
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
              e2 = sps%orb(i2)%e

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

  function GetTwoBodyChannel(this,ch) result(two)
    class(TwoBodySpace), intent(in) :: this
    integer, intent(in) :: ch
    type(TwoBodyChannel) :: two
    two = this%jpz(ch)
  end function GetTwoBodyChannel

  function GetTwoBodyChannelJPZ(this,J,P,Z) result(two)
    class(TwoBodySpace), intent(in) :: this
    integer, intent(in) :: J,P,Z
    integer :: ch
    type(TwoBodyChannel) :: two
    ch = this%jpz2ch(J,P,Z)
    two = this%jpz(ch)
  end function GetTwoBodyChannelJPZ


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
    integer :: a, b
    integer :: cnt

    this%j = j
    this%p = p
    this%z = z
    this%nst = n
    allocate(this%n2spi1(n))
    allocate(this%n2spi2(n))
    allocate(this%spis2n(sps%norbs,sps%norbs))
    allocate(this%iphase(sps%norbs,sps%norbs))
#ifdef ModelSpaceDebug
    write(*,'(a,i3,a,i3,a,i3)') "Two-body channel: J=", j, ", P=", p, ", Tz=", z
#endif
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
        if(triag(j1, j2, 2*j)) cycle
        if((-1) ** (l1+l2) /= p) cycle
        if(z1 + z2 /= 2*z) cycle
        if(i1 == i2 .and. mod(j,2) == 1) cycle

        a = i1; b = i2
        if( sps%orb(i1)%occ > sps%orb(i2)%occ ) then
          ! hp -> ph
          a = i2; b = i1
        end if

        cnt = cnt + 1
        this%n2spi1(cnt) = a
        this%n2spi2(cnt) = b
        this%spis2n(a,b) = cnt
        this%spis2n(b,a) = cnt
        this%iphase(a,b) = 1
        this%iphase(b,a) = -(-1) ** ((j1+j2)/2 - j)
#ifdef ModelSpaceDebug
        write(*,'(a,i3,a,i3,a,i6)') "i1=",a,", i2=",b,", Num=",cnt
#endif
      end do
    end do
#ifdef ModelSpaceDebug
    write(*,*)
#endif
  end subroutine InitTwoBodyChannel

  subroutine FinNonOrthIsospinThreeBodySpace(this)
    class(NonOrthIsospinThreeBodySpace), intent(inout) :: this
    integer :: ich
    do ich = 1, this%NChan
      call this%jpt(ich)%fin()
    end do
    deallocate(this%jpt)
    deallocate(this%jpt2ch)
  end subroutine FinNonOrthIsospinThreeBodySpace

  subroutine InitNonOrthIsospinThreeBodySpace(this, sps, e2max, e3max)
    use CommonLibrary, only: triag
    class(NonOrthIsospinThreeBodySpace), intent(inout) :: this
    type(OrbitsIsospin), intent(in) :: sps
    integer, intent(in) :: e2max, e3max
    integer :: j, p, t, ich, nidx, n
    integer :: i1, j1, l1, e1
    integer :: i2, j2, l2, e2
    integer :: i3, j3, l3, e3
    integer :: j12, t12, j12min, j12max, t12min, t12max
    integer :: nj
    integer, allocatable :: jj(:), pp(:), tt(:), nn(:), nnidx(:)
    type :: spis2n
      integer, allocatable :: spis2nidx(:,:,:)
    end type spis2n
    type(spis2n), allocatable :: ch(:)

    allocate(this%jpt2ch(1:2*min(e3max,3*sps%lmax)+3,-1:1,1:3))
    this%jpt2ch(:,:,:) = 0
    ich = 0
    do j = 1, 2*min(e3max,3*sps%lmax)+3, 2
      do p = 1, -1, -2
        do t = 1, 3, 2
          n = 0
          nidx = 0
          do i1 = 1, sps%norbs
            j1 = sps%orb(i1)%j
            l1 = sps%orb(i1)%l
            e1 = sps%orb(i1)%e
            do i2 = 1, i1
              j2 = sps%orb(i2)%j
              l2 = sps%orb(i2)%l
              e2 = sps%orb(i2)%e
              if(e1 + e2 > e2max) cycle
              do i3 = 1, i2
                j3 = sps%orb(i3)%j
                l3 = sps%orb(i3)%l
                e3 = sps%orb(i3)%e
                if(e1 + e3 > e2max) cycle
                if(e2 + e3 > e2max) cycle
                if(e1 + e2 + e3 > e3max) cycle
                if((-1) ** (l1+l2+l3) /= p) cycle

                nj = 0
                j12min = (j1+j2)/2
                j12max = abs(j1-j2)/2
                t12min = 1
                t12max = 0
                do j12 = abs(j1-j2)/2, (j1+j2)/2
                  if(triag(2*j12,j3,j)) cycle
                  do t12 = 0, 1
                    if(triag(2*t12, 1,t)) cycle
                    if(i1 == i2 .and. mod(j12+t12,2) == 0) cycle
                    n = n + 1
                    nj = nj + 1
                    j12max = max(j12, j12max)
                    j12min = min(j12, j12min)
                    t12max = max(t12, t12max)
                    t12min = min(t12, t12max)

                  end do
                end do
                if(nj /= 0) nidx = nidx + 1
              end do
            end do
          end do
          if(n /= 0) then
            ich = ich + 1
            this%jpt2ch(j,p,t) = ich
          end if

        end do
      end do
    end do
    this%NChan = ich
    this%emax = sps%emax
    this%e2max = e2max
    this%e3max = e3max

    allocate(this%jpt(ich))
    allocate(jj(ich))
    allocate(pp(ich))
    allocate(tt(ich))
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
        do t = 1, 3, 2
          n = 0
          nidx = 0
          do i1 = 1, sps%norbs
            j1 = sps%orb(i1)%j
            l1 = sps%orb(i1)%l
            e1 = sps%orb(i1)%e
            do i2 = 1, i1
              j2 = sps%orb(i2)%j
              l2 = sps%orb(i2)%l
              e2 = sps%orb(i2)%e
              if(e1 + e2 > e2max) cycle
              do i3 = 1, i2
                j3 = sps%orb(i3)%j
                l3 = sps%orb(i3)%l
                e3 = sps%orb(i3)%e
                if(e1 + e3 > e2max) cycle
                if(e2 + e3 > e2max) cycle
                if(e1 + e2 + e3 > e3max) cycle
                if((-1) ** (l1+l2+l3) /= p) cycle

                nj = 0
                j12min = (j1+j2)/2
                j12max = abs(j1-j2)/2
                t12min = 1
                t12max = 0
                do j12 = abs(j1-j2)/2, (j1+j2)/2
                  if(triag(2*j12,j3,j)) cycle
                  do t12 = 0, 1
                    if(triag(2*t12, 1,t)) cycle
                    if(i1 == i2 .and. mod(j12+t12,2) == 0) cycle
                    n = n + 1
                    nj = nj + 1
                    j12max = max(j12, j12max)
                    j12min = min(j12, j12min)
                    t12max = max(t12, t12max)
                    t12min = min(t12, t12max)

                  end do
                end do
                if(nj /= 0) then
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
            tt(ich) = t
            nn(ich) = n
            nnidx(ich) = nidx
          end if

        end do
      end do
    end do

    do ich = 1, this%NChan
      call this%jpt(ich)%init(jj(ich),pp(ich),tt(ich),nn(ich),nnidx(ich),&
          & ch(ich)%spis2nidx, sps, e2max, e3max)
      deallocate(ch(ich)%spis2nidx)
      call this%jpt(ich)%set( jj(ich),pp(ich),tt(ich), sps, e2max, e3max)
    end do
    deallocate(ch)
  end subroutine InitNonOrthIsospinThreeBodySpace

  subroutine FinNonOrthIsospinThreeBodyChannel(this)
    class(NonOrthIsospinThreeBodyChannel), intent(inout) :: this
    integer :: idx
    do idx = 1, this%n_idx
      deallocate(this%idxqn(idx)%JT2n)
    end do

    do idx = 1, this%n_sort
      call this%sort(idx)%fin()
    end do

    deallocate(this%sort)
    deallocate(this%idxqn)
    deallocate(this%n2spi1)
    deallocate(this%n2spi2)
    deallocate(this%n2spi3)
    deallocate(this%n2J12 )
    deallocate(this%n2T12 )
    deallocate(this%spis2idx)
    deallocate(this%sorting)
  end subroutine FinNonOrthIsospinThreeBodyChannel

  subroutine InitNonOrthIsospinThreeBodyChannel(this,j,p,t,n,nidx,spis2idx,sps,e2max,e3max)
    use CommonLibrary, only: triag
    class(NonOrthIsospinThreeBodyChannel), intent(inout) :: this
    integer, intent(in) :: j, p, t, n, nidx, spis2idx(:,:,:), e2max, e3max
    type(OrbitsIsospin), intent(in) :: sps
    integer :: j12max, j12min, t12max, t12min, j12, t12
    integer :: i1, j1, l1, e1
    integer :: i2, j2, l2, e2
    integer :: i3, j3, l3, e3
    integer :: nj, cnt, idx

    this%j = j
    this%p = p
    this%t = t
    this%nst = n
    this%n_idx = nidx
#ifdef ModelSpaceDebug
    write(*,'(a,i3,a,i3,a,i3,a)') "Three-body channel: J=", j, "/2, P=", p, ", T=", t, "/2"
#endif
    allocate(this%idxqn(this%n_idx))
    allocate(this%n2spi1(n))
    allocate(this%n2spi2(n))
    allocate(this%n2spi3(n))
    allocate(this%n2J12( n))
    allocate(this%n2T12( n))
    allocate(this%spis2idx(sps%norbs,sps%norbs,sps%norbs))
    this%spis2idx = spis2idx

    idx = 0
    do i1 = 1, sps%norbs
      j1 = sps%orb(i1)%j
      l1 = sps%orb(i1)%l
      e1 = sps%orb(i1)%e
      do i2 = 1, i1
        j2 = sps%orb(i2)%j
        l2 = sps%orb(i2)%l
        e2 = sps%orb(i2)%e
        if(e1 + e2 > e2max) cycle
        do i3 = 1, i2
          j3 = sps%orb(i3)%j
          l3 = sps%orb(i3)%l
          e3 = sps%orb(i3)%e
          if(e1 + e3 > e2max) cycle
          if(e2 + e3 > e2max) cycle
          if(e1 + e2 + e3 > e3max) cycle
          if((-1) ** (l1+l2+l3) /= p) cycle

          j12min = (j1+j2)/2
          j12max = abs(j1-j2)/2
          t12min = 1
          t12max = 0
          nj = 0
          do j12 = abs(j1-j2)/2, (j1+j2)/2
            if(triag(2*j12,j3,j)) cycle
            do t12 = 0, 1
              if(triag(2*t12, 1,t)) cycle
              if(i1 == i2 .and. mod(j12+t12,2) == 0) cycle
              nj = nj + 1
              j12max = max(j12, j12max)
              j12min = min(j12, j12min)
              t12max = max(t12, t12max)
              t12min = min(t12, t12min)

            end do
          end do

          if(nj /= 0) then
            idx = idx + 1
            this%idxqn(idx)%j12max = j12max
            this%idxqn(idx)%j12min = j12min
            this%idxqn(idx)%t12max = t12max
            this%idxqn(idx)%t12min = t12min
            allocate(this%idxqn(idx)%JT2n(j12min:j12max,t12min:t12max))
          end if
        end do
      end do
    end do

    cnt = 0
    idx = 0
    do i1 = 1, sps%norbs
      j1 = sps%orb(i1)%j
      l1 = sps%orb(i1)%l
      e1 = sps%orb(i1)%e
      do i2 = 1, i1
        j2 = sps%orb(i2)%j
        l2 = sps%orb(i2)%l
        e2 = sps%orb(i2)%e
        if(e1 + e2 > e2max) cycle
        do i3 = 1, i2
          j3 = sps%orb(i3)%j
          l3 = sps%orb(i3)%l
          e3 = sps%orb(i3)%e
          if(e1 + e3 > e2max) cycle
          if(e2 + e3 > e2max) cycle
          if(e1 + e2 + e3 > e3max) cycle
          if((-1) ** (l1+l2+l3) /= p) cycle

          j12min = (j1+j2)/2
          j12max = abs(j1-j2)/2
          t12min = 1
          t12max = 0
          nj = 0
          do j12 = abs(j1-j2)/2, (j1+j2)/2
            if(triag(2*j12,j3,j)) cycle
            do t12 = 0, 1
              if(triag(2*t12, 1,t)) cycle
              if(i1 == i2 .and. mod(j12+t12,2) == 0) cycle
              nj = nj + 1
              cnt = cnt + 1
              this%idxqn(idx+1)%JT2n(j12,t12) = cnt
            end do
          end do
          if(nj /= 0) then
            idx = idx + 1
          end if
        end do
      end do
    end do
  end subroutine InitNonOrthIsospinThreeBodyChannel

  subroutine SetSortingIndices(this,j,p,t,sps,e2max,e3max)
    use CommonLibrary, only: triag
    class(NonOrthIsospinThreeBodyChannel), intent(inout) :: this
    type(OrbitsIsospin), intent(in) :: sps
    integer, intent(in) :: j, p, t, e2max, e3max
    integer :: i1, l1, j1, e1
    integer :: i2, l2, j2, e2
    integer :: i3, l3, j3, e3
    integer :: a, b, c, n_recouple, idx_sorted
    integer :: cnt, idx
    integer :: j12min, j12, j12max
    integer :: t12min, t12, t12max

    allocate(this%sorting(sps%norbs,sps%norbs,sps%norbs))
    this%sorting(:,:,:) = 0
    idx = 0
    do i1 = 1, sps%norbs
      j1 = sps%orb(i1)%j
      l1 = sps%orb(i1)%l
      e1 = sps%orb(i1)%e
      do i2 = 1, sps%norbs
        j2 = sps%orb(i2)%j
        l2 = sps%orb(i2)%l
        e2 = sps%orb(i2)%e
        if(e1 + e2 > e2max) cycle
        do i3 = 1, sps%norbs
          j3 = sps%orb(i3)%j
          l3 = sps%orb(i3)%l
          e3 = sps%orb(i3)%e
          if(e1 + e3 > e2max) cycle
          if(e2 + e3 > e2max) cycle
          if(e1 + e2 + e3 > e3max) cycle
          if((-1) ** (l1+l2+l3) /= p) cycle
          idx = idx + 1
          this%sorting(i1,i2,i3) = idx
        end do
      end do
    end do

    this%n_sort = idx
    allocate(this%sort(this%n_sort))
    cnt = 0
    do i1 = 1, sps%norbs
      j1 = sps%orb(i1)%j
      l1 = sps%orb(i1)%l
      e1 = sps%orb(i1)%e
      do i2 = 1, sps%norbs
        j2 = sps%orb(i2)%j
        l2 = sps%orb(i2)%l
        e2 = sps%orb(i2)%e
        if(e1 + e2 > e2max) cycle
        do i3 = 1, sps%norbs
          j3 = sps%orb(i3)%j
          l3 = sps%orb(i3)%l
          e3 = sps%orb(i3)%e
          if(e1 + e3 > e2max) cycle
          if(e2 + e3 > e2max) cycle
          if(e1 + e2 + e3 > e3max) cycle
          if((-1) ** (l1+l2+l3) /= p) cycle

          call Sort123(i1,i2,i3,a,b,c,n_recouple)
          idx        = this%sorting(i1,i2,i3)
          idx_sorted = this%spis2idx(a,b,c)
          this%sort(idx)%idx_sorted = idx_sorted
          if(idx == 0) cycle
          if(idx_sorted == 0) cycle

          j12min = (j1+j2)/2
          j12max = abs(j1-j2)/2
          t12min = 1
          t12max = 0

          do j12 = iabs(j1-j2)/2, (j1+j2)/2
            if(triag(2*j12, j3, j)) cycle
            do t12 = 0, 1
              if(triag(2*t12,  1, t)) cycle
              if(i1 == i2 .and. mod(j12+t12, 2) == 0) cycle
              j12min = min(j12min, j12)
              j12max = max(j12max, j12)
              t12min = min(t12min, t12)
              t12max = max(t12max, t12)
            end do
          end do


          this%sort(idx)%j12min = j12min
          this%sort(idx)%j12max = j12max
          this%sort(idx)%t12min = t12min
          this%sort(idx)%t12max = t12max
          allocate(this%sort(idx)%jt(j12min:j12max,t12min:t12max))
          call this%sort(idx)%init(sps,i1,i2,i3,a,b,c,&
              & n_recouple,j,t,this%idxqn(idx_sorted))

        end do
      end do
    end do
  end subroutine SetSortingIndices

  subroutine Sort123(i1,i2,i3,ii1,ii2,ii3,n_recouple)
    integer, intent(in) :: i1, i2, i3
    integer, intent(out) :: ii1, ii2, ii3, n_recouple

    n_recouple = 0
    if(i1 >= i2 .and. i2 >= i3) n_recouple = 1                ! i1 >= i2 >= i3
    if(i1 >= i2 .and. i2 <  i3 .and. i1 >= i3) n_recouple = 2 ! i1 >= i3 >  i2
    if(i1 >= i2 .and. i2 <  i3 .and. i1 <  i3) n_recouple = 3 ! i3 >  i1 >= i2
    if(i1 <  i2 .and. i1 >= i3) n_recouple = 4                ! i2 >  i1 >= i3
    if(i1 <  i2 .and. i1 <  i3 .and. i2 >= i3) n_recouple = 5 ! i2 >= i3 >  i1
    if(i1 <  i2 .and. i1 <  i3 .and. i2 <  i3) n_recouple = 6 ! i3 >  i2 >  i1

    select case(n_recouple)
    case(1)
      ii1 = i1; ii2 = i2; ii3 = i3
    case(2)
      ii1 = i1; ii2 = i3; ii3 = i2
    case(3)
      ii1 = i3; ii2 = i1; ii3 = i2
    case(4)
      ii1 = i2; ii2 = i1; ii3 = i3
    case(5)
      ii1 = i2; ii2 = i3; ii3 = i1
    case(6)
      ii1 = i3; ii2 = i2; ii3 = i1
    case default
      write(*,*) "Error in Sort123"
      stop
    end select
  end subroutine Sort123

  subroutine fin_sort_index(this)
    class(sort_index), intent(inout) :: this
    integer :: j12, t12
    if(this%idx_sorted == 0) return
    do j12 = this%j12min, this%j12max
      do t12 = this%t12min, this%t12max
        if(this%JT(j12,t12)%n == 0) cycle
        deallocate(this%JT(j12,t12)%TrnsCoef)
        deallocate(this%JT(j12,t12)%Idx2num)
      end do
    end do
    deallocate(this%JT)
  end subroutine fin_sort_index

  subroutine init_sort_index(this, sps, i1, i2, i3, a, b, c, &
        &    n_recouple, j, t, NOIAQN)
    use CommonLibrary, only: triag, sjs, hat
    class(sort_index), intent(inout) :: this
    type(OrbitsIsospin), intent(in) :: sps
    integer, intent(in) :: i1, i2, i3, a, b, c, n_recouple, j, t
    type(NonOrthIsospinAdditionalQN), intent(in) :: NOIAQN
    integer :: n, loop
    integer :: j1, j2, j3, ja, jb, jc
    integer :: j12, t12, jab, tab

    j1 = sps%orb(i1)%j
    j2 = sps%orb(i2)%j
    j3 = sps%orb(i3)%j
    ja = sps%orb(a)%j
    jb = sps%orb(b)%j
    jc = sps%orb(c)%j
    do j12 = iabs(j1 - j2) / 2, (j1 + j2)/2
      do t12 = 0, 1
        if(i1 == i2 .and. mod(j12 + t12, 2) == 0) cycle
        if(triag(2 * t12,  1, t)) cycle
        if(triag(2 * j12, j3, j)) cycle
        do loop = 1, 2
          if(loop == 1) then
            n = 0
            if(n_recouple == 1 .or. n_recouple == 4) then
              n = 1
            else
              do jab = iabs(ja - jb)/2, (ja + jb)/2
                do tab = 0, 1
                  if(a == b .and. mod(jab + tab, 2) == 0) cycle
                  if(triag(2 * tab,  1, t)) cycle
                  if(triag(2 * jab, jc, j)) cycle
                  n = n + 1
                end do
              end do
            end if
            this%jt(j12, t12)%n = n
            allocate(this%jt(j12,t12)%idx2num( n))
            allocate(this%jt(j12,t12)%TrnsCoef(n))
            this%jt(j12,t12)%idx2num(:) = 0
#ifdef single_precision
            this%jt(j12,t12)%TrnsCoef(:) = 0.0
#else
            this%jt(j12,t12)%TrnsCoef(:) = 0.d0
#endif
          elseif(loop == 2) then
            n = 0
            select case(n_recouple)
            case(1)
              n = 1
#ifdef single_precision
              this%jt(j12,t12)%TrnsCoef(n) = 1.0
#else
              this%jt(j12,t12)%TrnsCoef(n) = 1.d0
#endif
              this%jt(j12,t12)%idx2num(n) = NOIAQN%JT2n(j12, t12)
            case(2)
              do jab = iabs(ja - jb)/2, (ja + jb)/2
                do tab = 0, 1
                  if(a == b .and. mod(jab + tab, 2) == 0) cycle
                  if(triag(2 * tab,  1, t)) cycle
                  if(triag(2 * jab, jc, j)) cycle
                  n = n + 1
#ifdef single_precision
                  this%jt(j12,t12)%TrnsCoef(n) = &
                      & real((-1.d0) ** ((jb + jc) / 2 + j12 + jab + t12 + tab) * &
                      & hat(2 * jab) * hat(2 * j12) * &
                      & sjs(ja, jb, 2 * jab, j, jc, 2 * j12) * &
                      & hat(2 * tab) * hat(2 * t12) * &
                      & sjs( 1,  1, 2 * tab, t,  1, 2 * t12))
#else
                  this%jt(j12,t12)%TrnsCoef(n) = &
                      & (-1.d0) ** ((jb + jc) / 2 + j12 + jab + t12 + tab) * &
                      & hat(2 * jab) * hat(2 * j12) * &
                      & sjs(ja, jb, 2 * jab, j, jc, 2 * j12) * &
                      & hat(2 * tab) * hat(2 * t12) * &
                      & sjs( 1,  1, 2 * tab, t,  1, 2 * t12)
#endif
                  this%jt(j12,t12)%idx2num(n) = NOIAQN%JT2n(jab, tab)
                end do
              end do
            case(3)
              do jab = iabs(ja - jb)/2, (ja + jb)/2
                do tab = 0, 1
                  if(a == b .and. mod(jab + tab, 2) == 0) cycle
                  if(triag(2 * tab,  1, t)) cycle
                  if(triag(2 * jab, jc, j)) cycle
                  n = n + 1
#ifdef single_precision
                  this%jt(j12,t12)%TrnsCoef(n) = real(&
                      & - (-1.d0) ** ((jb + jc) / 2 + j12 + t12) * &
                      & hat(2 * jab) * hat(2 * j12) * &
                      & sjs(ja, jb, 2 * jab, jc, j, 2 * j12) * &
                      & hat(2 * tab) * hat(2 * t12) * &
                      & sjs( 1,  1, 2 * tab,  1, t, 2 * t12))
#else
                  this%jt(j12,t12)%TrnsCoef(n) = &
                      & - (-1.d0) ** ((jb + jc) / 2 + j12 + t12) * &
                      & hat(2 * jab) * hat(2 * j12) * &
                      & sjs(ja, jb, 2 * jab, jc, j, 2 * j12) * &
                      & hat(2 * tab) * hat(2 * t12) * &
                      & sjs( 1,  1, 2 * tab,  1, t, 2 * t12)
#endif
                  this%jt(j12,t12)%idx2num(n) = NOIAQN%JT2n(jab, tab)
                end do
              end do
            case(4)
              n = n + 1
#ifdef single_precision
              this%jt(j12,t12)%TrnsCoef(n) = real(&
                  & (-1.d0) ** ((ja + jb) / 2 - j12 - t12))
#else
              this%jt(j12,t12)%TrnsCoef(n) = &
                  & (-1.d0) ** ((ja + jb) / 2 - j12 - t12)
#endif
              this%jt(j12,t12)%idx2num(n) = NOIAQN%JT2n(j12, t12)
            case(5)
              do jab = iabs(ja - jb)/2, (ja + jb)/2
                do tab = 0, 1
                  if(a == b .and. mod(jab + tab, 2) == 0) cycle
                  if(triag(2 * tab,  1, t)) cycle
                  if(triag(2 * jab, jc, j)) cycle
                  n = n + 1
#ifdef single_precision
                  this%jt(j12,t12)%TrnsCoef(n) = real(&
                      & - (-1.d0) ** ((ja + jb) / 2 + jab+ tab) * &
                      & hat(2 * jab) * hat(2 * j12) * &
                      & sjs(ja, jb, 2 * jab, j, jc, 2 * j12) * &
                      & hat(2 * tab) * hat(2 * t12) * &
                      & sjs( 1,  1, 2 * tab, t, 1, 2 * t12))
#else
                  this%jt(j12,t12)%TrnsCoef(n) = &
                      & - (-1.d0) ** ((ja + jb) / 2 + jab+ tab) * &
                      & hat(2 * jab) * hat(2 * j12) * &
                      & sjs(ja, jb, 2 * jab, j, jc, 2 * j12) * &
                      & hat(2 * tab) * hat(2 * t12) * &
                      & sjs( 1,  1, 2 * tab, t, 1, 2 * t12)
#endif
                  this%jt(j12,t12)%idx2num(n) = NOIAQN%JT2n(jab, tab)
                end do
              end do
            case(6)
              do jab = iabs(ja - jb)/2, (ja + jb)/2
                do tab = 0, 1
                  if(a == b .and. mod(jab + tab, 2) == 0) cycle
                  if(triag(2 * tab,  1, t)) cycle
                  if(triag(2 * jab, jc, j)) cycle
                  n = n + 1
#ifdef single_precision
                  this%jt(j12,t12)%TrnsCoef(n) = &
                      & real(- &
                      & hat(2 * jab) * hat(2 * j12) * &
                      & sjs(ja, jb, 2 * jab, jc, j, 2 * j12) * &
                      & hat(2 * tab) * hat(2 * t12) * &
                      & sjs( 1,  1, 2 * tab,  1, t, 2 * t12))
#else
                  this%jt(j12,t12)%TrnsCoef(n) = &
                      & - &
                      & hat(2 * jab) * hat(2 * j12) * &
                      & sjs(ja, jb, 2 * jab, jc, j, 2 * j12) * &
                      & hat(2 * tab) * hat(2 * t12) * &
                      & sjs( 1,  1, 2 * tab,  1, t, 2 * t12)
#endif
                  this%jt(j12,t12)%idx2num(n) = NOIAQN%JT2n(jab, tab)
                end do
              end do

            end select
          end if
        end do
        !write(*,'(5i3,a,3i3)') i1,i2,i3,j12,t12, ", ", a,b,c
        !write(*,'(10f10.4)') this%jt(j12,t12)%TrnsCoef

      end do
    end do
  end subroutine init_sort_index

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
    this%emax = sps%emax
    this%e2max = e2max
    this%e3max = e3max
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

  function GetThreeBodyChannel(this,ch) result(thr)
    class(ThreeBodySpace), intent(in) :: this
    integer, intent(in) :: ch
    type(ThreeBodyChannel) :: thr
    thr = this%jpz(ch)
  end function GetThreeBodyChannel

  function GetThreeBodyChannelJPZ(this,J,P,Z) result(thr)
    class(ThreeBodySpace), intent(in) :: this
    integer, intent(in) :: J,P,Z
    integer :: ch
    type(ThreeBodyChannel) :: thr
    ch = this%jpz2ch(J,P,Z)
    thr = this%jpz(ch)
  end function GetThreeBodyChannelJPZ

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
!program main
!  use Profiler, only: timer
!  use ModelSpace, only: MSpace
!  use CommonLibrary, only: &
!      &init_dbinomial_triangle, fin_dbinomial_triangle
!  type(MSpace) :: ms
!
!  call timer%init()
!  call init_dbinomial_triangle()
!  !call ms%init('O18', 20.d0, 4, 8, e3max=4, lmax=4)
!  call ms%init('O16', 20.d0, 14, 28, e3max=14)
!  call ms%fin()
!  call fin_dbinomial_triangle()
!  call timer%fin()
!end program main
