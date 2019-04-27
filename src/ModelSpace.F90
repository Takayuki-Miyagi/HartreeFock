! Define one-, two-, and three-body model space.
! The three-body space used in Hartree-Fock calculation is not orthonormalized
! Note that the orthonormalized three-body space is also defined
! with the coefficient of fractional parentage, which will not be used in Hartree-Fock calculation.
! However, this definition may be useful for beyond Hartree-Fock calculation methods.
!
module ModelSpace
  use omp_lib
  use SingleParticleState
  use OneBodyModelSpace
  use TwoBodyModelSpace
  use ThreeBodyModelSpace
  implicit none

  public :: MSpace

  private :: FinMSpace
  private :: InitMSpaceFromAZN
  private :: InitMSPaceFromReference
  private :: AssignCoreValence
  private :: GetParticleHoleOrbits
  private :: GetAZNFromReference
  private :: ReleaseThreeBody21
  private :: ReleaseThreeBody

  type :: MSpace
    type(Orbits) :: sps
    type(OrbitsIsospin) :: isps
    type(OneBodySpace) :: one
    type(TwoBodySpace) :: two
    type(ThreeBodySpace) :: thr
    type(NonOrthIsospinThreeBodySpace) :: thr21
    logical :: is_constructed=.false.
    logical :: is_three_body_jt =.false.
    logical :: is_three_body =.false.
    real(8), allocatable :: NOCoef(:)
    real(8) :: hw = 20.d0, beta = 0.d0
    integer, allocatable :: holes(:)
    integer, allocatable :: particles(:)
    character(:), allocatable :: Nucl
    character(:), allocatable :: Core
    integer :: A = 0, Z = 0, N = 0
    integer :: Ac= 0, Zc= 0, Nc= 0
    integer :: np, nh
    integer :: emax = -1
    integer :: e2max = -1
    integer :: e3max = -1
    integer :: lmax = -1
    integer :: e_fermi = -1
  contains
    procedure :: FinMSpace
    procedure :: InitMSpace
    procedure :: InitSubSpace
    procedure :: InitMSpaceFromReference
    procedure :: InitMSpaceFromAZN
    procedure :: InitMSpaceFromFile
    procedure :: GetConfFromFile
    procedure :: AssignCoreValence
    procedure :: GetParticleHoleOrbits
    procedure :: SetEFermi
    procedure :: GetEFermi
    procedure :: ReleaseThreeBody21
    procedure :: ReleaseThreeBody

    generic :: init => InitMspace
    generic :: fin => FinMSpace
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
    if(this%is_three_body_jt) then
      call this%thr21%fin()
    end if

    if(this%is_three_body) then
      call this%thr%fin()
    end if
  end subroutine FinMSpace

  subroutine InitMSpace(this, Nucl, Core, valence_orbits, filename, &
        & hw, emax, e2max, e3max, lmax, beta, is_three_body_jt, is_three_body)
    use ClassSys, only: sys
    class(MSpace), intent(inout) :: this
    character(*), intent(in), optional :: Core, Nucl, valence_orbits, filename
    real(8), intent(in), optional :: hw, beta
    integer, intent(in), optional :: emax, e2max, e3max, lmax
    logical, intent(in), optional :: is_three_body_jt, is_three_body

    if(present(filename)) then
      call this%InitMspaceFromFile(filename, hw, emax, e2max, &
          & e3max, lmax, beta, is_three_body_jt, is_three_body)
      return
    end if

    if(present(Nucl)) then
      call this%InitMSpaceFromReference(Nucl, Core, valence_orbits, &
          & hw, emax, e2max, e3max, lmax, beta, is_three_body_jt, &
          & is_three_body)
      return
    end if

    write(*,"(a)") "Error: input of InitMSpace"
    stop
  end subroutine InitMSpace

  subroutine InitMSpaceFromReference(this, Nucl, Core, valence_orbits, &
        & hw, emax, e2max, e3max, lmax, beta, is_three_body_jt, is_three_body)
    class(MSpace), intent(inout) :: this
    character(*), intent(in) :: Nucl
    character(*), intent(in), optional :: Core, valence_orbits
    real(8), intent(in), optional :: hw
    integer, intent(in), optional :: emax, e2max
    integer, intent(in), optional :: e3max, lmax
    real(8), intent(in), optional :: beta
    logical, intent(in), optional :: is_three_body_jt, is_three_body
    integer :: A=0, Z=0, N=0, Ac=0, Zc=0, Nc=0

    call GetAZNFromReference(Nucl,A,Z,N)
    if(.not. present(Core)) then
      Ac = A; Zc = Z; Nc = N
    end if

    if(present(Core)) then
      call GetAZNFromReference(Core,Ac,Zc,Nc)
    end if

    call this%InitMSpaceFromAZN([A,Z,N],[Ac,Zc,Nc],valence_orbits, &
        & hw,emax,e2max,e3max,lmax,beta,is_three_body_jt, is_three_body)
  end subroutine InitMSpaceFromReference

  subroutine InitMSpaceFromAZN(this, Nucl, Core, valence_orbits, hw, emax, e2max, e3max, &
        & lmax, beta, is_three_body_jt, is_three_body)
    use Profiler, only: timer
    use ClassSys, only: sys
    class(MSpace), intent(inout), target :: this
    character(*), intent(in), optional :: valence_orbits
    real(8), intent(in), optional :: hw
    integer, intent(in), optional :: Core(3), Nucl(3)
    integer, intent(in), optional :: emax, e2max
    integer, intent(in), optional :: lmax, e3max
    real(8), intent(in), optional :: beta
    logical, intent(in), optional :: is_three_body_jt, is_three_body
    type(sys) :: s
    integer :: i
    type(SingleParticleOrbit), pointer :: o
    real(8) :: ti

    ti = omp_get_wtime()
    call timer%cmemory()
    write(*,*)
    write(*,'(a)') "### Model-space constructor ###"


    this%is_constructed = .true.
    this%Ac = 0; this%Zc = 0; this%Nc = 0
    this%A = 0; this%Z = 0; this%N = 0

    if(present(Core)) then
      this%Ac = Core(1); this%Zc = Core(2); this%Nc = Core(3)
    end if

    if(present(Nucl)) then
      this%A = Nucl(1); this%Z = Nucl(2); this%N = Nucl(3)
    end if

    this%Core = adjustl(trim(elements(this%Zc+1)) // trim(s%str(this%Ac)))
    this%Nucl = adjustl(trim(elements(this%Z+1)) // trim(s%str(this%A)))
    if(present(hw)) this%hw = hw
    if(present(emax)) this%emax = emax
    if(present(e2max)) this%e2max = e2max
    this%lmax = emax
    if(present(e3max)) this%e3max = e3max
    if(present(lmax)) this%lmax = lmax
    if(present(beta)) this%beta = beta
    if(present(is_three_body)) this%is_three_body = is_three_body
    if(present(is_three_body_jt)) this%is_three_body_jt = is_three_body_jt
    write(*,'(a,i3,a,i3,a,i3)',advance='no') " emax=",this%emax,", e2max=",this%e2max
    if(present(e3max)) then
      write(*,'(a,i3)',advance='no') ", e3max=",this%e3max
    end if
    write(*,'(a,i3)') ", lmax=",this%lmax
    write(*,'(a,f6.2,a)') " hw = ",this%hw, " MeV"
    call this%sps%init(this%emax,this%lmax)
    call this%isps%init(this%emax,this%lmax)
    call this%AssignCoreValence(valence_orbits)
    call this%GetParticleHoleOrbits()
    call this%SetEFermi()

    write(*,'(2a)') " Target Nuclide is ", trim(this%Nucl)
    write(*,'(2a)') "   Core Nuclide is ", trim(this%Core)
    write(*,'(a)') "      p/h, idx,  n,  l,  j, tz,   occupation"
    do i = 1, this%sps%norbs ! print hole states
      o => this%sps%GetOrbit(i)
      if(o%ph /= 0) cycle
      write(*,'(a10,5i4,f14.6)') '     hole:', i, o%n, o%l, o%j, o%z, o%occ
    end do
    do i = 1, this%sps%norbs ! print valence states
      o => this%sps%GetOrbit(i)
      if(o%ph /= 2) cycle
      write(*,'(a10,5i4,f14.6)') '  valence:', i, o%n, o%l, o%j, o%z, o%occ
    end do

    write(*,*)
    call this%InitSubSpace()
    write(*,*)
    call timer%countup_memory("Model Space")
    call timer%Add('Construct Model Space', omp_get_wtime()-ti)
  end subroutine InitMSpaceFromAZN

  subroutine InitMSpaceFromFile(this, filename, hw, emax, e2max, e3max, &
        & lmax, beta, is_three_body_jt, is_three_body)
    use Profiler, only: timer
    use ClassSys, only: sys
    class(MSpace), intent(inout), target :: this
    character(*), intent(in) :: filename
    real(8), intent(in), optional :: hw
    integer, intent(in), optional :: emax, e2max
    integer, intent(in), optional :: lmax, e3max
    real(8), intent(in), optional :: beta
    logical, intent(in), optional :: is_three_body_jt, is_three_body
    type(sys) :: s
    integer :: Ac, Zc, Nc, A, Z, N, i, l
    type(SingleParticleOrbit), pointer :: o
    real(8) :: ti

    ti = omp_get_wtime()
    call timer%cmemory()
    write(*,*)
    write(*,'(a)') "### Model-space constructor ###"
    this%is_constructed = .true.

    if(present(hw)) this%hw = hw
    if(present(emax)) this%emax = emax
    if(present(e2max)) this%e2max = e2max
    this%lmax = emax
    if(present(e3max)) this%e3max = e3max
    if(present(lmax)) this%lmax = lmax
    if(present(beta)) this%beta = beta
    if(present(is_three_body)) this%is_three_body = is_three_body
    if(present(is_three_body_jt)) this%is_three_body_jt = is_three_body_jt
    write(*,'(a,i3,a,i3,a,i3)',advance='no') " emax=",this%emax,", e2max=",this%e2max
    if(present(e3max)) then
      write(*,'(a,i3)',advance='no') ", e3max=",this%e3max
    end if
    write(*,'(a,i3)') ", lmax=",this%lmax
    write(*,'(a,f6.2,a)') " hw = ",this%hw, " MeV"

    call this%sps%init(this%emax,this%lmax)
    call this%isps%init(this%emax,this%lmax)
    call this%GetConfFromFile(filename, Ac, Zc, Nc, A, Z, N)
    call this%GetParticleHoleOrbits()
    call this%SetEFermi()

    this%A = A
    this%Z = Z
    this%N = N
    this%Ac = Ac
    this%Zc = Zc
    this%Nc = Nc
    this%Core = trim(elements(this%Zc+1)) // trim(s%str(this%Ac))
    this%Nucl = trim(elements(this%Z+1)) // trim(s%str(this%A))
    write(*,'(2a)') " Target Nuclide is ", trim(this%Nucl)
    write(*,'(2a)') "   Core Nuclide is ", trim(this%Core)
    write(*,'(a)') "      p/h, idx,  n,  l,  j, tz,   occupation"
    do i = 1, this%sps%norbs ! print hole states
      o => this%sps%GetOrbit(i)
      if(o%ph /= 0) cycle
      write(*,'(a10,5i4,f14.6)') '     hole:', i, o%n, o%l, o%j, o%z, this%NOcoef(i)
    end do
    do i = 1, this%sps%norbs ! print valence states
      o => this%sps%GetOrbit(i)
      if(o%ph /= 2) cycle
      write(*,'(a10,5i4,f14.6)') '  valence:', i, o%n, o%l, o%j, o%z, this%NOcoef(i)
    end do
    write(*,*)
    call this%InitSubSpace()
    write(*,*)

    call timer%countup_memory("Model Space")
    call timer%Add('Construct Model Space', omp_get_wtime()-ti)
  end subroutine InitMSpaceFromFile

  subroutine InitSubSpace(this)
    class(MSpace), intent(inout) :: this
    call this%one%init(this%sps)
    write(*,'(a)') "  # of J, parity, and tz Channels:"
    write(*,'(a,i3)') "    OneBody: ", this%one%NChan
    call this%two%init(this%sps, this%e2max)
    write(*,'(a,i3)') "    TwoBody: ", this%two%NChan

    if(this%e3max > 0 .and. this%is_three_body_jt) then
      call this%thr21%init(this%sps, this%isps, this%e2max, this%e3max)
      write(*,'(a,i3)') "  ThreeBody (2+1): ", this%thr21%NChan
    end if

    if(this%e3max > 0 .and. this%is_three_body) then
      call this%thr%init(this%sps, this%e2max, this%e3max)
      write(*,'(a,i3)') "  ThreeBody: ", this%thr%NChan
    end if
  end subroutine InitSubSpace

  subroutine GetConfFromFile(this, filename, Ac, Zc, Nc, A, Z, N)
    use MyLibrary, only: skip_comment
    class(MSpace), intent(inout), target :: this
    character(*), intent(in) :: filename
    character(20) :: cvp
    integer, intent(out) :: Ac, Zc, Nc, A, Z, N
    integer :: nn, ll, jj, zz, idx, occ_num
    integer :: runit=20
    type(SingleParticleOrbit), pointer :: o

    Z = 0
    N = 0
    Zc = 0
    Nc = 0
    allocate(this%NOCoef(this%sps%norbs))
    this%NOCoef(:) = 0.d0
    open(runit, file=filename, status='old', action='read')
    call skip_comment(runit, '!')
    do
      read(runit,*,end=999) nn, ll, jj, zz, cvp, occ_num
      idx = this%sps%nljz2idx(nn,ll,jj,zz)
      this%NOcoef(idx) = dble(occ_num)/dble(jj+1)
      o => this%sps%GetOrbit(idx)
      if(zz == -1) Z = Z + occ_num
      if(zz ==  1) N = N + occ_num
      if(cvp(1:1) == "c" .and. zz == -1) Zc = Zc + occ_num
      if(cvp(1:1) == "c" .and. zz ==  1) Nc = Nc + occ_num
      if(cvp(1:1) == "c") call o%SetHoleParticleValence(0)
      if(cvp(1:1) == "v") call o%SetHoleParticleValence(2)
      if(cvp(1:1) == "p") call o%SetHoleParticleValence(1)
    end do
999 close(runit)
    A = Z + N
    Ac = Zc + Nc

    do idx = 1, this%sps%norbs
      o => this%sps%orb(idx)
      call o%SetOccupation(this%NOcoef(idx))
    end do

    do idx = 1, this%sps%norbs
      o => this%sps%orb(idx)
      if(o%ph == 0) cycle
      if(o%ph == 2) cycle
      call o%SetHoleParticleValence(1)
    end do
  end subroutine GetConfFromFile

  subroutine AssignCoreValence(this, valence_orbits)
    ! This should be called after obtaining A, Z, N
    use ClassSys, only: sys
    class(MSpace), intent(inout), target :: this
    character(*), intent(in), optional :: valence_orbits
    character(256), allocatable :: v_orbits(:)
    character(:), allocatable :: vlabel
    integer :: Z, N
    integer :: e, l, j, g, ns
    integer :: zz, nn, vz, vn
    integer :: idxp, idxn, idx
    type(SingleParticleOrbit), pointer :: o
    type(sys) :: s

    allocate(this%NOCoef(this%sps%norbs))
    Z = this%Zc
    N = this%Nc
    l = 0
    zz = 0
    do e = 0, this%sps%emax
      do g = 2*e+1, -2*e+1, -4
        j = abs(g)
        if(g > 0) l = (j-1)/2
        if(g < 0) l = (j+1)/2
        ns = (e-l)/2
        if(l > this%sps%lmax) cycle
        if(Z - zz == 0) cycle
        zz = zz + (j+1)
        idxp = this%sps%nljz2idx(ns,l,j,-1)
        o => this%sps%GetOrbit(idxp)

        call o%SetHoleParticleValence(0)
        if(Z - zz < 0) then
          write(*,"(a, i4)") "Error, proton core configuration has to be close:", Z
          stop
        end if
      end do
    end do

    nn = 0
    do e = 0, this%sps%emax
      do g = 2*e+1, -2*e+1, -4
        j = abs(g)
        if(g > 0) l = (j-1)/2
        if(g < 0) l = (j+1)/2
        ns = (e-l)/2
        if(l > this%sps%lmax) cycle
        if(N - nn == 0) cycle
        nn = nn + (j+1)
        idxn = this%sps%nljz2idx(ns,l,j, 1)
        o => this%sps%GetOrbit(idxn)

        call o%SetHoleParticleValence(0)
        if(N - nn < 0) then
          write(*,"(a, i4)") "Error, neutron core configuration has to be close:", N
          stop
        end if
      end do
    end do


    if(Z - zz > 0 .or. N - nn > 0) then
      write(*,'(a)', advance='no') "Error: emax is too small: "
      write(*,'(a,i4,a,i3,a,i3)') "Ac = ", this%Ac, ", Zc = ", this%Zc, ", Nc = ", this%Nc
      stop
    end if

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
      o => this%sps%orb(l)
      call o%SetOccupation(this%NOcoef(l))
    end do

    if(present(valence_orbits)) then
      call s%split(valence_orbits, ",", v_orbits)
      if(v_orbits(1) == "" .or. v_orbits(1) == "none") then
        do l = 1, this%sps%norbs
          o => this%sps%orb(l)
          if(o%ph == 0) cycle
          if(abs(1.d0 - o%occ) > 1.d-8 .and. abs(o%occ) > 1.d-8) then
            write(*,"(a)") "Error: no valence orbits and fractional occupation"
            stop
          end if
        end do
      end if
    end if

    if(allocated(v_orbits)) then
      do l = 1, size(v_orbits)
        vlabel = v_orbits(l)
        if(vlabel == "" .or. vlabel == "none") cycle
        idx = this%sps%GetIndexFromLabel(vlabel)
        o => this%sps%orb(idx)
        if(o%ph == 0) then
          write(*,"(2a)") "Error: conflict occurs core and valence orbit", trim(vlabel)
        end if
        call o%SetHoleParticleValence(2)
      end do
    end if

    do l = 1, this%sps%norbs
      o => this%sps%orb(l)
      if(o%ph == 0) cycle
      if(o%ph == 2) cycle
      call o%SetHoleParticleValence(1)
    end do

#ifdef ModelSpaceDebug
    write(*,'(a)') "In AssignCoreValence:"
    write(*,'(a)') "   n,  l,  j, tz, cpv,   occupation"
    do l = 1, this%sps%norbs
      o => this%sps%orb(l)
      if(o%ph == 0) vlabel = "c"
      if(o%ph == 2) vlabel = "v"
      if(o%ph == 1) vlabel = "p"
      write(*,'(4i4,a5,f14.6)') o%n, o%l, o%j, o%z, vlabel, this%NOcoef(l)
    end do
#endif
  end subroutine AssignCoreValence

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
    character(*), intent(in), optional :: Nucl
    integer, intent(out) :: A, Z, N
    integer :: i
    character(:), allocatable :: str
    character(20) :: str_A, element

    if(.not. present(Nucl) .or. Nucl == "") then
      A = 0
      Z = 0
      N = 0
      return
    end if

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

  subroutine SetEFermi(this)
    class(MSpace), intent(inout) :: this
    type(SingleParticleOrbit), pointer :: o
    integer :: i
    do i = 1, this%sps%norbs
      o => this%sps%GetOrbit(i)
      if(o%ph == 1) cycle
      this%e_fermi = max(this%e_fermi, o%e)
    end do
  end subroutine SetEFermi

  function GetEFermi(this) result(e_fermi)
    class(MSpace), intent(inout) :: this
    integer ::e_fermi
    e_fermi = this%e_fermi
  end function GetEFermi

  subroutine ReleaseThreeBody21(this)
    class(MSpace), intent(inout) :: this
    if(.not. this%is_three_body_jt) return
    call this%thr21%fin()
    this%is_three_body_jt = .false.
  end subroutine ReleaseThreeBody21

  subroutine ReleaseThreeBody(this)
    class(MSpace), intent(inout) :: this
    if(.not. this%is_three_body) return
    call this%thr%fin()
    this%is_three_body = .false.
  end subroutine ReleaseThreeBody
end module ModelSpace

! main for test
!program main
!  use Profiler, only: timer
!  use ModelSpace, only: MSpace
!  type(MSpace) :: ms
!
!  call timer%init()
!  !call ms%init('O18', 20.d0, 4, 8, e3max=4, lmax=4)
!  call ms%init('O16', 20.d0, 6, 6, e3max=6, is_three_body_jt=.true., is_three_body=.true.)
!  call ms%fin()
!  call timer%fin()
!end program main
