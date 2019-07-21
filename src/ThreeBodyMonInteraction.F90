module ThreeBodyMonInteraction
  use omp_lib
  use Profiler, only: timer
  use SingleParticleState
  implicit none

  public :: ThreeBodyMonForce

  private :: InitOneBodyChannels
  private :: FinOneBodyChannels
  private :: InitThreeBodyMonChannel
  private :: FinThreeBodyMonChannel
  private :: InitThreeBodyMonSpace
  private :: FinThreeBodyMonSpace

  type :: OneBodyChannels
    integer :: emax
    integer :: NChan
    integer, allocatable :: j(:)
    integer, allocatable :: p(:)
    integer, allocatable :: jp2ch(:,:)
  contains
    procedure :: init => InitOneBodyChannels
    procedure :: fin => FinOneBodyChannels
  end type OneBodyChannels

  type :: ThreeBodyMonChannel
    integer :: j, p, t
    integer :: ch1, ch2, ch3, j12
    integer :: n_state
    integer, allocatable :: n2abct(:,:)
    integer, allocatable :: nst2n(:,:,:,:)
  contains
    procedure :: init => InitThreeBodyMonChannel
    procedure :: fin => FinThreeBodyMonChannel
  end type ThreeBodyMonChannel

  type :: ThreeBodyMonSpace
    type(ThreeBodyMonChannel), allocatable :: Chan(:)
    type(OrbitsIsospin), pointer :: sps
    type(OneBodyChannels) :: one
    integer, allocatable :: TMon2Ch(:,:)
    integer, allocatable :: T2Ch(:)
    integer, allocatable :: Mon2Ch(:,:,:)
    integer :: emax, e2max, e3max
    integer :: NChan, NMon, NT
    logical :: is_Constructed=.false.
  contains
    procedure :: init => InitThreeBodyMonSpace
    procedure :: fin => FinThreeBodyMonSpace
  end type ThreeBodyMonSpace

  type :: ThreeBodyMonForceChannel
#ifdef single_precision_three_body_file
    real(4), allocatable :: v(:,:)
#else
    real(8), allocatable :: v(:,:)
#endif
  end type ThreeBodyMonForceChannel

  type :: ThreeBodyMonForce
    type(ThreeBodyMonForceChannel), allocatable :: MatCh(:)
    type(ThreeBodyMonSpace) :: thr
    type(Orbits), pointer :: sps
    type(OrbitsIsospin), pointer :: isps
    real(8), allocatable :: Cgs(:,:,:,:,:,:)
    integer :: emax, e2max, e3max
    logical :: zero = .true.
  contains
    procedure :: InitThreeBodyMonForce
    procedure :: FinThreeBodyMonForce
    procedure :: GetThBMEMon_isospin
    procedure :: GetThBMEMon_pn

    generic :: init => InitThreeBodyMonForce
    generic :: fin => FinThreeBodyMonForce
    generic :: GetMonThBME => GetThBMEMon_isospin, GetThBMEMon_pn
  end type ThreeBodyMonForce

  type :: Read3BodyMonopole
    character(:), allocatable :: file_3n
    integer :: emax3=-1, e2max3=-1, e3max3=-1, lmax3=-1
  contains
    procedure :: Set3BodyMonopoleFile          ! setter
    procedure :: Set3BodyMonopoleFileBoundaries! setter
    generic :: set => Set3BodyMonopoleFile, Set3BodyMonopoleFileBoundaries
    procedure :: ReadThreeBodyMonopole
    ! methods for three-body marix element
    procedure :: ReadFile
    procedure :: read_scalar_me3j_ascii_txt
    procedure :: read_scalar_me3j_gzip
    procedure :: read_scalar_me3j_ascii
    !procedure :: read_scalar_me3j_binary_stream
    !procedure :: read_scalar_me3j_binary_comp
    !procedure :: read_scalar_me3j_binary
  end type Read3BodyMonopole
contains
  subroutine FinOneBodyChannels(this)
    class(OneBodyChannels) :: this
    deallocate(this%j)
    deallocate(this%p)
    deallocate(this%jp2ch)
  end subroutine FinOneBodyChannels

  subroutine InitOneBodyChannels(this, sps)
    class(OneBodyChannels) :: this
    type(OrbitsIsospin), intent(in) :: sps
    integer :: la, ja, cnt
    this%emax = sps%emax
    allocate(this%jp2ch(2*sps%emax+1,-1:1))
    cnt = 0
    do la = 0, sps%emax
      do ja = abs(2*la-1), 2*la+1
        cnt = cnt + 1
      end do
    end do
    this%NChan = cnt
    allocate(this%j(this%NChan))
    allocate(this%p(this%NChan))

    cnt = 0
    do la = 0, sps%emax
      do ja = abs(2*la-1), 2*la+1
        cnt = cnt + 1
        this%jp2ch(ja,(-1)**la) = cnt
        this%j(cnt) = ja
        this%p(cnt) = (-1)**la
      end do
    end do
  end subroutine InitOneBodyChannels

  subroutine FinThreeBodyMonSpace(this)
    class(ThreeBodyMonSpace), intent(inout) :: this
    integer :: ich
    if(.not. this%is_Constructed) return
    do ich = 1, this%NChan
      call this%chan(ich)%fin()
    end do
    deallocate(this%TMon2Ch)
    deallocate(this%T2Ch)
    deallocate(this%Mon2Ch)
    deallocate(this%Chan)
    call this%one%fin()
    this%sps => null()
    this%is_Constructed=.false.
  end subroutine FinThreeBodyMonSpace

  subroutine InitThreeBodyMonSpace(this, sps, e2max, e3max)
    use MyLibrary, only: triag
    class(ThreeBodyMonSpace), intent(inout) :: this
    type(OrbitsIsospin), target, intent(in) :: sps
    integer, intent(in) :: e2max, e3max
    integer :: ttot
    integer :: cnt
    integer :: ch1, ch2, ch3
    integer :: j1, j2, j3, p1, p2, p3
    integer :: i1, i2, i3, t12
    type(SingleParticleOrbitIsospin), pointer :: o1, o2, o3
    integer :: ich, ich_T, ich_Mon
    integer, allocatable :: pp(:), tt(:)
    integer, allocatable :: ch1_(:), ch2_(:), ch3_(:)
    this%sps => sps
    this%emax = sps%emax
    this%e2max = e2max
    this%e3max = e3max

    call this%one%init(sps)
    this%NT = 2
    this%NMon = this%one%NChan**3
    allocate(this%Mon2Ch(this%one%NChan, this%one%NChan, this%one%NChan))
    allocate(this%T2Ch(1:3))
    allocate(this%TMon2Ch(this%NT, this%NMon))
    this%Mon2Ch(:,:,:) = 0
    this%T2Ch(:) = 0
    this%TMon2Ch(:,:) = 0
    ich = 0
    do ttot = 1, 3, 2

      do ch1 = 1, this%one%NChan
        j1 = this%one%j(ch1)
        p1 = this%one%p(ch1)
        do ch2 = 1, this%one%NChan
          j2 = this%one%j(ch2)
          p2 = this%one%p(ch2)
          do ch3 = 1, this%one%NChan
            j3 = this%one%j(ch3)
            p3 = this%one%p(ch3)

            cnt = 0
            do i1 = 1, sps%norbs
              o1 => sps%orb(i1)
              if(j1 /= o1%j) cycle
              if(p1 /= (-1)**o1%l) cycle
              do i2 = 1, sps%norbs
                o2 => sps%orb(i2)
                if(j2 /= o2%j) cycle
                if(p2 /= (-1)**o2%l) cycle
                if(o1%e + o2%e > e2max) cycle
                do i3 = 1, sps%norbs
                  o3 => sps%orb(i3)
                  if(j3 /= o3%j) cycle
                  if(p3 /= (-1)**o3%l) cycle
                  if(o1%e + o3%e > e2max) cycle
                  if(o2%e + o3%e > e2max) cycle
                  if(o1%e + o2%e + o3%e > e3max) cycle
                  do t12 = 0, 1
                    if(triag(2*t12, 1, ttot)) cycle
                    cnt = cnt + 1
                  end do
                end do
              end do
            end do

            if(cnt /= 0) ich = ich + 1
          end do
        end do
      end do
    end do

    this%NChan = ich
    allocate(pp(this%NChan))
    allocate(tt(this%NChan))
    allocate(ch1_(this%NChan))
    allocate(ch2_(this%NChan))
    allocate(ch3_(this%NChan))

    ich = 0
    ich_T = 0
    do ttot = 1, 3, 2
      ich_T = ich_T+1

      ich_Mon = 0
      do ch1 = 1, this%one%NChan
        j1 = this%one%j(ch1)
        p1 = this%one%p(ch1)
        do ch2 = 1, this%one%NChan
          j2 = this%one%j(ch2)
          p2 = this%one%p(ch2)
          do ch3 = 1, this%one%NChan
            j3 = this%one%j(ch3)
            p3 = this%one%p(ch3)
            ich_Mon = ich_Mon + 1

            cnt = 0
            do i1 = 1, sps%norbs
              o1 => sps%orb(i1)
              if(j1 /= o1%j) cycle
              if(p1 /= (-1)**o1%l) cycle
              do i2 = 1, sps%norbs
                o2 => sps%orb(i2)
                if(j2 /= o2%j) cycle
                if(p2 /= (-1)**o2%l) cycle
                if(o1%e + o2%e > e2max) cycle
                do i3 = 1, sps%norbs
                  o3 => sps%orb(i3)
                  if(j3 /= o3%j) cycle
                  if(p3 /= (-1)**o3%l) cycle
                  if(o1%e + o3%e > e2max) cycle
                  if(o2%e + o3%e > e2max) cycle
                  if(o1%e + o2%e + o3%e > e3max) cycle
                  do t12 = 0, 1
                    if(triag(2*t12, 1, ttot)) cycle
                    cnt = cnt + 1
                  end do
                end do
              end do
            end do
            if(cnt /= 0) then
              ich = ich + 1
              this%TMon2Ch(ich_T, ich_Mon) = ich
              this%T2Ch(ttot) = ich_T
              this%Mon2Ch(ch1, ch2, ch3) = ich_Mon
              tt(ich) = ttot
              ch1_(ich) = ch1
              ch2_(ich) = ch2
              ch3_(ich) = ch3
            end if
          end do
        end do
      end do
    end do

    allocate(this%chan(this%NChan))
    do ich = 1, this%NChan
      call this%chan(ich)%init(sps, this%one, tt(ich), &
          & ch1_(ich), ch2_(ich), ch3_(ich), e2max, e3max)
    end do
    deallocate(tt)
    deallocate(ch1_)
    deallocate(ch2_)
    deallocate(ch3_)
    this%is_Constructed=.true.
  end subroutine InitThreeBodyMonSpace

  subroutine FinThreeBodyMonChannel(this)
    class(ThreeBodyMonChannel), intent(inout) :: this
    deallocate(this%n2abct)
    deallocate(this%nst2n)
  end subroutine FinThreeBodyMonChannel

  subroutine InitThreeBodyMonChannel(this, sps, one, t, ch1, ch2, ch3, e2max, e3max)
    use MyLibrary, only: triag
    class(ThreeBodyMonChannel), intent(inout) :: this
    type(OrbitsIsospin), target, intent(in) :: sps
    type(OneBodyChannels), intent(in) :: one
    integer, intent(in) :: t, ch1, ch2, ch3, e2max, e3max
    type(SingleParticleOrbitIsospin), pointer :: o1, o2, o3
    integer :: i1, i2, i3, cnt
    integer :: j1, j2, j3, t12, p1, p2, p3

    this%t = t
    this%ch1 = ch1
    this%ch2 = ch2
    this%ch3 = ch3

    j1 = one%j(ch1)
    p1 = one%p(ch1)
    j2 = one%j(ch2)
    p2 = one%p(ch2)
    j3 = one%j(ch3)
    p3 = one%p(ch3)
    allocate(this%nst2n(0:sps%emax/2, 0:sps%emax/2, 0:sps%emax/2, 0:1))
    this%nst2n(:,:,:,:) = 0

    cnt = 0
    do i1 = 1, sps%norbs
      o1 => sps%orb(i1)
      if(j1 /= o1%j) cycle
      if(p1 /= (-1)**o1%l) cycle
      do i2 = 1, sps%norbs
        o2 => sps%orb(i2)
        if(j2 /= o2%j) cycle
        if(p2 /= (-1)**o2%l) cycle
        if(o1%e + o2%e > e2max) cycle
        do i3 = 1, sps%norbs
          o3 => sps%orb(i3)
          if(j3 /= o3%j) cycle
          if(p3 /= (-1)**o3%l) cycle
          if(o1%e + o3%e > e2max) cycle
          if(o2%e + o3%e > e2max) cycle
          if(o1%e + o2%e + o3%e > e3max) cycle
          do t12 = 0, 1
            if(triag(2*t12, 1, t)) cycle
            cnt = cnt + 1
          end do
        end do
      end do
    end do
    this%n_state = cnt
    allocate(this%n2abct(4,this%n_state))

    cnt = 0
    do i1 = 1, sps%norbs
      o1 => sps%orb(i1)
      if(j1 /= o1%j) cycle
      if(p1 /= (-1)**o1%l) cycle
      do i2 = 1, sps%norbs
        o2 => sps%orb(i2)
        if(j2 /= o2%j) cycle
        if(p2 /= (-1)**o2%l) cycle
        if(o1%e + o2%e > e2max) cycle
        do i3 = 1, sps%norbs
          o3 => sps%orb(i3)
          if(j3 /= o3%j) cycle
          if(p3 /= (-1)**o3%l) cycle
          if(o1%e + o3%e > e2max) cycle
          if(o2%e + o3%e > e2max) cycle
          if(o1%e + o2%e + o3%e > e3max) cycle
          do t12 = 0, 1
            if(triag(2*t12, 1, t)) cycle
            cnt = cnt + 1
            this%n2abct(:,cnt) = [i1,i2,i3,t12]
            this%nst2n(o1%n, o2%n, o3%n, t12) = cnt
          end do
        end do
      end do
    end do

  end subroutine InitThreeBodyMonChannel

  subroutine InitThreeBodyMonForce(this, sps, isps, e2max, e3max)
    use MyLibrary, only: triag, dcg
    class(ThreeBodyMonForce), intent(inout) :: this
    type(Orbits), intent(in), target :: sps
    type(OrbitsIsospin), intent(in), target :: isps
    integer, intent(in) :: e2max, e3max
    integer :: ch, n_state
    integer :: t1, t2, t3, z1, z2, z3

    this%sps => sps
    this%isps => isps
    this%emax = sps%emax
    this%e2max = e2max
    this%e3max = e3max

    call this%thr%init(isps, e2max, e3max)
    allocate(this%MatCh(this%thr%NChan))

    do ch = 1, this%thr%NChan
      n_state = this%thr%chan(ch)%n_state
      allocate(this%MatCh(ch)%v(n_state,n_state))
#ifdef single_precision_three_body_file
      this%MatCh(ch)%v(:,:) = 0.0
#else
      this%MatCh(ch)%v(:,:) = 0.d0
#endif
    end do

    allocate(this%CGs(0:3,-3:3, 0:3, -3:3, 0:3, -3:3))
    this%CGs(:,:,:,:,:,:) = 0.d0

    do t1 = 0,3
      do t2 = 0,3
        do t3 = 0,3
          do z1 = -t1, t1
            do z2 = -t2, t2
              do z3 = -t3, t3
                if(z1+z2 /= z3) cycle
                if(triag(t1,t2,t3)) cycle
                this%CGs(t1,z1,t2,z2,t3,z3) = dcg(t1,z1,t2,z2,t3,z3)
              end do
            end do
          end do
        end do
      end do
    end do
    this%zero = .false.
  end subroutine InitThreeBodyMonForce

  subroutine FinThreeBodyMonForce(this)
    class(ThreeBodyMonForce), intent(inout) :: this
    integer :: ch
    if(.not. allocated(this%MatCh)) return
    do ch = 1, this%thr%NChan
      deallocate(this%MatCh(ch)%v)
    end do
    deallocate(this%MatCh)
    deallocate(this%CGs)
    call this%thr%fin()
    this%zero = .true.
  end subroutine FinThreeBodyMonForce

  function GetThBMEMon_pn(this,i1,i2,i3,i4,i5,i6) result(r)
    use MyLibrary, only: dcg
    class(ThreeBodyMonForce), intent(in) :: this
    integer, intent(in) :: i1,i2,i3,i4,i5,i6
    type(Orbits), pointer :: sps
    type(OrbitsIsospin), pointer :: isps
    type(SingleParticleOrbit), pointer :: o1,o2,o3,o4,o5,o6
    integer :: z1, z2, z3, z4, z5, z6, T12, T45, T
    integer :: a, b, c, d, e, f, Z
    real(8) :: r
    !real(8) :: ti ! --- test

    !ti = omp_get_wtime() ! --- test (conflict with omp)
    sps => this%sps
    isps => this%isps
    r = 0.d0
    o1 => sps%GetOrbit(i1)
    o2 => sps%GetOrbit(i2)
    o3 => sps%GetOrbit(i3)
    o4 => sps%GetOrbit(i4)
    o5 => sps%GetOrbit(i5)
    o6 => sps%GetOrbit(i6)
    if(o1%j /= o4%j) return
    if(o2%j /= o5%j) return
    if(o3%j /= o6%j) return
    !if(o1%e + o2%e > this%thr%e2max) return
    !if(o4%e + o5%e > this%thr%e2max) return
    !if(o1%e + o3%e > this%thr%e2max) return
    !if(o4%e + o6%e > this%thr%e2max) return
    !if(o2%e + o3%e > this%thr%e2max) return
    !if(o5%e + o6%e > this%thr%e2max) return
    if(o1%e + o2%e + o3%e > this%thr%e3max) return
    if(o4%e + o5%e + o6%e > this%thr%e3max) return
    z1 = o1%z; z2 = o2%z; z3 = o3%z
    z4 = o4%z; z5 = o5%z; z6 = o6%z
    if(z1+z2+z3 /= z4+z5+z6) return
    Z = z1 + z2 + z3

    a = isps%nlj2idx( o1%n, o1%l, o1%j )
    b = isps%nlj2idx( o2%n, o2%l, o2%j )
    c = isps%nlj2idx( o3%n, o3%l, o3%j )
    d = isps%nlj2idx( o4%n, o4%l, o4%j )
    e = isps%nlj2idx( o5%n, o5%l, o5%j )
    f = isps%nlj2idx( o6%n, o6%l, o6%j )
    do T12 = 0, 1
      if(abs(z1+z2) > 2*T12) cycle
      do T45 = 0, 1
        if(abs(z4+z5) > 2*T45) cycle
        do T = max(abs(2*T12-1),abs(2*T45-1)), min(2*T12+1,2*T45+1), 2
          if(abs(Z) > T) cycle
          r = r + &
              & this%GetMonThBME(a,b,c,T12,d,e,f,T45,T) * &
              & this%CGs(1,z1,1,z2,2*T12,z1+z2) * this%CGs(2*T12,z1+z2,1,z3,T,Z) * &
              & this%CGs(1,z4,1,z5,2*T45,z4+z5) * this%CGs(2*T45,z4+z5,1,z6,T,Z)
        end do
      end do
    end do
    !call timer%add("GetThBME_pn", omp_get_wtime()-ti) ! --- test (conflict with omp)
  end function GetThBMEMon_pn

  function GetThBMEMon_Isospin(this,i1,i2,i3,T12,i4,i5,i6,T45,T) result(r)
    class(ThreeBodyMonForce), intent(in) :: this
    integer, intent(in) :: i1,i2,i3,i4,i5,i6
    integer, intent(in) :: T12,T45,T
    type(OrbitsIsospin), pointer :: isps
    type(SingleParticleOrbitIsospin), pointer :: o1,o2,o3,o4,o5,o6
    integer :: ch, bra, ket, ch_t, ch_mon
    integer :: ch1, ch2, ch3
    integer :: P123, P456
    real(8) :: r
    !real(8) :: ti ! --- test

    r = 0.d0
    !ti = omp_get_wtime() ! --- test (conflict with omp)
    isps => this%isps
    o1 => isps%GetOrbit(i1)
    o2 => isps%GetOrbit(i2)
    o3 => isps%GetOrbit(i3)
    o4 => isps%GetOrbit(i4)
    o5 => isps%GetOrbit(i5)
    o6 => isps%GetOrbit(i6)
    if(o1%j /= o4%j) return
    if(o2%j /= o5%j) return
    if(o3%j /= o6%j) return
    if(o1%l /= o4%l) return
    if(o2%l /= o5%l) return
    if(o3%l /= o6%l) return

    P123 = (-1) ** (o1%l+o2%l+o3%l)
    P456 = (-1) ** (o4%l+o5%l+o6%l)
    if(P123 * P456 /= 1) then
      write(*,*) 'Warning: in GetThBMEIso_scalar: P'
      return
    end if

    ch1 = this%thr%one%jp2ch(o1%j, (-1)**o1%l)
    ch2 = this%thr%one%jp2ch(o2%j, (-1)**o2%l)
    ch3 = this%thr%one%jp2ch(o3%j, (-1)**o3%l)

    ch_t = this%thr%t2ch(T)
    ch_mon = this%thr%mon2ch(ch1,ch2,ch3)
    if(ch_t * ch_mon == 0) return
    ch = this%thr%tmon2ch(ch_t,ch_mon)
    if(ch == 0) return

    bra = this%thr%chan(ch)%nst2n(o1%n,o2%n,o3%n,t12)
    ket = this%thr%chan(ch)%nst2n(o4%n,o5%n,o6%n,t45)
    if(bra * ket == 0) return
    r = dble(this%MatCh(ch)%v(bra,ket))
    !call timer%add("GetThBME_isospin", omp_get_wtime()-ti) ! --- test (conflict with omp)
  end function GetThBMEMon_Isospin

  !
  !
  ! reading three-body scalar
  !
  !
  subroutine Set3BodyMonopoleFile(this, file_3n)
    use ClassSys, only: sys
    class(Read3BodyMonopole), intent(inout) :: this
    character(*), intent(in), optional :: file_3n
    type(sys) :: s
    logical :: ex

    ! -- three-body file
    this%file_3n = 'none'
    if(present(file_3n)) this%file_3n = file_3n
    select case(this%file_3n)
    case("NONE", "none", "None")
    case default
      ex = s%isfile(this%file_3n, "SetReadFiles: three-body file")
    end select
  end subroutine Set3BodyMonopoleFile

  subroutine Set3BodyMonopoleFileBoundaries(this, emax, e2max, e3max, lmax)
    class(Read3BodyMonopole), intent(inout) :: this
    integer, intent(in) :: emax, e2max, e3max, lmax

    this%emax3 = emax
    this%e2max3= e2max
    this%e3max3= e3max
    this%lmax3 = lmax

  end subroutine Set3BodyMonopoleFileBoundaries

  subroutine ReadThreeBodyMonopole(this, V)
    class(Read3BodyMonopole), intent(in) :: this
    type(ThreeBodyMonForce), intent(inout) :: V
    real(8) :: ti

    ti = omp_get_wtime()

    select case(this%file_3n)
    case('None', 'NONE', 'none')
      write(*,*) "No three-body matrix element."
      return
    case default

      call this%ReadFile(V)
    end select

    call timer%Add('Read from file', omp_get_wtime()-ti)

  end subroutine ReadThreeBodyMonopole

  subroutine ReadFile(this,thr)
    use ClassSys, only: sys
    class(Read3BodyMonopole), intent(in) :: this
    type(ThreeBodyMonForce), intent(inout) :: thr
    type(sys) :: s

    if(s%find(this%file_3n,'.txt')) then
      call this%read_scalar_me3j_ascii_txt(thr)
      return
    end if

    if(s%find(this%file_3n,'.me3j.gz')) then
      call this%read_scalar_me3j_gzip(thr)
      return
    end if

    if(s%find(this%file_3n,'.me3j')) then
      call this%read_scalar_me3j_ascii(thr)
      return
    end if

    !if(s%find(this%file_3n,'stream.bin')) then
    !  call this%read_scalar_me3j_binary_stream(thr)
    !  return
    !end if

    !if(s%find(this%file_3n,'.bin') .and. s%find(this%file_3n,'_comp')) then
    !  call this%read_scalar_me3j_binary_comp(thr)
    !  return
    !end if

    !if(s%find(this%file_3n,'.bin')) then
    !  call this%read_scalar_me3j_binary(thr)
    !  return
    !end if

    !call this%read_scalar_me3j_binary_stream(thr)
  end subroutine ReadFile

  subroutine read_scalar_me3j_ascii_txt(this,thr)
    class(Read3BodyMonopole), intent(in) :: this
    type(ThreeBodyMonForce), intent(inout) :: thr
    type(OrbitsIsospin) :: spsf
#ifdef single_precision_three_body_file
    real(4), allocatable :: v(:)
#else
    real(8), allocatable :: v(:)
#endif
    integer(8) :: nelm, n
    integer :: runit = 22, io

    write(*,'(a)') "Reading three-body scalar line-by-line from human-readable file"

    call spsf%init(this%emax3, this%lmax3)
    nelm = count_scalar_3bme(spsf, this%e2max3, this%e3max3)
    allocate(v(nelm))
    open(runit, file=this%file_3n, action='read', iostat=io)
    if(io /= 0) then
      write(*,'(2a)') "File opening error: ", trim(this%file_3n)
      return
    end if
    do n = 1, nelm
      read(runit,*) v(n)
    end do
    close(runit)

    call store_scalar_3bme(thr,v,spsf,this%e2max3,this%e3max3)

    deallocate(v)
    call spsf%fin()
  end subroutine read_scalar_me3j_ascii_txt

  subroutine read_scalar_me3j_ascii(this,thr)
    class(Read3BodyMonopole), intent(in) :: this
    type(ThreeBodyMonForce), intent(inout) :: thr
    type(OrbitsIsospin) :: spsf
#ifdef single_precision_three_body_file
    real(4), allocatable :: v(:)
#else
    real(8), allocatable :: v(:)
#endif
    integer(8) :: nelm, n
    integer :: runit = 22

    write(*,'(a)') "Reading three-body scalar line-by-line from human-readable file"

    call spsf%init(this%emax3, this%lmax3)
    nelm = count_scalar_3bme(spsf, this%e2max3, this%e3max3)
    allocate(v(nelm))

    open(runit, file=this%file_3n, action='read')
    read(runit,*)
    do n = 1, nelm/10
      read(runit,*) v((n-1)*10+1 : n*10)
    end do

    ! basically, this is not needed
    if(nelm - (nelm/10) * 10 > 0) then
      read(runit,*) v((nelm/10)*10+1 : nelm)
    end if
    close(runit)

    call store_scalar_3bme(thr,v,spsf,this%e2max3,this%e3max3)

    deallocate(v)
    call spsf%fin()
  end subroutine read_scalar_me3j_ascii

  subroutine read_scalar_me3j_gzip(this,thr)
    use, intrinsic :: iso_c_binding
    use MyLibrary, only: gzip_open, gzip_readline, gzip_close
    class(Read3BodyMonopole), intent(in) :: this
    type(ThreeBodyMonForce), intent(inout) :: thr
    type(OrbitsIsospin) :: spsf
#ifdef single_precision_three_body_file
    real(4), allocatable :: v(:)
#else
    real(8), allocatable :: v(:)
#endif
    integer(8) :: nelm, n
    character(256) :: header, buffer
    type(c_ptr) :: fp, err

    write(*,'(a)') "Reading three-body scalar line-by-line from gzip file"

    call spsf%init(this%emax3, this%lmax3)
    nelm = count_scalar_3bme(spsf, this%e2max3, this%e3max3)
    allocate(v(nelm))

    fp = gzip_open(this%file_3n, "rt")
    err = gzip_readline(fp, header, len(header))
    do n = 1, nelm/10
      err = gzip_readline(fp, buffer, len(buffer))
      read(buffer,*) v((n-1)*10+1 : n*10)
    end do

    ! basically, this is not needed
    if(nelm - (nelm/10) * 10 > 0) then
      err = gzip_readline(fp, buffer, len(buffer))
      read(buffer,*) v((nelm/10)*10+1 : nelm)
    end if
    err = gzip_close(fp)

    call store_scalar_3bme(thr,v,spsf,this%e2max3,this%e3max3)

    deallocate(v)
    call spsf%fin()
  end subroutine read_scalar_me3j_gzip

  function count_scalar_3bme(spsf, e2max, e3max) result(r)
    type(OrbitsIsospin), intent(in) :: spsf
    integer, intent(in) :: e2max, e3max
    integer(8) :: r
    integer :: i1, l1, j1, e1
    integer :: i2, l2, j2, e2
    integer :: i3, l3, j3, e3
    integer :: i4, l4, j4, e4
    integer :: i5, l5, j5, e5, i5max
    integer :: i6, l6, j6, e6, i6max
    integer :: T12, T45, T, P123, P456

    r = 0
    do i1 = 1, spsf%norbs
      l1 = spsf%orb(i1)%l
      j1 = spsf%orb(i1)%j
      e1 = spsf%orb(i1)%e
      do i2 = 1, spsf%norbs
        l2 = spsf%orb(i2)%l
        j2 = spsf%orb(i2)%j
        e2 = spsf%orb(i2)%e
        if(e1 + e2 > e2max) cycle
        do i3 = 1, spsf%norbs
          l3 = spsf%orb(i3)%l
          j3 = spsf%orb(i3)%j
          e3 = spsf%orb(i3)%e
          if(e1 + e3 > e2max) cycle
          if(e2 + e3 > e2max) cycle
          if(e1 + e2 + e3 > e3max) cycle

          P123 = (-1) ** (l1+l2+l3)

          do i4 = 1, spsf%norbs
            l4 = spsf%orb(i4)%l
            j4 = spsf%orb(i4)%j
            e4 = spsf%orb(i4)%e

            if(j1 /= j4) cycle
            if(l1 /= l4) cycle

            i5max = i4
            if(i1 == i4) i5max = i2

            do i5 = 1, spsf%norbs
              l5 = spsf%orb(i5)%l
              j5 = spsf%orb(i5)%j
              e5 = spsf%orb(i5)%e
              if(j2 /= j5) cycle
              if(l2 /= l5) cycle
              if(e4 + e5 > e2max) cycle

              i6max = i5
              if(i1 == i4 .and. i2 == i5) i6max = i3

              do i6 = 1, spsf%norbs
                l6 = spsf%orb(i6)%l
                j6 = spsf%orb(i6)%j
                e6 = spsf%orb(i6)%e
                if(j3 /= j6) cycle
                if(l3 /= l6) cycle
                if(e4 + e6 > e2max) cycle
                if(e5 + e6 > e2max) cycle
                if(e4 + e5 + e6 > e3max) cycle

                P456 = (-1) ** (l4+l5+l6)

                if(P123 /= P456) cycle

                do T12 = 0, 1
                  do T45 = 0, 1
                    do T = max(abs(2*T12-1),abs(2*T45-1)),&
                          &min(   (2*T12+1),   (2*T45+1)), 2

                      r = r + 1

                    end do

                  end do
                end do
              end do
            end do
          end do


        end do
      end do
    end do
  end function count_scalar_3bme

  subroutine store_scalar_3bme(thr,v,spsf,e2max,e3max)
    type(ThreeBodyMonForce), intent(inout), target :: thr
    type(OrbitsIsospin), intent(in) :: spsf
#ifdef single_precision_three_body_file
    real(4), allocatable :: v(:)
#else
    real(8), allocatable :: v(:)
#endif
    integer, intent(in) :: e2max, e3max
    type(ThreeBodyMonSpace), pointer :: ms
    type(OrbitsIsospin), pointer :: sps
    integer(8) :: cnt
    integer :: i1, n1, l1, j1, e1, ch1
    integer :: i2, n2, l2, j2, e2, ch2
    integer :: i3, n3, l3, j3, e3, ch3
    integer :: i4, n4, l4, j4, e4
    integer :: i5, n5, l5, j5, e5, i5max
    integer :: i6, n6, l6, j6, e6, i6max
    integer :: T12, T45, T, P123, P456
    integer :: ch, ch_t, ch_mon, bra, ket

    ms => thr%thr
    sps => ms%sps

    cnt = 0
    do i1 = 1, spsf%norbs
      n1 = spsf%orb(i1)%n
      l1 = spsf%orb(i1)%l
      j1 = spsf%orb(i1)%j
      e1 = spsf%orb(i1)%e
      do i2 = 1, spsf%norbs
        n2 = spsf%orb(i2)%n
        l2 = spsf%orb(i2)%l
        j2 = spsf%orb(i2)%j
        e2 = spsf%orb(i2)%e
        if(e1 + e2 > e2max) cycle
        do i3 = 1, spsf%norbs
          n3 = spsf%orb(i3)%n
          l3 = spsf%orb(i3)%l
          j3 = spsf%orb(i3)%j
          e3 = spsf%orb(i3)%e
          if(e1 + e3 > e2max) cycle
          if(e2 + e3 > e2max) cycle
          if(e1 + e2 + e3 > e3max) cycle

          P123 = (-1) ** (l1+l2+l3)

          do i4 = 1, spsf%norbs
            n4 = spsf%orb(i4)%n
            l4 = spsf%orb(i4)%l
            j4 = spsf%orb(i4)%j
            e4 = spsf%orb(i4)%e
            if(j1 /= j4) cycle
            if(l1 /= l4) cycle

            i5max = i4
            if(i1 == i4) i5max = i2

            do i5 = 1, spsf%norbs
              n5 = spsf%orb(i5)%n
              l5 = spsf%orb(i5)%l
              j5 = spsf%orb(i5)%j
              e5 = spsf%orb(i5)%e
              if(j2 /= j5) cycle
              if(l2 /= l5) cycle
              if(e4 + e5 > e2max) cycle

              i6max = i5
              if(i1 == i4 .and. i2 == i5) i6max = i3

              do i6 = 1, spsf%norbs
                n6 = spsf%orb(i6)%n
                l6 = spsf%orb(i6)%l
                j6 = spsf%orb(i6)%j
                e6 = spsf%orb(i6)%e
                if(j3 /= j6) cycle
                if(l3 /= l6) cycle
                if(e4 + e6 > e2max) cycle
                if(e5 + e6 > e2max) cycle
                if(e4 + e5 + e6 > e3max) cycle

                P456 = (-1) ** (l4+l5+l6)

                if(P123 /= P456) cycle

                do T12 = 0, 1
                  do T45 = 0, 1
                    do T = max(abs(2*T12-1),abs(2*T45-1)),&
                          &min(   (2*T12+1),   (2*T45+1)), 2
                      cnt = cnt + 1

                      if(e1 > ms%emax) cycle
                      if(e2 > ms%emax) cycle
                      if(e3 > ms%emax) cycle

                      if(e4 > ms%emax) cycle
                      if(e5 > ms%emax) cycle
                      if(e6 > ms%emax) cycle

                      if(e1 + e2 > ms%e2max) cycle
                      if(e2 + e3 > ms%e2max) cycle
                      if(e3 + e1 > ms%e2max) cycle

                      if(e4 + e5 > ms%e2max) cycle
                      if(e5 + e6 > ms%e2max) cycle
                      if(e6 + e4 > ms%e2max) cycle

                      if(e1 + e2 + e3 > ms%e3max) cycle
                      if(e4 + e5 + e6 > ms%e3max) cycle

                      ch1 = ms%one%jp2ch(j1,(-1)**l1)
                      ch2 = ms%one%jp2ch(j2,(-1)**l2)
                      ch3 = ms%one%jp2ch(j3,(-1)**l3)
                      ch_t = ms%t2ch(T)
                      ch_mon = ms%mon2ch(ch1,ch2,ch3)
                      if(ch_t * ch_mon == 0) cycle
                      ch = ms%tmon2ch(ch_t,ch_mon)
                      if(ch == 0) cycle

                      bra = ms%chan(ch)%nst2n(n1,n2,n3,t12)
                      ket = ms%chan(ch)%nst2n(n4,n5,n6,t45)
                      if(bra * ket == 0) cycle
                      thr%MatCh(ch)%v(bra,ket) = v(cnt)
                      thr%MatCh(ch)%v(ket,bra) = v(cnt)
                    end do

                  end do
                end do


              end do
            end do
          end do

        end do
      end do
    end do
  end subroutine store_scalar_3bme


end module ThreeBodyMonInteraction
