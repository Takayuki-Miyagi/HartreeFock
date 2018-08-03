module ModelSpace
  use common_library
  use InputParameters, only: parameters
  use MPIFunction, only: myrank
  implicit none
  private :: Aop3, A3drct, A3exc1, A3exc2, SingleParticleOrbit, iter3, &
      & InitSPOIsospin, FinSPOIsospin, InitSPO, FinSPO, GetReferenceState, GetNOCoef, &
      & InitMSpace, FinMSpace, InitOneBodyChannel, FinOneBodyChannel, InitOneBodySpace, &
      & FinOneBodySpace, InitTwoBodySpace, FinTwoBodySpace, InitTwoBodyChannel, FinTwoBodyChannel, &
      & InitThreeBodySpace, FinThreeBodySpace, InitThreeBodyChannel, FinThreeBodyChannel, &
      & InitAdditionalQN, FinAdditionalQN, CntDim, GenIter
  public :: spo_isospin, spo_pn, SpinParityTz, OneBodyChannel, OneBodySpace, TwoBodyChannel, &
      & TwoBodySpace, ThreeBodyChannel, ThreeBodySpace, AdditionalQN, MSpace

  type :: SingleParticleOrbit
    integer :: n
    integer, allocatable :: nn(:), ll(:), jj(:), nshell(:)
  end type SingleParticleOrbit

  ! single particle state (isospin symmetry)
  type, extends(SingleParticleOrbit) :: spo_isospin
    integer, allocatable :: spo2n(:,:,:)
  contains
    procedure :: init => InitSPOIsospin
    procedure :: fin => FinSPOIsospin
  end type spo_isospin

  ! single particle state (proton neutron)
  type, extends(SingleParticleOrbit) :: spo_pn
    integer, allocatable :: itz(:)
    integer, allocatable :: spo2n(:,:,:,:)
    integer, allocatable :: h_orbits(:)
    integer, allocatable :: p_orbits(:)
    real(8), allocatable :: nocoef(:)
  contains
    procedure :: init => InitSPO
    procedure :: fin => FinSPO
    procedure :: GetReferenceState
    procedure :: GetNOCoef
  end type spo_pn

  type :: SpinParityTz
    integer :: n
    integer, allocatable :: j(:), p(:), tz(:)
    integer, allocatable :: jptz2n(:,:,:)
    integer, allocatable :: ndim(:)
  end type SpinParityTz

  type :: OneBodyChannel
    integer :: n
    integer, allocatable :: n2label(:), label2n(:)
  contains
    procedure :: init => InitOneBodyChannel
    procedure :: fin => FinOneBodyChannel
  end type OneBodyChannel

  type :: TwoBodyChannel
    integer :: n
    integer, allocatable :: n2label1(:), n2label2(:)
    integer, allocatable :: labels2n(:,:)
    integer, allocatable :: iphase(:,:)
  contains
    procedure :: init => InitTwoBodyChannel
    procedure :: fin => FinTwoBodyChannel
  end type TwoBodyChannel

  type :: AdditionalQN
    integer :: nphys, n
    integer, allocatable :: labels2n(:)
    integer, allocatable :: n2label1(:)
    integer, allocatable :: n2label2(:)
    integer, allocatable :: n2label3(:)
    integer, allocatable :: n2label4(:)
    real(8), allocatable :: cfp(:,:)
  contains
    procedure :: init => InitAdditionalQN
    procedure :: fin => finAdditionalQN
  end type AdditionalQN

  type :: ThreeBodyChannel
    integer :: n, nsub
    type(AdditionalQN), allocatable :: idx(:)
    integer, allocatable :: n2label1(:), n2label2(:), n2label3(:), n2label4(:)
    integer, allocatable :: labels2nsub(:,:,:)
  contains
    procedure :: init => InitThreeBodyChannel
    procedure :: fin => FinThreeBodyChannel
  end type ThreeBodyChannel

  type, extends(SpinParityTz) :: OneBodySpace
    type(OneBodyChannel), allocatable :: jptz(:)
    integer, allocatable :: label2jptz(:)
  contains
    procedure :: init => InitOneBodySpace
    procedure :: fin => FinOneBodySpace
  end type OneBodySpace

  type, extends(SpinParityTz) :: TwoBodySpace
    type(TwoBodyChannel), allocatable :: jptz(:)
  contains
    procedure :: init => InitTwoBodySpace
    procedure :: fin => FinTwoBodySpace
  end type TwoBodySpace

  type, extends(SpinParityTz) :: ThreeBodySpace
    type(ThreeBodyChannel), allocatable :: jptz(:)
  contains
    procedure :: init => InitThreeBodySpace
    procedure :: fin => FinThreeBodySpace
  end type ThreeBodySpace

  type :: MSpace
    type(OneBodySpace) :: one
    type(TwoBodySpace) :: two
    type(ThreeBodySpace) :: thr
  contains
    procedure :: init => InitMSpace
    procedure :: fin => FinMSpace
  end type MSpace

  type :: iter3
    integer, allocatable :: ii1(:), ii2(:), ii3(:)
  contains
    procedure :: GenIter
  end type iter3

contains
  subroutine InitSPOIsospin(this, params)
    class(spo_isospin), intent(inout) :: this
    type(parameters), intent(in) :: params
    integer :: emax, n, nl, nn, ll, is, j
    emax = params%emax_3nf
    allocate(this%spo2n(0:emax/2, 0:emax, 1:2*emax+1))
    this%spo2n = 0
    this%n = (emax + 1) * (emax + 2) / 2
    n = this%n
    allocate(this%nn(n))
    allocate(this%ll(n))
    allocate(this%jj(n))
    allocate(this%nshell(n))
    this%nn = 0; this%ll = 0; this%jj = 0
    this%nshell = 0
#ifdef debug
    write(*,'(a)') '##################################################'
    write(*,'(a)') ' Single-Particle State (Isospin formalism)'
    write(*,'(a)') '##################################################'
    write(*,'(a)') '      i,     n,     l,     j, shell'
#endif
    n = 0
    do nl = 0, emax

      do ll = 0, nl
        if(mod(nl - ll, 2) == 1) cycle
        nn = (nl - ll) / 2
        do is = -1, 1, 2
          j = 2 * ll + is
          if(j < 0) cycle
          n = n + 1
          this%nn(n) = nn
          this%ll(n) = ll
          this%jj(n) = j
          this%nshell(n) = 2 * nn + ll
          this%spo2n(nn, ll, j) = n
#ifdef debug
          write(*,'(5i7)') n, nn, ll, j, 2 * nn + ll
#endif
        end do
      end do
    end do
#ifdef debug
    write(*,*)
#endif
  end subroutine InitSPOIsospin

  subroutine FinSPOIsospin(this)
    class(spo_isospin), intent(inout) :: this
    deallocate(this%nn)
    deallocate(this%ll)
    deallocate(this%jj)
    deallocate(this%nshell)
    deallocate(this%spo2n)
  end subroutine FinSPOIsospin

  subroutine InitSPO(this, params)
    class(spo_pn), intent(inout) :: this
    type(parameters), intent(in) :: params
    integer :: emax, n, nl, nn, ll, is, j, itz
    emax = params%emax
    allocate(this%spo2n(0:emax/2, 0:emax, 1:2*emax+1, -1:1))
    this%spo2n = 0
    this%n = (emax + 1) * (emax + 2)
    n = this%n
    allocate(this%nn(n))
    allocate(this%ll(n))
    allocate(this%jj(n))
    allocate(this%itz(n))
    allocate(this%nshell(n))
    allocate(this%h_orbits(n))
    allocate(this%p_orbits(n))
    allocate(this%nocoef(n))
    this%nn = 0; this%ll = 0; this%jj = 0
    this%itz = 0; this%nshell = 0
    this%h_orbits = 0
    this%p_orbits = 0
    this%nocoef = 0.d0
    n = 0
#ifdef debug
    write(*,'(a)') '##################################################'
    write(*,'(a)') ' Single-Particle State (PN formalism)'
    write(*,'(a)') '##################################################'
    write(*,'(a)') '      i,     n,     l,     j,    iz, shell'
#endif
    do nl = 0, emax
      do ll = 0, nl
        if(mod(nl - ll, 2) == 1) cycle
        nn = (nl - ll) / 2
        do is = -1, 1, 2
          j = 2 * ll + is
          if(j < 0) cycle
          do itz = -1, 1, 2
            n = n + 1
            this%nn(n) = nn
            this%ll(n) = ll
            this%jj(n) = j
            this%itz(n) = itz
            this%nshell(n) = 2 * nn + ll
            this%spo2n(nn, ll, j, itz) = n
#ifdef debug
            write(*,'(6i7)') n, nn, ll, j, itz, 2 * nn + ll
#endif
          end do
        end do
      end do
    end do
#ifdef debug
    write(*,*)
#endif
    call this%GetReferenceState(params)
    call this%GetNOCoef(params)
  end subroutine InitSPO

  subroutine FinSPO(this)
    class(spo_pn), intent(inout) :: this
    deallocate(this%nn)
    deallocate(this%ll)
    deallocate(this%jj)
    deallocate(this%itz)
    deallocate(this%nshell)
    deallocate(this%spo2n)
    deallocate(this%h_orbits)
    deallocate(this%p_orbits)
    deallocate(this%nocoef)
  end subroutine FinSPO

  subroutine GetReferenceState(this, params)
    class(spo_pn), intent(inout) :: this
    type(parameters), intent(in) :: params
    integer :: num, i
    integer :: n, l, j, itz
    integer :: iunit = 40
    open(iunit, file = params%reference, status = 'old')
    read(iunit,*) num
    if(num < 1) return
    call skip_comment(iunit, '!')
    do i = 1, num
      read(iunit,*) n, l, j, itz
      this%h_orbits(this%spo2n(n,l,j,itz)) = 1
    end do
    if(myrank == 0) then
      write(*,'(2x, a)') 'Core Orbits: '
      do i = 1, this%n
        if(this%h_orbits(i) /= 0) then
          write(*,'(2x,a,i3,a,i3,a,i3,a,i3)') 'n = ', this%nn(i), &
              & ',  l = ', this%ll(i), ',  j = ', this%jj(i), &
              & ',  itz = ', this%itz(i)
        end if
      end do
      write(*,*)
    end if
  end subroutine GetReferenceState

  subroutine GetNOCoef(this, params)
    class(spo_pn), intent(inout) :: this
    type(parameters), intent(in) :: params
    integer :: num, i
    integer :: n, l, j, itz
    integer :: iunit = 40
    real(8) :: f
    open(iunit, file = params%nocoef, status = 'old')
    read(iunit,*) num
    if(num < 1) return
    call skip_comment(iunit, '!')
    do i = 1, num
      read(iunit,*) n, l, j, itz, f
      this%nocoef(this%spo2n(n,l,j,itz)) = dble(f)
    end do
    if(myrank == 0) then
      write(*,'(2x, a)') 'Coefficients for Normal Ordering: '
      do i = 1, this%n
        if(this%nocoef(i) /= 0) then
          write(*,'(2x,a,i3,a,i3,a,i3,a,i3,a,f7.3)') 'n = ', this%nn(i), &
              & ',  l = ', this%ll(i), ',  j = ', this%jj(i), &
              & ',  itz = ', this%itz(i), ',  Coef = ', this%nocoef(i)
        end if
      end do
      write(*,*)
    end if
  end subroutine GetNOCoef

 ! Model Space
  subroutine InitMSpace(this, sps, params)
    class(MSpace), intent(inout) :: this
    type(spo_pn), intent(in) :: sps
    type(parameters), intent(in) :: params
    call this%one%init(sps, params)
    call this%two%init(sps, params)
    !call this%thr%init(sps, params)
  end subroutine InitMSpace

  subroutine FinMSpace(this)
    class(MSpace), intent(inout) :: this
    call this%one%fin()
    call this%two%fin()
    !call this%thr%fin()
  end subroutine FinMSpace

  ! Model Space One Body
  subroutine InitOneBodySpace(this, sps, params)
    class(OneBodySpace), intent(inout) :: this
    type(spo_pn), intent(in) :: sps
    type(parameters), intent(in) :: params
    integer :: emax, ich, n, i1
    integer :: j, p, itz, loop
    emax = params%emax
    do loop = 1, 2
      ich = 0
      do j = 1, 2 * emax + 1, 2
        do p = 1, -1, -2
          do itz = -1, 1, 2
            n = 0
            do i1 = 1, sps%n
              if(sps%jj(i1) /= j) cycle
              if((-1) ** sps%ll(i1) /= p) cycle
              if(sps%itz(i1) /= itz) cycle
              n = n + 1
            end do
            if(n /= 0) ich = ich + 1
            if(loop == 2 .and. n /= 0) then
              this%j(ich) = j
              this%p(ich) = p
              this%tz(ich) = itz
              this%ndim(ich) = n
              this%jptz(ich)%n = n
              this%jptz2n(j, p, itz) = ich
            end if
          end do
        end do
      end do
      if(loop == 1) then
        this%n = ich
        allocate(this%jptz(ich))
        allocate(this%j(ich))
        allocate(this%p(ich))
        allocate(this%tz(ich))
        allocate(this%ndim(ich))
        allocate(this%jptz2n(2*emax+1,-1:1,-1:1))
        this%j(:) = 0
        this%p(:) = 0
        this%tz(:) = 0
        this%ndim(:) = 0
        this%jptz2n(:,:,:) = 0
      end if
    end do

    do ich = 1, this%n
      call this%jptz(ich)%init(this%j(ich), this%p(ich), this%tz(ich), sps)
#ifdef debug
      write(*,'(a, i3, a, i3, a, i3, a, i3)') 'J = ', this%j(ich), '/2,  P = ', &
          & this%p(ich), ',  Tz = ', this%tz(ich), '/2,  n = ', this%jptz(ich)%n
#endif
    end do

    allocate(this%label2jptz(sps%n))
    do i1 = 1, sps%n
      j = sps%jj(i1)
      p = (-1) ** sps%ll(i1)
      itz = sps%itz(i1)
      ich = this%jptz2n(j, p, itz)
      this%label2jptz(i1) = ich
    end do

  end subroutine InitOneBodySpace

  subroutine FinOneBodySpace(this)
    class(OneBodySpace), intent(inout) :: this
    integer :: ich
    do ich = 1, this%n
      call this%jptz(ich)%fin()
    end do
    deallocate(this%jptz)
    deallocate(this%j)
    deallocate(this%p)
    deallocate(this%tz)
    deallocate(this%jptz2n)
    deallocate(this%label2jptz)
  end subroutine FinOneBodySpace

  subroutine InitOneBodyChannel(this, j, p, tz, sps)
    class(OneBodyChannel), intent(inout) :: this
    type(spo_pn), intent(in) :: sps
    integer, intent(in) :: j, p, tz
    integer :: n, m, i1
    n = this%n
    m = sps%n
    allocate(this%n2label(n))
    allocate(this%label2n(m))
    this%n2label(:) = 0
    this%label2n(:) = 0
    n = 0
    do i1 = 1, m
      if(sps%jj(i1) /= j) cycle
      if((-1) ** sps%ll(i1) /= p) cycle
      if(sps%itz(i1) /= tz) cycle
      n = n + 1
      this%n2label(n) = i1
      this%label2n(i1) = n
    end do
  end subroutine InitOneBodyChannel

  subroutine FinOneBodyChannel(this)
    class(OneBodyChannel), intent(inout) :: this
    deallocate(this%n2label)
    deallocate(this%label2n)
  end subroutine FinOneBodyChannel

  ! Model Space Two Body
  subroutine InitTwoBodySpace(this, sps, params)
    class(TwoBodySpace), intent(inout) :: this
    type(spo_pn), intent(in) :: sps
    type(parameters), intent(in) :: params
    integer :: e2max, ich, n, i1, i2
    integer :: j1, l1, itz1
    integer :: j2, l2, itz2
    integer :: j, p, itz, loop
    e2max = params%e2max
    do loop = 1, 2
      ich = 0
      do j = 0, e2max + 1
        do p = 1, -1, -2
          do itz = -1, 1
            n = 0
            do i1 = 1, sps%n
              j1 = sps%jj(i1)
              l1 = sps%ll(i1)
              itz1 = sps%itz(i1)
              do i2 = 1, i1
                j2 = sps%jj(i2)
                l2 = sps%ll(i2)
                itz2 = sps%itz(i2)
                if(sps%nshell(i1) + sps%nshell(i2) > e2max) cycle
                if(triag(j1, j2, 2*j)) cycle
                if((-1) ** (l1 + l2) /= p) cycle
                if(itz1 + itz2 /= 2 * itz) cycle
                if(i1 == i2 .and. mod(j, 2) == 1) cycle
                n = n + 1
              end do
            end do
            if(n /= 0) ich = ich + 1
            if(loop == 2 .and. n /= 0) then
              this%j(ich) = j
              this%p(ich) = p
              this%tz(ich) = itz
              this%jptz(ich)%n = n
              this%ndim(ich) = n
              this%jptz2n(j, p, itz) = ich
            end if
          end do
        end do
      end do
      if(loop == 1) then
        this%n = ich
        allocate(this%jptz(ich))
        allocate(this%j(ich))
        allocate(this%p(ich))
        allocate(this%tz(ich))
        allocate(this%ndim(ich))
        allocate(this%jptz2n(0:e2max+1,-1:1,-1:1))
        this%j(:) = 0
        this%p(:) = 0
        this%tz(:) = 0
        this%ndim(:) = 0
        this%jptz2n(:,:,:) = 0
      end if
    end do

    do ich = 1, this%n
      call this%jptz(ich)%Init(params, sps, this%j(ich), this%p(ich), this%tz(ich))
#ifdef debug
      write(*,'(a, i3, a, i3, a, i3, a, i6)') 'J = ', this%j(ich), ',  P = ', &
          & this%p(ich), ',  Tz = ', this%tz(ich), ',  n = ', this%jptz(ich)%n
#endif
    end do
#ifdef debug
    write(*,*)
#endif

  end subroutine InitTwoBodySpace

  subroutine FinTwoBodySpace(this)
    class(TwoBodySpace), intent(inout) :: this
    integer :: ich
    do ich = 1, this%n
      call this%jptz(ich)%fin()
    end do
    deallocate(this%jptz)
    deallocate(this%j)
    deallocate(this%p)
    deallocate(this%tz)
    deallocate(this%jptz2n)
  end subroutine FinTwoBodySpace

  subroutine InitTwoBodyChannel(this, params, sps, j, p, tz)
    class(TwoBodyChannel), intent(inout) :: this
    type(parameters), intent(in) :: params
    type(spo_pn), intent(in) :: sps
    integer, intent(in) :: j, p, tz
    integer :: n, m, i1, i2
    integer :: j1, l1, tz1, j2, l2, tz2
    n = this%n
    m = sps%n
    allocate(this%n2label1(n))
    allocate(this%n2label2(n))
    allocate(this%labels2n(m,m))
    allocate(this%iphase(m,m))
    this%n2label1 = 0
    this%n2label2 = 0
    this%labels2n = 0
    this%iphase = 0
    n = 0
    do i1 = 1, m
      j1 = sps%jj(i1)
      l1 = sps%ll(i1)
      tz1 = sps%itz(i1)
      do i2 = 1, i1
        j2 = sps%jj(i2)
        l2 = sps%ll(i2)
        tz2 = sps%itz(i2)
        if(sps%nshell(i1) + sps%nshell(i2) > params%e2max) cycle
        if(triag(j1, j2, 2*j)) cycle
        if((-1) ** (l1 + l2) /= p) cycle
        if(tz1 + tz2 /= 2 * tz) cycle
        if(i1 == i2 .and. mod(j, 2) == 1) cycle
        n = n + 1
        this%n2label1(n) = i1
        this%n2label2(n) = i2
        this%labels2n(i1,i2) = n
        this%labels2n(i2,i1) = n
        this%iphase(i1,i2) = 1
        this%iphase(i2,i1) = -(-1) ** ((j1 + j2) / 2 - j)
      end do
    end do
  end subroutine InitTwoBodyChannel

  subroutine FinTwoBodyChannel(this)
    class(TwoBodyChannel), intent(inout) :: this
    deallocate(this%n2label1)
    deallocate(this%n2label2)
    deallocate(this%labels2n)
    deallocate(this%iphase)
  end subroutine FinTwoBodyChannel

  ! Model Space Three Body -- not used in Hartree Fock calculations --
  subroutine InitThreeBodySpace(this, sps, params)
    class(ThreeBodySpace), intent(inout) :: this
    type(spo_pn), intent(in) :: sps
    type(parameters), intent(in) :: params
    integer :: ich, n
    integer :: i1, j1, l1, itz1
    integer :: i2, j2, l2, itz2
    integer :: i3, j3, l3, itz3
    integer :: j, p, itz, loop, ni, nj
    integer, allocatable :: labels2idx(:,:,:)
    allocate(labels2idx(sps%n, sps%n, sps%n))
    do loop = 1, 2
      ich = 0
      do j = 1, 2 * params%e3max + 3, 2
        do p = 1, -1, -2
          do itz = -3, 3, 2
            n = 0
            labels2idx(:,:,:) = 0
            do i1 = 1, sps%n
              j1 = sps%jj(i1)
              l1 = sps%ll(i1)
              itz1 = sps%itz(i1)
              do i2 = 1, i1
                j2 = sps%jj(i2)
                l2 = sps%ll(i2)
                itz2 = sps%itz(i2)
                if(sps%nshell(i1) + sps%nshell(i2) > params%e2max) cycle
                do i3 = 1, i2
                  j3 = sps%jj(i3)
                  l3 = sps%ll(i3)
                  itz3 = sps%itz(i3)
                  if(sps%nshell(i1) + sps%nshell(i3) > params%e2max) cycle
                  if(sps%nshell(i2) + sps%nshell(i3) > params%e2max) cycle
                  if(sps%nshell(i1) + sps%nshell(i2) + sps%nshell(i3) > params%e3max) cycle
                  if(itz1 + itz2 + itz3 /= itz) cycle
                  if((-1) ** (l1 + l2 + l3) /= p) cycle
                  call CntDim(sps, i1, i2, i3, j, ni, nj)
                  if(ni /= 0) n = n + 1
                  if(loop == 2 .and. ni /= 0) then
                    labels2idx(i1,i2,i3) = n
                  end if
                end do
              end do
            end do
            if(n /= 0) ich = ich + 1
            if(loop == 2 .and. n /= 0) then
              this%j(ich) = j
              this%p(ich) = p
              this%tz(ich) = itz
              this%jptz(ich)%nsub = n
              this%jptz2n(j, p, itz) = ich
              this%jptz(ich)%labels2nsub(:,:,:) = labels2idx(:,:,:)
            end if
          end do
        end do
      end do
      if(loop == 1) then
        this%n = ich
        allocate(this%jptz(ich))
        allocate(this%j(ich))
        allocate(this%p(ich))
        allocate(this%tz(ich))
        allocate(this%ndim(ich))
        allocate(this%jptz2n(1:2*params%e3max+3,-1:1,-3:3))
        this%j(:) = 0
        this%p(:) = 0
        this%tz(:) = 0
        this%ndim(:) = 0
        this%jptz2n(:,:,:) = 0
        do ich = 1, this%n
          allocate(this%jptz(ich)%labels2nsub(sps%n, sps%n, sps%n))
          this%jptz(ich)%labels2nsub(:,:,:) = 0
        end do
      end if
    end do
    deallocate(labels2idx)

    do ich = 1, this%n
      call this%jptz(ich)%init(sps, params, this%j(ich), this%p(ich), this%tz(ich))
      this%ndim(ich) = this%jptz(ich)%n
#ifdef debug
      write(*,'(a, i3, a, i3, a, i3, a, i8)') 'J = ', this%j(ich), '/2,  P = ', &
          & this%p(ich), ',  Tz = ', this%tz(ich), '/2,  n = ', this%jptz(ich)%n
#endif
    end do
#ifdef debug
    write(*,*)
#endif

  end subroutine InitThreeBodySpace

  subroutine FinThreeBodySpace(this)
    class(ThreeBodySpace), intent(inout) :: this
    integer :: ich
    do ich = 1, this%n
      call this%jptz(ich)%fin()
    end do
    deallocate(this%jptz2n)
    deallocate(this%j)
    deallocate(this%p)
    deallocate(this%tz)
  end subroutine FinThreeBodySpace

  subroutine InitThreeBodyChannel(this, sps, params, j, p, tz)
    class(ThreeBodyChannel), intent(inout) :: this
    type(spo_pn), intent(in) :: sps
    type(parameters), intent(in) :: params
    integer, intent(in) :: j, p, tz
    integer :: m, n, idx, ni, nj, nn
    integer :: i1, i2, i3

    allocate(this%idx(this%nsub))
    n = 0
    m = sps%n
    do i1 = 1, m
      do i2 = 1, i1
        if(sps%nshell(i1) + sps%nshell(i2) > params%e2max) cycle
        do i3 = 1, i2
          if(sps%nshell(i1) + sps%nshell(i3) > params%e2max) cycle
          if(sps%nshell(i2) + sps%nshell(i3) > params%e2max) cycle
          if(sps%nshell(i1) + sps%nshell(i2) + sps%nshell(i3) > params%e3max) cycle
          if(sps%itz(i1) + sps%itz(i2) + sps%itz(i3) /= tz) cycle
          if((-1) ** (sps%ll(i1) + sps%ll(i2) + sps%ll(i3)) /= p) cycle
          idx = this%labels2nsub(i1,i2,i3)
          if(idx == 0) cycle
          call CntDim(sps, i1, i2, i3, j, ni, nj)
          this%idx(idx)%nphys = ni
          this%idx(idx)%n = nj
          n = n + ni
          if(ni > 0) then
            allocate(this%idx(idx)%labels2n(ni))
            this%idx(idx)%labels2n = 0
            if(nj > 0) then
              allocate(this%idx(idx)%cfp(nj, ni))
              this%idx(idx)%cfp = 0.d0
              allocate(this%idx(idx)%n2label1(nj))
              allocate(this%idx(idx)%n2label2(nj))
              allocate(this%idx(idx)%n2label3(nj))
              allocate(this%idx(idx)%n2label4(nj))
              this%idx(idx)%n2label1 = 0
              this%idx(idx)%n2label2 = 0
              this%idx(idx)%n2label3 = 0
              this%idx(idx)%n2label4 = 0
            end if
          end if
        end do
      end do
    end do
    this%n = n
    if(n < 1) return
    allocate(this%n2label1(n))
    allocate(this%n2label2(n))
    allocate(this%n2label3(n))
    allocate(this%n2label4(n))
    this%n2label1 = 0
    this%n2label2 = 0
    this%n2label3 = 0
    this%n2label4 = 0

    n = 0
    do i1 = 1, m
      do i2 = 1, i1
        do i3 = 1, i2
          idx = this%labels2nsub(i1,i2,i3)
          if(idx == 0) cycle
          if(this%idx(idx)%nphys /= 0) then
            do nn = 1, this%idx(idx)%nphys
              n = n + 1
              this%idx(idx)%labels2n(nn) = n
              this%n2label1(n) = i1
              this%n2label2(n) = i2
              this%n2label3(n) = i3
              this%n2label4(n) = nn
            end do
          end if
          call this%idx(idx)%Init(sps, i1, i2, i3, j)
        end do
      end do
    end do

  end subroutine InitThreeBodyChannel

  subroutine FinThreeBodyChannel(this)
    class(ThreeBodyChannel), intent(inout) :: this
    integer :: idx
    deallocate(this%n2label1)
    deallocate(this%n2label2)
    deallocate(this%n2label3)
    deallocate(this%n2label4)
    do idx = 1, this%nsub
      call this%idx(idx)%Fin()
    end do
    deallocate(this%idx)
    deallocate(this%labels2nsub)
  end subroutine FinThreeBodyChannel

  subroutine InitAdditionalQN(this, sps, i1, i2, i3, j)
    use MatrixDouble, only: DMat
    use VectorDouble, only: DVec
    use LinAlgLib, only: EigenSolSymD, assignment(=), operator(+), &
        & operator(*), operator(-), operator(/)
    class(AdditionalQN), intent(inout) :: this
    type(spo_pn), intent(in) :: sps
    integer, intent(in) :: i1, i2, i3, j
    type(iter3) :: ite
    integer :: nj, nch, i, k, kk
    integer :: bra, ket
    integer :: ii1, ii2, ii3, jj1, jj2, jj3, jj12
    integer :: ii4, ii5, ii6, jj45
    type(DMat) :: mat
    type(EigenSolSymD) :: sol
    call ite%GenIter(i1, i2, i3, nch)
    nj = 0
    do i = 1, nch
      ii1 = ite%ii1(i)
      ii2 = ite%ii2(i)
      ii3 = ite%ii3(i)
      jj1 = sps%jj(ii1)
      jj2 = sps%jj(ii2)
      jj3 = sps%jj(ii3)
      do jj12 = iabs(jj1 - jj2) / 2, (jj1 + jj2) / 2
        if(ii1 == ii2 .and. mod(jj12, 2) == 1) cycle
        if(triag(2*jj12, jj3, j)) cycle
        nj = nj + 1
        this%n2label1(nj) = ii1
        this%n2label2(nj) = ii2
        this%n2label3(nj) = ii3
        this%n2label4(nj) = jj12
      end do
    end do
    if(this%nphys < 1) return
    call mat%Ini(this%n, this%n)
    do bra = 1, this%n
      ii1  = this%n2label1(bra)
      ii2  = this%n2label2(bra)
      ii3  = this%n2label3(bra)
      jj12 = this%n2label4(bra)
      do ket = 1, this%n
        ii4  = this%n2label1(ket)
        ii5  = this%n2label2(ket)
        ii6  = this%n2label3(ket)
        jj45 = this%n2label4(ket)
        mat%m(bra, ket) = Aop3(sps, ii1, ii2, ii3, jj12, &
            & ii4, ii5, ii6, jj45, j)
      end do
    end do
    call sol%init(mat)
    call sol%DiagSym(mat)
    do i = 1, this%n
      kk = 0
      do k = 1, this%n
        if(abs(1.d0-sol%eig%v(k)).le.1.d-4) then
          kk=kk+1
          this%cfp(i,kk) = sol%vec%m(i, k)
        end if
      end do
    end do
    call mat%Fin()
    call sol%fin()
  end subroutine InitAdditionalQN

  subroutine FinAdditionalQN(this)
    class(AdditionalQN), intent(inout) :: this
    if(this%nphys > 0) then
      deallocate(this%labels2n)
      deallocate(this%cfp)
    end if
    if(this%n > 0) then
      deallocate(this%n2label1)
      deallocate(this%n2label2)
      deallocate(this%n2label3)
      deallocate(this%n2label4)
    end if
  end subroutine FinAdditionalQN

  subroutine CntDim(sps, i1, i2, i3, j, ni, nj)
    type(spo_pn), intent(in) :: sps
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
      jj1 = sps%jj(ii1)
      jj2 = sps%jj(ii2)
      jj3 = sps%jj(ii3)
      do jj12 = iabs(jj1 - jj2) / 2, (jj1 + jj2) / 2
        if(ii1 == ii2 .and. mod(jj12, 2) == 1) cycle
        if(triag(2*jj12, jj3, j)) cycle
        if(jmin > jj12) jmin = jj12
        if(jmax < jj12) jmax = jj12
        nj = nj + 1
        tr = tr + Aop3(sps, ii1, ii2, ii3, jj12, ii1, ii2, ii3, jj12, j)
      end do
    end do
    !#ifdef debug
    !      write(*,'(a, 4i3, a, f12.6)') 'Dimension from trace: (', &
    !          & i1, i2, i3, j, '),  ', tr
    !#endif
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

    if(i2 == i3) then
      if(i1 == i2) then
        num = 1
        icase = 1
      else
        num = 3
        icase = 2
      end if
    else
      if(i1 == i2) then
        num = 3
        icase = 3
      elseif(i1 == i3) then
        num = 3
        icase = 4
      else
        num = 6
        icase = 5
      end if
    end if
    allocate(ite%ii1(num), ite%ii2(num), ite%ii3(num))
    select case(icase)
    case(1)
      ite%ii1(1) = i1; ite%ii2(1) = i1; ite%ii3(1) = i1
    case(2)
      ite%ii1(1) = i1; ite%ii2(1) = i2; ite%ii3(1) = i2
      ite%ii1(2) = i2; ite%ii2(2) = i1; ite%ii3(2) = i2
      ite%ii1(3) = i2; ite%ii2(3) = i2; ite%ii3(3) = i1
    case(3)
      ite%ii1(1) = i1; ite%ii2(1) = i1; ite%ii3(1) = i3
      ite%ii1(2) = i3; ite%ii2(2) = i1; ite%ii3(2) = i1
      ite%ii1(3) = i1; ite%ii2(3) = i3; ite%ii3(3) = i1
    case(4)
      ite%ii1(1) = i1; ite%ii2(1) = i2; ite%ii3(1) = i1
      ite%ii1(2) = i1; ite%ii2(2) = i1; ite%ii3(2) = i2
      ite%ii1(3) = i2; ite%ii2(3) = i1; ite%ii3(3) = i1
    case(5)
      ite%ii1(1) = i1; ite%ii2(1) = i2; ite%ii3(1) = i3
      ite%ii1(2) = i3; ite%ii2(2) = i1; ite%ii3(2) = i2
      ite%ii1(3) = i2; ite%ii2(3) = i3; ite%ii3(3) = i1
      ite%ii1(4) = i2; ite%ii2(4) = i1; ite%ii3(4) = i3
      ite%ii1(5) = i3; ite%ii2(5) = i2; ite%ii3(5) = i1
      ite%ii1(6) = i1; ite%ii2(6) = i3; ite%ii3(6) = i2
    end select
  end subroutine GenIter

  real(8) function Aop3(sps, i1, i2, i3, j12, &
        & i4, i5, i6, j45, j) result(a)
    type(spo_pn), intent(in) :: sps
    integer, intent(in) :: i1, i2, i3, i4, i5, i6
    integer, intent(in) :: j12, j45, j
    real(8) :: phase12, phase23, phase31
    phase12 = (-1.d0) ** ((sps%jj(i1) + sps%jj(i2))/2 - j45)
    phase23 = (-1.d0) ** ((sps%jj(i2) + sps%jj(i3))/2 - j45)
    phase31 = (-1.d0) ** ((sps%jj(i3) + sps%jj(i1))/2 - j45)
    a = 0.d0
    a = a + A3drct(i1, i2, i3, j12, i4, i5, i6, j45)
    a = a - A3drct(i1, i2, i3, j12, i5, i4, i6, j45) * phase12
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
    use common_library, only: sjs
    type(spo_pn), intent(in) :: sps
    integer, intent(in) :: i1, i2, i3, i4, i5, i6
    integer, intent(in) :: j12, j45, j
    integer :: j1, j2, j3
    e = 0.d0
    if(i1 /= i5) return
    if(i2 /= i6) return
    if(i3 /= i4) return
    j1 = sps%jj(i1)
    j2 = sps%jj(i2)
    j3 = sps%jj(i3)
    e = - (-1.d0) ** ((j1 + j2) / 2 + j12) * &
        & dsqrt(dble(2 * j12 + 1)) * dsqrt(dble(2 * j45 + 1)) * &
        & sjs(j1, j2, 2*j12, j, j3, 2*j45)
  end function A3exc1

  real(8) function A3exc2(sps, i1, i2, i3, j12, i4, i5, i6, j45, j) result(e)
    use common_library, only: sjs
    type(spo_pn), intent(in) :: sps
    integer, intent(in) :: i1, i2, i3, i4, i5, i6
    integer, intent(in) :: j12, j45, j
    integer :: j1, j2, j3
    e = 0.d0
    if(i1 /= i6) return
    if(i2 /= i4) return
    if(i3 /= i5) return
    j1 = sps%jj(i1)
    j2 = sps%jj(i2)
    j3 = sps%jj(i3)
    e = - (-1.d0) ** ((j2 + j3) / 2 + j45) * &
        & dsqrt(dble(2 * j12 + 1)) * dsqrt(dble(2 * j45 + 1)) * &
        & sjs(j1, j2, 2*j12, j3, j, 2*j45)
  end function A3exc2

end module ModelSpace
