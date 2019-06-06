module ThreeBodyMonopoleSpace
  use omp_lib
  use SingleParticleState
  implicit none

  public :: ThreeBodyMonSpace
  public :: ThreeBodyMonChannel
  public :: ThreeBodyMonABC

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

  type :: coef
    integer :: n_sum = 0
    integer, allocatable :: idx2n(:)
#ifdef single_precision_three_body_file
    real(4), allocatable :: TrnsCoef(:)
#else
    real(8), allocatable :: TrnsCoef(:)
#endif
  end type coef

  type :: ThreeBodyMonABC
    integer :: j12min, j12max
    integer :: t12min, t12max
    type(coef), allocatable :: jt(:,:)
  contains
    procedure :: FinThreeBodyMonABC
    generic :: fin => FinThreeBodyMonABC
  end type ThreeBodyMonABC

  type :: ThreeBodyMonSorting
    integer :: j, p, t
    integer :: n_idx = 0
    integer, allocatable :: abc2n(:,:,:,:)
    type(ThreeBodyMonABC), allocatable :: idx(:)
  contains
    procedure :: init => InitThreeBodyMonSorting
    procedure :: fin => FinThreeBodyMonSorting
  end type ThreeBodyMonSorting

  type :: ThreeBodyMonChannel
    integer :: j, p, t
    integer :: ch1, ch2, ch3
    integer :: n_state
    integer, allocatable :: n2abcjt(:,:)
    integer, allocatable :: nsjt2n(:,:,:,:,:)
  contains
    procedure :: init => InitThreeBodyMonChannel
    procedure :: fin => FinThreeBodyMonChannel
  end type ThreeBodyMonChannel

  type :: ThreeBodyMonSpace
    integer :: emax, e2max, e3max
    type(Orbits), pointer :: sps
    type(OrbitsIsospin), pointer :: isps

    ! Three-body space for monopole interaction
    type(ThreeBodyMonChannel), allocatable :: chan(:)
    type(OneBodyChannels) :: one
    integer, allocatable :: JPTMon2Ch(:,:)
    integer, allocatable :: JPT2Ch(:,:,:)
    integer, allocatable :: Mon2Ch(:,:,:)
    integer :: NChan, NJPT, NMon

    ! Three-body space for sorting monopole interaction
    type(ThreeBodyMonSorting), allocatable :: sorted(:)
    integer, allocatable :: JPT2SortedCh(:,:,:)
    integer :: NChanSorted
  contains
    procedure :: init => InitThreeBodyMonSpace
    procedure :: fin => FinThreeBodyMonSpace
  end type ThreeBodyMonSpace

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
    do ich = 1, this%NChan
      call this%chan(ich)%fin()
    end do
    do ich = 1, this%NChanSorted
      call this%sorted(ich)%fin()
    end do
    deallocate(this%JPTMon2Ch)
    deallocate(this%JPT2Ch)
    deallocate(this%Mon2Ch)
    deallocate(this%JPT2SortedCh)
    deallocate(this%Chan)
    deallocate(this%sorted)
    call this%one%fin()
    this%sps => null()
    this%isps => null()
  end subroutine FinThreeBodyMonSpace

  subroutine InitThreeBodyMonSpace(this, sps, isps, e2max, e3max)
    use MyLibrary, only: triag
    class(ThreeBodyMonSpace), intent(inout) :: this
    type(Orbits), target, intent(in) :: sps
    type(OrbitsIsospin), target, intent(in) :: isps
    integer, intent(in) :: e2max, e3max
    integer :: jtot, ptot, ttot
    integer :: cnt
    integer :: ch1, ch2, ch3
    integer :: j1, j2, j3, j12, p1, p2, p3
    integer :: i1, i2, i3, t12
    type(SingleParticleOrbitIsospin), pointer :: o1, o2, o3
    integer :: ich, ich_JPT, ich_Mon
    integer, allocatable :: jj(:), pp(:), tt(:)
    integer, allocatable :: ch1_(:), ch2_(:), ch3_(:)
    this%sps => sps
    this%isps => isps
    this%emax = sps%emax
    this%e2max = e2max
    this%e3max = e3max

    call this%one%init(isps)

    this%NJPT = (2*e3max+3) * 4
    this%NMon = this%one%NChan**3
    allocate(this%JPTMon2Ch(this%NJPT, this%NMon))
    allocate(this%JPT2Ch(2 * e3max + 3, -1:1, 1:3))
    allocate(this%Mon2Ch(this%one%NChan, this%one%NChan, this%one%NChan))
    this%JPTMon2Ch(:,:) = 0
    this%JPT2Ch(:,:,:) = 0
    this%Mon2Ch(:,:,:) = 0
    ich = 0
    do ttot = 1, 3, 2
      do jtot = 1, 2 * e3max + 3, 2
        do ptot = 1, -1, -2

          do ch1 = 1, this%one%NChan
            j1 = this%one%j(ch1)
            p1 = this%one%p(ch1)
            do ch2 = 1, this%one%NChan
              j2 = this%one%j(ch2)
              p2 = this%one%p(ch2)
              do ch3 = 1, this%one%NChan
                j3 = this%one%j(ch3)
                p3 = this%one%p(ch3)
                if( p1*p2*p3 .ne. ptot) cycle

                cnt = 0
                do i1 = 1, sps%norbs
                  o1 => isps%orb(i1)
                  if(j1 /= o1%j) cycle
                  if(p1 /= (-1)**o1%l) cycle
                  do i2 = 1, i1
                    o2 => isps%orb(i2)
                    if(j2 /= o2%j) cycle
                    if(p2 /= (-1)**o2%l) cycle
                    if(o1%e + o2%e > e2max) cycle
                    do i3 = 1, i2
                      o3 => isps%orb(i3)
                      if(j3 /= o3%j) cycle
                      if(p3 /= (-1)**o3%l) cycle
                      if(o1%e + o3%e > e2max) cycle
                      if(o2%e + o3%e > e2max) cycle
                      if(o1%e + o2%e + o3%e > e3max) cycle
                      do j12 = abs(j1 - j2)/2, (j1 + j2)/2
                        if(triag(2*j12, j3, jtot)) cycle
                        do t12 = 0, 1
                          if(triag(2*t12, 1, ttot)) cycle
                          if(i1 == i2 .and. mod(j12 + t12, 2) == 0) cycle
                          cnt = cnt + 1
                        end do
                      end do
                    end do
                  end do
                end do
                if(cnt /= 0) ich = ich + 1

              end do
            end do
          end do
        end do
      end do
    end do
    this%NChan = ich
    allocate(jj(this%NChan))
    allocate(pp(this%NChan))
    allocate(tt(this%NChan))
    allocate(ch1_(this%NChan))
    allocate(ch2_(this%NChan))
    allocate(ch3_(this%NChan))

    ich = 0
    ich_JPT = 0
    do ttot = 1, 3, 2
      do jtot = 1, 2 * e3max + 3, 2
        do ptot = 1, -1, -2
          ich_JPT = ich_JPT+1

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
                if( p1*p2*p3 .ne. ptot) cycle

                cnt = 0
                do i1 = 1, sps%norbs
                  o1 => isps%orb(i1)
                  if(j1 /= o1%j) cycle
                  if(p1 /= (-1)**o1%l) cycle
                  do i2 = 1, i1
                    o2 => isps%orb(i2)
                    if(j2 /= o2%j) cycle
                    if(p2 /= (-1)**o2%l) cycle
                    if(o1%e + o2%e > e2max) cycle
                    do i3 = 1, i2
                      o3 => isps%orb(i3)
                      if(j3 /= o3%j) cycle
                      if(p3 /= (-1)**o3%l) cycle
                      if(o1%e + o3%e > e2max) cycle
                      if(o2%e + o3%e > e2max) cycle
                      if(o1%e + o2%e + o3%e > e3max) cycle
                      do j12 = abs(j1-j2)/2, (j1+j2)/2
                        if(triag(2*j12, j3, jtot)) cycle
                        do t12 = 0, 1
                          if(triag(2*t12, 1, ttot)) cycle
                          if(i1 == i2 .and. mod(j12 + t12, 2) == 0) cycle
                          cnt = cnt + 1
                        end do
                      end do
                    end do
                  end do
                end do
                if(cnt /= 0) then
                  ich = ich + 1
                  this%JPTMon2Ch(ich_JPT, ich_Mon) = ich
                  this%JPT2Ch(jtot, ptot, ttot) = ich_JPT
                  this%Mon2Ch(ch1, ch2, ch3) = ich_Mon
                  jj(ich) = jtot
                  pp(ich) = ptot
                  tt(ich) = ttot
                  ch1_(ich) = ch1
                  ch2_(ich) = ch2
                  ch3_(ich) = ch3
                end if
              end do
            end do
          end do

        end do
      end do
    end do

    allocate(this%chan(this%NChan))
    do ich = 1, this%NChan
      call this%chan(ich)%init(isps, this%one, jj(ich), pp(ich), tt(ich), &
          & ch1_(ich), ch2_(ich), ch3_(ich), e2max, e3max)
    end do
    deallocate(jj)
    deallocate(pp)
    deallocate(tt)
    deallocate(ch1_)
    deallocate(ch2_)
    deallocate(ch3_)

    ich = 0
    do ttot = 1, 3, 2
      do jtot = 1, 2 * e3max + 3, 2
        do ptot = 1, -1, -2

          cnt = 0
          do i1 = 1, sps%norbs
            o1 => isps%orb(i1)
            j1 = o1%j
            do i2 = 1, i1
              o2 => isps%orb(i2)
              j2 = o2%j
              if(o1%e + o2%e > e2max) cycle
              do i3 = 1, i2
                o3 => isps%orb(i3)
                j3 = o3%j
                if(j3 /= o3%j) cycle
                if(o1%e + o3%e > e2max) cycle
                if(o2%e + o3%e > e2max) cycle
                if(o1%e + o2%e + o3%e > e3max) cycle
                if( (-1)**(o1%l+o2%l+o3%l) .ne. ptot) cycle

                do j12 = abs(j1-j2)/2, (j1+j2)/2
                  if(triag(2*j12, j3, jtot)) cycle
                  do t12 = 0, 1
                    if(triag(2*t12, 1, ttot)) cycle
                    if(i1 == i2 .and. mod(j12 + t12, 2) == 0) cycle
                    cnt = cnt + 1
                  end do
                end do
              end do
            end do
          end do
          if(cnt /= 0) ich = ich + 1
        end do
      end do
    end do
    this%NChanSorted = ich
    allocate(jj(this%NChanSorted))
    allocate(pp(this%NChanSorted))
    allocate(tt(this%NChanSorted))

    ich = 0
    do ttot = 1, 3, 2
      do jtot = 1, 2 * e3max + 3, 2
        do ptot = 1, -1, -2

          cnt = 0
          do i1 = 1, sps%norbs
            o1 => isps%orb(i1)
            j1 = o1%j
            do i2 = 1, i1
              o2 => isps%orb(i2)
              j2 = o2%j
              if(o1%e + o2%e > e2max) cycle
              do i3 = 1, i2
                o3 => isps%orb(i3)
                j3 = o3%j
                if(j3 /= o3%j) cycle
                if(o1%e + o3%e > e2max) cycle
                if(o2%e + o3%e > e2max) cycle
                if(o1%e + o2%e + o3%e > e3max) cycle
                if( (-1)**(o1%l+o2%l+o3%l) .ne. ptot) cycle

                do j12 = abs(j1-j2)/2, (j1+j2)/2
                  if(triag(2*j12, j3, jtot)) cycle
                  do t12 = 0, 1
                    if(triag(2*t12, 1, ttot)) cycle
                    if(i1 == i2 .and. mod(j12 + t12, 2) == 0) cycle
                    cnt = cnt + 1
                    jj(ich) = jtot
                    pp(ich) = ptot
                    tt(ich) = ttot
                  end do
                end do
              end do
            end do
          end do
          if(cnt /= 0) ich = ich + 1
        end do
      end do
    end do

    allocate(this%sorted(this%NChanSorted))
    do ich = 1, this%NChanSorted
      call this%sorted(ich)%init(isps, this, jj(ich), pp(ich), tt(ich), e2max, e3max)
    end do
    deallocate(jj)
    deallocate(pp)
    deallocate(tt)

  end subroutine InitThreeBodyMonSpace

  subroutine FinThreeBodyMonChannel(this)
    class(ThreeBodyMonChannel), intent(inout) :: this
    deallocate(this%n2abcjt)
    deallocate(this%nsjt2n)
  end subroutine FinThreeBodyMonChannel

  subroutine InitThreeBodyMonChannel(this, sps, one, j, p, t, ch1, ch2, ch3, e2max, e3max)
    use MyLibrary, only: triag
    class(ThreeBodyMonChannel), intent(inout) :: this
    type(OrbitsIsospin), target, intent(in) :: sps
    type(OneBodyChannels), intent(in) :: one
    integer, intent(in) :: j, p, t, ch1, ch2, ch3, e2max, e3max
    type(SingleParticleOrbitIsospin), pointer :: o1, o2, o3
    integer :: i1, i2, i3, cnt
    integer :: j1, j2, j3, j12, t12, p1, p2, p3

    this%j = j
    this%p = p
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
    allocate(this%nsjt2n(0:sps%emax/2, 0:sps%emax/2, 0:sps%emax/2, 0:max(2*sps%emax, e2max)+1, 0:1))
    this%nsjt2n(:,:,:,:,:) = 0

    cnt = 0
    do i1 = 1, sps%norbs
      o1 => sps%orb(i1)
      if(j1 /= o1%j) cycle
      if(p1 /= (-1)**o1%l) cycle
      do i2 = 1, i1
        o2 => sps%orb(i2)
        if(j2 /= o2%j) cycle
        if(p2 /= (-1)**o2%l) cycle
        if(o1%e + o2%e > e2max) cycle
        do i3 = 1, i2
          o3 => sps%orb(i3)
          if(j3 /= o3%j) cycle
          if(p3 /= (-1)**o3%l) cycle
          if(o1%e + o3%e > e2max) cycle
          if(o2%e + o3%e > e2max) cycle
          if(o1%e + o2%e + o3%e > e3max) cycle
          do j12 = abs(j1-j2)/2, (j1+j2)/2
            if(triag(2*j12, j3, j)) cycle
            do t12 = 0, 1
              if(triag(2*t12, 1, t)) cycle
              if(i1 == i2 .and. mod(j12 + t12, 2) == 0) cycle
              cnt = cnt + 1
            end do
          end do
        end do
      end do
    end do
    this%n_state = cnt
    allocate(this%n2abcjt(5,this%n_state))

    cnt = 0
    do i1 = 1, sps%norbs
      o1 => sps%orb(i1)
      if(j1 /= o1%j) cycle
      if(p1 /= (-1)**o1%l) cycle
      do i2 = 1, i1
        o2 => sps%orb(i2)
        if(j2 /= o2%j) cycle
        if(p2 /= (-1)**o2%l) cycle
        if(o1%e + o2%e > e2max) cycle
        do i3 = 1, i2
          o3 => sps%orb(i3)
          if(j3 /= o3%j) cycle
          if(p3 /= (-1)**o3%l) cycle
          if(o1%e + o3%e > e2max) cycle
          if(o2%e + o3%e > e2max) cycle
          if(o1%e + o2%e + o3%e > e3max) cycle
          do j12 = abs(j1-j2)/2, (j1+j2)/2
            if(triag(2*j12, j3, j)) cycle
            do t12 = 0, 1
              if(triag(2*t12, 1, t)) cycle
              if(i1 == i2 .and. mod(j12 + t12, 2) == 0) cycle
              cnt = cnt + 1
              this%n2abcjt(:,cnt) = [i1,i2,i3,j12,t12]
              this%nsjt2n(o1%n, o2%n, o3%n, j12, t12) = cnt
            end do
          end do
        end do
      end do
    end do
  end subroutine InitThreeBodyMonChannel

  subroutine FinThreeBodyMonSorting(this)
    class(ThreeBodyMonSorting), intent(inout) :: this
    integer :: idx
    do idx = 1, this%n_idx
      call this%idx(idx)%fin()
    end do
    deallocate(this%abc2n)
    deallocate(this%idx)
  end subroutine FinThreeBodyMonSorting

  subroutine FinThreeBodyMonABC(this)
    class(ThreeBodyMonABC), intent(inout) :: this
    integer :: j12, t12
    do j12 = this%j12min, this%j12max
      do t12 = this%t12min, this%t12max
        if(allocated(this%jt)) then
          deallocate(this%jt(j12,t12)%TrnsCoef)
          deallocate(this%jt(j12,t12)%idx2n)
        end if
      end do
    end do
    deallocate(this%jt)
  end subroutine FinThreeBodyMonABC

  subroutine InitThreeBodyMonSorting(this, sps, space, j, p, t, e2max, e3max)
    use MyLibrary, only: triag
    class(ThreeBodyMonSorting), intent(inout) :: this
    type(ThreeBodyMonSpace), intent(in) :: space
    type(OrbitsIsospin), intent(in) :: sps
    integer, intent(in) :: j, p, t, e2max, e3max
    integer :: ich, ich_jpt, ich_mon
    integer :: i1, l1, j1, e1
    integer :: i2, l2, j2, e2
    integer :: i3, l3, j3, e3
    integer :: a, b, c, n_recouple, cha, chb, chc
    integer :: cnt, idx
    integer :: j12min, j12, j12max
    integer :: t12min, t12, t12max

    allocate(this%abc2n(4,sps%norbs,sps%norbs,sps%norbs))
    this%abc2n(:,:,:,:) = 0
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
          this%abc2n(1,i1,i2,i3) = idx
        end do
      end do
    end do

    this%n_idx = idx
    allocate(this%idx(this%n_idx))
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

          cha = space%one%jp2ch(sps%orb(a)%j, (-1)**sps%orb(a)%l)
          chb = space%one%jp2ch(sps%orb(b)%j, (-1)**sps%orb(b)%l)
          chc = space%one%jp2ch(sps%orb(c)%j, (-1)**sps%orb(c)%l)
          idx        = this%abc2n(1,i1,i2,i3)
          this%abc2n(2:4,i1,i2,i3) = [cha,chb,chc]
          if(idx == 0) cycle

          ich_jpt = space%JPT2Ch(j,p,t)
          ich_Mon = space%Mon2Ch(cha,chb,chc)
          if(ich_jpt * ich_Mon == 0) cycle
          ich = space%JPTMon2Ch(ich_jpt, ich_Mon)
          if(ich == 0) cycle

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


          this%idx(idx)%j12min = j12min
          this%idx(idx)%j12max = j12max
          this%idx(idx)%t12min = t12min
          this%idx(idx)%t12max = t12max
          allocate(this%idx(idx)%jt(j12min:j12max,t12min:t12max))
          call set_trans_coef_for_sorting(this%idx(idx),space%chan(ich), sps,&
              & i1,i2,i3,a,b,c,n_recouple,j,t)
        end do
      end do
    end do
  end subroutine InitThreeBodyMonSorting

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

  subroutine set_trans_coef_for_sorting(this, chan, sps, i1, i2, i3, a, b, c, &
        &    n_recouple, j, t)
    use MyLibrary, only: triag, sjs, hat
    type(ThreeBodyMonABC), intent(inout) :: this
    type(ThreeBodyMonChannel), intent(in) :: chan
    type(OrbitsIsospin), intent(in) :: sps
    integer, intent(in) :: i1, i2, i3, a, b, c, n_recouple, j, t
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
            this%jt(j12, t12)%n_sum = n
            allocate(this%jt(j12,t12)%idx2n( n))
            allocate(this%jt(j12,t12)%TrnsCoef(n))
            this%jt(j12,t12)%idx2n(:) = 0
#ifdef single_precision_three_body_file
            this%jt(j12,t12)%TrnsCoef(:) = 0.0
#else
            this%jt(j12,t12)%TrnsCoef(:) = 0.d0
#endif
          elseif(loop == 2) then
            n = 0
            select case(n_recouple)
            case(1)
              n = 1
#ifdef single_precision_three_body_file
              this%jt(j12,t12)%TrnsCoef(n) = 1.0
#else
              this%jt(j12,t12)%TrnsCoef(n) = 1.d0
#endif
              this%jt(j12, t12)%idx2n(n) = chan%nsjt2n(sps%orb(a)%n, sps%orb(b)%n, sps%orb(c)%n, jab, tab)
            case(2)
              do jab = iabs(ja - jb)/2, (ja + jb)/2
                do tab = 0, 1
                  if(a == b .and. mod(jab + tab, 2) == 0) cycle
                  if(triag(2 * tab,  1, t)) cycle
                  if(triag(2 * jab, jc, j)) cycle
                  n = n + 1
#ifdef single_precision_three_body_file
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
                  this%jt(j12, t12)%idx2n(n) = chan%nsjt2n(sps%orb(a)%n, sps%orb(b)%n, sps%orb(c)%n, jab, tab)
                end do
              end do
            case(3)
              do jab = iabs(ja - jb)/2, (ja + jb)/2
                do tab = 0, 1
                  if(a == b .and. mod(jab + tab, 2) == 0) cycle
                  if(triag(2 * tab,  1, t)) cycle
                  if(triag(2 * jab, jc, j)) cycle
                  n = n + 1
#ifdef single_precision_three_body_file
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
                  this%jt(j12, t12)%idx2n(n) = chan%nsjt2n(sps%orb(a)%n, sps%orb(b)%n, sps%orb(c)%n, jab, tab)
                end do
              end do
            case(4)
              n = n + 1
#ifdef single_precision_three_body_file
              this%jt(j12,t12)%TrnsCoef(n) = real(&
                  & (-1.d0) ** ((ja + jb) / 2 - j12 - t12))
#else
              this%jt(j12,t12)%TrnsCoef(n) = &
                  & (-1.d0) ** ((ja + jb) / 2 - j12 - t12)
#endif
              this%jt(j12, t12)%idx2n(n) = chan%nsjt2n(sps%orb(a)%n, sps%orb(b)%n, sps%orb(c)%n, jab, tab)
            case(5)
              do jab = iabs(ja - jb)/2, (ja + jb)/2
                do tab = 0, 1
                  if(a == b .and. mod(jab + tab, 2) == 0) cycle
                  if(triag(2 * tab,  1, t)) cycle
                  if(triag(2 * jab, jc, j)) cycle
                  n = n + 1
#ifdef single_precision_three_body_file
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
                  this%jt(j12, t12)%idx2n(n) = chan%nsjt2n(sps%orb(a)%n, sps%orb(b)%n, sps%orb(c)%n, jab, tab)
                end do
              end do
            case(6)
              do jab = iabs(ja - jb)/2, (ja + jb)/2
                do tab = 0, 1
                  if(a == b .and. mod(jab + tab, 2) == 0) cycle
                  if(triag(2 * tab,  1, t)) cycle
                  if(triag(2 * jab, jc, j)) cycle
                  n = n + 1
#ifdef single_precision_three_body_file
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
                  this%jt(j12, t12)%idx2n(n) = chan%nsjt2n(sps%orb(a)%n, sps%orb(b)%n, sps%orb(c)%n, jab, tab)
                end do
              end do

            end select
          end if
        end do
        !write(*,'(5i3,a,3i3)') i1,i2,i3,j12,t12, ", ", a,b,c
        !write(*,'(10f10.4)') this%jt(j12,t12)%TrnsCoef
      end do
    end do
  end subroutine set_trans_coef_for_sorting
end module ThreeBodyMonopoleSpace
