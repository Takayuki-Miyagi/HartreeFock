module ThreeBodyInteraction
  use omp_lib
  use SingleParticleState
  implicit none

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
    integer :: j = -1
    integer :: p = 0
    integer :: t = -1
    integer :: n_state = 0
    integer :: n_idx = 0
    integer :: n_sort = 0
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
    type(Orbits), pointer :: sps
    type(OrbitsIsospin), pointer :: isps
    integer, allocatable :: jpt2ch(:,:,:)
    integer :: NChan
    integer :: emax, e2max, e3max
  contains
    procedure :: init => InitNonOrthIsospinThreeBodySpace
    procedure :: fin => FinNonOrthIsospinThreeBodySpace
  end type NonOrthIsospinThreeBodySpace

  type :: ThreeBodyForceChannel
#ifdef single_precision_three_body_force
    real(4), allocatable :: v(:)
#else
    real(8), allocatable :: v(:)
#endif
  end type ThreeBodyForceChannel

  type :: ThreeBodyForce
    type(ThreeBodyForceChannel), allocatable :: jpt(:)
    type(NonOrthIsospinThreeBodySpace), pointer :: thr
  contains
    procedure :: InitThreeBodyForce
    procedure :: FinThreeBodyForce
    procedure :: GetThBME_isospin
    procedure :: GetThBME_pn

    generic :: init => InitThreeBodyForce
    generic :: fin => FinThreeBodyForce
    generic :: GetThBME => GetThBME_isospin, GetThBME_pn
  end type ThreeBodyForce

contains
  subroutine InitThreeBodyForce(this, thr)
    class(ThreeBodyForce), intent(inout) :: this
    type(NonOrthIsospinThreeBodySpace), target, intent(in) :: thr
    integer :: ch, n_state, n

    this%thr => thr
    allocate(this%jpt(thr%NChan))

    do ch = 1, thr%NChan
      n_state = thr%jpt(ch)%n_state
      n = n_state * (n_state+1) / 2
      allocate(this%jpt(ch)%v(n))
    end do
  end subroutine InitThreeBodyForce

  subroutine FinThreeBodyForce(this)
    class(ThreeBodyForce), intent(inout) :: this
    integer :: ch
    do ch = 1, this%thr%NChan
      deallocate(this%jpt(ch)%v)
    end do
    deallocate(this%jpt)
    this%thr => null()
  end subroutine FinThreeBodyForce

  function GetThBME_pn(this,i1,i2,i3,J12,&
        & i4,i5,i6,J45,J) result(r)
    use MyLibrary, only: dcg
    class(ThreeBodyForce), intent(in) :: this
    integer, intent(in) :: i1,i2,i3,i4,i5,i6
    integer, intent(in) :: J12,J45,J
    type(SingleParticleOrbit), pointer :: o1,o2,o3,o4,o5,o6
    integer :: z1, z2, z3, z4, z5, z6, T12, T45, T
    integer :: a, b, c, d, e, f
    integer :: P, Z, ch
    real(8) :: r
    r = 0.d0
    o1 => this%thr%sps%GetOrbit(i1)
    o2 => this%thr%sps%GetOrbit(i2)
    o3 => this%thr%sps%GetOrbit(i3)
    o4 => this%thr%sps%GetOrbit(i4)
    o5 => this%thr%sps%GetOrbit(i5)
    o6 => this%thr%sps%GetOrbit(i6)
    if(o1%e + o2%e + o3%e > this%thr%e3max) return
    if(o4%e + o5%e + o6%e > this%thr%e3max) return
    if(i1 == i2 .and. mod(J12, 2) == 1) return
    if(i4 == i5 .and. mod(J45, 2) == 1) return
    z1 = o1%z
    z2 = o2%z
    z3 = o3%z
    z4 = o4%z
    z5 = o5%z
    z6 = o6%z

    a = this%thr%isps%nlj2idx( o1%n, o1%l, o1%j )
    b = this%thr%isps%nlj2idx( o2%n, o2%l, o2%j )
    c = this%thr%isps%nlj2idx( o3%n, o3%l, o3%j )
    d = this%thr%isps%nlj2idx( o4%n, o4%l, o4%j )
    e = this%thr%isps%nlj2idx( o5%n, o5%l, o5%j )
    f = this%thr%isps%nlj2idx( o6%n, o6%l, o6%j )

    P = (-1) ** (o1%l+o2%l+o3%l)
    Z = z1 + z2 + z3
    do T12 = 0, 1
      if(abs(z1+z2) > 2*T12) cycle
      do T45 = 0, 1
        if(abs(z4+z5) > 2*T45) cycle
        do T = max(abs(2*T12-1),abs(2*T45-1)), min(2*T12+1,2*T45+1), 2
          if(abs(Z) > T) cycle
          ch = this%thr%jpt2ch(J,P,T)
          if(ch == 0) cycle
          r = r + &
              & this%GetThBME(a,b,c,J12,T12,&
              & d,e,f,J45,T45,J,T) * &
              & dcg(1,z1,1,z2,2*T12,z1+z2) * dcg(2*T12,z1+z2,1,z3,T,Z) * &
              & dcg(1,z4,1,z5,2*T45,z4+z5) * dcg(2*T45,z4+z5,1,z6,T,Z)
        end do
      end do
    end do
  end function GetThBME_pn

  function GetThBME_Isospin(this,i1,i2,i3,J12,T12,&
        & i4,i5,i6,J45,T45,J,T) result(r)
    class(ThreeBodyForce), intent(in) :: this
    integer, intent(in) :: i1,i2,i3,i4,i5,i6
    integer, intent(in) :: J12,T12,J45,T45,J,T
    type(SingleParticleOrbitIsospin), pointer :: o1,o2,o3,o4,o5,o6
    type(NonOrthIsospinThreeBodySpace), pointer :: tbs
    integer :: ch, idxbra, idxket, bra, ket
    integer :: P123, P456
    integer :: isorted_bra, isorted_ket
    integer :: ibra, iket, nmax, nmin, idx_vec
    real(8) :: r

    r = 0.d0
    if(i1 == i2 .and. mod(J12+T12,2) == 0) return
    if(i4 == i5 .and. mod(J45+T45,2) == 0) return
    o1 => this%thr%isps%GetOrbit(i1)
    o2 => this%thr%isps%GetOrbit(i2)
    o3 => this%thr%isps%GetOrbit(i3)
    o4 => this%thr%isps%GetOrbit(i4)
    o5 => this%thr%isps%GetOrbit(i5)
    o6 => this%thr%isps%GetOrbit(i6)
    tbs => this%thr
    P123 = (-1) ** (o1%l+o2%l+o3%l)
    P456 = (-1) ** (o4%l+o5%l+o6%l)
    if(P123 * P456 /= 1) then
      write(*,*) 'Warning: in GetThBMEIso_scalar: P'
      return
    end if
    ch = tbs%jpt2ch(J,P123,T)
    if(ch == 0) return
    idxbra = tbs%jpt(ch)%sorting(i1,i2,i3)
    idxket = tbs%jpt(ch)%sorting(i4,i5,i6)
    if(idxbra * idxket == 0) return
    if(i1 == i2 .and. mod(J12+T12,2) == 0) return
    if(i4 == i5 .and. mod(J45+T45,2) == 0) return
    isorted_bra = tbs%jpt(ch)%sort(idxbra)%idx_sorted
    isorted_ket = tbs%jpt(ch)%sort(idxket)%idx_sorted
    if(isorted_bra * isorted_ket == 0) return

    do ibra = 1, tbs%jpt(ch)%sort(idxbra)%JT(J12,T12)%n
      bra = tbs%jpt(ch)%sort(idxbra)%JT(J12,T12)%idx2num(ibra)
      do iket = 1, tbs%jpt(ch)%sort(idxket)%JT(J45,T45)%n
        ket = tbs%jpt(ch)%sort(idxket)%JT(J45,T45)%idx2num(iket)
        nmax = max(bra,ket)
        nmin = min(bra,ket)
        idx_vec = nmax * (nmax - 1) / 2 + nmin
        r = r + dble(this%jpt(ch)%v(idx_vec) * &
            & tbs%jpt(ch)%sort(idxbra)%JT(J12,T12)%TrnsCoef(ibra) * &
            & tbs%jpt(ch)%sort(idxket)%JT(J45,T45)%TrnsCoef(iket))
      end do
    end do
  end function GetThBME_Isospin

  subroutine FinNonOrthIsospinThreeBodySpace(this)
    class(NonOrthIsospinThreeBodySpace), intent(inout) :: this
    integer :: ich
    do ich = 1, this%NChan
      call this%jpt(ich)%fin()
    end do
    deallocate(this%jpt)
    deallocate(this%jpt2ch)
  end subroutine FinNonOrthIsospinThreeBodySpace

  subroutine InitNonOrthIsospinThreeBodySpace(this, isps, sps, e2max, e3max)
    use MyLibrary, only: triag
    class(NonOrthIsospinThreeBodySpace), intent(inout) :: this
    type(OrbitsIsospin), intent(in), target :: isps
    type(Orbits), intent(in), target :: sps
    integer, intent(in) :: e2max, e3max
    integer :: j, p, t, ich, nidx, n
    integer :: i1, i2, i3
    type(SingleParticleOrbitIsospin), pointer :: o1, o2, o3
    integer :: j12, t12, j12min, j12max, t12min, t12max
    integer :: nj
    integer, allocatable :: jj(:), pp(:), tt(:), nn(:), nnidx(:)
    type :: spis2n
      integer, allocatable :: spis2nidx(:,:,:)
    end type spis2n
    type(spis2n), allocatable :: ch(:)

    this%sps => sps
    this%isps => isps
    allocate(this%jpt2ch(1:2*min(e3max,3*sps%lmax)+3,-1:1,1:3))
    this%jpt2ch(:,:,:) = 0
    ich = 0
    do j = 1, 2*min(e3max,3*sps%lmax)+3, 2
      do p = 1, -1, -2
        do t = 1, 3, 2
          n = 0
          nidx = 0
          do i1 = 1, sps%norbs
            o1 => isps%GetOrbit(i1)
            do i2 = 1, i1
              o2 => isps%GetOrbit(i2)
              if(o1%e + o2%e > e2max) cycle
              do i3 = 1, i2
                o3 => isps%GetOrbit(i3)
                if(o1%e + o3%e > e2max) cycle
                if(o2%e + o3%e > e2max) cycle
                if(o1%e + o2%e + o3%e > e3max) cycle
                if((-1) ** (o1%l+o2%l+o3%l) /= p) cycle

                nj = 0
                j12min = (o1%j+o2%j)/2
                j12max = abs(o1%j-o2%j)/2
                t12min = 1
                t12max = 0
                do j12 = abs(o1%j-o2%j)/2, (o1%j+o2%j)/2
                  if(triag(2*j12,o3%j,j)) cycle
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
            o1 => isps%GetOrbit(i1)
            do i2 = 1, i1
              o2 => isps%GetOrbit(i2)
              if(o1%e + o2%e > e2max) cycle
              do i3 = 1, i2
                o3 => isps%GetOrbit(i3)
                if(o1%e + o3%e > e2max) cycle
                if(o2%e + o3%e > e2max) cycle
                if(o1%e + o2%e + o3%e > e3max) cycle
                if((-1) ** (o1%l+o2%l+o3%l) /= p) cycle

                nj = 0
                j12min = (o1%j+o2%j)/2
                j12max = abs(o1%j-o2%j)/2
                t12min = 1
                t12max = 0
                do j12 = abs(o1%j-o2%j)/2, (o1%j+o2%j)/2
                  if(triag(2*j12,o3%j,j)) cycle
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
          & ch(ich)%spis2nidx, isps, e2max, e3max)
      deallocate(ch(ich)%spis2nidx)
      call this%jpt(ich)%set( jj(ich),pp(ich),tt(ich), isps, e2max, e3max)
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
    use MyLibrary, only: triag
    class(NonOrthIsospinThreeBodyChannel), intent(inout) :: this
    integer, intent(in) :: j, p, t, n, nidx, spis2idx(:,:,:), e2max, e3max
    type(OrbitsIsospin), intent(in), target :: sps
    integer :: j12max, j12min, t12max, t12min, j12, t12
    integer :: i1, i2, i3
    type(SingleParticleOrbitIsospin), pointer :: o1, o2, o3
    integer :: nj, cnt, idx

    this%j = j
    this%p = p
    this%t = t
    this%n_state = n
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
      o1 => sps%GetOrbit(i1)
      do i2 = 1, i1
        o2 => sps%GetOrbit(i2)
        if(o1%e + o2%e > e2max) cycle
        do i3 = 1, i2
          o3 => sps%GetOrbit(i3)
          if(o1%e + o3%e > e2max) cycle
          if(o2%e + o3%e > e2max) cycle
          if(o1%e + o2%e + o3%e > e3max) cycle
          if((-1) ** (o1%l+o2%l+o3%l) /= p) cycle

          j12min = (o1%j+o2%j)/2
          j12max = abs(o1%j-o2%j)/2
          t12min = 1
          t12max = 0
          nj = 0
          do j12 = abs(o1%j-o2%j)/2, (o1%j+o2%j)/2
            if(triag(2*j12,o3%j,j)) cycle
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
      o1 => sps%orb(i1)
      do i2 = 1, i1
        o2 => sps%orb(i2)
        if(o1%e + o2%e > e2max) cycle
        do i3 = 1, i2
          o3 => sps%orb(i3)
          if(o1%e + o3%e > e2max) cycle
          if(o2%e + o3%e > e2max) cycle
          if(o1%e + o2%e + o3%e > e3max) cycle
          if((-1) ** (o1%l+o2%l+o3%l) /= p) cycle

          j12min = (o1%j+o2%j)/2
          j12max = abs(o1%j-o2%j)/2
          t12min = 1
          t12max = 0
          nj = 0
          do j12 = abs(o1%j-o2%j)/2, (o1%j+o2%j)/2
            if(triag(2*j12,o3%j,j)) cycle
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
    use MyLibrary, only: triag
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
    use MyLibrary, only: triag, sjs, hat
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

end module ThreeBodyInteraction
