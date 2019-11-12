module ThreeBodyModelSpace
  use omp_lib
  use SingleParticleState
  implicit none

  public :: ThreeBodySpace
  public :: ThreeBodyChannel

  ! Methods for orthonormal three-body space
  private :: FinThreeBodySpace
  private :: InitThreeBodySpace
  private :: FinThreeBodyChannel
  private :: InitThreeBodyChannel
  private :: FinAdditionalQN
  private :: InitAdditionalQN
  private :: CntDim
  private :: GenIter
  private :: Aop3
  private :: A3drct
  private :: A3exc1
  private :: A3exc2
  private :: FindIndexFromABCJ

  ! Methods for non-orthonormal three-body space
  private :: InitNonOrthIsospinThreeBodySpace
  private :: FinNonOrthIsospinThreeBodySpace
  private :: InitNonOrthIsospinThreeBodyChannel
  private :: FinNonOrthIsospinThreeBodyChannel
  private :: SetSortingIndices
  private :: Sort123
  private :: fin_sort_index
  private :: init_sort_index

  !
  ! Orthonormal three-body space
  ! |abci:JT>, here i is additional quantum number
  ! a, b, c: pn formalism
  type :: AdditionalQN
    integer :: north, nphys
    integer, allocatable :: idx2n(:)
    integer, allocatable :: n2spi1(:) ! permutations
    integer, allocatable :: n2spi2(:) ! permutations
    integer, allocatable :: n2spi3(:) ! permutations
    integer, allocatable :: n2J12(:)  ! Jab
    real(8), allocatable :: cfp(:,:)  ! (non-orth|orth)
  contains
    procedure :: init => InitAdditionalQN
    procedure :: fin => FinAdditionalQN
    procedure :: find => FindIndexFromABCJ
  end type AdditionalQN

  type :: ThreeBodyChannel
    integer :: j = -1
    integer :: p = 0
    integer :: z = 100
    integer :: n_state = 0
    integer :: n_cco_state = 0 ! hole-hole-particle     (not relevant for decoupling)
    integer :: n_cov_state = 0 ! hole-particle-valence  (not relevant for decoupling)
    integer :: n_coo_state = 0 ! hole-particle-particle (not relevant for decoupling)
    integer :: n_ccc_state = 0 ! hole-hole-hole
    integer :: n_ccv_state = 0 ! hole-hole-valence
    integer :: n_cvv_state = 0 ! hole-valence-valence
    integer :: n_vvv_state = 0 ! valence-valence-valence
    integer :: n_ovv_state = 0 ! particle-valence-valence
    integer :: n_oov_state = 0 ! particle-particle-valence
    integer :: n_ooo_state = 0 ! particle-particle-particle
    integer :: n_idx = 0
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
    type(Orbits), pointer :: sps
    integer, allocatable :: jpz2ch(:,:,:)
    integer :: NChan
    integer :: emax
    integer :: e2max
    integer :: e3max
  contains
    procedure :: init => InitThreeBodySpace
    procedure :: fin => FinThreeBodySpace
  end type ThreeBodySpace
  !
  ! end orthonormal three-body space
  !

  !
  ! Non-orthonomal three-body space
  ! |abcJ_{ab}T_{ab}:JT>
  ! a, b, c: isospin formalism
  type :: coef
    integer :: n = 0
    integer, allocatable :: idx2num(:)
#ifdef single_precision_three_body_file
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
  !
  ! end non-orthonomal three-body space
  !

  type, private :: iter3
    integer, allocatable :: ii1(:), ii2(:), ii3(:)
  contains
    procedure :: GenIter
  end type iter3
contains
  ! Methods for orthonormal three-body space
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
    type(Orbits), target, intent(in) :: sps
    integer, intent(in) :: e2max, e3max
    integer :: j, p, z, ich, nidx, n
    integer :: i1, i2, i3
    type(SingleParticleOrbit), pointer :: o1, o2, o3
    integer :: ni, nj
    integer, allocatable :: jj(:), pp(:), zz(:), nn(:), nnidx(:)
    type :: spis2n
      integer, allocatable :: spis2nidx(:,:,:)
    end type spis2n
    type(spis2n), allocatable :: ch(:)

    if(allocated(this%jpz)) call this%fin()

    this%sps => sps
    allocate(this%jpz2ch(1:2*min(e3max,3*sps%lmax)+3,-1:1,-3:3))
    this%jpz2ch(:,:,:) = 0
    ich = 0
    do j = 1, 2*min(e3max,3*sps%lmax)+3, 2
      do p = 1, -1, -2
        do z = -3, 3, 2
          n = 0
          nidx = 0
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
                if(o1%z + o2%z + o3%z /= z) cycle
                if((-1) ** (o1%l+o2%l+o3%l) /= p) cycle

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
            o1 => sps%GetOrbit(i1)
            do i2 = 1, i1
              o2 => sps%GetOrbit(i2)
              if(o1%e + o2%e > e2max) cycle
              do i3 = 1, i2
                o3 => sps%GetOrbit(i3)
                if(o1%e + o3%e > e2max) cycle
                if(o2%e + o3%e > e2max) cycle
                if(o1%e + o2%e + o3%e > e3max) cycle
                if(o1%z + o2%z + o3%z /= z) cycle
                if((-1) ** (o1%l+o2%l+o3%l) /= p) cycle

                call CntDim(sps,i1,i2,i3,j,ni,nj)
                n = n + ni
                if(ni /= 0) then
                  nidx = nidx + 1
                  ch(ich+1)%spis2nidx(i1,i2,i3) = nidx
                  ch(ich+1)%spis2nidx(i2,i3,i1) = nidx
                  ch(ich+1)%spis2nidx(i3,i1,i2) = nidx
                  ch(ich+1)%spis2nidx(i2,i1,i3) = nidx
                  ch(ich+1)%spis2nidx(i3,i2,i1) = nidx
                  ch(ich+1)%spis2nidx(i1,i3,i2) = nidx
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
    integer :: i1, i2, i3
    type(SingleParticleOrbit), pointer :: o1, o2, o3
    integer :: ni, nj, idx, i, cnt, loop, cnt_sub
    integer :: cvo1, cvo2, cvo3
    integer, allocatable :: nni(:), nnj(:)

    this%j = j
    this%p = p
    this%z = z
    this%n_state = n
    this%n_idx = nidx
    this%spis2idx = spis2idx

    allocate(this%idx(this%n_idx))
    allocate(nni(this%n_idx))
    allocate(nnj(this%n_idx))
    allocate(this%n2spi1(n))
    allocate(this%n2spi2(n))
    allocate(this%n2spi3(n))
    allocate(this%n2labl(n))

    cnt = 0
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
          if(o1%z + o2%z + o3%z /= z) cycle
          if((-1) ** (o1%l+o2%l+o3%l) /= p) cycle
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
    do loop = 1, 27
      cnt_sub = 0
      do i1 = 1, sps%norbs
        o1 => sps%GetOrbit(i1)
        do i2 = 1, i1
          o2 => sps%GetOrbit(i2)
          if(o1%e + o2%e > e2max) cycle
          do i3 = 1, i2
            o3 => sps%GetOrbit(i3)

            cvo1 = o1%GetCoreValenceOutside()
            cvo2 = o2%GetCoreValenceOutside()
            cvo3 = o3%GetCoreValenceOutside()

            if( loop == 1 .and. .not. (cvo1 == 0 .and. cvo2 == 0 .and. cvo3 == 2)) cycle ! cco
            if( loop == 2 .and. .not. (cvo1 == 0 .and. cvo2 == 2 .and. cvo3 == 0)) cycle ! coc
            if( loop == 3 .and. .not. (cvo1 == 2 .and. cvo2 == 0 .and. cvo3 == 0)) cycle ! occ

            if( loop == 4 .and. .not. (cvo1 == 0 .and. cvo2 == 2 .and. cvo3 == 1)) cycle ! cov
            if( loop == 5 .and. .not. (cvo1 == 2 .and. cvo2 == 1 .and. cvo3 == 0)) cycle ! ovc
            if( loop == 6 .and. .not. (cvo1 == 1 .and. cvo2 == 0 .and. cvo3 == 2)) cycle ! vco
            if( loop == 7 .and. .not. (cvo1 == 2 .and. cvo2 == 0 .and. cvo3 == 1)) cycle ! ocv
            if( loop == 8 .and. .not. (cvo1 == 1 .and. cvo2 == 2 .and. cvo3 == 0)) cycle ! voc
            if( loop == 9 .and. .not. (cvo1 == 0 .and. cvo2 == 1 .and. cvo3 == 2)) cycle ! cvo

            if( loop ==10 .and. .not. (cvo1 == 0 .and. cvo2 == 2 .and. cvo3 == 2)) cycle ! coo
            if( loop ==11 .and. .not. (cvo1 == 2 .and. cvo2 == 2 .and. cvo3 == 0)) cycle ! ooc
            if( loop ==12 .and. .not. (cvo1 == 2 .and. cvo2 == 0 .and. cvo3 == 2)) cycle ! oco

            if( loop ==13 .and. .not. (cvo1 == 0 .and. cvo2 == 0 .and. cvo3 == 0)) cycle ! ccc

            if( loop ==14 .and. .not. (cvo1 == 0 .and. cvo2 == 0 .and. cvo3 == 1)) cycle ! ccv
            if( loop ==15 .and. .not. (cvo1 == 1 .and. cvo2 == 0 .and. cvo3 == 0)) cycle ! vcc
            if( loop ==16 .and. .not. (cvo1 == 0 .and. cvo2 == 1 .and. cvo3 == 0)) cycle ! cvc

            if( loop ==17 .and. .not. (cvo1 == 0 .and. cvo2 == 1 .and. cvo3 == 1)) cycle ! cvv
            if( loop ==18 .and. .not. (cvo1 == 1 .and. cvo2 == 1 .and. cvo3 == 0)) cycle ! vvc
            if( loop ==19 .and. .not. (cvo1 == 1 .and. cvo2 == 0 .and. cvo3 == 1)) cycle ! vcv

            if( loop ==20 .and. .not. (cvo1 == 1 .and. cvo2 == 1 .and. cvo3 == 1)) cycle ! vvv

            if( loop ==21 .and. .not. (cvo1 == 2 .and. cvo2 == 1 .and. cvo3 == 1)) cycle ! ovv
            if( loop ==22 .and. .not. (cvo1 == 1 .and. cvo2 == 1 .and. cvo3 == 2)) cycle ! vvo
            if( loop ==23 .and. .not. (cvo1 == 1 .and. cvo2 == 2 .and. cvo3 == 1)) cycle ! vov

            if( loop ==24 .and. .not. (cvo1 == 2 .and. cvo2 == 2 .and. cvo3 == 1)) cycle ! oov
            if( loop ==25 .and. .not. (cvo1 == 1 .and. cvo2 == 2 .and. cvo3 == 2)) cycle ! voo
            if( loop ==26 .and. .not. (cvo1 == 2 .and. cvo2 == 1 .and. cvo3 == 2)) cycle ! ovo

            if( loop ==27 .and. .not. (cvo1 == 2 .and. cvo2 == 2 .and. cvo3 == 2)) cycle ! ooo
            if(o1%e + o3%e > e2max) cycle
            if(o2%e + o3%e > e2max) cycle
            if(o1%e + o2%e + o3%e > e3max) cycle
            if(o1%z + o2%z + o3%z /= z) cycle
            if((-1) ** (o1%l+o2%l+o3%l) /= p) cycle
            idx = this%spis2idx(i1,i2,i3)
            if(idx == 0) cycle
            ni = nni(idx)
            nj = nnj(idx)
            call this%idx(idx)%init(sps,i1,i2,i3,j,ni,nj)
            do i = 1, ni
              cnt = cnt + 1
              cnt_sub = cnt_sub + 1
              this%n2spi1(cnt) = i1
              this%n2spi2(cnt) = i2
              this%n2spi3(cnt) = i3
              this%n2labl(cnt) = i
              this%idx(idx)%idx2n(i) = cnt
            end do

          end do
        end do
      end do
      if(1 <= loop .and. loop <= 3) this%n_cco_state = this%n_cco_state + cnt_sub
      if(4 <= loop .and. loop <= 9) this%n_cov_state = this%n_cov_state + cnt_sub
      if(10<= loop .and. loop <=12) this%n_coo_state = this%n_coo_state + cnt_sub
      if(loop == 13) this%n_ccc_state = cnt_sub
      if(14<= loop .and. loop <=16) this%n_ccv_state = this%n_ccv_state + cnt_sub
      if(17<= loop .and. loop <=19) this%n_cvv_state = this%n_cvv_state + cnt_sub
      if(loop == 20) this%n_vvv_state = cnt_sub
      if(21<= loop .and. loop <=23) this%n_ovv_state = this%n_ovv_state + cnt_sub
      if(24<= loop .and. loop <=26) this%n_oov_state = this%n_oov_state + cnt_sub
      if(loop == 27) this%n_ooo_state = cnt_sub
    end do
#ifdef ModelSpaceDebug
    write(*,'(a,i3,a,i3,a,i3,a,i6)') "Three-body channel: J=", j, ", P=", p, ", Tz=", z, ", # of states=", this%n_state
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
    use MyLibrary, only: triag
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

  function FindIndexFromABCJ(this, a, b, c, Jab) result(r)
    class(AdditionalQN), intent(in) :: this
    integer, intent(in) :: a, b, c, Jab
    integer :: n
    integer :: r
    r = 0
    do n = 1, this%nphys
      if(a /= this%n2spi1(n)) cycle
      if(b /= this%n2spi2(n)) cycle
      if(c /= this%n2spi3(n)) cycle
      if(Jab/= this%n2J12(n)) cycle
      r = n
      return
    end do
  end function FindIndexFromABCJ

  subroutine CntDim(sps, i1, i2, i3, j, ni, nj)
    use MyLibrary, only: triag
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
    use MyLibrary, only: sjs
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
    use MyLibrary, only: sjs
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

  ! Methods for non-orthonormal three-body space
  subroutine FinNonOrthIsospinThreeBodySpace(this)
    class(NonOrthIsospinThreeBodySpace), intent(inout) :: this
    integer :: ich
    do ich = 1, this%NChan
      call this%jpt(ich)%fin()
    end do
    deallocate(this%jpt)
    deallocate(this%jpt2ch)
    this%sps => null()
    this%isps => null()
  end subroutine FinNonOrthIsospinThreeBodySpace

  subroutine InitNonOrthIsospinThreeBodySpace(this, sps, isps, e2max, e3max)
    use MyLibrary, only: triag
    class(NonOrthIsospinThreeBodySpace), intent(inout) :: this
    type(Orbits), target, intent(in) :: sps
    type(OrbitsIsospin), target, intent(in) :: isps
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

    if(allocated(this%jpt)) call this%fin()

    this%sps => sps
    this%isps => isps

    allocate(this%jpt2ch(1:2*min(e3max,3*isps%lmax)+3,-1:1,1:3))
    this%jpt2ch(:,:,:) = 0
    ich = 0
    do j = 1, 2*min(e3max,3*isps%lmax)+3, 2
      do p = 1, -1, -2
        do t = 1, 3, 2
          n = 0
          nidx = 0
          do i1 = 1, isps%norbs
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
      allocate(ch(ich)%spis2nidx(isps%norbs,isps%norbs,isps%norbs))
      ch(ich)%spis2nidx(:,:,:) = 0
    end do

    ich = 0
    do j = 1, 2*min(e3max,3*isps%lmax)+3, 2
      do p = 1, -1, -2
        do t = 1, 3, 2
          n = 0
          nidx = 0
          do i1 = 1, isps%norbs
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
              this%jt(j12,t12)%idx2num(n) = NOIAQN%JT2n(j12, t12)
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
                  this%jt(j12,t12)%idx2num(n) = NOIAQN%JT2n(jab, tab)
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
              this%jt(j12,t12)%idx2num(n) = NOIAQN%JT2n(j12, t12)
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
end module ThreeBodyModelSpace
