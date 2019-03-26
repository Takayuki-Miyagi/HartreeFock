module ThreeBodyModelSpace
  use omp_lib
  use SingleParticleState
  implicit none

  public :: ThreeBodySpace
  public :: ThreeBodyChannel

  private :: FinThreeBodySpace
  private :: InitThreeBodySpace
  private :: FinThreeBodyChannel
  private :: InitThreeBodyChannel
  private :: FinAdditionalQN
  private :: InitAdditionalQN
  private :: CntDim
  private :: GenIter
  private :: Aop3, A3drct, A3exc1, A3exc2

  type :: AdditionalQN
    integer :: north, nphys
    integer, allocatable :: idx2n(:)
    integer, allocatable :: n2spi1(:) ! permutations
    integer, allocatable :: n2spi2(:) ! permutations
    integer, allocatable :: n2spi3(:) ! permutations
    integer, allocatable :: n2J12(:)  ! Jab
    real(8), allocatable :: cfp(:,:)  ! (orth|non-orth)
  contains
    procedure :: init => InitAdditionalQN
    procedure :: fin => FinAdditionalQN
  end type AdditionalQN

  type :: ThreeBodyChannel
    integer :: j = -1
    integer :: p = 0
    integer :: z = 100
    integer :: n_state = 0
    integer :: n_hhp_state = 0 ! hole-hole-particle     (not relevant)
    integer :: n_hpv_state = 0 ! hole-particle-valence  (not relevant)
    integer :: n_hpp_state = 0 ! hole-particle-particle (not relevant)
    integer :: n_hhh_state = 0 ! hole-hole-hole
    integer :: n_hhv_state = 0 ! hole-hole-valence
    integer :: n_hvv_state = 0 ! hole-valence-valence
    integer :: n_vvv_state = 0 ! valence-valence-valence
    integer :: n_pvv_state = 0 ! particle-valence-valence
    integer :: n_ppv_state = 0 ! particle-particle-valence
    integer :: n_ppp_state = 0 ! particle-particle-particle
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
    integer, allocatable :: jpz2ch(:,:,:)
    integer :: NChan
    integer :: emax
    integer :: e2max
    integer :: e3max
  contains
    procedure :: init => InitThreeBodySpace
    procedure :: fin => FinThreeBodySpace
  end type ThreeBodySpace


  type, private :: iter3
    integer, allocatable :: ii1(:), ii2(:), ii3(:)
  contains
    procedure :: GenIter
  end type iter3
contains
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
    integer :: i1, i2, i3
    type(SingleParticleOrbit), pointer :: o1, o2, o3
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
    !integer :: i1, j1, l1, z1, e1
    !integer :: i2, j2, l2, z2, e2
    !integer :: i3, j3, l3, z3, e3
    integer :: ni, nj, idx, i, cnt, loop, cnt_sub
    integer, allocatable :: nni(:), nnj(:)

    this%j = j
    this%p = p
    this%z = z
    this%n_state = n
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
    do loop = 1, 27
      cnt_sub = 0
      do i1 = 1, sps%norbs
        o1 => sps%GetOrbit(i1)
        do i2 = 1, i1
          o2 => sps%GetOrbit(i2)
          if(o1%e + o2%e > e2max) cycle
          do i3 = 1, i2
            o3 => sps%GetOrbit(i3)
            if(.not. (loop == 1 .and. o1%ph == 0 .and. o2%ph == 0 .and. o3%ph == 1)) cycle ! hhp
            if(.not. (loop == 2 .and. o1%ph == 0 .and. o2%ph == 1 .and. o3%ph == 0)) cycle ! hph
            if(.not. (loop == 3 .and. o1%ph == 1 .and. o2%ph == 0 .and. o3%ph == 0)) cycle ! phh

            if(.not. (loop == 4 .and. o1%ph == 0 .and. o2%ph == 1 .and. o3%ph == 2)) cycle ! hpv
            if(.not. (loop == 5 .and. o1%ph == 1 .and. o2%ph == 2 .and. o3%ph == 0)) cycle ! pvh
            if(.not. (loop == 6 .and. o1%ph == 2 .and. o2%ph == 0 .and. o3%ph == 1)) cycle ! vhp
            if(.not. (loop == 7 .and. o1%ph == 1 .and. o2%ph == 0 .and. o3%ph == 2)) cycle ! phv
            if(.not. (loop == 8 .and. o1%ph == 2 .and. o2%ph == 1 .and. o3%ph == 0)) cycle ! vph
            if(.not. (loop == 9 .and. o1%ph == 0 .and. o2%ph == 2 .and. o3%ph == 1)) cycle ! hvp

            if(.not. (loop ==10 .and. o1%ph == 0 .and. o2%ph == 1 .and. o3%ph == 1)) cycle ! hpp
            if(.not. (loop ==11 .and. o1%ph == 1 .and. o2%ph == 1 .and. o3%ph == 0)) cycle ! pph
            if(.not. (loop ==12 .and. o1%ph == 1 .and. o2%ph == 0 .and. o3%ph == 1)) cycle ! php

            if(.not. (loop ==13 .and. o1%ph == 0 .and. o2%ph == 0 .and. o3%ph == 0)) cycle ! hhh

            if(.not. (loop ==14 .and. o1%ph == 0 .and. o2%ph == 0 .and. o3%ph == 2)) cycle ! hhv
            if(.not. (loop ==15 .and. o1%ph == 2 .and. o2%ph == 0 .and. o3%ph == 0)) cycle ! vhh
            if(.not. (loop ==16 .and. o1%ph == 0 .and. o2%ph == 2 .and. o3%ph == 0)) cycle ! hvh

            if(.not. (loop ==17 .and. o1%ph == 0 .and. o2%ph == 2 .and. o3%ph == 2)) cycle ! hvv
            if(.not. (loop ==18 .and. o1%ph == 2 .and. o2%ph == 2 .and. o3%ph == 0)) cycle ! vvh
            if(.not. (loop ==19 .and. o1%ph == 2 .and. o2%ph == 0 .and. o3%ph == 2)) cycle ! vhv

            if(.not. (loop ==20 .and. o1%ph == 2 .and. o2%ph == 2 .and. o3%ph == 2)) cycle ! vvv

            if(.not. (loop ==21 .and. o1%ph == 1 .and. o2%ph == 2 .and. o3%ph == 2)) cycle ! pvv
            if(.not. (loop ==22 .and. o1%ph == 2 .and. o2%ph == 2 .and. o3%ph == 1)) cycle ! vvp
            if(.not. (loop ==23 .and. o1%ph == 2 .and. o2%ph == 1 .and. o3%ph == 2)) cycle ! vpv

            if(.not. (loop ==24 .and. o1%ph == 1 .and. o2%ph == 1 .and. o3%ph == 2)) cycle ! ppv
            if(.not. (loop ==25 .and. o1%ph == 2 .and. o2%ph == 1 .and. o3%ph == 1)) cycle ! vpp
            if(.not. (loop ==26 .and. o1%ph == 1 .and. o2%ph == 2 .and. o3%ph == 1)) cycle ! pvp

            if(.not. (loop ==27 .and. o1%ph == 1 .and. o2%ph == 1 .and. o3%ph == 1)) cycle ! ppp

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

            if(1 <= loop .and. loop <= 3) this%n_hhp_state = this%n_hhp_state + cnt_sub
            if(4 <= loop .and. loop <= 9) this%n_hpv_state = this%n_hpv_state + cnt_sub
            if(10<= loop .and. loop <=12) this%n_hpp_state = this%n_hpp_state + cnt_sub
            if(loop == 13) this%n_hhh_state = cnt_sub
            if(14<= loop .and. loop <=16) this%n_hhv_state = this%n_hhv_state + cnt_sub
            if(17<= loop .and. loop <=19) this%n_hvv_state = this%n_hvv_state + cnt_sub
            if(loop == 20) this%n_vvv_state = cnt_sub
            if(21<= loop .and. loop <=23) this%n_pvv_state = this%n_pvv_state + cnt_sub
            if(24<= loop .and. loop <=26) this%n_ppv_state = this%n_ppv_state + cnt_sub
            if(loop == 27) this%n_ppp_state = cnt_sub
          end do
        end do
      end do
    end do

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
end module ThreeBodyModelSpace
