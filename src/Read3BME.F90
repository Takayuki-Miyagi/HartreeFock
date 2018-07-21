module read_3BME
  use InputParameters, only: parameters
  use MPIFunction, only: myrank
  use ModelSpace, only: spo_isospin
  implicit none
  private :: hat, triag
  type :: coef
    integer :: n = 0
    integer, allocatable :: idx2num(:)
#ifdef single_precision
    real(4), allocatable :: TrnsCoef(:)
#else
    real(8), allocatable :: TrnsCoef(:)
#endif
  end type coef

  type :: sorting
    integer :: idx, jmin, jmax, tmin, tmax
    type(coef), allocatable :: jt(:,:)
  contains
    procedure :: CalcCoef
  end type sorting

  type :: indices
    integer :: jmin, jmax, tmin, tmax
    integer, allocatable :: labels2n(:,:)
  contains
  end type indices

  type :: iThreeBodyChannel
    integer :: n, nsub
    type(indices), allocatable :: idx(:)
    type(sorting), allocatable :: idxsort(:)
    integer, allocatable :: labelssort(:,:,:)
    integer, allocatable :: labels2nsub(:,:,:)
#ifdef single_precision
    real(4), allocatable :: v3(:)
#else
    real(8), allocatable :: v3(:)
#endif
    real(8) :: usedmem = 0.d0
  contains
    procedure :: InitiThreeBodyChannel
    procedure :: FiniThreeBodyChannel
  end type iThreeBodyChannel

  type :: iThreeBodyScalar
    type(iThreeBodyChannel), allocatable :: jpt(:)
    integer, allocatable :: jpt2n(:,:,:)
    integer, allocatable :: j(:)
    integer, allocatable :: p(:)
    integer, allocatable :: t(:)
    integer :: n
    real(8) :: usedmem = 0.d0
  contains
    procedure :: init => InitiThreeBodyScalar
    procedure :: fin => FiniThreeBodyScalar
    procedure :: InitiThreeBodySpace
    procedure :: InitiReadScalar
    procedure :: Store3BME
  end type iThreeBodyScalar
#ifdef single_precision
  real(4), allocatable :: vtemp(:)
#else
  real(8), allocatable :: vtemp(:)
#endif

contains
  subroutine InitiThreeBodyScalar(this, sps, params, filename)
    class(iThreeBodyScalar), intent(inout) :: this
    type(spo_isospin), intent(in) :: sps
    type(parameters), intent(in) :: params
    character(len=*), intent(in) :: filename

    call this%InitiThreeBodySpace(sps, params)
    call this%InitiReadScalar(sps, params, filename)
  end subroutine InitiThreeBodyScalar

  subroutine FiniThreeBodyScalar(this, sps)
    class(iThreeBodyScalar), intent(inout) :: this
    type(spo_isospin), intent(in) :: sps
    integer :: ich
    do ich = 1, this%n
      call this%jpt(ich)%FiniThreeBodyChannel(sps)
    end do
    deallocate(this%j)
    deallocate(this%p)
    deallocate(this%t)
    deallocate(this%jpt2n)
  end subroutine FiniThreeBodyScalar

  subroutine InitiThreeBodySpace(this, sps, params)
    class(iThreeBodyScalar), intent(inout) :: this
    type(spo_isospin), intent(in) :: sps
    type(parameters), intent(in) :: params
    integer :: ich
    real(8) :: cnt
    integer :: e3max, n, m, nsp
    integer :: i1, j1, l1, j12, t12
    integer :: i2, j2, l2
    integer :: i3, j3, l3
    integer :: j, p, t, loop
    integer, allocatable :: labels2idx(:,:,:)

    e3max = params%e3max_3nf
    nsp = sps%n
    allocate(labels2idx(nsp, nsp, nsp))
    do loop = 1, 2
      ich = 0
      do j = 1, 2 * e3max + 3, 2
        do p = 1, -1, -2
          do t = 1, 3, 2

            labels2idx(:,:,:) = 0
            n = 0
            do i1 = 1, nsp
              j1 = sps%jj(i1)
              l1 = sps%ll(i1)
              do i2 = 1, i1
                j2 = sps%jj(i2)
                l2 = sps%ll(i2)
                do i3 = 1, i2
                  j3 = sps%jj(i3)
                  l3 = sps%ll(i3)
                  if((-1) ** (l1 + l2 + l3) /= p) cycle
                  if(sps%nshell(i1) + sps%nshell(i2) + sps%nshell(i3) > e3max) cycle
                  m = 0
                  do j12 = iabs(j1 - j2) / 2, (j1 + j2) / 2
                    do t12 = 0, 1
                      if(i1 == i2 .and. mod(j12 + t12, 2) == 0) cycle
                      if(triag(2 * t12,  1, t)) cycle
                      if(triag(2 * j12, j3, j)) cycle
                      m = m + 1
                    end do
                  end do
                  if(m /= 0) n = n + 1
                  if(loop == 2 .and. m /= 0) then
                    labels2idx(i1,i2,i3) = n
                  end if
                end do
              end do
            end do
            if(n /= 0) ich = ich + 1
            if(loop == 2 .and. n /= 0) then
              this%j(ich) = j
              this%p(ich) = p
              this%t(ich) = t
              this%jpt(ich)%nsub = n
              this%jpt2n(j, p, t) = ich
              this%jpt(ich)%labels2nsub(:,:,:) = labels2idx(:,:,:)

            end if
          end do
        end do
      end do
      if(loop == 1) then
        this%n = ich
        allocate(this%jpt(ich))
        allocate(this%j(ich))
        allocate(this%p(ich))
        allocate(this%t(ich))
        allocate(this%jpt2n(1:2*e3max+3,-1:1,1:3))
        this%j(:) = 0
        this%p(:) = 0
        this%t(:) = 0
        this%jpt2n(:,:,:) = 0
        do ich = 1, this%n
          allocate(this%jpt(ich)%labels2nsub(nsp,nsp,nsp))
          this%jpt(ich)%labels2nsub(:,:,:) = 0
        end do
      end if
    end do
    deallocate(labels2idx)

    cnt = 0.d0
    do ich = 1, this%n
      call this%jpt(ich)%InitiThreeBodyChannel(params, sps, this%j(ich), this%p(ich), this%t(ich))
      cnt = cnt + dble(this%jpt(ich)%n) * (dble(this%jpt(ich)%n) + 1) / 2
      this%usedmem = this%usedmem + this%jpt(ich)%usedmem
    end do
#ifdef single_precision
    this%usedmem = this%usedmem + cnt * 4.d0 / (1024.d0 ** 3)
#else
    this%usedmem = this%usedmem + cnt * 8.d0 / (1024.d0 ** 3)
#endif
    if(myrank == 0) then
      write(*,'(a, f9.4, a)') &
      & 'Estimated Memory for storing isospin-sym 3BME: ', &
      &  this%usedmem, " GB"
    end if
  end subroutine InitiThreeBodySpace

  subroutine InitiThreeBodyChannel(this, params, sps, j, p, t)
    class(iThreeBodyChannel), intent(inout) :: this
    type(spo_isospin), intent(in) :: sps
    type(parameters), intent(in) :: params
    integer, intent(in) :: j, p, t
    integer :: m, n, idx, loop, idxsort
    integer :: jmin, jmax, tmin, tmax
    integer :: i1, j1, l1, j12, t12
    integer :: i2, j2, l2
    integer :: i3, j3, l3
    integer :: a, b, c, sort
    real(8) :: cnt
    allocate(this%idx(this%nsub))
    n = 0
    m = sps%n
    do i1 = 1, m
      j1 = sps%jj(i1); l1 = sps%ll(i1)
      do i2 = 1, i1
        j2 = sps%jj(i2); l2 = sps%ll(i2)
        do i3 = 1, i2
          j3 = sps%jj(i3); l3 = sps%ll(i3)
          if((-1) ** (l1 + l2 + l3) /= p) cycle
          idx = this%labels2nsub(i1,i2,i3)
          if(idx == 0) cycle
          do loop = 1, 2
            jmin = 100; jmax = -100
            tmin = 100; tmax = -100
            do j12 = iabs(j1 - j2) / 2, (j1 + j2) / 2
              do t12 = 0, 1
                if(i1 == i2 .and. mod(j12 + t12, 2) == 0) cycle
                if(triag(2 * t12,  1, t)) cycle
                if(triag(2 * j12, j3, j)) cycle
                if(t12 > tmax) tmax = t12
                if(t12 < tmin) tmin = t12
                if(j12 > jmax) jmax = j12
                if(j12 < jmin) jmin = j12
                if(loop == 2) then
                  n = n + 1
                  this%idx(idx)%labels2n(j12, t12) = n
                end if
              end do
            end do
            if(loop == 1) then
              this%idx(idx)%jmin = jmin
              this%idx(idx)%jmax = jmax
              this%idx(idx)%tmin = tmin
              this%idx(idx)%tmax = tmax
              allocate(this%idx(idx)%labels2n(jmin:jmax, tmin:tmax))
            end if
          end do
        end do
      end do
    end do
    this%n = n
    allocate(this%v3(n*(n+1)/2))

    allocate(this%labelssort(m,m,m))
    this%labelssort(:,:,:) = 0
    idxsort = 0
    do i1 = 1, m
      j1 = sps%jj(i1)
      l1 = sps%ll(i1)
      do i2 = 1, m
        j2 = sps%jj(i2)
        l2 = sps%ll(i2)
        do i3 = 1, m
          j3 = sps%jj(i3)
          l3 = sps%ll(i3)
          if((-1) ** (l1 + l2 + l3) /= p) cycle
          if(sps%nshell(i1) + sps%nshell(i2) + sps%nshell(i3) > params%e3cut) cycle
          idxsort = idxsort + 1
          this%labelssort(i1,i2,i3) = idxsort
        end do
      end do
    end do
    allocate(this%idxsort(idxsort))
    cnt = 0.d0
    do i1 = 1, m
      j1 = sps%jj(i1)
      l1 = sps%ll(i1)
      do i2 = 1, m
        j2 = sps%jj(i2)
        l2 = sps%ll(i2)
        do i3 = 1, m
          j3 = sps%jj(i3)
          l3 = sps%ll(i3)
          if((-1) ** (l1 + l2 + l3) /= p) cycle
          if(sps%nshell(i1) + sps%nshell(i2) + sps%nshell(i3) > params%e3cut) cycle
          call SortABC(i1,i2,i3,a,b,c,sort)
          idx = this%labels2nsub(a,b,c)
          idxsort = this%labelssort(i1,i2,i3)
          this%idxsort(idxsort)%idx = idx
          if(idx < 1) cycle
          jmin = 100; jmax = -100
          tmin = 100; tmax = -100
          do j12 = iabs(j1 - j2) / 2, (j1 + j2) / 2
            do t12 = 0, 1
              if(i1 == i2 .and. mod(j12 + t12, 2) == 0) cycle
              if(triag(2 * t12,  1, t)) cycle
              if(triag(2 * j12, j3, j)) cycle
              if(t12 > tmax) tmax = t12
              if(t12 < tmin) tmin = t12
              if(j12 > jmax) jmax = j12
              if(j12 < jmin) jmin = j12
            end do
          end do
          this%idxsort(idxsort)%jmin = jmin
          this%idxsort(idxsort)%jmax = jmax
          this%idxsort(idxsort)%tmin = tmin
          this%idxsort(idxsort)%tmax = tmax
          allocate(this%idxsort(idxsort)%jt(jmin:jmax,tmin:tmax))
          call this%idxsort(idxsort)%CalcCoef(sps,i1,i2,i3,a,b,c,sort,j,t,this%idx(idx))
          do j12 = jmin, jmax
            do t12 = tmin, tmax
              cnt = cnt + this%idxsort(idxsort)%jt(j12,t12)%n
            end do
          end do
        end do
      end do
    end do
#ifdef single_precision
    this%usedmem = cnt * 8.d0 / (1024.d0 ** 3)
#else
    this%usedmem = cnt * 12.d0 / (1024.d0 ** 3)
#endif
  end subroutine InitiThreeBodyChannel

  subroutine CalcCoef(this, sps, i1, i2, i3, a, b, c, sort, j, t, idxt)
    use common_library, only: sjs
    class(sorting), intent(inout) :: this
    type(spo_isospin), intent(in) :: sps
    type(indices), intent(in) :: idxt
    integer, intent(in) :: i1, i2, i3, a, b, c, sort, j, t
    integer :: n, loop
    integer :: j1, j2, j3, ja, jb, jc
    integer :: j12, t12, jab, tab
    j1 = sps%jj(i1)
    j2 = sps%jj(i2)
    j3 = sps%jj(i3)
    ja = sps%jj(a)
    jb = sps%jj(b)
    jc = sps%jj(c)
    do j12 = iabs(j1 - j2) / 2, (j1 + j2)/2
      do t12 = 0, 1
        if(i1 == i2 .and. mod(j12 + t12, 2) == 0) cycle
        if(triag(2 * t12,  1, t)) cycle
        if(triag(2 * j12, j3, j)) cycle
        do loop = 1, 2
          if(loop == 1) then
            n = 0
            if(sort == 1 .or. sort == 4) then
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
            select case(sort)
            case(1)
              n = 1
#ifdef single_precision
              this%jt(j12,t12)%TrnsCoef(n) = 1.0
#else
              this%jt(j12,t12)%TrnsCoef(n) = 1.d0
#endif
              this%jt(j12,t12)%idx2num(n) = idxt%labels2n(j12, t12)
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
                  this%jt(j12,t12)%idx2num(n) = idxt%labels2n(jab, tab)
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
                  this%jt(j12,t12)%idx2num(n) = idxt%labels2n(jab, tab)
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
              this%jt(j12,t12)%idx2num(n) = idxt%labels2n(j12, t12)
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
                  this%jt(j12,t12)%idx2num(n) = idxt%labels2n(jab, tab)
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
                  this%jt(j12,t12)%idx2num(n) = idxt%labels2n(jab, tab)
                end do
              end do

            end select
          end if
        end do

      end do
    end do

  end subroutine CalcCoef

  subroutine FiniThreeBodyChannel(this, sps)
    class(iThreeBodyChannel), intent(inout) :: this
    type(spo_isospin), intent(in) :: sps
    integer :: i1, i2, i3, m, idx, j12, t12
    integer :: idxs
    m = sps%n
    do idx = 1, this%nsub
      deallocate(this%idx(idx)%labels2n)
    end do

    do i1 = 1, m
      do i2 = 1, m
        do i3 = 1, m
          idxs = this%labelssort(i1,i2,i3)
          if(idxs < 1) cycle
          idx = this%idxsort(idxs)%idx
          if(idx < 1) cycle
          do j12 = this%idxsort(idxs)%jmin, &
                & this%idxsort(idxs)%jmax
            do t12 = this%idxsort(idxs)%tmin, &
                  & this%idxsort(idxs)%tmax
              if(this%idxsort(idxs)%jt(j12,t12)%n < 1) cycle
              deallocate(this%idxsort(idxs)%jt(j12,t12)%idx2num)
              deallocate(this%idxsort(idxs)%jt(j12,t12)%TrnsCoef)
            end do
          end do
          deallocate(this%idxsort(idxs)%jt)
        end do
      end do
    end do

    deallocate(this%idx)
    deallocate(this%labels2nsub)
    deallocate(this%labelssort)
    deallocate(this%v3)
  end subroutine FiniThreeBodyChannel

  subroutine InitiReadScalar(this, sps, params, filename)
    use class_sys, only: sy
    class(iThreeBodyScalar), intent(inout) :: this
    type(spo_isospin), intent(in) :: sps
    type(parameters), intent(in) :: params
    character(len=*), intent(in) :: filename
    type(sy) :: sys
    integer :: iunit = 20
    integer :: i
    integer(8) :: num
    real(8) :: st
    logical :: ex
    inquire(file = filename, exist = ex)
    if(.not. ex) then
      write(*, '(a, " does not exist.")') trim(filename)
      stop
    end if

    call Cnt3BME(sps, params, num)
#ifdef single_precision
    st = dble(num) * 4.d0 / (1024.d0 ** 3)
#else
    st = dble(num) * 8.d0 / (1024.d0 ** 3)
#endif
    allocate(vtemp(num))

    if(sys%find(filename, '.txt')) then
      !call read_3bme_linebyline_txt(sps, params, filename)
      open(iunit, file = filename, status = 'old', access='stream', form='formatted')
      !read(iunit, *)  i ! old version
      read(iunit, *) vtemp
      close(iunit)
    elseif(sys%find(filename, '.bin')) then
      !call read_3bme_linebyline(sps, params, filename)
      open(iunit, form = 'unformatted', file = filename, status = 'old', access='stream')
      !open(iunit, form = 'unformatted', file = filename, status = 'old') ! old version
      !read(iunit)  i ! old version
      read(iunit) vtemp
      close(iunit)
    else
      write(*,'(2a)') 'The file format cannot be detected: ', trim(filename)
      stop
    end if

    if(myrank == 0) then
      write(*,'(a, f9.4, a)') &
          & 'Estimated Memory for temporal 3N file: ', &
          &  st, " GB"
    end if

    call this%Store3bme(sps, params)
    deallocate(vtemp)
  end subroutine InitiReadScalar

  function Get3BME(this, a, b, c, jab, tab, d, e, f, jde, tde) result(v)
#ifdef single_precision
    real(4) :: v
#else
    real(8) :: v
#endif
    type(iThreeBodyChannel), intent(in) :: this
    integer, intent(in) :: a, b, c, d, e, f
    integer, intent(in) :: jab, tab, jde, tde
    integer :: i, j, bra, ket
    integer :: idsb, idsk
#ifdef single_precision
    v = 0.0
#else
    v = 0.d0
#endif
    idsb = this%labelssort(a,b,c)
    idsk = this%labelssort(d,e,f)
    if(idsb < 1) return
    if(idsk < 1) return
    if(a == b .and. mod(jab + tab, 2) == 0) return
    if(d == e .and. mod(jde + tde, 2) == 0) return
    if(this%idxsort(idsb)%idx < 1) return
    if(this%idxsort(idsk)%idx < 1) return
    do i = 1, this%idxsort(idsb)%jt(jab,tab)%n
      bra = this%idxsort(idsb)%jt(jab,tab)%idx2num(i)
      do j = 1, this%idxsort(idsk)%jt(jde,tde)%n
        ket = this%idxsort(idsk)%jt(jde,tde)%idx2num(j)
        v = v + this%v3(max(bra,ket) * (max(bra,ket)-1) / 2 + min(bra,ket)) * &
            & this%idxsort(idsb)%jt(jab,tab)%TrnsCoef(i) * &
            & this%idxsort(idsk)%jt(jde,tde)%TrnsCoef(j)
      end do
    end do
  end function Get3BME

  subroutine SortABC(i1,i2,i3,ii1,ii2,ii3,n_recouple)
    integer, intent(in) :: i1, i2, i3
    integer, intent(out) :: ii1, ii2, ii3, n_recouple
    if(i1 .ge. i2) then
      if(i2 .ge. i3) then
        ! i1 >= i2 >= i3
        ii1 = i1; ii2 = i2; ii3 = i3
        n_recouple = 1
      elseif(i2 .lt. i3) then
        if(i1 .ge. i3) then
          ! i1 >= i3 > i2
          ii1 = i1; ii2 = i3; ii3 = i2
          n_recouple = 2
        elseif(i1 .lt. i3) then
          ! i3 > i1 >= i2
          ii1 = i3; ii2 = i1; ii3 = i2
          n_recouple = 3
        end if
      end if

    elseif(i1 .lt. i2) then
      if(i1 .ge. i3) then
        ! i2 > i1 >= i3
        ii1 = i2; ii2 = i1; ii3 = i3
        n_recouple = 4
      elseif(i1 .lt. i3) then

        if(i2 .ge. i3) then
          ! i2 >= i3 > i1
          ii1 = i2; ii2 = i3; ii3 = i1
          n_recouple = 5
        elseif(i2 .lt. i3) then
          ! i3 > i2 > i1
          ii1 = i3; ii2 = i2; ii3 = i1
          n_recouple = 6
        end if
      end if
    end if
  end subroutine SortABC

  subroutine Cnt3BME(sps, params, numtot)
    type(spo_isospin), intent(in) :: sps
    type(parameters), intent(in) :: params
    integer(8), intent(out) :: numtot
    integer :: i1, i2, i3, i4, i5, i6
    integer :: l1, l2, l3, l4, l5, l6
    integer :: j1, j2, j3, j4, j5, j6
    integer :: i5max, i6max
    integer :: p123, p456
    integer :: j12, j45, t12, t45
    integer :: j, t
    integer :: total_num
    integer :: emax, e2max, e3max
    emax = params%emax_3nf
    e2max = params%e2max_3nf
    e3max = params%e3max_3nf
    total_num = 0
    do i1 = 1, sps%n
      l1 = sps%ll(i1); j1 = sps%jj(i1)
      if(sps%nshell(i1) > emax) cycle
      do i2 = 1, i1
        l2 = sps%ll(i2); j2 = sps%jj(i2)
        if(sps%nshell(i2) > emax) cycle
        if(sps%nshell(i1) + sps%nshell(i2) > e2max) cycle
        do i3 = 1, i2
          l3 = sps%ll(i3); j3 = sps%jj(i3)
          if(sps%nshell(i3) > e3max) cycle
          if(sps%nshell(i2) + sps%nshell(i3) > e2max) cycle
          if(sps%nshell(i1) + sps%nshell(i3) > e2max) cycle
          if(sps%nshell(i1) + sps%nshell(i2) + sps%nshell(i3) > e3max) cycle
          p123 = (-1) ** (l1 + l2 + l3)

          do i4 = 1, i1
            l4 = sps%ll(i4); j4 = sps%jj(i4)
            if(sps%nshell(i4) > emax) cycle
            if(i1 == i4) then
              i5max = i2
            else
              i5max = i4
            end if
            do i5 = 1, i5max
              l5 = sps%ll(i5); j5 = sps%jj(i5)
              if(sps%nshell(i5) > emax) cycle
              if(sps%nshell(i4) + sps%nshell(i5) > e2max) cycle
              if(i1 == i4 .and. i2 == i5) then
                i6max = i3
              else
                i6max = i5
              end if
              do i6 = 1, i6max
                l6 = sps%ll(i6); j6 = sps%jj(i6)
                if(sps%nshell(i6) > emax) cycle
                if(sps%nshell(i5) + sps%nshell(i6) > e2max) cycle
                if(sps%nshell(i4) + sps%nshell(i6) > e2max) cycle
                if(sps%nshell(i4) + sps%nshell(i5) + sps%nshell(i6) > e3max) cycle
                p456 = (-1) ** (l4 + l5 +l6)
                if(p123 /= p456) cycle

                do j12 = iabs(j1 - j2)/2, (j1 + j2)/2
                  do j45 = iabs(j4 - j5)/2, (j4 + j5)/2
                    do j = max(iabs(2*j12 - j3), iabs(2*j45 - j6)), &
                          & min(2*j12 + j3, 2*j45 + j6), 2
                      do t12 = 0, 1
                        do t45 = 0, 1
                          do t = 1, min(2*t12 + 1, 2*t45 + 1), 2
                            total_num = total_num + 1
                          end do
                        end do
                      end do
                    end do
                  end do
                end do
              end do
            end do
          end do
        end do
      end do
    end do
    numtot = total_num
  end subroutine Cnt3BME

  subroutine read_3bme_linebyline_txt(sps, params, f)
    type(spo_isospin), intent(in) :: sps
    type(parameters), intent(in) :: params
    character(len=*), intent(in) :: f
    integer :: i1, i2, i3, i4, i5, i6
    integer :: l1, l2, l3, l4, l5, l6
    integer :: j1, j2, j3, j4, j5, j6
    integer :: i5max, i6max
    integer :: p123, p456
    integer :: j12, j45, t12, t45
    integer :: j, t
    integer :: total_num
    integer :: emax, e2max, e3max
    integer :: iunit = 20
    emax = params%emax_3nf
    e2max = params%e2max_3nf
    e3max = params%e3max_3nf
    open(iunit, file=f, form = 'formatted')
    total_num = 0
    do i1 = 1, sps%n
      l1 = sps%ll(i1); j1 = sps%jj(i1)
      if(sps%nshell(i1) > emax) cycle
      do i2 = 1, i1
        l2 = sps%ll(i2); j2 = sps%jj(i2)
        if(sps%nshell(i2) > emax) cycle
        if(sps%nshell(i1) + sps%nshell(i2) > e2max) cycle
        do i3 = 1, i2
          l3 = sps%ll(i3); j3 = sps%jj(i3)
          if(sps%nshell(i3) > e3max) cycle
          if(sps%nshell(i2) + sps%nshell(i3) > e2max) cycle
          if(sps%nshell(i1) + sps%nshell(i3) > e2max) cycle
          if(sps%nshell(i1) + sps%nshell(i2) + sps%nshell(i3) > e3max) cycle
          p123 = (-1) ** (l1 + l2 + l3)

          do i4 = 1, i1
            l4 = sps%ll(i4); j4 = sps%jj(i4)
            if(sps%nshell(i4) > emax) cycle
            if(i1 == i4) then
              i5max = i2
            else
              i5max = i4
            end if
            do i5 = 1, i5max
              l5 = sps%ll(i5); j5 = sps%jj(i5)
              if(sps%nshell(i5) > emax) cycle
              if(sps%nshell(i4) + sps%nshell(i5) > e2max) cycle
              if(i1 == i4 .and. i2 == i5) then
                i6max = i3
              else
                i6max = i5
              end if
              do i6 = 1, i6max
                l6 = sps%ll(i6); j6 = sps%jj(i6)
                if(sps%nshell(i6) > emax) cycle
                if(sps%nshell(i5) + sps%nshell(i6) > e2max) cycle
                if(sps%nshell(i4) + sps%nshell(i6) > e2max) cycle
                if(sps%nshell(i4) + sps%nshell(i5) + sps%nshell(i6) > e3max) cycle
                p456 = (-1) ** (l4 + l5 +l6)
                if(p123 /= p456) cycle

                do j12 = iabs(j1 - j2)/2, (j1 + j2)/2
                  do j45 = iabs(j4 - j5)/2, (j4 + j5)/2
                    do j = max(iabs(2*j12 - j3), iabs(2*j45 - j6)), &
                          & min(2*j12 + j3, 2*j45 + j6), 2
                      do t12 = 0, 1
                        do t45 = 0, 1
                          do t = 1, min(2*t12 + 1, 2*t45 + 1), 2
                            total_num = total_num + 1
                            read(iunit,*) vtemp(total_num)
                          end do
                        end do
                      end do
                    end do
                  end do
                end do
              end do
            end do
          end do
        end do
      end do
    end do
    close(iunit)
  end subroutine read_3bme_linebyline_txt

  subroutine read_3bme_linebyline(sps, params, f)
    type(spo_isospin), intent(in) :: sps
    type(parameters), intent(in) :: params
    character(len=*), intent(in) :: f
    integer :: i1, i2, i3, i4, i5, i6
    integer :: l1, l2, l3, l4, l5, l6
    integer :: j1, j2, j3, j4, j5, j6
    integer :: i5max, i6max
    integer :: p123, p456
    integer :: j12, j45, t12, t45
    integer :: j, t
    integer :: total_num
    integer :: emax, e2max, e3max
    integer :: iunit = 20
    emax = params%emax_3nf
    e2max = params%e2max_3nf
    e3max = params%e3max_3nf
    open(iunit, file=f, form = 'unformatted')
    total_num = 0
    do i1 = 1, sps%n
      l1 = sps%ll(i1); j1 = sps%jj(i1)
      if(sps%nshell(i1) > emax) cycle
      do i2 = 1, i1
        l2 = sps%ll(i2); j2 = sps%jj(i2)
        if(sps%nshell(i2) > emax) cycle
        if(sps%nshell(i1) + sps%nshell(i2) > e2max) cycle
        do i3 = 1, i2
          l3 = sps%ll(i3); j3 = sps%jj(i3)
          if(sps%nshell(i3) > e3max) cycle
          if(sps%nshell(i2) + sps%nshell(i3) > e2max) cycle
          if(sps%nshell(i1) + sps%nshell(i3) > e2max) cycle
          if(sps%nshell(i1) + sps%nshell(i2) + sps%nshell(i3) > e3max) cycle
          p123 = (-1) ** (l1 + l2 + l3)

          do i4 = 1, i1
            l4 = sps%ll(i4); j4 = sps%jj(i4)
            if(sps%nshell(i4) > emax) cycle
            if(i1 == i4) then
              i5max = i2
            else
              i5max = i4
            end if
            do i5 = 1, i5max
              l5 = sps%ll(i5); j5 = sps%jj(i5)
              if(sps%nshell(i5) > emax) cycle
              if(sps%nshell(i4) + sps%nshell(i5) > e2max) cycle
              if(i1 == i4 .and. i2 == i5) then
                i6max = i3
              else
                i6max = i5
              end if
              do i6 = 1, i6max
                l6 = sps%ll(i6); j6 = sps%jj(i6)
                if(sps%nshell(i6) > emax) cycle
                if(sps%nshell(i5) + sps%nshell(i6) > e2max) cycle
                if(sps%nshell(i4) + sps%nshell(i6) > e2max) cycle
                if(sps%nshell(i4) + sps%nshell(i5) + sps%nshell(i6) > e3max) cycle
                p456 = (-1) ** (l4 + l5 +l6)
                if(p123 /= p456) cycle

                do j12 = iabs(j1 - j2)/2, (j1 + j2)/2
                  do j45 = iabs(j4 - j5)/2, (j4 + j5)/2
                    do j = max(iabs(2*j12 - j3), iabs(2*j45 - j6)), &
                          & min(2*j12 + j3, 2*j45 + j6), 2
                      do t12 = 0, 1
                        do t45 = 0, 1
                          do t = 1, min(2*t12 + 1, 2*t45 + 1), 2
                            total_num = total_num + 1
                            read(iunit) vtemp(total_num)
                          end do
                        end do
                      end do
                    end do
                  end do
                end do
              end do
            end do
          end do
        end do
      end do
    end do
    close(iunit)
  end subroutine read_3bme_linebyline

  subroutine Store3BME(this, sps, params)
    class(iThreeBodyScalar), intent(inout) :: this
    type(spo_isospin), intent(in) :: sps
    type(parameters), intent(in) :: params
    integer :: i1, i2, i3, i4, i5, i6
    integer :: l1, l2, l3, l4, l5, l6
    integer :: j1, j2, j3, j4, j5, j6
    integer :: i5max, i6max
    integer :: p123, p456
    integer :: j12, j45, t12, t45
    integer :: j, t, ich
    integer :: idxb, idxk, bra, ket
    integer(8) :: total_num
    integer :: emax, e2max, e3max, e3cut
    emax = params%emax_3nf
    e2max = params%e2max_3nf
    e3max = params%e3max_3nf
    e3cut = params%e3cut

    total_num = 0
    do i1 = 1, sps%n
      l1 = sps%ll(i1); j1 = sps%jj(i1)
      if(sps%nshell(i1) > emax) cycle
      do i2 = 1, i1
        l2 = sps%ll(i2); j2 = sps%jj(i2)
        if(sps%nshell(i2) > emax) cycle
        if(sps%nshell(i1) + sps%nshell(i2) > e2max) cycle
        do i3 = 1, i2
          l3 = sps%ll(i3); j3 = sps%jj(i3)
          if(sps%nshell(i3) > e3max) cycle
          if(sps%nshell(i2) + sps%nshell(i3) > e2max) cycle
          if(sps%nshell(i1) + sps%nshell(i3) > e2max) cycle
          if(sps%nshell(i1) + sps%nshell(i2) + sps%nshell(i3) > e3max) cycle
          p123 = (-1) ** (l1 + l2 + l3)

          do i4 = 1, i1
            l4 = sps%ll(i4); j4 = sps%jj(i4)
            if(sps%nshell(i4) > emax) cycle
            if(i1 == i4) then
              i5max = i2
            else
              i5max = i4
            end if
            do i5 = 1, i5max
              l5 = sps%ll(i5); j5 = sps%jj(i5)
              if(sps%nshell(i5) > emax) cycle
              if(sps%nshell(i4) + sps%nshell(i5) > e2max) cycle
              if(i1 == i4 .and. i2 == i5) then
                i6max = i3
              else
                i6max = i5
              end if
              do i6 = 1, i6max
                l6 = sps%ll(i6); j6 = sps%jj(i6)
                if(sps%nshell(i6) > emax) cycle
                if(sps%nshell(i5) + sps%nshell(i6) > e2max) cycle
                if(sps%nshell(i4) + sps%nshell(i6) > e2max) cycle
                if(sps%nshell(i4) + sps%nshell(i5) + sps%nshell(i6) > e3max) cycle
                p456 = (-1) ** (l4 + l5 +l6)
                if(p123 /= p456) cycle

                do j12 = iabs(j1 - j2)/2, (j1 + j2)/2
                  do j45 = iabs(j4 - j5)/2, (j4 + j5)/2
                    do j = max(iabs(2*j12 - j3), iabs(2*j45 - j6)), &
                          & min(2*j12 + j3, 2*j45 + j6), 2
                      do t12 = 0, 1
                        do t45 = 0, 1
                          do t = 1, min(2*t12 + 1, 2*t45 + 1), 2
                            total_num = total_num + 1
                            if(sps%nshell(i1) + sps%nshell(i2) + sps%nshell(i3) > e3cut) cycle
                            if(sps%nshell(i4) + sps%nshell(i5) + sps%nshell(i6) > e3cut) cycle
                            if(i1 == i2 .and. mod(j12 + t12, 2) == 0) cycle
                            if(i4 == i5 .and. mod(j45 + t45, 2) == 0) cycle
                            ich = this%jpt2n(j, p123, t)
                            idxb = this%jpt(ich)%labels2nsub(i1,i2,i3)
                            idxk = this%jpt(ich)%labels2nsub(i4,i5,i6)
                            if(idxb * idxk == 0) cycle
                            bra = this%jpt(ich)%idx(idxb)%labels2n(j12,t12)
                            ket = this%jpt(ich)%idx(idxk)%labels2n(j45,t45)
                            if(bra * ket == 0) cycle
                            this%jpt(ich)%v3(max(bra,ket) * (max(bra,ket)-1)/2 + min(bra,ket)) = vtemp(total_num)
                          end do
                        end do
                      end do
                    end do
                  end do
                end do
              end do
            end do
          end do
        end do
      end do
    end do
  end subroutine Store3BME

  real(8) function hat(i)
    integer, intent(in) :: i
    hat = dsqrt(dble(i + 1))
  end function hat

  logical function triag(i,j,k)
    implicit none
    integer,intent(in)::i,j,k
    triag = ((i-(j+k))*(i-abs(j-k)) > 0)
  end function triag
end module read_3BME
