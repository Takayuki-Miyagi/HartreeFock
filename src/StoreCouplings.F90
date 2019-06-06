module StoreCouplings
  use omp_lib
  use Profiler, only: timer
  implicit none

  public :: CGsStore
  public :: SixJsStore
  public :: NineJsStore
  public :: TMbracketStore

  ! methods for CG coefficient

  ! methods for 6-j symbol
  private :: InitSixJsStore
  private :: FinSixJsStore
  private :: GetStoredSixJ
  private :: InitSixJs_intermediate
  private :: FinSixJs_intermediate
  private :: InitSixJs_final
  private :: FinSixJs_final

  ! methods for 9-j symbol
  private :: InitNineJsStore
  private :: FinNineJsStore
  private :: GetStoredNineJ
  private :: InitNineJs_intermediate
  private :: FinNineJs_intermediate
  private :: InitNineJs_final
  private :: FinNineJs_final

  ! methods for Talmi-Moshinsky bracket
  private :: InitTMbracketStore
  private :: FinTMbracketStore
  private :: GetStoredTMbracket
  private :: InitTMbk_intermediate
  private :: FinTMbk_intermediate

  ! CG coefficient storing
  !
  ! CG%v(j1i,j2i)%v(m1i,m2i,j12i)
  !   (  j1  j2 | j12 }
  ! = (  m1  m2 | m12 }
  !
  !   j1i = j1       (j1: integer)
  !   j1i = (j1+1)/2 (j1: half integer)
  type :: CGs_final
    real(8), allocatable :: v(:,:,:)
  contains
    procedure :: InitCGs_final
    procedure :: FinCGs_final

    generic :: init => InitCGs_final
    generic :: fin => FinCGs_final
  end type CGs_final

  type :: CGsStore
    type(CGs_final), allocatable :: v(:,:)
    integer :: j1min, j1max
    integer :: j2min, j2max
    logical :: half_j1
    logical :: half_j2
    logical :: half_j12
  contains
    procedure :: InitCGsStore
    procedure :: FinCGsStore
    procedure :: GetStoredCG

    generic :: init => InitCGsStore
    generic :: fin => FinCGsStore
    generic :: get => GetStoredCG
  end type CGsStore

  ! 6-j symbol storing
  !
  ! sixj%v(j1i,j2i,j3i)%v(j12i,j23i)%v(ji)
  !   {  j1  j2 j12 }
  ! = {  j3   J j23 }
  !
  !   j1i = j1       (j1: integer)
  !   j1i = (j1+1)/2 (j1: half integer)
  type :: SixJs_final
    real(8), allocatable :: v(:)
  contains
    procedure :: InitSixJs_final
    procedure :: FinSixJs_final

    generic :: init => InitSixJs_final
    generic :: fin => FinSixJs_final
  end type SixJs_final

  type :: SixJs_intermediate
    type(SixJs_final), allocatable :: v(:,:)
    integer :: j12min, j12max
    integer :: j23min, j23max
  contains
    procedure :: InitSixJs_intermediate
    procedure :: FinSixJs_intermediate

    generic :: init => InitSixJs_intermediate
    generic :: fin => FinSixJs_intermediate
  end type SixJs_intermediate

  type :: SixJsStore
    type(SixJs_intermediate), allocatable :: v(:,:,:)
    integer :: j1min, j1max
    integer :: j2min, j2max
    integer :: j3min, j3max
    logical :: half_j1
    logical :: half_j2
    logical :: half_j3
    logical :: half_j12
    logical :: half_j23
    logical :: half_j
  contains
    procedure :: InitSixJsStore
    procedure :: FinSixJsStore
    procedure :: GetStoredSixJ

    generic :: init => InitSixJsStore
    generic :: fin => FinSixJsStore
    generic :: get => GetStoredSixJ
  end type SixJsStore

  ! 9-j symbol storing
  !
  ! ninej%v(j1i,j2i,j3i,j4i)%v(j12i,j34i,j13i,j24i)%v(ji)
  !   {  j1  j2 j12 }
  ! = {  j3  j4 j34 }
  !   { j13 j24   J }
  !
  !   j1i = j1       (j1: integer)
  !   j1i = (j1+1)/2 (j1: half integer)
  type :: NineJs_final
    real(8), allocatable :: v(:)
  contains
    procedure :: InitNineJs_final
    procedure :: FinNineJs_final

    generic :: init => InitNineJs_final
    generic :: fin => FinNineJs_final
  end type NineJs_final

  type :: NineJs_intermediate
    type(NineJs_final), allocatable :: v(:,:)
    integer, allocatable :: jds2i_1234(:,:), jds2i_1324(:,:)
    integer, allocatable :: i2j12d(:), i2j34d(:), i2j13d(:), i2j24d(:)
    integer :: imax_1234, imax_1324
  contains
    procedure :: InitNineJs_intermediate
    procedure :: FinNineJs_intermediate

    generic :: init => InitNineJs_intermediate
    generic :: fin => FinNineJs_intermediate
  end type NineJs_intermediate

  type :: NineJsStore
    type(NineJs_intermediate), allocatable :: v(:,:,:,:)
    integer :: j1min, j1max
    integer :: j2min, j2max
    integer :: j3min, j3max
    integer :: j4min, j4max
    logical :: half_j1
    logical :: half_j2
    logical :: half_j3
    logical :: half_j4
    logical :: half_j
  contains
    procedure :: InitNineJsStore
    procedure :: FinNineJsStore
    procedure :: GetStoredNineJ

    generic :: init => InitNineJsStore
    generic :: fin => FinNineJsStore
    generic :: get => GetStoredNineJ
  end type NineJsStore

  ! Talmi-Moshinsky bracket storing
  !
  type :: TMbk_final
    real(8), allocatable :: v(:)
    real(8), allocatable :: vinv(:)
  end type TMbk_final

  type :: TMbk_intermediate
    type(TMbk_final), allocatable :: bk(:,:)
    integer, allocatable :: idx(:,:,:,:)
    integer, allocatable :: ll1(:)
    integer, allocatable :: ll2(:)
    integer, allocatable :: nn1(:)
    integer, allocatable :: nn2(:)
    integer :: n_idx
  contains
    procedure :: InitTMbk_intermediate
    procedure :: FinTMbk_intermediate

    generic :: init => InitTMbk_intermediate
    generic :: fin => FinTMbk_intermediate
  end type TMbk_intermediate

  type :: TMbracketStore
    type(TMbk_intermediate), allocatable :: N(:)
    real(8) :: mass_ratio
    integer :: Nmax
  contains
    procedure :: InitTMbracketStore
    procedure :: FinTMbracketStore
    procedure :: GetStoredTMbracket

    generic :: init => InitTMbracketStore
    generic :: fin => FinTMbracketStore
    generic :: get => GetStoredTMbracket
  end type TMbracketStore

contains
  !
  !
  !  CG storing
  !
  !
  subroutine FinCGsStore(this)
    class(CGsStore), intent(inout) :: this
    integer :: j1, j2
    do j1 = this%j1min, this%j1max
      do j2 = this%j2min, this%j2max
        call this%v(j1,j2)%fin()
      end do
    end do
    deallocate(this%v)
  end subroutine FinCGsStore

  subroutine InitCGsStore(this, j1dmin, j1dmax, half_j1, &
        & j2dmin, j2dmax, half_j2)
    class(CGsStore), intent(inout) :: this
    integer, intent(in) :: j1dmin, j1dmax, j2dmin, j2dmax
    logical, intent(in) :: half_j1, half_j2
    logical :: half_j12
    integer :: j1min, j1max, j2min, j2max
    integer :: j1, j2, j1d, j2d
    real(8) :: ti

    ti = omp_get_wtime()
    call timer%cmemory()

    j1min = j1dmin/2
    j1max = j1dmax/2
    j2min = j2dmin/2
    j2max = j2dmax/2
    if(half_j1) j1min = (j1dmin+1)/2
    if(half_j1) j1max = (j1dmax+1)/2
    if(half_j2) j2min = (j2dmin+1)/2
    if(half_j2) j2max = (j2dmax+1)/2

    this%j1min = j1min
    this%j1max = j1max
    this%j2min = j2min
    this%j2max = j2max
    this%half_j1 = half_j1
    this%half_j2 = half_j2

    half_j12 = half_j1 .neqv. half_j2
    this%half_j12 = half_j12

    allocate(this%v(j1min:j1max,j2min:j2max))

    do j1 = j1min, j1max
      j1d = 2*j1
      if(half_j1) j1d = 2*j1-1
      do j2 = j2min, j2max
        j2d = 2*j2
        if(half_j2) j2d = 2*j2-1
        call this%v(j1,j2)%init(j1d,j2d,half_j1,half_j2,half_j12)
      end do
    end do
    call timer%countup_memory("Storing CG coefs")
    call timer%add('Storing CG coefs', omp_get_wtime() - ti)
  end subroutine InitCGsStore

  function GetStoredCG(this, j1d, m1d, j2d, m2d, jd, md) result(r)
    ! inputs are double j
    class(CGsStore), intent(in) :: this
    integer, intent(in) :: j1d, m1d, j2d, m2d, jd, md
    integer :: j1, j2, m1, m2, j
    real(8) :: r

    r = 0.d0
    if(m1d + m2d /= md) then
      write(*,"(a)") "Warning, GetStoredCG: m1+m2 != m3"
      return
    end if
    j1 = j1d / 2
    j2 = j2d / 2
    m1 = m1d / 2
    m2 = m2d / 2
    j = jd / 2

    if(this%half_j1 ) then
      j1  = (j1d +1) / 2
      m1  = (m1d +1) / 2
    end if
    if(this%half_j2 ) then
      j2  = (j2d +1) / 2
      m2  = (m2d +1) / 2
    end if

    if(this%half_j12 ) then
      j  = (jd +1) / 2
    end if

    r = this%v(j1,j2)%v(m1,m2,j)
  end function GetStoredCG

  subroutine FinCGs_final(this)
    class(CGs_final), intent(inout) :: this
    if(.not. allocated(this%v)) return
    deallocate(this%v)
  end subroutine FinCGs_final

  subroutine InitCGs_final(this,j1d,j2d,half_j1,half_j2,half_j12)
    use MyLibrary, only: dcg
    class(CGs_final), intent(inout) :: this
    integer, intent(in) :: j1d, j2d
    logical, intent(in) :: half_j1, half_j2, half_j12
    integer :: jdmax, jdmin, jmax, jmin
    integer :: m1dmin, m1dmax, m1min, m1max
    integer :: m2dmin, m2dmax, m2min, m2max
    integer :: j, jd, m1, m1d, m2, m2d

    jdmin = abs(j1d-j2d)
    jdmax =     j1d+j2d
    m1dmin = -j1d; m1dmax = j1d
    m2dmin = -j2d; m2dmax = j2d

    jmin = jdmin / 2
    jmax = jdmax / 2
    if(half_j12) then
      jmin = (jdmin+1) / 2
      jmax = (jdmax+1) / 2
    end if

    m1min = m1dmin / 2
    m1max = m1dmax / 2
    if(half_j1) then
      m1min = (m1dmin+1)/2
      m1max = (m1dmax+1)/2
    end if

    m2min = m2dmin / 2
    m2max = m2dmax / 2
    if(half_j2) then
      m2min = (m2dmin+1)/2
      m2max = (m2dmax+1)/2
    end if

    allocate(this%v(m1min:m1max,m2min:m2max,jmin:jmax))
    do m1 = m1min, m1max
      m1d = 2*m1
      if(half_j1) m1d = 2*m1-1
      do m2 = m2min, m2max
        m2d = 2*m2
        if(half_j2) m2d = 2*m2-1
        do j = jmin, jmax
          jd = 2*j
          if(half_j12) jd = 2*j-1
          this%v(m1,m2,j) = dcg(j1d,m1d,j2d,m2d,jd,m1d+m2d)
        end do
      end do
    end do
  end subroutine InitCGs_final
  !
  !
  !  end CG storing
  !
  !
  !
  !
  !  6-j storing
  !
  !
  subroutine FinSixJsStore(this)
    class(SixJsStore), intent(inout) :: this
    integer :: j1, j2, j3
    do j1 = this%j1min, this%j1max
      do j2 = this%j2min, this%j2max
        do j3 = this%j3min, this%j3max
          call this%v(j1,j2,j3)%fin()
        end do
      end do
    end do
    deallocate(this%v)
  end subroutine FinSixJsStore

  subroutine InitSixJsStore(this, j1dmin, j1dmax, half_j1, &
        & j2dmin, j2dmax, half_j2, j3dmin, j3dmax, half_j3)
    class(SixJsStore), intent(inout) :: this
    integer, intent(in) :: j1dmin, j1dmax, j2dmin, j2dmax, j3dmin, j3dmax
    logical, intent(in) :: half_j1, half_j2, half_j3
    logical :: half_j12, half_j23, half_j123, half_j231
    integer :: j1min, j1max, j2min, j2max, j3min, j3max
    integer :: j1, j2, j3, j1d, j2d, j3d
    real(8) :: ti

    ti = omp_get_wtime()
    call timer%cmemory()

    j1min = j1dmin/2
    j1max = j1dmax/2
    j2min = j2dmin/2
    j2max = j2dmax/2
    j3min = j3dmin/2
    j3max = j3dmax/2
    if(half_j1) j1min = (j1dmin+1)/2
    if(half_j1) j1max = (j1dmax+1)/2
    if(half_j2) j2min = (j2dmin+1)/2
    if(half_j2) j2max = (j2dmax+1)/2
    if(half_j3) j3min = (j3dmin+1)/2
    if(half_j3) j3max = (j3dmax+1)/2

    this%j1min = j1min
    this%j1max = j1max
    this%j2min = j2min
    this%j2max = j2max
    this%j3min = j3min
    this%j3max = j3max
    this%half_j1 = half_j1
    this%half_j2 = half_j2
    this%half_j3 = half_j3

    half_j12 = half_j1 .neqv. half_j2
    half_j23 = half_j2 .neqv. half_j3
    half_j123 = half_j12 .neqv. half_j3
    half_j231 = half_j23 .neqv. half_j1
    if(half_j123 .neqv. half_j231) then
      write(*,*) "Error occures in 6-j symbol storing"
      return
    end if

    this%half_j12 = half_j12
    this%half_j23 = half_j23
    this%half_j = half_j123

    allocate(this%v(j1min:j1max,j2min:j2max,j3min:j3max))

    do j1 = j1min, j1max
      j1d = 2*j1
      if(half_j1) j1d = 2*j1-1
      do j2 = j2min, j2max
        j2d = 2*j2
        if(half_j2) j2d = 2*j2-1
        do j3 = j3min, j3max
          j3d = 2*j3
          if(half_j3) j3d = 2*j3-1
          call this%v(j1,j2,j3)%init(j1d,j2d,j3d, &
              & half_j1, half_j3, half_j12,half_j23)
        end do
      end do
    end do
    call timer%countup_memory("Storing 6-j couplings")
    call timer%add('Storing 6-j couplings', omp_get_wtime() - ti)
  end subroutine InitSixJsStore

  function GetStoredSixJ(this, j1d, j2d, j12d, j3d, jd, j23d) result(r)
    ! inputs are double j
    class(SixJsStore), intent(in) :: this
    integer, intent(in) :: j1d, j2d, j12d, j3d, jd, j23d
    integer :: j1, j2, j12, j3, j, j23
    real(8) :: r

    j1 = j1d / 2
    j2 = j2d / 2
    j3 = j3d / 2
    j12 = j12d / 2
    j23 = j23d / 2
    j = jd / 2
    if(this%half_j1 ) j1  = (j1d +1) / 2
    if(this%half_j2 ) j2  = (j2d +1) / 2
    if(this%half_j3 ) j3  = (j3d +1) / 2
    if(this%half_j12) j12 = (j12d+1) / 2
    if(this%half_j23) j23 = (j23d+1) / 2
    if(this%half_j  ) j   = (jd  +1) / 2

    r = this%v(j1,j2,j3)%v(j12,j23)%v(j)
  end function GetStoredSixJ

  subroutine FinSixJs_intermediate(this)
    class(SixJs_intermediate), intent(inout) :: this
    integer :: j12, j23
    do j12 = this%j12min, this%j12max
      do j23 = this%j23min, this%j23max
        call this%v(j12,j23)%fin()
      end do
    end do
    deallocate(this%v)
  end subroutine FinSixJs_intermediate

  subroutine InitSixJs_intermediate(this,j1d,j2d,j3d,half_j1,half_j3,&
        & half_j12,half_j23)
    class(Sixjs_intermediate), intent(inout) :: this
    integer, intent(in) :: j1d, j2d, j3d
    logical, intent(in) :: half_j12, half_j23, half_j1, half_j3
    integer :: j12, j23, j12d, j23d
    integer :: j12dmin, j12dmax, j12min, j12max
    integer :: j23dmin, j23dmax, j23min, j23max
    logical :: half_j, half_j123, half_j231

    j12dmin = abs(j1d - j2d)
    j12dmax =     j1d + j2d
    j23dmin = abs(j2d - j3d)
    j23dmax =     j2d + j3d

    j12min = j12dmin/2
    j12max = j12dmax/2
    j23min = j23dmin/2
    j23max = j23dmax/2
    if(half_j12) then
      j12min = (j12dmin+1)/2
      j12max = (j12dmax+1)/2
    end if

    if(half_j23) then
      j23min = (j23dmin+1)/2
      j23max = (j23dmax+1)/2
    end if

    half_j123 = half_j12 .neqv. half_j3
    half_j231 = half_j23 .neqv. half_j1
    if(half_j123 .neqv. half_j231) then
      write(*,*) "Error occures in 6-j symbol storing"
      return
    end if
    half_j = half_j123

    this%j12min = j12min
    this%j12max = j12max
    this%j23min = j23min
    this%j23max = j23max

    allocate(this%v(j12min:j12max,j23min:j23max))
    do j12 = this%j12min, this%j12max
      j12d = 2*j12
      if(half_j12) j12d = 2*j12-1
      do j23 = this%j23min, this%j23max
        j23d = 2*j23
        if(half_j23) j23d = 2*j23-1
        call this%v(j12,j23)%init(j1d,j2d,j3d,j12d,j23d,half_j)
      end do
    end do
  end subroutine InitSixJs_intermediate

  subroutine FinSixJs_final(this)
    class(Sixjs_final), intent(inout) :: this
    if(.not. allocated(this%v)) return
    deallocate(this%v)
  end subroutine FinSixJs_final

  subroutine InitSixJs_final(this,j1d,j2d,j3d,j12d,j23d,half_j)
    use MyLibrary, only: sjs
    class(Sixjs_final), intent(inout) :: this
    integer, intent(in) :: j1d, j2d, j3d, j12d, j23d
    logical, intent(in) :: half_j
    integer :: jdmax, jdmin, jd, jmax, jmin, j
    jdmin = max(abs(j12d-j3d), abs(j1d-j23d))
    jdmax = min(   (j12d+j3d),    (j1d+j23d))
    if(jdmin > jdmax) return
    jmin = jdmin / 2
    jmax = jdmax / 2
    if(half_j) then
      jmin = (jdmin+1) / 2
      jmax = (jdmax+1) / 2
    end if

    allocate(this%v(jmin:jmax))
    do j = jmin, jmax
      jd = 2*j
      if(half_j) jd = 2*j - 1
      this%v(j) = sjs(j1d,j2d,j12d,j3d,jd,j23d)
    end do
  end subroutine InitSixJs_final
  !
  !
  !  end 6-j storing
  !
  !

  !
  !
  !  9-j storing
  !
  !
  subroutine FinNineJsStore(this)
    class(NineJsStore), intent(inout) :: this
    integer :: j1, j2, j3, j4
    do j1 = this%j1min, this%j1max
      do j2 = this%j2min, this%j2max
        do j3 = this%j3min, this%j3max
          do j4 = this%j4min, this%j4max
            call this%v(j1,j2,j3,j4)%fin()
          end do
        end do
      end do
    end do
    deallocate(this%v)
  end subroutine FinNineJsStore

  subroutine InitNineJsStore(this, j1dmin, j1dmax, half_j1, j2dmin, j2dmax, half_j2, &
        & j3dmin, j3dmax, half_j3, j4dmin, j4dmax, half_j4, jd)
    class(NineJsStore), intent(inout) :: this
    integer, intent(in) :: j1dmin, j1dmax, j2dmin, j2dmax, j3dmin, j3dmax, j4dmin, j4dmax
    logical, intent(in) :: half_j1, half_j2, half_j3, half_j4
    integer, intent(in), optional :: Jd
    logical :: half_j12, half_j34, half_j13, half_j24, half_j1234, half_j1324
    integer :: j1min, j1max, j2min, j2max, j3min, j3max, j4min, j4max
    integer :: j1d, j2d, j3d, j4d
    integer :: j1, j2, j3, j4
    real(8) :: ti

    ti = omp_get_wtime()
    call timer%cmemory()

    j1min = j1dmin/2
    j1max = j1dmax/2
    j2min = j2dmin/2
    j2max = j2dmax/2
    j3min = j3dmin/2
    j3max = j3dmax/2
    j4min = j4dmin/2
    j4max = j4dmax/2
    if(half_j1) j1min = (j1dmin+1)/2
    if(half_j1) j1max = (j1dmax+1)/2
    if(half_j2) j2min = (j2dmin+1)/2
    if(half_j2) j2max = (j2dmax+1)/2
    if(half_j3) j3min = (j3dmin+1)/2
    if(half_j3) j3max = (j3dmax+1)/2
    if(half_j4) j4min = (j4dmin+1)/2
    if(half_j4) j4max = (j4dmax+1)/2

    this%j1min = j1min
    this%j1max = j1max
    this%j2min = j2min
    this%j2max = j2max
    this%j3min = j3min
    this%j3max = j3max
    this%j4min = j4min
    this%j4max = j4max
    this%half_j1 = half_j1
    this%half_j2 = half_j2
    this%half_j3 = half_j3
    this%half_j4 = half_j4

    half_j12 = half_j1 .neqv. half_j2
    half_j34 = half_j3 .neqv. half_j4
    half_j13 = half_j1 .neqv. half_j3
    half_j24 = half_j2 .neqv. half_j4

    half_j1234 = half_j12 .neqv. half_j34
    half_j1324 = half_j13 .neqv. half_j24
    if(half_j1234 .neqv. half_j1324) then
      write(*,*) "Error occures in 9-j symbol storing"
      return
    end if
    this%half_j = half_j1234

    allocate(this%v(j1min:j1max,j2min:j2max,j3min:j3max,j4min:j4max))

    do j1 = j1min, j1max
      j1d = 2*j1
      if(half_j1) j1d = 2*j1-1
      do j2 = j2min, j2max
        j2d = 2*j2
        if(half_j2) j2d = 2*j2-1
        do j3 = j3min, j3max
          j3d = 2*j3
          if(half_j3) j3d = 2*j3-1
          do j4 = j4min, j4max
            j4d = 2*j4
            if(half_j4) j4d = 2*j4-1
            call this%v(j1,j2,j3,j4)%init(j1d,j2d,j3d,j4d, &
                & half_j12,half_j34,half_j13,half_j24,jd)
          end do
        end do
      end do
    end do

    call timer%countup_memory("Storing 9-j couplings")
    call timer%add('Storing 9-j couplings', omp_get_wtime() - ti)
  end subroutine InitNineJsStore

  function GetStoredNineJ(this, j1d, j2d, j12d, j3d, j4d, j34d, j13d, j24d, jd) result(r)
    ! inputs are double j
    class(NineJsStore), intent(in) :: this
    integer, intent(in) :: j1d, j2d, j3d, j4d, j12d, j34d, j13d, j24d, jd
    integer :: j1, j2, j3, j4, j
    integer :: i1234, i1324
    real(8) :: r

    j1 = j1d / 2
    j2 = j2d / 2
    j3 = j3d / 2
    j4 = j4d / 2
    j = jd / 2

    if(this%half_j1 ) j1  = (j1d +1) / 2
    if(this%half_j2 ) j2  = (j2d +1) / 2
    if(this%half_j3 ) j3  = (j3d +1) / 2
    if(this%half_j4 ) j4  = (j4d +1) / 2
    if(this%half_j  ) j   = (jd  +1) / 2

    i1234 = this%v(j1,j2,j3,j4)%jds2i_1234(j12d,j34d)
    i1324 = this%v(j1,j2,j3,j4)%jds2i_1324(j13d,j24d)
    if(i1234 * i1324 == 0) then
      write(*,"(a)") "Warning: GetStoredNineJ"
      return
    end if

    r = this%v(j1,j2,j3,j4)%v(i1234,i1324)%v(j)
  end function GetStoredNineJ

  subroutine FinNineJs_intermediate(this)
    class(NineJs_intermediate), intent(inout) :: this
    integer :: i1234, i1324
    do i1234 = 1, this%imax_1234
      do i1324 = 1, this%imax_1324
        call this%v(i1234,i1324)%fin()
      end do
    end do
    deallocate(this%v)
  end subroutine FinNineJs_intermediate

  subroutine InitNineJs_intermediate(this,j1d,j2d,j3d,j4d,&
        & half_j12,half_j34,half_j13,half_j24,jd)
    use MyLibrary, only: triag
    class(NineJs_intermediate), intent(inout) :: this
    integer, intent(in) :: j1d,j2d,j3d,j4d
    logical, intent(in) :: half_j12, half_j34, half_j13, half_j24
    integer, intent(in), optional :: jd
    logical :: half_j
    integer :: j12,j34,j13,j24
    integer :: j12d,j34d,j13d,j24d
    integer :: j12dmin,j12dmax
    integer :: j34dmin,j34dmax
    integer :: j13dmin,j13dmax
    integer :: j24dmin,j24dmax
    integer :: cnt, i1234, i1324

    j12dmin = abs(j1d-j2d)
    j34dmin = abs(j3d-j4d)
    j13dmin = abs(j1d-j3d)
    j24dmin = abs(j2d-j4d)
    j12dmax =     j1d+j2d
    j34dmax =     j3d+j4d
    j13dmax =     j1d+j3d
    j24dmax =     j2d+j4d

    half_j = half_j12 .neqv. half_j34
    cnt = 0
    do j12d = j12dmin, j12dmax, 2
      do j34d = j34dmin, j34dmax, 2
        if(present(jd)) then
          if(triag(j12d, j34d, jd)) cycle
        end if
        cnt = cnt + 1
      end do
    end do
    this%imax_1234 = cnt
    allocate(this%jds2i_1234(j12dmin:j12dmax, j34dmin:j34dmax))
    allocate(this%i2j12d(this%imax_1234))
    allocate(this%i2j34d(this%imax_1234))
    this%jds2i_1234(:,:) = 0
    cnt = 0
    do j12d = j12dmin, j12dmax, 2
      do j34d = j34dmin, j34dmax, 2
        if(present(jd)) then
          if(triag(j12d, j34d, jd)) cycle
        end if
        cnt = cnt + 1
        this%jds2i_1234(j12d,j34d) = cnt
        this%i2j12d(cnt) = j12d
        this%i2j34d(cnt) = j34d
      end do
    end do

    cnt = 0
    do j13d = j13dmin, j13dmax, 2
      do j24d = j24dmin, j24dmax, 2
        if(present(jd)) then
          if(triag(j13d, j24d, jd)) cycle
        end if
        cnt = cnt + 1
      end do
    end do
    this%imax_1324 = cnt
    allocate(this%jds2i_1324(j13dmin:j13dmax, j24dmin:j24dmax))
    allocate(this%i2j13d(this%imax_1324))
    allocate(this%i2j24d(this%imax_1324))
    this%jds2i_1324(:,:) = 0
    cnt = 0
    do j13d = j13dmin, j13dmax, 2
      do j24d = j24dmin, j24dmax, 2
        if(present(jd)) then
          if(triag(j13d, j24d, jd)) cycle
        end if
        cnt = cnt + 1
        this%jds2i_1324(j13d,j24d) = cnt
        this%i2j13d(cnt) = j13d
        this%i2j24d(cnt) = j24d
      end do
    end do

    allocate(this%v(this%imax_1234,this%imax_1324))
    do i1234 = 1, this%imax_1234
      j12d = this%i2j12d(i1234)
      j34d = this%i2j34d(i1234)
      do i1324 = 1, this%imax_1324
        j13d = this%i2j13d(i1324)
        j24d = this%i2j24d(i1324)
        call this%v(i1234,i1324)%init(j1d,j2d,j3d,j4d,&
            & j12d,j34d,j13d,j24d,half_j,jd)
      end do
    end do
  end subroutine InitNineJs_intermediate

  subroutine FinNineJs_final(this)
    class(NineJs_final), intent(inout) :: this
    if(.not. allocated(this%v)) return
    deallocate(this%v)
  end subroutine FinNineJs_final

  subroutine InitNineJs_final(this,j1d,j2d,j3d,j4d,j12d,j34d,j13d,j24d,half_j,jd)
    use MyLibrary, only: snj
    class(NineJs_final), intent(inout) :: this
    integer, intent(in) :: j1d,j2d,j3d,j4d,j12d,j34d,j13d,j24d
    logical, intent(in) :: half_j
    integer, intent(in), optional :: jd
    integer :: j, jj
    integer :: jmin=-1, jmax=-1
    integer :: jdmin, jdmax


    jdmin = max(abs(j12d-j34d),abs(j13d-j24d))
    jdmax = min(   (j12d+j34d),   (j13d+j24d))
    if(present(jd)) then
      jdmin = jd
      jdmax = jd
    end if

    if(jdmin > jdmax) then
      return
    end if

    jmin = jdmin/2
    jmax = jdmax/2
    if(half_j) then
      jmin = (jdmin + 1)/2
      jmax = (jdmax + 1)/2
    end if

    allocate(this%v(jmin:jmax))
    do j = jmin, jmax
      jj = 2*j
      if(half_j) jj = 2*j-1
      this%v(j) = snj(j1d,j2d,j12d,j3d,j4d,j34d,j13d,j24d,jj)
    end do

  end subroutine InitNineJs_final

  !
  !
  !  end 9-j storing
  !
  !

  !
  !
  !  Talmi-Moshinsky bracket storing
  !
  !

  subroutine FinTMbracketStore(this)
    class(TMbracketStore), intent(inout) :: this
    integer :: N
    do N = 0, this%Nmax
      call this%N(N)%fin()
    end do
    deallocate(this%N)
  end subroutine FinTMbracketStore

  subroutine InitTMbracketStore(this, Nmax, d)
    class(TMbracketStore), intent(inout) :: this
    integer, intent(in) :: Nmax
    real(8), intent(in) :: d
    integer :: Ntot
    real(8) :: ti

    ti = omp_get_wtime()
    call timer%cmemory()

    this%mass_ratio = d
    this%Nmax = Nmax
    allocate(this%N(0:Nmax))
    do ntot = 0, Nmax
      call this%N(ntot)%init(ntot, d)
    end do

    call timer%countup_memory("Storing Talmi-Moshinsky brackets")
    call timer%add('Storing Talmi-Moshinsky brackets', omp_get_wtime() - ti)
  end subroutine InitTMbracketStore

  function GetStoredTMbracket(this, n1, l1, n2, l2, n3, l3, n4, l4, lam) result(r)
    class(TMbracketStore), intent(inout) :: this
    integer, intent(in) :: n1, l1, n2, l2, n3, l3, n4, l4, lam
    integer :: Nmax, idx_bra, idx_ket
    real(8) :: r

    r = 0.d0
    Nmax = 2*n1 + l1 + 2*n2 + l2
    if(Nmax /= 2*n3 + l3 + 2*n4 + l4) return

    idx_bra = this%N(Nmax)%idx(n1,l1,n2,l2)
    idx_ket = this%N(Nmax)%idx(n3,l3,n4,l4)
    if(idx_bra * idx_ket == 0) then
      write(*,*) "error : GetStoreadTMbracket !"
      return
    end if
    r = this%N(Nmax)%bk(idx_bra,idx_ket)%v(lam)
  end function GetStoredTMbracket

  subroutine FinTMbk_intermediate(this)
    class(TMbk_intermediate), intent(inout) :: this
    integer :: bra, ket
    do bra = 1, this%n_idx
      do ket = 1, this%n_idx
        deallocate(this%bk(bra,ket)%v)
        deallocate(this%bk(bra,ket)%vinv)
      end do
    end do
    deallocate(this%bk)
    deallocate(this%ll1)
    deallocate(this%ll2)
    deallocate(this%nn1)
    deallocate(this%nn2)
    deallocate(this%idx)
  end subroutine FinTMbk_intermediate

  subroutine InitTMbk_intermediate(this, Nmax, d)
    use MyLibrary, only: gmosh
    class(TMbk_intermediate), intent(inout) :: this
    integer, intent(in) :: Nmax
    real(8), intent(in) :: d
    integer :: idx, n1, l1, n2, l2, n
    integer :: bra, ket, n3, l3, n4, l4
    integer :: lammin, lammax, lam

    allocate(this%idx(0:Nmax/2, 0:Nmax, 0:Nmax/2, 0:Nmax))
    idx = 0
    do n1 = 0, Nmax / 2
      do l1 = 0, Nmax - 2 * n1
        do n2 = 0, (Nmax - 2 * n1 - l1) / 2
          l2 = Nmax - 2 * n1 - 2 * n2 - l1
          idx = idx + 1
        end do
      end do
    end do
    n = idx
    this%n_idx = n
    allocate(this%bk(n,n))
    allocate(this%ll1(n))
    allocate(this%ll2(n))
    allocate(this%nn1(n))
    allocate(this%nn2(n))

    idx = 0
    do n1 = 0, Nmax / 2
      do l1 = 0, Nmax - 2 * n1
        do n2 = 0, (Nmax - 2 * n1 - l1) / 2
          l2 = Nmax - 2 * n1 - 2 * n2 - l1
          idx = idx + 1
          this%ll1(idx) = l1
          this%ll2(idx) = l2
          this%nn1(idx) = n1
          this%nn2(idx) = n2
          this%idx(n1, l1, n2, l2) = idx
        end do
      end do
    end do

    n = this%n_idx
    !$omp parallel
    !$omp do private(bra, l1, l2, n1, n2, ket, l3, l4, n3, n4, &
    !$omp &          lammin, lammax, lam) schedule(dynamic)
    do bra = 1, n
      l1 = this%ll1(bra)
      l2 = this%ll2(bra)
      n1 = this%nn1(bra)
      n2 = this%nn2(bra)
      do ket = 1, n
        l3 = this%ll1(ket)
        l4 = this%ll2(ket)
        n3 = this%nn1(ket)
        n4 = this%nn2(ket)

        lammin = max(abs(l1 - l2), abs(l3 - l4))
        lammax = min(l1 + l2, l3 + l4)
        allocate(this%bk(bra, ket)%v(lammin:lammax))
        allocate(this%bk(bra, ket)%vinv(lammin:lammax))
        this%bk(bra, ket)%v(:) = 0.d0
        this%bk(bra, ket)%vinv(:) = 0.d0
        do lam = lammin, lammax
          this%bk(bra, ket)%v(lam)    = gmosh(n1, l1, n2, l2, n3, l3, n4, l4, lam, d)
          this%bk(bra, ket)%vinv(lam) = gmosh(n1, l1, n2, l2, n3, l3, n4, l4, lam, 1.d0/d)
        end do
      end do
    end do
    !$omp end do
    !$omp end parallel
  end subroutine InitTMbk_intermediate

end module StoreCouplings
