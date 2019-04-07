module Operators
  use omp_lib
  use ModelSpace
  use OneBodyOperator
  use TwoBodyOperator
  use ThreeBodyOperator
  use ThreeBodyInteraction
  implicit none

  public :: Ops

  private :: FinOps
  private :: InitOps
  private :: InitOpsFromString
  private :: SetOperatorFromFile
  private :: SetOperator
  private :: PrintOperator
  private :: CopyOps
  private :: SumOps
  private :: SubtractOps
  private :: ScaleOps
  private :: NormalOrdering
  private :: NO2BApprox
  !private :: UnNormalOrdering
  private :: DiscardThreeBodyForce
  private :: DiscardThreeBodyPart
  !private :: NO2BApproximation

  type :: Ops
    character(32) :: oprtr
    real(8) :: zero
    type(MSpace), pointer :: ms
    type(OneBodyPart) :: one
    type(TwoBodyPart) :: two
    type(ThreeBodyPart) :: thr
    type(ThreeBodyForce) :: thr21
    logical :: Scalar=.false.
    logical :: is_normal_ordered = .false.
    integer :: jr, pr, zr
    integer :: rank
  contains
    procedure :: FinOps
    procedure :: InitOps
    procedure :: InitOpsFromString
    procedure :: SetOperatorFromFile
    procedure :: SetOperator
    procedure :: PrintOperator
    procedure :: CopyOps
    procedure :: SumOps
    procedure :: SubtractOps
    procedure :: ScaleOps

    procedure :: NormalOrdering
    procedure :: NO2BApprox
    procedure :: UnNormalOrdering2B
    procedure :: DiscardThreeBodyForce
    procedure :: DiscardThreeBodyPart

    generic :: init => InitOps, InitOpsFromString
    generic :: fin => FinOps
    generic :: set => SetOperatorFromFile, SetOperator
    generic :: prt => PrintOperator
    generic :: assignment(=) => CopyOps
    generic :: operator(+) => SumOps
    generic :: operator(-) => SubtractOps
    generic :: operator(*) => ScaleOps
  end type Ops

contains
  subroutine FinOps(this)
    class(Ops), intent(inout) :: this
    this%zero = 0.d0
    call this%one%fin()
    call this%two%fin()
    if(this%rank == 3 .and. this%ms%is_three_body_jt) call this%thr21%fin()
    if(this%rank == 3 .and. this%ms%is_three_body) call this%thr%fin()
    this%ms => null()
    this%is_normal_ordered = .false.
  end subroutine FinOps

  subroutine InitOpsFromString(this, optr, ms, rank)
    use DefineOperators, only: GetOperatorRank
    class(Ops), intent(inout) :: this
    type(MSpace), intent(in) :: ms
    character(*), intent(in) :: optr
    integer, intent(in) :: rank
    integer :: jr, pr, zr

    call GetOperatorRank(optr,jr,pr,zr)
    call this%init(jr,pr,zr,optr,ms,rank)
  end subroutine InitOpsFromString

  subroutine InitOps(this, jr, pr, zr, optr, ms, rank)
    use Profiler, only: timer
    class(Ops), intent(inout) :: this
    integer, intent(in) :: jr, pr, zr, rank
    character(*), intent(in) :: optr
    type(MSpace), intent(in), target :: ms
    real(8) :: ti

    ti = omp_get_wtime()
    call timer%cmemory()

    this%is_normal_ordered = .false.
    this%ms => ms
    this%jr = jr
    this%pr = pr
    this%zr = zr
    this%oprtr = optr
    this%rank = rank
    if(this%jr == 0 .and. this%pr == 1 .and. this%zr == 0) this%Scalar = .true.
    if(this%jr /= 0 .or.  this%pr /= 1 .or.  this%zr /= 0) this%Scalar = .false.

    this%zero = 0.d0
    call this%one%init(ms%one, this%Scalar, optr, jr, pr, zr)
    call this%two%init(ms%two, this%Scalar, optr, jr, pr, zr)

    if(rank == 3 .and. this%ms%is_three_body_jt) then
      call this%thr21%init(ms%thr21) ! Three-body interaction
    end if

    if(rank == 3 .and. this%ms%is_three_body) then
      call this%thr%init(ms%thr, this%Scalar, optr, jr, pr, zr)
    end if

    call timer%countup_memory(trim(optr))
    call timer%Add('Construct '//trim(optr), omp_get_wtime()-ti)
  end subroutine InitOps

  subroutine CopyOps(a, b)
    class(Ops), intent(inout) :: a
    type(Ops), intent(in) :: b

    a%oprtr = b%oprtr
    a%Scalar = b%Scalar
    a%is_normal_ordered = b%is_normal_ordered
    a%jr = b%jr
    a%pr = b%pr
    a%zr = b%zr
    a%rank = b%rank

    a%ms => b%ms
    a%zero = b%zero
    a%one = b%one
    a%two = b%two

    if(a%rank == 3 .and. a%ms%is_three_body_jt) then
      a%thr21 = b%thr21
    end if

    if(a%rank == 3 .and. a%ms%is_three_body) then
      a%thr = b%thr
    end if
  end subroutine CopyOps

  function SumOps(a, b) result(c)
    class(Ops), intent(in) :: a, b
    type(Ops) :: c

    c = a
    c%zero = a%zero + b%zero
    c%one = a%one + b%one
    c%two = a%two + b%two
    if(a%rank == 3 .and. a%ms%is_three_body_jt) then
      c%thr21 = a%thr21 + b%thr21
    end if

    if(a%rank == 3 .and. a%ms%is_three_body) then
      c%thr = a%thr + b%thr
    end if
  end function SumOps

  function SubtractOps(a, b) result(c)
    class(Ops), intent(in) :: a, b
    type(Ops) :: c

    c = a
    c%zero = a%zero - b%zero
    c%one = a%one - b%one
    c%two = a%two - b%two

    if(a%rank == 3 .and. a%ms%is_three_body_jt) then
      c%thr21 = a%thr21 - b%thr21
    end if

    if(a%rank == 3 .and. a%ms%is_three_body) then
      c%thr = a%thr - b%thr
    end if
  end function SubtractOps

  function ScaleOps(a, b) result(c)
    class(Ops), intent(in) :: a
    real(8), intent(in) :: b
    type(Ops) :: c

    c = a
    c%zero = a%zero * b
    c%one = a%one * b
    c%two = a%two * b
    if(a%rank == 3 .and. a%ms%is_three_body_jt) then
      c%thr21 = a%thr21 * b
    end if

    if(a%rank == 3 .and. a%ms%is_three_body) then
      c%thr = a%thr * b
    end if
  end function ScaleOps

  subroutine SetOperatorFromFile(this, file_nn, file_3n, &
        & bound_2b_file, bound_3b_file)
    class(Ops), intent(inout), target :: this
    character(*), intent(in) :: file_nn, file_3n
    integer, intent(in), optional :: bound_2b_file(3), bound_3b_file(4)
    type(OneBodyPart) :: Tcm_one, Hcm_one
    type(TwoBodyPart) :: Tcm_two, Hcm_two
    type(Read2BodyFiles) :: rd2
    type(Read3BodyFiles) :: rd3
    type(MSpace), pointer :: ms

    if(file_nn == 'none' .and. file_3n == 'none') return
    ms => this%ms
    write(*,'(3a)') "Set ", trim(this%oprtr), " operator from files"
    call rd2%set(file_nn)
    ! -- boundary for two-body file
    call rd2%set(this%ms%emax, this%ms%e2max, this%ms%e3max)
    if(present(bound_2b_file)) then
      call rd2%set(bound_2b_file(1), bound_2b_file(2), bound_2b_file(3))
    end if

    ! -- boundary for three-body file
    call rd3%set(file_3n)
    call rd3%set(ms%emax, ms%e2max, ms%e3max, ms%lmax)
    if(present(bound_3b_file)) then
      call rd3%set(bound_3b_file(1), bound_3b_file(2), bound_3b_file(3), bound_3b_file(4))
    end if

    call this%one%set(ms%hw, ms%A, ms%Z, ms%N)
    write(*,'(2a)') "2B file: ", trim(file_nn)
    call rd2%ReadTwoBodyFile(this%two)
    if(this%rank == 3 .and. this%ms%is_three_body_jt) then
      write(*,'(2a)') "3B file: ", trim(file_3n)
      call rd3%ReadIsospinThreeBodyFile(this%thr21)
    end if

    if(this%rank == 3 .and. this%ms%is_three_body_jt .and. this%ms%is_three_body) then
      call this%thr%set(this%thr21)
      call this%DiscardThreeBodyForce()
    end if

    select case(this%oprtr)
    case('hamil', 'Hamil')
      call Tcm_one%init(ms%one, .true., 'Tcm', 0, 1, 0)
      call Tcm_two%init(ms%two, .true., 'Tcm', 0, 1, 0)

      call Tcm_one%set(ms%hw, ms%A, ms%Z, ms%N)
      call Tcm_two%set(ms%hw, ms%A, ms%Z, ms%N)

      this%one = this%one - Tcm_one
      this%two = this%two - Tcm_two

      if(ms%beta > 1.d-4) then
        call Hcm_one%init(ms%one, .true., 'Hcm', 0, 1, 0)
        call Hcm_two%init(ms%two, .true., 'Hcm', 0, 1, 0)

        call Hcm_one%set(ms%hw, ms%A, ms%Z, ms%N)
        call Hcm_two%set(ms%hw, ms%A, ms%Z, ms%N)

        this%zero = this%zero - (1.5d0 * ms%hw * ms%beta)
        this%one = this%one + (Hcm_one * ms%beta)
        this%two = this%two + (Hcm_two * ms%beta)
      end if
    case default
      return
    end select
  end subroutine SetOperatorFromFile

  subroutine SetOperator(this)
    use Profiler, only: timer
    class(Ops), intent(inout), target :: this
    real(8) :: ti

    ti = omp_get_wtime()
    select case(this%oprtr)
    case('Hcm', 'HCM')
      this%zero = - 1.5d0 * this%ms%hw
    case default
    end select

    call this%one%set(this%ms%hw, this%ms%A, this%ms%Z, this%ms%N)
    call this%two%set(this%ms%hw, this%ms%A, this%ms%Z, this%ms%N)
    call timer%Add("Set operator "//trim(this%oprtr), omp_get_wtime()-ti)
  end subroutine SetOperator

  subroutine PrintOperator(this, iunit)
    class(Ops), intent(in) :: this
    integer, intent(in), optional :: iunit
    integer :: ut
    ut = 6
    if(present(iunit)) ut = iunit

    write(ut,'(a)') "## Zero-body part"
    write(ut,'(f12.6)') this%zero

    write(ut,'(a)') "## One-body part"
    call this%one%prt(iunit)

    write(ut,'(a)') "## Two-body part"
    call this%two%prt(iunit)

    if(this%rank==3 .and. this%ms%is_three_body_jt) then
      write(ut,'(a)') "## Three-body part (2+1, isospin assumed) "
      call this%thr21%prt(iunit)
    end if

    if(this%rank==3 .and. this%ms%is_three_body) then
      write(ut,'(a)') "## Three-body part"
      call this%thr%prt(iunit)
    end if
  end subroutine PrintOperator

  subroutine NormalOrdering(this) ! Normal ordering w.r.t. refernece state
    class(Ops), intent(inout), target :: this
    type(TwoBodyPart) :: op2from3
    type(OneBodyPart) :: op1from3, op1from2
    real(8) :: op0from1=0.d0, op0from2=0.d0, op0from3=0.d0

    if(this%is_normal_ordered) then
      write(*,*) "Operator ", trim(this%oprtr), " is already normal ordered!"
      return
    end if


    op1from2 = this%two%NormalOrderingFrom2To1(this%ms%one)

    op0from2 = op1from2%NormalOrderingFrom1To0()
    op0from1 = this%one%NormalOrderingFrom1To0()

    this%one = this%one + op1from2
    this%zero = this%zero + op0from1 + op0from2 * 0.5d0

    if(this%rank == 3) then
      if(this%ms%is_three_body_jt) then
        op2from3 = this%thr21%NormalOrderingFrom3To2(this%ms%two)
        op1from3 = op2from3%NormalOrderingFrom2To1(this%ms%one)
        op0from3 = op1from3%NormalOrderingFrom1To0()
        this%two = this%two + op2from3
        this%one = this%one + op1from3 * 0.5d0
        this%zero = this%zero + op0from3 / 6.d0
      end if

      if(this%ms%is_three_body) then
        op2from3 = this%thr%NormalOrderingFrom3To2(this%ms%two)
        op1from3 = op2from3%NormalOrderingFrom2To1(this%ms%one)
        op0from3 = op1from3%NormalOrderingFrom1To0()
        this%two = this%two + op2from3
        this%one = this%one + op1from3 * 0.5d0
        this%zero = this%zero + op0from3 / 6.d0
      end if
    end if
    this%is_normal_ordered = .true.
  end subroutine NormalOrdering

  subroutine NO2BApprox(this) ! Normal ordering w.r.t. refernece state
    class(Ops), intent(inout), target :: this
    call this%NormalOrdering()
    if(this%rank==3) then
      call this%DiscardThreeBodyForce()
      call this%DiscardThreeBodyPart()
    end if
    this%rank = 2
    this%is_normal_ordered = .true.
  end subroutine NO2BApprox

  subroutine UnNormalOrdering2B(this) ! From NO2B approx. Op
    class(Ops), intent(inout), target :: this
    type(OneBodyPart) :: op1from2
    real(8) :: op0from1, op0from2

    if(this%rank /= 2) then
      write(*,'(a)') "Warning: UnNormalOrdering2B is only for the NO2B approximated Operators."
    end if

    if(.not. this%is_normal_ordered) then
      write(*,*) "Operator ", trim(this%oprtr), " is already unnormal ordered!"
      return
    end if

    op1from2 = this%two%NormalOrderingFrom2To1(this%ms%one)

    op0from2 = op1from2%NormalOrderingFrom1To0()
    op0from1 = this%one%NormalOrderingFrom1To0()

    ! this%two is unchanged
    this%one = this%one - op1from2
    this%zero = this%zero - op0from1 + op0from2 * 0.5d0

    call op1from2%fin()
    this%is_normal_ordered = .false.
  end subroutine UnNormalOrdering2B

  subroutine DiscardThreeBodyForce(this)
    class(Ops), intent(inout) :: this
    if(.not. this%ms%is_three_body_jt) return
    call this%thr21%fin()
    call this%ms%ReleaseThreeBody21()
  end subroutine DiscardThreeBodyFOrce

  subroutine DiscardThreeBodyPart(this)
    class(Ops), intent(inout) :: this
    if(.not. this%ms%is_three_body) return
    call this%thr%fin()
    call this%ms%ReleaseThreeBody()
  end subroutine DiscardThreeBodyPart

end module Operators

!program test
!  use Profiler, only: timer
!  use ModelSpace, only: MSpace
!  use Operators
!
!  implicit none
!
!  type(MSpace) :: ms
!  type(Ops) :: h
!  character(:), allocatable :: file_nn, file_3n
!
!  call timer%init()
!  call ms%init('O16', 35.d0, 4, 8, e3max=6, is_three_body_jt=.true.)
!  call h%init('hamil',ms, 3)
!
!  file_nn = '/home/takayuki/MtxElmnt/2BME/TwBME-HO_NN-only_N3LO_EM500_srg2.00_hw35_emax8_e2max16.me2j.gz'
!  file_3n = '/home/takayuki/MtxElmnt/3BME/ThBME_srg2.00_N3LO_EM500_ChEFT_N2LO_cD-0.20cE0.098_Local2_IS_hw35_ms6_6_6.me3j'
!
!  call h%set(file_nn,file_3n,[8,16,8],[6,6,6,6])
!  call h%NormalOrdering()
!  write(*,*) h%zero
!  call h%fin()
!  call ms%fin()
!
!  call timer%fin()
!
!end program test
