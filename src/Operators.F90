module Operators
  use omp_lib
  use ModelSpace
  use NOperators
  implicit none

  type :: Op
    character(32) :: optr
    real(8) :: zero
    type(NBodyPart) :: one, two
!#ifdef single_precision
    type(NBodyPartSp) :: thr
!#else
!    type(NBodyPart) :: thr
!#endif
    logical :: is_normal_ordered = .false.
    logical :: Scalar
    logical :: is_three_body = .false.
    integer :: jr, pr, zr
  contains
    procedure :: fin => FinOp
    procedure :: InitOp
    procedure :: InitOpFromString
    generic :: init => InitOp, InitOpFromString

    procedure :: SetOperatorFromFile
    procedure :: SetOperator
    generic :: set => SetOperatorFromFile, SetOperator

    procedure :: prt => PrintOperator

    procedure :: CopyOp
    procedure :: SumOp
    procedure :: SubtractOp
    procedure :: ScaleOp
    generic :: assignment(=) => CopyOp
    generic :: operator(+) => SumOp
    generic :: operator(-) => SubtractOp
    generic :: operator(*) => ScaleOp

    procedure :: NormalOrdering
    procedure :: UnNormalOrdering
    procedure :: DiscardThreeBodyPart
    procedure :: NO2BApproximation
  end type Op

contains
  subroutine FinOp(this)
    class(Op), intent(inout) :: this
    call this%one%fin()
    call this%two%fin()
    if(.not. this%is_three_body) return
    call this%thr%fin()
  end subroutine FinOp

  subroutine InitOpFromString(this, optr, ms, is_three)
    use DefineOperators, only: GetOperatorRank
    class(Op), intent(inout) :: this
    type(MSpace), intent(in) :: ms
    character(*), intent(in) :: optr
    logical, intent(in) :: is_three
    integer :: jr, pr, zr

    call GetOperatorRank(optr,jr,pr,zr)
    call this%init(jr,pr,zr,optr,ms,is_three)
  end subroutine InitOpFromString

  subroutine InitOp(this, jr, pr, zr, optr, ms, is_three)
    use Profiler, only: timer
    class(Op), intent(inout) :: this
    integer, intent(in) :: jr, pr, zr
    character(*), intent(in) :: optr
    type(MSpace), intent(in) :: ms
    logical, intent(in) :: is_three
    real(8) :: ti

    ti = omp_get_wtime()
    call timer%cmemory()

    this%jr = jr
    this%pr = pr
    this%zr = zr
    this%optr = optr
    if(this%jr == 0 .and. this%pr == 1 .and. this%zr == 0) this%Scalar = .true.
    if(this%jr /= 0 .or.  this%pr /= 1 .or.  this%zr /= 0) this%Scalar = .false.

    call this%one%init(ms%one, this%Scalar, optr, jr, pr, zr)
    call this%two%init(ms%two, this%Scalar, optr, jr, pr, zr)
    if(.not. is_three) then
      call timer%countup_memory(trim(optr))
      call timer%Add('Construct '//trim(optr), omp_get_wtime()-ti)
      return
    end if
    this%is_three_body = .true.
    call this%thr%init(ms%thr, this%Scalar, optr, jr, pr, zr)

    call timer%countup_memory(trim(optr))
    call timer%Add('Construct '//trim(optr), omp_get_wtime()-ti)
  end subroutine InitOp

  subroutine CopyOp(a, b)
    class(Op), intent(inout) :: a
    type(Op), intent(in) :: b

    a%optr = b%optr
    a%Scalar = b%Scalar
    a%is_three_body = b%is_three_body
    a%jr = b%jr
    a%pr = b%pr
    a%zr = b%zr
    a%zero = b%zero
    a%one = b%one
    a%two = b%two
    if(.not. a%is_three_body) return
    a%thr = b%thr
  end subroutine CopyOp

  function SumOp(a, b) result(c)
    class(Op), intent(in) :: a, b
    type(Op) :: c

    c = a
    c%zero = a%zero + b%zero
    c%one = a%one + b%one
    c%two = a%two + b%two
    if(.not. a%is_three_body) return
    c%thr = a%thr + b%thr
  end function SumOp

  function SubtractOp(a, b) result(c)
    class(Op), intent(in) :: a, b
    type(Op) :: c

    c = a
    c%zero = a%zero - b%zero
    c%one = a%one - b%one
    c%two = a%two - b%two
    if(.not. a%is_three_body) return
    c%thr = a%thr - b%thr
  end function SubtractOp

  function ScaleOp(a, b) result(c)
    class(Op), intent(in) :: a
    real(8), intent(in) :: b
    type(Op) :: c

    c = a
    c%zero = a%zero * b
    c%one = a%one * b
    c%two = a%two * b
    if(.not. a%is_three_body) return
    c%thr = a%thr * b
  end function ScaleOp

  subroutine SetOperatorFromFile(this, ms, file_nn, file_3n, &
        & bound_2b_file, bound_3b_file)
    class(Op), intent(inout) :: this
    type(MSpace), intent(in) :: ms
    character(*), intent(in) :: file_nn, file_3n
    integer, intent(in), optional :: bound_2b_file(3), bound_3b_file(4)
    type(NBodyPart) :: Tcm_one, Tcm_two
    type(NBodyPart) :: Hcm_one, Hcm_two
    type(ReadFiles) :: rd

    if(file_nn == 'none' .and. file_3n == 'none') return

    write(*,'(3a)') "Set ", trim(this%optr), " operator from files"
    ! -- boundary for two-body file
    rd%file_nn = file_nn
    rd%emax2 = ms%emax
    rd%e2max2 = ms%e2max
    rd%lmax2 = ms%lmax
    if(present(bound_2b_file)) then
      rd%emax2 = bound_2b_file(1)
      rd%e2max2= bound_2b_file(2)
      rd%lmax2 = bound_2b_file(3)
    end if

    ! -- boundary for three-body file
    rd%file_3n = file_3n
    rd%emax3 = ms%emax
    rd%e2max3 = ms%e2max
    rd%e3max3 = ms%e3max
    rd%lmax3 = ms%lmax
    if(present(bound_3b_file)) then
      rd%emax3 = bound_3b_file(1)
      rd%e2max3= bound_3b_file(2)
      rd%e3max3= bound_3b_file(3)
      rd%lmax3 = bound_3b_file(4)
    end if

    call this%one%set(ms%one, ms%sps, ms%hw, ms%A, ms%Z, ms%N)
    write(*,'(2a)') "2B file: ", trim(file_nn)
    call this%two%ReadFile(ms%sps, ms%two, rd)
    if(this%is_three_body) then
      write(*,'(2a)') "3B file: ", trim(file_3n)
      call this%thr%ReadFile(ms%isps,ms%thr, rd)
    end if

    select case(this%optr)
    case('hamil', 'Hamil')
      call Tcm_one%init(ms%one, .true., 'Tcm', 0, 1, 0)
      call Tcm_two%init(ms%two, .true., 'Tcm', 0, 1, 0)

      call Tcm_one%set(ms%one, ms%sps, ms%hw, ms%A, ms%Z, ms%N)
      call Tcm_two%set(ms%two, ms%sps, ms%hw, ms%A, ms%Z, ms%N)

      this%one = this%one - Tcm_one
      this%two = this%two - Tcm_two

      if(ms%beta > 1.d-4) then
        call Hcm_one%init(ms%one, .true., 'CMhamil', 0, 1, 0)
        call Hcm_two%init(ms%two, .true., 'CMhamil', 0, 1, 0)

        call Hcm_one%set(ms%one, ms%sps, ms%hw, ms%A, ms%Z, ms%N)
        call Hcm_two%set(ms%two, ms%sps, ms%hw, ms%A, ms%Z, ms%N)

        this%zero = this%zero - (1.5d0 * ms%hw * ms%beta)
        this%one = this%one + (Hcm_one * ms%beta)
        this%two = this%two + (Hcm_two * ms%beta)
      end if
    case default
      return
    end select
  end subroutine SetOperatorFromFile

  subroutine SetOperator(this, ms)
    class(Op), intent(inout) :: this
    type(MSpace), intent(in) :: ms

    call this%one%set(ms%one, ms%sps, ms%hw, ms%A, ms%Z, ms%N)
    call this%two%set(ms%two, ms%sps, ms%hw, ms%A, ms%Z, ms%N)
    this%is_three_body = .false.
  end subroutine SetOperator

  subroutine PrintOperator(this, iunit)
    class(Op), intent(in) :: this
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

    if(.not. this%is_three_body) return
    write(ut,'(a)') "## Three-body part"
    call this%thr%prt(iunit)
  end subroutine PrintOperator

  subroutine NormalOrdering(this, ms) ! NO2B approx
    class(Op), intent(inout) :: this
    type(MSpace), intent(in) :: ms
    type(NBodyPart) :: op2from3, op1from3, op1from2
    real(8) :: op0from1, op0from2, op0from3

    if(this%is_normal_ordered) then
      write(*,*) "Operator ", trim(this%optr), " is already normal ordered!"
      return
    end if

    if(this%is_three_body) then
!#ifdef single_precision
      op2from3 = this%thr%NormalOrderingFromSp3To2(ms)
!#else
!      op2from3 = this%thr%NormalOrderingFrom3To2(ms)
!#endif
    end if

    op1from3 = op2from3%NormalOrderingFrom2To1(ms)
    op1from2 = this%two%NormalOrderingFrom2To1(ms)

    op0from3 = op1from3%NormalOrderingFrom1To0(ms)
    op0from2 = op1from2%NormalOrderingFrom1To0(ms)
    op0from1 = this%one%NormalOrderingFrom1To0(ms)

    this%two = this%two + op2from3
    this%one = this%one + op1from2 + op1from3 * 0.5d0
    this%zero = this%zero + op0from1 + op0from2 * 0.5d0 + op0from3 / 6.d0

    this%is_normal_ordered = .true.

  end subroutine NormalOrdering

  subroutine UnNormalOrdering(this, ms) ! From NO2B approx. Op
    class(Op), intent(inout) :: this
    type(MSpace), intent(in) :: ms
    type(NBodyPart) :: op1from2
    real(8) :: op0from1, op0from2

    if(this%is_three_body) then
      write(*,'(a)') "Warning: this unnormal ordering process is only for the NO2B approximated Operators."
    end if

    if(.not. this%is_normal_ordered) then
      write(*,*) "Operator ", trim(this%optr), " is already unnormal ordered!"
      return
    end if

    op1from2 = this%two%NormalOrderingFrom2To1(ms)

    op0from2 = op1from2%NormalOrderingFrom1To0(ms)
    op0from1 = this%one%NormalOrderingFrom1To0(ms)

    ! this%two  is unchanged
    this%one = this%one - op1from2
    this%zero = this%zero - op0from1 + op0from2 * 0.5d0

    this%is_normal_ordered = .false.

  end subroutine UnNormalOrdering

  subroutine DiscardThreeBodyPart(this)
    class(Op), intent(inout) :: this
    if(.not. this%is_three_body) then
      write(*,'(a)') "This operator does not include three-body part."
      write(*,'(a)') "No need to discard three-body part."
      return
    end if
    call this%thr%fin()
    this%is_three_body = .false.
  end subroutine DiscardThreeBodyPart

  subroutine NO2BApproximation(this,ms)
    class(Op), intent(inout) :: this
    type(MSpace), intent(in) :: ms
    if(.not. this%is_three_body) then
      write(*,'(a)') "This operator does not include three-body part."
      write(*,'(a)') "No need to do NO2B approximation."
      return
    end if

    if(this%is_normal_ordered) then
      call this%DiscardThreeBodyPart()
      return
    end if
    call this%NormalOrdering(ms)
    call this%DiscardThreeBodyPart()

  end subroutine NO2BApproximation

end module Operators

!program test
!  use Profiler, only: timer
!  use CommonLibrary, only: &
!      &init_dbinomial_triangle, fin_dbinomial_triangle
!  use ModelSpace, only: MSpace
!  use Operators
!
!  implicit none
!
!  type(MSpace) :: ms
!  type(Op) :: h_me2j
!  type(Op) :: h_myg
!  type(Op) :: h_diff
!  character(:), allocatable :: file_nn, file_3n
!
!  call timer%init()
!  call init_dbinomial_triangle()
!
!  call ms%init('O16', 35.d0, 2, 4, e3max=2)
!  call h_me2j%init('hamil',ms,.false.)
!  call h_myg%init('hamil',ms,.false.)
!
!  file_nn = '/home/takayuki/TwBME-HO_NN-only_N3LO_EM500_bare_hw35_emax6_e2max12.txt.me2j'
!  file_nn = '/home/takayuki/TwBME-HO_NN-only_N3LO_EM500_bare_hw35_emax6_e2max12.bin.me2j'
!  file_3n = 'none'
!  call h_me2j%set(ms,file_nn,file_3n,[6,12,6],[2,2,2,2])
!
!  file_nn = '/home/takayuki/TwBME-HO_NN-only_N3LO_EM500_bare_hw35_emax6_e2max12.txt.myg'
!  file_3n = 'none'
!  call h_myg%set(ms,file_nn,file_3n,[6,12,6],[2,2,2,2])
!
!  h_diff = h_myg - h_me2j
!  call h_diff%prt()
!
!  call h_diff%fin()
!  call h_me2j%fin()
!  call h_myg%fin()
!  call ms%fin()
!
!  call fin_dbinomial_triangle()
!  call timer%fin()
!
!end program test
