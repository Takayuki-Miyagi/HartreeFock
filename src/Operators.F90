module Operators
  use omp_lib
  use myfort
  use ModelSpace
  use OneBodyOperator
  use TwoBodyOperator
  use ThreeBodyOperator
  use ThreeBodyInteraction
  use ThreeBodyMonInteraction
  use ThreeBodyNO2BInteraction
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
  private :: Truncate
  private :: NormalOrdering
  private :: ReNormalOrdering2B
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
    type(ThreeBodyNO2BForce) :: thr21_no2b
    type(ThreeBodyMonForce) :: thr21_mon
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

    procedure :: Truncate
    procedure :: NormalOrdering
    procedure :: ReNormalOrdering2B
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
    if(.not. associated(this%ms)) return
    this%zero = 0.d0
    call this%one%fin()
    call this%two%fin()
    call this%thr21%fin()
    call this%thr%fin()
    call this%thr21_no2b%fin()
    call this%thr21_mon%fin()
    this%ms => null()
    this%is_normal_ordered = .false.
  end subroutine FinOps

  subroutine InitOpsFromString(this, optr, ms, rank, type_3bme)
    use DefineOperators, only: GetOperatorRank
    class(Ops), intent(inout) :: this
    type(MSpace), intent(in) :: ms
    character(*), intent(in) :: optr
    integer, intent(in) :: rank
    character(*), intent(in), optional :: type_3bme
    integer :: jr, pr, zr

    call GetOperatorRank(optr,jr,pr,zr)
    call this%init(jr,pr,zr,optr,ms,rank, type_3bme)
  end subroutine InitOpsFromString

  subroutine InitOps(this, jr, pr, zr, optr, ms, rank, type_3bme)
    class(Ops), intent(inout) :: this
    integer, intent(in) :: jr, pr, zr, rank
    character(*), intent(in) :: optr
    type(MSpace), intent(in), target :: ms
    character(*), intent(in), optional :: type_3bme
    character(:), allocatable :: type_thbme
    real(8) :: ti

    ti = omp_get_wtime()

    type_thbme = 'full'
    if(present(type_3bme)) type_thbme = type_3bme

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

    select case(type_thbme)
    case("monopole", "mon", "Monopole", "Mon")
      if(rank == 3) call this%thr21_mon%init(ms%sps, ms%isps, ms%e2max, ms%e3max)
    case("NO2B", "no2b")
      if(rank == 3) call this%thr21_no2b%init(ms%sps, ms%isps, ms%e2max, ms%e3max)
    case("full", "FULL", "Full")
      if(rank == 3 .and. this%ms%is_three_body_jt) then
        call this%thr21%init(ms%thr21) ! Three-body interaction
      end if

      if(rank == 3 .and. this%ms%is_three_body) then
        call this%thr%init(ms%thr, this%Scalar, optr, jr, pr, zr)
      end if
    case default
      write(*,*) "Unknown three-body matrix treatment"
      return
    end select

    call timer%Add('Construct '//trim(optr), omp_get_wtime()-ti)
  end subroutine InitOps

  subroutine CopyOps(a, b)
    class(Ops), intent(inout) :: a
    type(Ops), intent(in) :: b
    real(8) :: ti

    ti = omp_get_wtime()

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
    call timer%Add('Copy Operator', omp_get_wtime()-ti)
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
    type(Read3BodyMonopole) :: rd3_mon
    type(Read3BodyNO2B) :: rd3_no2b
    type(MSpace), pointer :: ms
    type(sys) :: s
    real(8) :: ti

    if(file_nn == 'none' .and. file_3n == 'none') return
    ms => this%ms
    write(*,'(3a)') "Set ", trim(this%oprtr), " operator from files"
    call rd2%set(file_nn)
    ! -- boundary for two-body file
    call rd2%set(this%ms%emax, this%ms%e2max, this%ms%e3max)
    if(present(bound_2b_file)) then
      call rd2%set(bound_2b_file(1), bound_2b_file(2), bound_2b_file(3))
    end if

    if(.not.s%find(this%oprtr,"file_")) then
      call this%one%set(ms%hw, ms%A, ms%Z, ms%N)
    end if
    write(*,'(2a)') "2B file: ", trim(file_nn)
    call rd2%ReadTwoBodyFile(this%two)

    !------ Three-body file ------
    ! full case
    if(this%rank == 3 .and. this%ms%is_three_body_jt .and. .not. this%thr21%zero) then
      ! -- boundary for three-body file
      call rd3%set(file_3n)
      call rd3%set(ms%emax, ms%e2max, ms%e3max, ms%lmax)
      if(present(bound_3b_file)) then
        call rd3%set(bound_3b_file(1), bound_3b_file(2), bound_3b_file(3), bound_3b_file(4))
      end if
      write(*,'(2a)') "3B file: ", trim(file_3n)
      call rd3%ReadIsospinThreeBodyFile(this%thr21)
      if(this%ms%is_three_body) then
        call this%thr%set(this%thr21)
        call this%DiscardThreeBodyForce()
      end if
    end if

    ! NO2B case
    if(.not. this%thr21_no2b%zero) then
      ! -- boundary for three-body file
      call rd3_no2b%set(file_3n)
      call rd3_no2b%set(ms%emax, ms%e2max, ms%e3max, ms%lmax)
      if(present(bound_3b_file)) then
        call rd3_no2b%set(bound_3b_file(1), bound_3b_file(2), bound_3b_file(3), bound_3b_file(4))
      end if
      write(*,'(2a)') "3B file: ", trim(file_3n)
      call rd3_no2b%ReadThreeBodyNO2B(this%thr21_no2b)
      if(this%ms%is_three_body) then
        write(*,*) "Wrong option in operator from file"
        stop
      end if
    end if

    ! Monopole case
    if(.not. this%thr21_mon%zero) then
      ! -- boundary for three-body file
      call rd3_mon%set(file_3n)
      call rd3_mon%set(ms%emax, ms%e2max, ms%e3max, ms%lmax)
      if(present(bound_3b_file)) then
        call rd3_mon%set(bound_3b_file(1), bound_3b_file(2), bound_3b_file(3), bound_3b_file(4))
      end if
      write(*,'(2a)') "3B file: ", trim(file_3n)
      call rd3_mon%ReadThreeBodyMonopole(this%thr21_mon)
      if(this%ms%is_three_body) then
        write(*,*) "Wrong option in operator from file"
        stop
      end if
    end if

    select case(this%oprtr)
    case('hamil', 'Hamil')
      call Tcm_one%init(ms%one, .true., 'Tcm', 0, 1, 0)
      call Tcm_two%init(ms%two, .true., 'Tcm', 0, 1, 0)

      ti = omp_get_wtime()
      call Tcm_one%set(ms%hw, ms%A, ms%Z, ms%N)
      call Tcm_two%set(ms%hw, ms%A, ms%Z, ms%N)

      this%one = this%one - Tcm_one
      this%two = this%two - Tcm_two
      call timer%add("Set Tcm operator",omp_get_wtime() - ti)

      if(ms%beta > 1.d-4) then
        call Hcm_one%init(ms%one, .true., 'Hcm', 0, 1, 0)
        call Hcm_two%init(ms%two, .true., 'Hcm', 0, 1, 0)

        call Hcm_one%set(ms%hw, ms%A, ms%Z, ms%N)
        call Hcm_two%set(ms%hw, ms%A, ms%Z, ms%N)

        this%zero = this%zero - (1.5d0 * ms%beta)
        this%one = this%one + (Hcm_one * (ms%beta / ms%hw))
        this%two = this%two + (Hcm_two * (ms%beta / ms%hw))
      end if
    case default
      return
    end select
  end subroutine SetOperatorFromFile

  subroutine SetOperator(this)
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

  function Truncate(op, ms_new) result(op_new)
    class(Ops), intent(inout) :: op
    type(MSpace), intent(in), target :: ms_new
    type(Ops) :: op_new
    integer :: i, j, k, l, m, n
    integer :: ii, jj, kk, ll, mm, nn
    type(Orbits), pointer :: sps_new, sps
    type(SingleParticleOrbit), pointer :: oii, ojj, okk, oll, omm, onn
    integer :: chbra_new, chket_new
    integer :: chbra_old, chket_old
    integer :: jbra_new, pbra_new, zbra_new, jket_new, pket_new, zket_new
    integer :: bra, ket, ibra, iket
    real(8) :: ti
    ti = omp_get_wtime()

    sps => op%ms%sps
    sps_new => ms_new%sps
    call op_new%init(op%jr, op%pr, op%zr, op%oprtr, ms_new, op%rank)

    op_new%zero = op%zero

    do ii = 1, sps_new%norbs
      oii => sps_new%GetOrbit(ii)
      do jj = ii, sps_new%norbs
        ojj => sps_new%GetOrbit(jj)

        i = sps%nljz2idx(oii%n, oii%l, oii%j, oii%z)
        j = sps%nljz2idx(ojj%n, ojj%l, ojj%j, ojj%z)
        call op_new%one%SetOBME( ii, jj, op%One%GetOBME(i,j) )

      end do
    end do

    do chbra_new = 1, ms_new%two%NChan
      jbra_new = ms_new%two%jpz(chbra_new)%j
      pbra_new = ms_new%two%jpz(chbra_new)%p
      zbra_new = ms_new%two%jpz(chbra_new)%z
      chbra_old = op%ms%two%jpz2ch(jbra_new,pbra_new,zbra_new)
      do chket_new = 1, ms_new%two%NChan
        jket_new = ms_new%two%jpz(chket_new)%j
        pket_new = ms_new%two%jpz(chket_new)%p
        zket_new = ms_new%two%jpz(chket_new)%z
        chket_old = op%ms%two%jpz2ch(jket_new,pket_new,zket_new)
        if(.not. op%two%MatCh(chbra_old, chket_old)%is) cycle
        if(.not. op_new%two%MatCh(chbra_new, chket_new)%is) cycle

        !$omp parallel
        !$omp do private(bra, ii, jj, ket, kk, ll, i, j, k, l, &
        !$omp &          oii, ojj, okk, oll)
        do bra = 1, ms_new%two%jpz(chbra_new)%n_state
          ii = ms_new%two%jpz(chbra_new)%n2spi1(bra)
          jj = ms_new%two%jpz(chbra_new)%n2spi2(bra)
          oii => sps_new%GetOrbit(ii)
          ojj => sps_new%GetOrbit(jj)
          do ket = 1, ms_new%two%jpz(chket_new)%n_state
            kk = ms_new%two%jpz(chket_new)%n2spi1(ket)
            ll = ms_new%two%jpz(chket_new)%n2spi2(ket)
            okk => sps_new%GetOrbit(kk)
            oll => sps_new%GetOrbit(ll)

            i = sps%nljz2idx(oii%n, oii%l, oii%j, oii%z)
            j = sps%nljz2idx(ojj%n, ojj%l, ojj%j, ojj%z)
            k = sps%nljz2idx(okk%n, okk%l, okk%j, okk%z)
            l = sps%nljz2idx(oll%n, oll%l, oll%j, oll%z)
            op_new%two%MatCh(chbra_new, chket_new)%m(bra, ket) = &
                & op%two%GetTwBME(i, j, k, l, jbra_new, jket_new)
          end do
        end do
        !$omp end do
        !$omp end parallel

      end do
    end do

    do chbra_new = 1, ms_new%thr%NChan
      jbra_new = ms_new%thr%jpz(chbra_new)%j
      pbra_new = ms_new%thr%jpz(chbra_new)%p
      zbra_new = ms_new%thr%jpz(chbra_new)%z
      chbra_old = op%ms%thr%jpz2ch(jbra_new,pbra_new,zbra_new)
      do chket_new = 1, ms_new%thr%NChan
        jket_new = ms_new%thr%jpz(chket_new)%j
        pket_new = ms_new%thr%jpz(chket_new)%p
        zket_new = ms_new%thr%jpz(chket_new)%z
        chket_old = op%ms%thr%jpz2ch(jket_new,pket_new,zket_new)
        if(.not. op%thr%MatCh(chbra_old, chket_old)%is) cycle
        if(.not. op_new%thr%MatCh(chbra_new, chket_new)%is) cycle

        !$omp parallel
        !$omp do private(bra, ii, jj, kk, ket, ll, mm, nn, &
        !$omp &          i, j, k, l, m, n, ibra, iket, &
        !$omp &          oii, ojj, okk, oll, omm, onn)
        do bra = 1, ms_new%thr%jpz(chbra_new)%n_state
          ii =   ms_new%thr%jpz(chbra_new)%n2spi1(bra)
          jj =   ms_new%thr%jpz(chbra_new)%n2spi2(bra)
          kk =   ms_new%thr%jpz(chbra_new)%n2spi3(bra)
          ibra = ms_new%thr%jpz(chbra_new)%n2labl(bra)
          oii => sps_new%GetOrbit(ii)
          ojj => sps_new%GetOrbit(jj)
          okk => sps_new%GetOrbit(kk)

          do ket = 1, ms_new%thr%jpz(chket_new)%n_state
            ll =   ms_new%thr%jpz(chket_new)%n2spi1(ket)
            mm =   ms_new%thr%jpz(chket_new)%n2spi2(ket)
            nn =   ms_new%thr%jpz(chket_new)%n2spi3(ket)
            iket = ms_new%thr%jpz(chket_new)%n2labl(ket)
            oll => sps_new%GetOrbit(ll)
            omm => sps_new%GetOrbit(mm)
            onn => sps_new%GetOrbit(nn)

            i = sps%nljz2idx(oii%n, oii%l, oii%j, oii%z)
            j = sps%nljz2idx(ojj%n, ojj%l, ojj%j, ojj%z)
            k = sps%nljz2idx(okk%n, okk%l, okk%j, okk%z)
            l = sps%nljz2idx(oll%n, oll%l, oll%j, oll%z)
            m = sps%nljz2idx(omm%n, omm%l, omm%j, omm%z)
            n = sps%nljz2idx(onn%n, onn%l, onn%j, onn%z)
            op_new%thr%MatCh(chbra_new, chket_new)%m(bra, ket) = &
                & op%thr%GetThBMEi(i, j, k, ibra, l, m, n, iket, jbra_new, jket_new)
          end do
        end do
        !$omp end do
        !$omp end parallel

      end do
    end do

    call timer%Add('Truncate Operator '//trim(op_new%oprtr), omp_get_wtime()-ti)
  end function Truncate

  function NormalOrdering(this) result(op) ! Normal ordering w.r.t. refernece state
    class(Ops), intent(inout), target :: this
    type(Ops) :: op
    type(TwoBodyPart) :: op2from3
    type(OneBodyPart) :: op1from3, op1from2
    real(8) :: op0from1=0.d0, op0from2=0.d0, op0from3=0.d0

    if(this%is_normal_ordered) then
      write(*,*) "Operator ", trim(this%oprtr), " is already normal ordered!"
      return
    end if
    op = this
    op1from2 = this%two%NormalOrderingFrom2To1(this%ms%one)

    op0from2 = op1from2%NormalOrderingFrom1To0()
    op0from1 = this%one%NormalOrderingFrom1To0()

    op%one = this%one + op1from2
    op%zero = this%zero + op0from1 + op0from2 * 0.5d0

    if(this%rank == 3) then
      if( .not. op%thr21%zero ) then
        op2from3 = op%thr21%NormalOrderingFrom3To2(this%ms%two)
        op1from3 = op2from3%NormalOrderingFrom2To1(this%ms%one)
        op0from3 = op1from3%NormalOrderingFrom1To0()
        op%two = op%two + op2from3
        op%one = op%one + op1from3 * 0.5d0
        op%zero = op%zero + op0from3 / 6.d0
      end if

      if( .not. this%thr21_no2b%zero ) then
        op2from3 = this%thr21_no2b%NormalOrderingFrom3To2(this%ms%two)
        op1from3 = op2from3%NormalOrderingFrom2To1(this%ms%one)
        op0from3 = op1from3%NormalOrderingFrom1To0()
        op%two = op%two + op2from3
        op%one = op%one + op1from3 * 0.5d0
        op%zero = op%zero + op0from3 / 6.d0
      end if

      if(this%ms%is_three_body) then
        op2from3 = op%thr%NormalOrderingFrom3To2(this%ms%two)
        op1from3 = op2from3%NormalOrderingFrom2To1(this%ms%one)
        op0from3 = op1from3%NormalOrderingFrom1To0()
        op%two = op%two + op2from3
        op%one = op%one + op1from3 * 0.5d0
        op%zero = op%zero + op0from3 / 6.d0
      end if
    end if
    op%is_normal_ordered = .true.
  end function NormalOrdering

  function NO2BApprox(this) result(op) ! Normal ordering w.r.t. refernece state
    class(Ops), intent(inout), target :: this
    type(Ops) :: op
    op = this%NormalOrdering()
    if(this%rank==3) then
      call op%DiscardThreeBodyForce()
      call op%DiscardThreeBodyPart()
    end if
    op%rank = 2
    op%is_normal_ordered = .true.
  end function NO2BApprox

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

  subroutine ReNormalOrdering2B(op)
    class(Ops), intent(inout) :: op
    type(MSpace), pointer :: ms
    real(8), allocatable :: NOcoef_orig(:)
    type(SingleParticleOrbit), pointer :: o
    integer :: idx

    ms => op%ms
    allocate(NOcoef_orig(ms%sps%norbs))
    NOcoef_orig = ms%NOcoef
    call op%UnNormalOrdering2B()
    do idx = 1, ms%sps%norbs
      o => ms%sps%GetOrbit(idx)
      if(o%GetCoreValenceOutside() == 1) then
        call o%SetOccupation(0.d0)
      end if
    end do
    op = op%NormalOrdering()

    ms%NOcoef = NOcoef_orig
    do idx = 1, ms%sps%norbs
      o => ms%sps%GetOrbit(idx)
      call o%SetOccupation(ms%NOcoef(idx))
    end do
    deallocate(NOcoef_orig)

  end subroutine ReNormalOrdering2B

end module Operators

!program test
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
