module HartreeFock
  use omp_lib
  use LinAlgLib
  use Operators
  use ThreeBodyMonInteraction
  use ThreeBodyNO2BInteraction
  implicit none

  public :: HFSolver

  private :: FinHFSolver
  private :: InitHFSolver
  private :: SolveHFSolver
  private :: PrintSPEs
  private :: DiagonalizeFockMatrix
  private :: SetOccupationMatrix
  private :: UpdateDensityMatrix
  private :: UpdateDensityMatrixFromCoef
  private :: UpdateFockMatrix
  private :: CalcEnergy
  private :: BasisTransform
  private :: BasisTransNO2BHamiltonian
  private :: BasisTransNO2BHamiltonianFromNO2B
  private :: BasisTransNO2BScalar
  private :: BasisTransNO2BTensor
  private :: BasisTransScalar
  private :: BasisTransTensor


  private :: InitMonopole2
  private :: InitMonopole3
  private :: FinMonopole
  private :: GetIndex2
  private :: GetSpLabels2
  private :: GetIndex3
  private :: GetSpLabels3

  type, private :: Monopole
    integer :: nidx = 0
    real(8), allocatable :: v(:)
    integer(8), allocatable :: idx(:)
    logical :: constructed = .false.
  contains
    procedure :: InitMonopole2
    procedure :: InitMonopole3
    procedure :: InitMonopole3_i
    procedure :: InitMonopole3FromMon
    procedure :: InitMonopole3FromNO2B
    procedure :: FinMonopole
  end type Monopole

  ! Note:
  ! alpha is parameter controlling
  ! density matrix mixing ratio
  ! (has to be 0.d0 < alpha <= 1.d0, alpha = 1.d0 means direct iteration)
  type :: HFSolver
    logical :: is_three_body
    integer :: n_iter_max = 1000
    integer :: rank
    real(8) :: alpha = 1.d0
    real(8) :: tol = 1.d-8
    real(8) :: diff
    real(8) :: e0  ! cm term
    real(8) :: e1  ! kinetic term
    real(8) :: e2  ! nn int. term
    real(8) :: e3  ! 3n int. term
    real(8) :: ehf ! hf energy

    logical :: is_roothaan=.false.
    type(OneBodyPart) :: F   ! Fock opeartor
    type(OneBodyPart) :: T   ! kinetic term
    type(OneBodyPart) :: V   ! one-body filed from 2body interaction
    type(OneBOdyPart) :: W   ! one-body filed form 3body interaction
    type(OneBodyPart) :: rho ! density matrix
    type(OneBodyPart) :: C   ! F's diagonalization coefficient, (HO|HF)
    type(OneBodyPart) :: Occ ! diagonal Occupation matrix
    type(OneBodyPart) :: S   ! norm kernel
    type(Monopole) :: V2, V3
  contains
    procedure :: fin => FinHFSolver
    procedure :: init => InitHFSolver
    procedure :: solve => SolveHFSolver
    procedure :: PrintSPEs
    procedure :: DiagonalizeFockMatrix
    procedure :: SetOccupationMatrix
    procedure :: UpdateDensityMatrix
    procedure :: UpdateDensityMatrixFromCoef
    procedure :: UpdateFockMatrix
    procedure :: CalcEnergy
    procedure :: BasisTransform
    procedure :: BasisTransNO2BHamiltonian
    procedure :: BasisTransNO2BHamiltonianFromNO2B
    procedure :: BasisTransNO2BScalar
    procedure :: BasisTransNO2BTensor
    procedure :: BasisTransScalar
    procedure :: BasisTransTensor
  end type HFSolver

contains

  subroutine FinHFSolver(this)
    class(HFSolver), intent(inout) :: this
    call this%rho%fin()
    call this%F%fin()
    call this%T%fin()
    call this%V%fin()
    call this%W%fin()
    !call this%S%fin()
    call this%V2%FinMonopole()
    call this%V3%FinMonopole()
  end subroutine FinHFSolver

  subroutine InitHFSolver(this,hamil,n_iter_max,tol,alpha,norm_kernel)
    use Profiler, only: timer
    class(HFSolver), intent(inout) :: this
    type(Ops), intent(in) :: hamil
    integer, intent(in), optional :: n_iter_max
    real(8), intent(in), optional :: tol, alpha
    type(OneBodyPart), intent(in), optional :: norm_kernel
    type(MSpace), pointer :: ms
    real(8) :: ti
    integer :: ch, ndim

    ms => hamil%ms
    if(hamil%is_normal_ordered) then
      write(*,'(a)') "Hamiltonian has to be normal ordered w.r.t the vacuum."
      return
    end if

    this%e0 = hamil%zero
    this%e1 = 0.d0
    this%e2 = 0.d0
    this%e3 = 0.d0
    this%diff = 1.d2
    if(present(n_iter_max)) this%n_iter_max = n_iter_max
    if(present(tol)) this%tol = tol
    if(present(alpha)) this%alpha = alpha
    if(present(norm_kernel)) then
      this%is_roothaan = .true.
      this%S = norm_kernel
    end if
    this%is_three_body = hamil%ms%is_three_body_jt .or. hamil%ms%is_three_body
    this%rank = hamil%rank
    call this%C%init(  ms%one, .true., 'UT',      0, 1, 0)
    call this%Occ%init(ms%one, .true., 'Occ',     0, 1, 0)
    call this%rho%init(ms%one, .true., 'DenMat',  0, 1, 0)
    call this%F%init(  ms%one, .true., 'FockOp',  0, 1, 0)
    call this%T%init(  ms%one, .true., 'kinetic', 0, 1, 0)
    call this%V%init(  ms%one, .true., 'NNint',   0, 1, 0)
    call this%W%init(  ms%one, .true., 'NNNint',  0, 1, 0)

    this%T = hamil%one

    ti = omp_get_wtime()
    call timer%cmemory()
    call this%V2%InitMonopole2(hamil%two)
    call timer%countup_memory('Monople 2Body int.')
    call timer%Add("Construct Monopole 2Body int.", omp_get_wtime() - ti)

    if(this%rank == 3 .and. this%is_three_body) then
      ti = omp_get_wtime()
      call timer%cmemory()

      if(.not. hamil%thr21_mon%zero) call this%V3%InitMonopole3FromMon(hamil%thr21_mon)
      if(.not. hamil%thr21_no2b%zero) call this%V3%InitMonopole3FromNO2B(hamil%thr21_no2b)
      if(.not. hamil%thr21%zero) then
        if(hamil%ms%is_three_body_jt) call this%V3%InitMonopole3(hamil%thr21)
        if(hamil%ms%is_three_body) call this%V3%InitMonopole3_i(hamil%thr)
      end if
      call timer%countup_memory('Monople 3Body int.')
      call timer%Add("Construct Monopole 3Body int.", omp_get_wtime() - ti)
    end if

    do ch = 1, this%C%one%NChan
      ndim = this%C%one%jpz(ch)%n_state
      call this%C%MatCh(ch,ch)%eye( ndim )
    end do
    call this%SetOccupationMatrix(ms%NOcoef)
    call this%UpdateDensityMatrix()
    call this%UpdateFockMatrix()

  end subroutine InitHFSolver

  subroutine SolveHFSolver(this)
    use Profiler, only: timer
    class(HFSolver), intent(inout) :: this
    integer :: iter
    real(8) :: ti

    ti = omp_get_wtime()
    write(*,'(2x,a,9x,a,10x,a,9x,a,7x,a,14x,a)') "iter", "zero-body", &
        & "one-body ", "two-body ", &
        & "three-body ", "Ehf"

    call this%CalcEnergy()
    write(*,'(4x,a,5f18.6)') "HO", this%e0, this%e1, &
        &  this%e2, this%e3,this%ehf
    do iter = 1, this%n_iter_max

      call this%DiagonalizeFockMatrix()
      call this%UpdateDensityMatrix()
      call this%UpdateFockMatrix()

      call this%CalcEnergy()

      write(*,'(2x,i4,5f18.6)') iter, this%e0, this%e1, &
          &  this%e2, this%e3,this%ehf

      if(this%tol > this%diff) exit

      if(iter == this%n_iter_max) then
        write(*,'(a,i5,a)') "Hartree-Fock iteration does not converge after ", &
            & this%n_iter_max, " iterations"
        write(*,'(es14.6)') this%diff
      end if
    end do

    call timer%Add("Hartree-Fock iteration",omp_get_wtime() - ti)
  end subroutine SolveHFSolver

  function BasisTransform(HF,Optr,is_NO2B) result(op)
    !  Input: Operator is HO basis operator (not normal ordered)
    ! Output: Operator is HF basis operator {not normal ordered                    , is_NO2B is false
    !                                       {normal ordered with NO2B approximation, is_NO2B is true
    !
    class(HFSolver), intent(inout) :: HF
    type(Ops), intent(in) :: Optr
    type(Ops) :: op
    logical, intent(in), optional :: is_NO2B
    logical :: NO2B=.True.

    if(present(is_NO2B)) NO2B = is_NO2B

    if(Optr%is_normal_ordered) then
      write(*,'(a)') "Normal ordering status of input operator is wrong!"
      stop
    end if

    call HF%UpdateDensityMatrixFromCoef()
    if(NO2B) then
      if(Optr%oprtr=='hamil' .or. Optr%oprtr=='Hamil') then
        call HF%UpdateFockMatrix()
        call HF%CalcEnergy()
        if(Optr%rank < 3) op = HF%BasisTransNO2BHamiltonian(Optr) ! w/o three-body force
        if(.not. Optr%thr21%zero) op = HF%BasisTransNO2BHamiltonian(Optr) ! w/ full three-body force
        if(.not. Optr%thr21_no2b%zero) op = HF%BasisTransNO2BHamiltonianFromNO2B(Optr) ! w/ no2b relevant three-body force
        if(.not. Optr%thr21_mon%zero) then
          write(*,*) "Wrong option, taking normal ordering"
          stop
        end if
        op%is_normal_ordered = .true.
        return
      end if

      if(Optr%Scalar) then
        op = HF%BasisTransNO2BScalar(Optr)
        op%is_normal_ordered = .true.
        return
      end if

      if(.not. Optr%Scalar) then
        write(*,*) "Not tested yet."
        return
        op = HF%BasisTransNO2BTensor(Optr)
        op%is_normal_ordered = .true.
        return
      end if
    end if

    if(Optr%Scalar) then
      op = HF%BasisTransScalar(Optr)
      return
    end if

    if(.not. Optr%Scalar) then
      write(*,*) "Not tested yet."
      return
      op = HF%BasisTransTensor(Optr)
    end if

  end function BasisTransform

  function BasisTransNO2BHamiltonian(HF,H) result(op)
    use Profiler, only: timer
    class(HFSolver), intent(in) :: HF
    type(Ops), intent(in) :: H
    type(Ops) :: op
    type(MSpace), pointer :: ms
    type(TwoBodyChannel), pointer :: ch_two
    type(Orbits), pointer :: sps
    integer :: ch, J, n, bra, ket, a, b, c, d, e, f, JJJ
    integer :: ea, eb, ec, ed
    integer :: je, le, ze, ee
    integer :: jf, lf, zf, ef
    real(8) :: ph, ti
    type(DMat) :: UT, V2, V3

    ti = omp_get_wtime()
    ms => H%ms
    sps => ms%sps
    call op%init(0, 1, 0, "hamil", ms, 2)
    op%zero = HF%ehf
    do ch = 1, ms%one%NChan
      op%one%MatCh(ch,ch)%DMat = HF%C%MatCh(ch,ch)%DMat%T() * &
          &  HF%F%MatCh(ch,ch)%DMat * HF%C%MatCh(ch,ch)%DMat
    end do

    do ch = 1, ms%two%NChan
      ch_two => ms%two%jpz(ch)
      J = ch_two%j
      n = ch_two%n_state
      call UT%zeros(n,n)
      call V3%zeros(n,n)
      V2 = H%two%MatCh(ch,ch)%DMat

      !$omp parallel
      !$omp do private(bra,a,b,ea,eb,ph,ket,c,d,ec,ed,&
      !$omp &  e,je,le,ze,ee,f,jf,lf,zf,ef,JJJ)
      do bra = 1, n
        a = ch_two%n2spi1(bra)
        b = ch_two%n2spi2(bra)
        ea = sps%orb(a)%e
        eb = sps%orb(b)%e
        ph = (-1.d0)**((sps%orb(a)%j+sps%orb(b)%j)/2-J)
        do ket = 1, n
          c = ch_two%n2spi1(ket)
          d = ch_two%n2spi2(ket)
          ec = sps%orb(c)%e
          ed = sps%orb(d)%e

          UT%m(bra,ket) = HF%C%GetOBME(a,c) * HF%C%GetOBME(b,d)
          if(a/=b) UT%m(bra,ket) = UT%m(bra,ket) - ph * &
              & HF%C%GetOBME(a,d) * HF%C%GetOBME(b,c)
          if(a==b) UT%m(bra,ket) = UT%m(bra,ket) * dsqrt(2.d0)
          if(c==d) UT%m(bra,ket) = UT%m(bra,ket) / dsqrt(2.d0)

          if(ket > bra) cycle
          if(H%rank/=3 .or. .not. H%ms%is_three_body_jt) cycle
          do e = 1, sps%norbs
            je = sps%orb(e)%j
            le = sps%orb(e)%l
            ze = sps%orb(e)%z
            ee = sps%orb(e)%e
            if(ea+eb+ee > ms%e3max) cycle
            do f = 1, ms%sps%norbs
              jf = sps%orb(f)%j
              lf = sps%orb(f)%l
              zf = sps%orb(f)%z
              ef = sps%orb(f)%e
              if(je /= jf) cycle
              if(le /= lf) cycle
              if(ze /= zf) cycle
              if(ec+ed+ef > ms%e3max) cycle

              do JJJ = abs(2*J-je), (2*J+je), 2
                V3%m(bra,ket) = V3%m(bra,ket) + HF%rho%GetOBME(e,f) * &
                    & dble(JJJ+1) * H%thr21%GetThBME(a,b,e,J,c,d,f,J,JJJ)
              end do

            end do
          end do
          V3%m(bra,ket) = V3%m(bra,ket) / dble(2*J+1)
          if(a==b) V3%m(bra,ket) = V3%m(bra,ket) / dsqrt(2.d0)
          if(c==d) V3%m(bra,ket) = V3%m(bra,ket) / dsqrt(2.d0)
          V3%m(ket,bra) = V3%m(bra,ket)
        end do
      end do
      !$omp end do
      !$omp end parallel
      op%two%MatCh(ch,ch)%DMat = UT%T() * (V2+V3) * UT

      call UT%fin()
      call V2%fin()
      call V3%fin()
    end do
    call timer%Add("BasisTransNO2BHamltonian",omp_get_wtime() - ti)

  end function BasisTransNO2BHamiltonian

  function BasisTransNO2BHamiltonianFromNO2B(HF,H) result(op)
    use Profiler, only: timer
    class(HFSolver), intent(in) :: HF
    type(Ops), intent(in) :: H
    type(Ops) :: op
    type(MSpace), pointer :: ms
    type(TwoBodyChannel), pointer :: ch_two
    type(Orbits), pointer :: sps
    integer :: ch, J, n, bra, ket, a, b, c, d, e, f
    integer :: ea, eb, ec, ed
    integer :: je, le, ze, ee
    integer :: jf, lf, zf, ef
    real(8) :: ph, ti
    type(DMat) :: UT, V2, V3

    ti = omp_get_wtime()
    ms => H%ms
    sps => ms%sps
    call op%init(0, 1, 0, "hamil", ms, 2)
    op%zero = HF%ehf
    do ch = 1, ms%one%NChan
      op%one%MatCh(ch,ch)%DMat = HF%C%MatCh(ch,ch)%DMat%T() * &
          &  HF%F%MatCh(ch,ch)%DMat * HF%C%MatCh(ch,ch)%DMat
    end do

    do ch = 1, ms%two%NChan
      ch_two => ms%two%jpz(ch)
      J = ch_two%j
      n = ch_two%n_state
      call UT%zeros(n,n)
      call V3%zeros(n,n)
      V2 = H%two%MatCh(ch,ch)%DMat

      !$omp parallel
      !$omp do private(bra,a,b,ea,eb,ph,ket,c,d,ec,ed,&
      !$omp &  e,je,le,ze,ee,f,jf,lf,zf,ef)
      do bra = 1, n
        a = ch_two%n2spi1(bra)
        b = ch_two%n2spi2(bra)
        ea = sps%orb(a)%e
        eb = sps%orb(b)%e
        ph = (-1.d0)**((sps%orb(a)%j+sps%orb(b)%j)/2-J)
        do ket = 1, n
          c = ch_two%n2spi1(ket)
          d = ch_two%n2spi2(ket)
          ec = sps%orb(c)%e
          ed = sps%orb(d)%e

          UT%m(bra,ket) = HF%C%GetOBME(a,c) * HF%C%GetOBME(b,d)
          if(a/=b) UT%m(bra,ket) = UT%m(bra,ket) - ph * &
              & HF%C%GetOBME(a,d) * HF%C%GetOBME(b,c)
          if(a==b) UT%m(bra,ket) = UT%m(bra,ket) * dsqrt(2.d0)
          if(c==d) UT%m(bra,ket) = UT%m(bra,ket) / dsqrt(2.d0)

          if(ket > bra) cycle
          if(H%rank/=3 .or. .not. H%ms%is_three_body_jt) cycle
          do e = 1, sps%norbs
            je = sps%orb(e)%j
            le = sps%orb(e)%l
            ze = sps%orb(e)%z
            ee = sps%orb(e)%e
            if(ea+eb+ee > ms%e3max) cycle
            do f = 1, ms%sps%norbs
              jf = sps%orb(f)%j
              lf = sps%orb(f)%l
              zf = sps%orb(f)%z
              ef = sps%orb(f)%e
              if(je /= jf) cycle
              if(le /= lf) cycle
              if(ze /= zf) cycle
              if(ec+ed+ef > ms%e3max) cycle

              V3%m(bra,ket) = V3%m(bra,ket) + HF%rho%GetOBME(e,f) * &
                  & H%thr21_no2b%GetNO2BThBME(a,b,e,c,d,f,J)

            end do
          end do
          V3%m(bra,ket) = V3%m(bra,ket) / dble(2*J+1)
          if(a==b) V3%m(bra,ket) = V3%m(bra,ket) / dsqrt(2.d0)
          if(c==d) V3%m(bra,ket) = V3%m(bra,ket) / dsqrt(2.d0)
          V3%m(ket,bra) = V3%m(bra,ket)
        end do
      end do
      !$omp end do
      !$omp end parallel
      op%two%MatCh(ch,ch)%DMat = UT%T() * (V2+V3) * UT
      call UT%fin()
      call V2%fin()
      call V3%fin()
    end do
    call timer%Add("BasisTransNO2BHamltonianFromNO2B",omp_get_wtime() - ti)

  end function BasisTransNO2BHamiltonianFromNO2B

  function BasisTransNO2BScalar(HF,Opr) result(op)
    use Profiler, only: timer
    class(HFSolver), intent(in) :: HF
    type(Ops), intent(in) :: Opr
    type(Ops) :: op
    type(MSpace), pointer :: ms
    type(TwoBodyPart) :: o2from3
    type(OneBodyPart) :: o1from3, o1from2
    integer :: ch, J, n, bra, ket, a, b, c, d, e, f, JJJ
    integer :: ea, eb, ec, ed
    integer :: je, le, ze, ee
    integer :: jf, lf, zf, ef
    real(8) :: ph, ti, o0from1, o0from2, o0from3
    type(DMat) :: UT, V2, V3

    ms => Opr%ms
    ti = omp_get_wtime()
    call op%init(0, 1, 0, opr%oprtr, ms, 2)
    op%zero = Opr%zero
    do ch = 1, ms%one%NChan
      Op%one%MatCh(ch,ch)%DMat = HF%C%MatCh(ch,ch)%DMat%T() * &
          &  Opr%one%MatCh(ch,ch)%DMat * HF%C%MatCh(ch,ch)%DMat
    end do

    call o2from3%init(ms%two, .true., Opr%oprtr, 0, 1, 0)
    do ch = 1, ms%two%NChan
      J = ms%two%jpz(ch)%j
      n = ms%two%jpz(ch)%n_state
      call UT%zeros(n,n)
      call V3%zeros(n,n)

      !$omp parallel
      !$omp do private(bra,a,b,ea,eb,ph,ket,c,d,ec,ed,&
      !$omp &  e,je,le,ze,ee,f,jf,lf,zf,ef,JJJ)
      do bra = 1, n
        a = ms%two%jpz(ch)%n2spi1(bra)
        b = ms%two%jpz(ch)%n2spi2(bra)
        ea = ms%sps%orb(a)%e
        eb = ms%sps%orb(b)%e
        ph = (-1.d0)**((ms%sps%orb(a)%j+ms%sps%orb(b)%j)/2-J)
        do ket = 1, n
          c = ms%two%jpz(ch)%n2spi1(ket)
          d = ms%two%jpz(ch)%n2spi2(ket)
          ec = ms%sps%orb(c)%e
          ed = ms%sps%orb(d)%e

          UT%m(bra,ket) = HF%C%GetOBME(a,c) * HF%C%GetOBME(b,d)
          if(a/=b) UT%m(bra,ket) = UT%m(bra,ket) - ph * &
              & HF%C%GetOBME(a,d) * HF%C%GetOBME(b,c)
          if(a==b) UT%m(bra,ket) = UT%m(bra,ket) * dsqrt(2.d0)
          if(c==d) UT%m(bra,ket) = UT%m(bra,ket) / dsqrt(2.d0)

          if(ket > bra) cycle
          if(Opr%rank/=3 .or. .not. Opr%ms%is_three_body_jt) cycle
          do e = 1, ms%sps%norbs
            je = ms%sps%orb(e)%j
            le = ms%sps%orb(e)%l
            ze = ms%sps%orb(e)%z
            ee = ms%sps%orb(e)%e
            if(ea+eb+ee > ms%e3max) cycle
            do f = 1, ms%sps%norbs
              jf = ms%sps%orb(f)%j
              lf = ms%sps%orb(f)%l
              zf = ms%sps%orb(f)%z
              ef = ms%sps%orb(f)%e
              if(je /= jf) cycle
              if(le /= lf) cycle
              if(ze /= zf) cycle
              if(ec+ed+ef > ms%e3max) cycle

              do JJJ = abs(2*J-je), (2*J+je), 2
                o2from3%MatCh(ch,ch)%m(bra,ket) = &
                    & o2from3%MatCh(ch,ch)%m(bra,ket) + &
                    & HF%rho%GetOBME(e,f) * &
                    & dble(JJJ+1) * Opr%thr21%GetThBME(a,b,e,J,c,d,f,J,JJJ)
              end do

            end do
          end do
          o2from3%MatCh(ch,ch)%m(bra,ket) = &
              & o2from3%MatCh(ch,ch)%m(bra,ket) / dble(2*J+1)
          if(a==b) o2from3%MatCh(ch,ch)%m(bra,ket) = &
              &    o2from3%MatCh(ch,ch)%m(bra,ket) / dsqrt(2.d0)
          if(c==d) o2from3%MatCh(ch,ch)%m(bra,ket) = &
              &    o2from3%MatCh(ch,ch)%m(bra,ket) / dsqrt(2.d0)
          o2from3%MatCh(ch,ch)%m(ket,bra) = &
              & o2from3%MatCh(ch,ch)%m(bra,ket)
        end do
      end do
      !$omp end do
      !$omp end parallel
      Op%two%MatCh(ch,ch)%DMat = &
          & UT%T() * Opr%two%MatCh(ch,ch)%DMat * UT
      if(Opr%rank/=3 .or. .not. Opr%ms%is_three_body_jt) o2from3%MatCh(ch,ch)%DMat = &
          & UT%T() * o2from3%MatCh(ch,ch)%DMat * UT
      call UT%fin()
      call V2%fin()
    end do

    o1from3 = o2from3%NormalOrderingFrom2To1(ms%one)
    o1from2 = Op%two%NormalOrderingFrom2To1(ms%one)

    o0from3 = o1from3%NormalOrderingFrom1To0()
    o0from2 = o1from2%NormalOrderingFrom1To0()
    o0from1 = Op%one%NormalOrderingFrom1To0()

    Op%zero = Op%zero + o0from1 + o0from2 * 0.5d0 + o0from3 / 6.d0
    Op%one = Op%one + o1from2 + o1from3 * 0.5d0
    Op%two = Op%two + o2from3

    call timer%Add("BasisTransNO2BScalar",omp_get_wtime() - ti)

  end function BasisTransNO2BScalar

  function BasisTransNO2BTensor(HF,Opr) result(op)
    use Profiler, only: timer
    class(HFSolver), intent(in) :: HF
    type(Ops), intent(in) :: Opr
    type(Ops) :: op
    type(MSpace), pointer :: ms
    type(TwoBodyPart) :: o2from3
    type(OneBodyPart) :: o1from2, o1from3
    integer :: chbra, chket, Jbra, Jket, nbra, nket
    integer :: bra, ket, a, b, c, d
    real(8) :: ph, ti, o0from1, o0from2, o0from3
    type(DMat) :: UTbra, UTket, V2

    ms => Opr%ms
    ti = omp_get_wtime()
    call op%init(opr%jr, opr%pr, opr%zr, opr%oprtr, ms, 2)
    Op%zero = 0.d0
    do chbra = 1, ms%one%NChan
      do chket = 1, ms%one%NChan
        Op%one%MatCh(chbra,chket)%DMat = HF%C%MatCh(chbra,chbra)%DMat%T() * &
            &  Opr%one%MatCh(chbra,chket)%DMat * HF%C%MatCh(chket,chket)%DMat
      end do
    end do

    call o2from3%init(ms%two, .true., Opr%oprtr, 0, 1, 0)
    do chbra = 1, ms%two%NChan
      Jbra = ms%two%jpz(chbra)%j
      nbra = ms%two%jpz(chbra)%n_state
      call UTbra%zeros(nbra,nbra)

      !$omp parallel
      !$omp do private(bra,a,b,ph,ket,c,d)
      do bra = 1, nbra
        a = ms%two%jpz(chbra)%n2spi1(bra)
        b = ms%two%jpz(chbra)%n2spi2(bra)
        ph = (-1.d0)**((ms%sps%orb(a)%j+ms%sps%orb(b)%j)/2-Jbra)
        do ket = 1, nbra
          c = ms%two%jpz(chbra)%n2spi1(ket)
          d = ms%two%jpz(chbra)%n2spi2(ket)

          UTbra%m(bra,ket) = HF%C%GetOBME(a,c) * HF%C%GetOBME(b,d)
          if(a/=b) UTbra%m(bra,ket) = UTbra%m(bra,ket) - ph * &
              & HF%C%GetOBME(a,d) * HF%C%GetOBME(b,c)
          if(a==b) UTbra%m(bra,ket) = UTbra%m(bra,ket) * dsqrt(2.d0)
          if(c==d) UTbra%m(bra,ket) = UTbra%m(bra,ket) / dsqrt(2.d0)
        end do
      end do
      !$omp end do
      !$omp end parallel

      do chket = 1, ms%two%NChan
        Jket = ms%two%jpz(chket)%j
        nket = ms%two%jpz(chket)%n_state
        call UTket%zeros(nket,nket)
        !$omp parallel
        !$omp do private(bra,a,b,ph,ket,c,d)
        do bra = 1, nket
          a = ms%two%jpz(chket)%n2spi1(bra)
          b = ms%two%jpz(chket)%n2spi2(bra)
          ph = (-1.d0)**((ms%sps%orb(a)%j+ms%sps%orb(b)%j)/2-Jket)
          do ket = 1, nket
            c = ms%two%jpz(chket)%n2spi1(ket)
            d = ms%two%jpz(chket)%n2spi2(ket)

            UTket%m(bra,ket) = HF%C%GetOBME(a,c) * HF%C%GetOBME(b,d)
            if(a/=b) UTket%m(bra,ket) = UTket%m(bra,ket) - ph * &
                & HF%C%GetOBME(a,d) * HF%C%GetOBME(b,c)
            if(a==b) UTket%m(bra,ket) = UTket%m(bra,ket) * dsqrt(2.d0)
            if(c==d) UTket%m(bra,ket) = UTket%m(bra,ket) / dsqrt(2.d0)
          end do
        end do
        !$omp end do
        !$omp end parallel
        Op%two%MatCh(chbra,chket)%DMat = &
            & UTbra%T() * Opr%two%MatCh(chbra,chket)%DMat * UTket
        o2from3%MatCh(chbra,chket)%DMat = &
            & UTbra%T() * o2from3%MatCh(chbra,chket)%DMat * UTket
        call UTbra%fin()
        call UTket%fin()
        call V2%fin()
      end do
    end do

    o1from3 = o2from3%NormalOrderingFrom2To1(ms%one)
    o1from2 = Op%two%NormalOrderingFrom2To1(ms%one)

    o0from3 = o1from3%NormalOrderingFrom1To0()
    o0from2 = o1from2%NormalOrderingFrom1To0()
    o0from1 = Op%one%NormalOrderingFrom1To0()

    Op%one = Op%one + o1from2 + o1from3 * 0.5d0
    Op%two = Op%two + o2from3

    call timer%Add("BasisTransNO2BTensor",omp_get_wtime() - ti)
  end function BasisTransNO2BTensor

  function BasisTransScalar(HF,Op) result(Opnew)
    use Profiler, only: timer
    class(HFSolver), intent(in) :: HF
    type(Ops), intent(in) :: Op
    type(Ops) :: Opnew
    type(MSpace), pointer :: ms
    type(TwoBodyChannel), pointer :: ch_two
    type(ThreeBodyChannel), pointer :: ch_thr
    type(Orbits), pointer :: sps
    type(ThreeBodyPartChannel) :: o3
    integer :: ch, J, n, bra, ket, a, b, c, d, e, f, ibra, iket
    real(8) :: ph, ti
    type(DMat) :: UT, V2

    ti = omp_get_wtime()
    ms => Op%ms
    sps => ms%sps
    call opnew%init(0, 1, 0, op%oprtr, ms, 3)
    do ch = 1, ms%one%NChan
      Opnew%one%MatCh(ch,ch)%DMat = HF%C%MatCh(ch,ch)%DMat%T() * &
          &  Op%one%MatCh(ch,ch)%DMat * HF%C%MatCh(ch,ch)%DMat
    end do

    do ch = 1, ms%two%NChan
      ch_two => ms%two%jpz(ch)
      J = ch_two%j
      n = ch_two%n_state
      call UT%zeros(n,n)
      V2 = Op%two%MatCh(ch,ch)%DMat

      !$omp parallel
      !$omp do private(bra,a,b,ph,ket,c,d)
      do bra = 1, n
        a = ch_two%n2spi1(bra)
        b = ch_two%n2spi2(bra)
        ph = (-1.d0)**((sps%orb(a)%j+sps%orb(b)%j)/2-J)
        do ket = 1, n
          c = ch_two%n2spi1(ket)
          d = ch_two%n2spi2(ket)

          UT%m(bra,ket) = HF%C%GetOBME(a,c) * HF%C%GetOBME(b,d)
          if(a/=b) UT%m(bra,ket) = UT%m(bra,ket) - ph * &
              & HF%C%GetOBME(a,d) * HF%C%GetOBME(b,c)
          if(a==b) UT%m(bra,ket) = UT%m(bra,ket) * dsqrt(2.d0)
          if(c==d) UT%m(bra,ket) = UT%m(bra,ket) / dsqrt(2.d0)
        end do
      end do
      !$omp end do
      !$omp end parallel
      Opnew%two%MatCh(ch,ch)%DMat = UT%T() * V2 * UT
      call UT%fin()
      call V2%fin()
    end do

    if(Op%rank < 3) return
    return
    do ch = 1, ms%thr%NChan
      ch_thr => ms%thr%jpz(ch)
      J = ch_thr%j
      n = ch_thr%n_state
      call o3%init(ch_thr,ch_thr)
      !$omp parallel
      !$omp do private(bra,a,b,c,ibra,ket,d,e,f,iket)
      do bra = 1, n
        a = ch_thr%n2spi1(bra)
        b = ch_thr%n2spi2(bra)
        c = ch_thr%n2spi3(bra)
        ibra = ch_thr%n2labl(bra)
        do ket = 1, n
          d = ch_thr%n2spi1(ket)
          e = ch_thr%n2spi2(ket)
          f = ch_thr%n2spi3(ket)
          iket = ch_thr%n2labl(ket)
          o3%m(bra,ket) = HF_ThBME_scalar(HF%C,Op%thr, &
              & a, b, c, ibra, d, e, f, iket, J)
        end do
      end do
      !$omp end do
      !$omp end parallel
      Opnew%thr%MatCh(ch,ch)%DMat = o3%DMat
    end do

    call timer%Add("BasisTransScalar",omp_get_wtime() - ti)
  end function BasisTransScalar

  function BasisTransTensor(HF,Op) result(Opnew)
    use Profiler, only: timer
    class(HFSolver), intent(in) :: HF
    type(Ops), intent(in) :: Op
    type(Ops) :: opnew
    type(MSpace), pointer :: ms
    type(ThreeBodyPartChannel) :: o3
    type(ThreeBodyChannel), pointer :: bra3, ket3
    integer :: chbra, chket, Jbra, Jket, nbra, nket
    integer :: bra, ket, a, b, c, d, e, f, ibra, iket
    real(8) :: ph, ti
    type(DMat) :: UTbra, UTket

    ti = omp_get_wtime()
    ms => Op%ms
    call opnew%init(op%jr, op%pr, op%zr, op%oprtr, ms, 3)
    Opnew%zero = 0.d0
    do chbra = 1, ms%one%NChan
      do chket = 1, ms%one%NChan
        Opnew%one%MatCh(chbra,chket)%DMat = HF%C%MatCh(chbra,chbra)%DMat%T() * &
            &  Op%one%MatCh(chbra,chket)%DMat * HF%C%MatCh(chket,chket)%DMat
      end do
    end do

    do chbra = 1, ms%two%NChan
      Jbra = ms%two%jpz(chbra)%j
      nbra = ms%two%jpz(chbra)%n_state
      call UTbra%zeros(nbra,nbra)

      !$omp parallel
      !$omp do private(bra,a,b,ph,ket,c,d)
      do bra = 1, nbra
        a = ms%two%jpz(chbra)%n2spi1(bra)
        b = ms%two%jpz(chbra)%n2spi2(bra)
        ph = (-1.d0)**((ms%sps%orb(a)%j+ms%sps%orb(b)%j)/2-Jket)
        do ket = 1, nbra
          c = ms%two%jpz(chbra)%n2spi1(ket)
          d = ms%two%jpz(chbra)%n2spi2(ket)

          UTbra%m(bra,ket) = HF%C%GetOBME(a,c) * HF%C%GetOBME(b,d)
          if(a/=b) UTbra%m(bra,ket) = UTbra%m(bra,ket) - ph * &
              & HF%C%GetOBME(a,d) * HF%C%GetOBME(b,c)
          if(a==b) UTbra%m(bra,ket) = UTbra%m(bra,ket) * dsqrt(2.d0)
          if(c==d) UTbra%m(bra,ket) = UTbra%m(bra,ket) / dsqrt(2.d0)
        end do
      end do
      !$omp end do
      !$omp end parallel

      do chket = 1, ms%two%NChan
        Jket = ms%two%jpz(chket)%j
        nket = ms%two%jpz(chket)%n_state
        call UTket%zeros(nket,nket)
        !$omp parallel
        !$omp do private(bra,a,b,ph,ket,c,d)
        do bra = 1, nket
          a = ms%two%jpz(chket)%n2spi1(bra)
          b = ms%two%jpz(chket)%n2spi2(bra)
          ph = (-1.d0)**((ms%sps%orb(a)%j+ms%sps%orb(b)%j)/2-Jket)
          do ket = 1, nket
            c = ms%two%jpz(chket)%n2spi1(ket)
            d = ms%two%jpz(chket)%n2spi2(ket)

            UTket%m(bra,ket) = HF%C%GetOBME(a,c) * HF%C%GetOBME(b,d)
            if(a/=b) UTket%m(bra,ket) = UTket%m(bra,ket) - ph * &
                & HF%C%GetOBME(a,d) * HF%C%GetOBME(b,c)
            if(a==b) UTket%m(bra,ket) = UTket%m(bra,ket) * dsqrt(2.d0)
            if(c==d) UTket%m(bra,ket) = UTket%m(bra,ket) / dsqrt(2.d0)
          end do
        end do
        !$omp end do
        !$omp end parallel
        Opnew%two%MatCh(chbra,chket)%DMat = &
            & UTbra%T() * Op%two%MatCh(chbra,chket)%DMat * UTket
        call UTbra%fin()
        call UTket%fin()
      end do
    end do

    if(Op%rank < 3) return

    do chbra = 1, ms%thr%NChan
      bra3 => ms%thr%jpz(chbra)
      Jbra = bra3%j
      nbra = bra3%n_state
      do chket = 1, ms%thr%NChan
        ket3 => ms%thr%jpz(chket)
        Jket = ket3%j
        nket = ket3%n_state

        call o3%init(bra3, ket3)
        !$omp parallel
        !$omp do private(bra,a,b,c,ibra,ket,d,e,f,iket)
        do bra = 1, nbra
          a = bra3%n2spi1(bra)
          b = bra3%n2spi2(bra)
          c = bra3%n2spi3(bra)
          ibra = bra3%n2labl(bra)
          do ket = 1, nket
            d = ket3%n2spi1(ket)
            e = ket3%n2spi2(ket)
            f = ket3%n2spi3(ket)
            iket = ket3%n2labl(ket)
            o3%m(bra,ket) = HF_ThBME_tensor(HF%C,Op%thr, &
                & a, b, c, ibra, d, e, f, iket, Jbra, Jket)
          end do
        end do
        !$omp end do
        !$omp end parallel
        opnew%thr%MatCh(chbra,chket) = o3
      end do
    end do
    call timer%Add("BasisTransTensor",omp_get_wtime() - ti)

  end function BasisTransTensor

  function HF_ThBME_scalar(CC,thr,a,b,c,ibra,d,e,f,iket,J) result(me)
    type(OneBodyPart), intent(in) :: CC
    type(ThreeBodyPart), intent(in) :: thr
    integer, intent(in) :: a, b, c, d, e, f
    integer, intent(in) :: ibra, iket, J
    integer :: i1, i2, i3, i4, i5, i6
    type(Orbits), pointer :: sps
    real(8) :: me, tr
    sps => CC%one%sps
    me = 0.d0
    do i1 = 1, sps%norbs
      do i2 = 1, sps%norbs
        do i3 = 1, sps%norbs
          do i4 = 1, sps%norbs
            do i5 = 1, sps%norbs
              do i6 = 1, sps%norbs
                tr =  CC%GetOBME(i1,a) * CC%GetOBME(i2,b) * CC%GetOBME(i3,c) * &
                    & CC%GetOBME(i4,d) * CC%GetOBME(i5,e) * CC%GetOBME(i6,f)
                if(abs(tr) < 1.d-10) cycle
                me = me + tr * thr%GetThBMEi(i1,i2,i3,ibra,i4,i5,i6,iket,J)
              end do
            end do
          end do
        end do
      end do
    end do
  end function HF_ThBME_scalar

  ! not completed yet
  function HF_ThBME_tensor(CC,thr,a,b,c,ibra,d,e,f,iket,Jbra,Jket) result(me)
    type(OneBodyPart), intent(in) :: CC
    type(ThreeBodyPart), intent(in) :: thr
    integer, intent(in) :: a, b, c, d, e, f
    integer, intent(in) :: ibra, iket, Jbra, Jket
    integer :: i1, i2, i3, i4, i5, i6
    type(Orbits), pointer :: sps
    real(8) :: me, tr
    sps => CC%one%sps
    me = 0.d0
    do i1 = 1, sps%norbs
      do i2 = 1, sps%norbs
        do i3 = 1, sps%norbs
          do i4 = 1, sps%norbs
            do i5 = 1, sps%norbs
              do i6 = 1, sps%norbs
                tr =  CC%GetOBME(i1,a) * CC%GetOBME(i2,b) * CC%GetOBME(i3,c) * &
                    & CC%GetOBME(i4,d) * CC%GetOBME(i5,e) * CC%GetOBME(i6,f)
                if(abs(tr) < 1.d-10) cycle
                me = me + tr * thr%GetThBMEi(i1,i2,i3,ibra,i4,i5,i6,iket,Jbra,Jket)
              end do
            end do
          end do
        end do
      end do
    end do
  end function HF_ThBME_tensor

  subroutine CalcEnergy(this)
    class(HFSolver), intent(inout) :: this
    type(OneBodySpace), pointer :: one
    integer :: ch, bra, ket, j

    one => this%F%one
    this%e1 = 0.d0
    this%e2 = 0.d0
    this%e3 = 0.d0
    do ch = 1, one%NChan
      j = one%jpz(ch)%j
      do bra = 1, one%jpz(ch)%n_state
        do ket = 1, one%jpz(ch)%n_state
          this%e1 = this%e1 + this%rho%MatCh(ch,ch)%m(bra,ket) * this%T%MatCh(ch,ch)%m(bra,ket) * dble(j+1)
          this%e2 = this%e2 + this%rho%MatCh(ch,ch)%m(bra,ket) * this%V%MatCh(ch,ch)%m(bra,ket) * dble(j+1) * 0.5d0
          this%e3 = this%e3 + this%rho%MatCh(ch,ch)%m(bra,ket) * this%W%MatCh(ch,ch)%m(bra,ket) * dble(j+1) / 6.d0
        end do
      end do
    end do
    this%ehf = this%e0 + this%e1 + this%e2 + this%e3
  end subroutine CalcEnergy

  subroutine DiagonalizeFockMatrix(this)
    class(HFSolver), intent(inout) :: this
    type(EigenSolSymD) :: sol
    type(GenEigenSolSymD) :: sol_gen
    integer :: ch

    if(this%is_roothaan) then
      do ch = 1, this%F%one%NChan
        call sol_gen%init(this%F%MatCh(ch,ch)%DMat, this%S%MatCh(ch,ch)%DMat)
        call sol_gen%DiagSym(this%F%MatCh(ch,ch)%DMat, this%S%MatCh(ch,ch)%DMat)
        this%C%MatCh(ch,ch)%DMat = sol%vec
        call sol%fin()
      end do
      return
    end if

    do ch = 1, this%F%one%NChan
      call sol%init(this%F%MatCh(ch,ch)%DMat)
      call sol%DiagSym(this%F%MatCh(ch,ch)%DMat)
      this%C%MatCh(ch,ch)%DMat = sol%vec
      call sol%fin()
    end do
  end subroutine DiagonalizeFockMatrix

  subroutine SetOccupationMatrix(this,NOcoef)
    class(HFSolver), intent(inout) :: this
    real(8), intent(in) :: NOCoef(:)
    type(OneBodySpace), pointer :: one
    integer :: i, ich, io

    one => this%F%one
    do ich = 1, one%NChan
      do i = 1, one%jpz(ich)%n_state
        io = one%jpz(ich)%n2spi(i)
        this%Occ%MatCh(ich,ich)%m(i,i) = NOcoef(io)
      end do
    end do
  end subroutine SetOccupationMatrix

  subroutine UpdateDensityMatrix(this)
    class(HFSolver), intent(inout) :: this
    integer :: ich
    do ich = 1, this%rho%one%NChan
      this%rho%MatCh(ich,ich)%DMat = &
          & this%rho%MatCh(ich,ich)%DMat * (1.d0-this%alpha) + &
          & (this%C%MatCh(ich,ich)%DMat * this%Occ%MatCh(ich,ich)%DMat * &
          & this%C%MatCh(ich,ich)%DMat%T()) * this%alpha
    end do
  end subroutine UpdateDensityMatrix

  subroutine UpdateDensityMatrixFromCoef(this)
    class(HFSolver), intent(inout) :: this
    integer :: ich
    do ich = 1, this%rho%one%NChan
      this%rho%MatCh(ich,ich)%DMat = &
          & (this%C%MatCh(ich,ich)%DMat * this%Occ%MatCh(ich,ich)%DMat * &
          & this%C%MatCh(ich,ich)%DMat%T())
    end do
  end subroutine UpdateDensityMatrixFromCoef

  subroutine UpdateFockMatrix(this)
    class(HFSolver), intent(inout) :: this
    type(OneBodyPart) :: Fold
    type(Orbits), pointer :: sps
    type(OneBodySpace), pointer :: one
    integer(8) :: i1, i2, i3, i4, i5, i6, idx
    integer :: num, ch1, ch2, ch3
    integer :: ii1, ii2, ii3, ii4, ii5, ii6
    integer :: bra, ket

    one => this%F%one
    sps => one%sps
    Fold = this%F
    do ch1 = 1, one%NChan
      this%V%Match(ch1,ch1)%m = 0.d0
      this%W%Match(ch1,ch1)%m = 0.d0
    end do

    do num = 1, this%V2%nidx
      idx = this%V2%idx(num)
      call GetSpLabels2(idx,i1,i2,i3,i4)
      ch1 = one%jpz2ch(sps%orb(i1)%j, (-1)**sps%orb(i1)%l, sps%orb(i1)%z)
      ch2 = one%jpz2ch(sps%orb(i2)%j, (-1)**sps%orb(i2)%l, sps%orb(i2)%z)
      ii1 = one%jpz(ch1)%spi2n(i1)
      ii2 = one%jpz(ch2)%spi2n(i2)
      ii3 = one%jpz(ch1)%spi2n(i3)
      ii4 = one%jpz(ch2)%spi2n(i4)
      this%V%MatCh(ch1,ch1)%m(ii1,ii3) = this%V%MatCh(ch1,ch1)%m(ii1,ii3) + &
          & this%rho%MatCh(ch2,ch2)%m(ii2,ii4) * this%V2%v(num)
    end do

    if(this%rank == 3 .and. this%is_three_body) then
      do num = 1, this%V3%nidx
        idx = this%V3%idx(num)
        call GetSpLabels3(idx,i1,i2,i3,i4,i5,i6)
        ch1 = one%jpz2ch(sps%orb(i1)%j, (-1)**sps%orb(i1)%l, sps%orb(i1)%z)
        ch2 = one%jpz2ch(sps%orb(i2)%j, (-1)**sps%orb(i2)%l, sps%orb(i2)%z)
        ch3 = one%jpz2ch(sps%orb(i3)%j, (-1)**sps%orb(i3)%l, sps%orb(i3)%z)

        ii1 = one%jpz(ch1)%spi2n(i1)
        ii2 = one%jpz(ch2)%spi2n(i2)
        ii3 = one%jpz(ch3)%spi2n(i3)

        ii4 = one%jpz(ch1)%spi2n(i4)
        ii5 = one%jpz(ch2)%spi2n(i5)
        ii6 = one%jpz(ch3)%spi2n(i6)

        this%W%MatCh(ch1,ch1)%m(ii1,ii4) = this%W%MatCh(ch1,ch1)%m(ii1,ii4) + &
            & this%rho%MatCh(ch2,ch2)%m(ii2,ii5) * &
            & this%rho%MatCh(ch3,ch3)%m(ii3,ii6) * &
            & this%V3%v(num)
      end do
    end if

    do ch1 = 1, one%NChan
      do bra = 1, one%jpz(ch1)%n_state
        do ket = 1, bra
          this%V%MatCh(ch1,ch1)%m(ket,bra) = this%V%MatCh(ch1,ch1)%m(bra,ket)
          this%W%MatCh(ch1,ch1)%m(ket,bra) = this%W%MatCh(ch1,ch1)%m(bra,ket)
        end do
      end do
    end do

    this%diff = 0.d0
    this%F = this%T + this%V + this%W * 0.5d0
    do ch1 = 1, one%NChan
      this%diff = max(maxval(abs(this%F%MatCh(ch1,ch1)%m - Fold%MatCh(ch1,ch1)%m)), this%diff)
    end do
    call Fold%fin()

  end subroutine UpdateFockMatrix

  subroutine PrintSPEs(this,ms)
    class(HFSolver), intent(in) :: this
    type(MSpace), intent(in) :: ms
    integer :: i, io, ch
    type(OneBodyPart) :: F_HF

    F_HF = this%F
    do ch = 1, ms%one%NChan
      F_HF%MatCh(ch,ch)%DMat = this%C%MatCh(ch,ch)%DMat%T() * &
          &  this%F%MatCh(ch,ch)%DMat * this%C%MatCh(ch,ch)%DMat
    end do

    write(*,'(a)') "  Hartree-Fock single-particle energies"

    do i = 1, size(ms%holes)
      io = ms%holes(i)
      write(*,'(a,a10,i4,f12.6)') 'hole:     ', trim(ms%sps%GetLabelFromIndex(io)), io, &
          & F_HF%GetOBME(io,io)
    end do

    do i = 1, size(ms%particles)
      io = ms%particles(i)
      write(*,'(a,a10,i4,f12.6)') 'particle: ', trim(ms%sps%GetLabelFromIndex(io)), io, &
          & F_HF%GetOBME(io,io)
    end do
    call F_HF%fin()
  end subroutine PrintSPEs

  subroutine FinMonopole(this)
    class(Monopole), intent(inout) :: this
    if(.not. this%constructed) return
    deallocate(this%v)
    deallocate(this%idx)
    this%constructed = .false.
  end subroutine FinMonopole

  subroutine InitMonopole2(this, vnn)
    class(Monopole), intent(inout) :: this
    type(TwoBodyPart), intent(in) :: vnn
    type(TwoBodySpace), pointer :: ms
    integer :: n, idx
    integer(8) :: i1, num
    integer(8) :: i2
    integer(8) :: i3
    integer(8) :: i4
    integer :: l1, j1, z1, e1
    integer :: l2, j2, z2, e2
    integer :: l3, j3, z3, e3
    integer :: l4, j4, z4, e4
    integer :: JJ
    real(8) :: v, norm

    if(this%constructed) return
    ms => vnn%two
    n = 0
    do i1 = 1, ms%sps%norbs
      l1 = ms%sps%orb(i1)%l
      j1 = ms%sps%orb(i1)%j
      z1 = ms%sps%orb(i1)%z
      e1 = ms%sps%orb(i1)%e

      do i3 = 1, i1
        l3 = ms%sps%orb(i3)%l
        j3 = ms%sps%orb(i3)%j
        z3 = ms%sps%orb(i3)%z
        e3 = ms%sps%orb(i3)%e
        if(j3 /= j1) cycle
        if((-1)**l3 /= (-1)**l1) cycle
        if(z3 /= z1) cycle

        do i2 = 1, ms%sps%norbs
          l2 = ms%sps%orb(i2)%l
          j2 = ms%sps%orb(i2)%j
          z2 = ms%sps%orb(i2)%z
          e2 = ms%sps%orb(i2)%e

          do i4 = 1, ms%sps%norbs
            l4 = ms%sps%orb(i4)%l
            j4 = ms%sps%orb(i4)%j
            z4 = ms%sps%orb(i4)%z
            e4 = ms%sps%orb(i4)%e
            if(j2 /= j4) cycle
            if((-1)**l2 /= (-1)**l4) cycle
            if(z2 /= z4) cycle

            if(e1 + e2 > ms%e2max) cycle
            if(e3 + e4 > ms%e2max) cycle

            n = n + 1

          end do
        end do
      end do
    end do
    this%nidx = n
    allocate(this%idx(n))
    allocate(this%v(n))

    n = 0
    do i1 = 1, ms%sps%norbs
      l1 = ms%sps%orb(i1)%l
      j1 = ms%sps%orb(i1)%j
      z1 = ms%sps%orb(i1)%z
      e1 = ms%sps%orb(i1)%e

      do i3 = 1, i1
        l3 = ms%sps%orb(i3)%l
        j3 = ms%sps%orb(i3)%j
        z3 = ms%sps%orb(i3)%z
        e3 = ms%sps%orb(i3)%e
        if(j3 /= j1) cycle
        if((-1)**l3 /= (-1)**l1) cycle
        if(z3 /= z1) cycle

        do i2 = 1, ms%sps%norbs
          l2 = ms%sps%orb(i2)%l
          j2 = ms%sps%orb(i2)%j
          z2 = ms%sps%orb(i2)%z
          e2 = ms%sps%orb(i2)%e

          do i4 = 1, ms%sps%norbs
            l4 = ms%sps%orb(i4)%l
            j4 = ms%sps%orb(i4)%j
            z4 = ms%sps%orb(i4)%z
            e4 = ms%sps%orb(i4)%e
            if(j2 /= j4) cycle
            if((-1)**l2 /= (-1)**l4) cycle
            if(z2 /= z4) cycle

            if(e1 + e2 > ms%e2max) cycle
            if(e3 + e4 > ms%e2max) cycle

            n = n + 1
            this%idx(n) = GetIndex2(i1,i2,i3,i4)
          end do
        end do
      end do
    end do

    this%v(:) = 0.d0
    !$omp parallel
    !$omp do private(idx,num,i1,i2,i3,i4,j1,j2,norm,v,JJ)
    do idx = 1, this%nidx
      num = this%idx(idx)
      call GetSpLabels2(num,i1,i2,i3,i4)
      j1 = ms%sps%orb(i1)%j
      j2 = ms%sps%orb(i2)%j

      norm = 1.d0
      if(i1 == i2) norm = norm * dsqrt(2.d0)
      if(i3 == i4) norm = norm * dsqrt(2.d0)
      v = 0.d0
      do JJ = abs(j1-j2)/2, (j1+j2)/2
        if(i1 == i2 .and. mod(JJ,2) == 1) cycle
        if(i3 == i4 .and. mod(JJ,2) == 1) cycle
        v = v + dble(2*JJ+1) * &
            & vnn%GetTwBME(&
            & int(i1,kind(JJ)), int(i2,kind(JJ)), int(i3,kind(JJ)), int(i4,kind(JJ)),JJ)
        ! need to convert integer(8) -> integer(4)
      end do
      this%v(idx) = v * norm / dble(j1+1)
    end do
    !$omp end do
    !$omp end parallel
    this%constructed = .true.
  end subroutine InitMonopole2

  subroutine InitMonopole3(this, v3n)
    use MyLibrary, only: triag
    class(Monopole), intent(inout) :: this
    type(ThreeBodyForce), intent(in) :: v3n
    type(NonOrthIsospinThreeBodySpace), pointer :: ms
    integer :: n, idx, JJ, JJJ
    integer(8) :: i1, i2, i3, i4, i5, i6, num
    integer :: l1, j1, z1, e1
    integer :: l2, j2, z2, e2
    integer :: l3, j3, z3, e3
    integer :: l4, j4, z4, e4
    integer :: l5, j5, z5, e5
    integer :: l6, j6, z6, e6
    real(8) :: v

    if(this%constructed) return
    ms => v3n%thr
    n = 0
    do i1 = 1, ms%sps%norbs
      l1 = ms%sps%orb(i1)%l
      j1 = ms%sps%orb(i1)%j
      z1 = ms%sps%orb(i1)%z
      e1 = ms%sps%orb(i1)%e

      do i4 = 1, i1
        l4 = ms%sps%orb(i4)%l
        j4 = ms%sps%orb(i4)%j
        z4 = ms%sps%orb(i4)%z
        e4 = ms%sps%orb(i4)%e

        if(j1 /= j4) cycle
        if((-1)**l1 /= (-1)**l4) cycle
        if(z1 /= z4) cycle

        do i2 = 1, ms%sps%norbs
          l2 = ms%sps%orb(i2)%l
          j2 = ms%sps%orb(i2)%j
          z2 = ms%sps%orb(i2)%z
          e2 = ms%sps%orb(i2)%e

          do i5 = 1, ms%sps%norbs
            l5 = ms%sps%orb(i5)%l
            j5 = ms%sps%orb(i5)%j
            z5 = ms%sps%orb(i5)%z
            e5 = ms%sps%orb(i5)%e
            if(j2 /= j5) cycle
            if((-1)**l2 /= (-1)**l5) cycle
            if(z2 /= z5) cycle

            if(e1 + e2 > ms%e2max) cycle
            if(e4 + e5 > ms%e2max) cycle

            do i3 = 1, ms%sps%norbs
              l3 = ms%sps%orb(i3)%l
              j3 = ms%sps%orb(i3)%j
              z3 = ms%sps%orb(i3)%z
              e3 = ms%sps%orb(i3)%e
              do i6 = 1, ms%sps%norbs
                l6 = ms%sps%orb(i6)%l
                j6 = ms%sps%orb(i6)%j
                z6 = ms%sps%orb(i6)%z
                e6 = ms%sps%orb(i6)%e
                if(j3 /= j6) cycle
                if((-1)**l3 /= (-1)**l6) cycle
                if(z3 /= z6) cycle

                if(e1+e3 > ms%e2max) cycle
                if(e2+e3 > ms%e2max) cycle
                if(e4+e6 > ms%e2max) cycle
                if(e5+e6 > ms%e2max) cycle
                if(e1+e2+e3 > ms%e3max) cycle
                if(e4+e5+e6 > ms%e3max) cycle

                n = n + 1
              end do
            end do

          end do
        end do
      end do
    end do

    this%nidx = n
    allocate(this%idx(n))
    allocate(this%v(n))


    n = 0
    do i1 = 1, ms%sps%norbs
      l1 = ms%sps%orb(i1)%l
      j1 = ms%sps%orb(i1)%j
      z1 = ms%sps%orb(i1)%z
      e1 = ms%sps%orb(i1)%e

      do i4 = 1, i1
        l4 = ms%sps%orb(i4)%l
        j4 = ms%sps%orb(i4)%j
        z4 = ms%sps%orb(i4)%z
        e4 = ms%sps%orb(i4)%e

        if(j1 /= j4) cycle
        if((-1)**l1 /= (-1)**l4) cycle
        if(z1 /= z4) cycle

        do i2 = 1, ms%sps%norbs
          l2 = ms%sps%orb(i2)%l
          j2 = ms%sps%orb(i2)%j
          z2 = ms%sps%orb(i2)%z
          e2 = ms%sps%orb(i2)%e

          do i5 = 1, ms%sps%norbs
            l5 = ms%sps%orb(i5)%l
            j5 = ms%sps%orb(i5)%j
            z5 = ms%sps%orb(i5)%z
            e5 = ms%sps%orb(i5)%e
            if(j2 /= j5) cycle
            if((-1)**l2 /= (-1)**l5) cycle
            if(z2 /= z5) cycle

            if(e1 + e2 > ms%e2max) cycle
            if(e4 + e5> ms%e2max) cycle

            do i3 = 1, ms%sps%norbs
              l3 = ms%sps%orb(i3)%l
              j3 = ms%sps%orb(i3)%j
              z3 = ms%sps%orb(i3)%z
              e3 = ms%sps%orb(i3)%e
              do i6 = 1, ms%sps%norbs
                l6 = ms%sps%orb(i6)%l
                j6 = ms%sps%orb(i6)%j
                z6 = ms%sps%orb(i6)%z
                e6 = ms%sps%orb(i6)%e
                if(j3 /= j6) cycle
                if((-1)**l3 /= (-1)**l6) cycle
                if(z3 /= z6) cycle

                if(e1+e3 > ms%e2max) cycle
                if(e2+e3 > ms%e2max) cycle
                if(e4+e6 > ms%e2max) cycle
                if(e5+e6 > ms%e2max) cycle
                if(e1+e2+e3 > ms%e3max) cycle
                if(e4+e5+e6 > ms%e3max) cycle

                n = n + 1
                this%idx(n) = GetIndex3(i1,i2,i3,i4,i5,i6)
              end do
            end do

          end do
        end do

      end do
    end do

    !$omp parallel
    !$omp do private(idx,num,i1,i2,i3,i4,i5,i6,j1,j2,j3,v,JJ,JJJ)
    do idx = 1, this%nidx
      num = this%idx(idx)
      call GetSpLabels3(num,i1,i2,i3,i4,i5,i6)
      j1 = ms%sps%orb(i1)%j
      j2 = ms%sps%orb(i2)%j
      j3 = ms%sps%orb(i3)%j
      v = 0.d0
      do JJ = abs(j1-j2)/2, (j1+j2)/2
        if(i1 == i2 .and. mod(JJ,2) == 1) cycle
        if(i4 == i5 .and. mod(JJ,2) == 1) cycle
        do JJJ = abs(2*JJ-j3), (2*JJ+j3), 2
          v = v + dble(JJJ+1) * &
              & v3n%GetThBME(&
              & int(i1,kind(JJ)),int(i2,kind(JJ)),int(i3,kind(JJ)),JJ,&
              & int(i4,kind(JJ)),int(i5,kind(JJ)),int(i6,kind(JJ)),JJ,JJJ)
          ! need to convert integer(8) -> integer(4)
          !write(*,'(8i4,f12.6)') i1,i2,i3,i4,i5,i6,JJ,JJJ,v3n%GetThBME(&
          !    & int(i1,kind(JJ)),int(i2,kind(JJ)),int(i3,kind(JJ)),JJ,&
          !    & int(i4,kind(JJ)),int(i5,kind(JJ)),int(i6,kind(JJ)),JJ,JJJ)
        end do
      end do
      this%v(idx) = v / dble(j1+1)
    end do
    !$omp end do
    !$omp end parallel
    this%constructed = .true.
  end subroutine InitMonopole3

  subroutine InitMonopole3_i(this, v3n)
    use MyLibrary, only: triag
    class(Monopole), intent(inout) :: this
    type(ThreeBodyPart), intent(in) :: v3n
    type(ThreeBodySpace), pointer :: ms
    integer :: n, idx, JJ, JJJ
    integer(8) :: i1, i2, i3, i4, i5, i6, num
    integer :: l1, j1, z1, e1
    integer :: l2, j2, z2, e2
    integer :: l3, j3, z3, e3
    integer :: l4, j4, z4, e4
    integer :: l5, j5, z5, e5
    integer :: l6, j6, z6, e6
    real(8) :: v

    if(this%constructed) return
    ms => v3n%thr
    n = 0
    do i1 = 1, ms%sps%norbs
      l1 = ms%sps%orb(i1)%l
      j1 = ms%sps%orb(i1)%j
      z1 = ms%sps%orb(i1)%z
      e1 = ms%sps%orb(i1)%e

      do i4 = 1, i1
        l4 = ms%sps%orb(i4)%l
        j4 = ms%sps%orb(i4)%j
        z4 = ms%sps%orb(i4)%z
        e4 = ms%sps%orb(i4)%e

        if(j1 /= j4) cycle
        if((-1)**l1 /= (-1)**l4) cycle
        if(z1 /= z4) cycle

        do i2 = 1, ms%sps%norbs
          l2 = ms%sps%orb(i2)%l
          j2 = ms%sps%orb(i2)%j
          z2 = ms%sps%orb(i2)%z
          e2 = ms%sps%orb(i2)%e

          do i5 = 1, ms%sps%norbs
            l5 = ms%sps%orb(i5)%l
            j5 = ms%sps%orb(i5)%j
            z5 = ms%sps%orb(i5)%z
            e5 = ms%sps%orb(i5)%e
            if(j2 /= j5) cycle
            if((-1)**l2 /= (-1)**l5) cycle
            if(z2 /= z5) cycle

            if(e1 + e2 > ms%e2max) cycle
            if(e4 + e5 > ms%e2max) cycle

            do i3 = 1, ms%sps%norbs
              l3 = ms%sps%orb(i3)%l
              j3 = ms%sps%orb(i3)%j
              z3 = ms%sps%orb(i3)%z
              e3 = ms%sps%orb(i3)%e
              do i6 = 1, ms%sps%norbs
                l6 = ms%sps%orb(i6)%l
                j6 = ms%sps%orb(i6)%j
                z6 = ms%sps%orb(i6)%z
                e6 = ms%sps%orb(i6)%e
                if(j3 /= j6) cycle
                if((-1)**l3 /= (-1)**l6) cycle
                if(z3 /= z6) cycle

                if(e1+e3 > ms%e2max) cycle
                if(e2+e3 > ms%e2max) cycle
                if(e4+e6 > ms%e2max) cycle
                if(e5+e6 > ms%e2max) cycle
                if(e1+e2+e3 > ms%e3max) cycle
                if(e4+e5+e6 > ms%e3max) cycle

                n = n + 1
              end do
            end do

          end do
        end do
      end do
    end do

    this%nidx = n
    allocate(this%idx(n))
    allocate(this%v(n))


    n = 0
    do i1 = 1, ms%sps%norbs
      l1 = ms%sps%orb(i1)%l
      j1 = ms%sps%orb(i1)%j
      z1 = ms%sps%orb(i1)%z
      e1 = ms%sps%orb(i1)%e

      do i4 = 1, i1
        l4 = ms%sps%orb(i4)%l
        j4 = ms%sps%orb(i4)%j
        z4 = ms%sps%orb(i4)%z
        e4 = ms%sps%orb(i4)%e

        if(j1 /= j4) cycle
        if((-1)**l1 /= (-1)**l4) cycle
        if(z1 /= z4) cycle

        do i2 = 1, ms%sps%norbs
          l2 = ms%sps%orb(i2)%l
          j2 = ms%sps%orb(i2)%j
          z2 = ms%sps%orb(i2)%z
          e2 = ms%sps%orb(i2)%e

          do i5 = 1, ms%sps%norbs
            l5 = ms%sps%orb(i5)%l
            j5 = ms%sps%orb(i5)%j
            z5 = ms%sps%orb(i5)%z
            e5 = ms%sps%orb(i5)%e
            if(j2 /= j5) cycle
            if((-1)**l2 /= (-1)**l5) cycle
            if(z2 /= z5) cycle

            if(e1 + e2 > ms%e2max) cycle
            if(e4 + e5 > ms%e2max) cycle

            do i3 = 1, ms%sps%norbs
              l3 = ms%sps%orb(i3)%l
              j3 = ms%sps%orb(i3)%j
              z3 = ms%sps%orb(i3)%z
              e3 = ms%sps%orb(i3)%e
              do i6 = 1, ms%sps%norbs
                l6 = ms%sps%orb(i6)%l
                j6 = ms%sps%orb(i6)%j
                z6 = ms%sps%orb(i6)%z
                e6 = ms%sps%orb(i6)%e
                if(j3 /= j6) cycle
                if((-1)**l3 /= (-1)**l6) cycle
                if(z3 /= z6) cycle

                if(e1+e3 > ms%e2max) cycle
                if(e2+e3 > ms%e2max) cycle
                if(e4+e6 > ms%e2max) cycle
                if(e5+e6 > ms%e2max) cycle
                if(e1+e2+e3 > ms%e3max) cycle
                if(e4+e5+e6 > ms%e3max) cycle

                n = n + 1
                this%idx(n) = GetIndex3(i1,i2,i3,i4,i5,i6)
              end do
            end do

          end do
        end do

      end do
    end do

    !$omp parallel
    !$omp do private(idx,num,i1,i2,i3,i4,i5,i6,j1,j2,j3,v,JJ,JJJ)
    do idx = 1, this%nidx
      num = this%idx(idx)
      call GetSpLabels3(num,i1,i2,i3,i4,i5,i6)
      j1 = ms%sps%orb(i1)%j
      j2 = ms%sps%orb(i2)%j
      j3 = ms%sps%orb(i3)%j
      v = 0.d0
      do JJ = abs(j1-j2)/2, (j1+j2)/2
        if(i1 == i2 .and. mod(JJ,2) == 1) cycle
        if(i4 == i5 .and. mod(JJ,2) == 1) cycle
        do JJJ = abs(2*JJ-j3), (2*JJ+j3), 2
          v = v + dble(JJJ+1) * &
              & v3n%GetThBMEJ(&
              & int(i1,kind(JJ)),int(i2,kind(JJ)),int(i3,kind(JJ)),JJ,&
              & int(i4,kind(JJ)),int(i5,kind(JJ)),int(i6,kind(JJ)),JJ,JJJ)
          ! need to convert integer(8) -> integer(4)
          !write(*,'(8i4,f12.6)') i1,i2,i3,i4,i5,i6,JJ,JJJ,v
        end do
      end do
      this%v(idx) = v / dble(j1+1)
    end do
    !$omp end do
    !$omp end parallel
    this%constructed = .true.
  end subroutine InitMonopole3_i

  subroutine InitMonopole3FromMon(this, v3n)
    use MyLibrary, only: triag
    class(Monopole), intent(inout) :: this
    type(ThreeBodyMonForce), intent(in), target :: v3n
    type(ThreeBodyMonSpace), pointer :: ms
    type(Orbits), pointer :: sps
    integer :: n, idx
    integer(8) :: i1, i2, i3, i4, i5, i6, num
    integer :: l1, j1, z1, e1
    integer :: l2, j2, z2, e2
    integer :: l3, j3, z3, e3
    integer :: l4, j4, z4, e4
    integer :: l5, j5, z5, e5
    integer :: l6, j6, z6, e6
    real(8) :: v

    if(this%constructed) return
    write(*,"(a)") " From Monopole file"
    ms => v3n%thr
    sps => v3n%sps
    n = 0
    do i1 = 1, sps%norbs
      l1 = sps%orb(i1)%l
      j1 = sps%orb(i1)%j
      z1 = sps%orb(i1)%z
      e1 = sps%orb(i1)%e

      do i4 = 1, i1
        l4 = sps%orb(i4)%l
        j4 = sps%orb(i4)%j
        z4 = sps%orb(i4)%z
        e4 = sps%orb(i4)%e

        if(j1 /= j4) cycle
        if((-1)**l1 /= (-1)**l4) cycle
        if(z1 /= z4) cycle

        do i2 = 1, sps%norbs
          l2 = sps%orb(i2)%l
          j2 = sps%orb(i2)%j
          z2 = sps%orb(i2)%z
          e2 = sps%orb(i2)%e

          do i5 = 1, sps%norbs
            l5 = sps%orb(i5)%l
            j5 = sps%orb(i5)%j
            z5 = sps%orb(i5)%z
            e5 = sps%orb(i5)%e
            if(j2 /= j5) cycle
            if((-1)**l2 /= (-1)**l5) cycle
            if(z2 /= z5) cycle

            if(e1 + e2 > ms%e2max) cycle
            if(e4 + e5 > ms%e2max) cycle

            do i3 = 1, sps%norbs
              l3 = sps%orb(i3)%l
              j3 = sps%orb(i3)%j
              z3 = sps%orb(i3)%z
              e3 = sps%orb(i3)%e
              do i6 = 1, sps%norbs
                l6 = sps%orb(i6)%l
                j6 = sps%orb(i6)%j
                z6 = sps%orb(i6)%z
                e6 = sps%orb(i6)%e
                if(j3 /= j6) cycle
                if((-1)**l3 /= (-1)**l6) cycle
                if(z3 /= z6) cycle

                if(e1+e3 > ms%e2max) cycle
                if(e2+e3 > ms%e2max) cycle
                if(e4+e6 > ms%e2max) cycle
                if(e5+e6 > ms%e2max) cycle
                if(e1+e2+e3 > ms%e3max) cycle
                if(e4+e5+e6 > ms%e3max) cycle

                n = n + 1
              end do
            end do

          end do
        end do
      end do
    end do

    this%nidx = n
    allocate(this%idx(n))
    allocate(this%v(n))


    n = 0
    do i1 = 1, sps%norbs
      l1 = sps%orb(i1)%l
      j1 = sps%orb(i1)%j
      z1 = sps%orb(i1)%z
      e1 = sps%orb(i1)%e

      do i4 = 1, i1
        l4 = sps%orb(i4)%l
        j4 = sps%orb(i4)%j
        z4 = sps%orb(i4)%z
        e4 = sps%orb(i4)%e

        if(j1 /= j4) cycle
        if((-1)**l1 /= (-1)**l4) cycle
        if(z1 /= z4) cycle

        do i2 = 1, sps%norbs
          l2 = sps%orb(i2)%l
          j2 = sps%orb(i2)%j
          z2 = sps%orb(i2)%z
          e2 = sps%orb(i2)%e

          do i5 = 1, sps%norbs
            l5 = sps%orb(i5)%l
            j5 = sps%orb(i5)%j
            z5 = sps%orb(i5)%z
            e5 = sps%orb(i5)%e
            if(j2 /= j5) cycle
            if((-1)**l2 /= (-1)**l5) cycle
            if(z2 /= z5) cycle

            if(e1 + e2 > ms%e2max) cycle
            if(e4 + e5 > ms%e2max) cycle

            do i3 = 1, sps%norbs
              l3 = sps%orb(i3)%l
              j3 = sps%orb(i3)%j
              z3 = sps%orb(i3)%z
              e3 = sps%orb(i3)%e
              do i6 = 1, sps%norbs
                l6 = sps%orb(i6)%l
                j6 = sps%orb(i6)%j
                z6 = sps%orb(i6)%z
                e6 = sps%orb(i6)%e
                if(j3 /= j6) cycle
                if((-1)**l3 /= (-1)**l6) cycle
                if(z3 /= z6) cycle

                if(e1+e3 > ms%e2max) cycle
                if(e2+e3 > ms%e2max) cycle
                if(e4+e6 > ms%e2max) cycle
                if(e5+e6 > ms%e2max) cycle
                if(e1+e2+e3 > ms%e3max) cycle
                if(e4+e5+e6 > ms%e3max) cycle

                n = n + 1
                this%idx(n) = GetIndex3(i1,i2,i3,i4,i5,i6)
              end do
            end do

          end do
        end do

      end do
    end do

    !$omp parallel
    !$omp do private(idx,num,i1,i2,i3,i4,i5,i6,v)
    do idx = 1, this%nidx
      num = this%idx(idx)
      call GetSpLabels3(num,i1,i2,i3,i4,i5,i6)
      v = v3n%GetMonThBME(&
          & int(i1,kind(idx)),int(i2,kind(idx)),int(i3,kind(idx)), &
          & int(i4,kind(idx)),int(i5,kind(idx)),int(i6,kind(idx)))
      this%v(idx) = v / dble(sps%orb(i1)%j+1)
    end do
    !$omp end do
    !$omp end parallel
    this%constructed = .true.
  end subroutine InitMonopole3FromMon

  subroutine InitMonopole3FromNO2B(this, v3n)
    use MyLibrary, only: triag
    class(Monopole), intent(inout) :: this
    type(ThreeBodyNO2BForce), intent(in), target :: v3n
    type(ThreeBodyNO2BSpace), pointer :: ms
    type(Orbits), pointer :: sps
    integer :: n, idx
    integer(8) :: i1, i2, i3, i4, i5, i6, num
    integer :: l1, j1, z1, e1
    integer :: l2, j2, z2, e2
    integer :: l3, j3, z3, e3
    integer :: l4, j4, z4, e4
    integer :: l5, j5, z5, e5
    integer :: l6, j6, z6, e6, J
    real(8) :: v

    if(this%constructed) return
    write(*,"(a)") " From NO2B file"
    ms => v3n%thr
    sps => v3n%sps
    n = 0
    do i1 = 1, sps%norbs
      l1 = sps%orb(i1)%l
      j1 = sps%orb(i1)%j
      z1 = sps%orb(i1)%z
      e1 = sps%orb(i1)%e

      do i4 = 1, i1
        l4 = sps%orb(i4)%l
        j4 = sps%orb(i4)%j
        z4 = sps%orb(i4)%z
        e4 = sps%orb(i4)%e

        if(j1 /= j4) cycle
        if((-1)**l1 /= (-1)**l4) cycle
        if(z1 /= z4) cycle

        do i2 = 1, sps%norbs
          l2 = sps%orb(i2)%l
          j2 = sps%orb(i2)%j
          z2 = sps%orb(i2)%z
          e2 = sps%orb(i2)%e

          do i5 = 1, sps%norbs
            l5 = sps%orb(i5)%l
            j5 = sps%orb(i5)%j
            z5 = sps%orb(i5)%z
            e5 = sps%orb(i5)%e
            if(j2 /= j5) cycle
            if((-1)**l2 /= (-1)**l5) cycle
            if(z2 /= z5) cycle

            if(e1 + e2 > ms%e2max) cycle
            if(e4 + e5 > ms%e2max) cycle

            do i3 = 1, sps%norbs
              l3 = sps%orb(i3)%l
              j3 = sps%orb(i3)%j
              z3 = sps%orb(i3)%z
              e3 = sps%orb(i3)%e
              do i6 = 1, sps%norbs
                l6 = sps%orb(i6)%l
                j6 = sps%orb(i6)%j
                z6 = sps%orb(i6)%z
                e6 = sps%orb(i6)%e
                if(j3 /= j6) cycle
                if((-1)**l3 /= (-1)**l6) cycle
                if(z3 /= z6) cycle

                if(e1+e3 > ms%e2max) cycle
                if(e2+e3 > ms%e2max) cycle
                if(e4+e6 > ms%e2max) cycle
                if(e5+e6 > ms%e2max) cycle
                if(e1+e2+e3 > ms%e3max) cycle
                if(e4+e5+e6 > ms%e3max) cycle

                n = n + 1
              end do
            end do

          end do
        end do
      end do
    end do

    this%nidx = n
    allocate(this%idx(n))
    allocate(this%v(n))


    n = 0
    do i1 = 1, sps%norbs
      l1 = sps%orb(i1)%l
      j1 = sps%orb(i1)%j
      z1 = sps%orb(i1)%z
      e1 = sps%orb(i1)%e

      do i4 = 1, i1
        l4 = sps%orb(i4)%l
        j4 = sps%orb(i4)%j
        z4 = sps%orb(i4)%z
        e4 = sps%orb(i4)%e

        if(j1 /= j4) cycle
        if((-1)**l1 /= (-1)**l4) cycle
        if(z1 /= z4) cycle

        do i2 = 1, sps%norbs
          l2 = sps%orb(i2)%l
          j2 = sps%orb(i2)%j
          z2 = sps%orb(i2)%z
          e2 = sps%orb(i2)%e

          do i5 = 1, sps%norbs
            l5 = sps%orb(i5)%l
            j5 = sps%orb(i5)%j
            z5 = sps%orb(i5)%z
            e5 = sps%orb(i5)%e
            if(j2 /= j5) cycle
            if((-1)**l2 /= (-1)**l5) cycle
            if(z2 /= z5) cycle

            if(e1 + e2 > ms%e2max) cycle
            if(e4 + e5 > ms%e2max) cycle

            do i3 = 1, sps%norbs
              l3 = sps%orb(i3)%l
              j3 = sps%orb(i3)%j
              z3 = sps%orb(i3)%z
              e3 = sps%orb(i3)%e
              do i6 = 1, sps%norbs
                l6 = sps%orb(i6)%l
                j6 = sps%orb(i6)%j
                z6 = sps%orb(i6)%z
                e6 = sps%orb(i6)%e
                if(j3 /= j6) cycle
                if((-1)**l3 /= (-1)**l6) cycle
                if(z3 /= z6) cycle

                if(e1+e3 > ms%e2max) cycle
                if(e2+e3 > ms%e2max) cycle
                if(e4+e6 > ms%e2max) cycle
                if(e5+e6 > ms%e2max) cycle
                if(e1+e2+e3 > ms%e3max) cycle
                if(e4+e5+e6 > ms%e3max) cycle

                n = n + 1
                this%idx(n) = GetIndex3(i1,i2,i3,i4,i5,i6)
              end do
            end do

          end do
        end do

      end do
    end do

    !$omp parallel
    !$omp do private(idx,num,i1,i2,i3,i4,i5,i6,j1,j2,J,v)
    do idx = 1, this%nidx
      num = this%idx(idx)
      call GetSpLabels3(num,i1,i2,i3,i4,i5,i6)
      j1 = sps%orb(i1)%j
      j2 = sps%orb(i2)%j
      v = 0.d0
      do J = abs(j1-j2)/2, (j1+j2)/2
        if(i1 == i2 .and. mod(J,2) == 1) cycle
        if(i4 == i5 .and. mod(J,2) == 1) cycle
        v = v + v3n%GetNO2BThBME(&
            & int(i1,kind(idx)),int(i2,kind(idx)),int(i3,kind(idx)), &
            & int(i4,kind(idx)),int(i5,kind(idx)),int(i6,kind(idx)),J)
      end do
      this%v(idx) = v / dble(sps%orb(i1)%j+1)
    end do
    !$omp end do
    !$omp end parallel
    this%constructed = .true.
  end subroutine InitMonopole3FromNO2B

  function GetIndex2(i1,i2,i3,i4) result(r)
    integer(8), intent(in) :: i1, i2, i3, i4
    integer(8) :: r
    !r = i1 + lshift(i2,10) + lshift(i3,20) + lshift(i4,30)
    r = i1 + shiftl(i2,10) + shiftl(i3,20) + shiftl(i4,30)
  end function GetIndex2

  subroutine GetSpLabels2(idx,i1,i2,i3,i4)
    integer(8), intent(in)  :: idx
    integer(8), intent(out) :: i1, i2, i3, i4
    i1 = mod(idx, int(1024,kind=8))
    i2 = mod(shiftr(idx,10), int(1024,kind=8))
    i3 = mod(shiftr(idx,20), int(1024,kind=8))
    i4 = mod(shiftr(idx,30), int(1024,kind=8))
  end subroutine GetSpLabels2

  function GetIndex3(i1,i2,i3,i4,i5,i6) result(r)
    integer(8), intent(in) :: i1, i2, i3, i4, i5, i6
    integer(8) :: r
    r = i1 + shiftl(i2,10) + shiftl(i3,20) + shiftl(i4,30) + &
        &    shiftl(i5,40) + shiftl(i6,50)
  end function GetIndex3

  subroutine GetSpLabels3(idx,i1,i2,i3,i4,i5,i6)
    integer(8), intent(in)  :: idx
    integer(8), intent(out) :: i1, i2, i3, i4, i5, i6
    i1 = mod(idx, int(1024,kind=8))
    i2 = mod(shiftr(idx,10), int(1024,kind=8))
    i3 = mod(shiftr(idx,20), int(1024,kind=8))
    i4 = mod(shiftr(idx,30), int(1024,kind=8))
    i5 = mod(shiftr(idx,40), int(1024,kind=8))
    i6 = mod(shiftr(idx,50), int(1024,kind=8))
  end subroutine GetSpLabels3
end module HartreeFock

!program test
!  use Profiler, only: timer
!  use Operators
!  use HartreeFock
!
!  implicit none
!
!  type(MSpace) :: ms
!  type(Ops) :: h
!  type(HFSolver) :: hf
!  character(:), allocatable :: file_nn, file_3n
!
!  call timer%init()
!  call ms%init('O16', 35.d0, 8, 16, e3max=6, is_three_body_jt=.true.)
!  call h%init('hamil',ms, 3)
!
!  file_nn = '/home/takayuki/MtxElmnt/2BME/TwBME-HO_NN-only_N3LO_EM500_srg2.00_hw35_emax8_e2max16.me2j.gz'
!  file_3n = '/home/takayuki/MtxElmnt/3BME/ThBME_srg2.00_N3LO_EM500_ChEFT_N2LO_cD-0.20cE0.098_Local2_IS_hw35_ms6_6_6.me3j'
!
!  call h%set(file_nn,file_3n,[8,16,8],[6,6,6,6])
!  call hf%init(h)
!  call hf%solve()
!  !call h%NormalOrdering()
!  call hf%fin()
!  call h%fin()
!  call ms%fin()
!
!  call timer%fin()
!
!end program test
