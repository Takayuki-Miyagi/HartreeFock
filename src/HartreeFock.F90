module HartreeFock
  use LinAlgLib
  use Operators
  implicit none

  private :: InitMonopole2
  private :: InitMonopole3_sp
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
    procedure :: InitMonopole3_sp
    procedure :: InitMonopole3
    procedure :: FinMonopole
  end type Monopole

  ! Note:
  ! alpha is parameter controlling
  ! density matrix mixing ratio
  ! (has to be 0.d0 < alpha <= 1.d0, alpha = 1.d0 means direct iteration)
  type :: HFSolver
    logical :: is_three_body
    integer :: n_iter_max = 1000
    real(8) :: alpha = 1.d0
    real(8) :: tol = 1.d-8
    real(8) :: diff
    real(8) :: e0  ! cm term
    real(8) :: e1  ! kinetic term
    real(8) :: e2  ! nn int. term
    real(8) :: e3  ! 3n int. term
    real(8) :: ehf ! hf energy

    type(NBodyPart) :: F   ! Fock opeartor
    type(NBodyPart) :: T   ! kinetic term
    type(NBodyPart) :: V   ! one-body filed from 2body interaction
    type(NBOdyPart) :: W   ! one-body filed form 3body interaction
    type(NBodyPart) :: rho ! density matrix
    type(NBodyPart) :: C   ! F's diagonalization coefficient, (HO|HF)
    type(NBodyPart) :: Occ ! diagonal Occupation matrix
    type(Monopole) :: V2, V3
  contains
    procedure :: fin => FinHFSolver
    procedure :: init => InitHFSolver
    procedure :: solve => SolveHFSolver
    procedure :: PrintSPEs
    procedure :: DiagonalizeFockMatrix
    procedure :: SetOccupationMatrix
    procedure :: UpdateDensityMatrix
    procedure :: UpdateFockMatrix
    procedure :: CalcEnergy
    procedure :: TransformToHF
    procedure :: HFBasisHamiltonian
    procedure :: HFBasisScalar
    !procedure :: HFBasisTensor
  end type HFSolver

contains

  subroutine FinHFSolver(this)
    class(HFSolver), intent(inout) :: this
    call this%rho%fin()
    call this%F%fin()
    call this%T%fin()
    call this%V%fin()
    call this%W%fin()
    call this%V2%FinMonopole()
    call this%V3%FinMonopole()
  end subroutine FinHFSolver

  subroutine InitHFSolver(this,ms,hamil,n_iter_max,tol,alpha)
    use Profiler, only: timer
    class(HFSolver), intent(inout) :: this
    type(MSpace), intent(in) :: ms
    type(Op), intent(in) :: hamil
    integer, intent(in), optional :: n_iter_max
    real(8), intent(in), optional :: tol, alpha
    real(8) :: ti
    integer :: ich

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
    this%is_three_body = hamil%is_three_body
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
    call this%V2%InitMonopole2(ms, hamil%two)
    call timer%countup_memory('Monople 2Body int.')
    call timer%Add("Construct Monopole 2Body int.", omp_get_wtime() - ti)

    if(this%is_three_body) then

      ti = omp_get_wtime()
      call timer%cmemory()
!#ifdef single_precision
      call this%V3%InitMonopole3_sp(ms, hamil%thr)
!#else
!      call this%V3%InitMonopole3(ms, hamil%thr)
!#endif
      call timer%countup_memory('Monople 3Body int.')
      call timer%Add("Construct Monopole 3Body int.", omp_get_wtime() - ti)
    end if

    do ich = 1, this%C%NChan
      call this%C%MatCh(ich,ich)%eye( this%C%MatCh(ich,ich)%ndims(1) )
    end do
    call this%SetOccupationMatrix(ms%one, ms%NOcoef)
    call this%UpdateDensityMatrix()
    call this%UpdateFockMatrix(ms%sps, ms%one)

  end subroutine InitHFSolver

  subroutine SolveHFSolver(this, sps, one)
    use Profiler, only: timer
    class(HFSolver), intent(inout) :: this
    type(Orbits), intent(in) :: sps
    type(OneBodySpace), intent(in) :: one
    integer :: iter
    real(8) :: ti

    ti = omp_get_wtime()
    write(*,'(2x,a,9x,a,10x,a,9x,a,7x,a,14x,a)') "iter", "zero-body", &
        & "one-body ", "two-body ", &
        & "three-body ", "Ehf"

    call this%CalcEnergy(one)
    write(*,'(4x,a,5f18.6)') "HO", this%e0, this%e1, &
        &  this%e2, this%e3,this%ehf
    do iter = 1, this%n_iter_max

      call this%DiagonalizeFockMatrix()
      call this%UpdateDensityMatrix()
      call this%UpdateFockMatrix(sps, one)

      call this%CalcEnergy(one)

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

  subroutine TransformToHF(HF,ms,Optr)
    !  Input: Operator is HO basis operator (not normal ordered)
    ! Output: Operator is HF basis operator (normal ordred, NO2B approximated)
    class(HFSolver), intent(in) :: HF
    type(MSpace), intent(in) :: ms
    type(Op), intent(inout) :: Optr

    if(Optr%is_normal_ordered) then
      write(*,'(a)') "In TransformToHF, input operator should not be normal ordered wrt target nucleus."
      return
    end if

    if(Optr%optr=='hamil' .or. Optr%optr=='Hamil') then
      call HF%HFBasisHamiltonian(ms,Optr)
      return
    end if

    if(Optr%Scalar) then
      call HF%HFBasisScalar(ms,Optr)
      return
    end if

    if(.not. Optr%Scalar) then
      write(*,*) "Not implemented yet."
      return
      !call HF%HFBasisTensor(ms,Optr)
    end if

  end subroutine TransformToHF

  subroutine HFBasisHamiltonian(HF,ms,H)
    use Profiler, only: timer
    class(HFSolver), intent(in) :: HF
    type(MSpace), intent(in) :: ms
    type(Op), intent(inout) :: H
    integer :: ch, J, n, bra, ket, a, b, c, d, e, f, JJJ
    integer :: ea, eb, ec, ed
    integer :: je, le, ze, ee
    integer :: jf, lf, zf, ef
    real(8) :: ph, ti
    type(DMat) :: UT, V2, V3

    ti = omp_get_wtime()

    H%zero = HF%ehf
    do ch = 1, ms%one%NChan
      H%one%MatCh(ch,ch)%DMat = HF%C%MatCh(ch,ch)%DMat%T() * &
          &  HF%F%MatCh(ch,ch)%DMat * HF%C%MatCh(ch,ch)%DMat
    end do

    do ch = 1, ms%two%NChan
      J = ms%two%jpz(ch)%j
      n = ms%two%jpz(ch)%nst
      call UT%zeros(n,n)
      call V3%zeros(n,n)
      V2 = H%two%MatCh(ch,ch)%DMat

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

          UT%m(bra,ket) = HF%C%GetOBME(ms%sps,ms%one,a,c) * HF%C%GetOBME(ms%sps,ms%one,b,d)
          if(a/=b) UT%m(bra,ket) = UT%m(bra,ket) - ph * &
              & HF%C%GetOBME(ms%sps,ms%one,a,d) * HF%C%GetOBME(ms%sps,ms%one,b,c)
          if(a==b) UT%m(bra,ket) = UT%m(bra,ket) * dsqrt(2.d0)
          if(c==d) UT%m(bra,ket) = UT%m(bra,ket) / dsqrt(2.d0)

          if(ket > bra) cycle
          if(.not. H%is_three_body) cycle
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
                V3%m(bra,ket) = V3%m(bra,ket) + HF%rho%GetOBME(ms%sps,ms%one,e,f) * &
                    & dble(JJJ+1) * H%thr%GetThBME(ms,a,b,e,J,c,d,f,J,JJJ)
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
      H%two%MatCh(ch,ch)%DMat = UT%T() * (V2+V3) * UT
      call UT%fin()
      call V2%fin()
      call V3%fin()
    end do
    call H%DiscardThreeBodyPart()

    call timer%Add("HFBasisHamltonian",omp_get_wtime() - ti)

  end subroutine HFBasisHamiltonian

  subroutine HFBasisScalar(HF,ms,Opr)
    use Profiler, only: timer
    class(HFSolver), intent(in) :: HF
    type(MSpace), intent(in) :: ms
    type(Op), intent(inout) :: Opr
    type(NBodyPart) :: o2from3, o1from3, o1from2
    integer :: ch, J, n, bra, ket, a, b, c, d, e, f, JJJ
    integer :: ea, eb, ec, ed
    integer :: je, le, ze, ee
    integer :: jf, lf, zf, ef
    real(8) :: ph, ti, o0from1, o0from2, o0from3
    type(DMat) :: UT, V2, V3

    ti = omp_get_wtime()
    Opr%zero = 0.d0
    do ch = 1, ms%one%NChan
      Opr%one%MatCh(ch,ch)%DMat = HF%C%MatCh(ch,ch)%DMat%T() * &
          &  Opr%one%MatCh(ch,ch)%DMat * HF%C%MatCh(ch,ch)%DMat
    end do

    call o2from3%init(ms%two, .true., Opr%optr, 0, 1, 0)
    do ch = 1, ms%two%NChan
      J = ms%two%jpz(ch)%j
      n = ms%two%jpz(ch)%nst
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

          UT%m(bra,ket) = HF%C%GetOBME(ms%sps,ms%one,a,c) * HF%C%GetOBME(ms%sps,ms%one,b,d)
          if(a/=b) UT%m(bra,ket) = UT%m(bra,ket) - ph * &
              & HF%C%GetOBME(ms%sps,ms%one,a,d) * HF%C%GetOBME(ms%sps,ms%one,b,c)
          if(a==b) UT%m(bra,ket) = UT%m(bra,ket) * dsqrt(2.d0)
          if(c==d) UT%m(bra,ket) = UT%m(bra,ket) / dsqrt(2.d0)

          if(ket > bra) cycle
          if(.not. Opr%is_three_body) cycle
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
                    & HF%rho%GetOBME(ms%sps,ms%one,e,f) * &
                    & dble(JJJ+1) * Opr%thr%GetThBME(ms,a,b,e,J,c,d,f,J,JJJ)
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
      Opr%two%MatCh(ch,ch)%DMat = &
          & UT%T() * Opr%two%MatCh(ch,ch)%DMat * UT
      o2from3%MatCh(ch,ch)%DMat = &
          & UT%T() * o2from3%MatCh(ch,ch)%DMat * UT
      call UT%fin()
      call V2%fin()
    end do
    call Opr%DiscardThreeBodyPart()

    o1from3 = o2from3%NormalOrderingFrom2To1(ms)
    o1from2 = Opr%two%NormalOrderingFrom2To1(ms)

    o0from3 = o1from3%NormalOrderingFrom1To0(ms)
    o0from2 = o1from2%NormalOrderingFrom1To0(ms)
    o0from1 = Opr%one%NormalOrderingFrom1To0(ms)

    Opr%zero = o0from1 + o0from2 * 0.5d0 + o0from3 / 6.d0
    Opr%one = Opr%one + o1from2 + o1from3 * 0.5d0
    Opr%two = Opr%two + o2from3

    call timer%Add("HFBasisScalar",omp_get_wtime() - ti)

  end subroutine HFBasisScalar

  !subroutine HFBasisTensor(HF,ms,Opr)
  !  use Profiler, only: timer
  !  class(HFSolver), intent(in) :: HF
  !  type(MSpace), intent(in) :: ms
  !  type(Op), intent(inout) :: Opr
  !  type(NBodyPart) :: o2from3, o1from3, o1from2
  !  integer :: chbra, chket, Jbra, Jket, nbra, nket
  !  integer :: bra, ket, a, b, c, d, e, f, JJJ
  !  integer :: ea, eb, ec, ed
  !  integer :: je, le, ze, ee
  !  integer :: jf, lf, zf, ef
  !  real(8) :: ph, ti, o0from1, o0from2, o0from3
  !  type(DMat) :: UT, V2, V3

  !  ti = omp_get_wtime()
  !  Opr%zero = 0.d0
  !  do ch = 1, ms%one%NChan
  !    Opr%one%MatCh(ch,ch)%DMat = HF%C%MatCh(ch,ch)%DMat%T() * &
  !        &  Opr%one%MatCh(ch,ch)%DMat * HF%C%MatCh(ch,ch)%DMat
  !  end do

  !  call o2from3%init(ms%two, .true., Opr%optr, 0, 1, 0)
  !  do chbra = 1, ms%two%NChan
  !    Jbra = ms%two%jpz(chbra)%j
  !    nbra = ms%two%jpz(chbra)%nst
  !    do chket = 1, ms%two%NChan
  !      Jket = ms%two%jpz(chket)%j
  !      nket = ms%two%jpz(chket)%nst
  !      call UT%zeros(nbra,nket)
  !      call V3%zeros(nbra,nket)

  !      !$omp parallel
  !      !$omp do private(bra,a,b,ea,eb,ph,ket,c,d,ec,ed,&
  !      !$omp &  e,je,le,ze,ee,f,jf,lf,zf,ef,JJJ)
  !      do bra = 1, nbra
  !        a = ms%two%jpz(ch)%n2spi1(bra)
  !        b = ms%two%jpz(ch)%n2spi2(bra)
  !        ea = ms%sps%orb(a)%e
  !        eb = ms%sps%orb(b)%e
  !        ph = (-1.d0)**((ms%sps%orb(a)%j+ms%sps%orb(b)%j)/2-J)
  !        do ket = 1, nket
  !          c = ms%two%jpz(ch)%n2spi1(ket)
  !          d = ms%two%jpz(ch)%n2spi2(ket)
  !          ec = ms%sps%orb(c)%e
  !          ed = ms%sps%orb(d)%e

  !          UT%m(bra,ket) = HF%C%GetOBME(ms%sps,ms%one,a,c) * HF%C%GetOBME(ms%sps,ms%one,b,d)
  !          if(a/=b) UT%m(bra,ket) = UT%m(bra,ket) - ph * &
  !              & HF%C%GetOBME(ms%sps,ms%one,a,d) * HF%C%GetOBME(ms%sps,ms%one,b,c)
  !          if(a==b) UT%m(bra,ket) = UT%m(bra,ket) * dsqrt(2.d0)
  !          if(c==d) UT%m(bra,ket) = UT%m(bra,ket) / dsqrt(2.d0)

  !          if(ket > bra) cycle
  !          if(.not. H%is_three_body) cycle
  !          do e = 1, ms%sps%norbs
  !            je = ms%sps%orb(e)%j
  !            le = ms%sps%orb(e)%l
  !            ze = ms%sps%orb(e)%z
  !            ee = ms%sps%orb(e)%e
  !            if(ea+eb+ee > ms%e3max) cycle
  !            do f = 1, ms%sps%norbs
  !              jf = ms%sps%orb(f)%j
  !              lf = ms%sps%orb(f)%l
  !              zf = ms%sps%orb(f)%z
  !              ef = ms%sps%orb(f)%e
  !              if(je /= jf) cycle
  !              if(le /= lf) cycle
  !              if(ze /= zf) cycle
  !              if(ec+ed+ef > ms%e3max) cycle

  !              do JJJ = abs(2*J-je), (2*J+je), 2
  !                o2from3%MatCh(ch,ch)%m(bra,ket) = &
  !                    & o2from3%MatCh(ch,ch)%m(bra,ket) + &
  !                    & HF%rho%GetOBME(ms%sps,ms%one,e,f) * &
  !                    & dble(JJJ+1) * Opr%thr%GetThBME(ms,a,b,e,J,c,d,f,J,JJJ)
  !              end do

  !            end do
  !          end do
  !          o2from3%MatCh(ch,ch)%m(bra,ket) = &
  !              & o2from3%MatCh(ch,ch)%m(bra,ket) / dble(2*J+1)
  !          if(a==b) o2from3%MatCh(ch,ch)%m(bra,ket) = &
  !              &    o2from3%MatCh(ch,ch)%m(bra,ket) / dsqrt(2.d0)
  !          if(c==d) o2from3%MatCh(ch,ch)%m(bra,ket) = &
  !              &    o2from3%MatCh(ch,ch)%m(bra,ket) / dsqrt(2.d0)
  !          o2from3%MatCh(ch,ch)%m(ket,bra) = &
  !              & o2from3%MatCh(ch,ch)%m(bra,ket)
  !        end do
  !      end do
  !      !$omp end do
  !      !$omp end parallel
  !      Opr%two%MatCh(ch,ch)%DMat = &
  !          & UT%T() * Opr%two%MatCh(ch,ch)%DMat * UT
  !      o2from3%MatCh(ch,ch)%m(ket,bra) = &
  !          & UT%T() * o2from3%MatCh(ch,ch)%m(bra,ket) * UT
  !      call UT%fin()
  !      call V2%fin()
  !    end do
  !  end do
  !  call H%DiscardThreeBodyPart()

  !  o1from3 = o2from3%NormalOrdering2To1(ms)
  !  o1from2 = Opr%two%NormalOrdering2To1(ms)

  !  o0from3 = o1from3%NormalOrdering1To0(ms)
  !  o0from2 = o1from2%NormalOrdering1To0(ms)
  !  o0from1 = Opr%one%NormalOrdering1To0(ms)

  !  Opr%zero = o0from1 + o0from2 * 0.5d0 + o0from3 / 6.d0
  !  Opr%one = Opr%one + o1from2 + o1from3 * 0.5d0
  !  Opr%two = Opr%two + o2from3

  !  call timer%Add("HFBasisTensor",omp_get_wtime() - ti)

  !end subroutine HFBasisTensor

  subroutine CalcEnergy(this,one)
    class(HFSolver), intent(inout) :: this
    type(OneBodySpace), intent(in) :: one
    integer :: ch, bra, ket, j

    this%e1 = 0.d0
    this%e2 = 0.d0
    this%e3 = 0.d0
    do ch = 1, this%F%NChan
      j = one%jpz(ch)%j
      do bra = 1, this%F%MatCh(ch,ch)%ndims(1)
        do ket = 1, this%F%MatCh(ch,ch)%ndims(1)
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
    integer :: ch

    do ch = 1, this%F%NChan
      call sol%init(this%F%MatCh(ch,ch)%DMat)
      call sol%DiagSym(this%F%MatCh(ch,ch)%DMat)
      this%C%MatCh(ch,ch)%DMat = sol%vec
      call sol%fin()
    end do
  end subroutine DiagonalizeFockMatrix

  subroutine SetOccupationMatrix(this,one,NOcoef)
    class(HFSolver), intent(inout) :: this
    type(OneBodySpace), intent(in) :: one
    real(8), intent(in) :: NOCoef(:)
    integer :: i, ich, io

    do ich = 1, one%NChan
      do i = 1, one%jpz(ich)%nst
        io = one%jpz(ich)%n2spi(i)
        this%Occ%MatCh(ich,ich)%m(i,i) = NOcoef(io)
      end do
    end do
  end subroutine SetOccupationMatrix

  subroutine UpdateDensityMatrix(this)
    class(HFSolver), intent(inout) :: this
    integer :: ich

    do ich = 1, this%rho%NChan
      this%rho%MatCh(ich,ich)%DMat = &
          & this%rho%MatCh(ich,ich)%DMat * (1.d0-this%alpha) + &
          & (this%C%MatCh(ich,ich)%DMat * this%Occ%MatCh(ich,ich)%DMat * &
          & this%C%MatCh(ich,ich)%DMat%T()) * this%alpha
    end do

  end subroutine UpdateDensityMatrix

  subroutine UpdateFockMatrix(this, sps, one)
    class(HFSolver), intent(inout) :: this
    type(Orbits), intent(in) :: sps
    type(OneBodySpace), intent(in) :: one
    type(NBodyPart) :: Fold
    integer(8) :: i1, i2, i3, i4, i5, i6, idx
    integer :: num, ch1, ch2, ch3
    integer :: ii1, ii2, ii3, ii4, ii5, ii6
    integer :: bra, ket

    Fold = this%F
    do ch1 = 1, this%F%NChan
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

    if(this%is_three_body) then
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

    do ch1 = 1, this%F%NChan
      do bra = 1, this%F%MatCh(ch1,ch1)%ndims(1)
        do ket = 1, bra
          this%V%MatCh(ch1,ch1)%m(ket,bra) = this%V%MatCh(ch1,ch1)%m(bra,ket)
          this%W%MatCh(ch1,ch1)%m(ket,bra) = this%W%MatCh(ch1,ch1)%m(bra,ket)
        end do
      end do
    end do

    this%diff = 0.d0
    this%F = this%T + this%V + this%W * 0.5d0
    do ch1 = 1, this%F%NChan
      this%diff = max(maxval(this%F%MatCh(ch1,ch1)%m - Fold%MatCh(ch1,ch1)%m), this%diff)
    end do
    call Fold%fin()

  end subroutine UpdateFockMatrix

  subroutine PrintSPEs(this,ms)
    class(HFSolver), intent(in) :: this
    type(MSpace), intent(in) :: ms
    integer :: i, io, ch
    type(NBodyPart) :: F_HF

    F_HF = this%F
    do ch = 1, ms%one%NChan
      F_HF%MatCh(ch,ch)%DMat = this%C%MatCh(ch,ch)%DMat%T() * &
          &  this%F%MatCh(ch,ch)%DMat * this%C%MatCh(ch,ch)%DMat
    end do

    write(*,'(a)') "  Hartree-Fock single-particle energies"

    do i = 1, size(ms%holes)
      io = ms%holes(i)
      write(*,'(a,a10,i4,f12.6)') 'hole:     ', trim(ms%sps%GetLabelFromIndex(io)), io, &
          & F_HF%GetOBME(ms%sps,ms%one,io,io)
    end do

    do i = 1, size(ms%particles)
      io = ms%particles(i)
      write(*,'(a,a10,i4,f12.6)') 'particle: ', trim(ms%sps%GetLabelFromIndex(io)), io, &
          & F_HF%GetOBME(ms%sps,ms%one,io,io)
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

  subroutine InitMonopole2(this, ms, vnn)
    class(Monopole), intent(inout) :: this
    type(MSpace), intent(in) :: ms
    type(NBodyPart), intent(in) :: vnn
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
            & vnn%GetTwBME(ms%sps,ms%two,&
            & int(i1,kind(JJ)), int(i2,kind(JJ)), int(i3,kind(JJ)), int(i4,kind(JJ)),JJ)
        ! need to convert integer(8) -> integer(4)
      end do
      this%v(idx) = v * norm / dble(j1+1)
    end do
    !$omp end do
    !$omp end parallel
    this%constructed = .true.
  end subroutine InitMonopole2

  subroutine InitMonopole3(this, ms, v3n)
    use CommonLibrary, only: triag
    class(Monopole), intent(inout) :: this
    type(MSpace), intent(in) :: ms
    type(NBodyPart), intent(in) :: v3n
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
              & v3n%GetThBME(ms,&
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
  end subroutine InitMonopole3

  subroutine InitMonopole3_sp(this, ms, v3n)
    use CommonLibrary, only: triag
    class(Monopole), intent(inout) :: this
    type(MSpace), intent(in) :: ms
    type(NBodyPartSp), intent(in) :: v3n
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
              & v3n%GetThBME(ms,&
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
  end subroutine InitMonopole3_sp

  function GetIndex2(i1,i2,i3,i4) result(r)
    integer(8), intent(in) :: i1, i2, i3, i4
    integer(8) :: r
    r = i1 + lshift(i2,10) + lshift(i3,20) + lshift(i4,30)
  end function GetIndex2

  subroutine GetSpLabels2(idx,i1,i2,i3,i4)
    integer(8), intent(in)  :: idx
    integer(8), intent(out) :: i1, i2, i3, i4
    i1 = mod(idx,1024)
    i2 = mod(rshift(idx,10),1024)
    i3 = mod(rshift(idx,20),1024)
    i4 = mod(rshift(idx,30),1024)
  end subroutine GetSpLabels2

  function GetIndex3(i1,i2,i3,i4,i5,i6) result(r)
    integer(8), intent(in) :: i1, i2, i3, i4, i5, i6
    integer(8) :: r
    r = i1 + lshift(i2,10) + lshift(i3,20) + lshift(i4,30) + &
        &    lshift(i5,40) + lshift(i6,50)
  end function GetIndex3

  subroutine GetSpLabels3(idx,i1,i2,i3,i4,i5,i6)
    integer(8), intent(in)  :: idx
    integer(8), intent(out) :: i1, i2, i3, i4, i5, i6
    i1 = mod(idx,1024)
    i2 = mod(rshift(idx,10),1024)
    i3 = mod(rshift(idx,20),1024)
    i4 = mod(rshift(idx,30),1024)
    i5 = mod(rshift(idx,40),1024)
    i6 = mod(rshift(idx,50),1024)
  end subroutine GetSpLabels3
end module HartreeFock

!program test
!  use Profiler, only: timer
!  use ClassSys, only: sys
!  use CommonLibrary, only: &
!      &init_dbinomial_triangle, fin_dbinomial_triangle
!  use ModelSpace, only: MSpace
!  use Operators
!  use HartreeFock, only: HFSolver
!  use MBPT
!  implicit none
!
!  type(MSpace) :: ms
!  type(Op) :: h
!  character(:), allocatable :: file_nn, file_3n
!  type(HFsolver) :: HF
!  type(sys) :: s
!  type(MBPTEnergy) :: PT
!
!  call timer%init()
!  call init_dbinomial_triangle()
!
!  call ms%init('O16', 25.d0, 6, 12, e3max=6)
!  !call h%init('hamil',ms,.false.) ! nn-only
!  call h%init('hamil',ms,.true.)  ! nn+3n
!
!  call h%set(ms,file_nn,file_3n,[8,12,8],[8,8,8,8])
!
!  call HF%init(ms,h,alpha=1.0d0)
!  call HF%solve(ms%sps,ms%one)
!
!  call HF%TransformToHF(ms,H)
!
!  call PT%calc(ms,H)
!
!  call HF%fin()
!
!  call h%fin()
!  call ms%fin()
!
!  call fin_dbinomial_triangle()
!  call timer%fin()
!
!end program test
