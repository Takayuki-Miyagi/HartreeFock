module HFSolver
  use MPIFunction, only: myrank
  use InputParameters, only: parameters
  use ModelSpace, only: spo_pn, MSpace, OneBodySpace, TwoBodySpace, ThreeBodySpace, &
      & OneBodyChannel, TwoBodyChannel, ThreeBodyChannel
  use read_3bme, only: iThreeBodyScalar
  use ScalarOperator
  use NormalOrdering, only: NOThree2Two_HFbasis, NormOrd
  use MatrixDouble, only: DMat
  use VectorDouble, only: DVec
  use LinAlgLib
  use Optimizer, only: FixedPointProblem
  implicit none
  private :: InitHFSol, FinHFSol, HFSolLoop, InverseTransformation, HFOneBodyEquation, &
      & getVectorDim, getVectorOpt, getUpdate
  public :: HFSol, TwoBodyScalarHFBasis
  type :: HFSol
    type(NBodyScalars) :: HFBasis ! Transformation
    real(8) :: e1hf, e2hf, e3hf, ehf
    real(8) :: e1ho, e2ho, e3ho, eho
  contains
    procedure :: init => InitHFSol
    procedure :: fin => FinHFSol
    procedure :: solve => HFSolLoop
  end type HFSol
contains
  subroutine InitHFSol(this, ms)
    class(HFSol), intent(inout) :: this
    type(MSpace), intent(in) :: ms
    integer :: ich, n
    call this%HFBasis%init(ms%one%SpinParityTz)
    do ich = 1, ms%one%n
      n = ms%one%jptz(ich)%n
      call this%HFBasis%jptz(ich)%eye(n)
    end do
  end subroutine InitHFSol

  subroutine FinHFSol(this)
    class(HFSol), intent(inout) :: this
    call this%HFBasis%fin()
  end subroutine FinHFSol

  subroutine HFSolLoop(this, params, sps, ms, hamil, thbme)
    implicit none
    class(HFSol), intent(inout) :: this
    type(parameters), intent(in) :: params
    type(spo_pn), intent(in) :: sps
    type(MSPace), intent(in) :: ms
    type(ScalarOperators), intent(inout) :: hamil
    type(iThreeBodyScalar), optional, intent(in) :: thbme
    type(NBodyScalars) :: h1i, h1f, w1i, w1f, w12, t1f ! One-Body operators
    type(NBodyScalars) :: v2f, w2f, w2i ! Two-Body operators
    type(NOThree2Two_HFbasis) :: NO
    integer :: nite
    type(FixedPointProblem) :: opt
    write(*,*) '##################################################'
    write(*,*) '########   Hartree-Fock Iterations  ##############'
    write(*,*) '##################################################'

    call InitHFSolLoop()

    do nite = 1, 999
      opt%nite = nite
      h1i = hamil%one + w1i
      call HFOneBodyEquation(h1i, hamil%one, h1f, t1f, this%HFBasis)
      v2f = TwoBodyScalarHFBasis(ms%one, ms%two, hamil%two, this%HFBasis)
      w1f = NormOrd(params, sps, ms%one, ms%two, v2f)
      if(present(thbme)) then
        w2i = NormOrd(params, sps, this%HFBasis, ms%two, thbme, NO)
        !stop
        w2f = TwoBodyScalarHFBasis(ms%one, ms%two, w2i, this%HFBasis)
        w12 = NormOrd(params, sps, ms%one, ms%two, w2f)
      end if
      this%e1hf = NormOrd(sps, ms%one, t1f)
      this%e2hf = NormOrd(sps, ms%one, w1f) * 0.5d0
      this%e3hf = NormOrd(sps, ms%one, w12) / 6.d0
      this%ehf = this%e1hf + this%e2hf + this%e3hf
      if(myrank == 0) then
        write(*,'(a, i4, a, f15.6, a, f15.6, a, f15.6, a, f15.6, a, es10.3)') &
            & 'nite = ', nite, '  e1 = ', this%e1hf, '  e2 = ', this%e2hf, &
            & '  e3 = ', this%e3hf, '  ehf = ', this%ehf, '  error = ', opt%si
      end if
      w1f = w1f + 0.5d0 * w12
      w1i = InverseTransformation(w1f, this%HFbasis)
      call getVectorOpt(opt, w1i)
      call opt%GetVector()
      if(opt%si < params%conv) exit
      call getUpdate(opt, w1i)
    end do

    call FinHFSolLoop()
  contains
    subroutine InitHFSolLoop()
      integer :: n
      call h1i%init(ms%one%SpinParityTz)
      call h1f%init(ms%one%SpinParityTz)
      call w1i%init(ms%one%SpinParityTz)
      call w1f%init(ms%one%SpinParityTz)
      call w12%init(ms%one%SpinParityTz)
      call t1f%init(ms%one%SpinParityTz)

      call v2f%init(ms%two%SpinParityTz)
      call w2f%init(ms%two%SpinParityTz)
      call w2i%init(ms%two%SpinParityTz)

      n = GetVectorDim(w1i)
      !call opt%init(n, alpha = 0.7d0, method = 'direct', m = 10)
      !call opt%init(n, alpha = 0.7d0, method = 'broyden', m = 10)
      call opt%init(n, alpha = 0.7d0, method = 'mbroyden', m = 10)

      if(present(thbme)) call NO%init(params, sps, ms%one, ms%two)
      w1i = NormOrd(params, sps, ms%one, ms%two, hamil%two)
      if(present(thbme)) then
        w2i = NormOrd(params, sps, this%HFBasis, ms%two, thbme, NO)
        w12 = NormOrd(params, sps, ms%one, ms%two, w2i)
      end if
      this%e1ho = NormOrd(sps, ms%one, hamil%one)
      this%e2ho = NormOrd(sps, ms%one, w1i) * 0.5d0
      this%e3ho = NormOrd(sps, ms%one, w12) / 6.d0
      this%eho = this%e1ho + this%e2ho + this%e3ho
      w1i = w1i + 0.5d0 * w12
      if(myrank == 0) then
        write(*,'(a, f9.4, a)') &
            & 'Memory around : ', &
            &  3.d0 * hamil%two%usedmem + NO%usedmem, &
            & ' GB will be additionally used during Hartree-Fock iteractions'
      end if
      if(myrank == 0) then
        write(*,'(a)') 'initial normal ordering'
        write(*,'(a, f15.6, a, f15.6, a, f15.6, a, f15.6)') &
            & '  e1 = ', this%e1ho, '  e2 = ', this%e2ho, &
            & '  e3 = ', this%e3ho, '  ehf = ', this%eho
      end if
    end subroutine InitHFSolLoop

    subroutine FinHFSolLoop()
      if(present(thbme)) call NO%fin()

      if(params%vac == 'vacumm') then
        !-- NO2B w.r.t. vacuum
        hamil%zero = this%e3hf
        hamil%one = t1f - 0.5d0 * w12
        hamil%two = v2f + w2f
      elseif(params%vac == 'ref') then
        !-- NO2B w.r.t. refernce state
        hamil%zero = this%ehf
        hamil%one = h1f
        hamil%two = v2f + w2f
      end if

      call h1i%fin()
      call h1f%fin()
      call w1i%fin()
      call w1f%fin()
      call w12%fin()
      call t1f%fin()

      call v2f%fin()
      call w2f%fin()
      call w2i%fin()
    end subroutine FinHFSolLoop
  end subroutine HFSolLoop

  function TwoBodyScalarHFBasis(one, two, a, Trs) result(b)
    type(NBodyScalars) :: b
    type(OneBodySpace), intent(in) :: one
    type(TwoBodySpace), intent(in) :: two
    type(NbodyScalars), intent(in) :: a, Trs
    type(DMat) :: U
    integer :: ich, n
    b = a
    do ich = 1, two%n
      n = two%jptz(ich)%n
      call U%ini(n,n)
      call UEmbedded(two%jptz(ich), U, one, Trs)
      b%jptz(ich) = U%T() * a%jptz(ich) * U
      call U%fin()
    end do
  end function TwoBodyScalarHFBasis

  function InverseTransformation(a, Trs) result(b)
    type(NBodyScalars) :: b
    type(NBodyScalars), intent(in) :: a, Trs
    integer :: ich
    b = a
    do ich = 1, size(a%jptz)
      b%jptz(ich) = Trs%jptz(ich) * a%jptz(ich) * Trs%jptz(ich)%T()
    end do
  end function InverseTransformation

  subroutine HFOneBodyEquation(h1i, t1i, h1f, t1f, Trs)
    type(NBodyScalars), intent(in) :: h1i, t1i
    type(NBodyScalars), intent(inout) :: h1f, t1f, Trs
    integer :: ich
    type(EigenSolSymD) :: sol

    do ich = 1, size(t1i%jptz)
      call sol%init(h1i%jptz(ich))
      call sol%DiagSym(h1i%jptz(ich))
      Trs%jptz(ich) = sol%vec
      call h1f%jptz(ich)%DiagMat(sol%eig)
      t1f%jptz(ich) = sol%vec%T() * t1i%jptz(ich) * sol%vec
      call sol%fin()
    end do
  end subroutine HFOneBodyEquation

  function getVectorDim(sc) result(n)
    integer :: n
    type(NBodyScalars), intent(in) :: sc
    integer :: ich, nn
    n = 0
    do ich = 1, size(sc%jptz)
      nn = size(sc%jptz(ich)%m, 1)
      n = n + (nn + 1) * nn / 2
    end do
  end function getVectorDim

  subroutine getVectorOpt(opt, w1)
    type(FixedPointProblem), intent(inout) :: opt
    type(NBodyScalars), intent(in) :: w1
    integer :: ich, i, j, num, n
    num = 0
    do ich = 1, size(w1%jptz)
      n = size(w1%jptz(ich)%m, 1)
      do i = 1, n
        do j = 1, i
          num = num + 1
          opt%fj%v(num) = w1%jptz(ich)%m(i, j)
        end do
      end do
    end do
  end subroutine getVectorOpt

  subroutine getUpdate(opt, w1)
    type(FixedPointProblem), intent(in) :: opt
    type(NBodyScalars), intent(inout) :: w1
    integer :: ich, i, j, num, n
    num = 0
    do ich = 1, size(w1%jptz)
      n = size(w1%jptz(ich)%m, 1)
      do i = 1, n
        do j = 1, i
          num = num + 1
          w1%jptz(ich)%m(i,j) = opt%xj%v(num)
          w1%jptz(ich)%m(j,i) = opt%xj%v(num)
        end do
      end do
    end do
  end subroutine getUpdate
end module HFSolver
