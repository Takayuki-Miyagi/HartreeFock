program HartreeFockMain
  use InputParameters, only: parameters, set_parameters
  use ModelSpace, only: Mspace, spo_pn, spo_isospin
  use RotationGroup, only: init_dbinomial_triangle, fin_dbinomial_triangle
!  use read_3BME, only: ithbdy
!  use Operators, only: Scalar, Tensor, OneBodyScalar
!  use HFCalc, only: HartreeFock, NO2B
#ifdef MPI
  use mpi
  use MPIFunction, only: myrank, nprocs, ierr
#else
  use MPIFunction, only: myrank
#endif
  use class_stopwatch, only: time_total, start_stopwatch, &
      & stop_stopwatch, print_summary_stopwatch, time_calcop, &
      & time_HF, time_MS, time_set_hamil, time_set_ope
  !$ use omp_lib
  implicit none
  type(parameters) :: params
  type(spo_pn) :: sps
  type(spo_isospin) :: isps
  type(MSpace) :: ms
!  type(ithbdy) :: thbdyfrc
!  type(Scalar) :: hamil, s, x, sclop
!  type(Tensor) :: tnsop
!  type(OneBodyScalar) :: HFT
!  type(HartreeFock) :: HF
  integer :: iunite = 118
  integer :: iunito = 119
  integer :: iunits = 120
  integer :: iunith = 121
!  character(10) :: RankOperator
!  integer :: jr, pr, zr
  call start_stopwatch(time_total)
#ifdef MPI
  call mpi_init(ierr)
  call mpi_comm_size(mpi_comm_world, nprocs, ierr)
  call mpi_comm_rank(mpi_comm_world, myrank, ierr)
#endif

  call set_parameters(params)
  open(iunite, file=params%egs, status = 'replace')
  call params%PrtParams(iunite)

  call init_dbinomial_triangle()
  call sps%init(params)

!  ! Model Space
  call start_stopwatch(time_MS)
  call ms%init(sps, params)
  call stop_stopwatch(time_MS)

  call start_stopwatch(time_set_hamil)
  ! Three-Body Force
  if(params%thbmefile /= 'None') call isps%init(params)
  if(params%thbmefile /= 'None') call thbdyfrc%SetThrBdyScl(isps, params, params%thbmefile)
  ! Hamiltonian
!  call hamil%InitScalar(ms, params)
!  call hamil%SetHamil(sps, ms, params, thbdyfrc, params%twbmefile)
!  call stop_stopwatch(time_set_hamil)
!
!  ! Hartree-Fock calculation
!  call start_stopwatch(time_HF)
!  call HFT%InitOneBodyScalar(ms%one)
!  call HF%HFLoop(params, sps, ms, hamil, HFT, thbdyfrc, iunite)
!  if(params%sv_hf_rslt) then
!    call HF%PrintHFResult(params, sps, ms, hamil)
!  end if
!  call stop_stopwatch(time_HF)
!  if(myrank == 0) write(iunite,'(4f15.6, a)') HF%e1, HF%e2, HF%e3, HF%etot, ' HF'
!  if(params%thbme) call thbdyfrc%ReleaseJPT(isps)
!  if(params%thbme) call isps%ReleaseSPSIsospin()
!
!  call hamil%ReleaseScalar(ms, params)
!  call HFT%ReleaseOneBodyScalar(ms%one)
!  call ms%ReleaseModelSpace(sps, params)
!  call sps%ReleaseSPS()
!  call fin_dbinomial_triangle()
!
!#ifdef MPI
!  call mpi_finalize(ierr)
!#endif
!  call stop_stopwatch(time_total)
!  call print_summary_stopwatch()
end program HartreeFockMain
