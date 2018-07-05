program HartreeFockMain
  use InputParameters, only: parameters, set_parameters
  use RotationGroup, only: init_dbinomial_triangle, fin_dbinomial_triangle
  use ModelSpace, only: Mspace, spo_pn, spo_isospin
  use read_3BME, only: iThreeBodyScalar
  use ScalarOperator
  use NormalOrdering
  use HFSolver, only: HFSol
  use WriteOperator
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
  type(iThreeBodyScalar) :: thbme
  type(ScalarOperators) :: hamil
  type(HFSol) :: sol
  integer :: iunite = 118
  call start_stopwatch(time_total)
#ifdef MPI
  call mpi_init(ierr)
  call mpi_comm_size(mpi_comm_world, nprocs, ierr)
  call mpi_comm_rank(mpi_comm_world, myrank, ierr)
#endif

  call set_parameters(params)

  call init_dbinomial_triangle()
  call sps%init(params)

  ! Model Space
  call start_stopwatch(time_MS)
  call ms%init(sps, params)
  call stop_stopwatch(time_MS)

  call start_stopwatch(time_set_hamil)
  ! Three-Body Force
  if(params%thbmefile /= 'None') call isps%init(params)
  if(params%thbmefile /= 'None') call thbme%init(isps, params, params%thbmefile)
  ! Hamiltonian
  call hamil%init(ms)
  call hamil%set(params, sps, ms, 'hamil', params%twbmefile, thbme)
  call stop_stopwatch(time_set_hamil)


  ! Hartree-Fock calculation
  call start_stopwatch(time_HF)
  call sol%init(ms)
  if(params%thbmefile /= 'None') call sol%solve(params, sps, ms, hamil, thbme)
  if(params%thbmefile == 'None') call sol%solve(params, sps, ms, hamil)
  call stop_stopwatch(time_HF)

  if(params%thbmefile /= 'None') call thbme%fin(isps)
  if(params%thbmefile /= 'None') call isps%fin()
  open(iunite, file=params%egs, status = 'replace')
  call params%PrtParams(iunite)
  call params%PrtParams(iunite)
  if(myrank == 0) write(iunite,'(4f15.6, a)') sol%e1ho, sol%e2ho, sol%e3ho, sol%eho, ' HO'
  if(myrank == 0) write(iunite,'(4f15.6, a)') sol%e1hf, sol%e2hf, sol%e3hf, sol%ehf, ' HF'

  if(params%sv_hf_rslt) then
    call write_hamil(params, sps, ms, hamil)
  end if
  call hamil%fin()

  ! --- other scalar operators




  ! --- other scalar operators
  call sol%fin()
  call ms%fin()
  call sps%fin()

  call fin_dbinomial_triangle()
#ifdef MPI
  call mpi_finalize(ierr)
#endif
  call stop_stopwatch(time_total)
  call print_summary_stopwatch()
end program HartreeFockMain
