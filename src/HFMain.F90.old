program HartreeFockMain
  use common_library, only: init_dbinomial_triangle, fin_dbinomial_triangle, &
      & set_physics_constant
  use InputParameters, only: parameters, set_parameters
  use ModelSpace, only: Mspace, spo_pn, spo_isospin
  use read_3BME, only: iThreeBodyScalar
  use ScalarOperator
  use NormalOrdering
  use HFSolver, only: HFSol
  use WriteOperator
  use MBPT3, only: MBPT
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
  type(ScalarOperators) :: hamil, scalar
  type(HFSol) :: hf_sol
  type(MBPT) :: mbpt_sol
  integer :: iunite = 118
  call start_stopwatch(time_total)
#ifdef MPI
  call mpi_init(ierr)
  call mpi_comm_size(mpi_comm_world, nprocs, ierr)
  call mpi_comm_rank(mpi_comm_world, myrank, ierr)
#endif

  call set_physics_constant()
  call set_parameters(params)

  call init_dbinomial_triangle()
  call sps%init(params)

  ! Model Space
  call start_stopwatch(time_MS)
  call ms%init(sps, params)
  call stop_stopwatch(time_MS)

  call start_stopwatch(time_set_hamil)
  ! Three-Body Force
  if(params%thbmefile /= 'None') call isps%init(params%emax_3nf)
  if(params%thbmefile /= 'None') call thbme%init(isps, params, params%thbmefile)
  ! Hamiltonian
  call hamil%init(ms)
  call hamil%set(params, sps, ms, 'hamil', params%twbmefile, thbme)
  call stop_stopwatch(time_set_hamil)


  ! Hartree-Fock calculation
  call start_stopwatch(time_HF)
  call hf_sol%init(ms)
  if(params%thbmefile /= 'None') call hf_sol%solve(params, sps, ms, hamil, thbme)
  if(params%thbmefile == 'None') call hf_sol%solve(params, sps, ms, hamil)
  call stop_stopwatch(time_HF)

  if(params%thbmefile /= 'None') call thbme%fin(isps)
  if(params%thbmefile /= 'None') call isps%fin()
  open(iunite, file=params%egs, status = 'replace')
  call params%PrtParams(iunite)
  if(myrank == 0) write(iunite,'(4f15.6, a)') hf_sol%e1ho, hf_sol%e2ho, hf_sol%e3ho, hf_sol%eho, ' HO'
  if(myrank == 0) write(iunite,'(4f15.6, a)') hf_sol%e1hf, hf_sol%e2hf, hf_sol%e3hf, hf_sol%ehf, ' HF'


  if(params%sv_hf_rslt) then
    call write_hamil(params, sps, ms, hamil)
  end if

  if(params%vac == 'ref' .and. params%MBPT) then
    call mbpt_sol%calc(params, sps, ms, hamil)
    if(myrank == 0) write(iunite,'(4f15.6, a)') mbpt_sol%e_0, mbpt_sol%e_2, mbpt_sol%e_3, &
      & mbpt_sol%e_0 + mbpt_sol%e_2 + mbpt_sol%e_3, ' MBPT'
  end if

  call hamil%fin()

  ! --- other scalar operators

  ! --- bare operator
  call scalar%init(ms)
  call scalar%set(params, sps, ms, 'rm')
  call hf_sol%HFBasis(params, sps, ms, scalar)
  if(myrank == 0) write(iunite,'(4f15.6, a)') hf_sol%e1hf, hf_sol%e2hf, hf_sol%e3hf, hf_sol%ehf, ' HF'
  call scalar%fin()


  ! --- two-body effective operator
  call scalar%init(ms)
  call scalar%set(params, sps, ms, 'rm', f2 = params%scfile2)
  call hf_sol%HFBasis(params, sps, ms, scalar)
  if(myrank == 0) write(iunite,'(4f15.6, a)') hf_sol%e1hf, hf_sol%e2hf, hf_sol%e3hf, hf_sol%ehf, ' HF'
  call scalar%fin()

  ! --- three-body operator
  if(params%scfile3 /= 'None') call isps%init(params%emax_3nf)
  if(params%scfile3 /= 'None') call thbme%init(isps, params, params%scfile3)
  call scalar%init(ms)
  call scalar%set(params, sps, ms, 'rm', f2 = params%scfile2, thbme = thbme)
  if(params%scfile3 /= 'None') call hf_sol%HFBasis(params, sps, ms, scalar, thbme)
  if(params%scfile3 == 'None') call hf_sol%HFBasis(params, sps, ms, scalar)
  if(myrank == 0) write(iunite,'(4f15.6, a)') hf_sol%e1hf, hf_sol%e2hf, hf_sol%e3hf, hf_sol%ehf, ' HF'
  if(params%scfile3 /= 'None') call thbme%fin(isps)
  if(params%scfile3 /= 'None') call isps%fin()
  call scalar%fin()

  ! --- other scalar operators
  call hf_sol%fin()
  call ms%fin()
  call sps%fin()

  call fin_dbinomial_triangle()
#ifdef MPI
  call mpi_finalize(ierr)
#endif
  call stop_stopwatch(time_total)
  call print_summary_stopwatch()
end program HartreeFockMain
