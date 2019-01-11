module class_stopwatch
#ifdef MPI
  use mpi
#endif
  !$ use omp_lib
  implicit none
  type stopwatch
    real(8) :: time = 0.d0, tstart = 0.d0
    integer :: ncall = 0
    logical :: is_on = .false.
  end type stopwatch
  type(stopwatch) :: time_total, time_set_hamil, time_HF
  type(stopwatch) :: time_calcop, time_set_ope, time_MS
  type(stopwatch) :: time_io_read, time_io_write

contains
  subroutine start_stopwatch(this, is_reset, is_mpi_barrier)
#ifdef MPI
    use MPIFunctions, only: ierr
#endif
    type(stopwatch), intent(inout) :: this
    logical, intent(in), optional :: is_reset, is_mpi_barrier
    real(8) :: wclock
    if (this%is_on) stop "ERROR: start_stopwatch"
    this%is_on = .true.
    !    !$ if (omp_get_thread_num()==0) then
#ifdef MPI
    if (present(is_mpi_barrier)) then
      if (is_mpi_barrier) call mpi_barrier(mpi_comm_world, ierr)
    end if
#endif
    this%tstart = wclock()
    if (present(is_reset)) then
      if (is_reset) then
        this%time = 0.d0
        this%ncall = 0
      end if
    end if
    !    !$ end if
  end subroutine start_stopwatch

  subroutine stop_stopwatch(this, time_last, is_mpi_barrier)
#ifdef MPI
    use MPIFunctions, only: ierr
#endif
    type(stopwatch), intent(inout) :: this
    real(8), intent(out), optional :: time_last
    logical, intent(in), optional :: is_mpi_barrier
    real(8) :: wclock, t
    if (.not. this%is_on) stop "ERROR: stop_stopwatch"
    this%is_on = .false.
    !    !$ if (omp_get_thread_num()==0) then
#ifdef MPI
    if (present(is_mpi_barrier)) then
      if (is_mpi_barrier) call mpi_barrier(mpi_comm_world, ierr)
    end if
#endif
    t = wclock() - this%tstart
    if (present(time_last)) time_last = t
    this%time = this%time + t
    this%ncall = this%ncall + 1
    !    !$ end if
  end subroutine stop_stopwatch

  subroutine reset_stopwatch(this)
    type(stopwatch), intent(inout) :: this
    this%is_on = .false.
    this%time = 0.d0
    this%tstart = 0.d0
    this%ncall = 0
  end subroutine reset_stopwatch

  function time_stopwatch(this) result (r)
    type(stopwatch), intent(in) :: this
    real(8) :: r
    r = this%time
  end function time_stopwatch

  function get_ctime_stopwatch(this, is_mpi_sync) result (r)
#ifdef MPI
    use MPIFunctions, only: ierr
#endif
    ! get current time of running stopwatch
    type(stopwatch), intent(in) :: this
    logical, intent(in), optional :: is_mpi_sync
    real(8) :: r
    real(8) :: wclock
    r = wclock() - this%tstart
#ifdef MPI
    if (.not. present(is_mpi_sync)) return
    if (.not. is_mpi_sync) return
    call mpi_bcast(r, 1, mpi_real8, 0, mpi_comm_world, ierr)
#endif
  end function get_ctime_stopwatch

  subroutine print_summary_stopwatch()
    use MPIFunction
    real(8) :: r
    integer :: n
    if (myrank==0) write(*,*)
    n = int(time_total%time)
    if (myrank==0) write(*,'(a,i5,a,i2.2,a,i2.2)') &
        '    summary of time, total = ', n/3600, ':', &
        mod(n, 3600)/60, ':', mod(n, 60)
    if (myrank==0) write(*,*)
    if (myrank==0) write(*,'(32x, a)') "time,    ncall, time/ncall,   ratio "
    call print_each(time_total, time_total, "total")
    call print_each(time_total, time_MS, "Set Model Space")
    call print_each(time_total, time_set_hamil, "Set Hamiltonian")
    call print_each(time_total, time_HF, "HF calculation")
    r = time_total%time &
        - time_MS%time &
        - time_set_hamil%time &
        - time_calcop%time &
        - time_HF%time
    if (myrank==0) write(*,'(1a25, 1f12.3, 22x, 1f9.4)') &
        "misc", r, r/time_total%time
    if(myrank == 0) write(*,*)
    if(myrank == 0) write(*,*) "    --------------------------------------------------"
    call print_each(time_HF, time_HF, "HF calculation")
    if(myrank == 0) write(*,*)
    if(myrank == 0) write(*,*) "    --------------------------------------------------"

    if(myrank == 0) write(*,*)
    if(myrank == 0) write(*,*) "    --------------------------------------------------"
    call print_each_wo_r(time_io_read,      "I/O LV read ")
    call print_each_wo_r(time_io_write,     "I/O LV write")
    if (myrank==0) write(*,*)
  contains

    subroutine print_each(total, this, title)
      type(stopwatch), intent(in) :: this, total
      character(len=*), intent(in) :: title
      if (this%time == 0.d0) return
      if (myrank==0) &
          write(*,'(1a25, 1f12.3, 1i10, 1f12.5, 1f9.4)') title,  &
          this%time, this%ncall, this%time/this%ncall, &
          this%time/total%time
    end subroutine print_each

    subroutine print_each_wo_r(this, title)
      type(stopwatch), intent(in) :: this
      character(len=*), intent(in) :: title
      if (this%time == 0.d0) return
      if (myrank==0) &
          write(*,'(1a25, 1f12.3, 1i10, 1f12.5)') title,  &
          this%time, this%ncall, this%time/this%ncall
    end subroutine print_each_wo_r

  end subroutine print_summary_stopwatch
end module class_stopwatch
