module MPIFunction
#ifdef MPI
  use mpi
#endif
  implicit none
  integer :: myrank = 0, nprocs = 1, ierr = 0
#ifdef MPI
  integer :: idest, idummy
  integer :: istatus(mpi_status_size)
  integer :: procs
#endif
contains
  subroutine master_slave(Method, ntotal, slranks)
    interface
      subroutine Method(cnt)
        integer, intent(in) :: cnt
      end subroutine Method
    end interface
    integer, intent(in) :: ntotal
    integer, allocatable, intent(inout) :: slranks(:)
    integer :: num_loops, i
    if(.not.allocated(slranks)) allocate(slranks(ntotal))
    slranks = 0
#ifdef MPI
    if(myrank == 0)then
      write(*,*)
      write(*,'(a)')"Master-slave procedure"
      write(*,'(a, i4)')"Total # of sets = ", ntotal
    end if
    call mpi_barrier(mpi_comm_world, ierr)
#endif
    if(myrank == 0) then
      do num_loops = 1, ntotal
#ifdef MPI
        call mpi_recv(idummy,1,mpi_integer,mpi_any_source,1,mpi_comm_world,istatus,ierr)
        idest = istatus(mpi_source)
        slranks(num_loops) = idest
        call mpi_send(num_loops,1,mpi_integer,idest,1,mpi_comm_world,ierr)
      end do
      do i = 1, nprocs-1 ! finalize process
        call mpi_recv(idummy,1,mpi_integer,mpi_any_source,1,mpi_comm_world,istatus,ierr)
        idest = istatus(mpi_source)
        num_loops = 0
        call mpi_send(num_loops,1,mpi_integer,idest,1,mpi_comm_world,ierr)
      end do
    else
      do
        call mpi_send(idummy,1,mpi_integer,0,1,mpi_comm_world,ierr)
        call mpi_recv(num_loops,1,mpi_integer,0,1,mpi_comm_world,istatus,ierr)
        if(num_loops == 0) exit
#endif
        call Method(num_loops)
      end do
    end if
#ifdef MPI
    call mpi_barrier(mpi_comm_world, ierr)
    if(myrank == 0) then
      write(*,'(a)')"End: Master-slave procedure"
    end if
#endif
  end subroutine master_slave
end module MPIFunction
