module Profiler
  use omp_lib
  use ClassSys, only: imap, dmap
  implicit none

  public :: timer

  type, private :: prof
    type(imap) :: counter
    type(dmap) :: timer
    real(8) :: start_time
  contains
    procedure :: init => InitProf
    procedure :: start => StartProf
    procedure :: fin => FinProf
    procedure :: add => AddProf
    procedure :: prt => PrintSummary
  end type prof

  type(prof) :: timer

contains
  subroutine InitProf(this)
    class(prof), intent(inout) :: this
    this%start_time = omp_get_wtime()
    call this%counter%append('Total',1)
    call this%timer%append('Total',0.d0)
  end subroutine InitProf

  subroutine StartProf(this,key)
    class(prof), intent(inout) :: this
    character(*), intent(in) :: key
    call this%counter%append(key,0)
    call this%timer%append(key,0.d0)
  end subroutine StartProf

  subroutine AddProf(this,key,time_in)
    class(prof), intent(inout) :: this
    character(*), intent(in) :: key
    real(8), intent(in) :: time_in
    integer :: n
    real(8) :: time
    if(.not. this%counter%Search(key)) call this%start(key)
    n = this%counter%Get(key)
    time = this%timer%Get(key)
    call this%counter%append(key,n+1)
    call this%timer%append(key,time+time_in)
  end subroutine AddProf

  subroutine PrintSummary(this)
    use MPIFunction, only: myrank
    class(prof), intent(inout) :: this
    integer :: n, i
    real(8) :: r, ttotal
    call this%timer%append('Total',omp_get_wtime() - this%start_time)
    if (myrank==0) write(*,*)
    ttotal = this%timer%get('Total')
    n = int(ttotal)
    if (myrank==0) write(*,'(a,i5,a,i2.2,a,i2.2)') &
        '    summary of time, total = ', n/3600, ':', &
        mod(n, 3600)/60, ':', mod(n, 60)
    if(myrank == 0) write(*,*)
    if (myrank==0) write(*,'(37x,a)') "time,    ncall, time/ncall,   ratio "
    r = ttotal
    do i = 1, this%counter%GetLen()
      if(this%counter%key(i) == 'Total') cycle
      call PrintEach(this%counter%key(i), this%counter%val(i), this%timer%val(i))
      r = r - this%timer%val(i)
    end do
    if (myrank==0) write(*,'(1a30, 1f12.3, 22x, 1f9.4)') &
        "misc", r, r/ttotal
    if (myrank==0)  write(*,*)
  contains
    subroutine PrintEach(title, ncall, time)
      character(*), intent(in) :: title
      integer, intent(in) :: ncall
      real(8), intent(in) :: time
      if (myrank==0) &
          write(*,'(1a30, 1f12.3, 1i10, 1f12.5, 1f9.4)') title,  &
          time, ncall, time/ncall, &
          time/ttotal
    end subroutine PrintEach

  end subroutine PrintSummary

  subroutine FinProf(this)
    class(prof), intent(out) :: this
  end subroutine FinProf

end module Profiler

!program test_Profiler
!  use omp_lib
!  use Profiler, only: prof
!  implicit none
!  type(prof) :: timer
!  integer :: i, j, k
!  real(8) :: ti, tj, tk
!  real(8) :: wa
!  call timer%init()
!
!  call timer%Start('iloop')
!  call timer%Start('jloop')
!  call timer%Start('kloop')
!
!  wa = 0.d0
!  do i = 1, 10000
!    ti = omp_get_wtime()
!    wa = wa + dble(i)
!    call timer%Add('iloop',omp_get_wtime()-ti)
!  end do
!
!  ti = omp_get_wtime()
!  wa = 0.d0
!  do j = 1, 10000000
!    wa = wa + dble(j)
!  end do
!  call timer%Add('jloop',omp_get_wtime()-ti)
!
!  ti = omp_get_wtime()
!  wa = 0.d0
!  do k = 1, 100000000
!    wa = wa + dble(k)
!  end do
!  call timer%Add('kloop',omp_get_wtime()-ti)
!
!  ti = omp_get_wtime()
!  wa = 0.d0
!  do i = 1, 10000000
!    wa = wa + dble(i)
!  end do
!  call timer%Add('iloop',omp_get_wtime()-ti)
!  ti = omp_get_wtime()
!  wa = 0.d0
!  do i = 1, 10000000
!    wa = wa + dble(i)
!  end do
!  call timer%Add('iloop',omp_get_wtime()-ti)
!  ti = omp_get_wtime()
!  wa = 0.d0
!  do i = 1, 10000000
!    wa = wa + dble(i)
!  end do
!  call timer%Add('iloop',omp_get_wtime()-ti)
!
!  call timer%prt()
!
!end program test_Profiler
