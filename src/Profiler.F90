module Profiler
  use omp_lib
  use ClassSys, only: imap, dmap
  implicit none

  public :: timer

  type, private :: prof
    type(imap) :: counter
    type(dmap) :: timer
    type(imap) :: nmem
    type(dmap) :: memories
    real(8) :: start_time
    real(8) :: total_memory = 0.d0
    real(8) :: max_memory = 0.d0
    character(2) :: memory_unit = 'MB'
    character(:), allocatable :: OS
    character(9) :: tempfile = 'temp_prof'
  contains
    procedure :: init => InitProf
    procedure :: start => StartProf
    procedure :: fin => FinProf
    procedure :: add => AddProf
    procedure :: prt => PrintSummary
    procedure :: set_unit => SetMemoryUnit
    procedure :: cmemory => CurrentMemory
    procedure :: tmemory => TargetMemory
    procedure :: countup_memory => CountUpMemory
    procedure :: countdw_memory => CountDownMemory
    procedure :: prt_memory => PrintMemory
  end type prof

  type(prof) :: timer

contains
  subroutine InitProf(this, ut)
    class(prof), intent(inout) :: this
    character(*), intent(in), optional :: ut
    character(512) :: os_str
    integer :: io

    call system( 'uname > ' // trim(this%tempfile) )
    open(999,file=this%tempfile,iostat=io)
    if(io /= 0) then
      write(*,*) "Error opening file in InitProf"
      stop
    end if
    read(999,*,iostat=io) os_str
    if(io /= 0) then
      write(*,*) "Error reading file in InitProf"
      stop
    end if
    close(999)
    call system( 'rm ' // trim(this%tempfile) )

    select case(os_str)
    case("Linux", "linux")
      this%OS = "Linux"
    case("Darwin","darwin")
      this%OS = "OSX"
    case default
      write(*,'(2a)') "OS type is not clear: ", os_str
      this%OS = "unkown"
    end select

    this%start_time = omp_get_wtime()
    call this%counter%append('Total',1)
    call this%timer%append('Total',0.d0)
    if(present(ut)) call this%set_unit(ut)
    if(this%os == 'Linux') then
      call this%nmem%append('Total',1)
      call this%memories%append('Total',0.d0)
    end if
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
    if (myrank==0) write(*,*)
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

  subroutine SetMemoryUnit(this,ut)
    class(prof), intent(inout) :: this
    character(*), intent(in) :: ut
    this%memory_unit = ut
  end subroutine SetMemoryUnit

  subroutine CurrentMemory(this, msg)
    class(prof), intent(inout) :: this
    character(*), intent(in), optional :: msg
    integer :: io, getpid
    character(512) :: c_num, fn_proc, tmp1, tmp2, cmd
    real(8) :: mem

    if(this%os /= 'Linux') return
    write(c_num,'(i0)') getpid()
    fn_proc = '/proc/' // trim(c_num) // '/status'
    cmd = 'cat ' // trim(fn_proc) // ' | grep VmRSS > ' // trim(this%tempfile)
    call system( trim(cmd) )
    open(999,file=this%tempfile,iostat=io)
    if(io /= 0) then
      write(*,*) "File opening error in CurrentMemory"
      stop
    end if

    read(999,*,iostat=io) tmp1, mem, tmp2

    select case(this%memory_unit)
    case('kB')
      mem = mem
    case('MB')
      mem = mem / dble(1024)
    case('GB')
      mem = mem / dble(1024) ** 2
    end select

    this%max_memory = max(mem,this%max_memory)
    this%total_memory = mem

    if(io /= 0) then
      write(*,*) "File reading error in CurrentMemory"
      stop
    end if
    close(999)
    cmd = 'rm ' // trim(this%tempfile)
    call system( trim(cmd) )


    if(.not. present(msg) ) return
    write(*,'(a,f12.4,2a)') trim(msg), this%total_memory, &
        & " ", trim(this%memory_unit)
  end subroutine CurrentMemory

  ! call timer%cmemory before using
  ! msg is object name
  subroutine TargetMemory(this, msg)
    class(prof), intent(inout) :: this
    character(*), intent(in) :: msg
    real(8) :: mem_prev

    if(this%os /= 'Linux') return
    mem_prev = this%total_memory
    call this%cmemory()
    write(*,'(3a,f12.4,2a)') "Used memory for ", trim(msg), " is ", &
        & this%total_memory - mem_prev, " ", trim(this%memory_unit)
  end subroutine TargetMemory

  ! call timer%cmemory before using
  ! msg is object name
  subroutine CountUpMemory(this, key)
    class(prof), intent(inout) :: this
    character(*), intent(in) :: key
    integer :: n
    real(8) :: mem_prev, mem

    if(this%os /= 'Linux') return
    if(.not. this%memories%Search(key)) then
      call this%nmem%append(key,0)
      call this%memories%append(key,0.d0)
    end if
    n = this%nmem%Get(key)
    mem = this%memories%Get(key)

    mem_prev = this%total_memory
    call this%cmemory()
    call this%nmem%append(key,n+1)
    call this%memories%append(key,mem+this%total_memory-mem_prev)
  end subroutine CountUpMemory

  ! call timer%cmemory before using
  ! msg is object name
  subroutine CountDownMemory(this, key)
    class(prof), intent(inout) :: this
    character(*), intent(in) :: key
    integer :: n
    real(8) :: mem_prev, mem

    if(this%os /= 'Linux') return
    n = this%nmem%Get(key)
    mem = this%memories%Get(key)

    mem_prev = this%total_memory
    call this%cmemory()
    call this%nmem%append(key,n-1)
    call this%memories%append(key,mem-this%total_memory+mem_prev)
  end subroutine CountDownMemory

  subroutine PrintMemory(this)
    use MPIFunction, only: myrank
    class(prof), intent(inout) :: this
    integer :: i
    real(8) :: mtotal, a

    if(this%os /= 'Linux') return
    if (myrank==0) write(*,*)
    call this%cmemory()
    mtotal = 0.d0
    select case(this%memory_unit)
    case('kB')
      mtotal = this%max_memory / dble(1024) ** 2
    case('MB')
      mtotal = this%max_memory / dble(1024)
    case('GB')
      mtotal = this%max_memory
    end select

    if (myrank==0) write(*,'(a,f12.6,a)') &
        & '      Max used memory = ', mtotal, ' GB'
    if (myrank==0) write(*,*)
    if (myrank==0) write(*,'(30x,a)') &
        & "memory (GB),    ncall,   memory/ncall (GB)"
    do i = 1, this%nmem%GetLen()
      if(this%nmem%key(i) == 'Total') cycle

      select case(this%memory_unit)
      case('kB')
        a = this%memories%val(i) / dble(1024) ** 2
      case('MB')
        a = this%memories%val(i) / dble(1024)
      case('GB')
        a = this%memories%val(i)
      end select

      call PrintEach(this%nmem%key(i), this%nmem%val(i), a)
    end do
    if (myrank==0)  write(*,*)
  contains
    subroutine PrintEach(title, ncall, amnt)
      character(*), intent(in) :: title
      integer, intent(in) :: ncall
      real(8), intent(in) :: amnt


      if (myrank==0) &
          write(*,'(1a30, 1f12.3, 1i10, 1f20.5)') title,  &
          amnt, ncall, amnt/ncall
    end subroutine PrintEach

  end subroutine PrintMemory

  subroutine FinProf(this)
    class(prof), intent(inout) :: this
    call this%prt_memory()
    call this%prt()
  end subroutine FinProf

end module Profiler

!program test_Profiler
!  use omp_lib
!  use Profiler, only: timer
!  implicit none
!  integer :: i, j, k
!  real(8) :: ti, tj, tk
!  real(8) :: wa
!  real(8), allocatable :: a(:), b(:)
!
!  call timer%init()
!  allocate(a(10000000))
!  a = 0.d0
!  call timer%cmemory(.true.,'a is allocated : ')
!  ti = omp_get_wtime()
!  allocate(b(100000000))
!  call timer%add("Allocation ", omp_get_wtime() - ti)
!  ti = omp_get_wtime()
!  b = 0.d0
!  call timer%add("All zero ", omp_get_wtime() - ti)
!  call timer%tmemory('b')
!
!
!  deallocate(a,b)
!
!  call timer%fin()
!
!end program test_Profiler
