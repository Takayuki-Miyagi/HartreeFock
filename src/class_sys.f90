module class_sys
  type :: sy
  contains
    procedure :: mkdir
    generic :: str => i2str, real2str, str2str
    procedure :: i2str
    procedure :: real2str
    procedure :: str2str
    procedure :: find
  end type sy
contains
  subroutine mkdir(this, dir)
  class(sy) :: this
    character(*), intent(in) :: dir
    character(len=1000) :: comm
    integer :: is, ie, idxdir
    integer :: lnblnk

    idxdir = lnblnk(dir)
    is = 1; ie = is + len('if [ ! -d ')
    write(comm(is:ie), '(a)') 'if [ ! -d '
    is = ie + 1; ie = is + idxdir
    write(comm(is:ie), '(a)') dir(1:idxdir)
    is = ie + 1
    ie = is + len(' ]; then (echo "Creating directory, ')
    write(comm(is:ie), '(a)') ' ]; then (echo "Creating directory, '
    is = ie + 1; ie = is + idxdir
    write(comm(is:ie), '(a)') dir(1:idxdir)
    is = ie + 1; ie = is + len('"; mkdir -p ') - 1
    write(comm(is:ie), '(a)') '"; mkdir -p '
    is = ie + 1; ie = is + idxdir - 1
    write(comm(is:ie), '(a)') dir(1:idxdir)
    is = ie + 1; ie = is + 4
    write(comm(is:ie), '(a)') '); fi'
    !write(*,'(a)') comm(1:ie)
    call system(comm(1:ie))
  end subroutine mkdir

  function i2str(this, i)
  class(sy) :: this
    character(:), allocatable :: i2str
    character(256) :: ist
    integer, intent(in) :: i
    write(ist, *) i
    i2str=adjustl(trim(ist))
  end function i2str

  function real2str(this, r) result(str)
  class(sy) :: this
    character(:), allocatable :: str
    character(256) :: rst
    real(8), intent(in) :: r
    integer :: l
    l = int(log10(abs(r) + 1.d-8))
    if(l >= 1) then
      write(rst, *) int(r)
    elseif(l < 1 .and. l >= 0) then
      write(rst, '(f10.2)') r
    elseif(l < 0 .and. l >= -1) then
      write(rst, '(f10.3)') r
    elseif(l < -1 .and. l >= -2) then
      write(rst, '(f10.4)') r
    elseif(l < -2 .and. l >= -3) then
      write(rst, '(f10.5)') r
    elseif(l < -3 .and. l >= -4) then
      write(rst, '(f10.6)') r
    elseif(l < -4 .and. l >= -5) then
      write(rst, '(f10.7)') r
    elseif(l < -5 .and. l >= -6) then
      write(rst, '(f10.8)') r
    elseif(l < -6 .and. l >= -7) then
      write(rst, '(f10.9)') r
    end if
    str = adjustl(trim(rst))
  end function real2str

  function str2str(this, str)
  class(sy) :: this
    character(:), allocatable :: str2str
    character(*), intent(in) :: str
    str2str=adjustl(trim(str))
  end function str2str

  logical function find(this, str, key) result(r)
  class(sy) :: this
    character(*), intent(in) :: str, key
    integer :: i
    i = index(str, key)
    if(i == 0) r = .false.
    if(i /= 0) r = .true.
  end function find
end module class_sys

!program test
!  use class_sys
!  type(sy) :: sys
!  character(:), allocatable :: a
!  a = sys%str(100)
!  write(*,*) a
!  a = sys%str(10.d0)
!  write(*,*) a
!  a = sys%str(0.1d0)
!  write(*,*) a
!  a = sys%str('abc')
!  write(*,*) a
!end program test
