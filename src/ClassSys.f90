module ClassSys
  implicit none

  public :: sys, iList, cList, dList, imap, cmap, dmap

  private :: mkdir, i2str, real2str, str2str
  private :: find, split, isfile_stop, isfile_func
  private :: GetLength_iList, append_iList, print_iList
  private :: GetLength_cList, append_cList, print_cList
  private :: GetLength_dList, append_dList, print_dList
  private :: GetLength_imap, append_imap, print_imap, Get_imap, Search_imap
  private :: GetLength_cmap, append_cmap, print_cmap, Get_cmap, Search_cmap
  private :: GetLength_dmap, append_dmap, print_dmap, Get_dmap, Search_dmap
  private :: search_key

  type :: sys
  contains
    procedure :: mkdir
    procedure :: i2str
    procedure :: real2str
    procedure :: str2str
    generic :: str => i2str, real2str, str2str
    procedure :: find
    procedure :: split
    procedure :: isfile_stop
    procedure :: isfile_func
    generic :: isfile => isfile_stop, isfile_func
  end type sys

  type :: iList
    integer, allocatable :: i(:)
  contains
    procedure :: GetLen => GetLength_ilist
    procedure :: append => append_ilist
    procedure :: prt => print_ilist
  end type iList

  integer, private, parameter :: Lenc = 256
  type :: cList
    character(Lenc), allocatable :: i(:)
  contains
    procedure :: GetLen => GetLength_clist
    procedure :: append => append_clist
    procedure :: prt => print_clist
  end type cList

  type :: dList
    real(8), allocatable :: i(:)
  contains
    procedure :: GetLen => GetLength_dlist
    procedure :: append => append_dlist
    procedure :: prt => print_dlist
  end type dList

  type :: imap
    integer, allocatable :: val(:)
    character(Lenc), allocatable :: key(:)
  contains
    procedure :: GetLen => GetLength_imap
    procedure :: append => append_imap
    procedure :: prt => print_imap
    procedure :: Get => Get_imap
    procedure :: Search => Search_imap
  end type imap

  type :: cmap
    character(Lenc), allocatable :: val(:)
    character(Lenc), allocatable :: key(:)
  contains
    procedure :: GetLen => GetLength_cmap
    procedure :: append => append_cmap
    procedure :: prt => print_cmap
    procedure :: Get => Get_cmap
    procedure :: Search => Search_cmap
  end type cmap

  type :: dmap
    real(8), allocatable :: val(:)
    character(Lenc), allocatable :: key(:)
  contains
    procedure :: GetLen => GetLength_dmap
    procedure :: append => append_dmap
    procedure :: prt => print_dmap
    procedure :: Get => Get_dmap
    procedure :: Search => Search_dmap
  end type dmap

contains

  subroutine mkdir(this, dir)
    class(sys), intent(in) :: this
    character(*), intent(in) :: dir
    character(Lenc) :: comm
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
    class(sys), intent(in) :: this
    character(:), allocatable :: i2str
    character(Lenc) :: ist
    integer, intent(in) :: i
    write(ist,*) i
    i2str=adjustl(trim(ist))
  end function i2str

  function real2str(this, r) result(str)
    class(sys), intent(in) :: this
    character(:), allocatable :: str
    character(Lenc) ::  rst
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
    class(sys), intent(in) :: this
    character(:), allocatable :: str2str
    character(*), intent(in) :: str
    str2str=adjustl(trim(str))
  end function str2str

  logical function find(this, str, key) result(r)
    class(sys), intent(in) :: this
    character(*), intent(in) :: str, key
    integer :: i
    i = index(str, key)
    if(i == 0) r = .false.
    if(i /= 0) r = .true.
  end function find

  subroutine split(this, str, key, splitted)
    class(sys), intent(in) :: this
    character(*), intent(in) :: str, key
    character(Lenc), allocatable, intent(out) :: splitted(:)
    integer :: len_key, n, i, n_elm
    integer, allocatable :: spos(:), epos(:)
    len_key = len(key)
    n_elm = 0
    n = 0
    do i = 1, len(str)-len_key+1
      if(str(i:i+len_key-1) == key) n = n + 1
    end do
    n_elm = n + 1
    allocate(splitted(n_elm))
    allocate(spos(n_elm))
    allocate(epos(n_elm))
    spos(1) = 1
    epos(n_elm) = len(str)
    n = 0
    do i = 1, len(str)-len_key+1
      if(str(i:i+len_key-1) == key) then
        n = n + 1
        spos(n+1) = i+len_key
        epos(n) = i-1
      end if
    end do

    do i = 1, n_elm-1
      splitted(i) = str(spos(i):epos(i))
    end do
    splitted(n_elm) = trim(str(spos(n_elm):))
    deallocate(spos, epos)
  end subroutine split

  function isfile_stop(this, f, msg) result(ex)
    class(sys), intent(in) :: this
    character(*), intent(in) :: f, msg
    logical :: ex

    inquire(file=f, exist=ex)
    if(.not. ex) then
      write(*,'(3a)') trim(f), ' is not found! In ', trim(msg)
      stop
    end if

  end function isfile_stop

  function isfile_func(this, f) result(ex)
    class(sys), intent(in) :: this
    character(*), intent(in) :: f
    logical :: ex
    inquire(file=f, exist=ex)
  end function isfile_func

  function GetLength_ilist(this) result(n)
    class(iList), intent(in) :: this
    integer :: n
    if(allocated(this%i)) then
      n = size(this%i)
    else
      n = 0
    end if
  end function GetLength_iList

  subroutine append_iList(this, i)
    class(iList), intent(inout) :: this
    integer, intent(in) :: i
    integer, allocatable :: list(:)
    integer :: n
    if(allocated(this%i)) then
      n = this%GetLen()
      allocate(list(n+1))
      list(1:n) = this%i(:)
      list(n+1) = i
      deallocate(this%i)
      allocate(this%i(n+1))
      this%i = list
      deallocate(list)
    else
      allocate(this%i(1))
      this%i(1) = i
    end if
  end subroutine append_iList

  subroutine print_iList(this, iunit)
    class(iList), intent(inout) :: this
    integer, optional, intent(in) :: iunit
    integer :: n, unt = 6
    if(present(iunit)) unt = iunit
    n = this%GetLen()
    if(n == 0) then
      write(unt,*) ''
    else
      write(unt,*) this%i
    end if
  end subroutine print_iList

  function GetLength_clist(this) result(n)
    class(cList), intent(in) :: this
    integer :: n
    if(allocated(this%i)) then
      n = size(this%i)
    else
      n = 0
    end if
  end function GetLength_cList

  subroutine append_cList(this, i)
    class(cList), intent(inout) :: this
    character(*), intent(in) :: i
    character(Lenc), allocatable :: list(:)
    integer :: n
    if(allocated(this%i)) then
      n = this%GetLen()
      allocate(list(n+1))
      list(1:n) = this%i(:)
      list(n+1) = i
      deallocate(this%i)
      allocate(this%i(n+1))
      this%i = list
      deallocate(list)
    else
      allocate(this%i(1))
      this%i(1) = i
    end if
  end subroutine append_cList

  subroutine print_cList(this, iunit)
    class(cList), intent(inout) :: this
    integer, optional, intent(in) :: iunit
    integer :: n, unt = 6, i
    if(present(iunit)) unt = iunit
    n = this%GetLen()
    if(n == 0) then
      write(unt,*) ''
    else
      do i = 1, n
        write(unt,'(2a)',advance='no') trim(this%i(i)), ' '
      end do
      write(unt,*)
    end if
  end subroutine print_cList

  function GetLength_dlist(this) result(n)
    class(dList), intent(in) :: this
    integer :: n
    if(allocated(this%i)) then
      n = size(this%i)
    else
      n = 0
    end if
  end function GetLength_dList

  subroutine append_dList(this, i)
    class(dList), intent(inout) :: this
    real(8), intent(in) :: i
    real(8), allocatable :: list(:)
    integer :: n
    if(allocated(this%i)) then
      n = this%GetLen()
      allocate(list(n+1))
      list(1:n) = this%i(:)
      list(n+1) = i
      deallocate(this%i)
      allocate(this%i(n+1))
      this%i = list
      deallocate(list)
    else
      allocate(this%i(1))
      this%i(1) = i
    end if
  end subroutine append_dList

  subroutine print_dList(this, iunit)
    class(dList), intent(inout) :: this
    integer, optional, intent(in) :: iunit
    integer :: n, unt = 6
    if(present(iunit)) unt = iunit
    n = this%GetLen()
    if(n == 0) then
      write(unt,*) ''
    else
      write(unt,*) this%i
    end if
  end subroutine print_dList

  function GetLength_imap(this) result(n)
    class(imap), intent(in) :: this
    integer :: n
    if(allocated(this%key)) then
      n = size(this%key)
    else
      n = 0
    end if
  end function GetLength_imap

  subroutine append_imap(this, key, val)
    class(imap), intent(inout) :: this
    character(*), intent(in) :: key
    integer, intent(in) :: val
    integer :: n, pos
    character(Lenc), allocatable :: keys(:)
    integer, allocatable :: vals(:)
    n = this%GetLen()
    if(n == 0) then
      allocate(this%key(1))
      allocate(this%val(1))
      this%key = key
      this%val = val
    else
      pos = search_key(this%key, key)
      if(pos == 0) then
        allocate(keys(n+1))
        allocate(vals(n+1))
        keys(1:n) = this%key(:)
        vals(1:n) = this%val(:)
        keys(n+1) = key
        vals(n+1) = val
        deallocate(this%key, this%val)
        allocate(this%key(n+1))
        allocate(this%val(n+1))
        this%key = keys
        this%val = vals
        deallocate(keys)
        deallocate(vals)
      else
        this%val(pos) = val
      end if
    end if
  end subroutine append_imap

  function Get_imap(this, key) result(r)
    class(imap), intent(in) :: this
    character(*), intent(in) :: key
    integer :: r
    integer :: n, pos
    n = this%GetLen()
    if(n == 0) then
      write(*,*) 'imap is empty'
      stop
    else
      pos = search_key(this%key, key)
      if(pos == 0) then
        write(*,'(2a)') 'key is not found, key:', trim(key)
        stop
      else
        r = this%val(pos)
      end if
    end if
  end function Get_imap

  subroutine print_imap(this, iunit)
    class(imap), intent(inout) :: this
    integer, optional, intent(in) :: iunit
    integer :: n, unt = 6, i
    if(present(iunit)) unt = iunit
    n = this%GetLen()
    if(n == 0) then
      write(unt,*) ''
    else
      write(unt,'(a)',advance='no') '{'
      do i = 1, n
        write(unt,'(2a,i5,a)',advance='no') trim(this%key(i)), ':', this%val(i), ', '
      end do
      write(unt,'(a)') '}'
    end if
  end subroutine print_imap

  function Search_imap(this, key) result(pos)
    class(imap), intent(in) :: this
    character(*), intent(in) :: key
    logical :: pos
    if(search_key(this%key, key) == 0) pos = .false.
    if(search_key(this%key, key) /= 0) pos = .true.
  end function Search_imap

  function GetLength_cmap(this) result(n)
    class(cmap), intent(in) :: this
    integer :: n
    if(allocated(this%key)) then
      n = size(this%key)
    else
      n = 0
    end if
  end function GetLength_cmap

  subroutine append_cmap(this, key, val)
    class(cmap), intent(inout) :: this
    character(*), intent(in) :: key
    character(*), intent(in) :: val
    integer :: n, pos
    character(Lenc), allocatable :: keys(:)
    character(Lenc), allocatable :: vals(:)
    n = this%GetLen()
    if(n == 0) then
      allocate(this%key(1))
      allocate(this%val(1))
      this%key = key
      this%val = val
    else
      pos = search_key(this%key, key)
      if(pos == 0) then
        allocate(keys(n+1))
        allocate(vals(n+1))
        keys(1:n) = this%key(:)
        vals(1:n) = this%val(:)
        keys(n+1) = key
        vals(n+1) = val
        deallocate(this%key, this%val)
        allocate(this%key(n+1))
        allocate(this%val(n+1))
        this%key = keys
        this%val = vals
        deallocate(keys)
        deallocate(vals)
      else
        this%val(pos) = val
      end if
    end if
  end subroutine append_cmap

  subroutine print_cmap(this, iunit)
    class(cmap), intent(inout) :: this
    integer, optional, intent(in) :: iunit
    integer :: n, unt = 6, i
    if(present(iunit)) unt = iunit
    n = this%GetLen()
    if(n == 0) then
      write(unt,*) ''
    else
      write(unt,'(a)',advance='no') '{'
      do i = 1, n
        write(unt,'(4a)',advance='no') trim(this%key(i)), ':', trim(this%val(i)), ', '
      end do
      write(unt,'(a)') '}'
    end if
  end subroutine print_cmap

  function Search_cmap(this, key) result(pos)
    class(cmap), intent(in) :: this
    character(*), intent(in) :: key
    logical :: pos
    if(search_key(this%key, key) == 0) pos = .false.
    if(search_key(this%key, key) /= 0) pos = .true.
  end function Search_cmap

  function Get_cmap(this, key) result(r)
    class(cmap), intent(in) :: this
    character(*), intent(in) :: key
    character(Lenc) :: r
    integer :: n, pos
    n = this%GetLen()
    if(n == 0) then
      write(*,*) 'imap is empty'
      stop
    else
      pos = search_key(this%key, key)
      if(pos == 0) then
        write(*,'(2a)') 'key is not found, key:', trim(key)
        stop
      else
        r = this%val(pos)
      end if
    end if
  end function Get_cmap

  function GetLength_dmap(this) result(n)
    class(dmap), intent(in) :: this
    integer :: n
    if(allocated(this%key)) then
      n = size(this%key)
    else
      n = 0
    end if
  end function GetLength_dmap

  subroutine append_dmap(this, key, val)
    class(dmap), intent(inout) :: this
    character(*), intent(in) :: key
    real(8), intent(in) :: val
    integer :: n, pos
    character(Lenc), allocatable :: keys(:)
    real(8), allocatable :: vals(:)
    n = this%GetLen()
    if(n == 0) then
      allocate(this%key(1))
      allocate(this%val(1))
      this%key = key
      this%val = val
    else
      pos = search_key(this%key, key)
      if(pos == 0) then
        allocate(keys(n+1))
        allocate(vals(n+1))
        keys(1:n) = this%key(:)
        vals(1:n) = this%val(:)
        keys(n+1) = key
        vals(n+1) = val
        deallocate(this%key, this%val)
        allocate(this%key(n+1))
        allocate(this%val(n+1))
        this%key = keys
        this%val = vals
        deallocate(keys)
        deallocate(vals)
      else
        this%val(pos) = val
      end if
    end if
  end subroutine append_dmap

  subroutine print_dmap(this, iunit)
    class(dmap), intent(inout) :: this
    integer, optional, intent(in) :: iunit
    integer :: n, unt = 6, i
    if(present(iunit)) unt = iunit
    n = this%GetLen()
    if(n == 0) then
      write(unt,*) ''
    else
      write(unt,'(a)',advance='no') '{'
      do i = 1, n
        write(unt,'(2a,es12.4,a)',advance='no') trim(this%key(i)), ':', this%val(i), ', '
      end do
      write(unt,'(a)') '}'
    end if
  end subroutine print_dmap

  function Get_dmap(this, key) result(r)
    class(dmap), intent(in) :: this
    character(*), intent(in) :: key
    real(8) :: r
    integer :: n, pos
    n = this%GetLen()
    if(n == 0) then
      write(*,*) 'imap is empty'
      stop
    else
      pos = search_key(this%key, key)
      if(pos == 0) then
        write(*,'(2a)') 'key is not found, key:', trim(key)
        stop
      else
        r = this%val(pos)
      end if
    end if
  end function Get_dmap

  function Search_dmap(this, key) result(pos)
    class(dmap), intent(in) :: this
    character(*), intent(in) :: key
    logical :: pos
    if(search_key(this%key, key) == 0) pos = .false.
    if(search_key(this%key, key) /= 0) pos = .true.
  end function Search_dmap

  function search_key(keys, key) result(n)
    character(*), intent(in) :: keys(:)
    character(*), intent(in) :: key
    integer :: i, n
    n = 0
    do i = 1, size(keys)
      if(keys(i) == key) then
        n = i
        return
      end if
    end do
  end function search_key

end module ClassSys

!program test
!  use ClassSys
!end program test
