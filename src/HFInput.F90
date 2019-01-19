module HFInput
  implicit none
  type :: InputParameters
    integer :: emax
    integer :: lmax
    integer :: e2max
    integer :: e3max
    real(8) :: hw
    character(32), allocatable :: Ops(:)
    character(:), allocatable :: Nucl
    ! two-body file
    character(512), allocatable :: files_nn(:)
    integer, allocatable :: emaxs_nn(:)
    integer, allocatable :: e2maxs_nn(:)
    integer, allocatable :: lmaxs_nn(:)
    ! three-body file
    character(512), allocatable :: files_3n(:)
    integer, allocatable :: emaxs_3n(:)
    integer, allocatable :: e2maxs_3n(:)
    integer, allocatable :: e3maxs_3n(:)
    integer, allocatable :: lmaxs_3n(:)
  contains
    procedure :: init => InitInputParameters
  end type InputParameters
contains
  subroutine InitInputParameters(this)
    class(InputParameters), intent(inout) :: this
    integer :: emax=6
    integer :: e2max=12
    integer :: e3max=6
    integer :: lmax=6
    real(8) :: hw=20.d0
    character(20) :: Nucl='O16'
    character(:), allocatable :: Operators
    ! two-body file
    character(:), allocatable :: files_nn
    character(:), allocatable :: emaxs_nn
    character(:), allocatable :: e2maxs_nn
    character(:), allocatable :: lmaxs_nn
    ! three-body file
    character(:), allocatable :: files_3n(:)
    character(:), allocatable :: emaxs_3n(:)
    character(:), allocatable :: e2maxs_3n(:)
    character(:), allocatable :: e3maxs_3n(:)
    character(:), allocatable :: lmaxs_3n(:)

    character(:), allocatable :: inputfile
    integer :: io
    namelist /input/ emax, e2max, e3max, lmax, hw, &
        & Nucl, Operators, files_nn, emaxs_nn, &
        & e2maxs_nn, lmaxs_nn, files_3n, &
        & emaxs_3n, e2maxs_3n, e3maxs_3n, lmaxs_3n

    call getarg(1,inputfile)
    open(118, file=inputfile, action='read', iostat=io)
    if(io /= 0) then
      write(*,'(2a)') 'File opening error: ', trim(inputfile)
      stop
    end if
    read(118,nml=input)
    close(118)

    this%emax = emax
    this%e2max = e2max
    this%e3max = e3max
    this%lmax = lmax
    this%hw = hw
    this%Nucl = Nucl




  end subroutine InitInputParameters
end module HFInput
