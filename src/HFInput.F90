module HFInput
  implicit none

  public :: InputParameters

  private :: InitInputParameters
  private :: PrintInputParameters
  type :: InputParameters
    integer :: emax
    integer :: lmax
    integer :: e2max
    integer :: e3max
    real(8) :: hw
    real(8) :: beta_cm
    real(8) :: alpha
    character(256), allocatable :: Ops(:)
    character(:), allocatable :: Nucl
    ! two-body file
    character(256) :: int_nn_file
    character(256), allocatable :: files_nn(:)
    integer :: emax_nn
    integer :: e2max_nn
    integer :: lmax_nn
    ! three-body file
    character(256) :: int_3n_file
    character(256), allocatable :: files_3n(:)
    integer :: emax_3n
    integer :: e2max_3n
    integer :: e3max_3n
    integer :: lmax_3n
    ! output files
    character(256) :: out_dir
    character(256) :: summary_file

    !
    logical :: is_Op_out
    logical :: is_MBPTscalar_full
    logical :: is_MBPTScalar
    logical :: is_MBPTEnergy
    character(256) :: Op_file_format

    ! atomic mode
    logical :: is_Atomic
    integer :: electron_number
  contains
    procedure :: init => InitInputParameters
    procedure :: PrintInputParameters
  end type InputParameters
contains
  subroutine InitInputParameters(this, inputfile)
    use ClassSys, only: sys
    class(InputParameters), intent(inout) :: this
    character(*), intent(in) :: inputfile
    integer :: emax=6
    integer :: e2max=12
    integer :: e3max=6
    integer :: lmax=-1
    real(8) :: hw=20.d0
    real(8) :: beta_cm = 0.d0 ! lowson's cm parameter, H => H + beta_cm * (hw/A) * Hcm
    real(8) :: alpha=1.d0
    character(20) :: Nucl='O16'
    character(256) :: int_nn_file
    character(256) :: int_3n_file

    character(1024) :: optrs="none"
    character(1024) :: files_nn="none"
    integer :: emax_nn=6
    integer :: e2max_nn=12
    integer :: lmax_nn=-1
    ! three-body file
    character(1024) :: files_3n="none"
    integer :: emax_3n=6
    integer :: e2max_3n=6
    integer :: e3max_3n=6
    integer :: lmax_3n=-1

    character(256) :: out_dir = '.'
    character(256) :: summary_file = 'summary.out'
    logical :: is_Op_out = .false.
    logical :: is_MBPTscalar_full = .false.
    logical :: is_MBPTScalar = .true.
    logical :: is_MBPTEnergy = .true.
    character(256) :: Op_file_format = "snt"
    ! atomic mode
    logical :: is_Atomic=.false.

    type(sys) :: s
    integer :: io
    namelist /input/ emax, e2max, e3max, lmax, hw, &
        & Nucl, int_nn_file, files_nn, emax_nn, optrs, &
        & e2max_nn, lmax_nn, int_3n_file, files_3n, &
        & emax_3n, e2max_3n, e3max_3n, lmax_3n, alpha, &
        & summary_file, is_Op_out, is_MBPTscalar_full, &
        & is_MBPTScalar, is_MBPTEnergy, beta_cm, out_dir,&
        & Op_file_format, is_Atomic

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
    this%hw = hw
    this%alpha = alpha
    this%Nucl = Nucl
    this%int_nn_file = int_nn_file
    this%int_3n_file = int_3n_file

    this%emax_nn = emax_nn
    this%e2max_nn = e2max_nn

    this%emax_3n = emax_3n
    this%e2max_3n = e2max_3n
    this%e3max_3n = e3max_3n

    this%lmax = lmax
    this%lmax_nn = lmax_nn
    this%lmax_3n = lmax_3n

    this%out_dir = out_dir
    if(this%out_dir == '') this%out_dir = '.'
    call s%mkdir(this%out_dir)
    this%summary_file = trim(this%out_dir) // "/" // trim(summary_file)

    this%is_Op_out = is_Op_out
    this%is_MBPTscalar_full = is_MBPTscalar_full
    this%is_MBPTEnergy = is_MBPTEnergy
    this%is_MBPTScalar = is_MBPTScalar
    this%Op_file_format = Op_file_format
    this%beta_cm = beta_cm
    this%is_Atomic = is_Atomic

    if(lmax == -1) this%lmax = emax
    if(lmax_nn == -1) this%lmax_nn = emax_nn
    if(lmax_3n == -1) this%lmax_3n = emax_3n

    call s%split(optrs, ',', this%Ops)
    call s%split(files_nn, ',', this%files_nn)
    call s%split(files_3n, ',', this%files_3n)

  end subroutine InitInputParameters

  subroutine PrintInputParameters(this,iunit)
    class(InputParameters), intent(in) :: this
    integer, intent(in), optional :: iunit
    integer :: iut = 6
    integer :: n

    if(present(iunit)) iut = iunit

    write(iut,'(a)') "######  Input parameters  ####  "
    write(iut,'(2a)') "#  Target nuclide is ", trim(this%Nucl)
    write(iut,'(a,i3,a,i3,a,i3,a,i3)') "#  Model space: emax =", &
        & this%emax, ", e2max =", this%e2max, &
        & ", e3max =", this%e3max, ", lmax =", this%lmax
    write(iut,'(a,f8.3,a)') "#  HO basis parameter hw = ", this%hw, " MeV"

    write(iut,'(a)') "#"
    write(iut,'(a)') "#  NN files:"
    write(iut,'(2a)') "#  ", trim(this%int_nn_file)
    do n = 1, size(this%files_nn)
      if( this%files_nn(n) == 'none' .or. this%files_nn(n) == '') cycle
      write(iut,'(2a)') "#  ", trim(this%files_nn(n))
    end do
    write(iut,'(a,i3,a,i3,a,i3)') "#  File boundaries are emax =",this%emax_nn, &
        & ", e2max =",this%e2max_nn, ", lmax =",this%lmax_nn

    write(iut,'(a)') "#"
    write(iut,'(a)') "#  3N files:"
    write(iut,'(2a)') "#  ", trim(this%int_3n_file)
    do n = 1, size(this%files_3n)
      if( this%files_3n(n) == 'none' .or. this%files_3n(n) == '') cycle
      write(iut,'(2a)') "#  ", trim(this%files_3n(n))
    end do
    write(iut,'(a,i3,a,i3,a,i3,a,i3)') "#  File boundaries are emax =",this%emax_3n, &
        & ", e2max =",this%e2max_3n, ", e3max =", this%e3max_3n, &
        & ", lmax =",this%lmax_3n
    if(this%beta_cm > 1.d-4) write(iut,'(a,f6.3)') "#  Lawson's beta parameter is ", this%beta_cm
    do n = 1, size(this%Ops)
      if( this%Ops(n) == 'none' .or. this%Ops(n) == '') cycle
      write(iut,'(a,a)') "#  Operator is ", trim(this%Ops(n))
    end do
    if(this%is_MBPTEnergy) write(iut,'(a)') "#  MBPT calc g.s. energy is done up to 3rd order"
    if(this%is_MBPTScalar) then
      write(iut,'(a)') "#  MBPT calc g.s. scalar is done up to 2nd order"
      if(this%is_MBPTscalar_full) write(iut, '(a)') "#  MBPT for scalar operator is fully done up to 2nd order."
      if(.not. this%is_MBPTscalar_full) write(iut, '(a)') "#  MBPT for scalar operator is approximately done. "
    end if
  end subroutine PrintInputParameters
end module HFInput
