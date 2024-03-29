module HFInput
  use myfort
  implicit none

  public :: InputParameters

  private :: InitInputParameters
  private :: PrintInputParameters
  type :: InputParameters
    integer :: emax
    integer :: lmax
    integer :: e2max
    integer :: emax_mbpt
    integer :: lmax_mbpt
    integer :: e2max_mbpt
    integer :: e3max
    integer :: NOXB
    real(8) :: hw
    real(8) :: beta_cm
    real(8) :: alpha
    character(256), allocatable :: Ops(:)
    character(:), allocatable :: Nucl
    character(:), allocatable :: Core
    character(:), allocatable :: valence_list
    ! two-body file
    character(256) :: int_nn_file
    character(256), allocatable :: files_nn(:)
    integer :: emax_nn
    integer :: e2max_nn
    integer :: lmax_nn
    ! three-body file
    character(:), allocatable :: type_3n_file
    character(:), allocatable :: int_3n_file
    character(256), allocatable :: files_3n(:)
    integer :: emax_3n
    integer :: e2max_3n
    integer :: e3max_3n
    integer :: lmax_3n
    ! output files
    character(:), allocatable :: out_dir
    character(:), allocatable :: summary_file
    integer :: emax_1n
    integer :: lmax_1n
    character(:), allocatable :: density_matrix_file

    character(:), allocatable :: iter_method
    integer :: iter_n_history
    !
    logical :: is_Op_out
    logical :: is_MBPTscalar_full
    logical :: is_MBPT
    logical :: HO_reference
    logical :: is_NAT
    logical :: is_4th_order
    logical :: EN_denominator ! Epstein-Nesbet denominator (Moller-Plesset is default)
    logical :: dynamic_reference
    character(:), allocatable :: OpFileName
    character(:), allocatable :: TransFileName
    logical :: find_optimal_frequency

    ! atomic mode
    logical :: is_Atomic
    integer :: electron_number
  contains
    procedure :: init => InitInputParameters
    procedure :: PrintInputParameters
  end type InputParameters
contains
  subroutine InitInputParameters(this, inputfile)
    class(InputParameters), intent(inout) :: this
    character(*), intent(in) :: inputfile
    integer :: emax=6
    integer :: e2max=-1
    integer :: e3max=-1
    integer :: lmax=-1
    integer :: emax_mbpt=-1
    integer :: e2max_mbpt=-1
    integer :: lmax_mbpt=-1
    integer :: NOXB=2
    real(8) :: hw=20.d0
    real(8) :: beta_cm = 0.d0 ! lowson's cm parameter, H => H + beta_cm * Hcm in the unit of MeV
    real(8) :: alpha=0.7d0
    character(20) :: Nucl='O16'
    character(20) :: Core=""
    character(512) :: valence_list = ""
    character(256) :: int_nn_file="none"
    character(256) :: int_3n_file="none"
    character(256) :: type_3n_file="full"

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
    integer :: emax_1n=6
    integer :: lmax_1n=-1
    character(256) :: density_matrix_file = 'none'
    logical :: is_Op_out = .false.
    logical :: HO_reference = .false.
    logical :: is_MBPTscalar_full = .true.
    logical :: is_MBPT = .true.
    logical :: is_NAT = .false.
    logical :: is_4th_order = .false.
    logical :: EN_denominator =.false. ! Epstein-Nesbet denominator (Moller-Plesset is default)
    logical :: dynamic_reference=.false.
    character(256) :: OpFileName = "default"
    character(256) :: TransFileName = "none"
    logical :: find_optimal_frequency=.false.
    ! atomic mode
    logical :: is_Atomic=.false.
    character(512) :: iter_method = "linear"
    integer :: iter_n_history = 10

    type(sys) :: s
    integer :: io

    namelist /input/ emax, e2max, e3max, lmax, hw, &
        & Nucl, int_nn_file, files_nn, emax_nn, optrs, &
        & e2max_nn, lmax_nn, int_3n_file, files_3n, &
        & emax_3n, e2max_3n, e3max_3n, lmax_3n, alpha, &
        & summary_file, is_Op_out, is_MBPTscalar_full, &
        & is_MBPT, beta_cm, out_dir,&
        & OpFileName, is_Atomic, Core, valence_list, is_NAT, &
        & type_3n_file, is_4th_order, density_matrix_file, emax_1n, lmax_1n, &
        & EN_denominator, dynamic_reference, iter_method, iter_n_history, HO_reference, &
        & find_optimal_frequency, TransFileName, NOXB, emax_mbpt, e2max_mbpt, lmax_mbpt

    open(118, file=inputfile, action='read', iostat=io)
    if(io /= 0) then
      write(*,'(2a)') 'File opening error: ', trim(inputfile)
      stop
    end if
    read(118,nml=input)
    close(118)

    this%emax = emax
    this%e2max = e2max
    if(e2max == -1) this%e2max = 2*emax
    this%e3max = e3max
    if(e3max == -1) this%e3max = 3*emax
    this%emax_mbpt = emax_mbpt
    this%e2max_mbpt = e2max_mbpt
    this%lmax_mbpt = lmax_mbpt
    this%NOXB = NOXB
    this%hw = hw
    this%alpha = alpha
    this%Nucl = Nucl
    this%Core = Core
    if(this%Core=="") this%Core=this%Nucl
    this%valence_list = valence_list
    this%int_nn_file = int_nn_file
    this%int_3n_file = int_3n_file
    this%type_3n_file = type_3n_file

    this%emax_1n = emax_1n

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
    this%density_matrix_file = trim(density_matrix_file)

    this%HO_reference = HO_reference
    this%is_Op_out = is_Op_out
    this%is_MBPTscalar_full = is_MBPTscalar_full
    this%is_MBPT = is_MBPT
    this%is_NAT = is_NAT
    this%OpFileName = OpFileName
    this%TransFileName = TransFileName
    this%beta_cm = beta_cm
    this%is_Atomic = is_Atomic
    this%is_4th_order = is_4th_order
    this%EN_denominator = EN_denominator
    this%dynamic_reference = dynamic_reference
    this%iter_method = iter_method
    this%iter_n_history = iter_n_history
    this%find_optimal_frequency = find_optimal_frequency

    if(lmax == -1) this%lmax = emax
    if(lmax_1n == -1) this%lmax_1n = emax_1n
    if(lmax_nn == -1) this%lmax_nn = emax_nn
    if(lmax_3n == -1) this%lmax_3n = emax_3n
    if(emax_mbpt == -1) this%emax_mbpt = this%emax
    if(e2max_mbpt == -1) this%e2max_mbpt = 2*this%emax_mbpt
    if(lmax_mbpt == -1) this%lmax_mbpt = this%emax_mbpt

    call s%split(optrs, ',', this%Ops)
    call s%split(files_nn, ',', this%files_nn)
    call s%split(files_3n, ',', this%files_3n)

    if( size(this%Ops) /= size(this%files_nn) ) then
      write(*,*) "# Number of Op is not same as the number of NN files. Assuming all NN files are 'none'."
      deallocate(this%files_nn)
      allocate(this%files_nn(size(this%Ops)))
      this%files_nn(:) = 'none'
    end if

    if( size(this%Ops) /= size(this%files_3n) ) then
      write(*,*) "# Number of Op is not same as the number of NNN files. Assuming all NNN files are 'none'."
      deallocate(this%files_3n)
      allocate(this%files_3n(size(this%Ops)))
      this%files_3n(:) = 'none'
    end if
  end subroutine InitInputParameters

  subroutine PrintInputParameters(this,iunit)
    class(InputParameters), intent(in) :: this
    integer, intent(in), optional :: iunit
    integer :: iut = 6
    integer :: n

    if(present(iunit)) iut = iunit

    if(this%is_atomic) then
      write(iut,'(a)') "######  Input parameters  ####  "
      write(iut,'(a)') "#  Hamiltonian file:"
      write(iut,'(2a)') "#  ", trim(this%int_nn_file)
      return
    end if

    write(iut,'(a)') "######  Input parameters  ####  "
    write(iut,'(2a)') "#  Target nuclide is ", trim(this%Nucl)
    write(iut,'(a,i3,a,i3,a,i3,a,i3)') "#  Model space: emax =", &
        & this%emax, ", e2max =", this%e2max, &
        & ", e3max =", this%e3max, ", lmax =", this%lmax
    write(iut,'(a,f8.3,a)') "#  HO basis parameter hw = ", this%hw, " MeV"
    if(this%emax /= this%emax_mbpt) then
      write(iut,'(a)') "#"
      write(iut,'(a,i3,a,i3,a,i3)') "#  Model space for MBPT: emax =", &
          & this%emax_mbpt, ", e2max =", this%e2max_mbpt, &
          & ", lmax =", this%lmax_mbpt
    end if

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
    write(iut,'(4a)') "#  ", trim(this%int_3n_file), ", type: ", trim(this%type_3n_file)
    do n = 1, size(this%files_3n)
      if( this%files_3n(n) == 'none' .or. this%files_3n(n) == '') cycle
      write(iut,'(2a)') "#  ", trim(this%files_3n(n))
    end do
    write(iut,'(a,i1,a)') "#  3NMEs: NO",this%NOXB,"B approx. will be used."
    write(iut,'(a,i3,a,i3,a,i3,a,i3)') "#  File boundaries are emax =",this%emax_3n, &
        & ", e2max =",this%e2max_3n, ", e3max =", this%e3max_3n, &
        & ", lmax =",this%lmax_3n
    if(this%beta_cm > 1.d-4) write(iut,'(a,f6.3,a)') "#  Lawson's beta parameter: ", this%beta_cm, " MeV"
    do n = 1, size(this%Ops)
      if( this%Ops(n) == 'none' .or. this%Ops(n) == '') cycle
      write(iut,'(a,a)') "#  Operator is ", trim(this%Ops(n))
    end do

    if(this%is_MBPT) then
      write(iut,'(a)') "#  MBPT calc g.s. energy is done up to 3rd order"
      if(.not. this%EN_denominator) write(iut,'(a)') "#  Moller-Plesset (MP) denominator"
      if(      this%EN_denominator) write(iut,'(a)') "#  Epstein-Nesbet (EN) denominator"
      write(iut,'(a)') "#  MBPT calc g.s. scalar is done up to 2nd order"
      if(this%is_MBPTscalar_full) write(iut, '(a)') "#  MBPT for scalar operator is fully done up to 2nd order."
      if(.not. this%is_MBPTscalar_full) write(iut, '(a)') "#  MBPT for scalar operator is approximately done. "
    end if
  end subroutine PrintInputParameters
end module HFInput
