module InputParameters
  use MPIFunction, only: myrank
  use class_sys, only: sy
  implicit none
  type(sy) :: sys
  type :: parameters
    real(8) :: hw

    ! two-body matrix element file
    character(256) :: twbmefile = 'None', scfile2 = 'None'

    ! three-body matrix element file
    integer :: emax_3nf, e2max_3nf, e3max_3nf, e3cut
    character(256) :: thbmefile = 'None'
    character(256) :: scfile3 = 'None'

    integer :: pmass ! Z
    integer :: nmass ! N
    integer :: mass  ! A
    logical :: HFloop ! Hartree-Fock calculation
    logical :: sv_hf_rslt
    logical :: MBPT
    integer :: emax, e2max, e3max ! model space
    real(8) :: conv ! tolerance
    character(256) :: vac = 'ref'
    character(256) :: egs = 'None', no2bhfile = 'None', fmt_hf_snt = 'None'
    character(256) :: reference = 'None', nocoef = 'None'

  contains
    procedure :: GetFileName
    procedure :: PrtParams
  end type parameters
contains
  subroutine set_parameters(params, hw_in, emax_in)
    type(parameters) :: params
    real(8), optional :: hw_in
    integer, optional :: emax_in
    integer :: pmass, nmass, mass
    real(8) :: hw, conv
    integer :: emax_2nf, e2max_2nf
    integer :: emax_3nf, e2max_3nf, e3max_3nf, e3cut
    character(256) :: inputfile, fmt_hf_snt, vac

    logical :: sv_hf_rslt, HFloop, NO2B, MBPT
    integer :: emax, e2max, e3max
    integer :: narg
    namelist /input/ pmass, nmass, mass, &
      & hw, conv, emax_2nf, e2max_2nf, &
      & emax_3nf, e2max_3nf, e3max_3nf, &
      & sv_hf_rslt, &
      & HFloop, NO2B, emax, e2max, &
      & e3cut, fmt_hf_snt, vac, MBPT, fmt_hf_snt

    narg = command_argument_count()
    call getarg(1, inputfile)
    open(118, file = inputfile, status = 'old')
    read(118, nml=input)
    close(118)
    if(narg >= 2) call getarg(2, params%reference)
    if(narg >= 3) call getarg(3, params%nocoef)
    if(narg >= 4) call getarg(4, params%twbmefile)
    if(narg >= 5) call getarg(5, params%thbmefile)
    if(narg >= 6) call getarg(6, params%scfile2)
    if(narg >= 7) call getarg(7, params%scfile3)

    params%pmass = pmass
    params%nmass = nmass
    params%mass  =  mass
    params%hw = hw
    params%conv = conv
    params%fmt_hf_snt = fmt_hf_snt
    params%emax_3nf = emax_3nf
    params%e2max_3nf = e2max_3nf
    params%e3max_3nf = e3max_3nf
    params%e3cut = e3cut
    params%HFloop = HFloop
    params%emax = emax
    params%e2max = e2max
    params%e3max = 18
    params%sv_hf_rslt = sv_hf_rslt
    params%vac = vac
    params%MBPT = MBPT
    params%fmt_hf_snt = fmt_hf_snt
    if(present(emax_in)) then
      params%emax = emax_in
      params%e2max = 2 * emax_in
    end if
    if(present(hw_in)) params%hw = hw_in
    call params%GetFileName()
    call params%PrtParams(6)
  end subroutine set_parameters

  subroutine GetFileName(params)
  class(parameters), intent(inout) :: params
    params%no2bhfile = 'None'
    params%no2bhfile = 'Hamil.snt.' // trim(params%fmt_hf_snt)
    params%egs = 'Summary.out'
  end subroutine GetFileName

  subroutine PrtParams(params, iunit)
  class(parameters), intent(in) :: params
    integer, intent(in) :: iunit
    ! Show Parameters
    if(myrank == 0) then
      write(iunit,'(a)') '! Calculation Parameters:'
      write(iunit,'(a, f6.2, a)') '! hw = ', params%hw, ' MeV'
      write(iunit,'(a,i3,a,i3,a,i3)') '! Z = ', params%pmass, &
      & ',  N = ', params%nmass, ',  A = ', params%mass
      write(iunit,'(a,i3,a,i3)') '! emax = ', params%emax, &
      & ',  e2max = ', params%e2max
      write(iunit,'(a)') '! Input Files:'
      write(iunit,'(2a)') '! reference-state file: ', trim(params%reference)
      write(iunit,'(2a)') '! normal-ordered-state file: ', trim(params%nocoef)
      write(iunit,'(2a)') '! 2BME file: ', trim(params%twbmefile)
      write(iunit,'(2a)') '! 3BME file: ', trim(params%thbmefile)
      write(iunit,'(2a)') '! 2BME scalar file: ', trim(params%scfile2)
      write(iunit,'(2a)') '! 3BME scalar file: ', trim(params%scfile3)
    end if
    write(iunit,'(a)') '! Output Files:'
    write(iunit,'(2a)') '! Ground-state energy: ', trim(params%egs)
    if(params%sv_hf_rslt) then
      write(iunit,'(2a)') '! NO2B file: ', trim(params%no2bhfile)
    end if
  end subroutine PrtParams
end module InputParameters

