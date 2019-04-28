program HFMain
  use omp_lib
  use Profiler, only: timer
  use HFInput
  use ModelSpace
  use Operators
  use HartreeFock
  use MBPT
  use WriteOperator
  use Atomic
  implicit none
  type(InputParameters) :: p
  type(MSpace) :: ms
  type(Ops) :: h, opr
  type(HFSolver) :: HF
  type(MBPTEnergy) :: PT
  type(MBPTScalar) :: PTs
  type(WriteFiles) :: w
  character(256) :: inputfile='none', conffile='none'
  integer :: n, istatus, wunit=23

  call timer%init()

  select case(command_argument_count())
  case(0)
    write(*,'(a)') "This code needs input file!"
    stop
  case(1)
    call get_command_argument(1,inputfile,status=istatus)
    write(*,'(2a)') "Input file: ", trim(inputfile)
  case(2)
    call get_command_argument(1,inputfile,status=istatus)
    call get_command_argument(2,conffile,status=istatus)
    write(*,'(4a)') "Input files: ", trim(inputfile), ", ", trim(conffile)
  case default
    write(*,'(a)') "Too many arguments!"
    stop
  end select

  call p%init(inputfile)
  call p%PrintInputParameters()

  if(command_argument_count() == 1 .and. p%is_Atomic) then
    write(*,"(a)") "I need a configuration file too!"
    stop
  end if

  if(command_argument_count() == 2 .and. p%is_Atomic) then
    call atomic_case(inputfile, conffile)
    call timer%fin()
    stop
  end if

  select case(p%int_3n_file)
  case('none', 'None', 'NONE')
    if(conffile == 'none') call ms%init(Nucl=p%Nucl, Core=p%Core, valence_orbits=p%valence_list, &
        & hw=p%hw, emax=p%emax, e2max=p%e2max, lmax=p%lmax, beta=p%beta_cm)
    if(conffile /= 'none') call ms%init(filename=conffile, hw=p%hw, emax=p%emax, e2max=p%e2max, lmax=p%lmax, beta=p%beta_cm)
  case default
    if(conffile == 'none') call ms%init(Nucl=p%Nucl, Core=p%Core, valence_orbits=p%valence_list, &
        & hw=p%hw, emax=p%emax, e2max=p%e2max, e3max=p%e3max, lmax=p%lmax, &
        & beta=p%beta_cm, is_three_body_jt=.true.)
    if(conffile /= 'none') call ms%init(filename=conffile, hw=p%hw, emax=p%emax, e2max=p%e2max, &
        & e3max=p%e3max, lmax=p%lmax, beta=p%beta_cm, is_three_body_jt=.true.)
  end select
  call w%init(p%emax, p%e2max)

  ! Hamiltonian -----
  select case(p%int_3n_file)
  case('none', 'None', 'NONE')
    call h%init('hamil',ms, 2)
  case default
    call h%init('hamil',ms, 3)
  end select
  call h%set(p%int_nn_file,p%int_3n_file,&
      & [p%emax_nn,p%e2max_nn,p%lmax_nn],&
      & [p%emax_3n,p%e2max_3n,p%e3max_3n,p%lmax_3n])

  call HF%init(h,alpha=p%alpha)
  call HF%solve()

  ! print single-particle energies
  !call HF%PrintSPEs(ms)
  call HF%TransformToHF(h)

  if(p%is_MBPTEnergy) then
    call PT%calc(H)
    open(wunit, file = p%summary_file, action='write',status='replace')
    call p%PrintInputParameters(wunit)
    write(wunit,'(a,f12.6)') "# max(| h / (e_h1 - e_p1) |)               = ", PT%perturbativity1b
    write(wunit,'(a,f12.6)') "# max(| v / (e_h1 + e_h2 - e_p1 - e_p2) |) = ", PT%perturbativity2b
    write(wunit,'(a,f12.6)') "# min( e_p ) - max( e_h ) = ", PT%energy_gap
    if(max(PT%perturbativity1b, PT%perturbativity2b) > 1.d0 .or. PT%energy_gap < 0.d0) then
      write(wunit,'(a)') "# MBPT might be dengerous! "
    end if
    write(wunit,'(a,6x,a,9x,a,9x,a,13x,a)') &
        & "# Operator", "HF energy", "2nd order", "3rd order", "Total"
    write(wunit,'(a,4f18.8)') 'hamil: ', PT%e_0, PT%e_2, PT%e_3, &
        & PT%e_0+PT%e_2+PT%e_3
    close(wunit)
  end if
  if(p%is_Op_out) then
    call w%SetFileName(p%out_dir, p%Op_file_format, h)
    call w%writef(p,h)
  end if
  ! Hamiltonian -----

  if(p%Ops(1) /= 'none' .or. p%Ops(1) /= '') then
    open(wunit, file = p%summary_file, action='write',status='old',position='append')
    write(wunit,'(a,1x,a,9x,a,9x,a,13x,a)') &
        & "# Operator", "HF exp. val.", "1st order", "2nd order", "Total"
    close(wunit)
  end if

  ! -- bare Ops --
  do n = 1, size(p%Ops)
    if(p%Ops(n) == 'none' .or. p%Ops(n) == "") cycle
    write(*,'(3a)') "## Calculating bare ", trim(p%Ops(n)), " operator"
    call opr%init(p%Ops(n),ms,2)
    call opr%set()
    call HF%TransformToHF(opr)
    if(p%is_MBPTScalar) then
      call PTs%calc(h,opr,p%is_MBPTScalar_full)
      open(wunit, file = p%summary_file, action='write',status='old',position='append')
      !write(wunit,'(3a)') "# Expectation value : <HF| ", trim(opr%optr)," |HF> "
      write(wunit,'(2a,4f18.8)') trim(p%Ops(n)), ": ", PTs%s_0, PTs%s_1, PTs%s_2, PTs%s_0+PTs%s_1+PTs%s_2
      close(wunit)
    end if
    if(p%is_Op_out) then
      call w%SetFileName(p%out_dir, p%Op_file_format, opr)
      call w%writef(p,opr)
    end if
    call opr%fin()
  end do

  ! -- Ops from NN file (srg evolved or two-body current) --
  do n = 1, size(p%files_nn)
    if(p%files_nn(n) == 'none' .or. p%Ops(n) == "") cycle
    write(*,'(4a)') "## Calculating ", trim(p%Ops(n)), " operator using ", trim(p%files_nn(n))
    call opr%init(p%Ops(n),ms, 2)
    call opr%set(p%files_nn(n), 'none', &
        & [p%emax_nn, p%e2max_nn,p%lmax_nn])
    call HF%TransformToHF(opr)
    if(p%is_MBPTScalar) then
      call PTs%calc(h,opr,p%is_MBPTScalar_full)
      open(wunit, file = p%summary_file, action='write',status='old',position='append')
      !write(wunit,'(3a)') "# Expectation value : <HF| ", trim(opr%optr)," |HF>:"
      !write(wunit,'(2a)') "# 2B file is ", trim(p%files_nn(n))
      write(wunit,'(2a, 4f18.8)') trim(p%Ops(n)), ": ", PTs%s_0, PTs%s_1, PTs%s_2, PTs%s_0+PTs%s_1+PTs%s_2
      close(wunit)
    end if
    if(p%is_Op_out) then
      call w%SetFileName(p%out_dir, p%Op_file_format, opr)
      call w%writef(p,opr)
    end if
    call opr%fin()
  end do

  ! -- Ops from NN+3N file (srg evolved) --
  do n = 1, size(p%files_3n)
    if(p%files_3n(n) == 'none' .or. p%Ops(n) == "") cycle
    write(*,'(4a)') "## Calculating ", trim(p%Ops(n)), " operator using ", trim(p%files_nn(n)), &
        & " and ", trim(p%files_3n(n))
    call opr%init(p%Ops(n),ms, 3)
    call opr%set(p%files_nn(n), p%files_3n(n), &
        & [p%emax_nn, p%e2max_nn,p%lmax_nn], &
        & [p%emax_3n,p%e2max_3n,p%e3max_3n,p%lmax_3n])
    call HF%TransformToHF(opr)
    if(p%is_MBPTScalar) then
      call PTs%calc(h,opr,p%is_MBPTScalar_full)
      open(wunit, file = p%summary_file, action='write',status='old',position='append')
      !write(wunit,'(3a)') "# Expectation value : <HF| ", trim(opr%optr)," |HF>:"
      !write(wunit,'(2a)') "# 2B file is ", trim(p%files_nn(n))
      !write(wunit,'(2a)') "# 3B file is ", trim(p%files_3n(n))
      !write(wunit,'(4f18.8)') PTs%s_0, PTs%s_1, PTs%s_2, PTs%s_0+PTs%s_1+PTs%s_2
      write(wunit,'(2a, 4f18.8)') trim(p%Ops(n)), ": ", PTs%s_0, PTs%s_1, PTs%s_2, PTs%s_0+PTs%s_1+PTs%s_2
      close(wunit)
    end if
    if(p%is_Op_out) then
      call w%SetFileName(p%out_dir, p%Op_file_format, opr)
      call w%writef(p,opr)
    end if
    call opr%fin()
  end do

  call h%fin()
  call HF%fin()
  call ms%fin()

  call timer%fin()
end program HFMain