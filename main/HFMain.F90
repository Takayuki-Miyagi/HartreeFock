program HFMain
  use omp_lib
  use myfort
  use HFInput
  use ModelSpace
  use Operators
  use ThreeBodyMonInteraction
  use HartreeFock
  use HFMBPT
  use WriteOperator
  use Atomic
  implicit none
  type(InputParameters) :: p
  type(MSpace) :: ms, msNAT
  type(Ops) :: h, opr, htr, Rho
  type(HFSolver) :: HF
  type(MBPTEnergy) :: PT
  type(MBPTScalar) :: PTs
  type(MBPTDMat) :: PTd
  type(WriteFiles) :: w
  type(sys) :: s
  logical :: isfile
  character(256) :: inputfile='none', conffile='none'
  integer :: n, istatus, wunit=23
  integer :: rank

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

  ! input files checking
  if(p%int_nn_file /= 'none') isfile = s%isfile(p%int_nn_file, "main, 2N interaction")
  if(p%int_3n_file /= 'none') isfile = s%isfile(p%int_3n_file, "main, 3N interaction")
  do n = 1, size(p%Ops)
    if(p%files_nn(n) /= 'none') isfile = s%isfile(p%files_nn(n), "main, 2N file")
    if(p%files_3n(n) /= 'none') isfile = s%isfile(p%files_3n(n), "main, 3N file")
  end do

  call w%init(p%emax, p%e2max)

  ! Model Space & Hamiltonian -----
  select case(p%int_3n_file)
  case('none', 'None', 'NONE')
    if(conffile == 'none') call ms%init(Nucl=p%Nucl, Core=p%Core, valence_orbits=p%valence_list, &
        & hw=p%hw, emax=p%emax, e2max=p%e2max, lmax=p%lmax, beta=p%beta_cm)
    if(conffile /= 'none') call ms%init(filename=conffile, hw=p%hw, emax=p%emax, e2max=p%e2max, lmax=p%lmax, beta=p%beta_cm)
    rank = 2
  case default
    if(p%type_3n_file=="full" .or. p%type_3n_file=="FULL") then
      if(conffile == 'none') call ms%init(Nucl=p%Nucl, Core=p%Core, valence_orbits=p%valence_list, &
          & hw=p%hw, emax=p%emax, e2max=p%e2max, e3max=p%e3max, lmax=p%lmax, &
          & beta=p%beta_cm, is_three_body_jt=.true.)
      if(conffile /= 'none') call ms%init(filename=conffile, hw=p%hw, emax=p%emax, e2max=p%e2max, &
          & e3max=p%e3max, lmax=p%lmax, beta=p%beta_cm, is_three_body_jt=.true.)
      rank = 3
    else
      if(conffile == 'none') call ms%init(Nucl=p%Nucl, Core=p%Core, valence_orbits=p%valence_list, &
          & hw=p%hw, emax=p%emax, e2max=p%e2max, e3max=p%e3max, lmax=p%lmax, beta=p%beta_cm)
      if(conffile /= 'none') call ms%init(filename=conffile, hw=p%hw, &
          & emax=p%emax, e2max=p%e2max, e3max=p%e3max, lmax=p%lmax, beta=p%beta_cm)
      rank = 3
    end if
  end select

  call h%init('hamil',ms, rank, p%type_3n_file)

  call h%set(p%int_nn_file,p%int_3n_file,&
        & [p%emax_nn,p%e2max_nn,p%lmax_nn],&
        & [p%emax_3n,p%e2max_3n,p%e3max_3n,p%lmax_3n])

  if( p%HO_reference ) then
    open(wunit, file = p%summary_file, action='write',status='replace')
    if( p%TransFileName == "none") then
      call HF%init(h,alpha=p%alpha,NOXB_3NF=p%NOXB)
      htr = HF%BasisTransform(h,NOXB=p%NOXB)
    end if
    if( p%TransFileName /= "none") then ! Transformation file reading mode, not HO
      call HF%ReadTransformationMatrix(h,p%TransFileName)
      htr = HF%BasisTransform(h,NOXB=p%NOXB)
    end if
    if( p%is_MBPT) then
      call PT%calc(htr, p%is_4th_order, p%EN_denominator)
      write(wunit,'(a,11x,a,6x,a,9x,a,9x,a,13x,a)') &
          & "#", "Operator", "Ref. exp. val.", "MBPT LO", "MBPT NLO", "Total"
      write(wunit,'(a20,4f18.8)') 'hamil', PT%e_0, PT%e_2, PT%e_3, &
          & PT%e_0+PT%e_2+PT%e_3
      close(wunit)
    else
      write(wunit,'(a,11x,a,6x,a)') &
          & "#", "Operator", "Ref. exp. val."
      write(wunit,'(a20,f18.8)') 'hamil', htr%zero
      close(wunit)
    end if
  end if

  if( .not. p%HO_reference ) then

    call HF%init(h,alpha=p%alpha,NOXB_3NF=p%NOXB)
    call HF%SetDynamicReference(p%dynamic_reference)
    call HF%SetIterMethod(p%iter_method, p%alpha, p%iter_n_history)
    call HF%solve()
    if( p%TransFileName /= "none") call HF%WriteTransformationMatrix(p%TransFileName)

    open(wunit, file = p%summary_file, action='write',status='replace')
    if( p%find_optimal_frequency ) then
      call HF%FindOptimalFrequency()
    end if

    if( .not. p%find_optimal_frequency ) then
      call p%PrintInputParameters(wunit)
      call HF%PrintSPEs(wunit)
      if(HF%fail) then
        write(wunit,"(a)") "#*************************************************"
        write(wunit,"(a)") "#*************************************************"
        write(wunit,"(a)") "# Hartree-Fock doesn't converge"
        write(*,"(a)") "#*************************************************"
        write(*,"(a)") "#*************************************************"
        write(*,"(a)") "# Hartree-Fock doesn't converge"
        stop
      end if

      if(.not. p%is_MBPT) then
        write(wunit,'(a,6x,a)') "# Operator", "Ref. exp. val."
        write(wunit,'(a,1f18.8)') 'hamil: ', HF%ehf
      end if

      if(p%is_MBPT) then
        htr = HF%BasisTransform(h,NOXB=p%NOXB)
        call PT%calc(htr, p%is_4th_order, p%EN_denominator)
        write(wunit,'(a,f12.6)') "# max(| h / (e_h1 - e_p1) |)               = ", PT%perturbativity1b
        write(wunit,'(a,f12.6)') "# max(| v / (e_h1 + e_h2 - e_p1 - e_p2) |) = ", PT%perturbativity2b
        write(wunit,'(a,f12.6)') "# min( e_p ) - max( e_h ) = ", PT%energy_gap
        if(max(PT%perturbativity1b, PT%perturbativity2b) > 1.d0 .or. PT%energy_gap < 0.d0) then
          write(wunit,'(a)') "# MBPT might be dengerous! "
        end if
        write(wunit,'(a,11x,a,6x,a,9x,a,9x,a,13x,a)') &
            & "#", "Operator", "Ref. exp. val.", "MBPT LO", "MBPT NLO", "Total"
        write(wunit,'(a20,4f18.8)') 'hamil', PT%e_0, PT%e_2, PT%e_3, &
            & PT%e_0+PT%e_2+PT%e_3

        if(p%is_Op_out) then ! This can be an input of an NCSM calculation
          call w%SetFileName(p%OpFileName, htr)
          call w%writef(p,htr)
        end if
      end if
    end if
    close(wunit)
  end if

  if(p%is_NAT) then
    htr = HF%BasisTransform(h,NOXB=p%NOXB)
    call PTd%init(HF, htr, p%EN_denominator)
    HF%C = PTd%C_HO2NAT
    if( p%TransFileName /= "none") call HF%WriteTransformationMatrix(p%TransFileName)
    select case(p%density_matrix_file)
    case("", "none", "NONE", "None")
      if(p%is_Op_out) then
        call Rho%init("DenMat", ms, 2)
        Rho%one = PTd%rho_HO
        call w%SetFileName(p%OpFileName, Rho)
        call w%writef(p,Rho)
      end if
    case default
      ! reading density matri file from M. Gennari
      call PTd%rho_HO%ReadOneBodyFile(p%density_matrix_file, p%emax_1n, p%lmax_1n)
      call PTd%GetCoef()
      HF%C = PTd%C_HO2NAT
      if(p%is_Op_out) then
        call Rho%init("DenMat", ms, 2)
        Rho%one = PTd%rho_HO
        call w%SetFileName(p%OpFileName, Rho)
        call w%writef(p,Rho)
      end if
    end select
    htr = HF%BasisTransform(h,NOXB=p%NOXB)
    call htr%UnNormalOrdering2B()

    call msNAT%init(Nucl=ms%Nucl, Core=ms%Core, hw=ms%hw, emax=p%emax_mbpt, e2max=p%e2max_mbpt, e3max=0)
    h = htr%truncate(msNAT)
    if(p%is_Op_out) then ! This can be an input of an NCSM calculation
      call w%SetFileName(p%OpFileName, h)
      call w%writef(p,h)
    end if

    call HF%fin()
    call HF%init(h,alpha=p%alpha)
    call HF%solve()
    if(p%is_MBPT) then
        htr = HF%BasisTransform(h,NOXB=p%NOXB)
        call PT%calc(htr, p%is_4th_order, p%EN_denominator)
        write(wunit,'(a,f12.6)') "# max(| h / (e_h1 - e_p1) |)               = ", PT%perturbativity1b
        write(wunit,'(a,f12.6)') "# max(| v / (e_h1 + e_h2 - e_p1 - e_p2) |) = ", PT%perturbativity2b
        write(wunit,'(a,f12.6)') "# min( e_p ) - max( e_h ) = ", PT%energy_gap
        if(max(PT%perturbativity1b, PT%perturbativity2b) > 1.d0 .or. PT%energy_gap < 0.d0) then
          write(wunit,'(a)') "# MBPT might be dengerous! "
        end if
        write(wunit,'(a,11x,a,6x,a,9x,a,9x,a,13x,a)') &
            & "#", "Operator", "Ref. exp. val.", "MBPT LO", "MBPT NLO", "Total"
        write(wunit,'(a20,4f18.8)') 'hamil', PT%e_0, PT%e_2, PT%e_3, &
            & PT%e_0+PT%e_2+PT%e_3
      end if

  end if
  call h%fin()

  ! Hamiltonian -----

  ! -- bare Ops --
  do n = 1, size(p%Ops)
    if(p%Ops(n) == 'none' .or. p%Ops(n) == "" .or. s%find(p%Ops(n), "_file")) cycle
    write(*,*)
    write(*,'(3a)') "## Calculating bare ", trim(p%Ops(n)), " operator"
    call opr%init(p%Ops(n),ms,2)
    call opr%set()
    opr = HF%BasisTransform(opr,NOXB=p%NOXB)
    if(p%is_MBPT) then
      if(s%find(p%Ops(n), 'Hamil') .or. p%Ops(n) == "Tkin") then
        call PTs%calc(htr,opr,p%is_MBPTScalar_full, p%EN_denominator, part_of_hamil=.true.)
      else
        call PTs%calc(htr,opr,p%is_MBPTScalar_full, p%EN_denominator)
      end if
      open(wunit, file = p%summary_file, action='write',status='old',position='append')
      !write(wunit,'(3a)') "# Expectation value : <HF| ", trim(opr%optr)," |HF> "
      write(wunit,'(a20, 4f18.8)') trim(p%Ops(n)), PTs%s_0, PTs%s_1, PTs%s_2, PTs%s_0+PTs%s_1+PTs%s_2
      close(wunit)
    else
      open(wunit, file = p%summary_file, action='write',status='old',position='append')
      write(wunit,'(a20, f18.8)') trim(p%Ops(n)), opr%zero
      close(wunit)
    end if
    if(p%is_Op_out) then
      call w%SetFileName(p%OpFileName, opr)
      call w%writef(p,opr)
    end if
    call opr%fin()
  end do

  ! -- Ops from NN file (srg evolved or two-body current) --
  do n = 1, size(p%Ops)
    if(n > size(p%files_nn)) cycle
    if(p%files_nn(n) == 'none' .or. p%Ops(n) == "") cycle
    write(*,*)
    write(*,'(4a)') "## Calculating ", trim(p%Ops(n)), " operator using ", trim(p%files_nn(n))
    call opr%init(p%Ops(n),ms, 2)
    call opr%set(p%files_nn(n), 'none', &
        & [p%emax_nn, p%e2max_nn,p%lmax_nn])
    opr = HF%BasisTransform(opr,NOXB=p%NOXB)
    if(p%is_MBPT) then
      if(s%find(p%Ops(n), 'Hamil')) call PTs%calc(htr,opr,p%is_MBPTScalar_full, p%EN_denominator, part_of_hamil=.true.)
      if(.not. s%find(p%Ops(n), 'Hamil')) call PTs%calc(htr,opr,p%is_MBPTScalar_full, p%EN_denominator)
      open(wunit, file = p%summary_file, action='write',status='old',position='append')
      !write(wunit,'(3a)') "# Expectation value : <HF| ", trim(opr%optr)," |HF>:"
      !write(wunit,'(2a)') "# 2B file is ", trim(p%files_nn(n))
      write(wunit,'(a20, 4f18.8)') trim(p%Ops(n)), PTs%s_0, PTs%s_1, PTs%s_2, PTs%s_0+PTs%s_1+PTs%s_2
      close(wunit)
    else
      open(wunit, file = p%summary_file, action='write',status='old',position='append')
      write(wunit,'(a20, f18.8)') trim(p%Ops(n)), opr%zero
      close(wunit)
    end if
    if(p%is_Op_out) then
      call w%SetFileName(p%OpFileName, opr)
      call w%writef(p,opr)
    end if
    call opr%fin()
  end do

  ! -- Ops from NN+3N file (srg evolved) --
  do n = 1, size(p%Ops)
    if(n > size(p%files_nn)) cycle
    if(n > size(p%files_3n)) cycle
    if(p%files_3n(n) == 'none' .or. p%Ops(n) == "") cycle
    write(*,*)
    write(*,'(6a)') "## Calculating ", trim(p%Ops(n)), " operator using ", trim(p%files_nn(n)), &
        & " and ", trim(p%files_3n(n))
    call opr%init(p%Ops(n),ms, 3, p%type_3n_file)
    call opr%set(p%files_nn(n), p%files_3n(n), &
        & [p%emax_nn, p%e2max_nn,p%lmax_nn], &
        & [p%emax_3n,p%e2max_3n,p%e3max_3n,p%lmax_3n])
    opr = HF%BasisTransform(opr,NOXB=p%NOXB)
    if(p%is_MBPT) then
      if(s%find(p%Ops(n), 'Hamil')) call PTs%calc(htr,opr,p%is_MBPTScalar_full, p%EN_denominator, part_of_hamil=.true.)
      if(.not. s%find(p%Ops(n), 'Hamil')) call PTs%calc(htr,opr,p%is_MBPTScalar_full, p%EN_denominator)
      open(wunit, file = p%summary_file, action='write',status='old',position='append')
      !write(wunit,'(3a)') "# Expectation value : <HF| ", trim(opr%optr)," |HF>:"
      !write(wunit,'(2a)') "# 2B file is ", trim(p%files_nn(n))
      !write(wunit,'(2a)') "# 3B file is ", trim(p%files_3n(n))
      !write(wunit,'(4f18.8)') PTs%s_0, PTs%s_1, PTs%s_2, PTs%s_0+PTs%s_1+PTs%s_2
      write(wunit,'(a20, 4f18.8)') trim(p%Ops(n)), PTs%s_0, PTs%s_1, PTs%s_2, PTs%s_0+PTs%s_1+PTs%s_2
      close(wunit)
    else
      open(wunit, file = p%summary_file, action='write',status='old',position='append')
      write(wunit,'(a20, f18.8)') trim(p%Ops(n)), opr%zero
      close(wunit)
    end if
    if(p%is_Op_out) then
      call w%SetFileName(p%OpFileName, opr)
      call w%writef(p,opr)
    end if
    call opr%fin()
  end do

  call htr%fin()
  call HF%fin()
  call ms%fin()

  call timer%fin()
end program HFMain
