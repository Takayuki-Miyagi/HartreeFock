program HFMain
  use omp_lib
  use myfort
  use HFInput
  use ModelSpace
  use Operators
  use ThreeBodyMonInteraction
  use HartreeFock
  use gMBPT
  use WriteOperator
  use HartreeFock
  implicit none
  type(InputParameters) :: p
  type(MSpace) :: ms
  type(Ops) :: h, htr
  type(HFSolver) :: HF
  type(gMBPTEnergy) :: PT
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
  case default
    write(*,'(a)') "Too many arguments!"
    stop
  end select

  call p%init(inputfile)
  call p%PrintInputParameters()

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

  open(wunit, file = p%summary_file, action='write',status='replace')
  call HF%init(h,alpha=p%alpha,NOXB_3NF=p%NOXB)
  htr = HF%BasisTransform(h,NOXB=p%NOXB)
  call PT%calc(htr, p%EN_denominator)
  write(wunit,'(a,11x,a,6x,a,9x,a,9x,a,13x,a)') &
    & "#", "Operator", "Ref. exp. val.", "MBPT LO", "MBPT NLO", "Total"
  write(wunit,'(a20,4f18.8)') 'hamil', PT%e_0, sum(PT%e_2), sum(PT%e_3), &
    & PT%e_0+sum(PT%e_2)+sum(PT%e_3)
  close(wunit)

  call HF%fin()
  call ms%fin()
  call timer%fin()
end program HFMain
