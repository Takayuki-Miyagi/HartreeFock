program BasisTransform
  use omp_lib
  use myfort
  use HFInput
  use ModelSpace
  use Operators
  use HartreeFock
  use WriteOperator
  implicit none
  type(InputParameters) :: par
  type(MSpace) :: ms, msNAT
  type(Ops) :: h, op_tr, op
  type(HFSolver) :: HF
  type(WriteFiles) :: w
  type(sys) :: s
  type(str) :: opname
  logical :: isfile
  type(str), allocatable :: splt(:)
  character(256) :: inputfile='none'
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

  call par%init(inputfile)
  call par%PrintInputParameters()
  if(par%int_nn_file /= 'none') isfile = s%isfile(par%int_nn_file, "main, 2N interaction")
  if(par%int_3n_file /= 'none') isfile = s%isfile(par%int_3n_file, "main, 3N interaction")
  do n = 1, size(par%Ops)
    if(par%files_nn(n) /= 'none') isfile = s%isfile(par%files_nn(n), "main, 2N file")
    if(par%files_3n(n) /= 'none') isfile = s%isfile(par%files_3n(n), "main, 3N file")
  end do

  call ms%init(Nucl=par%Nucl, Core=par%Core, valence_orbits=par%valence_list, &
    & hw=par%hw, emax=par%emax, e2max=par%e2max, lmax=par%lmax, beta=par%beta_cm)
  rank = 2
  call h%init('hamil',ms, rank, par%type_3n_file)
  write(*,*) __FILE__, " ", par%TransFileName
  call HF%ReadTransformationMatrix(h,par%TransFileName)

  do n = 1, size(par%Ops)
    if(par%Ops(n) == "") cycle
    write(*,*)
    write(*,'(3a)') "Transforming ", trim(par%Ops(n)), " operator"
    call op%init(par%Ops(n),ms,2)
    call op%set()
    op_tr = HF%BasisTransform(op,NOXB=par%NOXB,is_NO=.false.)
    call msNAT%init(Nucl=ms%Nucl, Core=ms%Core, hw=ms%hw, emax=par%emax_mbpt, e2max=par%e2max_mbpt, e3max=0)
    op = op_tr%truncate(msNAT)
    call w%SetFileName(par%OpFileName, op)
    call w%writef(par,op)
    op = op%NormalOrdering()
    write(*,'(a,f16.8)') "NO 0-body contribution: ", op%Zero
    call msNAT%fin()
  end do

  do n = 1, size(par%files_nn) ! format: JPTz^(filename)
    if(par%files_nn(n) == 'none') cycle
    write(*,*)
    write(*,'(4a)') "Transforming an operator fom ", trim(par%files_nn(n))
    opname = trim(par%files_nn(n))//trim("_file_")//trim(par%jpt_of_files(n))
    call op%init(opname%val,ms,2)
    op%oprtr = trim("file_")//trim(par%jpt_of_files(n))
    if(size(par%files_n) > n) call op%one%ReadOneBodyFile(par%files_n(n), par%emax, par%lmax)
    call op%set(par%files_nn(n), 'none', [par%emax_nn, par%e2max_nn,par%lmax_nn])
    op_tr = HF%BasisTransform(op,NOXB=par%NOXB,is_NO=.false.)
    call msNAT%init(Nucl=ms%Nucl, Core=ms%Core, hw=ms%hw, emax=par%emax_mbpt, e2max=par%e2max_mbpt, e3max=0)
    op = op_tr%truncate(msNAT)
    call w%SetFileName(par%OpFileName, op)
    call w%writef(par,op)
    op = op%NormalOrdering()
    write(*,'(a,f16.8)') "NO 0-body contribution: ", op%Zero
    call msNAT%fin()
  end do

  call timer%fin()

end program BasisTransform

