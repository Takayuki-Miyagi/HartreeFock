program HFMain
  use omp_lib
  use Profiler, only: timer
  use CommonLibrary, only: &
      &init_dbinomial_triangle, fin_dbinomial_triangle
  use HFInput
  use ModelSpace
  use Operators
  use HartreeFock
  use MBPT
  use WriteOperator
  implicit none
  type(InputParameters) :: p
  type(MSpace) :: ms
  type(Op) :: h, opr
  type(HFSolver) :: HF
  type(MBPTEnergy) :: PT
  type(WriteFiles) :: w
  character(256) :: inputfile
  real(8) :: ti
  integer :: n

  call timer%init()

  call getarg(1,inputfile)
  call p%init(inputfile)
  call p%PrintInputParameters()

  ti = omp_get_wtime()
  call init_dbinomial_triangle()
  call timer%Add("Init for Rotation Group",omp_get_wtime()-ti)

  if(p%e3max == 0) then
    call ms%init(p%Nucl, p%hw, p%emax, p%e2max)
  elseif(p%e3max > 0) then
    call ms%init(p%Nucl, p%hw, p%emax, p%e2max, e3max=p%e3max)
  end if
  call w%init(p%emax, p%e2max)

  ! Hamiltonian -----
  call h%init('hamil',ms,.true.)
  call h%set(ms,p%int_nn_file,p%int_3n_file,&
      & [p%emax_nn,p%e2max_nn,p%lmax_nn],&
      & [p%emax_3n,p%e2max_3n,p%e3max_3n,p%lmax_3n])

  call HF%init(ms,h,alpha=p%alpha)
  call HF%solve(ms%sps,ms%one)
  call HF%TransformToHF(ms,h)


  call PT%calc(ms,H)
  call PT%write_summary(p)
  !call w%writef(p,ms,h)
  call h%fin()
  ! Hamiltonian -----

  ! -- bare Ops --
  do n = 1, size(p%Ops)
    if(p%Ops(n) == 'none' .or. p%Ops(n) == "") cycle
    call opr%init(p%Ops(n),ms,.false.)
    call opr%set(ms)
    call HF%TransformToHF(ms,opr)
    call w%writef(p,ms,opr)
  end do

  ! -- Ops from NN file (srg evolved or two-body current) --
  do n = 1, size(p%files_nn)
    if(p%files_nn(n) == 'none' .or. p%Ops(n) == "") cycle
    call opr%init(p%Ops(n),ms,.false.)
    call opr%set(ms, p%files_nn(n), 'none', &
        & [p%emax_nn, p%e2max_nn,p%lmax_nn])
    call HF%TransformToHF(ms,opr)
    call w%writef(p,ms,opr)
    call opr%fin()
  end do

  ! -- Ops from NN+3N file (srg evolved or two-body current) --
  do n = 1, size(p%files_3n)
    if(p%files_3n(n) == 'none' .or. p%Ops(n) == "") cycle
    call opr%init(p%Ops(n),ms,.true.)
    call opr%set(ms, p%files_nn(n), p%files_3n(n), &
        & [p%emax_nn, p%e2max_nn,p%lmax_nn], &
        & [p%emax_3n,p%e2max_3n,p%e3max_3n,p%lmax_3n])
    call HF%TransformToHF(ms,opr)
    call w%writef(p,ms,opr)
    call opr%fin()
  end do

  call HF%fin()
  call ms%fin()

  call fin_dbinomial_triangle()
  call timer%fin()
end program HFMain
