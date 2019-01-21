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
  implicit none
  type(InputParameters) :: p
  type(MSpace) :: ms
  type(Op) :: h
  type(HFSolver) :: HF
  type(MBPTEnergy) :: PT
  character(256) :: inputfile
  real(8) :: ti

  call timer%init()
  ti = omp_get_wtime()
  call init_dbinomial_triangle()
  call timer%Add("Init for Rotation Group",omp_get_wtime()-ti)

  call getarg(1,inputfile)
  call p%init(inputfile)
  call p%PrintInputParameters()

  if(p%e3max == 0) then
    call ms%init(p%Nucl, p%hw, p%emax, p%e2max)
  elseif(p%e3max > 0) then
    call ms%init(p%Nucl, p%hw, p%emax, p%e2max, e3max=p%e3max)
  end if

  call h%init('hamil',ms,.true.)
  call h%set(ms,p%files_nn(1),p%files_3n(1),&
      & [p%emax_nn,p%e2max_nn,p%lmax_nn],&
      & [p%emax_3n,p%e2max_3n,p%e3max_3n,p%lmax_3n])

  call HF%init(ms,h,alpha=p%alpha)
  call HF%solve(ms%sps,ms%one)
  call HF%TransformToHF(ms,h)

  call HF%fin()

  call PT%calc(ms,H)

  call h%fin()
  call ms%fin()

  call fin_dbinomial_triangle()
  call timer%fin()
end program HFMain
