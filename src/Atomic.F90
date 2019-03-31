module Atomic
  use HFInput
  use ModelSpace
  use Operators
  use HartreeFock
  use MBPT
  implicit none
contains
  subroutine atomic_case(inputfile, conffile)
    character(*), intent(in) :: inputfile, conffile
    type(InputParameters) :: p
    type(MSpace) :: ms
    type(Ops) :: h
    type(HFSolver) :: HF
    type(MBPTEnergy) :: PT
    integer :: wunit = 6

    call p%init(inputfile)
    call read_hamil_from_snt(ms, h,p%int_nn_file, conffile)
    call HF%init(h,alpha=p%alpha)
    call HF%solve()
    call HF%TransformToHF(h)
    if(p%is_MBPTEnergy) then
      call PT%calc(H)
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
  end subroutine atomic_case

  subroutine read_hamil_from_snt(ms, h, f, conffile)
    use MyLibrary, only: skip_comment
    type(Ops), intent(inout) :: h
    type(MSpace), intent(inout), target :: ms
    character(*), intent(in) :: f, conffile
    integer :: runit = 25
    integer :: a, b, c, d, J, lines, idx, n, l, z, e
    integer :: porbs, norbs, pc, nc, line
    integer :: io, i, N_e
    integer :: nmax, lmax, jmax
    integer :: nmin, lmin, jmin
    type(SingleParticleOrbit), pointer :: o
    real(8) :: me

    h%oprtr = "hamil"
    h%jr = 0
    h%pr = 1
    h%zr = 0
    ms%is_three_body_jt = .false.
    ms%is_three_body = .false.

    open(runit, file=f, action='read',iostat=io)
    if(io /= 0) then
      write(*,'(2a)') 'File open error: ', trim(f)
      return
    end if
    call skip_comment(runit,'#')
    read(runit,*) porbs, norbs, pc, nc
    lines = porbs+norbs
    ms%sps%norbs = lines
    allocate(ms%sps%orb(lines))
    do line = 1, lines
      read(runit,*) idx, n, l, j, e, z
      call ms%sps%orb(line)%set(n,l,j,z,idx,e)
    end do

    nmax = -1
    lmax = -1
    jmax = -1
    nmin = 100000
    lmin = 100000
    jmin = 100000
    do i = 1, ms%sps%norbs
      ms%sps%emax = max(ms%sps%emax, ms%sps%orb(i)%e)
      nmax = max(nmax, ms%sps%orb(i)%n)
      lmax = max(lmax, ms%sps%orb(i)%l)
      jmax = max(jmax, ms%sps%orb(i)%j)
      nmin = min(nmin, ms%sps%orb(i)%n)
      lmin = min(lmin, ms%sps%orb(i)%l)
      jmin = min(jmin, ms%sps%orb(i)%j)
    end do
    ms%sps%lmax = lmax
    allocate(ms%sps%nljz2idx(nmin:nmax, lmin:lmax, jmin:jmax, -1:1))
    ms%sps%nljz2idx(:,:,:,:) = 0
    do line = 1, ms%sps%norbs
      o => ms%sps%GetOrbit(line)
      ms%sps%nljz2idx(o%n, o%l, o%j, o%z) = o%idx
    end do


    call GetElectronConfFromFile(ms, conffile, N_e)
    call ms%GetParticleHoleOrbits()
    ms%A = N_e
    ms%Z = N_e
    ms%N = 0
    call ms%one%init(ms%sps)
    call ms%two%init(ms%sps, 2*ms%sps%emax)

    h%ms => ms
    call h%one%init(ms%one, .true., "hamil", 0, 1, 0)
    call h%two%init(ms%two, .true., "hamil", 0, 1, 0)

    call skip_comment(runit,'#')
    read(runit,*) lines
    call skip_comment(runit,'#')
    do line = 1, lines
      read(runit,*) a, b, me ! one-body part
      call h%one%SetOBME(a,b,me)
    end do

    call skip_comment(runit,'#')
    read(runit,*) lines
    call skip_comment(runit,'#')
    do line = 1, lines
      read(runit,*) a, b, c, d, J, me
#ifdef TwoBodyOperatorDebug
      write(*,'(5i3,f12.6)') a, b, c, d, J, me
#endif
      call h%two%SetTwBME(a,b,c,d,J,me)
    end do
    close(runit)
  end subroutine read_hamil_from_snt

  subroutine GetElectronConfFromFile(ms, f, N_e)
    use MyLibrary, only: skip_comment
    type(MSpace), intent(inout) :: ms
    character(*), intent(in) :: f
    integer, intent(inout) :: N_e
    integer :: runit=20
    integer :: nn, ll, jj, idx, occ_num
    type(SingleParticleOrbit), pointer :: o

    N_e = 0
    allocate(ms%NOCoef(ms%sps%norbs))
    ms%NOCoef(:) = 0.d0
    open(runit, file=f, status='old', action='read')
    call skip_comment(runit, '!')
    do
      read(runit,*,end=999) nn, ll, jj, occ_num
      idx = ms%sps%nljz2idx(nn,ll,jj,-1)
      if(idx == 0) stop "configuration file error"
      ms%NOcoef(idx) = dble(occ_num)/dble(jj+1)
      N_e = N_e + occ_num
    end do
999 close(runit)

    do idx = 1, ms%sps%norbs
      o => ms%sps%GetOrbit(idx)
      call o%SetOccupation(ms%NOcoef(idx))
      if(ms%NOcoef(idx) < 1.d-6) then
        call o%SetHoleParticleValence(1)
      elseif(abs(1.d0 - ms%NOCoef(idx)) < 1.d-6) then
        call o%SetHoleParticleValence(0)
      else
        call o%SetHoleParticleValence(2)
      end if
    end do
  end subroutine GetElectronConfFromFile

end module Atomic
