module Atomic
  use HFInput
  use ModelSpace
  use Operators
  use HartreeFock
  use HFMBPT
  implicit none
contains
  subroutine atomic_case(inputfile, conffile)
    character(*), intent(in) :: inputfile, conffile
    type(InputParameters) :: p
    type(MSpace) :: ms
    type(Ops) :: h
    type(HFSolver) :: HF
    type(MBPTEnergy) :: PT
    integer :: wunit = 15

    call p%init(inputfile)
    call read_hamil_from_snt(ms, h,p%int_nn_file, conffile)
    call HF%init(h,alpha=p%alpha)
    call HF%solve()

    if(.not. p%is_MBPTEnergy) then
      open(wunit, file = p%summary_file, action='write',status='replace')
      call p%PrintInputParameters(wunit)
      call print_single_particle_energies(HF,wunit)
      write(wunit,'(a,6x,a,9x,a,9x,a,13x,a)') &
          & "# Operator", "HF energy", "2nd order", "3rd order", "Total"
      write(wunit,'(a,4f18.8)') 'hamil: ', HF%ehf, 0.d0, 0.d0, HF%ehf
      close(wunit)
    end if

    if(p%is_MBPTEnergy) then
      h = HF%BasisTransform(h)
      call PT%calc(H)
      open(wunit, file = p%summary_file, action='write',status='replace')
      call p%PrintInputParameters(wunit)
      call print_single_particle_energies(HF,wunit)
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

    call HF%fin()
    call h%fin()
    call ms%fin()
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
    real(8) :: me, kin, pot, ti

    ti = omp_get_wtime()
    h%oprtr = "hamil"
    h%jr = 0
    h%pr = 1
    h%zr = 0
    ms%is_three_body_jt = .false.
    ms%is_three_body = .false.

    open(runit, file=f, action='read',iostat=io)
    if(io /= 0) then
      write(*,*)
      write(*,'(2a)') 'File open error: ', trim(f)
      write(*,*)
      return
    end if
    call skip_comment(runit,'#')
    read(runit,*) porbs, norbs, pc, nc
    call skip_comment(runit,'#')
    lines = porbs+norbs
    ms%sps%norbs = lines
    allocate(ms%sps%orb(lines))
    do line = 1, lines
      read(runit,*) idx, n, l, j, z, e
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
      read(runit,*) a, b, kin, pot ! one-body part
      call h%one%SetOBME(a, b, kin-dble(N_e)*pot)
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
    call timer%Add("Model space & Hamiltonian", omp_get_wtime()-ti)
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
        call o%SetCoreValenceOutside(2)
      elseif(abs(1.d0 - ms%NOCoef(idx)) < 1.d-6) then
        call o%SetCoreValenceOutside(0)
      else
        call o%SetCoreValenceOutside(1)
      end if
    end do
    write(*,"(a,i4,2a)") "# number of electrons ", N_e, " from ", trim(f)
  end subroutine GetElectronConfFromFile

  subroutine print_single_particle_energies(HF, iunit)
    type(HFSolver), intent(in) :: HF
    integer, intent(in) :: iunit
    type(MSpace), pointer :: ms
    type(OneBodyPart) :: F_HF
    type(SingleParticleOrbit), pointer :: oi
    character(40) :: char_spe, sp_label, ph
    character(256) :: line
    integer :: ch, i

    ms => HF%ms
    F_HF = HF%F
    do ch = 1, ms%one%NChan
      F_HF%MatCh(ch,ch)%DMat = HF%C%MatCh(ch,ch)%DMat%T() * &
          &  HF%F%MatCh(ch,ch)%DMat * HF%C%MatCh(ch,ch)%DMat
    end do
    write(iunit,'(a)') "#  Hartree-Fock single-particle energies (spe)"
    write(iunit,'(a)') "#    Orbit:                        spe     Occ"
    do i = 1, ms%sps%norbs
      oi => ms%sps%GetOrbit(i)
      sp_label = trim(ms%sps%GetLabelFromIndex(i))
      ph = "hole"
      if(oi%GetHoleParticle() == 1) ph = "particle"
      write(line,"(a,a8,a,a12,f14.6,f8.4)") "# ", trim(sp_label(2:)), ": ", &
          & trim(ph), F_HF%GetOBME(i,i), &
          & HF%Occ%GetOBME(i,i)
      write(iunit,"(a)") trim(line)
    end do
    call F_HF%fin()

  end subroutine print_single_particle_energies

end module Atomic
