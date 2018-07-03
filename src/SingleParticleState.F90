module SingleParticleModule
  use InputParameters
  implicit none
  type :: sps_array
    ! single particle state
    integer :: total_orbits
    integer, allocatable :: nn(:), ll(:), jj(:), itz(:), nshell(:)
    integer, allocatable :: num_orbit(:,:,:,:)
    integer, allocatable :: h_orbts(:)
    real(8), allocatable :: nocoef(:)
  contains
    procedure :: SetSPS
    procedure :: SetSPSIsospin
    procedure :: ReleaseSPS
    procedure :: ReleaseSPSIsospin
    procedure :: GetReferenceState
    procedure :: GetNOCoef
  end type sps_array
contains
  subroutine SetSPS(this, params)
    class(sps_array), intent(inout) :: this
      type(parameters), intent(in) :: params
      integer :: emax, n, nl, nn, ll, is, j, itz
      emax = params%emax
      allocate(this%num_orbit(0:emax/2, 0:emax, 1:2*emax+1, -1:1))
      this%num_orbit = 0
      this%total_orbits = (emax + 1) * (emax + 2)
      n = this%total_orbits
      allocate(this%nn(n))
      allocate(this%ll(n))
      allocate(this%jj(n))
      allocate(this%itz(n))
      allocate(this%nshell(n))
      allocate(this%h_orbts(n))
      allocate(this%nocoef(n))
      this%nn = 0; this%ll = 0; this%jj = 0
      this%itz = 0; this%nshell = 0
      this%h_orbts = 0
      this%nocoef = 0.d0
      n = 0
#ifdef debug
      write(*,'(a)') '##################################################'
      write(*,'(a)') ' Single-Particle State (PN formalism)'
      write(*,'(a)') '##################################################'
      write(*,'(a)') '      i,     n,     l,     j,    iz, shell'
#endif
      do nl = 0, emax
        do ll = 0, nl
          if(mod(nl - ll, 2) == 1) cycle
          nn = (nl - ll) / 2
          do is = -1, 1, 2
            j = 2 * ll + is
            if(j < 0) cycle
            do itz = -1, 1, 2
              n = n + 1
              this%nn(n) = nn
              this%ll(n) = ll
              this%jj(n) = j
              this%itz(n) = itz
              this%nshell(n) = 2 * nn + ll
              this%num_orbit(nn, ll, j, itz) = n
#ifdef debug
              write(*,'(6i7)') n, nn, ll, j, itz, 2 * nn + ll
#endif
            end do
          end do
        end do
      end do
#ifdef debug
      write(*,*)
#endif
      call this%GetReferenceState()
      call this%GetNOCoef()
    end subroutine SetSPS

  subroutine SetSPSIsospin(this, params)
    class(sps_array), intent(inout) :: this
      type(parameters), intent(in) :: params
      integer :: emax, n, nl, nn, ll, is, j
      emax = params%emax_3nf
      allocate(this%num_orbit(0:emax/2, 0:emax, 1:2*emax+1, 1))
      this%num_orbit = 0
      this%total_orbits = (emax + 1) * (emax + 2) / 2
      n = this%total_orbits
      allocate(this%nn(n))
      allocate(this%ll(n))
      allocate(this%jj(n))
      allocate(this%nshell(n))
      this%nn = 0; this%ll = 0; this%jj = 0
      this%nshell = 0
#ifdef debug
      write(*,'(a)') '##################################################'
      write(*,'(a)') ' Single-Particle State (Isospin formalism)'
      write(*,'(a)') '##################################################'
      write(*,'(a)') '      i,     n,     l,     j, shell'
#endif
      n = 0
      do nl = 0, emax
        do ll = 0, nl
          if(mod(nl - ll, 2) == 1) cycle
          nn = (nl - ll) / 2
          do is = -1, 1, 2
            j = 2 * ll + is
            if(j < 0) cycle
            n = n + 1
            this%nn(n) = nn
            this%ll(n) = ll
            this%jj(n) = j
            this%nshell(n) = 2 * nn + ll
            this%num_orbit(nn, ll, j, 1) = n
#ifdef debug
            write(*,'(5i7)') n, nn, ll, j, 2 * nn + ll
#endif
          end do
        end do
      end do
#ifdef debug
      write(*,*)
#endif
  end subroutine SetSPSIsospin

  subroutine GetReferenceState(this)
    class(sps_array), intent(inout) :: this
      character(256) :: ref
      integer :: num, i
      integer :: n, l, j, itz
      integer :: iunit = 40
      call getarg(2, ref)
      open(iunit, file = ref, status = 'old')
      read(iunit,*) num
      if(num < 1) return
      call skip_comment(iunit)
      do i = 1, num
        read(iunit,*) n, l, j, itz
        this%h_orbts(this%num_orbit(n,l,j,itz)) = 1
      end do
      if(myrank == 0) then
        write(*,'(2x, a)') 'Core Orbits: '
        do i = 1, this%total_orbits
          if(this%h_orbts(i) /= 0) then
            write(*,'(2x,a,i3,a,i3,a,i3,a,i3)') 'n = ', this%nn(i), &
                & ',  l = ', this%ll(i), ',  j = ', this%jj(i), &
                & ',  itz = ', this%itz(i)
          end if
        end do
        write(*,*)
      end if
  end subroutine GetReferenceState

  subroutine GetNOCoef(this)
    class(sps_array), intent(inout) :: this
      character(256) :: no
      integer :: num, i
      integer :: n, l, j, itz
      integer :: iunit = 40
      real(8) :: f
      call getarg(3, no)
      open(iunit, file = no, status = 'old')
      read(iunit,*) num
      if(num < 1) return
      call skip_comment(iunit)
      do i = 1, num
        read(iunit,*) n, l, j, itz, f
        this%nocoef(this%num_orbit(n,l,j,itz)) = dble(f)
      end do
      if(myrank == 0) then
        write(*,'(2x, a)') 'Coefficients for Normal Ordering: '
        do i = 1, this%total_orbits
          if(this%nocoef(i) /= 0) then
            write(*,'(2x,a,i3,a,i3,a,i3,a,i3,a,f7.3)') 'n = ', this%nn(i), &
                & ',  l = ', this%ll(i), ',  j = ', this%jj(i), &
                & ',  itz = ', this%itz(i), ',  Coef = ', this%nocoef(i)
          end if
        end do
        write(*,*)
      end if
  end subroutine GetNOCoef

  subroutine ReleaseSPS(this)
    class(sps_array), intent(inout) :: this
      deallocate(this%nn)
      deallocate(this%ll)
      deallocate(this%jj)
      deallocate(this%itz)
      deallocate(this%nshell)
      deallocate(this%num_orbit)
      deallocate(this%h_orbts)
      deallocate(this%nocoef)
  end subroutine ReleaseSPS

  subroutine ReleaseSPSIsospin(this)
    class(sps_array), intent(inout) :: this
      deallocate(this%nn)
      deallocate(this%ll)
      deallocate(this%jj)
      deallocate(this%nshell)
      deallocate(this%num_orbit)
  end subroutine ReleaseSPSIsospin

  subroutine skip_comment(nfile)
      implicit none
      integer,intent(in)::nfile
      character(1),parameter::com1='!', com2='#'
      character(1)::c1,c2
      read(nfile,'(2a1)') c1, c2
      do while  (c1 == com1 .or. c1 == com2 .or. c2 == com1 .or. c2 == com2)
        read(nfile,'(2a1)') c1, c2
      end do
      backspace(nfile)
  end subroutine skip_comment

end module SingleParticleModule
