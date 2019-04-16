! Write operators in snt (Tokyo) format
module WriteOperator
  use omp_lib
  use HFInput
  use Operators
  implicit none

  public :: WriteFiles
  private :: InitWriteFiles
  private :: WriteFile
  private :: SetFileName

  type :: WriteFiles
    integer :: emax
    integer :: e2max
    character(256) :: filename = 'none'
  contains
    procedure :: init => InitWriteFiles
    procedure :: writef => WriteFile
    procedure :: SetFileName
  end type WriteFiles
contains

  subroutine InitWriteFiles(this, emax, e2max)
    class(WriteFiles), intent(inout) :: this
    integer, intent(in) :: emax, e2max
    this%emax = emax
    this%e2max = e2max
    this%filename = "op.out"
  end subroutine InitWriteFiles

  subroutine SetFileName(this, outdir, file_format, op)
    use ClassSys, only: sys
    class(WriteFiles), intent(inout) :: this
    character(*), intent(in) :: outdir
    character(*), intent(in) :: file_format
    type(Ops), intent(in) :: op
    type(MSpace), pointer :: ms
    type(sys) :: s

    ms => op%ms
    this%filename = trim(outdir) // "/" // trim(op%oprtr) // '_' // trim(ms%Nucl) // &
        & '_HF_hw' // trim(s%str(ms%hw)) // &
        & '_e' // trim(s%str(this%emax)) // &
        & '_2e' // trim(s%str(this%e2max)) // '.' // trim(file_format)
  end subroutine SetFileName

  subroutine WriteFile(this,p,op)
    use ClassSys, only: sys
    use Profiler, only: timer
    class(WriteFiles), intent(inout) :: this
    type(InputParameters), intent(in) :: p
    type(Ops), intent(in) :: op
    real(8) :: ti
    type(sys) :: s

    ti = omp_get_wtime()

    if(s%find(this%filename, "kshell.snt")) then
      call write_operator_kshell_snt(this%filename, p, op)
      return
    end if

    if(s%find(this%filename, "snt")) then
      call write_operator_ascii_snt(this%filename, p, op)
      return
    end if

    call timer%Add("Write to file", omp_get_wtime()-ti)
  end subroutine WriteFile

  subroutine write_operator_ascii_snt(f, p, op)
    character(*), intent(in) :: f
    type(InputParameters), intent(in) :: p
    type(Ops), intent(in) :: op
    type(SingleParticleOrbit) :: o
    type(MSpace), pointer :: ms
    integer :: wunit = 33
    integer :: i, cnt
    integer :: chbra, chket, bra, ket, max_ket
    integer :: a, b, c, d, Jbra, Jket

    ms => op%ms
    open(wunit, file=f, action = 'write')
    call p%PrintInputParameters(wunit)
    write(wunit,'(a,f18.8)') "# zero-body term: ", op%zero
    write(wunit,'(a)') "# model space"
    write(wunit,'(4i4)') (ms%emax+1)*(ms%emax+2)/2,&
        & (ms%emax+1)*(ms%emax+2)/2, 0, 0
    cnt = 0
    do i = 1, ms%sps%norbs
      if(ms%sps%orb(i)%e > ms%emax) cycle
      cnt = cnt + 1
      o = ms%sps%orb(i)
      write(wunit,'(5i4)') cnt, o%n, o%l, o%j, o%z
    end do
    write(wunit,'(2a)') "# ", trim(op%oprtr)

    cnt = 0
    do chbra = 1, ms%one%NChan
      do chket = 1, chbra
        if(.not. op%one%MatCh(chbra,chket)%is) cycle
        do bra = 1, ms%one%jpz(chbra)%n_state
          max_ket = ms%one%jpz(chket)%n_state
          if(chbra == chket) max_ket = bra
          do ket = 1, max_ket
            if(abs(op%one%MatCh(chbra,chket)%m(bra,ket)) < 1.d-8) cycle
            a = ms%one%jpz(chbra)%n2spi(bra)
            b = ms%one%jpz(chket)%n2spi(ket)
            if(ms%sps%orb(a)%e>ms%emax) cycle
            if(ms%sps%orb(b)%e>ms%emax) cycle
            cnt = cnt + 1
          end do
        end do
      end do
    end do

    write(wunit,'(i6,i2)') cnt, 0
    do chbra = 1, ms%one%NChan
      do chket = 1, chbra
        if(.not. op%one%MatCh(chbra,chket)%is) cycle
        do bra = 1, ms%one%jpz(chbra)%n_state
          max_ket = ms%one%jpz(chket)%n_state
          if(chbra == chket) max_ket = bra
          do ket = 1, max_ket
            if(abs(op%one%MatCh(chbra,chket)%m(bra,ket)) < 1.d-8) cycle
            a = ms%one%jpz(chbra)%n2spi(bra)
            b = ms%one%jpz(chket)%n2spi(ket)
            if(ms%sps%orb(a)%e>ms%emax) cycle
            if(ms%sps%orb(b)%e>ms%emax) cycle
            write(wunit,'(2i4,f18.8)') a,b,op%one%MatCh(chbra,chket)%m(bra,ket)
          end do
        end do
      end do
    end do

    cnt = 0
    do chbra = 1, ms%two%NChan
      do chket = 1, chbra
        if(.not. op%two%MatCh(chbra,chket)%is) cycle
        do bra = 1, ms%two%jpz(chbra)%n_state
          max_ket = ms%two%jpz(chket)%n_state
          if(chbra == chket) max_ket = bra
          do ket = 1, max_ket
            if(abs(op%two%MatCh(chbra,chket)%m(bra,ket)) < 1.d-8) cycle
            cnt = cnt + 1
          end do
        end do
      end do
    end do

    write(wunit,'(i12,i2 )') cnt,0
    do chbra = 1, ms%two%NChan
      Jbra = ms%two%jpz(chbra)%j
      do chket = 1, chbra
        Jket = ms%two%jpz(chket)%j
        if(.not. op%two%MatCh(chbra,chket)%is) cycle
        do bra = 1, ms%two%jpz(chbra)%n_state
          max_ket = ms%two%jpz(chket)%n_state
          if(chbra == chket) max_ket = bra
          do ket = 1, max_ket
            if(abs(op%two%MatCh(chbra,chket)%m(bra,ket)) < 1.d-8) cycle
            a = ms%two%jpz(chbra)%n2spi1(bra)
            b = ms%two%jpz(chbra)%n2spi2(bra)
            c = ms%two%jpz(chket)%n2spi1(ket)
            d = ms%two%jpz(chket)%n2spi2(ket)
            if(ms%sps%orb(a)%e>ms%emax) cycle
            if(ms%sps%orb(b)%e>ms%emax) cycle
            if(ms%sps%orb(c)%e>ms%emax) cycle
            if(ms%sps%orb(d)%e>ms%emax) cycle
            if(ms%sps%orb(a)%e+ms%sps%orb(b)%e>ms%e2max) cycle
            if(ms%sps%orb(c)%e+ms%sps%orb(d)%e>ms%e2max) cycle
            if(op%Scalar) write(wunit,'(5i4,f18.8)') &
                & a,b,c,d,Jbra,op%two%MatCh(chbra,chket)%m(bra,ket)
            if(.not. op%Scalar) write(wunit,'(6i4,f18.8)') &
                & a,b,c,d,Jbra,Jket,op%two%MatCh(chbra,chket)%m(bra,ket)
          end do
        end do
      end do
    end do

    close(wunit)
  end subroutine write_operator_ascii_snt

  subroutine write_operator_kshell_snt(f, p, op)
    character(*), intent(in) :: f
    type(InputParameters), intent(in) :: p
    type(Ops), intent(in) :: op
    type(SingleParticleOrbit) :: o
    type(SingleParticleOrbit), pointer :: oa, ob, oc, od
    type(MSpace), pointer :: ms
    type(Orbits) :: ksps
    integer :: wunit = 33
    integer :: i, cnt
    integer :: chbra, chket, bra, ket, max_ket
    integer :: a, b, c, d, Jbra, Jket
    integer :: ak, bk, ck, dk

    ms => op%ms
    call ksps%init(ms%emax, ms%lmax, "kshell")
    open(wunit, file=f, action = 'write')
    call p%PrintInputParameters(wunit)
    write(wunit,'(a,f18.8)') "# zero-body term: ", op%zero
    write(wunit,'(a)') "# model space"
    write(wunit,'(4i4)') (ms%emax+1)*(ms%emax+2)/2,&
        & (ms%emax+1)*(ms%emax+2)/2, 0, 0
    cnt = 0
    do i = 1, ksps%norbs
      if(ksps%orb(i)%e > ms%emax) cycle
      cnt = cnt + 1
      o = ksps%orb(i)
      write(wunit,'(5i4)') cnt, o%n, o%l, o%j, o%z
    end do
    write(wunit,'(2a)') "# ", trim(op%oprtr)

    cnt = 0
    do chbra = 1, ms%one%NChan
      do chket = 1, chbra
        if(.not. op%one%MatCh(chbra,chket)%is) cycle
        do bra = 1, ms%one%jpz(chbra)%n_state
          max_ket = ms%one%jpz(chket)%n_state
          if(chbra == chket) max_ket = bra
          do ket = 1, max_ket
            if(abs(op%one%MatCh(chbra,chket)%m(bra,ket)) < 1.d-8) cycle
            a = ms%one%jpz(chbra)%n2spi(bra)
            b = ms%one%jpz(chket)%n2spi(ket)
            if(ms%sps%orb(a)%e>ms%emax) cycle
            if(ms%sps%orb(b)%e>ms%emax) cycle
            cnt = cnt + 1
          end do
        end do
      end do
    end do

    write(wunit,'(i6,i2)') cnt, 0
    do chbra = 1, ms%one%NChan
      do chket = 1, chbra
        if(.not. op%one%MatCh(chbra,chket)%is) cycle
        do bra = 1, ms%one%jpz(chbra)%n_state
          max_ket = ms%one%jpz(chket)%n_state
          if(chbra == chket) max_ket = bra
          do ket = 1, max_ket
            if(abs(op%one%MatCh(chbra,chket)%m(bra,ket)) < 1.d-8) cycle
            a = ms%one%jpz(chbra)%n2spi(bra)
            b = ms%one%jpz(chket)%n2spi(ket)
            if(ms%sps%orb(a)%e>ms%emax) cycle
            if(ms%sps%orb(b)%e>ms%emax) cycle
            oa => ms%sps%GetOrbit(a)
            ob => ms%sps%GetOrbit(b)
            ak = ksps%nljz2idx(oa%n, oa%l, oa%j, oa%z)
            bk = ksps%nljz2idx(ob%n, ob%l, ob%j, ob%z)
            write(wunit,'(2i4,f18.8)') ak,bk,op%one%MatCh(chbra,chket)%m(bra,ket)
          end do
        end do
      end do
    end do

    cnt = 0
    do chbra = 1, ms%two%NChan
      do chket = 1, chbra
        if(.not. op%two%MatCh(chbra,chket)%is) cycle
        do bra = 1, ms%two%jpz(chbra)%n_state
          max_ket = ms%two%jpz(chket)%n_state
          if(chbra == chket) max_ket = bra
          do ket = 1, max_ket
            if(abs(op%two%MatCh(chbra,chket)%m(bra,ket)) < 1.d-8) cycle
            cnt = cnt + 1
          end do
        end do
      end do
    end do

    write(wunit,'(i12,i2)') cnt, 0
    do chbra = 1, ms%two%NChan
      Jbra = ms%two%jpz(chbra)%j
      do chket = 1, chbra
        Jket = ms%two%jpz(chket)%j
        if(.not. op%two%MatCh(chbra,chket)%is) cycle
        do bra = 1, ms%two%jpz(chbra)%n_state
          max_ket = ms%two%jpz(chket)%n_state
          if(chbra == chket) max_ket = bra
          do ket = 1, max_ket
            if(abs(op%two%MatCh(chbra,chket)%m(bra,ket)) < 1.d-8) cycle
            a = ms%two%jpz(chbra)%n2spi1(bra)
            b = ms%two%jpz(chbra)%n2spi2(bra)
            c = ms%two%jpz(chket)%n2spi1(ket)
            d = ms%two%jpz(chket)%n2spi2(ket)
            if(ms%sps%orb(a)%e>ms%emax) cycle
            if(ms%sps%orb(b)%e>ms%emax) cycle
            if(ms%sps%orb(c)%e>ms%emax) cycle
            if(ms%sps%orb(d)%e>ms%emax) cycle
            if(ms%sps%orb(a)%e+ms%sps%orb(b)%e>ms%e2max) cycle
            if(ms%sps%orb(c)%e+ms%sps%orb(d)%e>ms%e2max) cycle
            oa => ms%sps%GetOrbit(a)
            ob => ms%sps%GetOrbit(b)
            oc => ms%sps%GetOrbit(c)
            od => ms%sps%GetOrbit(d)
            ak = ksps%nljz2idx(oa%n, oa%l, oa%j, oa%z)
            bk = ksps%nljz2idx(ob%n, ob%l, ob%j, ob%z)
            ck = ksps%nljz2idx(oc%n, oc%l, oc%j, oc%z)
            dk = ksps%nljz2idx(od%n, od%l, od%j, od%z)
            if(op%Scalar) write(wunit,'(5i4,f18.8)') &
                & ak,bk,ck,dk,Jbra,op%two%MatCh(chbra,chket)%m(bra,ket)
            if(.not. op%Scalar) write(wunit,'(6i4,f18.8)') &
                & ak,bk,ck,dk,Jbra,Jket,op%two%MatCh(chbra,chket)%m(bra,ket)
          end do
        end do
      end do
    end do

    close(wunit)
  end subroutine write_operator_kshell_snt

end module WriteOperator
