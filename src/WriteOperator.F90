! Write operators in snt (Tokyo) format
module WriteOperator
  use omp_lib
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
  end subroutine InitWriteFiles

  subroutine SetFileName(this, ms, opr, filename)
    use ClassSys, only: sys
    class(WriteFiles), intent(inout) :: this
    type(MSpace), intent(in) :: ms
    type(Op), intent(in) :: opr
    character(*), intent(in), optional :: filename
    type(sys) :: s

    if(present(filename)) then
      this%filename = filename
      return
    end if

    this%filename = trim(opr%optr) // '_' // trim(ms%Nucl) // &
        & '_HF_hw' // trim(s%str(ms%hw)) // &
        & '_e' // trim(s%str(this%emax)) // &
        & '_2e' // trim(s%str(this%e2max)) // '.txt.snt'
  end subroutine SetFileName

  subroutine WriteFile(this,p,ms,opr,filename)
    use Profiler, only: timer
    use HFInput
    class(WriteFiles), intent(inout) :: this
    type(InputParameters), intent(in) :: p
    type(MSpace), intent(in) :: ms
    type(Op), intent(in) :: opr
    character(*), intent(in), optional :: filename
    type(SingleParticleOrbit) :: o
    integer :: wunit = 33
    integer :: i, cnt
    integer :: chbra, chket, bra, ket, max_ket
    integer :: a, b, c, d, Jbra, Jket
    real(8) :: ti

    ti = omp_get_wtime()
    call this%SetFileName(ms,opr,filename)

    open(wunit, file=this%filename, action = 'write')
    call p%PrintInputParameters(wunit)
    write(wunit,'(a,f18.8)') "# zero-body term: ", opr%zero
    write(wunit,'(a)') "# model space"
    write(wunit,'(4i4)') (this%emax+1)*(this%emax+2)/2,&
        & (this%emax+1)*(this%emax+2)/2, 0, 0
    cnt = 0
    do i = 1, ms%sps%norbs
      if(ms%sps%orb(i)%e > this%emax) cycle
      cnt = cnt + 1
      !o = ms%sps%get(i)
      o = ms%sps%orb(i)
      write(wunit,'(5i4)') cnt, o%n, o%l, o%j, o%z
    end do
    write(wunit,'(2a)') "# ", trim(opr%optr)

    cnt = 0
    do chbra = 1, ms%one%NChan
      do chket = 1, chbra
        if(.not. opr%one%MatCh(chbra,chket)%is) cycle
        do bra = 1, ms%one%jpz(chbra)%nst
          max_ket = ms%one%jpz(chket)%nst
          if(chbra == chket) max_ket = bra
          do ket = 1, max_ket
            if(abs(opr%one%MatCh(chbra,chket)%m(bra,ket)) < 1.d-8) cycle
            a = ms%one%jpz(chbra)%n2spi(bra)
            b = ms%one%jpz(chket)%n2spi(ket)
            if(ms%sps%orb(a)%e>this%emax) cycle
            if(ms%sps%orb(b)%e>this%emax) cycle
            cnt = cnt + 1
          end do
        end do
      end do
    end do

    write(wunit,'(i6,i2,f6.2)') cnt, 0, ms%hw
    do chbra = 1, ms%one%NChan
      do chket = 1, chbra
        if(.not. opr%one%MatCh(chbra,chket)%is) cycle
        do bra = 1, ms%one%jpz(chbra)%nst
          max_ket = ms%one%jpz(chket)%nst
          if(chbra == chket) max_ket = bra
          do ket = 1, max_ket
            if(abs(opr%one%MatCh(chbra,chket)%m(bra,ket)) < 1.d-8) cycle
            a = ms%one%jpz(chbra)%n2spi(bra)
            b = ms%one%jpz(chket)%n2spi(ket)
            if(ms%sps%orb(a)%e>this%emax) cycle
            if(ms%sps%orb(b)%e>this%emax) cycle
            write(wunit,'(2i4,f18.8)') a,b,opr%one%MatCh(chbra,chket)%m(bra,ket)
          end do
        end do
      end do
    end do

    cnt = 0
    do chbra = 1, ms%two%NChan
      do chket = 1, chbra
        if(.not. opr%two%MatCh(chbra,chket)%is) cycle
        do bra = 1, ms%two%jpz(chbra)%nst
          max_ket = ms%two%jpz(chket)%nst
          if(chbra == chket) max_ket = bra
          do ket = 1, max_ket
            if(abs(opr%two%MatCh(chbra,chket)%m(bra,ket)) < 1.d-8) cycle
            cnt = cnt + 1
          end do
        end do
      end do
    end do

    write(wunit,'(i12,i2,f6.2)') cnt, 0, ms%hw
    do chbra = 1, ms%two%NChan
      Jbra = ms%two%jpz(chbra)%j
      do chket = 1, chbra
        Jket = ms%two%jpz(chket)%j
        if(.not. opr%two%MatCh(chbra,chket)%is) cycle
        do bra = 1, ms%two%jpz(chbra)%nst
          max_ket = ms%two%jpz(chket)%nst
          if(chbra == chket) max_ket = bra
          do ket = 1, max_ket
            if(abs(opr%two%MatCh(chbra,chket)%m(bra,ket)) < 1.d-8) cycle
            a = ms%two%jpz(chbra)%n2spi1(bra)
            b = ms%two%jpz(chbra)%n2spi2(bra)
            c = ms%two%jpz(chket)%n2spi1(ket)
            d = ms%two%jpz(chket)%n2spi2(ket)
            if(ms%sps%orb(a)%e>this%emax) cycle
            if(ms%sps%orb(b)%e>this%emax) cycle
            if(ms%sps%orb(c)%e>this%emax) cycle
            if(ms%sps%orb(d)%e>this%emax) cycle
            if(ms%sps%orb(a)%e+ms%sps%orb(b)%e>this%e2max) cycle
            if(ms%sps%orb(c)%e+ms%sps%orb(d)%e>this%e2max) cycle
            if(opr%Scalar) write(wunit,'(5i4,f18.8)') &
                & a,b,c,d,Jbra,opr%two%MatCh(chbra,chket)%m(bra,ket)
            if(.not. opr%Scalar) write(wunit,'(6i4,f18.8)') &
                & a,b,c,d,Jbra,Jket,opr%two%MatCh(chbra,chket)%m(bra,ket)
          end do
        end do
      end do
    end do

    close(wunit)
    call timer%Add("Write to file", omp_get_wtime()-ti)
  end subroutine WriteFile

end module WriteOperator
