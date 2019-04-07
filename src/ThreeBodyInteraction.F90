module ThreeBodyInteraction
  use omp_lib
  use Profiler, only: timer
  use SingleParticleState
  use ThreeBodyModelSpace
  implicit none

  public :: ThreeBodyForce

  private :: InitThreeBodyForce
  private :: FinThreeBodyForce
  private :: GetThBME_isospin
  private :: GetThBME_pn
  private :: ReadIsospinThreeBodyFile
  private :: CopyThreeBodyForce
  private :: SumThreeBodyForce
  private :: SubtractThreeBodyForce
  private :: ScaleThreeBodyForce
  private :: PrintThreeBodyForce
  private :: NormalOrderingFrom3To2

  private :: Set3BodyReadFile        ! setter
  private :: Set3BodyFileBoundaries  ! setter
  private :: ReadScalar3BFile
  private :: ReadTensor3BFile
  private :: read_scalar_me3j_ascii_txt
  private :: read_scalar_me3j_gzip
  private :: read_scalar_me3j_ascii
  private :: read_scalar_me3j_binary_stream
  private :: read_scalar_me3j_binary_comp
  private :: read_scalar_me3j_binary

  type :: ThreeBodyForceChannel
#ifdef single_precision_three_body_file
    real(4), allocatable :: v(:,:)
#else
    real(8), allocatable :: v(:,:)
#endif
  end type ThreeBodyForceChannel

  type :: ThreeBodyForce
    type(ThreeBodyForceChannel), allocatable :: MatCh(:,:)
    type(NonOrthIsospinThreeBodySpace), pointer :: thr
    logical :: Scalar
  contains
    procedure :: InitThreeBodyForce
    procedure :: FinThreeBodyForce
    procedure :: GetThBME_isospin
    procedure :: GetThBME_pn
    procedure :: CopyThreeBodyForce
    procedure :: SumThreeBodyForce
    procedure :: SubtractThreeBodyForce
    procedure :: ScaleThreeBodyForce
    procedure :: PrintThreeBodyForce
    procedure :: NormalOrderingFrom3To2

    generic :: init => InitThreeBodyForce
    generic :: fin => FinThreeBodyForce
    generic :: GetThBME => GetThBME_isospin, GetThBME_pn
    generic :: prt => PrintThreeBodyForce
    generic :: assignment(=) => CopyThreeBodyForce
    generic :: operator(+) => SumThreeBodyForce
    generic :: operator(-) => SubtractThreeBodyForce
    generic :: operator(*) => ScaleThreeBodyForce
  end type ThreeBodyForce

  type :: Read3BodyFiles
    character(:), allocatable :: file_3n
    integer :: emax3=-1, e2max3=-1, e3max3=-1, lmax3=-1
  contains
    procedure :: Set3BodyReadFile        ! setter
    procedure :: Set3BodyFileBoundaries  ! setter
    generic :: set => Set3BodyReadFile, Set3BodyFileBoundaries
    procedure :: ReadIsospinThreeBodyFile
    ! methods for three-body marix element
    procedure :: ReadScalar3BFile
    procedure :: ReadTensor3BFile
    procedure :: read_scalar_me3j_ascii_txt
    procedure :: read_scalar_me3j_gzip
    procedure :: read_scalar_me3j_ascii
    procedure :: read_scalar_me3j_binary_stream
    procedure :: read_scalar_me3j_binary_comp
    procedure :: read_scalar_me3j_binary
  end type Read3BodyFiles

contains
  subroutine InitThreeBodyForce(this, thr)
    class(ThreeBodyForce), intent(inout) :: this
    type(NonOrthIsospinThreeBodySpace), target, intent(in) :: thr
    integer :: ch, n_state
    ! so far only scalar
    this%Scalar = .true.
    this%thr => thr
    allocate(this%MatCh(thr%NChan,thr%NChan))

    do ch = 1, thr%NChan
      n_state = thr%jpt(ch)%n_state
      allocate(this%MatCh(ch,ch)%v(n_state,n_state))
#ifdef single_precision_three_body_file
      this%MatCh(ch,ch)%v(:,:) = 0.0
#else
      this%MatCh(ch,ch)%v(:,:) = 0.d0
#endif
    end do
  end subroutine InitThreeBodyForce

  subroutine FinThreeBodyForce(this)
    class(ThreeBodyForce), intent(inout) :: this
    integer :: ch
    do ch = 1, this%thr%NChan
      deallocate(this%MatCh(ch,ch)%v)
    end do
    deallocate(this%MatCh)
    this%thr => null()
  end subroutine FinThreeBodyForce

  subroutine CopyThreeBodyForce(a, b)
    class(ThreeBodyForce), intent(inout) :: a
    type(ThreeBodyForce), intent(in) :: b
    integer :: ch, ndim
    if(allocated(a%MatCH)) call a%fin()
    a%thr => b%thr
    allocate(a%MatCh(a%thr%NChan,a%thr%NChan))
    do ch = 1, a%thr%NChan
      ndim = a%thr%jpt(ch)%n_state
      allocate(a%MatCh(ch,ch)%v(ndim,ndim))
      a%MatCh(ch,ch)%v(:,:) = b%MatCh(ch,ch)%v(:,:)
    end do
  end subroutine CopyThreeBodyForce

  function SumThreeBodyForce(a, b) result(c)
    class(ThreeBodyForce), intent(in) :: a
    type(ThreeBodyForce), intent(in) :: b
    type(ThreeBodyForce) :: c
    integer :: ch
    c = a
    do ch = 1, a%thr%NChan
      c%MatCh(ch,ch)%v(:,:) = a%MatCh(ch,ch)%v(:,:) + b%MatCh(ch,ch)%v(:,:)
    end do
  end function SumThreeBodyForce

  function SubtractThreeBodyForce(a, b) result(c)
    class(ThreeBodyForce), intent(in) :: a
    type(ThreeBodyForce), intent(in) :: b
    type(ThreeBodyForce) :: c
    integer :: ch
    c = a
    do ch = 1, a%thr%NChan
      c%MatCh(ch,ch)%v(:,:) = a%MatCh(ch,ch)%v(:,:) - b%MatCh(ch,ch)%v(:,:)
    end do
  end function SubtractThreeBodyForce

  function ScaleThreeBodyForce(a, b) result(c)
    class(ThreeBodyForce), intent(in) :: a
    real(8), intent(in) :: b
    type(ThreeBodyForce) :: c
    integer :: ch
    c = a
    do ch = 1, a%thr%NChan
#ifdef single_precision_three_body_file
      c%MatCh(ch,ch)%v(:,:) = real(b) * a%MatCh(ch,ch)%v(:,:)
#else
      c%MatCh(ch,ch)%v(:,:) = b * a%MatCh(ch,ch)%v(:,:)
#endif
    end do
  end function ScaleThreeBodyForce

  function GetThBME_pn(this,i1,i2,i3,J12,&
        & i4,i5,i6,J45,J) result(r)
    use MyLibrary, only: dcg
    class(ThreeBodyForce), intent(in) :: this
    integer, intent(in) :: i1,i2,i3,i4,i5,i6
    integer, intent(in) :: J12,J45,J
    type(Orbits), pointer :: sps
    type(OrbitsIsospin), pointer :: isps
    type(SingleParticleOrbit), pointer :: o1,o2,o3,o4,o5,o6
    integer :: z1, z2, z3, z4, z5, z6, T12, T45, T
    integer :: a, b, c, d, e, f, Z
    real(8) :: r
    !real(8) :: ti ! --- test

    !ti = omp_get_wtime() ! --- test (conflict with omp)
    sps => this%thr%sps
    isps => this%thr%isps
    r = 0.d0
    if(i1 == i2 .and. mod(J12, 2) == 1) return
    if(i4 == i5 .and. mod(J45, 2) == 1) return
    o1 => sps%GetOrbit(i1)
    o2 => sps%GetOrbit(i2)
    o3 => sps%GetOrbit(i3)
    o4 => sps%GetOrbit(i4)
    o5 => sps%GetOrbit(i5)
    o6 => sps%GetOrbit(i6)
    if(o1%e + o2%e + o3%e > this%thr%e3max) return
    if(o4%e + o5%e + o6%e > this%thr%e3max) return
    z1 = o1%z; z2 = o2%z; z3 = o3%z
    z4 = o4%z; z5 = o5%z; z6 = o6%z
    if(z1+z2+z3 /= z4+z5+z6) return
    Z = z1 + z2 + z3

    a = isps%nlj2idx( o1%n, o1%l, o1%j )
    b = isps%nlj2idx( o2%n, o2%l, o2%j )
    c = isps%nlj2idx( o3%n, o3%l, o3%j )
    d = isps%nlj2idx( o4%n, o4%l, o4%j )
    e = isps%nlj2idx( o5%n, o5%l, o5%j )
    f = isps%nlj2idx( o6%n, o6%l, o6%j )
    do T12 = 0, 1
      if(abs(z1+z2) > 2*T12) cycle
      do T45 = 0, 1
        if(abs(z4+z5) > 2*T45) cycle
        do T = max(abs(2*T12-1),abs(2*T45-1)), min(2*T12+1,2*T45+1), 2
          if(abs(Z) > T) cycle
          r = r + &
              & this%GetThBME(a,b,c,J12,T12,&
              & d,e,f,J45,T45,J,T) * &
              & dcg(1,z1,1,z2,2*T12,z1+z2) * dcg(2*T12,z1+z2,1,z3,T,Z) * &
              & dcg(1,z4,1,z5,2*T45,z4+z5) * dcg(2*T45,z4+z5,1,z6,T,Z)
        end do
      end do
    end do
    !call timer%add("GetThBME_pn", omp_get_wtime()-ti) ! --- test (conflict with omp)
  end function GetThBME_pn

  function GetThBME_Isospin(this,i1,i2,i3,J12,T12,&
        & i4,i5,i6,J45,T45,J,T) result(r)
    class(ThreeBodyForce), intent(in) :: this
    integer, intent(in) :: i1,i2,i3,i4,i5,i6
    integer, intent(in) :: J12,T12,J45,T45,J,T
    type(OrbitsIsospin), pointer :: isps
    type(SingleParticleOrbitIsospin), pointer :: o1,o2,o3,o4,o5,o6
    type(NonOrthIsospinThreeBodySpace), pointer :: tbs
    type(coef), pointer :: bra_ch, ket_ch
    integer :: ch, idxbra, idxket, bra, ket
    integer :: P123, P456
    integer :: isorted_bra, isorted_ket
    integer :: ibra, iket
    real(8) :: r
    !real(8) :: ti ! --- test

    r = 0.d0
    !ti = omp_get_wtime() ! --- test (conflict with omp)
    if(i1 == i2 .and. mod(J12+T12,2) == 0) return
    if(i4 == i5 .and. mod(J45+T45,2) == 0) return
    isps => this%thr%isps
    o1 => isps%GetOrbit(i1)
    o2 => isps%GetOrbit(i2)
    o3 => isps%GetOrbit(i3)
    o4 => isps%GetOrbit(i4)
    o5 => isps%GetOrbit(i5)
    o6 => isps%GetOrbit(i6)
    tbs => this%thr
    P123 = (-1) ** (o1%l+o2%l+o3%l)
    P456 = (-1) ** (o4%l+o5%l+o6%l)
    if(P123 * P456 /= 1) then
      write(*,*) 'Warning: in GetThBMEIso_scalar: P'
      return
    end if
    ch = tbs%jpt2ch(J,P123,T)
    if(ch == 0) return
    idxbra = tbs%jpt(ch)%sorting(i1,i2,i3)
    idxket = tbs%jpt(ch)%sorting(i4,i5,i6)
    if(idxbra * idxket == 0) return
    isorted_bra = tbs%jpt(ch)%sort(idxbra)%idx_sorted
    isorted_ket = tbs%jpt(ch)%sort(idxket)%idx_sorted
    if(isorted_bra * isorted_ket == 0) return

    bra_ch => tbs%jpt(ch)%sort(idxbra)%JT(J12,T12)
    ket_ch => tbs%jpt(ch)%sort(idxket)%JT(J45,T45)

    do ibra = 1, bra_ch%n
      bra = bra_ch%idx2num(ibra)
      do iket = 1, ket_ch%n
        ket = ket_ch%idx2num(iket)
        r = r + dble(this%MatCh(ch,ch)%v(bra,ket) * &
            & bra_ch%TrnsCoef(ibra) * ket_ch%TrnsCoef(iket))
      end do
    end do
    !call timer%add("GetThBME_isospin", omp_get_wtime()-ti) ! --- test (conflict with omp)
  end function GetThBME_Isospin

  subroutine PrintThreeBodyForce(this, wunit)
    use ClassSys, only: sys
    use LinAlgLib
    class(ThreeBodyForce), intent(in) :: this
    integer, intent(in), optional :: wunit
    integer :: chbra, chket
    type(sys) :: s
    character(:), allocatable :: msg
    type(DMat) :: m

    do chbra = 1, this%thr%NChan
      do chket = 1, this%thr%NChan
        if(chbra /= chket) cycle
        msg = "hamil " // trim(s%str(chbra)) &
            &  // " " // trim(s%str(chket))
        call m%ini(this%thr%jpt(chbra)%n_state, this%thr%jpt(chket)%n_state)
        m%m(:,:) = this%MatCh(chbra,chket)%v(:,:)
        call m%prt(msg=msg,iunit=wunit)
        call m%fin()
      end do
    end do
  end subroutine PrintThreeBodyForce

  ! only Three-body interaction
  function NormalOrderingFrom3To2(this, two) result(r)
    use MyLibrary, only: triag, sjs
    use TwoBodyOperator
    class(ThreeBodyForce), intent(in) :: this
    type(TwoBodySpace), intent(in) :: two
    type(TwoBodyPart) :: r
    type(SingleParticleOrbit), pointer :: oh
    integer :: i1, i2, i3, i4, ih, J12, J34, jh, Jbra, Jket
    integer :: e1, e2, e3, e4, eh
    integer :: chbra, chket, bra, ket, bra_max, ket_max
    real(8) :: v, fact, vsum, tfact

    !call r%init(two,this%Scalar,this%oprtr,this%jr,this%pr,this%zr)
    call r%init(two,this%Scalar,"hamil",0,1,0)

    do chbra = 1, two%NChan
      do chket = 1, chbra
        if( .not. r%MatCh(chbra,chket)%is) cycle

        bra_max = two%jpz(chbra)%n_state
        ket_max = two%jpz(chket)%n_state

        J12 = two%jpz(chbra)%j
        J34 = two%jpz(chket)%j

        do bra = 1, bra_max
          i1 = two%jpz(chbra)%n2spi1(bra)
          i2 = two%jpz(chbra)%n2spi2(bra)
          e1 = two%sps%orb(i1)%e
          e2 = two%sps%orb(i2)%e
          if(chbra == chket) ket_max = bra
          do ket = 1, ket_max
            i3 = two%jpz(chket)%n2spi1(ket)
            i4 = two%jpz(chket)%n2spi2(ket)
            e3 = two%sps%orb(i3)%e
            e4 = two%sps%orb(i4)%e

            fact = 1.d0
            if(i1==i2) fact = fact / dsqrt(2.d0)
            if(i3==i4) fact = fact / dsqrt(2.d0)

            vsum = 0.d0
            do ih = 1, two%sps%norbs
              oh => two%sps%GetOrbit(ih)
              if(oh%occ < 1.d-6) cycle
              jh = oh%j
              eh = oh%e
              if(e1+e2+eh > this%thr%e3max) cycle
              if(e3+e4+eh > this%thr%e3max) cycle

              v = 0.d0
              do Jbra = abs(2*J12-jh), (2*J12+jh), 2
                do Jket = abs(2*J34-jh), (2*J34+jh), 2
                  if(triag(Jbra,Jket,0)) cycle
                  tfact = sqrt(dble( (Jbra+1) * (Jket+1) ))
                  if(.not. this%Scalar) then
                    tfact = sqrt( dble( (Jbra+1) * (Jket+1) ) ) * &
                        & (-1.d0) ** ( J12 + (jh+Jket)/2 + 0) * &
                        & sjs(2*J12,Jbra,jh,Jket,2*J34,0)
                  end if
                  v = v + tfact * oh%occ * &
                      & this%GetThBME(i1,i2,ih,J12,i3,i4,ih,J34,Jket)
                end do
              end do
            end do
            r%MatCh(chbra,chket)%m(bra,ket) = vsum * fact / sqrt( dble( (2*J12+1) * (2*J34+1) ) )
            r%MatCh(chket,chbra)%m(ket,bra) = r%MatCh(chbra,chket)%m(bra,ket) * &
                & (-1.d0) ** (J12-J34)
          end do
        end do

      end do
    end do

  end function NormalOrderingFrom3To2


  !
  !
  ! reading three-body scalar
  !
  !
  subroutine Set3BodyReadFile(this, file_3n)
    use ClassSys, only: sys
    class(Read3BodyFiles), intent(inout) :: this
    character(*), intent(in), optional :: file_3n
    type(sys) :: s
    logical :: ex

    ! -- three-body file
    this%file_3n = 'none'
    if(present(file_3n)) this%file_3n = file_3n
    select case(this%file_3n)
    case("NONE", "none", "None")
    case default
      ex = s%isfile(this%file_3n, "SetReadFiles: three-body file")
    end select
  end subroutine Set3BodyReadFile

  subroutine Set3BodyFileBoundaries(this, emax, e2max, e3max, lmax)
    class(Read3BodyFiles), intent(inout) :: this
    integer, intent(in) :: emax, e2max, e3max, lmax

    this%emax3 = emax
    this%e2max3= e2max
    this%e3max3= e3max
    this%lmax3 = lmax

  end subroutine Set3BodyFileBoundaries

  subroutine ReadIsospinThreeBodyFile(this, V)
    class(Read3BodyFiles), intent(in) :: this
    type(ThreeBodyForce), intent(inout) :: V
    real(8) :: ti

    ti = omp_get_wtime()

    select case(this%file_3n)
    case('None', 'NONE', 'none')
      write(*,*) "No three-body matrix element."
      return
    case default

      if(V%Scalar) call this%ReadScalar3BFile(V)
      if(.not. V%Scalar) call this%ReadTensor3BFile(V)
    end select

    call timer%Add('Read from file', omp_get_wtime()-ti)

  end subroutine ReadIsospinThreeBodyFile

  subroutine ReadScalar3BFile(this,thr)
    use ClassSys, only: sys
    class(Read3BodyFiles), intent(in) :: this
    type(ThreeBodyForce), intent(inout) :: thr
    type(sys) :: s

    if(s%find(this%file_3n,'.txt')) then
      call this%read_scalar_me3j_ascii_txt(thr)
      return
    end if

    if(s%find(this%file_3n,'.me3j.gz')) then
      call this%read_scalar_me3j_gzip(thr)
      return
    end if

    if(s%find(this%file_3n,'.me3j')) then
      call this%read_scalar_me3j_ascii(thr)
      return
    end if

    if(s%find(this%file_3n,'stream.bin')) then
      call this%read_scalar_me3j_binary_stream(thr)
      return
    end if

    if(s%find(this%file_3n,'.bin') .and. s%find(this%file_3n,'_comp')) then
      call this%read_scalar_me3j_binary_comp(thr)
      return
    end if

    if(s%find(this%file_3n,'.bin')) then
      call this%read_scalar_me3j_binary(thr)
      return
    end if

    call this%read_scalar_me3j_binary_stream(thr)
  end subroutine ReadScalar3BFile

  subroutine read_scalar_me3j_ascii_txt(this,thr)
    class(Read3BodyFiles), intent(in) :: this
    type(ThreeBodyForce), intent(inout) :: thr
    type(OrbitsIsospin) :: spsf
#ifdef single_precision_three_body_file
    real(4), allocatable :: v(:)
#else
    real(8), allocatable :: v(:)
#endif
    integer(8) :: nelm, n
    integer :: runit = 22, io

    write(*,'(a)') "Reading three-body scalar line-by-line from human-readable file"

    call spsf%init(this%emax3, this%lmax3)
    nelm = count_scalar_3bme(spsf, this%e2max3, this%e3max3)
    allocate(v(nelm))
    open(runit, file=this%file_3n, action='read', iostat=io)
    if(io /= 0) then
      write(*,'(2a)') "File opening error: ", trim(this%file_3n)
      return
    end if
    do n = 1, nelm
      read(runit,*) v(n)
    end do
    close(runit)

    call store_scalar_3bme(thr,v,spsf,this%e2max3,this%e3max3)

    deallocate(v)
    call spsf%fin()
  end subroutine read_scalar_me3j_ascii_txt

  subroutine read_scalar_me3j_ascii(this,thr)
    class(Read3BodyFiles), intent(in) :: this
    type(ThreeBodyForce), intent(inout) :: thr
    type(OrbitsIsospin) :: spsf
#ifdef single_precision_three_body_file
    real(4), allocatable :: v(:)
#else
    real(8), allocatable :: v(:)
#endif
    integer(8) :: nelm, n
    integer :: runit = 22

    write(*,'(a)') "Reading three-body scalar line-by-line from human-readable file"

    call spsf%init(this%emax3, this%lmax3)
    nelm = count_scalar_3bme(spsf, this%e2max3, this%e3max3)
    allocate(v(nelm))

    open(runit, file=this%file_3n, action='read')
    read(runit,*)
    do n = 1, nelm/10
      read(runit,*) v((n-1)*10+1 : n*10)
    end do

    ! basically, this is not needed
    if(nelm - (nelm/10) * 10 > 0) then
      read(runit,*) v((nelm/10)*10+1 : nelm)
    end if
    close(runit)

    call store_scalar_3bme(thr,v,spsf,this%e2max3,this%e3max3)

    deallocate(v)
    call spsf%fin()
  end subroutine read_scalar_me3j_ascii

  subroutine read_scalar_me3j_gzip(this,thr)
    use, intrinsic :: iso_c_binding
    use MyLibrary, only: gzip_open, gzip_readline, gzip_close
    class(Read3BodyFiles), intent(in) :: this
    type(ThreeBodyForce), intent(inout) :: thr
    type(OrbitsIsospin) :: spsf
#ifdef single_precision_three_body_file
    real(4), allocatable :: v(:)
#else
    real(8), allocatable :: v(:)
#endif
    integer(8) :: nelm, n
    character(256) :: header, buffer
    type(c_ptr) :: fp, err

    write(*,'(a)') "Reading three-body scalar line-by-line from gzip file"

    call spsf%init(this%emax3, this%lmax3)
    nelm = count_scalar_3bme(spsf, this%e2max3, this%e3max3)
    allocate(v(nelm))

    fp = gzip_open(this%file_3n, "rt")
    err = gzip_readline(fp, header, len(header))
    do n = 1, nelm/10
      err = gzip_readline(fp, buffer, len(buffer))
      read(buffer,*) v((n-1)*10+1 : n*10)
    end do

    ! basically, this is not needed
    if(nelm - (nelm/10) * 10 > 0) then
      err = gzip_readline(fp, buffer, len(buffer))
      read(buffer,*) v((nelm/10)*10+1 : nelm)
    end if
    err = gzip_close(fp)

    call store_scalar_3bme(thr,v,spsf,this%e2max3,this%e3max3)

    deallocate(v)
    call spsf%fin()
  end subroutine read_scalar_me3j_gzip

  subroutine read_scalar_me3j_binary(this,thr)
    class(Read3BodyFiles), intent(in) :: this
    type(ThreeBodyForce), intent(inout) :: thr
    type(OrbitsIsospin) :: spsf
#ifdef single_precision_three_body_file
    real(4), allocatable :: v(:)
#else
    real(8), allocatable :: v(:)
#endif
    integer(8) :: nelm, i
    integer :: runit = 22, io

    write(*,'(a)') "Reading three-body scalar line-by-line from binary file"
    call spsf%init(this%emax3, this%lmax3)
    nelm = count_scalar_3bme(spsf, this%e2max3, this%e3max3)
    write(*,*) "Number of three-body matrix element: ", nelm
    allocate(v(nelm))
    open(runit, file=this%file_3n, action='read', iostat=io, &
        & form='unformatted')
    if(io /= 0) then
      write(*,'(2a)') "File opening error: ", trim(this%file_3n)
      return
    end if
    do i = 1, nelm
      read(runit) v(i)
    end do
    close(runit)

    call store_scalar_3bme(thr,v,spsf,this%e2max3,this%e3max3)

    deallocate(v)
    call spsf%fin()
  end subroutine read_scalar_me3j_binary

  subroutine read_scalar_me3j_binary_comp(this,thr)
    class(Read3BodyFiles), intent(in) :: this
    type(ThreeBodyForce), intent(inout) :: thr
    type(OrbitsIsospin) :: spsf
#ifdef single_precision_three_body_file
    real(4), allocatable :: v(:)
#else
    real(8), allocatable :: v(:)
#endif
    integer(8) :: nelm
    integer :: runit = 22, io

    write(*,'(a)') "Reading three-body scalar from comp format binary file"
    call spsf%init(this%emax3, this%lmax3)
    nelm = count_scalar_3bme(spsf, this%e2max3, this%e3max3)
    allocate(v(nelm))
    open(runit, file=this%file_3n, action='read', iostat=io, &
        & form='unformatted')
    if(io /= 0) then
      write(*,'(2a)') "File opening error: ", trim(this%file_3n)
      return
    end if
    read(runit) v
    close(runit)

    call store_scalar_3bme(thr,v,spsf,this%e2max3,this%e3max3)

    deallocate(v)
    call spsf%fin()
  end subroutine read_scalar_me3j_binary_comp

  subroutine read_scalar_me3j_binary_stream(this,thr)
    class(Read3BodyFiles), intent(in) :: this
    type(ThreeBodyForce), intent(inout) :: thr
    type(OrbitsIsospin) :: spsf
#ifdef single_precision_three_body_file
    real(4), allocatable :: v(:)
#else
    real(8), allocatable :: v(:)
#endif
    integer(8) :: nelm
    integer :: runit = 22, io

    write(*,'(a)') "Reading three-body scalar from stream i/o format binary file"
    call spsf%init(this%emax3, this%lmax3)
    nelm = count_scalar_3bme(spsf, this%e2max3, this%e3max3)
    allocate(v(nelm))
    open(runit, file=this%file_3n, action='read', iostat=io, &
        & form='unformatted', access='stream')
    if(io /= 0) then
      write(*,'(2a)') "File opening error: ", trim(this%file_3n)
      return
    end if
    read(runit) v
    close(runit)

    call store_scalar_3bme(thr,v,spsf,this%e2max3,this%e3max3)

    deallocate(v)
    call spsf%fin()
  end subroutine read_scalar_me3j_binary_stream

  function count_scalar_3bme(spsf, e2max, e3max) result(r)
    type(OrbitsIsospin), intent(in) :: spsf
    integer, intent(in) :: e2max, e3max
    integer(8) :: r
    integer :: i1, l1, j1, e1
    integer :: i2, l2, j2, e2
    integer :: i3, l3, j3, e3
    integer :: i4, l4, j4, e4
    integer :: i5, l5, j5, e5, i5max
    integer :: i6, l6, j6, e6, i6max
    integer :: J12, T12, J45, T45, J, T, P123, P456

    r = 0
    do i1 = 1, spsf%norbs
      l1 = spsf%orb(i1)%l
      j1 = spsf%orb(i1)%j
      e1 = spsf%orb(i1)%e
      do i2 = 1, i1
        l2 = spsf%orb(i2)%l
        j2 = spsf%orb(i2)%j
        e2 = spsf%orb(i2)%e
        if(e1 + e2 > e2max) cycle
        do i3 = 1, i2
          l3 = spsf%orb(i3)%l
          j3 = spsf%orb(i3)%j
          e3 = spsf%orb(i3)%e
          if(e1 + e3 > e2max) cycle
          if(e2 + e3 > e2max) cycle
          if(e1 + e2 + e3 > e3max) cycle

          P123 = (-1) ** (l1+l2+l3)

          do i4 = 1, i1
            l4 = spsf%orb(i4)%l
            j4 = spsf%orb(i4)%j
            e4 = spsf%orb(i4)%e

            i5max = i4
            if(i1 == i4) i5max = i2

            do i5 = 1, i5max
              l5 = spsf%orb(i5)%l
              j5 = spsf%orb(i5)%j
              e5 = spsf%orb(i5)%e
              if(e4 + e5 > e2max) cycle

              i6max = i5
              if(i1 == i4 .and. i2 == i5) i6max = i3

              do i6 = 1, i6max
                l6 = spsf%orb(i6)%l
                j6 = spsf%orb(i6)%j
                e6 = spsf%orb(i6)%e
                if(e4 + e6 > e2max) cycle
                if(e5 + e6 > e2max) cycle
                if(e4 + e5 + e6 > e3max) cycle

                P456 = (-1) ** (l4+l5+l6)

                if(P123 /= P456) cycle
                do J12 = abs(j1-j2)/2, (j1+j2)/2
                  do J45 = abs(j4-j5)/2, (j4+j5)/2
                    do J = max(abs(2*J12-j3),abs(2*J45-j6)),&
                          &min(   (2*J12+j3),   (2*J45+j6)), 2

                      do T12 = 0, 1
                        do T45 = 0, 1
                          do T = max(abs(2*T12-1),abs(2*T45-1)),&
                                &min(   (2*T12+1),   (2*T45+1)), 2

                            r = r + 1

                          end do

                        end do
                      end do
                    end do
                  end do
                end do


              end do
            end do
          end do


        end do
      end do
    end do
  end function count_scalar_3bme

  subroutine store_scalar_3bme(thr,v,spsf,e2max,e3max)
    type(ThreeBodyForce), intent(inout) :: thr
    type(OrbitsIsospin), intent(in) :: spsf
#ifdef single_precision_three_body_file
    real(4), allocatable :: v(:)
#else
    real(8), allocatable :: v(:)
#endif
    integer, intent(in) :: e2max, e3max
    type(NonOrthIsospinThreeBodySpace), pointer :: ms
    type(OrbitsIsospin), pointer :: sps
    integer(8) :: cnt
    integer :: i1, l1, j1, e1
    integer :: i2, l2, j2, e2
    integer :: i3, l3, j3, e3
    integer :: i4, l4, j4, e4
    integer :: i5, l5, j5, e5, i5max
    integer :: i6, l6, j6, e6, i6max
    integer :: J12, T12, J45, T45, J, T, P123, P456
    integer :: ch, idxb, idxk, bra, ket

    ms => thr%thr
    sps => ms%isps

    cnt = 0
    do i1 = 1, spsf%norbs
      l1 = spsf%orb(i1)%l
      j1 = spsf%orb(i1)%j
      e1 = spsf%orb(i1)%e
      do i2 = 1, i1
        l2 = spsf%orb(i2)%l
        j2 = spsf%orb(i2)%j
        e2 = spsf%orb(i2)%e
        if(e1 + e2 > e2max) cycle
        do i3 = 1, i2
          l3 = spsf%orb(i3)%l
          j3 = spsf%orb(i3)%j
          e3 = spsf%orb(i3)%e
          if(e1 + e3 > e2max) cycle
          if(e2 + e3 > e2max) cycle
          if(e1 + e2 + e3 > e3max) cycle

          P123 = (-1) ** (l1+l2+l3)

          do i4 = 1, i1
            l4 = spsf%orb(i4)%l
            j4 = spsf%orb(i4)%j
            e4 = spsf%orb(i4)%e

            i5max = i4
            if(i1 == i4) i5max = i2

            do i5 = 1, i5max
              l5 = spsf%orb(i5)%l
              j5 = spsf%orb(i5)%j
              e5 = spsf%orb(i5)%e
              if(e4 + e5 > e2max) cycle

              i6max = i5
              if(i1 == i4 .and. i2 == i5) i6max = i3

              do i6 = 1, i6max
                l6 = spsf%orb(i6)%l
                j6 = spsf%orb(i6)%j
                e6 = spsf%orb(i6)%e
                if(e4 + e6 > e2max) cycle
                if(e5 + e6 > e2max) cycle
                if(e4 + e5 + e6 > e3max) cycle

                P456 = (-1) ** (l4+l5+l6)

                if(P123 /= P456) cycle
                do J12 = abs(j1-j2)/2, (j1+j2)/2
                  do J45 = abs(j4-j5)/2, (j4+j5)/2
                    do J = max(abs(2*J12-j3),abs(2*J45-j6)),&
                          &min(   (2*J12+j3),   (2*J45+j6)), 2

                      do T12 = 0, 1
                        do T45 = 0, 1
                          do T = max(abs(2*T12-1),abs(2*T45-1)),&
                                &min(   (2*T12+1),   (2*T45+1)), 2
                            cnt = cnt + 1

                            if(e1 > ms%emax) cycle
                            if(e2 > ms%emax) cycle
                            if(e3 > ms%emax) cycle

                            if(e4 > ms%emax) cycle
                            if(e5 > ms%emax) cycle
                            if(e6 > ms%emax) cycle

                            if(e1 + e2 > ms%e2max) cycle
                            if(e2 + e3 > ms%e2max) cycle
                            if(e3 + e1 > ms%e2max) cycle

                            if(e4 + e5 > ms%e2max) cycle
                            if(e5 + e6 > ms%e2max) cycle
                            if(e6 + e4 > ms%e2max) cycle

                            if(e1 + e2 + e3 > ms%e3max) cycle
                            if(e4 + e5 + e6 > ms%e3max) cycle

                            if(i1==i2 .and. mod(J12+T12,2)==0) then
                              if(abs(v(cnt)) > 1.d-6) then
                                write(*,*) "Warning: something wrong, this three-body matrix element has to be zero."
                              end if
                              cycle
                            end if

                            if(i4==i5 .and. mod(J45+T45,2)==0) then
                              if(abs(v(cnt)) > 1.d-6) then
                                write(*,*) "Warning: something wrong, this three-body matrix element has to be zero."
                              end if
                              cycle
                            end if

                            ch = ms%jpt2ch(J,P123,T)
                            idxb = ms%jpt(ch)%spis2idx(i1,i2,i3)
                            idxk = ms%jpt(ch)%spis2idx(i4,i5,i6)
                            if(idxb*idxk == 0) cycle
                            bra = ms%jpt(ch)%idxqn(idxb)%JT2n(J12,T12)
                            ket = ms%jpt(ch)%idxqn(idxk)%JT2n(J45,T45)
                            if(bra*ket == 0) cycle
                            thr%MatCh(ch,ch)%v(bra,ket) = v(cnt)
                            thr%MatCh(ch,ch)%v(ket,bra) = v(cnt)
                          end do

                        end do
                      end do
                    end do
                  end do
                end do


              end do
            end do
          end do


        end do
      end do
    end do
  end subroutine store_scalar_3bme

  subroutine store_scalar_3bme_sp(thr,v,spsf,e2max,e3max)
    type(ThreeBodyForce), intent(inout) :: thr
    type(OrbitsIsospin), intent(in) :: spsf
#ifdef single_precision_three_body_file
    real(4), intent(in) :: v(:)
#else
    real(8), intent(in) :: v(:)
#endif
    integer, intent(in) :: e2max, e3max
    type(NonOrthIsospinThreeBodySpace), pointer :: ms
    type(OrbitsIsospin), pointer :: sps
    integer(8) :: cnt
    integer :: i1, l1, j1, e1
    integer :: i2, l2, j2, e2
    integer :: i3, l3, j3, e3
    integer :: i4, l4, j4, e4
    integer :: i5, l5, j5, e5, i5max
    integer :: i6, l6, j6, e6, i6max
    integer :: J12, T12, J45, T45, J, T, P123, P456
    integer :: ch, idxb, idxk, bra, ket

    ms => thr%thr
    sps => ms%isps
    cnt = 0
    do i1 = 1, spsf%norbs
      l1 = spsf%orb(i1)%l
      j1 = spsf%orb(i1)%j
      e1 = spsf%orb(i1)%e
      do i2 = 1, i1
        l2 = spsf%orb(i2)%l
        j2 = spsf%orb(i2)%j
        e2 = spsf%orb(i2)%e
        if(e1 + e2 > e2max) cycle
        do i3 = 1, i2
          l3 = spsf%orb(i3)%l
          j3 = spsf%orb(i3)%j
          e3 = spsf%orb(i3)%e
          if(e1 + e3 > e2max) cycle
          if(e2 + e3 > e2max) cycle
          if(e1 + e2 + e3 > e3max) cycle

          P123 = (-1) ** (l1+l2+l3)

          do i4 = 1, i1
            l4 = spsf%orb(i4)%l
            j4 = spsf%orb(i4)%j
            e4 = spsf%orb(i4)%e

            i5max = i4
            if(i1 == i4) i5max = i2

            do i5 = 1, i5max
              l5 = spsf%orb(i5)%l
              j5 = spsf%orb(i5)%j
              e5 = spsf%orb(i5)%e
              if(e4 + e5 > e2max) cycle

              i6max = i5
              if(i1 == i4 .and. i2 == i5) i6max = i3

              do i6 = 1, i6max
                l6 = spsf%orb(i6)%l
                j6 = spsf%orb(i6)%j
                e6 = spsf%orb(i6)%e
                if(e4 + e6 > e2max) cycle
                if(e5 + e6 > e2max) cycle
                if(e4 + e5 + e6 > e3max) cycle

                P456 = (-1) ** (l4+l5+l6)

                if(P123 /= P456) cycle
                do J12 = abs(j1-j2)/2, (j1+j2)/2
                  do J45 = abs(j4-j5)/2, (j4+j5)/2
                    do J = max(abs(2*J12-j3),abs(2*J45-j6)),&
                          &min(   (2*J12+j3),   (2*J45+j6)), 2

                      do T12 = 0, 1
                        do T45 = 0, 1
                          do T = max(abs(2*T12-1),abs(2*T45-1)),&
                                &min(   (2*T12+1),   (2*T45+1)), 2
                            cnt = cnt + 1

                            if(e1 > ms%emax) cycle
                            if(e2 > ms%emax) cycle
                            if(e3 > ms%emax) cycle

                            if(e4 > ms%emax) cycle
                            if(e5 > ms%emax) cycle
                            if(e6 > ms%emax) cycle

                            if(e1 + e2 > ms%e2max) cycle
                            if(e2 + e3 > ms%e2max) cycle
                            if(e3 + e1 > ms%e2max) cycle

                            if(e4 + e5 > ms%e2max) cycle
                            if(e5 + e6 > ms%e2max) cycle
                            if(e6 + e4 > ms%e2max) cycle

                            if(e1 + e2 + e3 > ms%e3max) cycle
                            if(e4 + e5 + e6 > ms%e3max) cycle

                            if(i1==i2 .and. mod(J12+T12,2)==0) then
                              if(abs(v(cnt)) > 1.d-6) then
                                write(*,*) "Warning: something wrong, this three-body matrix element has to be zero."
                              end if
                              cycle
                            end if

                            if(i4==i5 .and. mod(J45+T45,2)==0) then
                              if(abs(v(cnt)) > 1.d-6) then
                                write(*,*) "Warning: something wrong, this three-body matrix element has to be zero."
                              end if
                              cycle
                            end if
                            ch = ms%jpt2ch(J,P123,T)
                            idxb = ms%jpt(ch)%spis2idx(i1,i2,i3)
                            idxk = ms%jpt(ch)%spis2idx(i4,i5,i6)
                            if(idxb*idxk == 0) cycle
                            bra = ms%jpt(ch)%idxqn(idxb)%JT2n(J12,T12)
                            ket = ms%jpt(ch)%idxqn(idxk)%JT2n(J45,T45)
                            if(bra*ket == 0) cycle
                            thr%MatCh(ch,ch)%v(bra,ket) = v(cnt)
                            thr%MatCh(ch,ch)%v(ket,bra) = v(cnt)
                          end do

                        end do
                      end do
                    end do
                  end do
                end do


              end do
            end do
          end do


        end do
      end do
    end do
  end subroutine store_scalar_3bme_sp

  !
  !
  ! reading three-body tensor
  !
  !

  subroutine ReadTensor3BFile(this,thr)
    use ClassSys, only: sys
    class(Read3BodyFiles), intent(in) :: this
    type(ThreeBodyForce), intent(inout) :: thr
    type(sys) :: s

    if(s%find(this%file_3n,'.txt')) then
      !call this%read_tensor_3bme_ascii(thr)
      return
    end if

    if(s%find(this%file_3n,'.bin')) then
      !call this%read_tensor_3bme_bin(thr)
      return
    end if
  end subroutine ReadTensor3BFile

  subroutine read_tensor_3bme_ascii(this,thr)
    class(Read3BodyFiles), intent(in) :: this
    type(ThreeBodyForce), intent(inout) :: thr
    write(*,*) "reading tensor is not ready"
    return
  end subroutine read_tensor_3bme_ascii

  subroutine read_tensor_3bme_bin(this,thr)
    class(Read3BodyFiles), intent(in) :: this
    type(ThreeBodyForce), intent(inout) :: thr
    write(*,*) "reading tensor is not ready"
    return
  end subroutine read_tensor_3bme_bin
end module ThreeBodyInteraction
