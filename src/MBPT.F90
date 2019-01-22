module MBPT
  use omp_lib
  use ModelSpace
  use Operators
  implicit none

  public :: MBPTEnergy

  private :: CalcEnergyCorr
  private :: energy_second
  private :: denom
  private :: Eket
  private :: Pari
  private :: Tz
  private :: energy_third_pp
  private :: energy_third_hh
  private :: energy_third_ph
  private :: cross_couple

  type :: MBPTEnergy
    real(8) :: e_0 = 0.d0
    real(8) :: e_2 = 0.d0
    real(8) :: e_3 = 0.d0
    real(8) :: e_3_pp = 0.d0
    real(8) :: e_3_hh = 0.d0
    real(8) :: e_3_ph = 0.d0
  contains
    procedure :: calc => CalcEnergyCorr
    procedure :: energy_second
    procedure :: energy_third
    procedure :: energy_third_pp
    procedure :: energy_third_hh
    procedure :: energy_third_ph
    procedure :: write_summary
  end type MBPTEnergy
contains

  subroutine write_summary(this,p)
    use HFInput
    class(MBPTEnergy), intent(inout) :: this
    type(InputParameters), intent(in) :: p
    integer :: wunit = 33

    open(wunit, file = 'SummaryHFMBPT.txt', action='write')
    call p%PrintInputParameters(wunit)
    write(wunit,'(4f18.8)') this%e_0, this%e_2, this%e_3, &
        &           this%e_0+this%e_2+this%e_3
    close(wunit)
  end subroutine write_summary

  subroutine CalcEnergyCorr(this,ms,hamil)
    class(MBPTEnergy), intent(inout) :: this
    type(MSpace), intent(in) :: ms
    type(Op), intent(in) :: hamil

    write(*,*)
    write(*,'(a)') " Many-body perturbation calculation up to 3rd order"
    write(*,*)

    this%e_0 = hamil%zero
    call this%energy_second(ms,hamil)
    write(*,'(a,f16.8)') "Second order correction: ", this%e_2

    call this%energy_third(ms,hamil)
    write(*,'(a,4f16.8)') "Third order corrections pp, hh, ph, and total: ", &
        & this%e_3_pp, this%e_3_hh, this%e_3_ph, this%e_3
    write(*,'(a,f16.8)') "Third order MBPT energy: ", this%e_0 + this%e_2 + this%e_3

  end subroutine CalcEnergyCorr

  subroutine energy_second(this,ms,h)
    ! a, b : particle
    ! i, j : hole
    !     _____________
    !    /\           /\
    !   /  \         /  \
    !   |  |         |  |
    ! a |  | i     b |  | j
    !   |  |         |  |
    !   \  /         \  /
    !    \/___________\/
    !
    ! \sum_{i>j,a>b} <ij||ab> <ab||ij> / denominator
    use Profiler, only: timer
    class(MBPTEnergy), intent(inout) :: this
    type(MSPace), intent(in) :: ms
    type(Op), intent(in) :: h
    integer :: ch, ab, ij, J2, n
    integer :: a, b, i, j
    real(8) :: vsum, v, ti

    ti = omp_get_wtime()

    vsum = 0.d0
    do ch = 1, ms%two%NChan
      J2 = ms%two%jpz(ch)%j
      n = ms%two%jpz(ch)%nst

      v = 0.d0
      !$omp parallel
      !$omp do private(ab,a,b,ij,i,j) reduction(+:v)
      do ab = 1, n
        a = ms%two%jpz(ch)%n2spi1(ab)
        b = ms%two%jpz(ch)%n2spi2(ab)
        if( ms%sps%orb(a)%ph /= 1 ) cycle
        if( ms%sps%orb(b)%ph /= 1 ) cycle
        do ij = 1, n
          i = ms%two%jpz(ch)%n2spi1(ij)
          j = ms%two%jpz(ch)%n2spi2(ij)
          if( ms%sps%orb(i)%ph /= 0 ) cycle
          if( ms%sps%orb(j)%ph /= 0 ) cycle

          v = v + h%two%GetTwBME(ms%sps,ms%two,i,j,a,b,J2)**2 / &
              & denom(ms%sps,ms%one,h%one,i,j,a,b)
        end do
      end do
      !$omp end do
      !$omp end parallel
      vsum = vsum + v * dble(2*J2+1)
    end do

    this%e_2 = vsum

    call timer%Add("Second order MBPT",omp_get_wtime()-ti)
  end subroutine energy_second

  function denom(sps,one,f,h1,h2,p1,p2) result(r)
    type(Orbits), intent(in) :: sps
    type(OneBodySpace), intent(in) :: one
    type(NBodyPart), intent(in) :: f
    integer, intent(in) :: h1, h2, p1, p2
    real(8) :: r

    r = f%GetOBME(sps,one,h1,h1) + f%GetOBME(sps,one,h2,h2) - &
        & f%GetOBME(sps,one,p1,p1) - f%GetOBME(sps,one,p2,p2)
  end function denom

  function Eket(sps,a,b)
    type(Orbits), intent(in) :: sps
    integer, intent(in) :: a, b
    integer :: Eket
    Eket = sps%orb(a)%e + sps%orb(b)%e
  end function Eket

  function Pari(sps,a,b) result(p)
    type(Orbits), intent(in) :: sps
    integer, intent(in) :: a, b
    integer :: p
    p = (-1) ** (sps%orb(a)%l + sps%orb(b)%l)
  end function Pari

  function Tz(sps,a,b)
    type(Orbits), intent(in) :: sps
    integer, intent(in) :: a, b
    integer :: Tz
    Tz = sps%orb(a)%z + sps%orb(b)%z
  end function Tz

  subroutine energy_third(this,ms,h)
    class(MBPTEnergy), intent(inout) :: this
    type(MSPace), intent(in) :: ms
    type(Op), intent(in) :: h

    call this%energy_third_hh(ms,h)
    call this%energy_third_pp(ms,h)
    call this%energy_third_ph(ms,h)
    this%e_3 = this%e_3_hh + this%e_3_pp + this%e_3_ph
  end subroutine energy_third

  subroutine energy_third_pp(this,ms,h)
    ! a, b, c, d: particle
    ! i, j      : hole
    !     _____________
    !    /\           /\
    !   /  \         /  \
    !   |  | a     b |  |
    ! i |  |_________|  | j
    !   |  |         |  |
    !   |  | c     d |  |
    !   \  /         \  /
    !    \/___________\/
    !
    ! <ij||ab> <ab||cd> <cd||ij> / 8 denominator
    use Profiler, only: timer
    class(MBPTEnergy), intent(inout) :: this
    type(MSPace), intent(in) :: ms
    type(Op), intent(in) :: h
    integer :: ch, J2, n, ab, cd, ij
    integer :: i, j, a, b, c, d
    real(8) :: v, vsum, ti

    ti = omp_get_wtime()
    vsum = 0.d0
    do ch = 1, ms%two%NChan
      J2 = ms%two%jpz(ch)%j
      n = ms%two%jpz(ch)%nst

      v = 0.d0
      !$omp parallel
      !$omp do private(ij,i,j,ab,a,b,cd,c,d) reduction(+:v)
      do ab = 1, n
        a = ms%two%jpz(ch)%n2spi1(ab)
        b = ms%two%jpz(ch)%n2spi2(ab)
        if( ms%sps%orb(a)%ph /= 1 ) cycle
        if( ms%sps%orb(b)%ph /= 1 ) cycle
        do cd = 1, n
          c = ms%two%jpz(ch)%n2spi1(cd)
          d = ms%two%jpz(ch)%n2spi2(cd)
          if( ms%sps%orb(c)%ph /= 1 ) cycle
          if( ms%sps%orb(d)%ph /= 1 ) cycle
          do ij = 1, n
            i = ms%two%jpz(ch)%n2spi1(ij)
            j = ms%two%jpz(ch)%n2spi2(ij)
            if( ms%sps%orb(i)%ph /= 0 ) cycle
            if( ms%sps%orb(j)%ph /= 0 ) cycle

            v = v + h%two%GetTwBME(ms%sps,ms%two,i,j,a,b,J2) * &
                & h%two%GetTwBME(ms%sps,ms%two,a,b,c,d,J2) * &
                & h%two%GetTwBME(ms%sps,ms%two,c,d,i,j,J2) / &
                & ( denom(ms%sps,ms%one,h%one,i,j,a,b) * &
                & denom(ms%sps,ms%one,h%one,i,j,c,d) )
          end do
        end do
      end do
      !$omp end do
      !$omp end parallel

      vsum = vsum + v * dble(2*J2+1)
    end do
    this%e_3_pp = vsum
    call timer%Add("Third order MBPT pp ladder",omp_get_wtime()-ti)
  end subroutine energy_third_pp

  subroutine energy_third_hh(this,ms,h)
    ! a, b      : particle
    ! i, j, k, l: hole
    !     _____________
    !    /\           /\
    !   /  \         /  \
    !   |  | i     j |  |
    ! a |  |_________|  | b
    !   |  |         |  |
    !   |  | k     l |  |
    !   \  /         \  /
    !    \/___________\/
    !
    ! <ij||ab> <kl||ij> <ab||kl> / 8 denominator
    use Profiler, only: timer
    class(MBPTEnergy), intent(inout) :: this
    type(MSPace), intent(in) :: ms
    type(Op), intent(in) :: h
    integer :: ch, J2, n, ab, ij, kl
    integer :: a, b, i, j, k, l
    real(8) :: v, vsum, ti

    ti = omp_get_wtime()

    vsum = 0.d0
    do ch = 1, ms%two%NChan
      J2 = ms%two%jpz(ch)%j
      n = ms%two%jpz(ch)%nst

      v = 0.d0
      !$omp parallel
      !$omp do private(ab,a,b,ij,i,j,kl,k,l) reduction(+:v)
      do ab = 1, n
        a = ms%two%jpz(ch)%n2spi1(ab)
        b = ms%two%jpz(ch)%n2spi2(ab)
        if( ms%sps%orb(a)%ph /= 1 ) cycle
        if( ms%sps%orb(b)%ph /= 1 ) cycle
        do ij = 1, n
          i = ms%two%jpz(ch)%n2spi1(ij)
          j = ms%two%jpz(ch)%n2spi2(ij)
          if( ms%sps%orb(i)%ph /= 0 ) cycle
          if( ms%sps%orb(j)%ph /= 0 ) cycle
          do kl = 1, n
            k = ms%two%jpz(ch)%n2spi1(kl)
            l = ms%two%jpz(ch)%n2spi2(kl)
            if( ms%sps%orb(k)%ph /= 0 ) cycle
            if( ms%sps%orb(l)%ph /= 0 ) cycle
            v = v + h%two%GetTwBME(ms%sps,ms%two,a,b,i,j,J2) * &
                & h%two%GetTwBME(ms%sps,ms%two,i,j,k,l,J2) * &
                & h%two%GetTwBME(ms%sps,ms%two,k,l,a,b,J2) / &
                & ( denom(ms%sps,ms%one,h%one,i,j,a,b) * &
                & denom(ms%sps,ms%one,h%one,k,l,a,b) )
          end do
        end do
      end do
      !$omp end do
      !$omp end parallel

      vsum = vsum + dble(2*J2+1) * v
    end do

    this%e_3_hh = vsum
    call timer%Add("Third order MBPT hh ladder",omp_get_wtime()-ti)
  end subroutine energy_third_hh

  subroutine energy_third_ph(this,ms,h)
    ! a, b, c : particle
    ! i, j, k : hole
    !     _____________
    !    /\           /\
    !   /  \         /  \
    !   |  | i     b |  |
    ! a |  |_________|  | j
    !   |  |         |  |
    !   |  | k     c |  |
    !   \  /         \  /
    !    \/___________\/
    !
    ! - <ij||ab> <kb||ic> <ac||kj> / denominator
    ! = - \sum_{L} [L] <ij|X|ab>_{L} <kb|X|ic>_{L} <ac|X|kj>_{L} / denominator
    ! <ij|X|ab>_{L} = \sum_{A} [A] {i j A} <ij:A|X|ab:A>
    !                              {a b L}
    use Profiler, only: timer
    use CommonLibrary, only: triag
    class(MBPTEnergy), intent(inout) :: this
    type(MSPace), intent(in) :: ms
    type(Op), intent(in) :: h
    integer :: J2
    integer :: a, b, c, i, j, k
    integer :: ia, ib, ic, ii, ij, ik
    integer :: ja, jb, jc, ji, jj, jk
    real(8) :: v, vsum, norm, ti

    ti = omp_get_wtime()

    vsum = 0.d0
    !$omp parallel
    !$omp do private(ia,ib,ic,ii,ij,ik, &
    !$omp &  a,b,c,i,j,k,ja,jb,jc,ji,jj,jk,norm,v) reduction(+:vsum)
    do ia = 1, size(ms%particles)
      do ib = 1, size(ms%particles)
        a = ms%particles(ia)
        b = ms%particles(ib)
        if(Eket(ms%sps,a,b) > ms%two%e2max) cycle
        do ii = 1, size(ms%holes)
          do ij = 1, size(ms%holes)
            i = ms%holes(ii)
            j = ms%holes(ij)
            if(Eket(ms%sps,i,j) > ms%two%e2max) cycle
            if(Pari(ms%sps,i,j) /= Pari(ms%sps,a,b)) cycle
            if(  Tz(ms%sps,i,j) /=   Tz(ms%sps,a,b)) cycle

            do ic = 1, size(ms%particles)
              do ik = 1, size(ms%holes)
                c = ms%particles(ic)
                k = ms%holes(ik)


                if(Eket(ms%sps,k,b) > ms%two%e2max) cycle
                if(Eket(ms%sps,i,c) > ms%two%e2max) cycle
                if(Eket(ms%sps,a,c) > ms%two%e2max) cycle
                if(Eket(ms%sps,k,j) > ms%two%e2max) cycle

                if(Pari(ms%sps,k,b) /= Pari(ms%sps,i,c)) cycle
                if(Pari(ms%sps,a,c) /= Pari(ms%sps,k,j)) cycle

                if(  Tz(ms%sps,k,b) /=   Tz(ms%sps,i,c)) cycle
                if(  Tz(ms%sps,a,c) /=   Tz(ms%sps,k,j)) cycle


                norm = 1.d0
                if(i==j) norm = norm * sqrt(2.d0)
                if(j==k) norm = norm * sqrt(2.d0)
                if(a==b) norm = norm * sqrt(2.d0)
                if(a==c) norm = norm * sqrt(2.d0)

                ja = ms%sps%orb(a)%j
                jb = ms%sps%orb(b)%j
                jc = ms%sps%orb(c)%j
                ji = ms%sps%orb(i)%j
                jj = ms%sps%orb(j)%j
                jk = ms%sps%orb(k)%j

                v = 0.d0
                do J2 = max(abs(ja-jj),abs(jb-ji),abs(jc-jk))/2, &
                      & min(   (ja+jj),   (jb+ji),   (jc+jk))/2
                  v = v + dble(2*J2+1) * &
                    & cross_couple(ms%sps,ms%two,h%two,i,j,a,b,J2) * &
                    & cross_couple(ms%sps,ms%two,h%two,k,b,i,c,J2) * &
                    & cross_couple(ms%sps,ms%two,h%two,a,c,k,j,J2)

                end do

                vsum = vsum + norm * v / &
                    & ( denom(ms%sps,ms%one,h%one,k,j,a,c) * &
                    &   denom(ms%sps,ms%one,h%one,i,j,a,b) )

              end do
            end do
          end do

        end do
      end do
    end do
    !$omp end do
    !$omp end parallel

    this%e_3_ph = - vsum
    call timer%Add("Third order MBPT ph ladder",omp_get_wtime()-ti)
  end subroutine energy_third_ph

  function cross_couple(sps,two,v,i,j,a,b,L) result(r)
    use CommonLibrary, only: triag,sjs
    type(Orbits), intent(in) :: sps
    type(TwoBodySpace), intent(in) :: two
    type(NBodyPart), intent(in) :: v
    integer, intent(in) :: i, j, a, b, L
    integer :: ji, jj, ja, jb
    integer :: J2
    real(8) :: r

    r = 0.d0

    if(Eket(sps,i,j) > two%e2max) return
    if(Eket(sps,a,b) > two%e2max) return
    if(Pari(sps,i,j) /= Pari(sps,a,b)) return
    if(  Tz(sps,i,j) /=   Tz(sps,a,b)) return

    ji = sps%orb(i)%j
    jj = sps%orb(j)%j
    ja = sps%orb(a)%j
    jb = sps%orb(b)%j

    if(triag(ji,jb,2*L)) return
    if(triag(ja,jj,2*L)) return

    do J2 = max( abs(ji-jj),abs(ja-jb) )/2, min( (ji+jj),(ja+jb) )/2
      if(i==j .and. mod(J2,2)==1) cycle
      if(a==b .and. mod(J2,2)==1) cycle
      r = r + dble(2*J2+1) * &
          & sjs(ji, jj, 2*J2, ja, jb, 2*L) * &
          & v%GetTwBME(sps,two,i,j,a,b,J2)
    end do
  end function cross_couple

end module MBPT
