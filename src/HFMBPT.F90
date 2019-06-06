module HFMBPT
  use omp_lib
  use ModelSpace
  use Operators
  use StoreCouplings
  use HartreeFock
  implicit none

  public :: MBPTEnergy, MBPTScalar, MBPTDMat

  private :: CalcEnergyCorr
  private :: energy_second
  private :: denom1b
  private :: denom2b
  private :: Eket
  private :: Pari
  private :: Tz
  private :: energy_third_pp
  private :: energy_third_hh
  private :: energy_third_ph
  private :: cross_couple
  private :: MBPTCriteria

  private :: CalcScalarCorr
  private :: scalar_first
  private :: scalar_second
  private :: scalar_second_s1p
  private :: scalar_second_s1h
  private :: scalar_second_s1ph
  private :: scalar_second_s2pp
  private :: scalar_second_s2hh
  private :: scalar_second_s2ph
  private :: scalar_second_v2pp
  private :: scalar_second_v2hh
  private :: scalar_second_v2ph

  type :: MBPTEnergy
    real(8) :: e_0 = 0.d0
    real(8) :: e_2 = 0.d0
    real(8) :: e_3 = 0.d0
    real(8) :: e_3_pp = 0.d0
    real(8) :: e_3_hh = 0.d0
    real(8) :: e_3_ph = 0.d0
    real(8) :: perturbativity1b = 0.d0
    real(8) :: perturbativity2b = 0.d0
    real(8) :: energy_gap = 0.d0
  contains
    procedure :: calc => CalcEnergyCorr
    procedure :: energy_second
    procedure :: energy_third
    procedure :: energy_third_pp
    procedure :: energy_third_hh
    procedure :: energy_third_ph
    procedure :: MBPTCriteria
  end type MBPTEnergy

  type :: MBPTDMat
    type(OneBodyPart) :: rho
    type(OneBodyPart) :: C_HO2HF
    type(OneBodyPart) :: C_HF2NAT
    type(OneBodyPart) :: C_HO2NAT
    type(OneBodyPart) :: Occ
    type(MSpace), pointer :: ms
  contains
    procedure :: init => InitMBPTDMat
    procedure :: fin => FinMBPTDMat
    procedure :: density_matrix_pp
    procedure :: density_matrix_hh
    procedure :: density_matrix_ph
    procedure :: GetCoef
  end type MBPTDMat

  type :: MBPTScalar
    real(8) :: s_0 = 0.d0
    real(8) :: s_1 = 0.d0
    real(8) :: s_2 = 0.d0
    real(8) :: s_2_s2pp = 0.d0
    real(8) :: s_2_s2hh = 0.d0
    real(8) :: s_2_s2ph = 0.d0
    real(8) :: s_2_s1h  = 0.d0
    real(8) :: s_2_s1p  = 0.d0
    real(8) :: s_2_s1ph = 0.d0
    real(8) :: s_2_v2pp = 0.d0
    real(8) :: s_2_v2hh = 0.d0
    real(8) :: s_2_v2ph = 0.d0
  contains
    procedure :: calc => CalcScalarCorr
    procedure :: scalar_first
    procedure :: scalar_second
    procedure :: scalar_second_s1p
    procedure :: scalar_second_s1h
    procedure :: scalar_second_s1ph
    procedure :: scalar_second_s2pp
    procedure :: scalar_second_s2hh
    procedure :: scalar_second_s2ph
    procedure :: scalar_second_v2pp
    procedure :: scalar_second_v2hh
    procedure :: scalar_second_v2ph
  end type MBPTScalar
  type(SixJsStore), private :: sixjs
contains

  subroutine CalcEnergyCorr(this,hamil)
    class(MBPTEnergy), intent(inout) :: this
    type(Ops), intent(in) :: hamil
    integer :: jmax

    write(*,*)
    write(*,'(a)') " Many-body perturbation calculation up to 3rd order"
    write(*,*)

    if(.not. hamil%is_normal_ordered) then
      write(*,"(a)") "In CalcEnergyCorr: "
      write(*,"(a)") " Hamiltonian has to be normal ordered"
      return
    end if

    jmax = 2*hamil%ms%sps%lmax+1
    call sixjs%init(1,jmax,.true., 1,jmax,.true., 1,jmax,.true.)
    call this%MBPTCriteria(hamil)

    this%e_0 = hamil%zero
    call this%energy_second(hamil)
    write(*,'(a,f16.8)') "Second order correction: ", this%e_2

    call this%energy_third(hamil)
    write(*,'(a,4f16.8)') "Third order corrections pp, hh, ph, and total: ", &
        & this%e_3_pp, this%e_3_hh, this%e_3_ph, this%e_3
    write(*,'(a,f16.8)') "Third order MBPT energy: ", this%e_0 + this%e_2 + this%e_3
    call sixjs%fin()

  end subroutine CalcEnergyCorr

  subroutine energy_second(this,h)
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
    type(Ops), intent(in) :: h
    type(MSPace), pointer :: ms
    type(TwoBodyChannel), pointer :: ch_two
    integer :: ch, ab, ij, J2, n
    integer :: a, b, i, j
    real(8) :: vsum, v, ti

    ti = omp_get_wtime()
    ms => h%ms
    vsum = 0.d0
    do ch = 1, ms%two%NChan
      ch_two => ms%two%jpz(ch)
      J2 = ch_two%j
      n = ch_two%n_state
      if(ch_two%n_hh_state < 1) cycle
      if(ch_two%n_pp_state < 1) cycle

      v = 0.d0
      !$omp parallel
      !$omp do private(ab,a,b,ij,i,j) reduction(+:v)
      do ab = 1, ch_two%n_pp_state
        a = ch_two%n2spi1(ch_two%pps(ab))
        b = ch_two%n2spi2(ch_two%pps(ab))
        do ij = 1, ch_two%n_hh_state
          i = ch_two%n2spi1(ch_two%hhs(ij))
          j = ch_two%n2spi2(ch_two%hhs(ij))

          v = v + h%two%GetTwBME(i,j,a,b,J2)**2 / &
              & denom2b(h%one,i,j,a,b)
        end do
      end do
      !$omp end do
      !$omp end parallel
      vsum = vsum + v * dble(2*J2+1)
    end do

    this%e_2 = vsum

    call timer%Add("Second order MBPT",omp_get_wtime()-ti)
  end subroutine energy_second

  function denom1b(f,h1,p1) result(r)
    type(OneBodyPart), intent(in) :: f
    integer, intent(in) :: h1, p1
    real(8) :: r

    r = f%GetOBME(h1,h1) - f%GetOBME(p1,p1)
  end function denom1b

  function denom2b(f,h1,h2,p1,p2) result(r)
    type(OneBodyPart), intent(in) :: f
    integer, intent(in) :: h1, h2, p1, p2
    real(8) :: r

    r = f%GetOBME(h1,h1) + f%GetOBME(h2,h2) - &
        & f%GetOBME(p1,p1) - f%GetOBME(p2,p2)
  end function denom2b

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

  subroutine energy_third(this,h)
    class(MBPTEnergy), intent(inout) :: this
    type(Ops), intent(in) :: h

    call this%energy_third_hh(h)
    call this%energy_third_pp(h)
    call this%energy_third_ph(h)
    this%e_3 = this%e_3_hh + this%e_3_pp + this%e_3_ph
  end subroutine energy_third

  subroutine energy_third_pp(this,h)
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
    type(Ops), intent(in) :: h
    type(MSPace), pointer :: ms
    type(TwoBodyChannel), pointer :: ch_two
    integer :: ch, J2, n, ab, cd, ij
    integer :: i, j, a, b, c, d
    real(8) :: v, vsum, ti

    ms => h%ms
    ti = omp_get_wtime()
    vsum = 0.d0
    do ch = 1, ms%two%NChan
      J2 = ms%two%jpz(ch)%j
      n = ms%two%jpz(ch)%n_state

      ch_two => ms%two%jpz(ch)
      if(ch_two%n_hh_state < 1) cycle
      if(ch_two%n_pp_state < 1) cycle

      v = 0.d0
      !$omp parallel
      !$omp do private(ab,a,b,cd,c,d,ij,i,j) reduction(+:v)
      do ab = 1, ch_two%n_pp_state
        a = ch_two%n2spi1(ch_two%pps(ab))
        b = ch_two%n2spi2(ch_two%pps(ab))
        do cd = 1, ch_two%n_pp_state
          c = ch_two%n2spi1(ch_two%pps(cd))
          d = ch_two%n2spi2(ch_two%pps(cd))
          do ij = 1, ch_two%n_hh_state
            i = ch_two%n2spi1(ch_two%hhs(ij))
            j = ch_two%n2spi2(ch_two%hhs(ij))

            v = v + h%two%GetTwBME(i,j,a,b,J2) * &
                & h%two%GetTwBME(a,b,c,d,J2) * &
                & h%two%GetTwBME(c,d,i,j,J2) / &
                & ( denom2b(h%one,i,j,a,b) * &
                &   denom2b(h%one,i,j,c,d) )
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

  subroutine energy_third_hh(this,h)
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
    type(Ops), intent(in) :: h
    type(MSPace), pointer :: ms
    type(TwoBodyChannel), pointer :: ch_two
    integer :: ch, J2, n, ab, ij, kl
    integer :: a, b, i, j, k, l
    real(8) :: v, vsum, ti

    ti = omp_get_wtime()

    ms => h%ms
    vsum = 0.d0
    do ch = 1, ms%two%NChan
      J2 = ms%two%jpz(ch)%j
      n = ms%two%jpz(ch)%n_state

      ch_two => ms%two%jpz(ch)
      if(ch_two%n_hh_state < 1) cycle
      if(ch_two%n_pp_state < 1) cycle

      v = 0.d0
      !$omp parallel
      !$omp do private(ab,a,b,ij,i,j,kl,k,l) reduction(+:v)
      do ab = 1, ch_two%n_pp_state
        a = ch_two%n2spi1(ch_two%pps(ab))
        b = ch_two%n2spi2(ch_two%pps(ab))
        do ij = 1, ch_two%n_hh_state
          i = ch_two%n2spi1(ch_two%hhs(ij))
          j = ch_two%n2spi2(ch_two%hhs(ij))
          do kl = 1, ch_two%n_hh_state
            k = ch_two%n2spi1(ch_two%hhs(kl))
            l = ch_two%n2spi2(ch_two%hhs(kl))

            v = v + h%two%GetTwBME(a,b,i,j,J2) * &
                & h%two%GetTwBME(i,j,k,l,J2) * &
                & h%two%GetTwBME(k,l,a,b,J2) / &
                & ( denom2b(h%one,i,j,a,b) * &
                &   denom2b(h%one,k,l,a,b) )
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

  subroutine energy_third_ph(this,h)
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
    use MyLibrary, only: triag
    class(MBPTEnergy), intent(inout) :: this
    type(Ops), intent(in) :: h
    type(MSPace), pointer :: ms
    integer :: J2
    integer :: a, b, c, i, j, k
    integer :: ia, ib, ic, ii, ij, ik
    integer :: ja, jb, jc, ji, jj, jk
    real(8) :: v, vsum, norm, ti

    ti = omp_get_wtime()
    ms => h%ms
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
                    & cross_couple(h%two,i,j,a,b,J2) * &
                    & cross_couple(h%two,k,b,i,c,J2) * &
                    & cross_couple(h%two,a,c,k,j,J2)

                end do

                vsum = vsum + norm * v / &
                    & ( denom2b(h%one,k,j,a,c) * &
                    &   denom2b(h%one,i,j,a,b) )

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

  function cross_couple(v,i,j,a,b,L) result(r)
    use MyLibrary, only: triag,sjs
    type(TwoBodyPart), intent(in) :: v
    integer, intent(in) :: i, j, a, b, L
    type(TwoBodySpace), pointer :: two
    type(Orbits), pointer :: sps
    integer :: ji, jj, ja, jb
    integer :: J2
    real(8) :: r

    r = 0.d0
    two => v%two
    sps => two%sps

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
      !r = r + dble(2*J2+1) * &
      !    & sjs(ji, jj, 2*J2, ja, jb, 2*L) * &
      !    & v%GetTwBME(i,j,a,b,J2)
      r = r + dble(2*J2+1) * &
          & sixjs%get(ji,jj,2*J2,ja,jb,2*L) * &
          & v%GetTwBME(i,j,a,b,J2)
    end do
  end function cross_couple

  subroutine MBPTCriteria(this,h)
    !
    ! MBPT criteria is not clear, but it definitely gets worse in the proton-neutron unbalance system.
    !
    class(MBPTEnergy), intent(inout) :: this
    type(Ops), intent(in) :: h
    type(MSpace), pointer :: ms
    integer :: a, b, i, j, norbs, ch, J2, n, ab, ij
    type(SingleParticleOrbit), pointer :: oa, oi
    real(8) :: max_denom1b, max_denom2b, max_val, lowest_particle, highest_hole

    ms => h%ms
    norbs = ms%sps%norbs
    lowest_particle = 100.d0
    do a = 1, norbs
      oa => ms%sps%orb(a)
      if( oa%ph /= 1 ) cycle
      lowest_particle = min(lowest_particle, h%one%GetOBME(a,a))
    end do

    highest_hole = -100.d0
    do a = 1, norbs
      oa => ms%sps%orb(a)
      if( oa%ph /= 0 ) cycle
      highest_hole = max(highest_hole, h%one%GetOBME(a,a))
    end do
    this%energy_gap = lowest_particle - highest_hole

    ! -- check HF condition, max_denom1b should be 0
    max_denom1b = 0.d0
    do a = 1, norbs
      oa => ms%sps%orb(a)
      if( oa%ph /= 1 ) cycle
      do i = 1, norbs
      oi => ms%sps%orb(i)
      if( oi%ph /= 0 ) cycle
        max_denom1b = max(max_denom1b, h%one%GetOBME(i,a) / abs(denom1b(h%one,i,a)))
      end do
    end do

    max_denom2b = 0.d0
    do ch = 1, ms%two%NChan
      J2 = ms%two%jpz(ch)%j
      n = ms%two%jpz(ch)%n_state

      max_val = 0.d0
      !$omp parallel
      !$omp do private(ab, a, b, ij, i, j) reduction(max:max_val) schedule(dynamic)
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
          max_val = max(max_val, &
              & abs(h%two%GetTwBME(i,j,a,b,J2)) / abs(denom2b(h%one,i,j,a,b)))

        end do
      end do
      !$omp end do
      !$omp end parallel
      max_denom2b = max(max_denom2b, max_val)
    end do

    write(*,'(a,f12.6)') " max(h / |e_h1 - e_p1|)               = ", max_denom1b
    write(*,'(a,f12.6)') " max(v / |e_h1 + e_h2 - e_p1 - e_p2|) = ", max_denom2b
    write(*,'(a,f12.6)') " min( e_p ) - max( e_h ) = ", this%energy_gap
    if(max(max_denom2b,max_denom1b) > 1.d0 .or. this%energy_gap < 0.d0) then
      write(*,'(a)') "Warning: MBPT might be dengerous!"
    end if
    this%perturbativity1b = max_denom1b
    this%perturbativity2b = max_denom2b
  end subroutine MBPTCriteria

  subroutine CalcScalarCorr(this,hamil,opr,is_MBPT_full)
    class(MBPTScalar), intent(inout) :: this
    type(Ops), intent(in) :: hamil, opr
    logical, intent(in) :: is_MBPT_full
    type(MSpace), pointer :: ms
    integer :: jmax
    write(*,*)
    write(*,'(a)') " Many-body perturbation calculation up to 2nd order (scalar)"
    write(*,*)

    if(.not. hamil%is_normal_ordered) then
      write(*,"(a)") "In CalcScalarCorr: "
      write(*,"(a)") " Hamiltonian has to be normal ordered"
      return
    end if

    if(.not. opr%is_normal_ordered) then
      write(*,"(a)") "In CalcScalarCorr: "
      write(*,"(a)") " Operator has to be normal ordered"
      return
    end if

    jmax = 2*hamil%ms%sps%lmax+1
    call sixjs%init(1,jmax,.true., 1,jmax,.true., 1,jmax,.true.)

    ms => hamil%ms
    this%s_0 = opr%zero
    if(is_MBPT_full) then
      call this%scalar_first(hamil,opr)
      write(*,'(a,f16.8)') "First order correction: ", this%s_1
    end if

    call this%scalar_second(hamil,opr,is_MBPT_full)
    write(*,'(a)') "Second order corrections: "
    write(*,'(a,f16.8)') "s1 p ladder  = ", this%s_2_s1p
    write(*,'(a,f16.8)') "s1 h ladder  = ", this%s_2_s1h
    write(*,'(a,f16.8)') "s1 ph bubble = ", this%s_2_s1ph
    if(is_MBPT_full) then
      write(*,'(a,f16.8)') "s2 pp ladder = ", this%s_2_s2pp
      write(*,'(a,f16.8)') "s2 hh ladder = ", this%s_2_s2hh
      write(*,'(a,f16.8)') "s2 ph ladder = ", this%s_2_s2ph
      write(*,'(a,f16.8)') "v2 pp ladder = ", this%s_2_v2pp
      write(*,'(a,f16.8)') "v2 hh ladder = ", this%s_2_v2hh
      write(*,'(a,f16.8)') "v2 ph ladder = ", this%s_2_v2ph
    end if
    write(*,'(a,f16.8)') "Total        = ", this%s_2

    call sixjs%fin()
  end subroutine CalcScalarCorr

  subroutine scalar_first(this,h,s)
    ! a, b : particle
    ! i, j : hole
    !
    !    /\===========/\
    !   /  \         /  \
    !   |  |         |  |
    ! a |  | i     b |  | j
    !   |  |         |  |
    !   \  /         \  /
    !    \/___________\/
    !
    ! 2 \sum_{i>j,a>b} <ij||ab> <ab||ij> / denominator
    use Profiler, only: timer
    class(MBPTScalar), intent(inout) :: this
    type(Ops), intent(in) :: h, s
    type(MSPace), pointer :: ms
    type(TwoBodyChannel), pointer :: ch_two
    integer :: ch, ab, ij, J2, n
    integer :: a, b, i, j
    real(8) :: vsum, v, ti

    ti = omp_get_wtime()
    ms => h%ms

    vsum = 0.d0
    do ch = 1, ms%two%NChan
      ch_two => ms%two%jpz(ch)
      J2 = ch_two%j
      n = ch_two%n_state

      v = 0.d0
      !$omp parallel
      !$omp do private(ab,a,b,ij,i,j) reduction(+:v)
      do ab = 1, ch_two%n_pp_state
        a = ch_two%n2spi1(ch_two%pps(ab))
        b = ch_two%n2spi2(ch_two%pps(ab))
        do ij = 1, ch_two%n_hh_state
          i = ch_two%n2spi1(ch_two%hhs(ij))
          j = ch_two%n2spi2(ch_two%hhs(ij))

          v = v + h%two%GetTwBME(i,j,a,b,J2) * &
              &   s%two%GetTwBME(a,b,i,j,J2) / &
              &   denom2b(h%one,i,j,a,b)
        end do
      end do
      !$omp end do
      !$omp end parallel
      vsum = vsum + v * dble(2*J2+1)
    end do

    this%s_1 = 2.d0 * vsum
    call timer%Add("First order MBPT for Scalar",omp_get_wtime()-ti)
  end subroutine scalar_first

  subroutine scalar_second(this,h,s, is_MBPT_full)
    class(MBPTScalar), intent(inout) :: this
    type(Ops), intent(in) :: h, s
    logical, intent(in) :: is_MBPT_full
    call this%scalar_second_s1p( h,s)
    call this%scalar_second_s1h( h,s)
    call this%scalar_second_s1ph(h,s)
    if(is_MBPT_full) then
      call this%scalar_second_s2pp(h,s)
      call this%scalar_second_s2hh(h,s)
      call this%scalar_second_s2ph(h,s)
      call this%scalar_second_v2pp(h,s)
      call this%scalar_second_v2hh(h,s)
      call this%scalar_second_v2ph(h,s)
    end if

    this%s_2 = this%s_2_s1h + this%s_2_s1p + this%s_2_s1ph + &
        & this%s_2_s2pp + this%s_2_s2hh + this%s_2_s2ph + &
        & this%s_2_v2pp + this%s_2_v2hh + this%s_2_v2ph
  end subroutine scalar_second

  subroutine scalar_second_s1p(this,h,s)
    ! a, b, c : particle
    ! i, j    : hole
    !     ____________
    !    /\           /\
    !   /  \         /  \
    !   |  |       c |  | j
    ! a |  | i       |=======x
    !   |  |       b |  |
    !   \  /         \  /
    !    \/___________\/
    !
    ! <ab||ij> <ij||ac> <b|x|c> / 2 denominator
    use Profiler, only: timer
    class(MBPTScalar), intent(inout) :: this
    type(Ops), intent(in) :: h, s
    type(MSPace), pointer :: ms
    integer :: ia, ib, ic, ii, ij
    integer :: a, b, c, i, j
    integer :: ja, jb, ji, jj, J2
    real(8) :: vsum, v, norm, ti

    ti = omp_get_wtime()
    ms => h%ms
    vsum = 0.d0
    !$omp parallel
    !$omp do private(ia,ib,a,b,ii,ij,i,j,ic,c,ja,jb,ji,jj,norm,v,J2) reduction(+:vsum)
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
              c = ms%particles(ic)
              if(Eket(ms%sps,a,c) > ms%two%e2max) cycle
              if(ms%sps%orb(b)%j /= ms%sps%orb(c)%j) cycle
              if(ms%sps%orb(b)%l /= ms%sps%orb(c)%l) cycle
              if(ms%sps%orb(b)%z /= ms%sps%orb(c)%z) cycle

              ja = ms%sps%orb(a)%j
              jb = ms%sps%orb(b)%j
              ji = ms%sps%orb(i)%j
              jj = ms%sps%orb(j)%j

              norm = 1.d0
              if(a==b) norm = norm * sqrt(2.d0)
              if(a==c) norm = norm * sqrt(2.d0)
              if(i==j) norm = norm * 2.d0

              v = 0.d0
              do J2 = max(abs(ja-jb), abs(ji-jj))/2, min(ja+jb, ji+jj)/2
                v = v + dble(2*J2+1) * &
                    & h%two%GetTwBME(a,b,i,j,J2) * &
                    & h%two%GetTwBME(i,j,a,c,J2) * &
                    & s%one%GetOBME(b,c)
              end do

              vsum = vsum + norm * v / &
                  & ( denom2b(h%one,i,j,a,b) * &
                  &   denom2b(h%one,i,j,a,c) )
            end do
          end do
        end do
      end do
    end do
    !$omp end do
    !$omp end parallel
    this%s_2_s1p = vsum * 0.5d0
    call timer%Add("Second order MBPT s1 p ladder",omp_get_wtime()-ti)
  end subroutine scalar_second_s1p

  subroutine scalar_second_s1h(this,h,s)
    ! a, b    : particle
    ! i, j, k : hole
    !     ____________
    !    /\           /\
    !   /  \         /  \
    !   |  |         |  | k
    ! a |  | i     b |  |====x
    !   |  |         |  | j
    !   \  /         \  /
    !    \/___________\/
    !
    ! - <ab||ij> <ik||ac> <j|x|k> / 2 denominator
    use Profiler, only: timer
    class(MBPTScalar), intent(inout) :: this
    type(Ops), intent(in) :: h, s
    type(MSpace), pointer :: ms
    integer :: ia, ib, ii, ij, ik
    integer :: a, b, i, j, k
    integer :: ja, jb, ji, jj, J2
    real(8) :: vsum, v, norm, ti

    ti = omp_get_wtime()
    ms => h%ms

    vsum = 0.d0
    !$omp parallel
    !$omp do private(ia,ib,a,b,ii,ij,i,j,ik,k,ja,jb,ji,jj,norm,v,J2) reduction(+:vsum)
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

            do ik = 1, size(ms%holes)
              k = ms%holes(ik)
              if(Eket(ms%sps,i,k) > ms%two%e2max) cycle
              if(ms%sps%orb(j)%j /= ms%sps%orb(k)%j) cycle
              if(ms%sps%orb(j)%l /= ms%sps%orb(k)%l) cycle
              if(ms%sps%orb(j)%z /= ms%sps%orb(k)%z) cycle

              ja = ms%sps%orb(a)%j
              jb = ms%sps%orb(b)%j
              ji = ms%sps%orb(i)%j
              jj = ms%sps%orb(j)%j

              norm = 1.d0
              if(i==j) norm = norm * sqrt(2.d0)
              if(i==k) norm = norm * sqrt(2.d0)
              if(a==b) norm = norm * 2.d0

              v = 0.d0
              do J2 = max(abs(ja-jb), abs(ji-jj))/2, min(ja+jb, ji+jj)/2
                v = v + dble(2*J2+1) * &
                    & h%two%GetTwBME(a,b,i,j,J2) * &
                    & h%two%GetTwBME(i,k,a,b,J2) * &
                    & s%one%GetOBME(j,k)
              end do

              vsum = vsum - norm * v / &
                  & ( denom2b(h%one,i,j,a,b) * &
                  &   denom2b(h%one,i,k,a,b) )
            end do
          end do
        end do
      end do
    end do
    !$omp end do
    !$omp end parallel
    this%s_2_s1h = vsum * 0.5d0
    call timer%Add("Second order MBPT s1 h ladder",omp_get_wtime()-ti)
  end subroutine scalar_second_s1h

  subroutine scalar_second_s1ph(this,h,s)
    ! a, b, c : particle
    ! i, j, k : hole
    !
    !    /\============x              /\=============x
    !   /  \ k                       /  \
    !   |  |___________            c |__|____________
    ! a |  |          /\   x 2 +     |  |           /\    x 2
    !   |  |         /  \          a |  |          /  \
    !   \  / i     b \  / j          \  / i      b \  / j
    !    \/___________\/              \/____________\/
    !
    ! - <ab||ij> <ij||kb> <a|x|k> / denominator
    ! + <ab||ij> <cj||ab> <c|x|i> / denominator
    use Profiler, only: timer
    class(MBPTScalar), intent(inout) :: this
    type(Ops), intent(in) :: h, s
    type(MSPace), pointer :: ms
    integer :: ia, ib, ic, ii, ij, ik
    integer :: a, b, c, i, j, k
    integer :: ja, jb, ji, jj, J2
    real(8) :: vsum, v, norm, ti

    ti = omp_get_wtime()
    ms => h%ms
    vsum = 0.d0
    !$omp parallel
    !$omp do private(ia,ib,a,b,ii,ij,i,j,ik,k,ja,jb,ji,jj,norm,v,J2,ic,c) reduction(+:vsum)
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

            do ik = 1, size(ms%holes)
              k = ms%holes(ik)
              if(Eket(ms%sps,k,b) > ms%two%e2max) cycle
              if(ms%sps%orb(a)%j /= ms%sps%orb(k)%j) cycle
              if(ms%sps%orb(a)%l /= ms%sps%orb(k)%l) cycle
              if(ms%sps%orb(a)%z /= ms%sps%orb(k)%z) cycle

              ja = ms%sps%orb(a)%j
              jb = ms%sps%orb(b)%j
              ji = ms%sps%orb(i)%j
              jj = ms%sps%orb(j)%j

              norm = 1.d0
              if(i==j) norm = norm * 2.d0
              if(a==b) norm = norm * sqrt(2.d0)

              v = 0.d0
              do J2 = max(abs(ja-jb), abs(ji-jj))/2, min(ja+jb, ji+jj)/2
                v = v + dble(2*J2+1) * &
                    & h%two%GetTwBME(a,b,i,j,J2) * &
                    & h%two%GetTwBME(i,j,k,b,J2) * &
                    & s%one%GetOBME(a,k)
              end do

              vsum = vsum - norm * v / &
                  & ( denom2b(h%one,i,j,a,b) * &
                  &   denom1b(h%one,k,a) )
            end do

            do ic = 1, size(ms%particles)
              c = ms%particles(ic)
              if(Eket(ms%sps,c,j) > ms%two%e2max) cycle
              if(ms%sps%orb(i)%j /= ms%sps%orb(c)%j) cycle
              if(ms%sps%orb(i)%l /= ms%sps%orb(c)%l) cycle
              if(ms%sps%orb(i)%z /= ms%sps%orb(c)%z) cycle

              ja = ms%sps%orb(a)%j
              jb = ms%sps%orb(b)%j
              ji = ms%sps%orb(i)%j
              jj = ms%sps%orb(j)%j

              norm = 1.d0
              if(i==j) norm = norm * sqrt(2.d0)
              if(a==b) norm = norm * 2.d0

              v = 0.d0
              do J2 = max(abs(ja-jb), abs(ji-jj))/2, min(ja+jb, ji+jj)/2
                v = v + dble(2*J2+1) * &
                    & h%two%GetTwBME(a,b,i,j,J2) * &
                    & h%two%GetTwBME(c,j,a,b,J2) * &
                    & s%one%GetOBME(c,i)
              end do

              vsum = vsum + norm * v / &
                  & ( denom2b(h%one,i,j,a,b) * &
                  &   denom1b(h%one,i,c) )
            end do


          end do
        end do
      end do
    end do
    !$omp end do
    !$omp end parallel
    this%s_2_s1ph = vsum
    call timer%Add("Second order MBPT s1 ph bubble",omp_get_wtime()-ti)
  end subroutine scalar_second_s1ph

  subroutine scalar_second_s2pp(this,h,s)
    ! a, b, c, d: particle
    ! i, j      : hole
    !     _____________
    !    /\           /\
    !   /  \         /  \
    !   |  | a     b |  |
    ! i |  |=========|  | j
    !   |  |         |  |
    !   |  | c     d |  |
    !   \  /         \  /
    !    \/___________\/
    !
    use Profiler, only: timer
    class(MBPTScalar), intent(inout) :: this
    type(Ops), intent(in) :: h, s
    type(MSpace), pointer :: ms
    type(TwoBodyChannel), pointer :: ch_two
    integer :: ch, J2, n, ab, cd, ij
    integer :: i, j, a, b, c, d
    real(8) :: v, vsum, ti

    ti = omp_get_wtime()
    ms => h%ms
    vsum = 0.d0
    do ch = 1, ms%two%NChan
      J2 = ms%two%jpz(ch)%j
      n = ms%two%jpz(ch)%n_state

      ch_two => ms%two%jpz(ch)
      if(ch_two%n_hh_state < 1) cycle
      if(ch_two%n_pp_state < 1) cycle

      v = 0.d0
      !$omp parallel
      !$omp do private(a,b,cd,c,d,ij,i,j) reduction(+:v)
      do ab = 1, ch_two%n_pp_state
        a = ch_two%n2spi1(ch_two%pps(ab))
        b = ch_two%n2spi2(ch_two%pps(ab))
        do cd = 1, ch_two%n_pp_state
          c = ch_two%n2spi1(ch_two%pps(cd))
          d = ch_two%n2spi2(ch_two%pps(cd))
          do ij = 1, ch_two%n_hh_state
            i = ch_two%n2spi1(ch_two%hhs(ij))
            j = ch_two%n2spi2(ch_two%hhs(ij))

            v = v + h%two%GetTwBME(i,j,a,b,J2) * &
                & s%two%GetTwBME(a,b,c,d,J2) * &
                & h%two%GetTwBME(c,d,i,j,J2) / &
                & ( denom2b(h%one,i,j,a,b) * &
                &   denom2b(h%one,i,j,c,d) )
          end do
        end do
      end do
      !$omp end do
      !$omp end parallel
      vsum = vsum + v * dble(2*J2+1)
    end do
    this%s_2_s2pp = vsum
    call timer%Add("Second order MBPT s2 pp ladder",omp_get_wtime()-ti)
  end subroutine scalar_second_s2pp

  subroutine scalar_second_s2hh(this,h,s)
    ! a, b      : particle
    ! i, j, k, l: hole
    !     _____________
    !    /\           /\
    !   /  \         /  \
    !   |  | i     j |  |
    ! a |  |=========|  | b
    !   |  |         |  |
    !   |  | k     l |  |
    !   \  /         \  /
    !    \/___________\/
    !
    use Profiler, only: timer
    class(MBPTScalar), intent(inout) :: this
    type(Ops), intent(in) :: h, s
    type(MSpace), pointer :: ms
    type(TwoBodyChannel), pointer :: ch_two
    integer :: ch, J2, n, ab, ij, kl
    integer :: a, b, i, j, k, l
    real(8) :: v, vsum, ti

    ti = omp_get_wtime()
    ms => h%ms

    vsum = 0.d0
    do ch = 1, ms%two%NChan
      J2 = ms%two%jpz(ch)%j
      n = ms%two%jpz(ch)%n_state

      ch_two => ms%two%jpz(ch)
      if(ch_two%n_hh_state < 1) cycle
      if(ch_two%n_pp_state < 1) cycle

      v = 0.d0
      !$omp parallel
      !$omp do private(a,b,ij,i,j,kl,k,l) reduction(+:v)
      do ab = 1, ch_two%n_pp_state
        a = ch_two%n2spi1(ch_two%pps(ab))
        b = ch_two%n2spi2(ch_two%pps(ab))
        do ij = 1, ch_two%n_hh_state
          i = ch_two%n2spi1(ch_two%hhs(ij))
          j = ch_two%n2spi2(ch_two%hhs(ij))
          do kl = 1, ch_two%n_hh_state
            k = ch_two%n2spi1(ch_two%hhs(kl))
            l = ch_two%n2spi2(ch_two%hhs(kl))

            v = v + h%two%GetTwBME(a,b,i,j,J2) * &
                & s%two%GetTwBME(i,j,k,l,J2) * &
                & h%two%GetTwBME(k,l,a,b,J2) / &
                & ( denom2b(h%one,i,j,a,b) * &
                &   denom2b(h%one,k,l,a,b) )
          end do
        end do
      end do
      !$omp end do
      !$omp end parallel

      vsum = vsum + dble(2*J2+1) * v
    end do

    this%s_2_s2hh = vsum
    call timer%Add("Second order MBPT s2 hh ladder",omp_get_wtime()-ti)
  end subroutine scalar_second_s2hh

  subroutine scalar_second_s2ph(this,h,s)
    ! a, b, c : particle
    ! i, j, k : hole
    !     _____________
    !    /\           /\
    !   /  \         /  \
    !   |  | i     b |  |
    ! a |  |=========|  | j
    !   |  |         |  |
    !   |  | k     c |  |
    !   \  /         \  /
    !    \/___________\/
    !
    use Profiler, only: timer
    class(MBPTScalar), intent(inout) :: this
    type(Ops), intent(in) :: h, s
    type(MSpace), pointer :: ms
    integer :: J2
    integer :: a, b, c, i, j, k
    integer :: ia, ib, ic, ii, ij, ik
    integer :: ja, jb, jc, ji, jj, jk
    real(8) :: v, vsum, norm, ti

    ti = omp_get_wtime()
    ms => h%ms

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
                    & cross_couple(h%two,i,j,a,b,J2) * &
                    & cross_couple(s%two,k,b,i,c,J2) * &
                    & cross_couple(h%two,a,c,k,j,J2)

                end do

                vsum = vsum + norm * v / &
                    & ( denom2b(h%one,k,j,a,c) * &
                    &   denom2b(h%one,i,j,a,b) )

              end do
            end do
          end do

        end do
      end do
    end do
    !$omp end do
    !$omp end parallel

    this%s_2_s2ph = - vsum
    call timer%Add("Second order MBPT s2 ph ladder",omp_get_wtime()-ti)
  end subroutine scalar_second_s2ph

  subroutine scalar_second_v2pp(this,h,s)
    ! a, b, c, d: particle
    ! i, j      : hole
    !
    !    /\===========/\
    !   /  \         /  \
    !   |  | a     b |  |
    ! i |  |_________|  | j
    !   |  |         |  |    x 2
    !   |  | c     d |  |
    !   \  /         \  /
    !    \/___________\/
    !
    use Profiler, only: timer
    class(MBPTScalar), intent(inout) :: this
    type(Ops), intent(in) :: h, s
    type(MSpace), pointer :: ms
    type(TwoBodyChannel), pointer :: ch_two
    integer :: ch, J2, n, ab, cd, ij
    integer :: i, j, a, b, c, d
    real(8) :: v, vsum, ti

    ti = omp_get_wtime()
    ms => h%ms
    vsum = 0.d0
    do ch = 1, ms%two%NChan
      J2 = ms%two%jpz(ch)%j
      n = ms%two%jpz(ch)%n_state

      ch_two => ms%two%jpz(ch)
      if(ch_two%n_hh_state < 1) cycle
      if(ch_two%n_pp_state < 1) cycle

      v = 0.d0
      !$omp parallel
      !$omp do private(a,b,cd,c,d,ij,i,j) reduction(+:v)
      do ab = 1, ch_two%n_pp_state
        a = ch_two%n2spi1(ch_two%pps(ab))
        b = ch_two%n2spi2(ch_two%pps(ab))
        do cd = 1, ch_two%n_pp_state
          c = ch_two%n2spi1(ch_two%pps(cd))
          d = ch_two%n2spi2(ch_two%pps(cd))
          do ij = 1, ch_two%n_hh_state
            i = ch_two%n2spi1(ch_two%hhs(ij))
            j = ch_two%n2spi2(ch_two%hhs(ij))

            v = v + s%two%GetTwBME(i,j,a,b,J2) * &
                & h%two%GetTwBME(a,b,c,d,J2) * &
                & h%two%GetTwBME(c,d,i,j,J2) / &
                & ( denom2b(h%one,i,j,a,b) * &
                &   denom2b(h%one,i,j,c,d) )
          end do
        end do
      end do
      !$omp end do
      !$omp end parallel
      vsum = vsum + v * dble(2*J2+1)
    end do
    this%s_2_v2pp = 2.d0 * vsum
    call timer%Add("Second order MBPT v2 pp ladder",omp_get_wtime()-ti)
  end subroutine scalar_second_v2pp

  subroutine scalar_second_v2hh(this,h,s)
    ! a, b      : particle
    ! i, j, k, l: hole
    !
    !    /\===========/\
    !   /  \         /  \
    !   |  | i     j |  |
    ! a |  |_________|  | b
    !   |  |         |  |   x 2
    !   |  | k     l |  |
    !   \  /         \  /
    !    \/___________\/
    !
    use Profiler, only: timer
    class(MBPTScalar), intent(inout) :: this
    type(Ops), intent(in) :: h, s
    type(MSpace), pointer :: ms
    type(TwoBodyChannel), pointer :: ch_two
    integer :: ch, J2, n, ab, ij, kl
    integer :: a, b, i, j, k, l
    real(8) :: v, vsum, ti

    ti = omp_get_wtime()
    ms => h%ms

    vsum = 0.d0
    do ch = 1, ms%two%NChan
      J2 = ms%two%jpz(ch)%j
      n = ms%two%jpz(ch)%n_state

      ch_two => ms%two%jpz(ch)
      if(ch_two%n_hh_state < 1) cycle
      if(ch_two%n_pp_state < 1) cycle

      v = 0.d0
      !$omp parallel
      !$omp do private(a,b,ij,i,j,kl,k,l) reduction(+:v)
      do ab = 1, ch_two%n_pp_state
        a = ch_two%n2spi1(ch_two%pps(ab))
        b = ch_two%n2spi2(ch_two%pps(ab))
        do ij = 1, ch_two%n_hh_state
          i = ch_two%n2spi1(ch_two%hhs(ij))
          j = ch_two%n2spi2(ch_two%hhs(ij))
          do kl = 1, ch_two%n_hh_state
            k = ch_two%n2spi1(ch_two%hhs(kl))
            l = ch_two%n2spi2(ch_two%hhs(kl))

            v = v + h%two%GetTwBME(a,b,i,j,J2) * &
                & s%two%GetTwBME(i,j,k,l,J2) * &
                & h%two%GetTwBME(k,l,a,b,J2) / &
                & ( denom2b(h%one,i,j,a,b) * &
                &   denom2b(h%one,k,l,a,b) )
          end do
        end do
      end do
      !$omp end do
      !$omp end parallel

      vsum = vsum + dble(2*J2+1) * v
    end do

    this%s_2_v2hh = 2.d0 * vsum
    call timer%Add("Second order MBPT v2 hh ladder",omp_get_wtime()-ti)
  end subroutine scalar_second_v2hh

  subroutine scalar_second_v2ph(this,h,s)
    ! a, b, c : particle
    ! i, j, k : hole
    !
    !    /\===========/\
    !   /  \         /  \
    !   |  | i     b |  |
    ! a |  |_________|  | j
    !   |  |         |  |   x 2
    !   |  | k     c |  |
    !   \  /         \  /
    !    \/___________\/
    !
    use Profiler, only: timer
    class(MBPTScalar), intent(inout) :: this
    type(Ops), intent(in) :: h, s
    type(MSpace), pointer :: ms
    integer :: J2
    integer :: a, b, c, i, j, k
    integer :: ia, ib, ic, ii, ij, ik
    integer :: ja, jb, jc, ji, jj, jk
    real(8) :: v, vsum, norm, ti

    ti = omp_get_wtime()
    ms => h%ms

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
                    & cross_couple(s%two,i,j,a,b,J2) * &
                    & cross_couple(h%two,k,b,i,c,J2) * &
                    & cross_couple(h%two,a,c,k,j,J2)

                end do

                vsum = vsum + norm * v / &
                    & ( denom2b(h%one,k,j,a,c) * &
                    &   denom2b(h%one,i,j,a,b) )

              end do
            end do
          end do

        end do
      end do
    end do
    !$omp end do
    !$omp end parallel

    this%s_2_v2ph = - 2.d0 * vsum
    call timer%Add("Second order MBPT v2 ph ladder",omp_get_wtime()-ti)
  end subroutine scalar_second_v2ph

  subroutine FinMBPTDMat(this)
    class(MBPTDMat), intent(inout) :: this
    call this%rho%fin()
    call this%Occ%fin()
    call this%C_HO2HF%fin()
    call this%C_HF2NAT%fin()
    call this%C_HO2NAT%fin()
    this%ms => null()
  end subroutine FinMBPTDMat

  subroutine InitMBPTDMat(this, HF, hamil)
    class(MBPTDMat), intent(inout) :: this
    type(HFSolver), intent(in) :: HF
    type(Ops), intent(in) :: hamil
    type(MSPace), pointer :: ms
    ms => hamil%ms
    this%ms => hamil%ms
    write(*,"(a)") "# Calculating density matrix w/ MBPT"
    call this%rho%init(ms%one, .true., 'DenMat',  0, 1, 0)
    call this%Occ%init(ms%one, .true., 'Occupation',  0, 1, 0)
    call this%C_HF2NAT%init(ms%one, .true., 'UT',  0, 1, 0)
    call this%C_HO2NAT%init(ms%one, .true., 'UT',  0, 1, 0)

    this%C_HO2HF = HF%C
    call this%density_matrix_pp(hamil)
    call this%density_matrix_hh(hamil)
    call this%density_matrix_ph(hamil)
    call this%GetCoef()
  end subroutine InitMBPTDMat

  subroutine density_matrix_pp(this, hamil)
    use Profiler, only: timer
    class(MBPTDMat), intent(inout) :: this
    type(Ops), intent(in) :: hamil
    type(MSPace), pointer :: ms
    type(Orbits), pointer :: sps
    integer :: a, b, c, i, j
    integer :: cnt, Nab, Jmin, Jmax, JJ
    type(SingleParticleOrbit), pointer :: oa, ob, oc, oi, oj
    real(8) :: r, e_acij, e_bcij, tbme, norm_ac, norm_ij, norm_bc
    integer, allocatable :: temp_a(:), temp_b(:)
    real(8) :: ti

    ti = omp_get_wtime()
    ms => hamil%ms
    sps => ms%sps

    cnt = 0
    do a = 1, sps%norbs
      oa => sps%GetOrbit(a)
      if(abs(oa%occ) > 1.d-6) cycle
      do b = a, sps%norbs
        ob => sps%GetOrbit(b)
        if(abs(ob%occ) > 1.d-6) cycle
        if(oa%l /= ob%l) cycle
        if(oa%j /= ob%j) cycle
        if(oa%z /= ob%z) cycle
        cnt = cnt + 1
      end do
    end do
    Nab = cnt
    if(Nab == 0) return
    allocate(temp_a(Nab), temp_b(Nab))

    cnt = 0
    do a = 1, sps%norbs
      oa => sps%GetOrbit(a)
      if(abs(oa%occ) > 1.d-6) cycle
      do b = a, sps%norbs
        ob => sps%GetOrbit(b)
        if(abs(ob%occ) > 1.d-6) cycle
        if(oa%l /= ob%l) cycle
        if(oa%j /= ob%j) cycle
        if(oa%z /= ob%z) cycle
        cnt = cnt + 1
        temp_a(cnt) = a
        temp_b(cnt) = b
      end do
    end do


    !$omp parallel
    !$omp do private(cnt, a, b, oa, ob, r, c, oc, i, oi, j, oj, &
    !$omp &          e_acij, e_bcij, Jmin, Jmax, tbme, JJ, norm_ac, norm_bc, norm_ij)
    do cnt = 1, Nab
      a = temp_a(cnt)
      b = temp_b(cnt)
      oa => sps%GetOrbit(a)
      ob => sps%GetOrbit(b)

      r = 0.d0
      do c = 1, sps%norbs
        oc => sps%GetOrbit(c)
        if(abs(oc%occ) > 1.d-6) cycle
        do i = 1, sps%norbs
          oi => sps%GetOrbit(i)
          if(abs(oi%occ) < 1.d-6) cycle
          do j = 1, sps%norbs
            oj => sps%GetOrbit(j)
            if(abs(oj%occ) < 1.d-6) cycle

            if((-1)**(oa%l+oc%l+oi%l+oj%l) == -1) cycle
            if(oa%z+oc%z /= oi%z+oj%z) cycle

            e_acij = denom2b(hamil%one,i,j,a,c)
            e_bcij = denom2b(hamil%one,i,j,b,c)
            if(abs(e_acij) < 1.d-8) cycle
            if(abs(e_bcij) < 1.d-8) cycle

            norm_ac = 1.d0
            norm_bc = 1.d0
            norm_ij = 1.d0
            if(a==c) norm_ac = sqrt(2.d0)
            if(b==c) norm_bc = sqrt(2.d0)
            if(i==j) norm_ij = 2.d0

            Jmin = max(abs(oa%j-oc%j), abs(oi%j-oj%j))/2
            Jmax = min(    oa%j+oc%j ,     oi%j+oj%j )/2

            tbme = 0.d0
            do JJ = Jmin, Jmax
              tbme = tbme + dble(2*JJ+1) * &
                  & hamil%two%GetTwBME(a, c, i, j, JJ) * &
                  & hamil%two%GetTwBME(i, j, b, c, JJ)
            end do
            r = r + tbme * norm_ac * norm_bc * norm_ij / (e_acij * e_bcij)
          end do
        end do
      end do
      r = 0.5d0 * r / dble(oa%j+1)
      call this%rho%SetOBME(a,b,r)
    end do
    !$omp end do
    !$omp end parallel
    deallocate(temp_a, temp_b)
    call timer%Add("density_matrix_pp",omp_get_wtime()-ti)
  end subroutine density_matrix_pp

  subroutine density_matrix_hh(this, hamil)
    use Profiler, only: timer
    class(MBPTDMat), intent(inout) :: this
    type(Ops), intent(in) :: hamil
    type(MSPace), pointer :: ms
    type(Orbits), pointer :: sps
    integer :: a, b, i, j, k
    integer :: cnt, Nij, Jmin, Jmax, JJ
    type(SingleParticleOrbit), pointer :: oa, ob, oi, oj, ok
    real(8) :: r, e_abik, e_abjk, tbme, norm_ik, norm_jk, norm_ab
    integer, allocatable :: temp_i(:), temp_j(:)
    real(8) :: ti

    ti = omp_get_wtime()
    ms => hamil%ms
    sps => ms%sps

    cnt = 0
    do i = 1, sps%norbs
      oi => sps%GetOrbit(i)
      if(abs(oi%occ) < 1.d-6) cycle
      do j = i, sps%norbs
        oj => sps%GetOrbit(j)
        if(abs(oj%occ) < 1.d-6) cycle
        if(oi%l /= oj%l) cycle
        if(oi%j /= oj%j) cycle
        if(oi%z /= oj%z) cycle
        cnt = cnt + 1
      end do
    end do
    Nij = cnt
    if(Nij == 0) return
    allocate(temp_i(Nij), temp_j(Nij))

    cnt = 0
    do i = 1, sps%norbs
      oi => sps%GetOrbit(i)
      if(abs(oi%occ) < 1.d-6) cycle
      do j = i, sps%norbs
        oj => sps%GetOrbit(j)
        if(abs(oj%occ) < 1.d-6) cycle
        if(oi%l /= oj%l) cycle
        if(oi%j /= oj%j) cycle
        if(oi%z /= oj%z) cycle
        cnt = cnt + 1
        temp_i(cnt) = i
        temp_j(cnt) = j
      end do
    end do


    !$omp parallel
    !$omp do private(cnt, i, j, oi, oj, r, a, oa, b, ob, k, ok, &
    !$omp &          e_abik, e_abjk, Jmin, Jmax, tbme, JJ, norm_ab, norm_ik, norm_jk)
    do cnt = 1, Nij
      i = temp_i(cnt)
      j = temp_j(cnt)
      oi => sps%GetOrbit(i)
      oj => sps%GetOrbit(j)

      r = 0.d0
      do a = 1, sps%norbs
        oa => sps%GetOrbit(a)
        if(abs(oa%occ) > 1.d-6) cycle

        do b = 1, sps%norbs
          ob => sps%GetOrbit(b)
          if(abs(ob%occ) > 1.d-6) cycle

          do k = 1, sps%norbs
            ok => sps%GetOrbit(k)
            if(abs(ok%occ) < 1.d-6) cycle

            if((-1)**(oa%l+ob%l+oi%l+ok%l) == -1) cycle
            if(oa%z + ob%z /= oi%z + ok%z) cycle

            e_abik = denom2b(hamil%one,i,k,a,b)
            e_abjk = denom2b(hamil%one,j,k,a,b)
            if(abs(e_abik) < 1.d-8) cycle
            if(abs(e_abjk) < 1.d-8) cycle

            Jmin = max(abs(oa%j-ob%j), abs(oi%j-ok%j))/2
            Jmax = min(    oa%j+ob%j ,     oi%j+ok%j )/2

            norm_ab = 1.d0
            norm_ik = 1.d0
            norm_jk = 1.d0
            if(a==b) norm_ab = 2.d0
            if(i==k) norm_ik = sqrt(2.d0)
            if(j==k) norm_jk = sqrt(2.d0)

            tbme = 0.d0
            do JJ = Jmin, Jmax
              tbme = tbme + dble(2*JJ+1) * &
                  & hamil%two%GetTwBME(i, k, a, b, JJ) * &
                  & hamil%two%GetTwBME(j, k, a, b, JJ)
            end do
            r = r + tbme * norm_ab * norm_ik * norm_jk / (e_abik * e_abjk)
          end do
        end do
      end do
      r = - 0.5d0 * r / dble(oi%j+1)
      if(i == j) r = r + oi%occ
      call this%rho%SetOBME(i,j,r)
    end do
    !$omp end do
    !$omp end parallel
    deallocate(temp_i, temp_j)
    call timer%Add("density_matrix_hh",omp_get_wtime()-ti)
  end subroutine density_matrix_hh

  subroutine density_matrix_ph(this, hamil)
    use Profiler, only: timer
    class(MBPTDMat), intent(inout) :: this
    type(Ops), intent(in) :: hamil
    type(MSPace), pointer :: ms
    type(Orbits), pointer :: sps
    integer :: a, b, c, i, j, k
    integer :: cnt, Nai, Jmin, Jmax, JJ
    type(SingleParticleOrbit), pointer :: oa, ob, oc, oi, oj, ok
    real(8) :: r, e_ai, e_bcij, e_abjk, tbme, norm_ab, norm_bc, norm_ij, norm_jk
    integer, allocatable :: temp_a(:), temp_i(:)
    real(8) :: ti

    ti = omp_get_wtime()
    ms => hamil%ms
    sps => ms%sps

    cnt = 0
    do a = 1, sps%norbs
      oa => sps%GetOrbit(a)
      if(abs(oa%occ) > 1.d-6) cycle
      do i = 1, sps%norbs
        oi => sps%GetOrbit(i)
        if(abs(oi%occ) < 1.d-6) cycle
        if(oa%l /= oi%l) cycle
        if(oa%j /= oi%j) cycle
        if(oa%z /= oi%z) cycle
        cnt = cnt + 1
      end do
    end do
    Nai = cnt
    if(Nai == 0) return
    allocate(temp_a(Nai), temp_i(Nai))

    cnt = 0
    do a = 1, sps%norbs
      oa => sps%GetOrbit(a)
      if(abs(oa%occ) > 1.d-6) cycle
      do i = 1, sps%norbs
        oi => sps%GetOrbit(i)
        if(abs(oi%occ) < 1.d-6) cycle
        if(oa%l /= oi%l) cycle
        if(oa%j /= oi%j) cycle
        if(oa%z /= oi%z) cycle
        cnt = cnt + 1
        temp_a(cnt) = a
        temp_i(cnt) = i
      end do
    end do

    !$omp parallel
    !$omp do private(cnt, a, i, oa, oi, r, b, ob, c, oc, j, oj, &
    !$omp &          e_ai, e_bcij, Jmin, Jmax, tbme, JJ, k, ok, &
    !$omp &          e_abjk, norm_bc, norm_ij, norm_ab, norm_jk)
    do cnt = 1, Nai
      a = temp_a(cnt)
      i = temp_i(cnt)
      oa => sps%GetOrbit(a)
      oi => sps%GetOrbit(i)

      r = 0.d0
      do b = 1, sps%norbs
        ob => sps%GetOrbit(b)
        if(abs(ob%occ) > 1.d-6) cycle
        do c = 1, sps%norbs
          oc => sps%GetOrbit(c)
          if(abs(oc%occ) > 1.d-6) cycle
          do j = 1, sps%norbs
            oj => sps%GetOrbit(j)
            if(abs(oj%occ) < 1.d-6) cycle
            if((-1)**(oa%l+oj%l+ob%l+oc%l) == -1) cycle
            if(oa%z+oj%z /= ob%z+oc%z) cycle
            e_ai = denom1b(hamil%one, i, a)
            e_bcij = denom2b(hamil%one, i, j, b, c)
            if(abs(e_ai) < 1.d-8) cycle
            if(abs(e_bcij) < 1.d-8) cycle

            Jmin = max(abs(oa%j-oj%j), abs(ob%j-oc%j))/2
            Jmax = min(    oa%j+oj%j ,     ob%j+oc%j )/2

            norm_bc = 1.d0
            norm_ij = 1.d0
            if(b==c) norm_bc = 2.d0
            if(i==j) norm_ij = sqrt(2.d0)
            tbme = 0.d0
            do JJ = Jmin, Jmax
              tbme = tbme + dble(2*JJ+1) * &
                  &  hamil%two%GetTwBME(a, j, b, c, JJ) * &
                  &  hamil%two%GetTwBME(b, c, i, j, JJ)
            end do
            r = r + tbme / (e_ai * e_bcij) * norm_bc * norm_ij
          end do
        end do
      end do

      do b = 1, sps%norbs
        ob => sps%GetOrbit(b)
        if(abs(ob%occ) > 1.d-6) cycle
        do j = 1, sps%norbs
          oj => sps%GetOrbit(j)
          if(abs(oj%occ) < 1.d-6) cycle
          do k = 1, sps%norbs
            ok => sps%GetOrbit(k)
            if(abs(ok%occ) < 1.d-6) cycle

            if((-1)**(oa%l+ob%l+oj%l+ok%l) == -1) cycle
            if(oa%z + ob%z /= oj%z + ok%z) cycle

            e_ai = denom1b(hamil%one, i, a)
            e_abjk = denom2b(hamil%one, j, k, a, b)
            if(abs(e_ai) < 1.d-8) cycle
            if(abs(e_abjk) < 1.d-8) cycle
            norm_ab = 1.d0
            norm_jk = 1.d0
            if(a==b) norm_ab = sqrt(2.d0)
            if(j==k) norm_jk = 2.d0
            Jmin = max(abs(oa%j-ob%j), abs(oj%j-ok%j))/2
            Jmax = min(    oa%j+ob%j ,     oj%j+ok%j )/2
            tbme = 0.d0
            do JJ = Jmin, Jmax
              tbme = tbme + dble(2*JJ+1) * &
                  &  hamil%two%GetTwBME(k, j, i, b, JJ) * &
                  &  hamil%two%GetTwBME(a, b, k, j, JJ)
            end do
            r = r - tbme / (e_ai * e_abjk) * norm_jk * norm_ab
          end do
        end do
      end do
      r = 0.5d0 * r / dble(oa%j+1)
      call this%rho%SetOBME(a,i,r)
    end do
    !$omp end do
    !$omp end parallel
    deallocate(temp_a, temp_i)
    call timer%Add("density_matrix_ph",omp_get_wtime()-ti)
  end subroutine density_matrix_ph

  subroutine GetCoef(this)
    class(MBPTDMat), intent(inout) :: this
    type(EigenSolSymD) :: sol
    integer :: ch, i, m, jj, iz
    real(8) :: N, Z
    N = 0.d0
    Z = 0.d0
    do ch = 1, this%rho%one%NChan
      jj = this%rho%one%jpz(ch)%j
      iz = this%rho%one%jpz(ch)%z
      call sol%init(this%rho%MatCh(ch,ch)%DMat)
      call sol%DiagSym(this%rho%MatCh(ch,ch)%DMat)
      m = size(sol%eig%v)
      do i = 1, m
        this%C_HF2NAT%MatCh(ch,ch)%m(i,:) = sol%vec%m(m-i+1,:)
        this%Occ%MatCh(ch,ch)%m(i,i) = sol%eig%v(m-i+1)
      end do
      this%C_HO2NAT%MatCh(ch,ch)%DMat = &
          & this%C_HO2HF%MatCh(ch,ch)%DMat * &
          & this%C_HF2NAT%MatCh(ch,ch)%DMat
      !call this%rho%MatCh(ch,ch)%DMat%prt("rho")
      !call this%Occ%MatCh(ch,ch)%prt("Occupation Number")
      if(iz == -1) Z = Z + sum(sol%eig%v) * dble(jj+1)
      if(iz ==  1) N = N + sum(sol%eig%v) * dble(jj+1)
      call sol%fin()
    end do
    write(*,"(a,i4,a,f6.2)") " Actual Z: ", this%ms%Z, ", Z from tr(rho): ", Z
    write(*,"(a,i4,a,f6.2)") " Actual N: ", this%ms%N, ", N from tr(rho): ", N
  end subroutine GetCoef

end module HFMBPT
