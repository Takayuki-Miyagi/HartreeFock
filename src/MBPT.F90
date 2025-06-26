!
! This is for a general MBPT case (F is not diagonal)
!
module gMBPT
  use omp_lib
  use myfort
  use ModelSpace
  use Operators
  implicit none

  ! Methods for general
  private :: get_denominator1
  private :: get_denominator2
  private :: get_denominator3
  private :: get_denominator4
  private :: get_pphh_part
  private :: get_xc_pphh_to_phph

  type :: gMBPTEnergy
    real(8) :: e_0 = 0.d0
    real(8) :: e_2(2) = 0.d0
    real(8) :: e_3(14) = 0.d0
  contains
    procedure :: calc => CalcEnergyCorr
  end type gMBPTEnergy

  type(SixJsStore), private :: sixjs
  logical, private :: EN_denominator=.false.  ! Epstein-Nesbet denominator (Moller-Plesset is default)

contains

  subroutine CalcEnergyCorr(this, hamil, EN_denominator_in)
    class(gMBPTEnergy), intent(inout) :: this
    type(Ops), intent(in) :: hamil
    logical, intent(in), optional :: EN_denominator_in
    real(8) :: s_1
    integer :: i, jmax
    if(present(EN_denominator_in)) EN_denominator = EN_denominator_in
    if(.not. EN_denominator) write(*,"(a)") "# Moller-Plesset (MP) denominator"
    if(EN_denominator)       write(*,"(a)") "# Epstein-Nesbet (EN) denominator"
    write(*,*)
    if(.not. hamil%is_normal_ordered) then
      write(*,"(a)") "In CalcEnergyCorr: "
      write(*,"(a)") " Hamiltonian has to be normal ordered"
      return
    end if
    jmax = 2*hamil%ms%sps%lmax+1
    call sixjs%init(1,jmax,.true., 1,jmax,.true., 1,jmax,.true.)
    this%e_0 = hamil%zero
    write(*,'(a,2f16.8)') "Reference energy: ", this%e_0

    ! Second order
    this%e_2(1) = energy_second_1(hamil)
    this%e_2(2) = energy_second_2(hamil)
    write(*,'(a,2f16.8)') "Second order corrections: ", this%e_2
    write(*,'(a,f16.8)') "Total 2nd order correction: ", sum(this%e_2)

    ! Third order
    write(*,*)
    this%e_3( 1) = energy_third_1(hamil)
    this%e_3( 2) = energy_third_2(hamil)
    this%e_3( 3) = energy_third_3(hamil)
    this%e_3( 4) = energy_third_4(hamil)
    this%e_3( 5) = energy_third_5(hamil)
    this%e_3( 6) = energy_third_6(hamil)
    this%e_3( 7) = energy_third_7(hamil)
    this%e_3( 8) = this%e_3(4)
    this%e_3( 9) = this%e_3(5)
    this%e_3(10) = energy_third_10(hamil)
    this%e_3(11) = energy_third_11(hamil)
    this%e_3(12) = this%e_3(10)
    this%e_3(13) = energy_third_13(hamil)
    this%e_3(14) = energy_third_14(hamil)
    do i = 1, size(this%e_3)
      write(*,'(a,i2,f16.8)') "Third order corrections: ", i, this%e_3(i)
    end do
    write(*,'(a,f16.8)') "Total 3rd order correction: ", sum(this%e_3)
    !write(*,*) energy_third_4_M(hamil), energy_third_4(hamil)
    !write(*,*) energy_third_5_M(hamil), energy_third_5(hamil)
    !write(*,*) energy_third_6_M(hamil), energy_third_6(hamil)
    !write(*,*) energy_third_7_M(hamil), energy_third_7(hamil)
    !write(*,*) energy_third_10_M(hamil), energy_third_10(hamil)
    !write(*,*) energy_third_11_M(hamil), energy_third_11(hamil)
    !write(*,*) energy_third_13_M(hamil), energy_third_13(hamil)
    !write(*,*) energy_third_14_M(hamil), energy_third_14(hamil)

  end subroutine CalcEnergyCorr

  function energy_second_1(h) result(r)
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
    type(Ops), intent(in) :: h
    type(MSPace), pointer :: ms
    type(TwoBodyChannel), pointer :: ch_two
    type(ThreeBodyChannel), pointer :: ch_three
    type(Orbits), pointer :: sps
    integer :: ch, ab, ij, J2, n
    integer :: a, b, c, i, j, k, abc, ijk, Jab, Jij
    type(SingleParticleOrbit), pointer :: oa, ob, oc, oi, oj, ok
    real(8) :: vsum, v, ti, r

    ti = omp_get_wtime()
    ms => h%ms
    sps => ms%sps
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
        a = ch_two%n2spi1(ch_two%pp_s(ab))
        b = ch_two%n2spi2(ch_two%pp_s(ab))
        do ij = 1, ch_two%n_hh_state
          i = ch_two%n2spi1(ch_two%hh_s(ij))
          j = ch_two%n2spi2(ch_two%hh_s(ij))

          v = v + h%two%GetTwBME(i,j,a,b,J2)**2 / &
              & get_denominator2(h,i,j,a,b)
        end do
      end do
      !$omp end do
      !$omp end parallel
      vsum = vsum + v * dble(2*J2+1)
    end do
    r = vsum
    call timer%Add("Second order MBPT",omp_get_wtime()-ti)
  end function energy_second_1

  function energy_second_2(h) result(r)
    ! a : particle
    ! i : hole
    !     ____________\/
    !    /\           /\
    !   /  \
    !   |  |
    ! a |  | i
    !   |  |
    !   \  /
    !    \/___________\/
    !                 /\
    ! \sum_{i a} <i||a> <a||i> / denominator
    type(Ops), intent(in) :: h
    type(MSPace), pointer :: ms
    type(Orbits), pointer :: sps
    integer :: a, i
    type(SingleParticleOrbit), pointer :: oa, oi
    real(8) :: vsum, ti, r, denom

    ti = omp_get_wtime()
    ms => h%ms
    sps => ms%sps
    vsum = 0.d0
    do a = 1, sps%norbs
      oa => sps%GetOrbit(a)
      if(oa%GetOccupation() > 1.d-6) cycle
      do i = 1, sps%norbs
        oi => sps%GetOrbit(i)
        if(oi%GetOccupation() < 1.d-6) cycle
        denom = 1.d0 / get_denominator1(h,i,a)
        vsum = vsum + h%one%GetOBME(a,i)**2 * denom

      end do
    end do
    r = vsum
    call timer%Add("Second order MBPT",omp_get_wtime()-ti)
  end function energy_second_2

  function energy_third_1(h) result(r)
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
    type(Ops), intent(in) :: h
    type(MSPace), pointer :: ms
    type(TwoBodyChannel), pointer :: ch_two
    integer :: ch, J2, i
    real(8) :: v, vsum, ti, r
    type(DMat) :: m1, m2, m3

    ms => h%ms
    ti = omp_get_wtime()
    vsum = 0.d0
    do ch = 1, ms%two%NChan
      J2 = ms%two%jpz(ch)%j

      ch_two => ms%two%jpz(ch)
      if(ch_two%n_hh_state < 1) cycle
      if(ch_two%n_pp_state < 1) cycle

      m1 = get_pphh_part(h%two%MatCh(ch,ch), h)
      m2 = h%two%MatCh(ch,ch)%get_pppp(ms%sps)
      m3 = m1%t() * (m2 * m1)
      v = 0.d0
      do i = 1, size(m3%m,1)
        v = v + m3%m(i,i)
      end do
      vsum = vsum + v * dble(2*J2+1)
    end do
    r = vsum
    call timer%Add("Third order MBPT pp ladder",omp_get_wtime()-ti)
  end function energy_third_1

  function energy_third_2(h) result(r)
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
    type(Ops), intent(in) :: h
    type(MSPace), pointer :: ms
    type(TwoBodyChannel), pointer :: ch_two
    integer :: ch, J2, i
    real(8) :: v, vsum, ti, r
    type(DMat) :: m1, m2, m3

    ti = omp_get_wtime()

    ms => h%ms
    vsum = 0.d0
    do ch = 1, ms%two%NChan
      J2 = ms%two%jpz(ch)%j

      ch_two => ms%two%jpz(ch)
      if(ch_two%n_hh_state < 1) cycle
      if(ch_two%n_pp_state < 1) cycle

      m1 = get_pphh_part(h%two%MatCh(ch,ch), h)
      m2 = h%two%MatCh(ch,ch)%get_hhhh(ms%sps)
      m3 = m1%t() * (m1 * m2)
      v = 0.d0
      do i = 1, size(m3%m,1)
        v = v + m3%m(i,i)
      end do
      vsum = vsum + dble(2*J2+1) * v
    end do
    r = vsum
    call timer%Add("Third order MBPT hh ladder",omp_get_wtime()-ti)
  end function energy_third_2

  function energy_third_3(h) result(r)
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
    type(Ops), intent(in) :: h
    type(MSPace), pointer :: ms
    integer :: ch, J2
    type(CrossCoupledTwoBodyChannel), pointer :: ch_cc
    integer :: i
    real(8) :: v, vsum, ti, r
    type(DMat) :: m1, m2, m3

    ti = omp_get_wtime()
    ms => h%ms
    vsum = 0.d0
    do ch = 1,ms%cc_two%NChan
      ch_cc => ms%cc_two%jpz(ch)
      J2 = ch_cc%j
      m1 = get_xc_pphh_to_phph(h%two, ch_cc, h)
      if(m1%n_row<1 .or. m1%n_col<1) cycle
      m2 = h%two%get_xc_hphp2phph(ch_cc)
      if(m2%n_row<1 .or. m2%n_col<1) cycle
      m3 = m1 * (m2 * m1%t())
      v = 0.d0
      do i = 1, m3%n_row
        v = v + m3%m(i,i)
      end do
      vsum = vsum + dble(2*J2+1) * v
      call m1%fin()
      call m2%fin()
      call m3%fin()
    end do
    r = - vsum
    call timer%Add("Third order MBPT ph ladder",omp_get_wtime()-ti)
  end function energy_third_3

  function energy_third_4(h) result(r)
    ! a, b, c : particle
    ! i, j, k : hole
    !     _____________\/
    !    /\            /\
    !   /  \
    ! c |__|____________
    !   |  |           /\
    ! a |  |          /  \
    !   \  / i      b \  / j
    !    \/____________\/
    !
    ! (1/2) <ab||ij> <cj||ab> <c||i> / denominator
    type(Ops), intent(in) :: h
    type(MSPace), pointer :: ms
    integer :: ch, J2, bra, ket, iph
    type(TwoBodyChannel), pointer :: ch_two
    type(TwoBodyPartChannel) :: tmp_phhh
    type(DMat) :: m1, m2
    real(8) :: vsum, v, norm, ti, r
    integer :: a, b, c, i, j, k
    type(SingleParticleOrbit), pointer :: oa, ob, oc, oi, oj, ok

    ti = omp_get_wtime()
    ms => h%ms
    vsum = 0.d0
    do ch = 1, ms%two%NChan
      ch_two => ms%two%jpz(ch)
      J2 = ch_two%j
      if(ch_two%n_hh_state < 1) cycle
      if(ch_two%n_pp_state < 1) cycle
      call tmp_phhh%init(ch_two, ch_two)
      m1 = get_pphh_part(h%two%MatCh(ch,ch), h)
      m2 = h%two%MatCh(ch,ch)%get_phpp(ms%sps)
      if(m2%n_row < 1 .or. m2%n_col < 1) cycle
      call tmp_phhh%set_phhh(ms%sps, m2*m1)

      v = 0.d0
      do c = 1, ms%sps%norbs
        oc => ms%sps%GetOrbit(c)
        if(abs(oc%GetOccupation()) > 1.d-6) cycle
        do i = 1, ms%sps%norbs
          oi => ms%sps%GetOrbit(i)
          if(abs(oi%GetOccupation()) < 1.d-6) cycle
          if(oc%j /= oi%j) cycle
          if(oc%l /= oi%l) cycle
          if(oc%z /= oi%z) cycle
          do j = 1, ms%sps%norbs
            oj => ms%sps%GetOrbit(i)
            if(abs(oj%GetOccupation()) < 1.d-6) cycle

            bra = ch_two%spis2n(i,j)
            ket = ch_two%spis2n(c,j)
            if(bra*ket == 0) cycle
            iph = ch_two%iphase(i,j) * ch_two%iphase(c,j)
            norm = 2.d0
            if(i==j) norm = norm*sqrt(2.d0)
            v = v + norm * dble(2*J2+1) * dble(iph) * &
                & tmp_phhh%m(bra,ket) * h%one%GetOBME(c,i) / &
                & get_denominator1(h,i,c)
          end do
        end do
      end do
      vsum = vsum + v
      call tmp_phhh%fin()
      call m1%fin()
      call m2%fin()
    end do
    r = vsum * 0.5d0
  end function energy_third_4

  function energy_third_5(h) result(r)
    ! a, b, c : particle
    ! i, j, k : hole
    !     ____________\/
    !    /\           /\
    !   /  \
    !   |  | k
    ! a |  |___________
    !   |  |          /\
    !   |  | i     j /  \ b
    !   \  /         \  /
    !    \/___________\/
    !
    ! - (1/2) <ab||ij> <ij||kb> <a||k> / denominator
    type(Ops), intent(in) :: h
    type(MSPace), pointer :: ms
    integer :: ch, J2, bra, ket, iph
    type(TwoBodyChannel), pointer :: ch_two
    type(TwoBodyPartChannel) :: tmp_phpp
    type(DMat) :: m1, m2
    real(8) :: vsum, v, norm, ti, r
    integer :: a, b, c, i, j, k
    type(SingleParticleOrbit), pointer :: oa, ob, oc, oi, oj, ok

    ti = omp_get_wtime()
    ms => h%ms
    vsum = 0.d0
    do ch = 1, ms%two%NChan
      ch_two => ms%two%jpz(ch)
      J2 = ch_two%j
      if(ch_two%n_hh_state < 1) cycle
      if(ch_two%n_pp_state < 1) cycle
      call tmp_phpp%init(ch_two, ch_two)
      m1 = get_pphh_part(h%two%MatCh(ch,ch), h)
      m2 = h%two%MatCh(ch,ch)%get_phhh(ms%sps)
      if(m2%n_row < 1 .or. m2%n_col < 1) cycle
      call tmp_phpp%set_phpp(ms%sps, m2*m1%t())

      v = 0.d0
      do a = 1, ms%sps%norbs
        oa => ms%sps%GetOrbit(a)
        if(abs(oa%GetOccupation()) > 1.d-6) cycle
        do k = 1, ms%sps%norbs
          ok => ms%sps%GetOrbit(k)
          if(abs(ok%GetOccupation()) < 1.d-6) cycle
          if(oa%j /= ok%j) cycle
          if(oa%l /= ok%l) cycle
          if(oa%z /= ok%z) cycle

          do b = 1, ms%sps%norbs
            ob => ms%sps%GetOrbit(b)
            if(abs(ob%GetOccupation()) > 1.d-6) cycle

            bra = ch_two%spis2n(a,b)
            ket = ch_two%spis2n(k,b)
            if(bra*ket == 0) cycle
            iph = ch_two%iphase(a,b) * ch_two%iphase(k,b)
            norm = 2.d0
            if(a==b) norm = norm*sqrt(2.d0)
            v = v + norm * dble(2*J2+1) * dble(iph) * &
                & tmp_phpp%m(bra,ket) * h%one%GetOBME(a,k) / &
                & get_denominator1(h,k,a)
          end do
        end do
      end do
      vsum = vsum + v
      call tmp_phpp%fin()
      call m1%fin()
      call m2%fin()
    end do
    r = vsum * (-0.5d0)
  end function energy_third_5

  function energy_third_6(h) result(r)
    ! a, b, c : particle
    ! i, j    : hole
    !     ____________
    !    /\           /\
    !   /  \         /  \
    !   |  |       c |  | j
    ! a |  | i       |_______\/
    !   |  |       b |  |    /\
    !   \  /         \  /
    !    \/___________\/
    !
    ! <ab||ij> <ij||ac> <b||c> / 2 denominator
    type(Ops), intent(in) :: h
    type(MSPace), pointer :: ms
    integer :: ch, J2, bra, ket
    type(TwoBodyChannel), pointer :: ch_two
    type(TwoBodyPartChannel) :: tmp
    type(DMat) :: m1, m2
    integer :: a, b, c, iph
    type(SingleParticleOrbit), pointer :: oa, ob, oc
    real(8) :: vsum, v, ti, r, norm

    ti = omp_get_wtime()
    ms => h%ms
    vsum = 0.d0

    do ch = 1, ms%two%NChan
      ch_two => ms%two%jpz(ch)
      J2 = ch_two%j
      if(ch_two%n_hh_state < 1) cycle
      if(ch_two%n_pp_state < 1) cycle
      m1 = get_pphh_part(h%two%MatCh(ch,ch), h)
      m2 = m1 * m1%t()
      call tmp%init(ch_two, ch_two)
      call tmp%set_pppp(ms%sps, m2)

      v = 0.d0
      do b = 1, ms%sps%norbs
        ob => ms%sps%GetOrbit(b)
        if(abs(ob%GetOccupation()) > 1.d-6) cycle
        do c = 1, ms%sps%norbs
          oc => ms%sps%GetOrbit(c)
          if(abs(oc%GetOccupation()) > 1.d-6) cycle
          if(ob%j /= oc%j) cycle
          if(ob%l /= oc%l) cycle
          if(ob%z /= oc%z) cycle
          do a = 1, ms%sps%norbs
            oa => ms%sps%GetOrbit(a)
            if(abs(oa%GetOccupation()) > 1.d-6) cycle

            bra = ch_two%spis2n(a,b)
            ket = ch_two%spis2n(a,c)
            if(bra*ket == 0) cycle
            iph = ch_two%iphase(a,b) * ch_two%iphase(a,c)
            norm = 2.d0
            if(a==b) norm = norm * sqrt(2.d0)
            if(a==c) norm = norm * sqrt(2.d0)
            v = v + norm * dble(2*J2+1) * tmp%m(bra,ket) * h%one%GetOBME(b,c) * dble(iph)
          end do
        end do
      end do
      vsum = vsum + v
      call tmp%fin()
      call m1%fin()
      call m2%fin()
    end do
    r = vsum * 0.5d0
  end function energy_third_6

  function energy_third_7(h) result(r)
    ! a, b, c : particle
    ! i, j    : hole
    !     ____________
    !    /\           /\
    !   /  \         /  \
    !   |  |         |  | k
    ! a |  | i     b |  |____\/
    !   |  |         |  | j  /\
    !   \  /         \  /
    !    \/___________\/
    !
    ! - <ab||ij> <ik||ab> <j||k> / 2 denominator
    type(Ops), intent(in) :: h
    type(MSpace), pointer :: ms
    integer :: ch, J2, bra, ket
    type(TwoBodyChannel), pointer :: ch_two
    type(TwoBodyPartChannel) :: tmp
    type(DMat) :: m1, m2
    integer :: i, j, k, iph
    type(SingleParticleOrbit), pointer :: oi, oj, ok
    real(8) :: vsum, v, norm, ti, r

    ti = omp_get_wtime()
    ms => h%ms

    vsum = 0.d0
    do ch = 1, ms%two%NChan
      ch_two => ms%two%jpz(ch)
      J2 = ch_two%j
      if(ch_two%n_hh_state < 1) cycle
      if(ch_two%n_pp_state < 1) cycle
      m1 = get_pphh_part(h%two%MatCh(ch,ch), h)
      m2 = m1%t() * m1
      call tmp%init(ch_two, ch_two)
      call tmp%set_hhhh(ms%sps, m2)

      v = 0.d0
      do j = 1, ms%sps%norbs
        oj => ms%sps%GetOrbit(j)
        if(abs(oj%GetOccupation()) < 1.d-6) cycle
        do k = 1, ms%sps%norbs
          ok => ms%sps%GetOrbit(k)
          if(abs(ok%GetOccupation()) < 1.d-6) cycle
          if(oj%j /= ok%j) cycle
          if(oj%l /= ok%l) cycle
          if(oj%z /= ok%z) cycle
          do i = 1, ms%sps%norbs
            oi => ms%sps%GetOrbit(i)
            if(abs(oi%GetOccupation()) < 1.d-6) cycle

            bra = ch_two%spis2n(i,j)
            ket = ch_two%spis2n(i,k)
            if(bra*ket == 0) cycle
            iph = ch_two%iphase(i,j) * ch_two%iphase(i,k)
            norm = 2.d0
            if(i==j) norm = norm * sqrt(2.d0)
            if(i==k) norm = norm * sqrt(2.d0)
            v = v - norm * dble(2*J2+1) * tmp%m(bra,ket) * h%one%GetOBME(j,k) * dble(iph)
          end do
        end do
      end do
      vsum = vsum + v
      call tmp%fin()
      call m1%fin()
      call m2%fin()
    end do
    r = vsum * 0.5d0
  end function energy_third_7

  function energy_third_10(h) result(r)
    ! a, b, c : particle
    ! i, j    : hole
    !     ____________\/
    !    /\           /\
    !   /  \         /  \
    !   |  |           ____\/
    ! a |  | i        /\   /\
    !   |  |         /  \
    !   \  /       j \  / b
    !    \/___________\/
    !
    !  <ab||ij> <j||b> <i||a> / denominator
    type(Ops), intent(in) :: h
    type(MSpace), pointer :: ms
    integer :: i, j, a, b
    type(SingleParticleOrbit), pointer :: oi, oj, oa, ob
    integer :: ch, J2, bra, ket
    type(TwoBodyChannel), pointer :: ch_two
    type(TwoBodyPartChannel) :: tmp
    type(DMat) :: m1, m2
    real(8) :: vsum, v, norm, ti, r, denom

    ms => h%ms

    vsum = 0.d0
    do i = 1, ms%sps%norbs
      oi => ms%sps%GetOrbit(i)
      if(abs(oi%GetOccupation()) < 1.d-6) cycle
      do j = 1, ms%sps%norbs
        oj => ms%sps%GetOrbit(j)
        if(abs(oj%GetOccupation()) < 1.d-6) cycle

        do a = 1, ms%sps%norbs
          oa => ms%sps%GetOrbit(a)
          if(abs(oa%GetOccupation()) > 1.d-6) cycle
          do b = 1, ms%sps%norbs
            ob => ms%sps%GetOrbit(b)
            if(abs(ob%GetOccupation()) > 1.d-6) cycle

            denom = 1.d0 / (get_denominator2(h,i,j,a,b) * get_denominator1(h,i,a))
            if(oi%l /= oa%l) cycle
            if(oi%j /= oa%j) cycle
            if(oi%z /= oa%z) cycle
            if(oj%l /= ob%l) cycle
            if(oj%j /= ob%j) cycle
            if(oj%z /= ob%z) cycle
            norm = 1.d0
            if(a==b) norm = norm * sqrt(2.d0)
            if(i==j) norm = norm * sqrt(2.d0)
            do J2 = abs(oa%j-ob%j)/2, (oa%j+ob%j)/2
              vsum = vsum + dble(2*J2+1) * norm * h%two%GetTwBME(a,b,i,j,J2) * h%one%GetOBME(j,b) * h%one%GetOBME(i,a) * denom
            end do

          end do
        end do

      end do
    end do
    r = vsum
  end function energy_third_10

  function energy_third_11(h) result(r)
    ! a, b, c : particle
    ! i, j    : hole
    !     ____________\/
    !    /\           /\
    !   /  \
    ! b \  / i
    !    \/____________
    !                 /\
    !                /  \
    !              j \  / a
    !    \/___________\/
    !    /\
    !  <j||a> <bj||ia> <b||i> / denominator
    type(Ops), intent(in) :: h
    type(MSpace), pointer :: ms
    integer :: i, j, a, b
    type(SingleParticleOrbit), pointer :: oi, oj, oa, ob
    integer :: ch, J2, bra, ket
    type(TwoBodyChannel), pointer :: ch_two
    type(TwoBodyPartChannel) :: tmp
    type(DMat) :: m1, m2
    real(8) :: vsum, v, norm, ti, r, denom

    ms => h%ms
    vsum = 0.d0
    do i = 1, ms%sps%norbs
      oi => ms%sps%GetOrbit(i)
      if(abs(oi%GetOccupation()) < 1.d-6) cycle
      do j = 1, ms%sps%norbs
        oj => ms%sps%GetOrbit(j)
        if(abs(oj%GetOccupation()) < 1.d-6) cycle

        do a = 1, ms%sps%norbs
          oa => ms%sps%GetOrbit(a)
          if(abs(oa%GetOccupation()) > 1.d-6) cycle
          do b = 1, ms%sps%norbs
            ob => ms%sps%GetOrbit(b)
            if(abs(ob%GetOccupation()) > 1.d-6) cycle

            denom = 1.d0 / (get_denominator1(h,j,a) * get_denominator1(h,i,b))
            if(oi%l /= ob%l) cycle
            if(oi%j /= ob%j) cycle
            if(oi%z /= ob%z) cycle
            if(oj%l /= oa%l) cycle
            if(oj%j /= oa%j) cycle
            if(oj%z /= oa%z) cycle
            do J2 = abs(oa%j-ob%j)/2, (oa%j+ob%j)/2
              vsum = vsum + dble(2*J2+1) * h%two%GetTwBME(b,j,i,a,J2) * h%one%GetOBME(j,a) * h%one%GetOBME(i,b) * denom
            end do

          end do
        end do

      end do
    end do
    r = vsum
  end function energy_third_11

  function energy_third_13(h) result(r)
    ! a, b : particle
    ! i, j : hole
    !       ______\/
    !     /|      /\
    !    / |
    !   /  | b
    !   |  |
    ! i |  |______\/
    !   |  |      /\
    !   \  | a
    !    \ |
    !     \|______\/
    !             /\
    !
    type(Ops), intent(in) :: h
    type(MSPace), pointer :: ms
    type(OneBodyChannel), pointer :: ch_one
    integer :: ch, a, b, i, j, ia, ib, ii, ij, J1, n, bra, ket
    type(Orbits), pointer :: sps
    type(SingleParticleOrbit), pointer :: oa, ob, oi, oj
    real(8) :: spp, r, ti

    ti = omp_get_wtime()
    ms => h%ms
    sps => ms%sps
    spp = 0.d0
    do i = 1, ms%nh
      ii = ms%holes(i)
      oi => sps%GetOrbit(ii)
      do a = 1, ms%np
        ia = ms%particles(a)
        oa => sps%GetOrbit(ia)
        if(oi%l /= oa%l) cycle
        if(oi%j /= oa%j) cycle
        if(oi%z /= oa%z) cycle
        do b = 1, ms%np
          ib = ms%particles(b)
          ob => sps%GetOrbit(ib)
          if(oa%l /= ob%l) cycle
          if(oa%j /= ob%j) cycle
          if(oa%z /= ob%z) cycle
          spp = spp + &
              & dble(oi%j+1) * &
              &  h%one%GetOBME(ii,ia) * h%one%GetOBME(ia,ib) * h%one%GetOBME(ib,ii) / &
              & ( get_denominator1(h,ii,ia) * get_denominator1(h,ii,ib) )
        end do
      end do
    end do
    r = spp
  end function energy_third_13

  function energy_third_14(h) result(r)
    ! a, b : particle
    ! i, j : hole
    !       ______\/
    !     /|      /\
    !    / |
    !   /  | j
    !   |  |
    ! a |  |______\/
    !   |  |      /\
    !   \  | i
    !    \ |
    !     \|______\/
    !             /\
    !
    type(Ops), intent(in) :: h
    type(MSPace), pointer :: ms
    type(OneBodyChannel), pointer :: ch_one
    integer :: ch, a, b, i, j, ia, ib, ii, ij, J1, n, bra, ket
    type(Orbits), pointer :: sps
    type(SingleParticleOrbit), pointer :: oa, ob, oi, oj
    real(8) :: shh, r, ti

    ti = omp_get_wtime()
    ms => h%ms
    sps => ms%sps
    shh = 0.d0
    do i = 1, ms%nh
      ii = ms%holes(i)
      oi => sps%GetOrbit(ii)
      do j = 1, ms%nh
        ij = ms%holes(j)
        oj => sps%GetOrbit(ij)
        if(oi%l /= oj%l) cycle
        if(oi%j /= oj%j) cycle
        if(oi%z /= oj%z) cycle
        do a = 1, ms%np
          ia = ms%particles(a)
          oa => sps%GetOrbit(ia)
          if(oi%l /= oa%l) cycle
          if(oi%j /= oa%j) cycle
          if(oi%z /= oa%z) cycle
          shh = shh - &
              & dble(oi%j+1) * &
              & h%one%GetOBME(ia,ii) * h%one%GetOBME(ii,ij) * h%one%GetOBME(ij,ia) / &
              & ( get_denominator1(h,ii,ia) * get_denominator1(h,ij,ia) )
        end do
      end do
    end do
    r = shh
  end function energy_third_14

  function get_denominator1(h,h1,p1) result(r)
    ! E_{0} - <p1h1| H | p1h1 >
    type(Ops), intent(in) :: h
    integer, intent(in) :: h1, p1
    real(8) :: r

    r = h%one%GetDenominator1(h1,p1)
    if(.not. EN_denominator) return
    r = r - h%two%GetTwBMEMon(p1,h1,p1,h1)
  end function get_denominator1

  function get_denominator2(h,h1,h2,p1,p2) result(r)
    ! E_{0} - <p1p2h1h2| H | p1p2h1h2 >
    type(Ops), intent(in) :: h
    integer, intent(in) :: h1, h2, p1, p2
    real(8) :: r

    r = h%one%GetDenominator2(h1,h2,p1,p2)
    if(.not. EN_denominator) return
    r = r - h%two%GetTwBMEMon(p1,p2,p1,p2) - h%two%GetTwBMEMon(h1,h2,h1,h2) + &
      & h%two%GetTwBMEMon(p1,h1,p1,h1) + h%two%GetTwBMEMon(p2,h2,p2,h2) + &
      & h%two%GetTwBMEMon(p2,h1,p2,h1) + h%two%GetTwBMEMon(p1,h2,p1,h2)
  end function get_denominator2

  function get_denominator3(h,h1,h2,h3,p1,p2,p3) result(r)
    ! E_{0} - <p1p2p3h1h2h3| H | p1p2p3h1h2h3 >
    type(Ops), intent(in) :: h
    integer, intent(in) :: h1, h2, h3, p1, p2, p3
    real(8) :: r

    r = h%one%GetDenominator3(h1,h2,h3,p1,p2,p3)
    if(.not. EN_denominator) return
    write(*,*) "EN denominator for 3p3h: Not implemented yet"
  end function get_denominator3

  function get_denominator4(h,h1,h2,h3,h4,p1,p2,p3,p4) result(r)
    ! E_{0} - <p1p2p3p4h1h2h3h4| H | p1p2p3p4h1h2h3h4 >
    type(Ops), intent(in) :: h
    integer, intent(in) :: h1, h2, h3, h4, p1, p2, p3, p4
    real(8) :: r

    r = h%one%GetDenominator4(h1,h2,h3,h4,p1,p2,p3,p4)
    if(.not. EN_denominator) return
    write(*,*) "EN denominator for 4p4h: Not implemented yet"
  end function get_denominator4

  function get_pphh_part(opch, h) result(m)
    type(TwoBodyPartChannel), intent(in) :: opch
    type(Ops), intent(in) :: h
    type(TwoBodyChannel), pointer :: ch_bra, ch_ket
    type(Orbits), pointer :: sps
    type(DMat) :: m
    integer :: ibra, iket, bra, ket
    integer :: a, b, c, d
    type(SingleParticleOrbit), pointer :: oa, ob, oc, od

    ch_bra => opch%ch_bra
    ch_ket => opch%ch_ket
    sps => h%ms%sps

    ibra = 0
    do bra = 1, ch_bra%n_state
      a = ch_bra%n2spi1(bra)
      b = ch_bra%n2spi2(bra)
      oa => sps%GetOrbit(a)
      ob => sps%GetOrbit(b)
      if(abs(oa%GetOccupation())+abs(ob%GetOccupation()) < 1.d-6) ibra = ibra+1
    end do

    iket = 0
    do bra = 1, ch_ket%n_state
      a = ch_ket%n2spi1(bra)
      b = ch_ket%n2spi2(bra)
      oa => sps%GetOrbit(a)
      ob => sps%GetOrbit(b)
      if(abs(oa%GetOccupation())*abs(ob%GetOccupation()) > 1.d-6) iket = iket+1
    end do

    call M%ini(ibra,iket)
    ibra = 0
    do bra = 1, ch_bra%n_state
      a = ch_bra%n2spi1(bra)
      b = ch_bra%n2spi2(bra)
      oa => sps%GetOrbit(a)
      ob => sps%GetOrbit(b)
      if(abs(oa%GetOccupation())+abs(ob%GetOccupation()) > 1.d-6) cycle
      ibra = ibra + 1
      iket = 0
      do ket = 1, ch_ket%n_state
        c = ch_ket%n2spi1(ket)
        d = ch_ket%n2spi2(ket)
        oc => sps%GetOrbit(c)
        od => sps%GetOrbit(d)
        if(abs(oc%GetOccupation())*abs(od%GetOccupation()) < 1.d-6) cycle
        iket = iket + 1
        M%m(ibra,iket) = opch%m(bra,ket) / get_denominator2(h,c,d,a,b)
      end do
    end do
  end function get_pphh_part

  function get_xc_pphh_to_phph(op, ch_cc, h) result(Mat)
    !  only for scalar
    !  _________________
    !  <ph:J| V |p'h':J> = \sum_{J'} [J'] {jp  jp' J'} <pp':J'|V|h'h:J'> / denominator
    !                                     {jh' jh  J }
    type(TwoBodyPart), intent(in) :: op
    type(Ops), intent(in) :: h
    type(CrossCoupledTwoBodyChannel), intent(in) :: ch_cc
    type(Orbits), pointer :: sps
    type(DMat) :: Mat
    integer :: a, b, c, d, K
    integer :: ibra, iket, bra, ket
    type(SingleParticleOrbit), pointer :: oa, ob, oc, od
    real(8) :: v, norm

    sps => op%two%sps
    K = ch_cc%j
    ibra = 0
    do bra = 1, ch_cc%n_state
      a = ch_cc%n2spi1(bra)
      b = ch_cc%n2spi2(bra)
      oa => sps%GetOrbit(a)
      ob => sps%GetOrbit(b)
      if(abs(oa%GetOccupation())+abs(ob%GetOccupation()) > 1.d-6 .and. &
          & abs(oa%GetOccupation())*abs(ob%GetOccupation()) < 1.d-6) ibra = ibra+1
    end do
    call Mat%zeros(ibra,ibra)
    if(ibra < 1) return
    ibra = 0
    do bra = 1, ch_cc%n_state
      a = ch_cc%n2spi1(bra) ! p
      b = ch_cc%n2spi2(bra) ! h
      oa => sps%GetOrbit(a)
      ob => sps%GetOrbit(b)
      if(abs(oa%GetOccupation())+abs(ob%GetOccupation()) < 1.d-6 .or. &
          & abs(oa%GetOccupation())*abs(oa%GetOccupation()) > 1.d-6) cycle
      ibra = ibra+1
      iket = 0
      do ket = 1, ch_cc%n_state
        c = ch_cc%n2spi1(ket) ! p
        d = ch_cc%n2spi2(ket) ! h
        oc => sps%GetOrbit(c)
        od => sps%GetOrbit(d)
        if(abs(oc%GetOccupation())+abs(od%GetOccupation()) < 1.d-6 .or. &
            & abs(oc%GetOccupation())*abs(od%GetOccupation()) > 1.d-6) cycle
        iket = iket+1
        if(oa%z+oc%z /= ob%z+od%z) cycle
        if(abs(oa%GetOccupation())+abs(oc%GetOccupation()) > 1.d-6) cycle
        if(abs(ob%GetOccupation())*abs(od%GetOccupation()) < 1.d-6) cycle
        norm = 1.d0
        if(a==c) norm = norm*sqrt(2.d0)
        if(b==d) norm = norm*sqrt(2.d0)
        v = op%get_xc1423(a,c,d,b,K)
        Mat%m(ibra,iket) = v * norm / get_denominator2(h,b,d,a,c)
      end do
    end do
  end function get_xc_pphh_to_phph

  function CheckParityZ(o1, o2, o3, o4) result(r)
    type(SingleParticleOrbit), intent(in) :: o1, o2, o3, o4
    logical :: r
    r = .false.
    if(mod(o1%l+o2%l+o3%l+o4%l,2)==1 .or. o1%z+o2%z-o3%z-o4%z/=0) r = .true.
  end function CheckParityZ

  function energy_third_4_M(h) result(r)
    type(Ops), intent(in) :: h
    type(MSPace), pointer :: ms
    type(TwoBodyChannel), pointer :: ch_two
    type(Orbits), pointer :: sps
    integer :: p1, p2, p3, h1, h2
    type(SingleParticleOrbit), pointer :: op1, op2, op3, oh1, oh2
    integer :: mp1, mp2, mp3, mh1, mh2
    real(8) :: vsum, v, ti, r

    ti = omp_get_wtime()
    ms => h%ms
    sps => ms%sps
    vsum = 0.d0
    do p1 = 1, sps%norbs
      op1 => sps%GetOrbit(p1)
      if(op1%GetOccupation() > 1.d-6) cycle
      do p2 = 1, sps%norbs
        op2 => sps%GetOrbit(p2)
        if(op2%GetOccupation() > 1.d-6) cycle
        do p3 = 1, sps%norbs
          op3 => sps%GetOrbit(p3)
          if(op3%GetOccupation() > 1.d-6) cycle
          do h1 = 1, sps%norbs
            oh1 => sps%GetOrbit(h1)
            if(oh1%GetOccupation() < 1.d-6) cycle
            do h2 = 1, sps%norbs
              oh2 => sps%GetOrbit(h2)
              if(oh2%GetOccupation() < 1.d-6) cycle

              if(CheckParityZ(op1,op2,oh1,oh2)) cycle
              if(CheckParityZ(op3,oh2,op1,op2)) cycle
              if(op3%l /= oh1%l) cycle
              if(op3%j /= oh1%j) cycle
              if(op3%z /= oh1%z) cycle

              do mp1 = -op1%j, op1%j, 2
                do mp2 = -op2%j, op2%j, 2
                  do mp3 = -op3%j, op3%j, 2
                    do mh1 = -oh1%j, oh1%j, 2
                      do mh2 = -oh2%j, oh2%j, 2

                        if(mp1+mp2 /= mh1+mh2) cycle
                        if(mp1+mp2 /= mp3+mh2) cycle
                        if(mh1 /= mp3) cycle

                        vsum = vsum + h%two%GetTwBME_mscheme(p1,mp1,p2,mp2,h1,mh1,h2,mh2,0) * &
                            & h%two%GetTwBME_mscheme(p3,mp3,h2,mh2,p1,mp1,p2,mp2,0) * &
                            & h%one%GetOBME(h1,p3) / &
                            & (get_denominator2(h,h1,h2,p1,p2) * get_denominator1(h,h1,p3) )

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
    r = 0.5d0 * vsum
  end function energy_third_4_M

  function energy_third_5_M(h) result(r)
    type(Ops), intent(in) :: h
    type(MSPace), pointer :: ms
    type(TwoBodyChannel), pointer :: ch_two
    type(Orbits), pointer :: sps
    integer :: p1, p2, h1, h2, h3
    type(SingleParticleOrbit), pointer :: op1, op2, oh1, oh2, oh3
    integer :: mp1, mp2, mh1, mh2, mh3
    real(8) :: vsum, v, ti, r

    ti = omp_get_wtime()
    ms => h%ms
    sps => ms%sps
    vsum = 0.d0
    do p1 = 1, sps%norbs
      op1 => sps%GetOrbit(p1)
      if(op1%GetOccupation() > 1.d-6) cycle

      do p2 = 1, sps%norbs
        op2 => sps%GetOrbit(p2)
        if(op2%GetOccupation() > 1.d-6) cycle

        do h1 = 1, sps%norbs
          oh1 => sps%GetOrbit(h1)
          if(oh1%GetOccupation() < 1.d-6) cycle

          do h2 = 1, sps%norbs
            oh2 => sps%GetOrbit(h2)
            if(oh2%GetOccupation() < 1.d-6) cycle

            do h3 = 1, sps%norbs
              oh3 => sps%GetOrbit(h3)
              if(oh3%GetOccupation() < 1.d-6) cycle

              if(CheckParityZ(op1,op2,oh1,oh2)) cycle
              if(CheckParityZ(oh1,oh2,oh3,op2)) cycle
              if(op1%l /= oh3%l) cycle
              if(op1%j /= oh3%j) cycle
              if(op1%z /= oh3%z) cycle

              do mp1 = -op1%j, op1%j, 2
                do mp2 = -op2%j, op2%j, 2
                  do mh1 = -oh1%j, oh1%j, 2
                    do mh2 = -oh2%j, oh2%j, 2
                      do mh3 = -oh3%j, oh3%j, 2

                        if(mp1+mp2 /= mh1+mh2) cycle
                        if(mh1+mh2 /= mh3+mp2) cycle
                        if(mp1 /= mh3) cycle

                        vsum = vsum + h%two%GetTwBME_mscheme(p1,mp1,p2,mp2,h1,mh1,h2,mh2,0) * &
                            & h%two%GetTwBME_mscheme(h1,mh1,h2,mh2,h3,mh3,p2,mp2,0) * &
                            & h%one%GetOBME(p1,h3) / &
                            & (get_denominator2(h,h1,h2,p1,p2) * get_denominator1(h,h3,p1) )

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
    r = (-0.5d0) * vsum
  end function energy_third_5_M

  function energy_third_6_M(h) result(r)
    type(Ops), intent(in) :: h
    type(MSPace), pointer :: ms
    type(TwoBodyChannel), pointer :: ch_two
    type(Orbits), pointer :: sps
    integer :: p1, p2, p3, h1, h2
    type(SingleParticleOrbit), pointer :: op1, op2, op3, oh1, oh2
    integer :: mp1, mp2, mp3, mh1, mh2
    real(8) :: vsum, v, ti, r

    ti = omp_get_wtime()
    ms => h%ms
    sps => ms%sps
    vsum = 0.d0
    do p1 = 1, sps%norbs
      op1 => sps%GetOrbit(p1)
      if(op1%GetOccupation() > 1.d-6) cycle

      do p2 = 1, sps%norbs
        op2 => sps%GetOrbit(p2)
        if(op2%GetOccupation() > 1.d-6) cycle

        do p3 = 1, sps%norbs
          op3 => sps%GetOrbit(p3)
          if(op3%GetOccupation() > 1.d-6) cycle

          do h1 = 1, sps%norbs
            oh1 => sps%GetOrbit(h1)
            if(oh1%GetOccupation() < 1.d-6) cycle

            do h2 = 1, sps%norbs
              oh2 => sps%GetOrbit(h2)
              if(oh2%GetOccupation() < 1.d-6) cycle

              if(CheckParityZ(op1,op2,oh1,oh2)) cycle
              if(CheckParityZ(oh1,oh2,op1,op3)) cycle
              if(op2%l /= op3%l) cycle
              if(op2%j /= op3%j) cycle
              if(op2%z /= op3%z) cycle

              do mp1 = -op1%j, op1%j, 2
                do mp2 = -op2%j, op2%j, 2
                  do mp3 = -op3%j, op3%j, 2
                    do mh1 = -oh1%j, oh1%j, 2
                      do mh2 = -oh2%j, oh2%j, 2

                        if(mp1+mp2 /= mh1+mh2) cycle
                        if(mp1+mp3 /= mh1+mh2) cycle
                        if(mp2 /= mp3) cycle

                        vsum = vsum + h%two%GetTwBME_mscheme(p1,mp1,p2,mp2,h1,mh1,h2,mh2,0) * &
                            & h%two%GetTwBME_mscheme(h1,mh1,h2,mh2,p1,mp1,p3,mp3,0) * &
                            & h%one%GetOBME(p2,p3) / &
                            & (get_denominator2(h,h1,h2,p1,p2) * get_denominator2(h,h1,h2,p1,p3) )

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
    r = 0.5d0 * vsum
  end function energy_third_6_M

  function energy_third_7_M(h) result(r)
    type(Ops), intent(in) :: h
    type(MSPace), pointer :: ms
    type(TwoBodyChannel), pointer :: ch_two
    type(Orbits), pointer :: sps
    integer :: p1, p2, h1, h2, h3
    type(SingleParticleOrbit), pointer :: op1, op2, oh1, oh2, oh3
    integer :: mp1, mp2, mh1, mh2, mh3
    real(8) :: vsum, v, ti, r

    ti = omp_get_wtime()
    ms => h%ms
    sps => ms%sps
    vsum = 0.d0
    do p1 = 1, sps%norbs
      op1 => sps%GetOrbit(p1)
      if(op1%GetOccupation() > 1.d-6) cycle

      do p2 = 1, sps%norbs
        op2 => sps%GetOrbit(p2)
        if(op2%GetOccupation() > 1.d-6) cycle

        do h1 = 1, sps%norbs
          oh1 => sps%GetOrbit(h1)
          if(oh1%GetOccupation() < 1.d-6) cycle

          do h2 = 1, sps%norbs
            oh2 => sps%GetOrbit(h2)
            if(oh2%GetOccupation() < 1.d-6) cycle

            do h3 = 1, sps%norbs
              oh3 => sps%GetOrbit(h3)
              if(oh3%GetOccupation() < 1.d-6) cycle

              if(CheckParityZ(op1,op2,oh1,oh2)) cycle
              if(CheckParityZ(op1,op2,oh3,oh2)) cycle
              if(oh1%l /= oh3%l) cycle
              if(oh1%j /= oh3%j) cycle
              if(oh1%z /= oh3%z) cycle

              do mp1 = -op1%j, op1%j, 2
                do mp2 = -op2%j, op2%j, 2
                  do mh1 = -oh1%j, oh1%j, 2
                    do mh2 = -oh2%j, oh2%j, 2
                      do mh3 = -oh3%j, oh3%j, 2

                        if(mp1+mp2 /= mh1+mh2) cycle
                        if(mp1+mp2 /= mh3+mh2) cycle
                        if(mh1 /= mh3) cycle

                        vsum = vsum + h%two%GetTwBME_mscheme(p1,mp1,p2,mp2,h1,mh1,h2,mh2,0) * &
                            & h%two%GetTwBME_mscheme(p1,mp1,p2,mp2,h3,mh3,h2,mh2,0) * &
                            & h%one%GetOBME(h1,h3) / &
                            & (get_denominator2(h,h1,h2,p1,p2) * get_denominator2(h,h2,h3,p1,p2) )

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
    r = (-0.5d0) * vsum
  end function energy_third_7_M

  function energy_third_10_M(h) result(r)
    type(Ops), intent(in) :: h
    type(MSPace), pointer :: ms
    type(TwoBodyChannel), pointer :: ch_two
    type(Orbits), pointer :: sps
    integer :: p1, p2, h1, h2
    type(SingleParticleOrbit), pointer :: op1, op2, oh1, oh2
    integer :: mp1, mp2, mh1, mh2
    real(8) :: vsum, v, ti, r

    ti = omp_get_wtime()
    ms => h%ms
    sps => ms%sps
    vsum = 0.d0
    do p1 = 1, sps%norbs
      op1 => sps%GetOrbit(p1)
      if(op1%GetOccupation() > 1.d-6) cycle

      do p2 = 1, sps%norbs
        op2 => sps%GetOrbit(p2)
        if(op2%GetOccupation() > 1.d-6) cycle

        do h1 = 1, sps%norbs
          oh1 => sps%GetOrbit(h1)
          if(oh1%GetOccupation() < 1.d-6) cycle

          do h2 = 1, sps%norbs
            oh2 => sps%GetOrbit(h2)
            if(oh2%GetOccupation() < 1.d-6) cycle

            if(CheckParityZ(op1,op2,oh1,oh2)) cycle
            if(oh1%l /= op1%l) cycle
            if(oh1%j /= op1%j) cycle
            if(oh1%z /= op1%z) cycle
            if(oh2%l /= op2%l) cycle
            if(oh2%j /= op2%j) cycle
            if(oh2%z /= op2%z) cycle

            do mp1 = -op1%j, op1%j, 2
              do mp2 = -op2%j, op2%j, 2
                do mh1 = -oh1%j, oh1%j, 2
                  do mh2 = -oh2%j, oh2%j, 2

                    if(mp1+mp2 /= mh1+mh2) cycle
                    if(mh1 /= mp1) cycle
                    if(mh2 /= mp2) cycle

                    vsum = vsum + h%two%GetTwBME_mscheme(p1,mp1,p2,mp2,h1,mh1,h2,mh2,0) * &
                        & h%one%GetOBME(h1,p1) * h%one%GetOBME(h2,p2) / &
                        & (get_denominator2(h,h1,h2,p1,p2) * get_denominator1(h,h1,p1) )

                  end do
                end do
              end do
            end do

          end do
        end do
      end do
    end do
    r = vsum
  end function energy_third_10_M

  function energy_third_11_M(h) result(r)
    type(Ops), intent(in) :: h
    type(MSPace), pointer :: ms
    type(TwoBodyChannel), pointer :: ch_two
    type(Orbits), pointer :: sps
    integer :: p1, p2, h1, h2
    type(SingleParticleOrbit), pointer :: op1, op2, oh1, oh2
    integer :: mp1, mp2, mh1, mh2
    real(8) :: vsum, v, ti, r

    ti = omp_get_wtime()
    ms => h%ms
    sps => ms%sps
    vsum = 0.d0
    do p1 = 1, sps%norbs
      op1 => sps%GetOrbit(p1)
      if(op1%GetOccupation() > 1.d-6) cycle

      do p2 = 1, sps%norbs
        op2 => sps%GetOrbit(p2)
        if(op2%GetOccupation() > 1.d-6) cycle

        do h1 = 1, sps%norbs
          oh1 => sps%GetOrbit(h1)
          if(oh1%GetOccupation() < 1.d-6) cycle

          do h2 = 1, sps%norbs
            oh2 => sps%GetOrbit(h2)
            if(oh2%GetOccupation() < 1.d-6) cycle

            if(CheckParityZ(op2,oh2,oh1,op1)) cycle
            if(oh2%l /= op1%l) cycle
            if(oh2%j /= op1%j) cycle
            if(oh2%z /= op1%z) cycle
            if(oh1%l /= op2%l) cycle
            if(oh1%j /= op2%j) cycle
            if(oh1%z /= op2%z) cycle

            do mp1 = -op1%j, op1%j, 2
              do mp2 = -op2%j, op2%j, 2
                do mh1 = -oh1%j, oh1%j, 2
                  do mh2 = -oh2%j, oh2%j, 2

                    if(mp2+mh2 /= mh1+mp1) cycle
                    if(mh2 /= mp1) cycle
                    if(mh1 /= mp2) cycle

                    vsum = vsum + h%two%GetTwBME_mscheme(p2,mp2,h2,mh2,h1,mh1,p1,mp1,0) * &
                        & h%one%GetOBME(h2,p1) * h%one%GetOBME(p2,h1) / &
                        & (get_denominator1(h,h2,p1) * get_denominator1(h,h1,p2) )

                  end do
                end do
              end do
            end do

          end do
        end do
      end do
    end do
    r = vsum
  end function energy_third_11_M

  function energy_third_13_M(h) result(r)
    type(Ops), intent(in) :: h
    type(MSPace), pointer :: ms
    type(TwoBodyChannel), pointer :: ch_two
    type(Orbits), pointer :: sps
    integer :: p1, p2, h1
    type(SingleParticleOrbit), pointer :: op1, op2, oh1
    integer :: mp1, mp2, mh1
    real(8) :: vsum, v, ti, r

    ti = omp_get_wtime()
    ms => h%ms
    sps => ms%sps
    vsum = 0.d0
    do p1 = 1, sps%norbs
      op1 => sps%GetOrbit(p1)
      if(op1%GetOccupation() > 1.d-6) cycle

      do p2 = 1, sps%norbs
        op2 => sps%GetOrbit(p2)
        if(op2%GetOccupation() > 1.d-6) cycle

        do h1 = 1, sps%norbs
          oh1 => sps%GetOrbit(h1)
          if(oh1%GetOccupation() < 1.d-6) cycle

          if(oh1%l /= op1%l) cycle
          if(oh1%j /= op1%j) cycle
          if(oh1%z /= op1%z) cycle
          if(oh1%l /= op2%l) cycle
          if(oh1%j /= op2%j) cycle
          if(oh1%z /= op2%z) cycle

          do mp1 = -op1%j, op1%j, 2
            do mp2 = -op2%j, op2%j, 2
              do mh1 = -oh1%j, oh1%j, 2

                if(mh1 /= mp1) cycle
                if(mh1 /= mp2) cycle

                vsum = vsum + h%one%GetOBME(h1,p1) * &
                    & h%one%GetOBME(p1,p2) * h%one%GetOBME(p2,h1) / &
                    & (get_denominator1(h,h1,p1) * get_denominator1(h,h1,p2) )

              end do
            end do
          end do

        end do
      end do
    end do
    r = vsum
  end function energy_third_13_M

  function energy_third_14_M(h) result(r)
    type(Ops), intent(in) :: h
    type(MSPace), pointer :: ms
    type(TwoBodyChannel), pointer :: ch_two
    type(Orbits), pointer :: sps
    integer :: p1, h1, h2
    type(SingleParticleOrbit), pointer :: op1, oh1, oh2
    integer :: mp1, mh1, mh2
    real(8) :: vsum, v, ti, r

    ti = omp_get_wtime()
    ms => h%ms
    sps => ms%sps
    vsum = 0.d0
    do p1 = 1, sps%norbs
      op1 => sps%GetOrbit(p1)
      if(op1%GetOccupation() > 1.d-6) cycle

      do h1 = 1, sps%norbs
        oh1 => sps%GetOrbit(h1)
        if(oh1%GetOccupation() < 1.d-6) cycle

        do h2 = 1, sps%norbs
          oh2 => sps%GetOrbit(h2)
          if(oh2%GetOccupation() < 1.d-6) cycle

          if(oh2%l /= op1%l) cycle
          if(oh2%j /= op1%j) cycle
          if(oh2%z /= op1%z) cycle
          if(oh1%l /= op1%l) cycle
          if(oh1%j /= op1%j) cycle
          if(oh1%z /= op1%z) cycle

            do mp1 = -op1%j, op1%j, 2
              do mh1 = -oh1%j, oh1%j, 2
                do mh2 = -oh2%j, oh2%j, 2

                if(mh1 /= mp1) cycle
                if(mh2 /= mp1) cycle

                vsum = vsum + h%one%GetOBME(p1,h1) * &
                    & h%one%GetOBME(h1,h2) * h%one%GetOBME(h2,p1) / &
                    & (get_denominator1(h,h1,p1) * get_denominator1(h,h2,p1) )

              end do
            end do
          end do

        end do
      end do
    end do
    r = vsum * (-1.d0)
  end function energy_third_14_M

end module gMBPT

