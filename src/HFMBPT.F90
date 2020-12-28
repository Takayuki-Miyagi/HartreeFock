module HFMBPT
  use omp_lib
  use myfort
  use ModelSpace
  use Operators
  use HartreeFock
  implicit none

  public :: MBPTEnergy, MBPTScalar, MBPTDMat

  ! Methods for MBPTEnergy
  private :: CalcEnergyCorr
  private :: MBPTCriteria
  private :: energy_second
  private :: energy_third_hh
  private :: energy_third_pp
  private :: energy_third_ph
  private :: energy_fourth_F1
  private :: energy_fourth_F3
  private :: energy_fourth_F4
  private :: energy_fourth_F5
  !private :: energy_fourth_F7
  !private :: energy_fourth_F8
  !private :: energy_fourth_F10
  !private :: energy_fourth_F12
  !private :: energy_fourth_F13
  !private :: energy_fourth_F14
  !private :: energy_fourth_F15
  !private :: energy_fourth_F16
  !private :: energy_fourth_F17
  !private :: energy_fourth_F19
  !private :: energy_fourth_F20
  !private :: energy_fourth_F21
  !private :: energy_fourth_F23
  !private :: energy_fourth_F25
  !private :: energy_fourth_F26
  !private :: energy_fourth_F27
  !private :: energy_fourth_F28
  !private :: energy_fourth_F29
  !private :: energy_fourth_F31
  !private :: energy_fourth_F32
  !private :: energy_fourth_F33
  !private :: energy_fourth_F34
  !private :: energy_fourth_F35
  !private :: energy_fourth_F36
  !private :: energy_fourth_F37
  !private :: energy_fourth_F38
  !private :: energy_fourth_F39

  ! Methods for Scalar Operator
  private :: CalcScalarCorr
  private :: scalar_first
  private :: scalar_first1
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

  ! Methods for density matrix
  private :: InitMBPTDMat
  private :: FinMBPTDMat
  private :: density_matrix_element_pp
  private :: density_matrix_element_hh
  private :: density_matrix_element_ph
  private :: GetCoef

  ! Methods for general
  private :: get_denominator1
  private :: get_denominator2
  private :: get_denominator3
  private :: get_denominator4
  private :: get_pphh_part
  private :: get_xc_pphh_to_phph

  type :: MBPTEnergy
    real(8) :: e_0 = 0.d0
    real(8) :: e_2 = 0.d0
    real(8) :: e_3 = 0.d0
    real(8) :: e_3_pp = 0.d0
    real(8) :: e_3_hh = 0.d0
    real(8) :: e_3_ph = 0.d0
    real(8) :: e_4(39)= 0.d0
    real(8) :: perturbativity1b = 0.d0
    real(8) :: perturbativity2b = 0.d0
    real(8) :: energy_gap = 0.d0
  contains
    procedure :: calc => CalcEnergyCorr
    procedure :: energy_third
    procedure :: energy_fourth
    procedure :: MBPTCriteria
  end type MBPTEnergy

  type :: MBPTDMat
    type(OneBodyPart) :: rho_HF
    type(OneBodyPart) :: rho_HO
    type(OneBodyPart) :: C_HO2HF
    type(OneBodyPart) :: C_HF2NAT
    type(OneBodyPart) :: C_HO2NAT
    type(OneBodyPart) :: Occ_HF
    type(OneBodyPart) :: Occ_HO
    type(MSpace), pointer :: ms
  contains
    procedure :: init => InitMBPTDMat
    procedure :: fin => FinMBPTDMat
    procedure :: GetCoef
  end type MBPTDMat

  type :: MBPTScalar
    logical :: part_of_hamil = .false.
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
    procedure :: scalar_second
  end type MBPTScalar
  type(SixJsStore), private :: sixjs
  logical, private :: EN_denominator=.false.  ! Epstein-Nesbet denominator (Moller-Plesset is default)

contains

  subroutine CalcEnergyCorr(this,hamil,fourth_order, EN_denominator_in)
    class(MBPTEnergy), intent(inout) :: this
    type(Ops), intent(in) :: hamil
    logical, intent(in), optional :: fourth_order
    logical, intent(in), optional :: EN_denominator_in
    real(8) :: s_1
    logical :: is_4th_order=.false.
    integer :: i, jmax

    if(present(fourth_order)) is_4th_order=fourth_order
    if(present(EN_denominator_in)) EN_denominator = EN_denominator_in
    write(*,*)

    if(.not. is_4th_order) then
      write(*,'(a)') "# Many-body perturbation calculation up to 3rd order"
    end if

    if(is_4th_order) then
      write(*,'(a)') "# Many-body perturbation calculation up to 4th order"
    end if
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
    call this%MBPTCriteria(hamil)

    this%e_0 = hamil%zero
    this%e_2 = energy_second(hamil)
    write(*,'(a,f16.8)') "Second order correction: ", this%e_2

    call this%energy_third(hamil)
    write(*,'(a,4f16.8)') "Third order corrections pp, hh, ph, and total: ", &
        & this%e_3_pp, this%e_3_hh, this%e_3_ph, this%e_3
    if(.not. is_4th_order) then
      write(*,'(a,f16.8)') "Third order MBPT energy: ", this%e_0 + this%e_2 + this%e_3
    end if

    if(is_4th_order) then
      call this%energy_fourth(hamil)
      do i = 1, 39
        write(*,"(a,i2,f16.8)") "Fourth order correction F", i, this%e_4(i)
      end do
      write(*,'(a,f16.8)') "Total Fourth order correction: ", sum(this%e_4)
      write(*,'(a,f16.8)') "Fourth order MBPT energy: ", this%e_0 + &
          & this%e_2 + this%e_3 + sum(this%e_4)
    end if
    call sixjs%fin()
    EN_denominator=.false.
  end subroutine CalcEnergyCorr

  function energy_second(h) result(r)
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
    integer :: ch, ab, ij, J2, n
    integer :: a, b, i, j
    real(8) :: vsum, v, ti, r

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
  end function energy_second

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

  subroutine energy_third(this,h)
    class(MBPTEnergy), intent(inout) :: this
    type(Ops), intent(in) :: h

    this%e_3_hh = energy_third_hh(h)
    this%e_3_pp = energy_third_pp(h)
    this%e_3_ph = energy_third_ph(h)
    this%e_3 = this%e_3_hh + this%e_3_pp + this%e_3_ph
  end subroutine energy_third

  subroutine energy_fourth(this,h)
    class(MBPTEnergy), intent(inout) :: this
    type(Ops), intent(in) :: h
    real(8) :: e

    this%e_4(1) = energy_fourth_F1(h)
    this%e_4(2) = this%e_4(1)
    this%e_4(3) = energy_fourth_F3(h)
    this%e_4(4) = energy_fourth_F4(h)
    this%e_4(5) = energy_fourth_F5(h)
    this%e_4(6) = this%e_4(5)
    !call this%energy_fourth_F7(h)
    !call this%energy_fourth_F8(h)
    !this%e_4(9) = this%e_4(8)
    !call this%energy_fourth_F10(h)
    !this%e_4(11) = this%e_4(10)
    !call this%energy_fourth_F12(h)
    !call this%energy_fourth_F13(h)
    !call this%energy_fourth_F14(h)
    !call this%energy_fourth_F15(h)
    !call this%energy_fourth_F16(h)
    !call this%energy_fourth_F17(h)
    !this%e_4(18) = this%e_4(17)
    !call this%energy_fourth_F19(h)
    !call this%energy_fourth_F20(h)
    !call this%energy_fourth_F21(h)
    !this%e_4(22) = this%e_4(21)
    !call this%energy_fourth_F23(h)
    !this%e_4(24) = this%e_4(23)
    !call this%energy_fourth_F25(h)
    !call this%energy_fourth_F26(h)
    !call this%energy_fourth_F27(h)
    !call this%energy_fourth_F28(h)
    !call this%energy_fourth_F29(h)
    !this%e_4(30) = this%e_4(31)
    !call this%energy_fourth_F31(h)
    !call this%energy_fourth_F32(h)
    !call this%energy_fourth_F33(h)
    !call this%energy_fourth_F34(h)
    !call this%energy_fourth_F35(h)
    !call this%energy_fourth_F36(h)
    !call this%energy_fourth_F37(h)
    !call this%energy_fourth_F38(h)
    !call this%energy_fourth_F39(h)
  end subroutine energy_fourth

  function energy_third_pp(h) result(r)
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
  end function energy_third_pp

  function energy_third_hh(h) result(r)
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
  end function energy_third_hh

  function energy_third_ph(h) result(r)
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
  end function energy_third_ph

  function energy_fourth_F1(h) result(r)
    ! p1, p2, p3, p4: particle
    ! h1, h2, h3, h4: hole
    !     _____________
    !    /\           /\
    !   /  \ h1   p2 /  \ h2
    !p1 |  |         \  /
    !   |__|__________\/
    !   |  |
    !   |  |
    !   |  |__________
    !p3 |  |          /\
    !   |  | h3   p4 /  \ h4
    !   \  /         \  /
    !    \/___________\/
    !
    ! F1 = (1/4) [J1][J2][K][L] (-1)^(J1+J2+K+L)
    ! { jh1 jh2 J1 } { jh1 jp4 J2 } { jh2 jh2 K }
    ! { jh2 jp3 K  } { jp4 jp3 K  } { jp4 jp4 L }
    ! <p1p2:J1||h1h2:J1> <p3h2:J1||p1p2:J1>
    ! <h1p4:J2||h3h4:J2> <h3h4:J2||p3p4:J2> / denominator
    !
    type(Ops), intent(in) :: h
    type(MSPace), pointer :: ms
    type(TwoBodyChannel), pointer :: tbs
    type(Orbits), pointer :: sps
    type(Ops) :: v_h1h2p3h2, v_h1p4p3p4
    type(DMat) :: m1, m2
    integer :: ch
    integer :: h1, h2, p3, p4, L, Kmin, Kmax, K
    type(SingleParticleOrbit), pointer :: oh1, oh2, op3, op4
    real(8) :: r, ti, vsum, v1, v2, denom, norm

    ti = omp_get_wtime()
    ms => h%ms
    sps => ms%sps

    ! intermediate stuffs
    call v_h1h2p3h2%init(0,1,0,"v for MBPT",ms,2)
    call v_h1p4p3p4%init(0,1,0,"v for MBPT",ms,2)
    do ch = 1, ms%two%NChan
      tbs => ms%two%jpz(ch)
      m1 = get_pphh_part(h%two%MatCh(ch,ch), h) * (-1.d0) ** tbs%j
      if(m1%n_row<1 .or. m1%n_col<1) cycle

      m2 = h%two%MatCh(ch,ch)%get_phpp(ms%sps)
      if(m2%n_row>0 .and. m2%n_col>0) call v_h1h2p3h2%two%MatCh(ch,ch)%set_phhh(ms%sps, m2 * m1)

      m2 = h%two%MatCh(ch,ch)%get_phhh(ms%sps)
      if(m2%n_row>0 .and. m2%n_col>0) call v_h1p4p3p4%two%MatCh(ch,ch)%set_phpp(ms%sps, m2 * m1%t())
      call m1%fin()
      call m2%fin()
    end do

    vsum = 0.d0
    do h2 = 1, sps%norbs
      oh2 => sps%GetOrbit(h2)
      if(oh2%GetOccupation() < 1.d-6) cycle
      do p4 = 1, sps%norbs
        op4 => sps%GetOrbit(p4)
        if(op4%GetOccupation() > 1.d-6) cycle

        do h1 = 1, sps%norbs
          oh1 => sps%GetOrbit(h1)
          if(oh1%GetOccupation() < 1.d-6) cycle
          do p3 = 1, sps%norbs
            op3 => sps%GetOrbit(p3)
            if(op3%GetOccupation() > 1.d-6) cycle
            norm = 1.d0
            if(h1==h2) norm = norm * sqrt(2.d0)
            if(p3==p4) norm = norm * sqrt(2.d0)
            denom = 1.d0 / get_denominator1(h,h1,p3)
            Kmin = abs(oh1%j-op3%j)/2
            Kmax = min(oh1%j+op3%j, 2*oh2%j, 2*op4%j)/2
            do K = Kmin, Kmax
              v1 = v_h1h2p3h2%two%get_xc1324(h1,h2,p3,h2,K)
              v2 = v_h1p4p3p4%two%get_xc1324(h1,p4,p3,p4,K)
              if(abs(v1*v2) < 1.d-12) cycle
              do L = abs(oh2%j-op4%j)/2, (oh2%j+op4%j)/2
                vsum = vsum + dble(2*K+1) * dble(2*L+1) * v1 * v2 * &
                    & sixjs%get(oh2%j, oh2%j, 2*K, op4%j, op4%j, 2*L) * &
                    & (-1.d0) ** (K+L) * denom * norm
              end do

            end do
          end do

        end do
      end do
    end do

    call v_h1h2p3h2%fin()
    call v_h1p4p3p4%fin()
    r = vsum
    call timer%Add("Fourth order MBPT F1",omp_get_wtime()-ti)
  end function energy_fourth_F1

  function energy_fourth_F3(h) result(r)
    ! p1, p2, p3, p4, p5: particle
    ! h1, h2, h3:         hole
    !     _____________
    !    /\           /\
    !   /  \ h1   p2 /  \ h2
    !p1 |  |         \  /
    !   |__|__________\/
    !   |  |
    !p3 |  |
    !   |__|__________
    !   |  |          /\
    !p4 |  |      p5 /  \ h3
    !   \  /         \  /
    !    \/___________\/
    !
    ! F3 = - (1/4) [J1][J2][K][L] (-1)^(J1+J2+K+L)
    ! { jp3 jh2 J1 } { jp3 jh3 J2 } { jh2 jh2 K }
    ! { jh2 jh1 K  } { jh3 jh1 K  } { jh3 jh3 L }
    ! <p1p2:J1||h1h2:J1> <p3h2:J1||p1p2:J1>
    ! <p4p5:J2||p3h3:J2> <h1h3:J2||p4p5:J2> / denominator
    !
    type(Ops), intent(in) :: h
    type(MSPace), pointer :: ms
    type(TwoBodyChannel), pointer :: tbs
    type(Orbits), pointer :: sps
    type(Ops) :: v_phhh
    type(DMat) :: m1, m2
    integer :: ch
    integer :: h1, h2, h3, p3, L, Kmin, Kmax, K
    type(SingleParticleOrbit), pointer :: oh1, oh2, oh3, op3
    real(8) :: r, ti, vsum, v1, v2, denom, norm

    ti = omp_get_wtime()
    ms => h%ms
    sps => ms%sps

    ! intermediate stuffs
    call v_phhh%init(0,1,0,"v for MBPT",ms,2)
    do ch = 1, ms%two%NChan
      tbs => ms%two%jpz(ch)
      m1 = get_pphh_part(h%two%MatCh(ch,ch), h) * (-1.d0) ** tbs%j
      if(m1%n_row<1 .or. m1%n_col<1) cycle
      m2 = h%two%MatCh(ch,ch)%get_phpp(ms%sps)
      if(m2%n_row<1 .or. m2%n_col<1) cycle
      call v_phhh%two%MatCh(ch,ch)%set_phhh(ms%sps, m2 * m1)
      call m1%fin()
      call m2%fin()
    end do

    vsum = 0.d0
    do h1 = 1, sps%norbs
      oh1 => sps%GetOrbit(h1)
      if(oh1%GetOccupation() < 1.d-6) cycle
      do h2 = 1, sps%norbs
        oh2 => sps%GetOrbit(h2)
        if(oh2%GetOccupation() < 1.d-6) cycle
        do h3 = 1, sps%norbs
          oh3 => sps%GetOrbit(h3)
          if(oh3%GetOccupation() < 1.d-6) cycle

          do p3 = 1, sps%norbs
            op3 => sps%GetOrbit(p3)
            if(op3%GetOccupation() > 1.d-6) cycle

            norm = 1.d0
            if(h1==h2) norm = norm * sqrt(2.d0)
            if(h1==h3) norm = norm * sqrt(2.d0)
            denom = 1.d0 / get_denominator1(h,h1,p3)
            Kmin = abs(oh1%j-op3%j)/2
            Kmax = min(oh1%j+op3%j, 2*oh2%j, 2*oh3%j)/2
            do K = Kmin, Kmax
              v1 = v_phhh%two%get_xc1324(p3,h2,h1,h2,K)
              v2 = v_phhh%two%get_xc1324(p3,h3,h1,h3,K)
              if(abs(v1*v2) < 1.d-12) cycle
              do L = abs(oh2%j-oh3%j)/2, (oh2%j+oh3%j)/2
                vsum = vsum + dble(2*K+1) * dble(2*L+1) * v1 * v2 * &
                    & sixjs%get(oh2%j, oh2%j, 2*K, oh3%j, oh3%j, 2*L) * &
                    & (-1.d0) ** (K+L) * denom * norm
              end do

            end do

          end do
        end do
      end do
    end do

    call v_phhh%fin()
    r = - vsum
    call timer%Add("Fourth order MBPT F3",omp_get_wtime()-ti)
  end function energy_fourth_F3

  function energy_fourth_F4(h) result(r)
    ! p1, p2, p3:         particle
    ! h1, h2, h3, h4, h5: hole
    !     _____________
    !    /\           /\
    !   /  \ h1   p2 /  \ h2
    !p1 |  |         \  /
    !   |  |__________\/
    !   |  |
    !   |  | h3
    !   |  |__________
    !   |  |          /\
    !   |  | h4   p3 /  \ h5
    !   \  /         \  /
    !    \/___________\/
    !
    ! F4 = - (1/4) [J1][J2][K][L] (-1)^(J1+J2+K+L)
    ! { jp1 jp2 J1 } { jp1 jp3 J2 } { jp2 jp2 K }
    ! { jp2 jh3 K  } { jp3 jh3 K  } { jp3 jp3 L }
    ! <p1p2:J1||h1h2:J1> <h1h2:J1||h3p2:J1>
    ! <h3p3:J2||h4h5:J2> <h4h5:J2||p1p3:J2> / denominator
    !
    type(Ops), intent(in) :: h
    type(MSPace), pointer :: ms
    type(TwoBodyChannel), pointer :: tbs
    type(Orbits), pointer :: sps
    type(Ops) :: v_phpp
    type(DMat) :: m1, m2
    integer :: ch
    integer :: p1, p2, p3, h3, L, Kmin, Kmax, K
    type(SingleParticleOrbit), pointer :: op1, op2, op3, oh3
    real(8) :: r, ti, vsum, v1, v2, denom, norm

    ti = omp_get_wtime()
    ms => h%ms
    sps => ms%sps

    ! intermediate stuffs
    call v_phpp%init(0,1,0,"v for MBPT",ms,2)
    do ch = 1, ms%two%NChan
      tbs => ms%two%jpz(ch)
      m1 = get_pphh_part(h%two%MatCh(ch,ch), h) * (-1.d0) ** tbs%j
      if(m1%n_row<1 .or. m1%n_col<1) cycle
      m2 = h%two%MatCh(ch,ch)%get_phhh(ms%sps)
      if(m2%n_row<1 .or. m2%n_col<1) cycle
      call v_phpp%two%MatCh(ch,ch)%set_phpp(ms%sps, m2*m1%t())
      call m1%fin()
      call m2%fin()
    end do

    vsum = 0.d0
    do h3 = 1, sps%norbs
      oh3 => sps%GetOrbit(h3)
      if(oh3%GetOccupation() < 1.d-6) cycle

      do p1 = 1, sps%norbs
        op1 => sps%GetOrbit(p1)
        if(op1%GetOccupation() > 1.d-6) cycle

        do p2 = 1, sps%norbs
          op2 => sps%GetOrbit(p2)
          if(op2%GetOccupation() > 1.d-6) cycle

          do p3 = 1, sps%norbs
            op3 => sps%GetOrbit(p3)
            if(op3%GetOccupation() > 1.d-6) cycle

            norm = 1.d0
            if(p1==p2) norm = norm * sqrt(2.d0)
            if(p1==p3) norm = norm * sqrt(2.d0)
            denom = 1.d0 / get_denominator1(h,h3,p1)
            Kmin = abs(op1%j-oh3%j)/2
            Kmax = min(op1%j+oh3%j, 2*op2%j, 2*op3%j)/2
            do K = Kmin, Kmax
              v1 = v_phpp%two%get_xc1324(p1,p2,h3,p2,K)
              v2 = v_phpp%two%get_xc1324(p1,p3,h3,p3,K)
              if(abs(v1*v2) < 1.d-12) cycle
              do L = abs(op2%j-op3%j)/2, (op2%j+op3%j)/2
                vsum = vsum + dble(2*K+1) * dble(2*L+1) * v1 * v2 * &
                    & sixjs%get(op2%j, op2%j, 2*K, op3%j, op3%j, 2*L) * &
                    & (-1.d0) ** (K+L) * denom * norm
              end do
            end do

          end do
        end do
      end do
    end do

    call v_phpp%fin()
    r = - vsum
    call timer%Add("Fourth order MBPT F4",omp_get_wtime()-ti)
  end function energy_fourth_F4

  function energy_fourth_F5(h) result(r)
    ! p1, p2, p3, p4: particle
    ! h1, h2, h3, h4: hole
    !     _____________
    !    /\           /\
    !p1 /  \ h1   p2 /  \ h2
    !   |  |         |  |
    !   |__|_________|  |
    !   |  |         |  |
    !p3 |  |      p4 |  |
    !   |  |_________|__|
    !   |  |         |  |
    !   |  | h3      |  | h4
    !   \  /         \  /
    !    \/___________\/
    !
    ! F5 = (1/16) [J]
    ! <h1h2:J||p1p2:J> <p1p2:J||p3p4:J>
    ! <p3p4:J||h3h4:J> <h3h4:J||h1h2:J> / denominator
    !
    type(Ops), intent(in) :: h
    type(MSPace), pointer :: ms
    type(TwoBodyChannel), pointer :: tbs
    type(TwoBodyChannel), pointer :: ch_two
    type(Orbits), pointer :: sps
    type(Ops) :: v_pphh_1, v_pphh_2
    type(DMat) :: m1, m2
    integer :: ch, ab, a, b, ij, i, j, J2
    real(8) :: r, ti, v, vsum, denom

    ti = omp_get_wtime()
    ms => h%ms
    sps => ms%sps

    ! intermediate stuffs
    call v_pphh_1%init(0,1,0,"v for MBPT",ms,2)
    call v_pphh_2%init(0,1,0,"v for MBPT",ms,2)
    do ch = 1, ms%two%NChan
      tbs => ms%two%jpz(ch)
      m1 = get_pphh_part(h%two%MatCh(ch,ch), h)
      if(m1%n_row<1 .or. m1%n_col<1) cycle
      m2 = h%two%MatCh(ch,ch)%get_pppp(ms%sps)
      if(m2%n_row>0 .and. m2%n_col>0) then
        call v_pphh_1%two%MatCh(ch,ch)%set_pphh(ms%sps, m2*m1)
      end if

      m2 = h%two%MatCh(ch,ch)%get_hhhh(ms%sps)
      if(m2%n_row>0 .and. m2%n_col>0) then
        call v_pphh_2%two%MatCh(ch,ch)%set_pphh(ms%sps, m1*m2)
      end if
      call m1%fin()
      call m2%fin()
    end do

    vsum = 0.d0
    do ch = 1, ms%two%NChan
      ch_two => ms%two%jpz(ch)
      J2 = ch_two%j
      if(ch_two%n_hh_state < 1) cycle
      if(ch_two%n_pp_state < 1) cycle
      v = 0.d0
      !$omp parallel
      !$omp do private(ab,a,b,ij,i,j,denom) reduction(+:v)
      do ab = 1, ch_two%n_pp_state
        a = ch_two%n2spi1(ch_two%pp_s(ab))
        b = ch_two%n2spi2(ch_two%pp_s(ab))
        do ij = 1, ch_two%n_hh_state
          i = ch_two%n2spi1(ch_two%hh_s(ij))
          j = ch_two%n2spi2(ch_two%hh_s(ij))
          denom = 1.d0 / get_denominator2(h,i,j,a,b)

          v = v + v_pphh_1%two%GetTwBME(i,j,a,b,J2) * &
              &   v_pphh_2%two%GetTwBME(a,b,i,j,J2) * denom
        end do
      end do
      !$omp end do
      !$omp end parallel
      vsum = vsum + v * dble(2*J2+1)
    end do
    call v_pphh_1%fin()
    call v_pphh_2%fin()
    r = vsum
    call timer%Add("Fourth order MBPT F5",omp_get_wtime()-ti)
  end function energy_fourth_F5

  function energy_fourth_F7(h) result(r)
    ! p1, p2, p3, p4: particle
    ! h1, h2, h3, h4: hole
    !     _____________
    !    /\           /\
    !p1 /  \ h1   p2 /  \ h2
    !   |  |         |  |
    !   |__|_________|__|
    !   |  |         |  |
    !p3 |  |         |  | h4
    !   |  |_________|  |
    !   |  |         |  |
    !   |  | h3   p4 |  |
    !   \  /         \  /
    !    \/___________\/
    !
    ! F7 = [K][J1][J2][J3][J4]
    ! { jh1 jh2 J1 } { jp3 jh2 J2 } { jh1 jp4 J3 } { jp3 jp4 J4 }
    ! { jp1 jp2 K  } { jp1 jh4 K  } { jh3 jp2 K  } { jh3 jh4 K  }
    ! <h1h2:J1||p1p2:J1> <p3h2:J2||p1h4:J2>
    ! <h1p4:J3||h3p2:J3> <p3p4:J4||h3h4:J4> / denominator
    !
    type(Ops), intent(in) :: h
    type(MSPace), pointer :: ms
    type(Orbits), pointer :: sps
    integer :: ch, J2
    type(CrossCoupledTwoBodyChannel), pointer :: ch_cc
    type(DMat) :: m1, m2, m3
    real(8) :: v, vsum, ti, r

    r = 0.d0
    ti = omp_get_wtime()
    ms => h%ms
    sps => ms%sps
    do ch = 1,ms%cc_two%NChan
      ch_cc => ms%cc_two%jpz(ch)
      J2 = ch_cc%j
      m1 = get_xc_pphh_to_phph(h%two, ch_cc, h)
      m2 = h%two%get_xc_hphp2phph(ch_cc)
      m3 = m2 * m1%t()
    end do

    !r = vsum
    call timer%Add("Fourth order MBPT F7",omp_get_wtime()-ti)
  end function energy_fourth_F7

!  function energy_fourth_F7(h) result(r)
!    ! p1, p2, p3, p4: particle
!    ! h1, h2, h3, h4: hole
!    !     _____________
!    !    /\           /\
!    !p1 /  \ h1   p2 /  \ h2
!    !   |  |         |  |
!    !   |__|_________|__|
!    !   |  |         |  |
!    !p3 |  |         |  | h4
!    !   |  |_________|  |
!    !   |  |         |  |
!    !   |  | h3   p4 |  |
!    !   \  /         \  /
!    !    \/___________\/
!    !
!    ! F7 = (1/16) [J]
!    ! <h1h2:J||p1p2:J> <p1p2:J||p3p4:J>
!    ! <p3p4:J||h3h4:J> <h3h4:J||h1h2:J> / denominator
!    !
!    type(Ops), intent(in) :: h
!    type(MSPace), pointer :: ms
!    type(TwoBodyChannel), pointer :: tbs
!    type(TwoBodyChannel), pointer :: ch_two
!    type(Orbits), pointer :: sps
!    type(Ops) :: v_pphh_1, v_pphh_2
!    type(DMat) :: m1, m2
!    integer :: ch, ab, a, b, ij, i, j, J2
!    real(8) :: r, ti, v, vsum, denom
!
!    ti = omp_get_wtime()
!    ms => h%ms
!    sps => ms%sps
!
!    ! intermediate stuffs
!    call v_pphh_1%init(0,1,0,"v for MBPT",ms,2)
!    call v_pphh_2%init(0,1,0,"v for MBPT",ms,2)
!    do ch = 1, ms%two%NChan
!      tbs => ms%two%jpz(ch)
!      m1 = h%two%MatCh(ch,ch)%get_pphh(ms%sps, h%one)
!      if(m1%n_row<1 .or. m1%n_col<1) cycle
!      m2 = h%two%MatCh(ch,ch)%get_pppp(ms%sps)
!      if(m2%n_row>0 .and. m2%n_col>0) then
!        call v_pphh_1%two%MatCh(ch,ch)%set_pphh(ms%sps, m2*m1)
!      end if
!
!      m2 = h%two%MatCh(ch,ch)%get_hhhh(ms%sps)
!      if(m2%n_row>0 .and. m2%n_col>0) then
!        call v_pphh_2%two%MatCh(ch,ch)%set_pphh(ms%sps, m1*m2)
!      end if
!      call m1%fin()
!      call m2%fin()
!    end do
!
!    vsum = 0.d0
!    do ch = 1, ms%two%NChan
!      ch_two => ms%two%jpz(ch)
!      J2 = ch_two%j
!      if(ch_two%n_hh_state < 1) cycle
!      if(ch_two%n_pp_state < 1) cycle
!      v = 0.d0
!      !$omp parallel
!      !$omp do private(ab,a,b,ij,i,j,denom) reduction(+:v)
!      do ab = 1, ch_two%n_pp_state
!        a = ch_two%n2spi1(ch_two%pps(ab))
!        b = ch_two%n2spi2(ch_two%pps(ab))
!        do ij = 1, ch_two%n_hh_state
!          i = ch_two%n2spi1(ch_two%hhs(ij))
!          j = ch_two%n2spi2(ch_two%hhs(ij))
!          denom = 1.d0 / get_denominator2(h,i,j,a,b)
!
!          v = v + v_pphh_1%two%GetTwBME(i,j,a,b,J2) * &
!              &   v_pphh_2%two%GetTwBME(a,b,i,j,J2) * denom
!        end do
!      end do
!      !$omp end do
!      !$omp end parallel
!      vsum = vsum + v * dble(2*J2+1)
!    end do
!    call v_pphh_1%fin()
!    call v_pphh_2%fin()
!    r = vsum
!    call timer%Add("Fourth order MBPT F7",omp_get_wtime()-ti)
!  end function energy_fourth_F7

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

    ch_bra => opch%ch_bra
    ch_ket => opch%ch_ket
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
      if( oa%GetHoleParticle() /= 1 ) cycle
      lowest_particle = min(lowest_particle, h%one%GetOBME(a,a))
    end do

    highest_hole = -100.d0
    do a = 1, norbs
      oa => ms%sps%orb(a)
      if( oa%GetHoleParticle() /= 0 ) cycle
      highest_hole = max(highest_hole, h%one%GetOBME(a,a))
    end do
    this%energy_gap = lowest_particle - highest_hole

    ! -- check HF condition, max_denom1b should be 0
    max_denom1b = 0.d0
    do a = 1, norbs
      oa => ms%sps%orb(a)
      if( oa%GetHoleParticle() /= 1 ) cycle
      do i = 1, norbs
      oi => ms%sps%orb(i)
      if( oi%GetHoleParticle() /= 0 ) cycle
        max_denom1b = max(max_denom1b, h%one%GetOBME(i,a) / abs(get_denominator1(h,i,a)))
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
        if( ms%sps%orb(a)%GetHoleParticle() /= 1 ) cycle
        if( ms%sps%orb(b)%GetHoleParticle() /= 1 ) cycle
        do ij = 1, n
          i = ms%two%jpz(ch)%n2spi1(ij)
          j = ms%two%jpz(ch)%n2spi2(ij)
          if( ms%sps%orb(i)%GetHoleParticle() /= 0 ) cycle
          if( ms%sps%orb(j)%GetHoleParticle() /= 0 ) cycle
          max_val = max(max_val, &
              & abs(h%two%GetTwBME(i,j,a,b,J2)) / abs(get_denominator2(h,i,j,a,b)))
          !!$omp critical
          !write(*,*) ch,a,b,i,j,max_val
          !!$omp end critical
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

  subroutine CalcScalarCorr(this,hamil,opr,is_MBPT_full,EN_denominator_in,part_of_hamil)
    class(MBPTScalar), intent(inout) :: this
    type(Ops), intent(in) :: hamil, opr
    logical, intent(in) :: is_MBPT_full
    logical, intent(in), optional :: EN_denominator_in
    logical, intent(in), optional :: part_of_hamil
    type(MSpace), pointer :: ms
    integer :: jmax
    real(8) :: s_1

    write(*,'(a)') " Many-body perturbation calculation up to 2nd order (scalar)"
    write(*,*)

    if(present(EN_denominator_in)) EN_denominator = EN_denominator_in
    if(present(part_of_hamil)) this%part_of_hamil = part_of_hamil
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

    if(.not. EN_denominator) write(*,"(a)") "# Moller-Plesset (MP) denominator"
    if(EN_denominator)       write(*,"(a)") "# Epstein-Nesbet (EN) denominator"
    if(this%part_of_hamil) write(*,"(a)") "# Calculating as a part of Hamiltonian"

    jmax = 2*hamil%ms%sps%lmax+1
    call sixjs%init(1,jmax,.true., 1,jmax,.true., 1,jmax,.true.)

    ms => hamil%ms
    this%s_0 = opr%zero
    if(is_MBPT_full) then
      s_1 = scalar_first1(hamil,opr)
      write(*,'(a,f16.8)') "First order correction from f: ", s_1
      s_1 = s_1 + scalar_second1(hamil,opr)
      write(*,'(a,f16.8)') "Second order correction from f: ", s_1
      s_1 = s_1 + scalar_first(hamil,opr)
      this%s_1 = s_1
      if(this%part_of_hamil) this%s_1 = 0.5d0 * this%s_1
      write(*,'(a,f16.8)') "First order correction: ", this%s_1
    end if

    call this%scalar_second(hamil,opr,is_MBPT_full,this%part_of_hamil)
    write(*,'(a)') "# Second order corrections: "
    write(*,'(a,f16.8)') "# s1 p ladder  = ", this%s_2_s1p
    write(*,'(a,f16.8)') "# s1 h ladder  = ", this%s_2_s1h
    write(*,'(a,f16.8)') "# s1 ph bubble = ", this%s_2_s1ph
    if(is_MBPT_full) then
      write(*,'(a,f16.8)') "# s2 pp ladder = ", this%s_2_s2pp
      write(*,'(a,f16.8)') "# s2 hh ladder = ", this%s_2_s2hh
      write(*,'(a,f16.8)') "# s2 ph ladder = ", this%s_2_s2ph
      write(*,'(a,f16.8)') "# v2 pp ladder = ", this%s_2_v2pp
      write(*,'(a,f16.8)') "# v2 hh ladder = ", this%s_2_v2hh
      write(*,'(a,f16.8)') "# v2 ph ladder = ", this%s_2_v2ph
    end if
    write(*,'(a,f16.8)') "# Total        = ", this%s_2

    call sixjs%fin()
    EN_denominator=.false.
    this%part_of_hamil=.false.
  end subroutine CalcScalarCorr

  function scalar_first1(h,s) result(r)
    ! a : particle
    ! i : hole
    !
    !    /\===========x
    !   /  \
    !   |  |
    ! a |  | i               x 2
    !   |  |
    !   \  /
    !    \/___________x
    !
    ! 2 \sum_{i,a} <i||a> <a||i> / denominator
    type(Ops), intent(in) :: h, s
    type(MSPace), pointer :: ms
    type(OneBodyChannel), pointer :: ch_one
    integer :: ch, a, i, J1, n, bra, ket
    type(Orbits), pointer :: sps
    type(SingleParticleOrbit), pointer :: oa, oi
    real(8) :: vsum, v, ti, r

    ti = omp_get_wtime()
    ms => h%ms
    sps => ms%sps

    vsum = 0.d0
    do bra = 1, ms%nh
      i = ms%holes(bra)
      oi => sps%GetOrbit(i)
      do ket = 1, ms%np
        a = ms%particles(ket)
        oa => sps%GetOrbit(a)
        if(oi%l /= oa%l) cycle
        if(oi%j /= oa%j) cycle
        if(oi%z /= oa%z) cycle
        vsum = vsum + dble(oi%j+1) * s%one%GetOBME(a,i) * h%one%GetOBME(a,i) / get_denominator1(h,i,a)
        !vsum = vsum + dble(oi%j+1) * s%one%GetOBME(a,i) * h%one%GetOBME(a,i) / (oi%e - oa%e) / ms%hw ! degenerate gap denominator
        !vsum = vsum + dble(oi%j+1) * s%one%GetOBME(a,i) * h%one%GetOBME(a,i) / (oi%e - oa%e) / 12.d0
      end do
    end do
    r = 2.d0 * vsum
    call timer%Add("First order MBPT for Scalar",omp_get_wtime()-ti)
  end function scalar_first1

  function scalar_first(h,s) result(r)
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
    type(Ops), intent(in) :: h, s
    type(MSPace), pointer :: ms
    type(TwoBodyChannel), pointer :: ch_two
    integer :: ch, ab, ij, J2, n
    integer :: a, b, i, j
    real(8) :: vsum, v, ti, r

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
        a = ch_two%n2spi1(ch_two%pp_s(ab))
        b = ch_two%n2spi2(ch_two%pp_s(ab))
        do ij = 1, ch_two%n_hh_state
          i = ch_two%n2spi1(ch_two%hh_s(ij))
          j = ch_two%n2spi2(ch_two%hh_s(ij))

          v = v + h%two%GetTwBME(i,j,a,b,J2) * &
              &   s%two%GetTwBME(a,b,i,j,J2) / &
              &   get_denominator2(h,i,j,a,b)
        end do
      end do
      !$omp end do
      !$omp end parallel
      vsum = vsum + v * dble(2*J2+1)
    end do
    r = 2.d0 * vsum
    call timer%Add("First order MBPT for Scalar",omp_get_wtime()-ti)
  end function scalar_first

  subroutine scalar_second(this,h,s, is_MBPT_full, part_of_hamil)
    class(MBPTScalar), intent(inout) :: this
    type(Ops), intent(in) :: h, s
    logical, intent(in) :: is_MBPT_full
    logical, intent(in) :: part_of_hamil

    if(.not. part_of_hamil) then
      this%s_2_s1p = scalar_second_s1p( h,s)
      this%s_2_s1h = scalar_second_s1h( h,s)
      this%s_2_s1ph= scalar_second_s1ph(h,s)
    end if

    if(is_MBPT_full) then
      this%s_2_s2pp = scalar_second_s2pp(h,s)
      this%s_2_s2hh = scalar_second_s2hh(h,s)
      this%s_2_s2ph = scalar_second_s2ph(h,s)
      if(.not. part_of_hamil) then
        this%s_2_v2pp = scalar_second_v2pp(h,s)
        this%s_2_v2hh = scalar_second_v2hh(h,s)
        this%s_2_v2ph = scalar_second_v2ph(h,s)
      end if
    end if

    this%s_2 = this%s_2_s1h + this%s_2_s1p + this%s_2_s1ph + &
        & this%s_2_s2pp + this%s_2_s2hh + this%s_2_s2ph + &
        & this%s_2_v2pp + this%s_2_v2hh + this%s_2_v2ph
  end subroutine scalar_second

  function scalar_second1(h,s) result(r)
    !
    !       ______x           _______x
    !     /|                /|               /|=======x          /|=======x
    !    / |               / |              / |                 / |
    !   /  | i            /  | b           /  | b              /  | j
    !   |  |              |  |             |  |                |  |
    ! a |  |======x     i |  |=======x   i |  |_______x 2    a |  |_______x 2
    !   |  |              |  |             |  |                |  |
    !   \  | j            \  | a           \  | a              \  | i
    !    \ |               \ |              \ |                 \ |
    !     \|______x         \|_______x       \|_______x          \|_______x
    !
    type(Ops), intent(in) :: h, s
    type(MSPace), pointer :: ms
    type(OneBodyChannel), pointer :: ch_one
    integer :: ch, a, b, i, j, ia, ib, ii, ij, J1, n, bra, ket
    type(Orbits), pointer :: sps
    type(SingleParticleOrbit), pointer :: oa, ob, oi, oj
    real(8) :: spp, shh, r, ti

    ti = omp_get_wtime()
    ms => h%ms
    sps => ms%sps
    spp = 0.d0; shh = 0.d0
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
              & (      h%one%GetOBME(ii,ia) * s%one%GetOBME(ia,ib) * h%one%GetOBME(ib,ii) + &
              & 2.d0 * s%one%GetOBME(ii,ia) * h%one%GetOBME(ia,ib) * h%one%GetOBME(ib,ii)) / &
              & ( get_denominator1(h,ii,ia) * get_denominator1(h,ii,ib) )
          !spp = spp + &
          !    & dble(oi%j+1) * &
          !    & (      h%one%GetOBME(ii,ia) * s%one%GetOBME(ia,ib) * h%one%GetOBME(ib,ii) + &
          !    & 2.d0 * s%one%GetOBME(ii,ia) * h%one%GetOBME(ia,ib) * h%one%GetOBME(ib,ii)) / &
          !    & ( (oi%e-oa%e) * (oi%e-ob%e) * 12.d0 )
        end do
      end do
    end do

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
              & (      h%one%GetOBME(ia,ii) * s%one%GetOBME(ii,ij) * h%one%GetOBME(ij,ia) + &
              & 2.d0 * s%one%GetOBME(ia,ii) * h%one%GetOBME(ii,ij) * h%one%GetOBME(ij,ia)) / &
              & ( get_denominator1(h,ii,ia) * get_denominator1(h,ij,ia) )
          !shh = shh - &
          !    & dble(oi%j+1) * &
          !    & (      h%one%GetOBME(ia,ii) * s%one%GetOBME(ii,ij) * h%one%GetOBME(ij,ia) + &
          !    & 2.d0 * s%one%GetOBME(ia,ii) * h%one%GetOBME(ii,ij) * h%one%GetOBME(ij,ia)) / &
          !    & ( (oi%e-oa%e) * (oj%e-oa%e) * 12.d0 )
        end do
      end do
    end do
    r = spp + shh
    call timer%Add("Second order MBPT for Scalar",omp_get_wtime()-ti)
  end function scalar_second1

  function scalar_second_s1p(h,s) result(r)
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
    type(Ops), intent(in) :: h, s
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
            norm = 1.d0
            if(a==b) norm = norm * sqrt(2.d0)
            if(a==c) norm = norm * sqrt(2.d0)
            v = v + norm * dble(2*J2+1) * tmp%m(bra,ket) * s%one%GetOBME(b,c) * dble(iph)
          end do
        end do
      end do
      vsum = vsum + v
      call tmp%fin()
      call m1%fin()
      call m2%fin()
    end do
    r = vsum
    call timer%Add("Second order MBPT s1 p ladder",omp_get_wtime()-ti)
  end function scalar_second_s1p

  function scalar_second_s1h(h,s) result(r)
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
    ! - <ab||ij> <ik||ab> <j|x|k> / 2 denominator
    type(Ops), intent(in) :: h, s
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
            norm = 1.d0
            if(i==j) norm = norm * sqrt(2.d0)
            if(i==k) norm = norm * sqrt(2.d0)
            v = v - norm * dble(2*J2+1) * tmp%m(bra,ket) * s%one%GetOBME(j,k) * dble(iph)
          end do
        end do
      end do
      vsum = vsum + v
      call tmp%fin()
      call m1%fin()
      call m2%fin()
    end do
    r = vsum
    call timer%Add("Second order MBPT s1 h ladder",omp_get_wtime()-ti)
  end function scalar_second_s1h

  function scalar_second_s1ph(h,s) result(r)
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
    type(Ops), intent(in) :: h, s
    type(MSPace), pointer :: ms
    integer :: ch, J2, bra, ket, iph
    type(TwoBodyChannel), pointer :: ch_two
    type(TwoBodyPartChannel) :: tmp_phpp, tmp_phhh
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
      call tmp_phhh%init(ch_two, ch_two)
      m1 = get_pphh_part(h%two%MatCh(ch,ch), h)
      m2 = h%two%MatCh(ch,ch)%get_phpp(ms%sps)
      if(m2%n_row < 1 .or. m2%n_col < 1) cycle
      call tmp_phhh%set_phhh(ms%sps, m2*m1)

      m2 = h%two%MatCh(ch,ch)%get_phhh(ms%sps)
      if(m2%n_row < 1 .or. m2%n_col < 1) cycle
      call tmp_phpp%set_phpp(ms%sps, m2*m1%t())

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
                & tmp_phhh%m(bra,ket) * s%one%GetOBME(c,i) / &
                & get_denominator1(h,i,c)
          end do
        end do
      end do

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
            v = v - norm * dble(2*J2+1) * dble(iph) * &
                & tmp_phpp%m(bra,ket) * s%one%GetOBME(a,k) / &
                & get_denominator1(h,k,a)
          end do
        end do
      end do
      vsum = vsum + v
      call tmp_phpp%fin()
      call tmp_phhh%fin()
      call m1%fin()
      call m2%fin()
    end do
    r = vsum
    call timer%Add("Second order MBPT s1 ph bubble",omp_get_wtime()-ti)
  end function scalar_second_s1ph

  function scalar_second_s2pp(h,s) result(r)
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
    type(Ops), intent(in) :: h, s
    type(MSpace), pointer :: ms
    type(TwoBodyChannel), pointer :: ch_two
    type(DMat) :: m1, m2, m3
    integer :: ch, J2, i
    real(8) :: v, vsum, ti, r

    ti = omp_get_wtime()
    ms => h%ms
    vsum = 0.d0
    do ch = 1, ms%two%NChan

      ch_two => ms%two%jpz(ch)
      J2 = ch_two%j
      if(ch_two%n_hh_state < 1) cycle
      if(ch_two%n_pp_state < 1) cycle

      m1 = get_pphh_part(h%two%MatCh(ch,ch), h)
      m2 = s%two%MatCH(ch,ch)%get_pppp(ms%sps)
      m3 = m1%t() * (m2 * m1)

      v = 0.d0
      do i = 1, size(m3%m, 1)
        v = v + m3%m(i,i)
      end do
      vsum = vsum + v * dble(2*J2+1)
    end do
    r = vsum
    call timer%Add("Second order MBPT s2 pp ladder",omp_get_wtime()-ti)
  end function scalar_second_s2pp

  function scalar_second_s2hh(h,s) result(r)
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
    type(Ops), intent(in) :: h, s
    type(MSpace), pointer :: ms
    type(TwoBodyChannel), pointer :: ch_two
    type(DMat) :: m1, m2, m3
    integer :: ch, J2, i
    real(8) :: v, vsum, ti, r

    ti = omp_get_wtime()
    ms => h%ms

    vsum = 0.d0
    do ch = 1, ms%two%NChan

      ch_two => ms%two%jpz(ch)
      J2 = ch_two%j
      if(ch_two%n_hh_state < 1) cycle
      if(ch_two%n_pp_state < 1) cycle

      m1 = get_pphh_part(h%two%MatCh(ch,ch), h)
      m2 = s%two%MatCH(ch,ch)%get_hhhh(ms%sps)
      m3 = m1 * (m2 * m1%t())

      v = 0.d0
      do i = 1, size(m3%m, 1)
        v = v + m3%m(i,i)
      end do
      vsum = vsum + v * dble(2*J2+1)
    end do
    r = vsum

    call timer%Add("Second order MBPT s2 hh ladder",omp_get_wtime()-ti)
  end function scalar_second_s2hh

  function scalar_second_s2ph(h,s) result(r)
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
    type(Ops), intent(in) :: h, s
    type(MSpace), pointer :: ms
    integer :: ch, J2, i
    type(CrossCoupledTwoBodyChannel), pointer :: ch_cc
    type(DMat) :: m1, m2, m3
    real(8) :: v, vsum, ti, r

    ti = omp_get_wtime()
    ms => h%ms

    vsum = 0.d0
    do ch = 1,ms%cc_two%NChan
      ch_cc => ms%cc_two%jpz(ch)
      J2 = ch_cc%j
      m1 = get_xc_pphh_to_phph(h%two, ch_cc, h)
      if(m1%n_row<1 .or. m1%n_col<1) cycle
      m2 = s%two%get_xc_hphp2phph(ch_cc)
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
    call timer%Add("Second order MBPT s2 ph ladder",omp_get_wtime()-ti)
  end function scalar_second_s2ph

  function scalar_second_v2pp(h,s) result(r)
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
    type(Ops), intent(in) :: h, s
    type(MSpace), pointer :: ms
    type(TwoBodyChannel), pointer :: ch_two
    type(DMat) :: m1, m2, m3
    integer :: ch, J2, i
    real(8) :: v, vsum, ti, r

    ti = omp_get_wtime()
    ms => h%ms
    vsum = 0.d0
    do ch = 1, ms%two%NChan

      ch_two => ms%two%jpz(ch)
      J2 = ch_two%j
      if(ch_two%n_hh_state < 1) cycle
      if(ch_two%n_pp_state < 1) cycle

      m1 = get_pphh_part(h%two%MatCh(ch,ch), h)
      m2 = h%two%MatCH(ch,ch)%get_pppp(ms%sps)
      m3 = get_pphh_part(s%two%MatCh(ch,ch), h)
      m3 = m3%t() * (m2 * m1)

      v = 0.d0
      do i = 1, size(m3%m, 1)
        v = v + m3%m(i,i)
      end do
      vsum = vsum + v * dble(2*J2+1)
    end do
    r = 2.d0 * vsum
    call timer%Add("Second order MBPT v2 pp ladder",omp_get_wtime()-ti)
  end function scalar_second_v2pp

  function scalar_second_v2hh(h,s) result(r)
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
    type(Ops), intent(in) :: h, s
    type(MSpace), pointer :: ms
    type(TwoBodyChannel), pointer :: ch_two
    type(DMat) :: m1, m2, m3
    integer :: ch, J2, i
    real(8) :: v, vsum, ti, r

    ti = omp_get_wtime()
    ms => h%ms

    vsum = 0.d0
    do ch = 1, ms%two%NChan

      ch_two => ms%two%jpz(ch)
      J2 = ch_two%j
      if(ch_two%n_hh_state < 1) cycle
      if(ch_two%n_pp_state < 1) cycle

      m1 = get_pphh_part(h%two%MatCh(ch,ch), h)
      m2 = h%two%MatCH(ch,ch)%get_hhhh(ms%sps)
      m3 = get_pphh_part(s%two%MatCh(ch,ch), h)
      m3 = m1 * (m2 * m3%t())

      v = 0.d0
      do i = 1, size(m3%m, 1)
        v = v + m3%m(i,i)
      end do
      vsum = vsum + v * dble(2*J2+1)
    end do
    r = 2.d0 * vsum
    call timer%Add("Second order MBPT v2 hh ladder",omp_get_wtime()-ti)
  end function scalar_second_v2hh

  function scalar_second_v2ph(h,s) result(r)
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
    type(Ops), intent(in) :: h, s
    type(MSpace), pointer :: ms
    integer :: ch, J2, i
    type(CrossCoupledTwoBodyChannel), pointer :: ch_cc
    type(DMat) :: m1, m2, m3
    real(8) :: v, vsum, ti, r

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
      m3 = get_xc_pphh_to_phph(s%two, ch_cc, h)
      if(m3%n_row<1 .or. m3%n_col<1) cycle
      m3 = m3 * (m2 * m1%t())
      v = 0.d0
      do i = 1, m3%n_row
        v = v + m3%m(i,i)
      end do
      vsum = vsum + dble(2*J2+1) * v
      call m1%fin()
      call m2%fin()
      call m3%fin()
    end do

    r = - 2.d0 * vsum
    call timer%Add("Second order MBPT v2 ph ladder",omp_get_wtime()-ti)
  end function scalar_second_v2ph

  subroutine FinMBPTDMat(this)
    class(MBPTDMat), intent(inout) :: this
    call this%rho_HF%fin()
    call this%rho_HO%fin()
    call this%Occ_HF%fin()
    call this%Occ_HO%fin()
    call this%C_HO2HF%fin()
    call this%C_HF2NAT%fin()
    call this%C_HO2NAT%fin()
    EN_denominator = .false.
    this%ms => null()
  end subroutine FinMBPTDMat

  subroutine InitMBPTDMat(this, HF, hamil, EN_denominator_in)
    class(MBPTDMat), intent(inout) :: this
    type(HFSolver), intent(in) :: HF
    type(Ops), intent(in) :: hamil
    logical, intent(in), optional :: EN_denominator_in
    type(MSPace), pointer :: ms
    type(Orbits), pointer :: sps
    integer :: idx, a, b, ch
    integer, allocatable :: aa(:), bb(:)
    type(SingleParticleOrbit), pointer :: oa, ob
    real(8) :: me, ti
    ms => hamil%ms
    this%ms => hamil%ms
    sps => ms%sps
    if(present(EN_denominator_in)) EN_denominator = EN_denominator_in
    write(*,"(a)") "# Calculating density matrix w/ MBPT"
    if(.not. EN_denominator) write(*,"(a)") "# Moller-Plesset (MP) denominator"
    if(EN_denominator)       write(*,"(a)") "# Epstein-Nesbet (EN) denominator"

    call this%rho_HF%init(ms%one, .true., 'DenMat',  0, 1, 0)
    call this%rho_HO%init(ms%one, .true., 'DenMat',  0, 1, 0)
    call this%Occ_HF%init(ms%one, .true., 'Occupation',  0, 1, 0)
    call this%Occ_HO%init(ms%one, .true., 'Occupation',  0, 1, 0)
    call this%C_HF2NAT%init(ms%one, .true., 'UT',  0, 1, 0)
    call this%C_HO2NAT%init(ms%one, .true., 'UT',  0, 1, 0)

    this%C_HO2HF = HF%C

    allocate(aa(sps%norbs*(sps%norbs+1)/2))
    allocate(bb(sps%norbs*(sps%norbs+1)/2))
    idx = 0
    do a = 1, sps%norbs
    do b = 1, a
    idx = idx + 1
    aa(idx) = a
    bb(idx) = b
    end do
    end do
    ti = omp_get_wtime()
    !$omp parallel
    !$omp do private(idx, a, b, oa, ob, me) schedule(dynamic)
    do idx = 1, sps%norbs*(sps%norbs+1)/2
    a = aa(idx)
    b = bb(idx)
    oa => sps%GetOrbit(a)
    ob => sps%GetOrbit(b)
    if(oa%j /= ob%j) cycle
    if(oa%l /= ob%l) cycle
    if(oa%z /= ob%z) cycle
    me = 0.d0
    if(abs(oa%GetOccupation()) > 1.d-6 .and. abs(ob%GetOccupation()) > 1.d-6) me = me + density_matrix_element_hh(a, b, hamil)
    if(abs(oa%GetOccupation()) < 1.d-6 .and. abs(ob%GetOccupation()) > 1.d-6) me = me + density_matrix_element_ph(a, b, hamil)
    if(abs(oa%GetOccupation()) > 1.d-6 .and. abs(ob%GetOccupation()) < 1.d-6) me = me + density_matrix_element_ph(b, a, hamil)
    if(abs(oa%GetOccupation()) < 1.d-6 .and. abs(ob%GetOccupation()) < 1.d-6) me = me + density_matrix_element_pp(a, b, hamil)
    call this%rho_HF%SetOBME(a,b,me)
    call this%rho_HF%SetOBME(b,a,me)
    end do
    !$omp end do
    !$omp end parallel
    call timer%Add("Density matrix @ 2nd order MBPT",omp_get_wtime()-ti)
    deallocate(aa,bb)

    do ch = 1, this%rho_HF%one%NChan
      this%rho_HO%MatCh(ch,ch)%DMat = &
        & this%C_HO2HF%MatCh(ch,ch)%DMat * &
        & this%rho_HF%MatCh(ch,ch)%DMat * &
        & this%C_HO2HF%MatCh(Ch,ch)%DMat%t()
    end do

    call this%GetCoef()
  end subroutine InitMBPTDMat

  function density_matrix_element_pp(a, b, hamil) result(r)
    integer, intent(in) :: a, b
    type(Ops), intent(in) :: hamil
    type(MSPace), pointer :: ms
    type(Orbits), pointer :: sps
    real(8) :: r
    integer :: c, i, j, JJ, Jmin, Jmax
    real(8) :: e_acij, e_bcij
    real(8) :: norm_ac, norm_bc, norm_ij, tbme
    type(SingleParticleOrbit), pointer :: oa, ob, oc, oi, oj

    ms => hamil%ms
    sps => ms%sps
    oa => sps%GetOrbit(a)
    ob => sps%GetOrbit(b)
    r = 0.d0
    do c = 1, sps%norbs
      oc => sps%GetOrbit(c)
      if(abs(oc%GetOccupation()) > 1.d-6) cycle
      do i = 1, sps%norbs
        oi => sps%GetOrbit(i)
        if(abs(oi%GetOccupation()) < 1.d-6) cycle
        do j = 1, sps%norbs
          oj => sps%GetOrbit(j)
          if(abs(oj%GetOccupation()) < 1.d-6) cycle

          if((-1)**(oa%l+oc%l+oi%l+oj%l) == -1) cycle
          if(oa%z+oc%z /= oi%z+oj%z) cycle

          e_acij = get_denominator2(hamil,i,j,a,c)
          e_bcij = get_denominator2(hamil,i,j,b,c)

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
  end function density_matrix_element_pp

  function density_matrix_element_hh(i, j, hamil) result(r)
    integer, intent(in) :: i, j
    type(Ops), intent(in) :: hamil
    type(MSPace), pointer :: ms
    type(Orbits), pointer :: sps
    real(8) :: r
    integer :: a, b, k, JJ, Jmin, Jmax
    real(8) :: e_abik, e_abjk
    real(8) :: norm_ab, norm_ik, norm_jk, tbme
    type(SingleParticleOrbit), pointer :: oa, ob, oi, oj, ok
    ms => hamil%ms
    sps => ms%sps
    oi => sps%GetOrbit(i)
    oj => sps%GetOrbit(j)
    r = 0.d0
    do a = 1, sps%norbs
      oa => sps%GetOrbit(a)
      if(abs(oa%GetOccupation()) > 1.d-6) cycle

      do b = 1, sps%norbs
        ob => sps%GetOrbit(b)
        if(abs(ob%GetOccupation()) > 1.d-6) cycle

        do k = 1, sps%norbs
          ok => sps%GetOrbit(k)
          if(abs(ok%GetOccupation()) < 1.d-6) cycle

          if((-1)**(oa%l+ob%l+oi%l+ok%l) == -1) cycle
          if(oa%z + ob%z /= oi%z + ok%z) cycle

          e_abik = get_denominator2(hamil,i,k,a,b)
          e_abjk = get_denominator2(hamil,j,k,a,b)

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
    if(i == j) r = r + oi%GetOccupation()
  end function density_matrix_element_hh

  function density_matrix_element_ph(a, i, hamil) result(r)
    integer, intent(in) :: a, i
    type(Ops), intent(in) :: hamil
    type(MSPace), pointer :: ms
    type(Orbits), pointer :: sps
    real(8) :: r
    integer :: b, c, j, k, JJ, Jmin, Jmax
    real(8) :: e_ai, e_bcij, e_abjk
    real(8) :: norm_bc, norm_ij, norm_ab, norm_jk, tbme
    type(SingleParticleOrbit), pointer :: oa, ob, oc, oi, oj, ok
    ms => hamil%ms
    sps => ms%sps
    oa => sps%GetOrbit(a)
    oi => sps%GetOrbit(i)
    r = 0.d0
    do b = 1, sps%norbs
      ob => sps%GetOrbit(b)
      if(abs(ob%GetOccupation()) > 1.d-6) cycle
      do c = 1, sps%norbs
        oc => sps%GetOrbit(c)
        if(abs(oc%GetOccupation()) > 1.d-6) cycle
        do j = 1, sps%norbs
          oj => sps%GetOrbit(j)
          if(abs(oj%GetOccupation()) < 1.d-6) cycle
          if((-1)**(oa%l+oj%l+ob%l+oc%l) == -1) cycle
          if(oa%z+oj%z /= ob%z+oc%z) cycle
          e_ai = get_denominator1(hamil, i, a)
          e_bcij = get_denominator2(hamil, i, j, b, c)

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
      if(abs(ob%GetOccupation()) > 1.d-6) cycle
      do j = 1, sps%norbs
        oj => sps%GetOrbit(j)
        if(abs(oj%GetOccupation()) < 1.d-6) cycle
        do k = 1, sps%norbs
          ok => sps%GetOrbit(k)
          if(abs(ok%GetOccupation()) < 1.d-6) cycle

          if((-1)**(oa%l+ob%l+oj%l+ok%l) == -1) cycle
          if(oa%z + ob%z /= oj%z + ok%z) cycle

          e_ai = get_denominator1(hamil, i, a)
          e_abjk = get_denominator2(hamil, j, k, a, b)

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
  end function density_matrix_element_ph

  subroutine GetCoef(this)
    class(MBPTDMat), intent(inout) :: this
    type(EigenSolSymD) :: sol
    integer :: ch, i, m, jj, iz
    real(8) :: N, Z
    N = 0.d0
    Z = 0.d0
    do ch = 1, this%rho_HF%one%NChan
      jj = this%rho_HF%one%jpz(ch)%j
      iz = this%rho_HF%one%jpz(ch)%z

      ! for HF => NAT
      call sol%init(this%rho_HF%MatCh(ch,ch)%DMat)
      call sol%DiagSym(this%rho_HF%MatCh(ch,ch)%DMat)
      m = size(sol%eig%v)
      do i = 1, m
        this%C_HF2NAT%MatCh(ch,ch)%m(:,i) = sol%vec%m(:,m-i+1)
        this%Occ_HF%MatCh(ch,ch)%m(i,i) = sol%eig%v(m-i+1)
      end do
      !call this%rho_HF%MatCh(ch,ch)%DMat%prt("rho HF")
      !call this%Occ_HF%MatCh(ch,ch)%DMat%prt("Occupation number HF")
      if(iz == -1) Z = Z + sum(sol%eig%v) * dble(jj+1)
      if(iz ==  1) N = N + sum(sol%eig%v) * dble(jj+1)
      call sol%fin()

      ! for HO => NAT
      call sol%init(this%rho_HO%MatCh(ch,ch)%DMat)
      call sol%DiagSym(this%rho_HO%MatCh(ch,ch)%DMat)
      m = size(sol%eig%v)
      do i = 1, m
        this%C_HO2NAT%MatCh(ch,ch)%m(:,i) = sol%vec%m(:,m-i+1)
        this%Occ_HO%MatCh(ch,ch)%m(i,i) = sol%eig%v(m-i+1)
      end do
      !call this%rho_HO%MatCh(ch,ch)%DMat%prt("rho HO")
      !call this%Occ_HO%MatCh(ch,ch)%DMat%prt("Occupation number HO")
      !if(iz == -1) Z = Z + sum(sol%eig%v) * dble(jj+1)
      !if(iz ==  1) N = N + sum(sol%eig%v) * dble(jj+1)
      call sol%fin()

    end do
    write(*,"(a,i4,a,f6.2)") " Actual Z: ", this%ms%Z, ", Z from tr(rho): ", Z
    write(*,"(a,i4,a,f6.2)") " Actual N: ", this%ms%N, ", N from tr(rho): ", N
  end subroutine GetCoef

! used to check the implementation of the matrix multiplication
!  function energy_fourth_F1_(h) result(r)
!    ! p1, p2, p3, p4: particle
!    ! h1, h2, h3, h4: hole
!    !     _____________
!    !    /\           /\
!    !   /  \ h1   p2 /  \ h2
!    !p1 |  |         \  /
!    !   |__|__________\/
!    !   |  |
!    !   |  |
!    !   |  |__________
!    !p3 |  |          /\
!    !   |  | h3   p4 /  \ h4
!    !   \  /         \  /
!    !    \/___________\/
!    !
!    ! F1 = - (1/4) [J1][J2][K][L]
!    ! { jh1 jh2 J1 } { jh1 jp4 J2 } { jh2 jh2 K }
!    ! { jh2 jp3 K  } { jp4 jp3 K  } { jp4 jp4 L }
!    ! <p1p2:J1||h1h2:J1> <p3h2:J1||p1p2:J1>
!    ! <h1p4:J2||h3h4:J2> <h3h4:J2||p3p4:J2> / denominator
!    !
!    type(Ops), intent(in) :: h
!    type(MSPace), pointer :: ms
!    type(Orbits), pointer :: sps
!    integer :: p1, p2, p3, p4
!    integer :: h1, h2, h3, h4
!    type(SingleParticleOrbit), pointer :: op1, op2, op3, op4
!    type(SingleParticleOrbit), pointer :: oh1, oh2, oh3, oh4
!    integer :: J1min, J1max, J2min, J2max, Kmin, Kmax, Lmin, Lmax
!    integer :: J1, J2, K, L
!    real(8) :: delp12, delh12, delh34, delp34
!    real(8) :: v1, v2, v, vsum, ti, r
!
!    ti = omp_get_wtime()
!    ms => h%ms
!    sps => ms%sps
!    vsum = 0.d0
!    do p1 = 1, sps%norbs
!      op1 => sps%GetOrbit(p1)
!      if(abs(op1%occ) > 1.d-6) cycle
!      do p2 = 1, sps%norbs
!        op2 => sps%GetOrbit(p2)
!        if(abs(op2%occ) > 1.d-6) cycle
!        if(op1%e + op2%e > ms%e2max) cycle
!        do p3 = 1, sps%norbs
!          op3 => sps%GetOrbit(p3)
!          if(abs(op3%occ) > 1.d-6) cycle
!          do p4 = 1, sps%norbs
!            op4 => sps%GetOrbit(p4)
!            if(abs(op4%occ) > 1.d-6) cycle
!            if(op3%e + op4%e > ms%e2max) cycle
!
!            do h1 = 1, sps%norbs
!              oh1 => sps%GetOrbit(h1)
!              if(abs(oh1%occ) < 1.d-6) cycle
!              if(oh1%e + op4%e > ms%e2max) cycle
!              do h2 = 1, sps%norbs
!                oh2 => sps%GetOrbit(h2)
!                if(abs(oh2%occ) < 1.d-6) cycle
!                if(oh1%e + oh2%e > ms%e2max) cycle
!                if(op3%e + oh2%e > ms%e2max) cycle
!                if(mod(op1%l + op2%l + oh1%l + oh2%l, 2)==1) cycle
!                if(mod(op3%l + oh2%l + op1%l + op2%l, 2)==1) cycle
!                if(op1%z + op2%z /= oh1%z + oh2%z) cycle
!                if(op3%z + oh2%z /= op1%z + op2%z) cycle
!                do h3 = 1, sps%norbs
!                  oh3 => sps%GetOrbit(h3)
!                  if(abs(oh3%occ) < 1.d-6) cycle
!                  do h4 = 1, sps%norbs
!                    oh4 => sps%GetOrbit(h4)
!                    if(abs(oh4%occ) < 1.d-6) cycle
!                    if(oh3%e + oh4%e > ms%e2max) cycle
!                    if(mod(oh1%l + op4%l + oh3%l + oh4%l, 2)==1) cycle
!                    if(mod(op3%l + op4%l + oh3%l + oh4%l, 2)==1) cycle
!                    if(oh1%z + op4%z /= oh3%z + oh4%z) cycle
!                    if(op3%z + op4%z /= oh3%z + oh4%z) cycle
!
!                    J1min = max(abs(op1%j-op2%j), abs(oh1%j-oh2%j), abs(op3%j-oh2%j))/2
!                    J1max = min(   (op1%j+op2%j),    (oh1%j+oh2%j),    (op3%j+oh2%j))/2
!                    J2min = max(abs(oh1%j-op4%j), abs(oh3%j-oh4%j), abs(op3%j-op4%j))/2
!                    J2max = min(   (oh1%j+op4%j),    (oh3%j+oh4%j),    (op3%j+op4%j))/2
!                    Kmin = max(abs(oh1%j-op3%j), 0)/2
!                    Kmax = min(   (oh1%j+op3%j), 2*oh2%j, 2*op4%j)/2
!                    Lmin = abs(oh2%j-op4%j)/2
!                    Lmax = (oh2%j+op4%j)/2
!                    if(J1min>J1max) cycle
!                    if(J2min>J2max) cycle
!                    if(Kmin>Kmax) cycle
!                    if(Lmin>Lmax) cycle
!
!                    v = 0.d0
!                    do L = Lmin, Lmax
!                      do K = Kmin, Kmax
!
!                        v1 = 0.d0
!                        v2 = 0.d0
!                        do J1 = J1min, J1max
!                          if(p1==p2 .and. mod(J1,2)==1) cycle
!                          if(h1==h2 .and. mod(J1,2)==1) cycle
!                          v1 = v1 + h%two%GetTwBME(p1,p2,h1,h2,J1) * &
!                              &     h%two%GetTwBME(p3,h2,p1,p2,J1) * &
!                              &     sixjs%get(oh1%j,oh2%j,2*J1,oh2%j,op3%j,2*K) * &
!                              &     dble(2*J1+1) * (-1.d0)**J1
!                        end do
!
!                        do J2 = J2min, J2max
!                          if(h3==h4 .and. mod(J2,2)==1) cycle
!                          if(p3==p4 .and. mod(J2,2)==1) cycle
!                          v2 = v2 + h%two%GetTwBME(h1,p4,h3,h4,J2) * &
!                              &     h%two%GetTwBME(h3,h4,p3,p4,J2) * &
!                              &     sixjs%get(oh1%j,op4%j,2*J2,op4%j,op3%j,2*K) * &
!                              &     dble(2*J2+1) * (-1.d0)**J2
!                        end do
!                        v = v + v1 * v2 * dble(2*K+1)*dble(2*L+1)* &
!                            &   sixjs%get(oh2%j,oh2%j,2*K,op4%j,op4%j,2*L) * (-1.d0)**(K+L)
!                      end do
!                    end do
!                    delp12 = 1.d0
!                    delh12 = 1.d0
!                    delh34 = 1.d0
!                    delp34 = 1.d0
!                    if(p1==p2) delp12 = 2.d0
!                    if(h1==h2) delh12 = sqrt(2.d0)
!                    if(p3==p4) delp34 = sqrt(2.d0)
!                    if(h3==h4) delh34 = 2.d0
!                    vsum = vsum + v * delp12*delh12*delp34*delh34 / ( &
!                        & denom2b(h%one,h1,h2,p1,p2) * &
!                        & denom1b(h%one,h1,p3) * denom2b(h%one,h3,h4,p3,p4))
!                  end do
!                end do
!              end do
!            end do
!
!          end do
!        end do
!      end do
!    end do
!    r = 0.25d0 * vsum
!    call timer%Add("Fourth order MBPT F1",omp_get_wtime()-ti)
!  end function energy_fourth_F1_
!
!  function energy_fourth_F3_(h) result(r)
!    ! p1, p2, p3, p4, p5: particle
!    ! h1, h2, h3:         hole
!    !     _____________
!    !    /\           /\
!    !   /  \ h1   p2 /  \ h2
!    !p1 |  |         \  /
!    !   |__|__________\/
!    !   |  |
!    !p3 |  |
!    !   |__|__________
!    !   |  |          /\
!    !p4 |  |      p5 /  \ h3
!    !   \  /         \  /
!    !    \/___________\/
!    !
!    ! F3 = - (1/4) [J1][J2][K][L] (-1)^(J1+J2+K+L)
!    ! { jp3 jh2 J1 } { jp3 jh3 J2 } { jh2 jh2 K }
!    ! { jh2 jh1 K  } { jh3 jh1 K  } { jh3 jh3 L }
!    ! <p1p2:J1||h1h2:J1> <p3h2:J1||p1p2:J1>
!    ! <p4p5:J2||p3h3:J2> <h1h3:J2||p4p5:J2> / denominator
!    !
!    type(Ops), intent(in) :: h
!    type(MSPace), pointer :: ms
!    type(Orbits), pointer :: sps
!    integer :: p1, p2, p3, p4, p5
!    integer :: h1, h2, h3
!    type(SingleParticleOrbit), pointer :: op1, op2, op3, op4, op5
!    type(SingleParticleOrbit), pointer :: oh1, oh2, oh3
!    integer :: J1min, J1max, J2min, J2max, Kmin, Kmax, Lmin, Lmax
!    integer :: J1, J2, K, L
!    real(8) :: delp12, delh12, delp45, delh13
!    real(8) :: v1, v2, v, vsum, ti, r
!
!
!    ti = omp_get_wtime()
!    ms => h%ms
!    sps => ms%sps
!    vsum = 0.d0
!    do p1 = 1, sps%norbs
!      op1 => sps%GetOrbit(p1)
!      if(abs(op1%occ) > 1.d-6) cycle
!      do p2 = 1, sps%norbs
!        op2 => sps%GetOrbit(p2)
!        if(abs(op2%occ) > 1.d-6) cycle
!        do p3 = 1, sps%norbs
!          op3 => sps%GetOrbit(p3)
!          if(abs(op3%occ) > 1.d-6) cycle
!          do p4 = 1, sps%norbs
!            op4 => sps%GetOrbit(p4)
!            if(abs(op4%occ) > 1.d-6) cycle
!            do p5 = 1, sps%norbs
!              op5 => sps%GetOrbit(p5)
!              if(abs(op5%occ) > 1.d-6) cycle
!
!              do h1 = 1, sps%norbs
!                oh1 => sps%GetOrbit(h1)
!                if(abs(oh1%occ) < 1.d-6) cycle
!                do h2 = 1, sps%norbs
!                  oh2 => sps%GetOrbit(h2)
!                  if(abs(oh2%occ) < 1.d-6) cycle
!                  do h3 = 1, sps%norbs
!                    oh3 => sps%GetOrbit(h3)
!                    if(abs(oh3%occ) < 1.d-6) cycle
!
!                    if((-1)**(op1%l+op2%l+oh1%l+oh2%l) == -1) cycle
!                    if((-1)**(op3%l+oh2%l+op1%l+op2%l) == -1) cycle
!                    if((-1)**(op4%l+op5%l+op3%l+oh3%l) == -1) cycle
!                    if((-1)**(oh1%l+oh3%l+op4%l+op5%l) == -1) cycle
!                    if(op1%z+op2%z /= oh1%z+oh2%z) cycle
!                    if(op3%z+oh2%z /= op1%z+op2%z) cycle
!                    if(op4%z+op5%z /= op3%z+oh3%z) cycle
!                    if(oh1%z+oh3%z /= op4%z+op5%z) cycle
!                    J1min = max(abs(op1%j-op2%j), abs(oh1%j-oh2%j), abs(op3%j-oh2%j))/2
!                    J1max = min(   (op1%j+op2%j),    (oh1%j+oh2%j),    (op3%j+oh2%j))/2
!                    J2min = max(abs(oh3%j-op3%j), abs(op4%j-op5%j), abs(oh1%j-oh3%j))/2
!                    J2max = min(   (oh3%j+op3%j),    (op4%j+op5%j),    (oh1%j+oh3%j))/2
!                    Kmin = max(abs(op3%j-oh1%j), 0)/2
!                    Kmax = min( (op3%j+oh1%j), 2*oh2%j, 2*oh3%j)/2
!                    Lmin = abs(oh2%j-oh3%j)/2
!                    Lmax = (oh2%j+oh3%j)/2
!                    if(J1min>J1max) cycle
!                    if(J2min>J2max) cycle
!                    if(Kmin>Kmax) cycle
!                    if(Lmin>Lmax) cycle
!
!                    v = 0.d0
!                    do L = Lmin, Lmax
!                      do K = Kmin, Kmax
!
!                        v1 = 0.d0
!                        v2 = 0.d0
!                        do J1 = J1min, J1max
!                          if(p1==p2 .and. mod(J1,2)==1) cycle
!                          if(h1==h2 .and. mod(J1,2)==1) cycle
!                          v1 = v1 + h%two%GetTwBME(p1,p2,h1,h2,J1) * &
!                              &     h%two%GetTwBME(p3,h2,p1,p2,J1) * &
!                              &     sixjs%get(op3%j,oh2%j,2*J1,oh2%j,oh1%j,2*K) * &
!                              &     dble(2*J1+1) * (-1.d0)**J1
!                        end do
!
!                        do J2 = J2min, J2max
!                          if(h1==h3 .and. mod(J2,2)==1) cycle
!                          if(p4==p5 .and. mod(J2,2)==1) cycle
!                          v2 = v2 + h%two%GetTwBME(p4,p5,p3,h3,J2) * &
!                              &     h%two%GetTwBME(h1,h3,p4,p5,J2) * &
!                              &     sixjs%get(op3%j,oh3%j,2*J2,oh3%j,oh1%j,2*K) * &
!                              &     dble(2*J2+1) * (-1.d0)**J2
!                        end do
!                        v = v + v1 * v2 * dble(2*K+1) * dble(2*L+1) * &
!                            &   sixjs%get(oh2%j,oh2%j,2*K,oh3%j,oh3%j,2*L) * (-1.d0)**(K+L)
!                      end do
!                    end do
!                    delh12 = 1.d0
!                    delh13 = 1.d0
!                    delp12 = 1.d0
!                    delp45 = 1.d0
!                    if(h1==h2) delh12 = sqrt(2.d0)
!                    if(h1==h3) delh13 = sqrt(2.d0)
!                    if(p1==p2) delp12 = 2.d0
!                    if(p4==p5) delp45 = 2.d0
!                    vsum = vsum + v * delp12*delh12*delh13*delp45 / ( &
!                        & denom2b(h%one,h1,h2,p1,p2) * &
!                        & denom1b(h%one,h1,p3) * denom2b(h%one,h1,h3,p4,p5))
!                  end do
!                end do
!              end do
!            end do
!          end do
!
!        end do
!      end do
!    end do
!    r = -0.25d0 * vsum
!    call timer%Add("Fourth order MBPT F3",omp_get_wtime()-ti)
!  end function energy_fourth_F3_
!
!  function energy_fourth_F4_(h) result(r)
!    ! p1, p2, p3:         particle
!    ! h1, h2, h3, h4, h5: hole
!    !     _____________
!    !    /\           /\
!    !   /  \ h1   p2 /  \ h2
!    !p1 |  |         \  /
!    !   |  |__________\/
!    !   |  |
!    !   |  | h3
!    !   |  |__________
!    !   |  |          /\
!    !   |  | h4   p3 /  \ h5
!    !   \  /         \  /
!    !    \/___________\/
!    !
!    ! F4 = - (1/4) [J1][J2][K][L] (-1)^(J1+J2+K+L)
!    ! { jp1 jp2 J1 } { jp1 jp3 J2 } { jp2 jp2 K }
!    ! { jp2 jh3 K  } { jp3 jh3 K  } { jp3 jp3 L }
!    ! <p1p2:J1||h1h2:J1> <h1h2:J1||h3p2:J1>
!    ! <h3p3:J2||h4h5:J2> <h4h5:J2||p1p3:J2> / denominator
!    !
!    type(Ops), intent(in) :: h
!    type(MSPace), pointer :: ms
!    type(Orbits), pointer :: sps
!    integer :: p1, p2, p3
!    integer :: h1, h2, h3, h4, h5
!    type(SingleParticleOrbit), pointer :: op1, op2, op3
!    type(SingleParticleOrbit), pointer :: oh1, oh2, oh3, oh4, oh5
!    integer :: J1min, J1max, J2min, J2max, Kmin, Kmax, Lmin, Lmax
!    integer :: J1, J2, K, L
!    real(8) :: delp12, delh12, delh45, delp13
!    real(8) :: v1, v2, v, vsum, ti, r
!
!
!    ti = omp_get_wtime()
!    ms => h%ms
!    sps => ms%sps
!    vsum = 0.d0
!    do p1 = 1, sps%norbs
!      do p2 = 1, sps%norbs
!        do p3 = 1, sps%norbs
!          op1 => sps%GetOrbit(p1)
!          op2 => sps%GetOrbit(p2)
!          op3 => sps%GetOrbit(p3)
!          if(abs(op1%occ) > 1.d-6) cycle
!          if(abs(op2%occ) > 1.d-6) cycle
!          if(abs(op3%occ) > 1.d-6) cycle
!
!          do h1 = 1, sps%norbs
!            do h2 = 1, sps%norbs
!              do h3 = 1, sps%norbs
!                do h4 = 1, sps%norbs
!                  do h5 = 1, sps%norbs
!                    oh1 => sps%GetOrbit(h1)
!                    oh2 => sps%GetOrbit(h2)
!                    oh3 => sps%GetOrbit(h3)
!                    oh4 => sps%GetOrbit(h4)
!                    oh5 => sps%GetOrbit(h5)
!                    if(abs(oh1%occ) < 1.d-6) cycle
!                    if(abs(oh2%occ) < 1.d-6) cycle
!                    if(abs(oh3%occ) < 1.d-6) cycle
!                    if(abs(oh4%occ) < 1.d-6) cycle
!                    if(abs(oh5%occ) < 1.d-6) cycle
!
!                    if((-1)**(op1%l+op2%l+oh1%l+oh2%l) == -1) cycle
!                    if((-1)**(oh3%l+op2%l+oh1%l+oh2%l) == -1) cycle
!                    if((-1)**(oh4%l+oh5%l+op3%l+oh3%l) == -1) cycle
!                    if((-1)**(op1%l+op3%l+oh4%l+oh5%l) == -1) cycle
!                    if(op1%z+op2%z /= oh1%z+oh2%z) cycle
!                    if(oh3%z+op2%z /= oh1%z+oh2%z) cycle
!                    if(oh4%z+oh5%z /= op3%z+oh3%z) cycle
!                    if(op1%z+op3%z /= oh4%z+oh5%z) cycle
!                    J1min = max(abs(op1%j-op2%j), abs(oh1%j-oh2%j), abs(oh3%j-op2%j))/2
!                    J1max = min(   (op1%j+op2%j),    (oh1%j+oh2%j),    (oh3%j+op2%j))/2
!                    J2min = max(abs(oh3%j-op3%j), abs(oh4%j-oh5%j), abs(op1%j-op3%j))/2
!                    J2max = min(   (oh3%j+op3%j),    (oh4%j+oh5%j),    (op1%j+op3%j))/2
!                    Kmin = max(abs(oh3%j-op1%j), 0)/2
!                    Kmax = min( (oh3%j+op1%j), 2*op2%j, 2*op3%j)/2
!                    Lmin = abs(op2%j-op3%j)/2
!                    Lmax = (op2%j+op3%j)/2
!                    if(J1min>J1max) cycle
!                    if(J2min>J2max) cycle
!                    if(Kmin>Kmax) cycle
!                    if(Lmin>Lmax) cycle
!
!                    v = 0.d0
!                    do L = Lmin, Lmax
!                      do K = Kmin, Kmax
!
!                        v1 = 0.d0
!                        do J1 = J1min, J1max
!                          if(p1==p2 .and. mod(J1,2)==1) cycle
!                          if(h1==h2 .and. mod(J1,2)==1) cycle
!                          v1 = v1 + h%two%GetTwBME(p1,p2,h1,h2,J1) * &
!                              &     h%two%GetTwBME(h1,h2,h3,p2,J1) * &
!                              &     sixjs%get(op1%j,op2%j,2*J1,op2%j,oh3%j,2*K) * &
!                              &     dble(2*J1+1) * (-1.d0)**J1
!                        end do
!
!                        v2 = 0.d0
!                        do J2 = J2min, J2max
!                          if(p1==p3 .and. mod(J2,2)==1) cycle
!                          if(h4==h5 .and. mod(J2,2)==1) cycle
!                          v2 = v2 + h%two%GetTwBME(h4,h5,h3,p3,J2) * &
!                              &     h%two%GetTwBME(p1,p3,h4,h5,J2) * &
!                              &     sixjs%get(op1%j,op3%j,2*J2,op3%j,oh3%j,2*K) * &
!                              &     dble(2*J2+1) * (-1.d0)**J2
!                        end do
!                        v = v + v1 * v2 * dble(2*K+1) * dble(2*L+1) * &
!                            &   sixjs%get(op2%j,op2%j,2*K,op3%j,op3%j,2*L) * (-1.d0)**(K+L)
!                      end do
!                    end do
!                    delp12 = 1.d0
!                    delp13 = 1.d0
!                    delh12 = 1.d0
!                    delh45 = 1.d0
!                    if(p1==p2) delp12 = sqrt(2.d0)
!                    if(p1==p3) delp13 = sqrt(2.d0)
!                    if(h1==h2) delh12 = 2.d0
!                    if(h4==h5) delh45 = 2.d0
!                    vsum = vsum + v * delp12*delh12*delp13*delh45 / ( &
!                        & denom2b(h%one,h1,h2,p1,p2) * &
!                        & denom1b(h%one,h3,p1) * denom2b(h%one,h4,h5,p1,p3))
!                  end do
!                end do
!              end do
!            end do
!          end do
!
!        end do
!      end do
!    end do
!    r = -0.25d0 * vsum
!    call timer%Add("Fourth order MBPT F4",omp_get_wtime()-ti)
!  end function energy_fourth_F4_
!
!  function energy_fourth_F5_(h) result(r)
!    ! p1, p2, p3, p4: particle
!    ! h1, h2, h3, h4: hole
!    !     _____________
!    !    /\           /\
!    !p1 /  \ h1   p2 /  \ h2
!    !   |  |         |  |
!    !   |__|_________|  |
!    !   |  |         |  |
!    !p3 |  |      p4 |  |
!    !   |  |_________|__|
!    !   |  |         |  |
!    !   |  | h3      |  | h4
!    !   \  /         \  /
!    !    \/___________\/
!    !
!    ! F5 = (1/16) [J]
!    ! <h1h2:J||p1p2:J> <p1p2:J||p3p4:J>
!    ! <p3p4:J||h3h4:J> <h3h4:J||h1h2:J> / denominator
!    !
!    type(Ops), intent(in) :: h
!    type(MSPace), pointer :: ms
!    type(Orbits), pointer :: sps
!    integer :: p1, p2, p3, p4
!    integer :: h1, h2, h3, h4
!    type(SingleParticleOrbit), pointer :: op1, op2, op3, op4
!    type(SingleParticleOrbit), pointer :: oh1, oh2, oh3, oh4
!    integer :: J, Jmin, Jmax
!    real(8) :: delp12, delh12, delh34, delp34
!    real(8) :: v, vsum, ti, r
!
!
!    ti = omp_get_wtime()
!    ms => h%ms
!    sps => ms%sps
!    vsum = 0.d0
!    do p1 = 1, sps%norbs
!      do p2 = 1, sps%norbs
!        do p3 = 1, sps%norbs
!          do p4 = 1, sps%norbs
!            op1 => sps%GetOrbit(p1)
!            op2 => sps%GetOrbit(p2)
!            op3 => sps%GetOrbit(p3)
!            op4 => sps%GetOrbit(p4)
!            if(abs(op1%occ) > 1.d-6) cycle
!            if(abs(op2%occ) > 1.d-6) cycle
!            if(abs(op3%occ) > 1.d-6) cycle
!            if(abs(op4%occ) > 1.d-6) cycle
!
!            do h1 = 1, sps%norbs
!              do h2 = 1, sps%norbs
!                do h3 = 1, sps%norbs
!                  do h4 = 1, sps%norbs
!                    oh1 => sps%GetOrbit(h1)
!                    oh2 => sps%GetOrbit(h2)
!                    oh3 => sps%GetOrbit(h3)
!                    oh4 => sps%GetOrbit(h4)
!                    if(abs(oh1%occ) < 1.d-6) cycle
!                    if(abs(oh2%occ) < 1.d-6) cycle
!                    if(abs(oh3%occ) < 1.d-6) cycle
!                    if(abs(oh4%occ) < 1.d-6) cycle
!
!                    if((-1)**(op1%l+op2%l+oh1%l+oh2%l) == -1) cycle
!                    if((-1)**(op1%l+op2%l+op3%l+op4%l) == -1) cycle
!                    if((-1)**(op3%l+op4%l+oh3%l+oh4%l) == -1) cycle
!                    if((-1)**(oh3%l+oh4%l+oh1%l+oh2%l) == -1) cycle
!                    if( op1%z+op2%z /= oh1%z+oh2%z) cycle
!                    if( op1%z+op2%z /= op3%z+op4%z) cycle
!                    if( op3%z+op4%z /= oh3%z+oh4%z) cycle
!                    if( oh3%z+oh4%z /= oh1%z+oh2%z) cycle
!                    Jmin = max(abs(op1%j-op2%j), abs(oh1%j-oh2%j), abs(op3%j-op4%j), abs(oh3%j-oh4%j))/2
!                    Jmax = min(   (op1%j+op2%j),    (oh1%j+oh2%j),    (op3%j+op4%j),    (oh3%j+oh4%j))/2
!                    if(Jmin>Jmax) cycle
!
!                    v = 0.d0
!                    do J = Jmin, Jmax
!                      if(p1==p2 .and. mod(J,2)==1) cycle
!                      if(p3==p4 .and. mod(J,2)==1) cycle
!                      if(h1==h2 .and. mod(J,2)==1) cycle
!                      if(h3==h4 .and. mod(J,2)==1) cycle
!                      v = v + dble(2*J+1) * &
!                          & h%two%GetTwBME(h1,h2,p1,p2,J) * &
!                          & h%two%GetTwBME(p1,p2,p3,p4,J) * &
!                          & h%two%GetTwBME(p3,p4,h3,h4,J) * &
!                          & h%two%GetTwBME(h3,h4,h1,h2,J)
!                    end do
!
!                    delp12 = 1.d0
!                    delp34 = 1.d0
!                    delh12 = 1.d0
!                    delh34 = 1.d0
!                    if(p1==p2) delp12 = 2.d0
!                    if(p3==p4) delp34 = 2.d0
!                    if(h1==h2) delh12 = 2.d0
!                    if(h3==h4) delh34 = 2.d0
!                    vsum = vsum + v * delp12*delh12*delp34*delh34 / ( &
!                        & get_denominator2(h,h1,h2,p1,p2) * &
!                        & get_denominator2(h,h3,h4,p3,p4) * &
!                        & get_denominator2(h,h1,h2,p3,p4))
!                  end do
!                end do
!              end do
!            end do
!
!          end do
!        end do
!      end do
!    end do
!    r =  vsum / 16.d0
!    call timer%Add("Fourth order MBPT F5",omp_get_wtime()-ti)
!  end function energy_fourth_F5_

  function energy_fourth_F7_(h) result(r)
    ! p1, p2, p3, p4: particle
    ! h1, h2, h3, h4: hole
    !     _____________
    !    /\           /\
    !p1 /  \ h1   p2 /  \ h2
    !   |  |         |  |
    !   |__|_________|__|
    !   |  |         |  |
    !p3 |  |         |  | h4
    !   |  |_________|  |
    !   |  |         |  |
    !   |  | h3   p4 |  |
    !   \  /         \  /
    !    \/___________\/
    !
    ! F7 = [K][J1][J2][J3][J4]
    ! { jh1 jh2 J1 } { jp3 jh2 J2 } { jh1 jp4 J3 } { jp3 jp4 J4 }
    ! { jp1 jp2 K  } { jp1 jh4 K  } { jh3 jp2 K  } { jh3 jh4 K  }
    ! <h1h2:J1||p1p2:J1> <p3h2:J2||p1h4:J2>
    ! <h1p4:J3||h3p2:J3> <p3p4:J4||h3h4:J4> / denominator
    !
    type(Ops), intent(in) :: h
    type(MSPace), pointer :: ms
    type(Orbits), pointer :: sps
    integer :: p1, p2, p3, p4
    integer :: h1, h2, h3, h4
    type(SingleParticleOrbit), pointer :: op1, op2, op3, op4
    type(SingleParticleOrbit), pointer :: oh1, oh2, oh3, oh4
    integer :: J1min, J1max, J1
    integer :: J2min, J2max, J2
    integer :: J3min, J3max, J3
    integer :: J4min, J4max, J4
    integer :: Kmin, Kmax, K
    real(8) :: v1, v2, v3, v4, r, ti, v, vsum, denom

    ti = omp_get_wtime()
    ms => h%ms
    sps => ms%sps
    vsum = 0.d0
    do p1 = 1, sps%norbs
      do p2 = 1, sps%norbs
        do p3 = 1, sps%norbs
          do p4 = 1, sps%norbs
            op1 => sps%GetOrbit(p1)
            op2 => sps%GetOrbit(p2)
            op3 => sps%GetOrbit(p3)
            op4 => sps%GetOrbit(p4)
            if(abs(op1%GetOccupation()) > 1.d-6) cycle
            if(abs(op2%GetOccupation()) > 1.d-6) cycle
            if(abs(op3%GetOccupation()) > 1.d-6) cycle
            if(abs(op4%GetOccupation()) > 1.d-6) cycle

            do h1 = 1, sps%norbs
              do h2 = 1, sps%norbs
                do h3 = 1, sps%norbs
                  do h4 = 1, sps%norbs
                    oh1 => sps%GetOrbit(h1)
                    oh2 => sps%GetOrbit(h2)
                    oh3 => sps%GetOrbit(h3)
                    oh4 => sps%GetOrbit(h4)
                    if(abs(oh1%GetOccupation()) < 1.d-6) cycle
                    if(abs(oh2%GetOccupation()) < 1.d-6) cycle
                    if(abs(oh3%GetOccupation()) < 1.d-6) cycle
                    if(abs(oh4%GetOccupation()) < 1.d-6) cycle

                    if((-1)**(op1%l+op2%l+oh1%l+oh2%l) == -1) cycle
                    if((-1)**(op3%l+oh2%l+op1%l+oh4%l) == -1) cycle
                    if((-1)**(oh1%l+op4%l+oh3%l+op2%l) == -1) cycle
                    if((-1)**(op3%l+op4%l+oh3%l+oh4%l) == -1) cycle

                    if( op1%z+op2%z /= oh1%z+oh2%z ) cycle
                    if( op3%z+oh2%z /= op1%z+oh4%z ) cycle
                    if( oh1%z+op4%z /= oh3%z+op2%z ) cycle
                    if( op3%z+op4%z /= oh3%z+oh4%z ) cycle

                    J1min = max(abs(op1%j-op2%j), abs(oh1%j-oh2%j))/2
                    J1max = min(   (op1%j+op2%j),    (oh1%j+oh2%j))/2
                    J2min = max(abs(op3%j-oh2%j), abs(op1%j-oh4%j))/2
                    J2max = min(   (op3%j+oh2%j),    (op1%j+oh4%j))/2
                    J3min = max(abs(oh1%j-op4%j), abs(oh3%j-op2%j))/2
                    J3max = min(   (oh1%j+op4%j),    (oh3%j+op2%j))/2
                    J4min = max(abs(op3%j-op4%j), abs(oh3%j-oh4%j))/2
                    J4max = min(   (op3%j+op4%j),    (oh3%j+oh4%j))/2
                    Kmin = max(abs(oh1%j-op2%j), abs(oh2%j-op1%j), abs(op3%j-oh4%j), abs(op4%j-oh3%j))/2
                    Kmax = min(   (oh1%j+op2%j),    (oh2%j+op1%j),    (op3%j+oh4%j),    (op4%j+oh3%j))/2
                    if(J1min > J1max) cycle
                    if(J2min > J2max) cycle
                    if(J3min > J3max) cycle
                    if(J4min > J4max) cycle
                    if(Kmin > Kmax) cycle

                    do K = Kmin, Kmax
                      v1 = 0.d0
                      do J1 = J1min, J1max
                        if(p1==p2 .and. mod(J1,2)==1) cycle
                        if(h1==h2 .and. mod(J1,2)==1) cycle
                        v1 = v1 + h%two%GetTwBME(p1,p2,h1,h2,J1) * &
                            & sixjs%get(op1%j,op2%j,2*J1,oh1%j,oh2%j,2*K) * dble(2*J1+1)
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
    r = vsum
    call timer%Add("Fourth order MBPT F7",omp_get_wtime()-ti)
  end function energy_fourth_F7_
end module HFMBPT
