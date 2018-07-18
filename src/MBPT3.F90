module MBPT3
  use common_library, only: Del, sjs
  use InputParameters, only: parameters
  use ModelSpace, only: spo_pn, MSpace, OneBodySpace, TwoBodySpace, ThreeBodySpace, &
      & OneBodyChannel, TwoBodyChannel, ThreeBodyChannel
  use ScalarOperator
  implicit none
  type :: MBPT
    real(8) :: e_0 = 0.d0, e_2 = 0.d0, e_3 = 0.d0
  contains
    procedure :: calc => CalcEnergyCorr
  end type MBPT
contains
  subroutine CalcEnergyCorr(this, params, sps, ms, hamil)
  class(MBPT), intent(inout) :: this
    type(parameters), intent(in) :: params
    type(spo_pn) :: sps
    type(MSpace) :: ms
    type(ScalarOperators), intent(in) :: hamil

    this%e_0 = hamil%zero
    call energy_second()
    call energy_third()
  contains

    subroutine energy_second()
      integer :: p1, p2, h1, h2
      integer :: lp1, lp2, lh1, lh2
      integer :: jp1, jp2, jh1, jh2
      integer :: zp1, zp2, zh1, zh2
      integer :: np1, np2, nh1, nh2
      integer :: j, jmax, jmin
      integer :: bra, ket, iphase
      integer :: pari, itz
      integer :: ich, ichp1, ichp2, ichh1, ichh2
      real(8) :: v, d
      real(8) :: ep1, ep2, eh1, eh2

      this%e_2 = 0.d0
      do p1 = 1, sps%n
        if(sps%h_orbits(p1) == 1) cycle
        jp1 = sps%jj(p1)
        lp1 = sps%ll(p1)
        zp1 = sps%itz(p1)
        ichp1 = ms%one%label2jptz(p1)
        np1 = ms%one%jptz(ichp1)%label2n(p1)
        ep1 = hamil%one%jptz(ichp1)%m(np1, np1)

        do p2 = 1, sps%n
          if(sps%h_orbits(p2) == 1) cycle
          jp2 = sps%jj(p2)
          lp2 = sps%ll(p2)
          zp2 = sps%itz(p2)
          ichp2 = ms%one%label2jptz(p2)
          np2 = ms%one%jptz(ichp2)%label2n(p2)
          ep2 = hamil%one%jptz(ichp2)%m(np2, np2)

          if(sps%nshell(p1) + sps%nshell(p2) > params%e2max) cycle

          do h1 = 1, sps%n
            if(sps%h_orbits(h1) /= 1) cycle
            jh1 = sps%jj(h1)
            lh1 = sps%ll(h1)
            zh1 = sps%itz(h1)
            ichh1 = ms%one%label2jptz(h1)
            nh1 = ms%one%jptz(ichh1)%label2n(h1)
            eh1 = hamil%one%jptz(ichh1)%m(nh1, nh1)

            do h2 = 1, sps%n
              if(sps%h_orbits(h2) /= 1) cycle
              jh2 = sps%jj(h2)
              lh2 = sps%ll(h2)
              zh2 = sps%itz(h2)
              ichh2 = ms%one%label2jptz(h2)
              nh2 = ms%one%jptz(ichh2)%label2n(h2)
              eh2 = hamil%one%jptz(ichh2)%m(nh2, nh2)

              if(sps%nshell(h1) + sps%nshell(h2) > params%e2max) cycle

              pari = (-1) ** (lp1 + lp2)
              if(pari /= (-1) ** (lh1 + lh2)) cycle
              itz = (zp1 + zp2)/2
              if(itz /= (zh1 + zh2)/2) cycle
              jmax = min(jp1 + jp2, jh1 + jh2)
              jmin = max(iabs(jp1 - jp2), iabs(jh1 - jh2))
              d = 1.d0 / (eh1 + eh2 - ep1 - ep2)

              v = 0.d0
              do j = jmin, jmax
                if(p1 == p2 .and. mod(j, 2) /= 0) cycle
                if(h1 == h2 .and. mod(j, 2) /= 0) cycle
                ich = ms%two%jptz2n(j, pari, itz)
                bra = ms%two%jptz(ich)%labels2n(p1, p2)
                ket = ms%two%jptz(ich)%labels2n(h1, h2)
                if(bra * ket == 0) cycle
                iphase = ms%two%jptz(ich)%iphase(p1, p2) * &
                  &      ms%two%jptz(ich)%iphase(h1, h2)
                v = v + dble(2 * j + 1) * hamil%two%jptz(ich)%m(bra,ket) ** 2
              end do

              this%e_2 = this%e_2 + 0.25d0 * v * d * (Del(p1,p2) * Del(h1,h2)) ** 2

            end do
          end do
        end do
      end do

    end subroutine energy_second

    subroutine energy_third()
      real(8) :: vpp, vhh, vph
      vhh = third_hh()
      vpp = third_pp()
      vph = third_ph()

      this%e_3 = vhh + vpp + vph

    end subroutine energy_third

    function third_hh() result(ehh)
      real(8) :: ehh
      integer :: p1, p2, h1, h2, h3, h4
      integer :: np1, np2, nh1, nh2, nh3, nh4
      integer :: jmax, jmin, j
      integer :: nh12, nh34, np12
      integer :: phh12, phh34, php12
      integer :: pari, itz, ich
      integer :: ichp1, ichp2, ichh1, ichh2, ichh3, ichh4
      real(8) :: v, d

      ehh = 0.d0
      do h1 = 1, sps%n
        if(sps%h_orbits(h1) /= 1) cycle
        do h2 = 1, sps%n
          if(sps%h_orbits(h2) /= 1) cycle
          pari = (-1) ** (sps%ll(h1) + sps%ll(h2))
          itz = (sps%itz(h1) + sps%itz(h2)) / 2
          if(sps%nshell(h1) + sps%nshell(h2) > params%e2max) cycle

          do h3 = 1, sps%n
            if(sps%h_orbits(h3) /= 1) cycle
            do h4 = 1, sps%n
              if(sps%h_orbits(h4) /= 1) cycle

              if(itz /= (sps%itz(h3) + sps%itz(h4)) / 2) cycle
              if(pari /= (-1) ** (sps%ll(h3) + sps%ll(h4)) ) cycle
              if(sps%nshell(h3) + sps%nshell(h4) > params%e2max) cycle

              do p1 = 1, sps%n
                if(sps%h_orbits(p1) == 1) cycle
                do p2 = 1, sps%n
                  if(sps%h_orbits(p2) == 1) cycle

                  if(itz /= (sps%itz(p1) + sps%itz(p2)) / 2) cycle
                  if(pari /= (-1) ** (sps%ll(p1) + sps%ll(p2))) cycle
                  if(sps%nshell(p1) + sps%nshell(p2) > params%e2max) cycle


                  jmax = min(sps%jj(h1) + sps%jj(h2), &
                    &        sps%jj(h3) + sps%jj(h4), &
                    &        sps%jj(p1) + sps%jj(p2)) / 2
                  jmin = max(abs(sps%jj(h1) - sps%jj(h2)), &
                    &        abs(sps%jj(h3) - sps%jj(h4)), &
                    &        abs(sps%jj(p1) - sps%jj(p2))) / 2

                  v = 0.d0
                  do j = jmin, jmax
                    if(h1 == h2 .and. mod(j, 2) == 1) cycle
                    if(h3 == h4 .and. mod(j, 2) == 1) cycle
                    if(p1 == p2 .and. mod(j, 2) == 1) cycle
                    ich = ms%two%jptz2n(j, pari, itz)
                    nh12 = ms%two%jptz(ich)%labels2n(h1, h2)
                    nh34 = ms%two%jptz(ich)%labels2n(h3, h4)
                    np12 = ms%two%jptz(ich)%labels2n(p1, p2)
                    phh12 = ms%two%jptz(ich)%iphase(h1, h2)
                    phh34 = ms%two%jptz(ich)%iphase(h3, h4)
                    php12 = ms%two%jptz(ich)%iphase(p1, p2)

                    v = v + dble(2 * j + 1) * &
                      & hamil%two%jptz(ich)%m(nh12, nh34) * dble(phh12 * phh34) * &
                      & hamil%two%jptz(ich)%m(nh34, np12) * dble(phh34 * php12) * &
                      & hamil%two%jptz(ich)%m(np12, nh12) * dble(php12 * phh12)

                  end do
                  ichh1 = ms%one%label2jptz(h1)
                  ichh2 = ms%one%label2jptz(h2)
                  ichh3 = ms%one%label2jptz(h3)
                  ichh4 = ms%one%label2jptz(h4)
                  ichp1 = ms%one%label2jptz(p1)
                  ichp2 = ms%one%label2jptz(p2)

                  nh1 = ms%one%jptz(ichh1)%label2n(h1)
                  nh2 = ms%one%jptz(ichh2)%label2n(h2)
                  nh3 = ms%one%jptz(ichh3)%label2n(h3)
                  nh4 = ms%one%jptz(ichh4)%label2n(h4)
                  np1 = ms%one%jptz(ichp1)%label2n(p1)
                  np2 = ms%one%jptz(ichp2)%label2n(p2)

                  d = 1.d0 / ((hamil%one%jptz(ichh1)%m(nh1,nh1) + &
                    &  hamil%one%jptz(ichh2)%m(nh2,nh2) - &
                    &  hamil%one%jptz(ichp1)%m(np1,np1) - &
                    &  hamil%one%jptz(ichp2)%m(np2,np2)) * &
                    & (hamil%one%jptz(ichh3)%m(nh3,nh3) + &
                    &  hamil%one%jptz(ichh4)%m(nh4,nh4) - &
                    &  hamil%one%jptz(ichp1)%m(np1,np1) - &
                    &  hamil%one%jptz(ichp2)%m(np2,np2)))

                  ehh = ehh * v * d * &
                    & (Del(p1,p2) * Del(h1,h2) * Del(h3,h4)) ** 2 / 8.d0


                end do
              end do
            end do
          end do
        end do
      end do
    end function third_hh

    function third_pp() result(epp)
      real(8) :: epp

    integer :: p1, p2, p3, p4, h1, h2
    integer :: np1, np2, nh1, nh2, np3, np4
    integer :: jmax, jmin, j
    integer :: np12, np34, nh12
    integer :: phh12, php34, php12
    integer :: pari, itz, ich
    integer :: ichp1, ichp2, ichh1, ichh2, ichp3, ichp4
    real(8) :: v, d
    epp = 0.d0

    do h1 = 1, sps%n
      if(sps%h_orbits(h1) /= 1) cycle

      do h2 = 1, sps%n
        if(sps%h_orbits(h2) /= 1) cycle

        if(sps%nshell(h1) + sps%nshell(h2) > params%e2max) cycle
        pari = (-1) ** (sps%ll(h1) + sps%ll(h2))
        itz = (sps%itz(h1) + sps%itz(h2)) / 2

        do p1 = 1, sps%n
          if(sps%h_orbits(p1) == 1) cycle

          do p2 = 1, sps%n
            if(sps%h_orbits(p2) == 1) cycle

            if(sps%nshell(p1) + sps%nshell(p2) > params%e2max) cycle
            if(pari /= (-1) ** (sps%ll(p1) + sps%ll(p2))) cycle
            if(itz /= (sps%itz(p1) + sps%itz(p2)) / 2) cycle

            do p3 = 1, sps%n
              if(sps%h_orbits(p3) == 1) cycle

              do p4 = 1, sps%n
                if(sps%h_orbits(p4) == 1) cycle

                if(sps%nshell(p3) + sps%nshell(p4) > params%e2max) cycle
                if(pari /= (-1) ** (sps%ll(p3) + sps%ll(p4))) cycle
                if(itz /= (sps%itz(p3) + sps%itz(p4)) / 2) cycle

                jmax = min(sps%jj(h1) + sps%jj(h2), &
                  &        sps%jj(p3) + sps%jj(p4), &
                  &        sps%jj(p1) + sps%jj(p2)) / 2
                jmin = max(abs(sps%jj(h1) - sps%jj(h2)), &
                  &        abs(sps%jj(p3) - sps%jj(p4)), &
                  &        abs(sps%jj(p1) - sps%jj(p2))) / 2

                v = 0.d0
                do j = jmin, jmax
                  if(h1 == h2 .and. mod(j, 2) == 1) cycle
                  if(p3 == p4 .and. mod(j, 2) == 1) cycle
                  if(p1 == p2 .and. mod(j, 2) == 1) cycle
                  ich = ms%two%jptz2n(j, pari, itz)
                  nh12 = ms%two%jptz(ich)%labels2n(h1, h2)
                  np34 = ms%two%jptz(ich)%labels2n(p3, p4)
                  np12 = ms%two%jptz(ich)%labels2n(p1, p2)
                  phh12 = ms%two%jptz(ich)%iphase(h1, h2)
                  php34 = ms%two%jptz(ich)%iphase(p3, p4)
                  php12 = ms%two%jptz(ich)%iphase(p1, p2)

                  v = v + dble(2 * j + 1) * &
                    & hamil%two%jptz(ich)%m(np12, np34) * dble(php12 * php34) * &
                    & hamil%two%jptz(ich)%m(np34, nh12) * dble(php34 * phh12) * &
                    & hamil%two%jptz(ich)%m(nh12, np12) * dble(phh12 * php12)
                end do

                ichp1 = ms%one%label2jptz(p1)
                ichp2 = ms%one%label2jptz(p2)
                ichh1 = ms%one%label2jptz(h1)
                ichh2 = ms%one%label2jptz(h2)
                ichp3 = ms%one%label2jptz(p3)
                ichp4 = ms%one%label2jptz(p4)
                np1 = ms%one%jptz(ichp1)%label2n(p1)
                np2 = ms%one%jptz(ichp2)%label2n(p2)
                nh1 = ms%one%jptz(ichh1)%label2n(h1)
                nh2 = ms%one%jptz(ichh2)%label2n(h2)
                np3 = ms%one%jptz(ichp3)%label2n(p3)
                np4 = ms%one%jptz(ichp4)%label2n(p4)

                d = 1.d0 / ((hamil%one%jptz(ichh1)%m(nh1,nh1) + &
                  &     hamil%one%jptz(ichh2)%m(nh2,nh2) - &
                  &     hamil%one%jptz(ichp1)%m(np1,np1) - &
                  &     hamil%one%jptz(ichp2)%m(np2,np2)) * &
                  &    (hamil%one%jptz(ichh1)%m(nh1,nh1) + &
                  &     hamil%one%jptz(ichh2)%m(nh2,nh2) - &
                  &     hamil%one%jptz(ichp3)%m(np3,np3) - &
                  &     hamil%one%jptz(ichp4)%m(np4,np4)))

                epp = epp + v * d * &
                  & (Del(p1,p2) * Del(p3,p4) * Del(h1,h2)) ** 2 / 8.d0

              end do
            end do

          end do
        end do

      end do
    end do

    end function third_pp

    function third_ph() result(eph)
      real(8) :: eph
      integer :: p1, p2, p3, h1, h2, h3
      integer :: pp1, pp2, pp3, hh1, hh2, hh3
      integer :: ichp1, ichp2, ichp3, ichh1, ichh2, ichh3
      integer :: jh3p2, jh3p2min, jh3p2max
      integer :: jh1h3, jh1h3min, jh1h3max
      integer :: jh1h2, jh1h2min, jh1h2max
      integer :: j, jmax, jmin, pha, ich
      integer :: bra, ket
      integer :: itzh3p2, iparh3p2
      integer :: itzh1h3, iparh1h3
      integer :: itzh1h2, iparh1h2
      real(8) :: delh1h3, delp1p2, delp1p3, delh1h2
      real(8) :: crv1, crv2, crv3, vsum, deno

      eph = 0.d0
      do h1 = 1, sps%n
        if(sps%h_orbits(h1) /= 1) cycle

        do h2 = 1, sps%n
          if(sps%h_orbits(h2) /= 1) cycle
          itzh1h2 = (sps%itz(h1) + sps%itz(h2)) / 2
          iparh1h2 = (-1) ** (sps%ll(h1) + sps%ll(h2))
          delh1h2 = 1.d0
          if(h1 == h2) delh1h2 = dsqrt(2.d0)

          do h3 = 1, sps%n
            if(sps%h_orbits(h3) /= 1) cycle
            itzh1h3 = (sps%itz(h1) + sps%itz(h3)) / 2
            iparh1h3 = (-1) ** (sps%ll(h1) + sps%ll(h3))
            delh1h3 = 1.d0
            if(h1 == h3) delh1h3 = dsqrt(2.d0)

            do p1 = 1, sps%n
              if(sps%h_orbits(p1) == 1) cycle

              do p2 = 1, sps%n
                if(sps%h_orbits(p2) == 1) cycle
                itzh3p2 = (sps%itz(h3) + sps%itz(p2)) / 2
                iparh3p2 = (-1) ** (sps%ll(h3) + sps%ll(p2))
                delp1p2 = 1.d0
                if(p1 == p2) delp1p2 = dsqrt(2.d0)
                if((sps%itz(p1) + sps%itz(p2)) / 2 /= itzh1h2) cycle
                if((-1) ** (sps%ll(p1) + sps%ll(p2)) /= iparh1h2) cycle

                do p3 = 1, sps%n
                  if(sps%h_orbits(p3) == 1) cycle
                  if((sps%itz(p1) + sps%itz(p3)) / 2 /= itzh1h3) cycle
                  if((-1) ** (sps%ll(p1) + sps%ll(p3)) /= iparh1h3) cycle
                  if((sps%itz(h2) + sps%itz(p3)) / 2 /= itzh3p2) cycle
                  if((-1) ** (sps%ll(h2) + sps%ll(p3)) /= iparh3p2) cycle
                  delp1p3 = 1.d0
                  if(p1 == p3) delp1p3 = dsqrt(2.d0)
                  if(sps%nshell(p1) + sps%nshell(p2) > params%e2max) cycle
                  if(sps%nshell(h1) + sps%nshell(h2) > params%e2max) cycle
                  if(sps%nshell(p3) + sps%nshell(h2) > params%e2max) cycle
                  if(sps%nshell(h3) + sps%nshell(p2) > params%e2max) cycle
                  if(sps%nshell(p1) + sps%nshell(p3) > params%e2max) cycle
                  if(sps%nshell(h1) + sps%nshell(h3) > params%e2max) cycle

                  jh1h2max = min(sps%jj(h1) + sps%jj(h2), &
                    &            sps%jj(p1) + sps%jj(p2)) / 2
                  jh1h2min = max(abs(sps%jj(h1) - sps%jj(h2)), &
                    &            abs(sps%jj(p1) - sps%jj(p2))) / 2
                  jh1h3max = min(sps%jj(h1) + sps%jj(h3), &
                    &           sps%jj(p1) + sps%jj(p3)) / 2
                  jh1h3min = max(abs(sps%jj(h1) - sps%jj(h3)), &
                    &            abs(sps%jj(p1) - sps%jj(p3))) / 2
                  jh3p2max = min(sps%jj(h3) + sps%jj(p2), &
                    &            sps%jj(h2) + sps%jj(p3)) / 2
                  jh3p2min = max(abs(sps%jj(h3) - sps%jj(p2)), &
                    &            abs(sps%jj(h2) - sps%jj(p3))) / 2

                  jmax = min(sps%jj(p1) + sps%jj(h1), &
                    &        sps%jj(p2) + sps%jj(h2), &
                    &        sps%jj(p3) + sps%jj(h3)) / 2
                  jmin = max(abs(sps%jj(p1) + sps%jj(h1)), &
                    &        abs(sps%jj(p2) + sps%jj(h2)), &
                    &        abs(sps%jj(p2) + sps%jj(h3))) / 2

                  vsum = 0.d0
                  do j = jmin, jmax
                    crv1 = 0.d0; crv2 = 0.d0; crv3 = 0.d0
                    do jh3p2 = jh3p2min, jh3p2max
                      ich = ms%two%jptz2n(jh3p2, iparh3p2, itzh3p2)
                      bra = ms%two%jptz(ich)%labels2n(p3,h2)
                      ket = ms%two%jptz(ich)%labels2n(h3,p2)
                      pha = ms%two%jptz(ich)%iphase(p3,h2) * &
                        &   ms%two%jptz(ich)%iphase(h3,p2)
                      if(bra * ket == 0) cycle
                      crv1 = crv1 + dble(2 * jh3p2 + 1) * &
                        &    sjs(sps%jj(p3), sps%jj(h2), 2 * jh3p2, &
                        &        sps%jj(p2), sps%jj(h3), 2 * j) * &
                        &    hamil%two%jptz(ich)%m(bra, ket) * dble(pha) * (-1.d0) ** (jh3p2 - j)
                    end do

                    do jh1h3 = jh1h3min, jh1h3max
                      ich = ms%two%jptz2n(jh1h3, iparh1h3, itzh1h3)
                      bra = ms%two%jptz(ich)%labels2n(h1,h3)
                      ket = ms%two%jptz(ich)%labels2n(p1,p3)
                      pha = ms%two%jptz(ich)%iphase(h1,h3) * &
                        &   ms%two%jptz(ich)%iphase(p1,p3)
                      if(bra * ket == 0) cycle
                      crv2 = crv2 + dble(2 * jh1h3 + 1) * &
                        &    sjs(sps%jj(h1), sps%jj(h3), 2 * jh1h3, &
                        &        sps%jj(p3), sps%jj(p1), 2 * j) * &
                        &    hamil%two%jptz(ich)%m(bra, ket) * dble(pha) * (-1.d0) ** (jh1h3 - j)
                    end do

                    do jh1h2 = jh1h2min, jh1h2max
                      ich = ms%two%jptz2n(jh1h2, iparh1h2, itzh1h2)
                      bra = ms%two%jptz(ich)%labels2n(p1,p2)
                      ket = ms%two%jptz(ich)%labels2n(h1,h2)
                      pha = ms%two%jptz(ich)%iphase(p1,p2) * &
                        &   ms%two%jptz(ich)%iphase(h1,h2)
                      if(bra * ket == 0) cycle
                      crv3 = crv3 + dble(2 * jh1h2 + 1) * &
                        &    sjs(sps%jj(p1), sps%jj(p3), 2 * jh1h2, &
                        &        sps%jj(h1), sps%jj(h2), 2 * j) * &
                        &    hamil%two%jptz(ich)%m(bra, ket) * dble(pha) * (-1.d0) ** (jh1h2 - j)
                    end do
                    vsum = vsum + dble(2 * j + 1) * crv1 * crv2 * crv3

                  end do
                  ichp1 = ms%one%label2jptz(p1)
                  ichp2 = ms%one%label2jptz(p2)
                  ichp3 = ms%one%label2jptz(p3)
                  ichh1 = ms%one%label2jptz(h1)
                  ichh2 = ms%one%label2jptz(h2)
                  ichh3 = ms%one%label2jptz(h3)
                  pp1 = ms%one%jptz(ichp1)%label2n(p1)
                  pp2 = ms%one%jptz(ichp2)%label2n(p2)
                  pp3 = ms%one%jptz(ichp3)%label2n(p3)
                  hh1 = ms%one%jptz(ichh1)%label2n(h1)
                  hh2 = ms%one%jptz(ichh2)%label2n(h2)
                  hh3 = ms%one%jptz(ichh3)%label2n(h3)
                  deno = (hamil%one%jptz(ichh1)%m(hh1,hh1) + &
                    &     hamil%one%jptz(ichh2)%m(hh2,hh2) - &
                    &     hamil%one%jptz(ichp1)%m(pp1,pp1) - &
                    &     hamil%one%jptz(ichp2)%m(pp2,pp2)) * &
                    &    (hamil%one%jptz(ichh1)%m(hh1,hh1) + &
                    &     hamil%one%jptz(ichh3)%m(hh3,hh3) - &
                    &     hamil%one%jptz(ichp1)%m(pp1,pp1) - &
                    &     hamil%one%jptz(ichp3)%m(pp3,pp3))

                  eph = eph - vsum * delh1h2 * delh1h3 * delp1p2 * delp1p3 / deno

                end do
              end do
            end do
          end do
        end do
      end do
    end function third_ph

  end subroutine CalcEnergyCorr

end module MBPT3
