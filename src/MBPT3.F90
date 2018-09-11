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
              jmax = min(jp1 + jp2, jh1 + jh2)/2
              jmin = max(iabs(jp1 - jp2), iabs(jh1 - jh2))/2
              d = 0.25d0 / (eh1 + eh2 - ep1 - ep2)

              v = 0.d0
              do j = jmin, jmax
                if(p1 == p2 .and. mod(j, 2) == 1) cycle
                if(h1 == h2 .and. mod(j, 2) == 1) cycle
                ich = ms%two%jptz2n(j, pari, itz)
                bra = ms%two%jptz(ich)%labels2n(p1,p2)
                ket = ms%two%jptz(ich)%labels2n(h1,h2)
                v = v + dble(2 * j + 1) * hamil%two%jptz(ich)%m(bra,ket) ** 2
              end do
              this%e_2 = this%e_2 + v * d * (Del(p1,p2) * Del(h1,h2)) ** 2

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
      write(*,*) '# MBPT third order contributions:'
      write(*,*) '# hh ladder, pp ladder, ph '
      write(*,'(3f15.6)') vhh, vpp, vph

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

                  ehh = ehh + v * d * &
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
      integer :: a, b, c, i, j, k
      integer :: aa, bb, cc, ii, jj, kk
      integer :: icha, ichb, ichc, ichi, ichj, ichk
      integer :: jab, jabmin, jabmax
      integer :: jac, jacmin, jacmax
      integer :: jkb, jkbmin, jkbmax
      integer :: jtot, jmax, jmin, pha, ich
      integer :: bra, ket
      integer :: itzab, ipab
      integer :: itzac, ipac
      integer :: itzkb, ipkb
      real(8) :: delab, delij, delac, delkj
      real(8) :: crv1, crv2, crv3, vsum, deno

      eph = 0.d0
      do i = 1, sps%n
        if(sps%h_orbits(i) /= 1) cycle

        do j = 1, sps%n
          if(sps%h_orbits(j) /= 1) cycle
          delij = 1.d0
          if(i == j) delij = dsqrt(2.d0)

          do k = 1, sps%n
            if(sps%h_orbits(k) /= 1) cycle
            delkj = 1.d0
            if(j == k) delkj = dsqrt(2.d0)

            do a = 1, sps%n
              if(sps%h_orbits(a) == 1) cycle

              do b = 1, sps%n
                if(sps%h_orbits(b) == 1) cycle
                itzab = (sps%itz(a) + sps%itz(b)) / 2
                ipab = (-1) ** (sps%ll(a) + sps%ll(b))
                itzkb = (sps%itz(k) + sps%itz(b)) / 2
                ipkb = (-1) ** (sps%ll(k) + sps%ll(b))
                delab = 1.d0
                if(a == b) delab = dsqrt(2.d0)
                if(sps%itz(i) + sps%itz(j) /= 2*itzab) cycle
                if((-1)** (sps%ll(i) + sps%ll(j)) /= ipab) cycle

                do c = 1, sps%n
                  if(sps%h_orbits(c) == 1) cycle
                  itzac = (sps%itz(a) + sps%itz(c)) / 2
                  ipac = (-1) ** (sps%ll(a) + sps%ll(c))
                  if(itzkb /= (sps%itz(i) + sps%itz(c))/2) cycle
                  if(ipkb /= (-1) ** (sps%ll(i) + sps%ll(c))) cycle
                  if(itzac /= (sps%itz(k) + sps%itz(j))/2) cycle
                  if(ipac /= (-1) ** (sps%ll(k) + sps%ll(j))) cycle
                  delac = 1.d0
                  if(a == c) delac = dsqrt(2.d0)

                  if(sps%nshell(a) + sps%nshell(b) > params%e2max) cycle
                  if(sps%nshell(i) + sps%nshell(j) > params%e2max) cycle
                  if(sps%nshell(k) + sps%nshell(b) > params%e2max) cycle
                  if(sps%nshell(i) + sps%nshell(c) > params%e2max) cycle
                  if(sps%nshell(a) + sps%nshell(c) > params%e2max) cycle
                  if(sps%nshell(k) + sps%nshell(j) > params%e2max) cycle

                  jabmax = min(sps%jj(i) + sps%jj(j), &
                      &        sps%jj(a) + sps%jj(b)) / 2
                  jabmin = max(abs(sps%jj(i) - sps%jj(j)), &
                      &        abs(sps%jj(a) - sps%jj(b))) / 2

                  jkbmax = min(sps%jj(k) + sps%jj(b), &
                      &        sps%jj(i) + sps%jj(c)) / 2
                  jkbmin = max(abs(sps%jj(k) - sps%jj(b)), &
                      &        abs(sps%jj(i) - sps%jj(c))) / 2

                  jacmax = min(sps%jj(a) + sps%jj(c), &
                      &        sps%jj(k) + sps%jj(j)) / 2
                  jacmin = max(abs(sps%jj(a) - sps%jj(c)), &
                      &        abs(sps%jj(k) - sps%jj(j))) / 2

                  jmax = min(sps%jj(a) + sps%jj(j), &
                      &      sps%jj(b) + sps%jj(i), &
                      &      sps%jj(c) + sps%jj(k)) / 2
                  jmin = max(abs(sps%jj(a) - sps%jj(j)), &
                      &      abs(sps%jj(b) - sps%jj(i)), &
                      &      abs(sps%jj(c) - sps%jj(k))) / 2

                  vsum = 0.d0
                  do jtot = jmin, jmax
                    crv1 = 0.d0; crv2 = 0.d0; crv3 = 0.d0
                    do jab = jabmin, jabmax
                      if(a == b .and. mod(jab, 2) == 1) cycle
                      if(i == j .and. mod(jab, 2) == 1) cycle
                      ich = ms%two%jptz2n(jab, ipab, itzab)
                      bra = ms%two%jptz(ich)%labels2n(i,j)
                      ket = ms%two%jptz(ich)%labels2n(a,b)
                      pha = ms%two%jptz(ich)%iphase(i,j) * &
                          & ms%two%jptz(ich)%iphase(a,b)
                      if(bra * ket == 0) cycle
                      crv1 = crv1 + dble(2*jab+1) * &
                          &  sjs(sps%jj(i), sps%jj(j), 2*jab, &
                          &      sps%jj(a), sps%jj(b), 2*jtot) * &
                          &  hamil%two%jptz(ich)%m(bra,ket) * dble(pha)
                    end do

                    do jkb = jkbmin, jkbmax
                      ich = ms%two%jptz2n(jkb, ipkb, itzkb)
                      bra = ms%two%jptz(ich)%labels2n(k,b)
                      ket = ms%two%jptz(ich)%labels2n(i,c)
                      pha = ms%two%jptz(ich)%iphase(k,b) * &
                          & ms%two%jptz(ich)%iphase(i,c)
                      if(bra * ket == 0) cycle
                      crv2 = crv2 + dble(2*jkb+1) * &
                          &  sjs(sps%jj(k), sps%jj(b), 2*jkb, &
                          &      sps%jj(i), sps%jj(c), 2*jtot) * &
                          &  hamil%two%jptz(ich)%m(bra,ket) * dble(pha)
                    end do

                    do jac = jacmin, jacmax
                      if(a == c .and. mod(jac, 2) == 1) cycle
                      if(k == j .and. mod(jac, 2) == 1) cycle
                      ich = ms%two%jptz2n(jac, ipac, itzac)
                      bra = ms%two%jptz(ich)%labels2n(a,c)
                      ket = ms%two%jptz(ich)%labels2n(k,j)
                      pha = ms%two%jptz(ich)%iphase(a,c) * &
                          & ms%two%jptz(ich)%iphase(k,j)
                      if(bra * ket == 0) cycle
                      crv3 = crv3 + dble(2*jac+1) * &
                          &  sjs(sps%jj(a), sps%jj(c), 2*jac, &
                          &      sps%jj(k), sps%jj(j), 2*jtot) * &
                          &  hamil%two%jptz(ich)%m(bra,ket) * dble(pha)
                    end do
                    !vsum = vsum + dble(2*jtot+1) * crv1 * crv2 * crv3
                    vsum = vsum + dble(2*jtot+1) ** 2 * crv1 * crv2 * crv3
                  end do

                  icha = ms%one%label2jptz(a)
                  ichb = ms%one%label2jptz(b)
                  ichc = ms%one%label2jptz(c)
                  ichi = ms%one%label2jptz(i)
                  ichj = ms%one%label2jptz(j)
                  ichk = ms%one%label2jptz(k)
                  aa = ms%one%jptz(icha)%label2n(a)
                  bb = ms%one%jptz(ichb)%label2n(b)
                  cc = ms%one%jptz(ichc)%label2n(c)
                  ii = ms%one%jptz(ichi)%label2n(i)
                  jj = ms%one%jptz(ichj)%label2n(j)
                  kk = ms%one%jptz(ichk)%label2n(k)
                  deno = (hamil%one%jptz(icha)%m(aa,aa) + &
                      &   hamil%one%jptz(ichb)%m(bb,bb) - &
                      &   hamil%one%jptz(ichi)%m(ii,ii) - &
                      &   hamil%one%jptz(ichj)%m(jj,jj)) * &
                      &  (hamil%one%jptz(icha)%m(aa,aa) + &
                      &   hamil%one%jptz(ichc)%m(cc,cc) - &
                      &   hamil%one%jptz(ichj)%m(jj,jj) - &
                      &   hamil%one%jptz(ichk)%m(kk,kk))

                  eph = eph - vsum * delij * delab * delac * delkj / deno

                end do
              end do
            end do

          end do
        end do
      end do
    end function third_ph
  end subroutine CalcEnergyCorr

end module MBPT3
