module DefineOperators
  use myfort
  implicit none
contains
  subroutine GetOperatorRank(optr, jr, pr, zr)
    character(*), intent(in) :: optr
    integer, intent(out) :: jr, pr, zr
    character(256), allocatable :: splt(:)
    character(:), allocatable :: str
    type(sys) :: s

    jr = 0
    pr = 0
    zr = 0

    if(s%find(optr, "file_")) then
      call s%split(optr, "file_", splt)
      str = splt(2)
      call s%split(str, "+", splt)
      if(size(splt) > 1) then
        read(splt(1),*) jr
        pr = 1
        read(splt(2),*) zr
        return
      end if

      call s%split(str, "-", splt)
      if(size(splt) > 1) then
        read(splt(1),*) jr
        pr = -1
        read(splt(2),*) zr
        return
      end if

      write(*,"(2a)") "operator from file: ", trim(optr)
    end if

    select case(optr)
    case('hamil', 'Hamil', "HOHamil", "Hohamil", 'hohamil' ,&
          & 'Hcm','HCM','RM2', 'Rm2', "rm2", &
          & 'Tcm', 'tcm', 'Rp2', 'RP2', 'rp2', 'Rn2', 'RN2', 'rn2',"DenMat",&
          & "Tkin", "tkin", 'Rp2so', 'rp2so', 'Rp2_1b', 'Rp2_2b', 'Rp4_LO', 'Rp4_NLO')
      jr = 0
      pr = 1
      zr = 0
    case default
      write(*,'(2a)') 'Unknown operator: ', trim(optr)
    end select
  end subroutine GetOperatorRank

  recursive function one_body_element(optr, ia, ib, hw, A, Z, N) result(r)
    real(8) :: r
    character(*), intent(in) :: optr
    real(8), intent(in) :: hw
    real(8) :: amnucl
    integer, intent(in) :: ia(4), ib(4), A, Z, N
    integer :: na, la, ja, za
    integer :: nb, lb, jb, zb

    r = 0.d0
    na = ia(1)
    la = ia(2)
    ja = ia(3)
    za = ia(4)

    nb = ib(1)
    lb = ib(2)
    jb = ib(3)
    zb = ib(4)

    select case(optr)
    case("Hamil", "hamil", 'kinetic', 'Kinetic')
      if(la /= lb) return
      if(ja /= jb) return
      if(za /= zb) return
      if(abs(na - nb) > 1) return
      if(na == nb) r = dble(2*na + la) + 1.5d0
      if(na == nb + 1) r = dsqrt(dble(na) * (dble(na+la)+0.5d0))
      if(na == nb - 1) r = dsqrt(dble(nb) * (dble(nb+lb)+0.5d0))
      r = r * 0.5d0 * hw
      return

    case("Tkin", "tkin")
      r = one_body_element("kinetic", ia, ib, hw, A, Z, N) - &
          & one_body_element("tcm", ia, ib, hw, A, Z, N)
      return

    case("Tcm", "tcm")
      if(la /= lb) return
      if(ja /= jb) return
      if(za /= zb) return
      if(abs(na - nb) > 1) return
      if(na == nb) r = dble(2*na + la) + 1.5d0
      if(na == nb + 1) r = dsqrt(dble(na) * (dble(na+la)+0.5d0))
      if(na == nb - 1) r = dsqrt(dble(nb) * (dble(nb+lb)+0.5d0))
      r = r * 0.5d0 * hw / dble(A)
      return

    case("HOHamil","Hohamil", "hohamil")
      if(na /= nb) return
      if(la /= lb) return
      if(ja /= jb) return
      if(za /= zb) return
      r = (dble(2*na+la)+1.5d0) * hw
      return

    case("Hcm","HCM")
      if(na /= nb) return
      if(la /= lb) return
      if(ja /= jb) return
      if(za /= zb) return
      r = (dble(2*na+la)+1.5d0) * hw / dble(A)
      return

    case("RM2","Rm2", "rm2")
      if(la /= lb) return
      if(ja /= jb) return
      if(za /= zb) return
      if(abs(na - nb) > 1) return
      if(na == nb) r = dble(2*na + la) + 1.5d0
      if(na == nb + 1) r = -dsqrt(dble(na) * (dble(na+la)+0.5d0))
      if(na == nb - 1) r = -dsqrt(dble(nb) * (dble(nb+lb)+0.5d0))
      amnucl = (amp + amn) * 0.5d0
      r = r * (1.d0/dble(A) - 1.d0/dble(A)**2) * hc ** 2 / (amnucl * hw)
      return

    case("Rp2","RP2", "rp2", "Rp2_1b")
      if(la /= lb) return
      if(ja /= jb) return
      if(za /= zb) return
      if(abs(na - nb) > 1) return
      if(na == nb) r = dble(2*na + la) + 1.5d0
      if(na == nb + 1) r = -dsqrt(dble(na) * (dble(na+la)+0.5d0))
      if(na == nb - 1) r = -dsqrt(dble(nb) * (dble(nb+lb)+0.5d0))
      if(za == -1) r = r * (1.d0 / dble(Z) - 2.d0 / dble(A*Z) + 1.d0 / dble(A)**2)
      if(za ==  1) r = r / dble(A)**2
      amnucl = (amp + amn) * 0.5d0
      r = r * hc ** 2 / (amnucl * hw)
      return

    case("Rp2_2b")
      return

    case("Rp2so", "rp2so")
      if(la /= lb) return
      if(ja /= jb) return
      if(za /= zb) return
      if(na /= nb) return
      if(ja == 2*la-1) r = dble(la)
      if(ja == 2*la+1) r =-dble(la+1)
      if(za == -1) r = (r+1.d0) * (gs+gv-1.d0)*0.5d0
      if(za ==  1) r = (r+1.d0) * (gs-gv) * 0.5d0
      r = (-1.d0) * r * hc**2 / ((amp+amn)*0.5d0)**2 / dble(Z)


    case("Rn2","RN2", "rn2")
      if(la /= lb) return
      if(ja /= jb) return
      if(za /= zb) return
      if(abs(na - nb) > 1) return
      if(na == nb) r = dble(2*na + la) + 1.5d0
      if(na == nb + 1) r = -dsqrt(dble(na) * (dble(na+la)+0.5d0))
      if(na == nb - 1) r = -dsqrt(dble(nb) * (dble(nb+lb)+0.5d0))
      if(za == -1) r = r / dble(A)**2
      if(za ==  1) r = r * (1.d0 / dble(N) - 2.d0 / dble(A*N) + 1.d0 / dble(A)**2)
      amnucl = (amp + amn) * 0.5d0
      r = r * hc ** 2 / (amnucl * hw)
      return

    case("Rp4_LO")
      if(la /= lb) return
      if(ja /= jb) return
      if(za /= zb) return
      r = radius_power(4, na, la, nb, lb, sqrt(2.d0 * hc**2 / ((amp + amn) * hw))) * (1.d0 - dble(za))
      r = 0.5d0 * r / dble(Z)

    case("Rp4_NLO")
      if(la /= lb) return
      if(ja /= jb) return
      if(za /= zb) return
      r = radius_power(4, na, la, nb, lb, sqrt(2.d0 * hc**2 / ((amp + amn) * hw))) * (1.d0 - dble(za))
      r = 0.5d0 * r / dble(Z) * (1.d0 - 4.d0 / dble(A))

    case default
      write(*,'(3a)') "In one_body_element, ", &
          & trim(optr), " is not implemented"
      return
    end select
  end function one_body_element

  recursive function two_body_element(optr, ia, ib, ic, id, Jab, Jcd, hw, A, Z, N) result(r)
    real(8) :: r
    character(*), intent(in) :: optr
    real(8), intent(in) :: hw
    real(8) :: amnucl
    integer, intent(in) :: ia(4), ib(4), ic(4), id(4), A, Z, N, Jab, Jcd
    integer :: Zab, Zcd, Pab, Pcd

    Pab = (-1)**(ia(2)+ib(2))
    Pcd = (-1)**(ic(2)+id(2))
    Zab = (ia(4)+ib(4))/2
    Zcd = (ic(4)+id(4))/2

    r = 0.d0
    select case(optr)
    case("Hamil", "hamil")
      return

    case('Tcm', 'tcm')
      if(Jab /= Jcd .or. Pab /= Pcd .or. Zab /= Zcd) then
        write(*,'(a,6i3)') "Error in SetTwoBodyChannel: ", Jab, Jcd, Pab, Pcd, Zab, Zcd
        return
      end if
      r = p_dot_p(ia,ib,ic,id,Jab) * hw / dble(A)
      return

    case("Tkin", "tkin")
      r = - two_body_element("tcm", ia, ib, ic, id, Jab, Jcd, hw, A, Z, N)
      return

    case("Hcm","HCM")
      if(Jab /= Jcd .or. Pab /= Pcd .or. Zab /= Zcd) then
        write(*,'(a,2i3)') "Error in SetTwoBodyChannel: ", Jab, Jcd
        return
      end if
      r = (p_dot_p(ia,ib,ic,id,Jab)+r_dot_r(ia,ib,ic,id,Jab)) * hw/dble(A)
      return

    case('RM2', 'Rm2', 'rm2')
      if(Jab /= Jcd .or. Pab /= Pcd .or. Zab /= Zcd) then
        write(*,'(a,2i3)') "Error in SetTwoBodyChannel: ", Jab, Jcd
        return
      end if
      amnucl = (amp + amn) * 0.5d0
      r = - 2.d0 * r_dot_r(ia,ib,ic,id,Jab) * hc**2 / (amnucl*hw*dble(A)**2)
      return

    case('Rp2', 'RP2', 'rp2', 'Rp2_2b')
      if(Jab /= Jcd .or. Pab /= Pcd .or. Zab /= Zcd) then
        write(*,'(a,2i3)') "Error in SetTwoBodyChannel: ", Jab, Jcd
        return
      end if
      if(Zab == -1) r = ( 2.d0 / dble(A)**2 - 4.d0 / dble(A*Z) ) * r_dot_r(ia,ib,ic,id,Jab)
      if(Zab ==  0) r = ( 2.d0 / dble(A)**2 - 2.d0 / dble(A*Z) ) * r_dot_r(ia,ib,ic,id,Jab)
      if(Zab ==  1) r =  2.d0 * r_dot_r(ia,ib,ic,id,Jab) / dble(A)**2
      amnucl = (amp + amn) * 0.5d0
      r = r * hc**2 / (amnucl*hw)
      return

    case("Rp2_1b")
      return
    case("Rp2so", "rp2so")
      return

    case('Rn2', 'RN2', 'rn2')
      if(Jab /= Jcd .or. Pab /= Pcd .or. Zab /= Zcd) then
        write(*,'(a,2i3)') "Error in SetTwoBodyChannel: ", Jab, Jcd
        return
      end if
      if(Zab == -1) r = 2.d0 * r_dot_r(ia,ib,ic,id,Jab) / dble(A)**2
      if(Zab ==  0) r = ( 2.d0 / dble(A)**2 - 2.d0 / dble(A*N) ) * r_dot_r(ia,ib,ic,id,Jab)
      if(Zab ==  1) r = ( 2.d0 / dble(A)**2 - 4.d0 / dble(A*N) ) * r_dot_r(ia,ib,ic,id,Jab)
      amnucl = (amp + amn) * 0.5d0
      r = r * hc**2 / (amnucl*hw)
      return

    case("Rp4_LO")
      return

    case("Rp4_NLO")
      if(Jab /= Jcd .or. Pab /= Pcd .or. Zab /= Zcd) then
        write(*,'(a,2i3)') "Error in SetTwoBodyChannel: ", Jab, Jcd
        return
      end if
      r = matel_r4(ia, ib, ic, id, Jab, 0, hw) - matel_r4(ia, ib, ic, id, Jab, 1, hw)
      r = r * (-2.d0) / dble(Z * A)

    case default
      write(*,'(3a)') "In two_body_element, ", &
          & trim(optr), " is not implemented"
      return
    end select
  end function two_body_element

  !
  !  (1/2) * m * omega**2 * r_{i} \cdot r_{j} / hw
  !
  function r_dot_r(a, b, c, d, J) result(r)
    integer, intent(in) :: a(4), b(4), c(4), d(4), J
    real(8) :: r
    integer :: ja, jb, jc, jd
    integer :: ea, eb, ec, ed
    real(8) :: fact

    r = 0.d0
    ea = 2*a(1) + a(2)
    eb = 2*b(1) + b(2)
    ec = 2*c(1) + c(2)
    ed = 2*d(1) + d(2)
    if(abs(ea+eb-ec-ed) > 2 .or. mod(abs(ea+eb-ec-ed),2)==1) return

    ja = a(3)
    jb = b(3)
    jc = c(3)
    jd = d(3)

    fact = 1.d0
    if(a(1) == b(1) .and. a(2) == b(2) .and. a(3) == b(3) .and. a(4) == b(4)) fact = fact / dsqrt(2.d0)
    if(c(1) == d(1) .and. c(2) == d(2) .and. c(3) == d(3) .and. c(4) == d(4)) fact = fact / dsqrt(2.d0)

    r = (-1.d0) ** ((jb+jc)/2+J) * &
        & sjs(ja,jb,2*J,jd,jc,2) * red_r_j(a,c) * red_r_j(b,d) + &
        & (-1.d0) ** ((jb+jc)/2) * &
        & sjs(ja,jb,2*J,jc,jd,2) * red_r_j(a,d) * red_r_j(b,c)

    r = r * fact
  end function r_dot_r

  !
  !  (1/2m) * p_{i} \cdot p_{j} / hw
  !
  function p_dot_p(a, b, c, d, J) result(r)
    integer, intent(in) :: a(4), b(4), c(4), d(4), J
    real(8) :: r
    integer :: ea, eb, ec, ed
    integer :: ja, jb, jc, jd
    real(8) :: fact

    r = 0.d0

    ea = 2*a(1) + a(2)
    eb = 2*b(1) + b(2)
    ec = 2*c(1) + c(2)
    ed = 2*d(1) + d(2)
    if(abs(ea+eb-ec-ed) > 2 .or. mod(abs(ea+eb-ec-ed),2)==1) return

    ja = a(3)
    jb = b(3)
    jc = c(3)
    jd = d(3)

    fact = 1.d0
    if(a(1) == b(1) .and. a(2) == b(2) .and. a(3) == b(3) .and. a(4) == b(4)) fact = fact / dsqrt(2.d0)
    if(c(1) == d(1) .and. c(2) == d(2) .and. c(3) == d(3) .and. c(4) == d(4)) fact = fact / dsqrt(2.d0)

    r = - (-1.d0) ** ((jb+jc)/2+J) * &
        & sjs(ja,jb,2*J,jd,jc,2) * red_nab_j(a,c) * red_nab_j(b,d) - &
        & (-1.d0) ** ((jb+jc)/2) * &
        & sjs(ja,jb,2*J,jc,jd,2) * red_nab_j(a,d) * red_nab_j(b,c)

    r = r * fact
  end function p_dot_p

  ! < a || r || b >
  function red_r_j(a, b) result(r)
    real(8) :: r
    integer, intent(in) :: a(4), b(4)
    integer :: na, la, ja, za
    integer :: nb, lb, jb, zb

    r = 0.d0

    na = a(1)
    la = a(2)
    ja = a(3)
    za = a(4)

    nb = b(1)
    lb = b(2)
    jb = b(3)
    zb = b(4)

    if(za /= zb) return

    r = (-1.d0) ** ((3+2*la+jb)/2) * dsqrt(dble(ja + 1) * dble(jb + 1)) * &
        &  sjs(ja, 2, jb, 2 * lb, 1, 2 * la) * red_r_l(na, la, nb, lb)
  end function red_r_j

  ! < a || \nabra || b >
  function red_nab_j(a, b) result(r)
    real(8) :: r
    integer, intent(in) :: a(4), b(4)
    integer :: na, la, ja, za
    integer :: nb, lb, jb, zb

    r = 0.d0

    na = a(1)
    la = a(2)
    ja = a(3)
    za = a(4)

    nb = b(1)
    lb = b(2)
    jb = b(3)
    zb = b(4)

    if(za /= zb) return

    r = (-1.d0) ** ((3+2*la+jb)/2) * dsqrt(dble(ja + 1) * dble(jb + 1)) * &
        &  sjs(ja, 2, jb, 2 * lb, 1, 2 * la) * red_nab_l(na, la, nb, lb)
  end function red_nab_j

    real(8) function red_r_l(n1, l1, n2, l2) result(rl)
    integer, intent(in) :: n1, l1, n2, l2
    if (n1 == n2 .and. l1 == l2-1) then
      rl = -dsqrt(dble(l2)*(dble(n2 + l2) + 0.5d0))
    elseif (n1 == n2-1 .and. l1 == l2+1) then
      rl = -dsqrt(dble(l2 + 1)*dble(n2))
    elseif (n1 == n2+1 .and. l1 == l2-1) then
      rl = dsqrt(dble(l2)*(dble(n2  + 1)))
    elseif (n1 == n2 .and. l1==l2+1) then
      rl = dsqrt(dble(l2+1)*(dble(n2 +l2)+1.5d0))
    else
      rl = 0.d0
    end if
  end function red_r_l

  real(8) function red_nab_l(n1, l1, n2, l2) result(nl)
    integer, intent(in) :: n1, l1, n2, l2
    if(n1 == n2 .and. l1 == l2+1) then
      nl = -dsqrt(dble(l2 + 1)*(dble(n2 + l2) + 1.5d0))
    elseif(n1 == n2-1 .and. l1 == l2+1) then
      nl = -dsqrt(dble(l2 + 1)*dble(n2))
    elseif(n1 == n2 .and. l1 == l2-1) then
      nl = -dsqrt(dble(l2)*(dble(n2 + l2) + 0.5d0))
    elseif(n1 == n2+1 .and. l1==l2-1) then
      nl = -dsqrt(dble(l2)*dble(n2 + 1))
    else
      nl = 0.d0
    end if
  end function red_nab_l

  function matel_r4(a, b, c, d, J, isospin, hw) result(res)
    integer, intent(in) :: a(4), b(4), c(4), d(4), J, isospin
    real(8), intent(in) :: hw
    real(8) :: res
    integer :: na, la, ja, za
    integer :: nb, lb, jb, zb
    integer :: nc, lc, jc, zc
    integer :: nd, ld, jd, zd
    real(8) :: norm
    na = a(1); la = a(2); ja = a(3); za = a(4);
    nb = b(1); lb = b(2); jb = b(3); zb = b(4);
    nc = c(1); lc = c(2); jc = c(3); zc = c(4);
    nd = d(1); ld = d(2); jd = d(3); zd = d(4);
    res = 0.d0
    if(triag(ja, jb, 2*J)) return
    if(triag(jc, jd, 2*J)) return
    norm = 1.d0
    if(na==nb .and. la==lb .and. ja==jb .and. za==zb) norm = norm / dsqrt(2.d0)
    if(nc==nd .and. lc==ld .and. jc==jd .and. zc==zd) norm = norm / dsqrt(2.d0)
    res = matel_not_asym(a, b, c, d, J) - matel_not_asym(a, b, d, c, J) * (-1.d0)**((jc+jd)/2-J)
  contains
    function matel_not_asym(i, j, k, l, J2b) result(me)
    integer, intent(in) :: i(4), j(4), k(4), l(4), J2b
    integer :: ni, li, ji, zi
    integer :: nj, lj, jj, zj
    integer :: nk, lk, jk, zk
    integer :: nl, ll, jl, zl
    real(8) :: me, iso1, iso2, bpar
    ni = i(1); li = i(2); ji = i(3); zi = i(4);
    nj = j(1); lj = j(2); jj = j(3); zj = j(4);
    nk = k(1); lk = k(2); jk = k(3); zk = k(4);
    nl = l(1); ll = l(2); jl = l(3); zl = l(4);
    bpar = sqrt(2.d0 * hc**2 / ((amp + amn) * hw))
    me = 0.d0
    if(zi /= zk) return
    if(zj /= zl) return
    if(mod(li + lk + 1, 2) == 1) return
    if(mod(lj + ll + 1, 2) == 1) return
    if(triag(ji, jk, 2)) return
    if(triag(jj, jl, 2)) return
    iso1 = 1.d0; iso2 = 1.d0
    if(isospin==0) then
      iso1 = 1.d0; iso2 = 1.d0
    else if(isospin==1) then
      iso1 = dble(zi); iso2 = dble(zj)
    else
      write(*,*) "Error: ", __FILE__, __LINE__
    end if
    me = (-1.d0)**((jj + jl)/2 + J2b) * sjs(ji, jj, 2*J2b, jl, jk, 2) * &
        & sqrt(dble(ji+1)*dble(jj+1)*dble(jk+1)*dble(jl+1)) * &
        & tjs(ji, 2, jk, -1, 0, 1) * tjs(jj, 2, jl, -1, 0, 1) * &
        & (iso1 * radius_power(3, ni, li, nk, lk, bpar) * radius_power(1, nj, lj, nl, ll, bpar) + &
        & iso2 * radius_power(1, ni, li, nk, lk, bpar) * radius_power(3, nj, lj, nl, ll, bpar))
    end function matel_not_asym
  end function matel_r4

  function radius_power(k, n1, l1, n2, l2, bpar) result(s)
    !
    !  radial integral for the harmonic oscillator wave function
    !
    !  radius_power(k, n1, l1, n2, l2) =  <n1, l1|r^k|n2, l2>
    !
    !  n1, n2: the number of nodes (n1, n2 = 0, 1, 2, ...)
    !  integral of r only
    !
    integer, intent(in) :: k, n1, l1, n2, l2
    real(8), intent(in) :: bpar
    real(8) :: s
    integer :: ll1, ll2, ll, i, imin, imax

    ll = l1 + l2 + k
    ll1 = l2 - l1 + k
    ll2 = l1 - l2 + k
    s = 0.0d0
    if (mod(ll, 2) == 1 .or. ll1 < 0 .or. ll2 < 0) then
       ! direct integral should be done instead
       write(*,*) 'error [radius_power]: input is out of range'
       write(*,'(1a,1i3,2x,1a,1i3,2x,1a,1i3)') &
            'l1 =', l1, 'l2 =', l2, 'k =', k
       stop
    end if


    ll1 = ll1/2
    ll2 = ll2/2

    imin = max(0, n1-ll1, n2-ll2)
    imax = min(n1, n2)
    do i = imin, imax
       s = s + dbinomial(ll1, n1-i) * dbinomial(ll2, n2-i) &
            & * (double_factorial(ll+2*i+1)/double_factorial(2*i))
    end do

    s = s * sqrt(double_factorial(2*n1)/double_factorial(2*n1+2*l1+1)) &
         & * sqrt(double_factorial(2*n2)/double_factorial(2*n2+2*l2+1)) &
         & * sqrt(1.0d0/2**k) * (-1)**(n1-n2) * bpar**k

  end function radius_power

  function dbinomial(n, m) result(s)
    !
    !  binomial coefficient: n_C_m
    !  s: double precision
    !
    integer, intent(in) :: n, m
    real(8) :: s, s1, s2
    integer :: i, m1

    s = 1.0d0
    m1 = min(m, n-m)
    if (m1 == 0) return
    if (n > 1000) then
       write(*,'(1a, 1i6, 1a)') '[dbinomial]: n =', n, ' is too large'
       stop
    end if

    if (n < 250) then
       s1 = 1.0d0
       s2 = 1.0d0
       do i = 1, m1
          s1 = s1 * (n-i+1)
          s2 = s2 * (m1-i+1)
       end do
       s = s1 / s2
    else
       do i = 1, m1
          s = (s * (n-i+1)) / (m1-i+1)
       end do
    endif

  end function dbinomial

end module DefineOperators
