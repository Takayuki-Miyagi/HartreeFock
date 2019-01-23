module DefineOperators
  implicit none
contains
  subroutine GetOperatorRank(optr, jr, pr, zr)
    character(*), intent(in) :: optr
    integer, intent(out) :: jr, pr, zr

    jr = 0
    pr = 0
    zr = 0

    select case(optr)
    case('hamil', 'Hamil', "HOHamil", "Hohamil", 'hohamil' ,&
          & 'CMHamil', 'CMhamil', 'cmhamil', 'RM2', 'Rm2', "rm2", &
          & 'Tcm', 'tcm')
      jr = 0
      pr = 1
      zr = 0
    case default
      write(*,'(2a)') 'Unknown operator: ', trim(optr)
    end select
  end subroutine GetOperatorRank

  function one_body_element(optr, ia, ib, hw, A, Z, N) result(r)
    use CommonLibrary, only: hc, amnucl
    real(8) :: r
    character(*), intent(in) :: optr
    real(8), intent(in) :: hw
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

    case("CMHamil","CMhamil")
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
      r = r * (1.d0 - 1.d0/dble(A)) * hc ** 2 / (amnucl * hw * A)
      return

    case default
      write(*,'(3a)') "In one_body_element, ", &
          & trim(optr), " is not implemented"
      return
    end select
  end function one_body_element

  function two_body_element(optr, ia, ib, ic, id, Jab, Jcd, hw, A, Z, N) result(r)
    use CommonLibrary, only: hc, amnucl
    real(8) :: r
    character(*), intent(in) :: optr
    real(8), intent(in) :: hw
    integer, intent(in) :: ia(4), ib(4), ic(4), id(4), A, Z, N, Jab, Jcd

    r = 0.d0
    select case(optr)
    case("Hamil", "hamil")
      return

    case('Tcm', 'tcm')
      if(Jab /= Jcd) then
        write(*,'(a,2i3)') "Error in SetTwoBodyChannel: ", Jab, Jcd
      end if
      r = p_dot_p(ia,ib,ic,id,Jab) * hw / dble(A)
      return

    case("CMHamil","CMhamil")
      if(Jab /= Jcd) then
        write(*,'(a,2i3)') "Error in SetTwoBodyChannel: ", Jab, Jcd
      end if
      r = (p_dot_p(ia,ib,ic,id,Jab)+r_dot_r(ia,ib,ic,id,Jab)) * hw/dble(A)
      return

    case('RM2', 'Rm2', 'rm2')
      if(Jab /= Jcd) then
        write(*,'(a,2i3)') "Error in SetTwoBodyChannel: ", Jab, Jcd
      end if
      r = - 4.d0 * r_dot_r(ia,ib,ic,id,Jab) * hc**2 / (amnucl*hw*dble(A)**2)
      return

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
    use CommonLibrary, only: sjs
    integer, intent(in) :: a(4), b(4), c(4), d(4), J
    real(8) :: r
    integer :: ja, jb, jc, jd
    real(8) :: fact

    r = 0.d0
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
    use CommonLibrary, only: sjs
    integer, intent(in) :: a(4), b(4), c(4), d(4), J
    real(8) :: r
    integer :: ja, jb, jc, jd
    real(8) :: fact

    r = 0.d0
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
    use CommonLibrary, only: red_r_l, sjs
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
    use CommonLibrary, only: red_nab_l, sjs
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

end module DefineOperators
