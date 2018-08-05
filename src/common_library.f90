module common_library
  use class_sys, only: sy
  implicit none
  type(sy), private :: sys
  integer, private, parameter :: n_trinomial = 100
  real(8), private, allocatable  :: dtrinomial(:,:,:)
  real(8), parameter :: pi = 3.141592741012573d0 ! \pi
  real(8), parameter:: hc = 197.32705d0         ! \hbar c [MeV fm]
  real(8), parameter :: alpha = 137.035999d0     ! electric fine structure constant
  real(8) :: amp      ! proton mass [MeV] can be changed in LQCD calc.
  real(8) :: amn      ! neutron mass [MeV] can be changed in LQCD calc.
  real(8) :: rmass    ! reduced mass
  real(8) :: amnucl   ! averaged nucleon mass

contains
  subroutine set_physics_constant
    amp = 938.27231d0        ! proton mass [MeV] can be changed in LQCD calc.
    amn = 939.56563d0        ! neutron mass [MeV] can be changed in LQCD calc.
    rmass = (amp * amn) / (amp + amn)
    amnucl = (amp + amn) * 0.5d0
  end subroutine set_physics_constant

  subroutine skip_comment(nfile, comment)
    implicit none
    integer,intent(in)::nfile
    character(*), intent(in) :: comment
    character(20) :: line
    read(nfile,'(a)') line
    do while  (sys%find(line, comment))
      read(nfile,'(a)') line
    end do
    backspace(nfile)
  end subroutine skip_comment

  real(8) function Del(i1, i2)
    integer, intent(in) :: i1, i2
    Del = 1.d0
    if(i1 == i2) Del = dsqrt(2.d0)
  end function Del

  real(8) function hat(j)
    integer, intent(in) :: j
    hat = sqrt(dble(j+1))
  end function hat

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

  logical function triag(i,j,k)
    implicit none
    integer,intent(in)::i,j,k
    triag = ((i-(j+k))*(i-abs(j-k)) > 0)
  end function triag

  subroutine init_dbinomial_triangle()
    ! cache for tuning 6-j coefficient
    integer :: n, m, i, j, k, info
    real(8) :: d

    allocate(dtrinomial(0:n_trinomial, 0:n_trinomial, 0:n_trinomial))
    dtrinomial(:,:,:) = 0.d0
    dtrinomial(0,0,0) = 1.d0
    dtrinomial(1,1,1) = 1.d0
    do i = 2, n_trinomial
      m = i - i/2 * 2
      dtrinomial(i,i,m) = 1.d0
      dtrinomial(i,m,i) = 1.d0
      n = m + 2
      do j = i, n, -2
        do k = n, j, 2
          dtrinomial(i, j, k) = dtrinomial(i, j, k-2) / dble(k)
          dtrinomial(i, k, j) = dtrinomial(i, j, k)
        end do
        dtrinomial(i, j-2, m) = dtrinomial(i, j, m) * dble(j)
        dtrinomial(i, m, j-2) = dtrinomial(i, j-2, m)
      end do
    end do

  end subroutine init_dbinomial_triangle

  subroutine fin_dbinomial_triangle()
    deallocate(dtrinomial)
  end subroutine fin_dbinomial_triangle

  !
  ! factorials for 3j,6j and 9j symbols
  ! for moshinsky trans brackets and for
  ! vector brackets
  !

  ! ------------ rotation group  ----------------------------

  function dcg(j1, m1, j2, m2, j3, m3) result(s)
    !
    !  Clebsch-Gordan coefficient
    !
    !  dcg(j1, m1, j2, m2, j3, m3)
    !  = ((j1)/2, (m1)/2, (j2)/2, (m2)/2 | (j3)/2, (m3)/2)
    implicit none
    integer, intent(in) :: j1, j2, j3, m1, m2, m3
    real(8) :: cg, s

    s = cg(j1, j2, j3, m1, m2, m3)
  end function dcg

  function sjs(j1, j2, j3, l1, l2, l3) result(s)
    !
    !  6j coefficient
    !
    !  d6j(j1, j2, j3, l1, l2, l3) = {(j1)/2 (j2)/2 (j3)/2}
    !                                {(l1)/2 (l2)/3 (l3)/2}
    implicit none
    integer, intent(in) :: j1, j2, j3, l1, l2, l3
    real(8) :: s, sixj

    s = sixj(j1, j2, j3, l1, l2, l3)

  end function sjs

  function snj(j11, j12, j13, j21, j22, j23, j31, j32, j33) result(s)
    !
    !  9j coefficient
    !
    !  d9j(j11, j12, j13, j21, j22, j23, j31, j32, j33)
    !
    !    {(j11)/2 (j12)/2 (j13)/2}
    !  = {(j21)/2 (j22)/2 (j23)/2}
    !    {(j31)/2 (j32)/2 (j33)/2}
    !
    implicit none
    integer, intent(in) :: j11, j12, j13, j21, j22, j23, j31, j32, j33
    real(8) :: s, ninej

    s = ninej(j11, j12, j13, j21, j22, j23, j31, j32, j33)

  end function snj

  function gmosh(nl, ll, nr, lr, n1, l1, n2, l2, lm, d) result(r)
    real(8) :: r
    integer, intent(in) :: nl, ll, nr, lr, n1, l1, n2, l2, lm
    real(8), intent(in) :: d
    integer :: ee, er, e1, e2, m, ed, eb, ec, ea, ld, lb, lc, la
    real(8) :: t, s

    r = 0.d0
    ee = 2*nl + ll
    er = 2*nr + lr
    e1 = 2*n1 + l1
    e2 = 2*n2 + l2
    if(ee + er /= e1 + e2) return
    if(triag(ll, lr, lm)) return
    if(triag(l1, l2, lm)) return
    t = dsqrt((d ** (e1 - er)) / ((1.d0 + d) ** (e1 + e2)))
    m = min(er, e2)
    s = 1.d0
    do ed = 0, m
      eb = er - ed
      ec = e2 - ed
      ea = e1 - er + ed

      do ld = ed, 0, -2
        do lb = eb, 0, -2
          if(triag(ld,lb,lr)) cycle
          do lc = ec, 0, -2
            if(triag(ld,lc,l2)) cycle
            do la = ea, 0, -2
              if(triag(la,lb,l1)) cycle
              if(triag(la,ll,lc)) cycle

              r = r + s * t * &
                  & snj(2*la, 2*lb, 2*l1, 2*lc, 2*ld, 2*l2, 2*ll, 2*lr, 2*lm) * &
                  & g(e1, l1, ea, la, eb, lb) * g(e2, l2, ec, lc, ed, ld) * &
                  & g(ee, ll, ea, la, ec, lc) * g(er, lr, eb, lb, ed, ld)

            end do
          end do
        end do
      end do
      s = s * (-d)
    end do
    r = r * (-1.d0) ** (n1 + n2 + nr + nl)
  contains
    function g(e1, l1, ea, la, eb, lb) result(r)
      real(8) :: r
      integer, intent(in) :: e1, l1, ea, la, eb, lb

      r = dcg(2*la, 0, 2*lb, 0, 2*l1, 0) * dsqrt((2*la + 1) * (2*lb + 1) * &
          & dtrinomial(e1 - l1, ea - la, eb - lb) * &
          & dtrinomial(e1 + l1 + 1, ea + la + 1, eb + lb + 1))

    end function g
  end function gmosh

end module common_library
