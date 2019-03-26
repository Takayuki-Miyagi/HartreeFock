module MyLibrary
  use, intrinsic :: iso_c_binding
  use ClassSys, only: sys
  implicit none

  real(8), parameter, public :: pi = 3.141592741012573d0 ! \pi
  real(8), parameter, public :: hc = 197.32705d0         ! \hbar c [MeV fm]
  real(8), public :: amp = 938.27231d0        ! proton mass [MeV] can be changed in LQCD calc.
  real(8), public :: amn = 939.56563d0        ! neutron mass [MeV] can be changed in LQCD calc.
  real(8), parameter, public :: alpha = 137.035999d0     ! electric fine structure constant

  ! cache for Talmi-Moshinsky bracket
  integer, private, parameter :: n_trinomial = 100
  real(8), private, allocatable  :: dtrinomial(:,:,:)

  ! C interfaces
  interface
    ! 3-j symbol
    function coupling_3j(j1,j2,j3,m1,m2,m3) bind(c,name='gsl_sf_coupling_3j')
      import c_int, c_double
      real(c_double) :: coupling_3j
      integer(c_int), value, intent(in) :: j1,j2,j3,m1,m2,m3
    end function coupling_3j

    ! 6-j symbol
    function coupling_6j(j1,j2,j3,j4,j5,j6) bind(c,name='gsl_sf_coupling_6j')
      import c_int, c_double
      real(c_double) :: coupling_6j
      integer(c_int), value, intent(in) :: j1,j2,j3,j4,j5,j6
    end function coupling_6j

    ! 9-j symbol
    function coupling_9j(j1,j2,j3,j4,j5,j6,j7,j8,j9) bind(c,name='gsl_sf_coupling_9j')
      import c_int, c_double
      real(c_double) :: coupling_9j
      integer(c_int), value, intent(in) :: j1,j2,j3,j4,j5,j6,j7,j8,j9
    end function coupling_9j

    ! factorial n!
    function factorial(n) bind(c,name='gsl_sf_fact')
      import c_double, c_int
      integer(c_int), value, intent(in) :: n
      real(c_double) :: factorial
    end function factorial

    ! double factorial n!!
    function double_factorial(n) bind(c,name='gsl_sf_doublefact')
      import c_double, c_int
      integer(c_int), value, intent(in) :: n
      real(c_double) :: double_factorial
    end function double_factorial

    ! Gamma function \Gamma(x)
    function gamma_function(x) bind(c,name='gsl_sf_gamma')
      import c_double
      real(c_double), value, intent(in) :: x
      real(c_double) :: gamma_function
    end function gamma_function

    ! log of Gamma function ln \Gamma(x)
    function ln_gamma(x) bind(c,name='gsl_sf_lngamma')
      import c_double
      real(c_double), value, intent(in) :: x
      real(c_double) :: ln_gamma
    end function ln_gamma

    ! spherical bessel function j_l(x)
    function spherical_bessel(l,x) bind(c,name='gsl_sf_bessel_jl')
      import c_int, c_double
      integer(c_int), value, intent(in) :: l
      real(c_double), value, intent(in) :: x
      real(c_double) :: spherical_bessel
    end function spherical_bessel

    ! Legendre polynomial P_l(x)
    function legendre_polynomial(l,x) bind(c,name='gsl_sf_legendre_Pl')
      import c_int, c_double
      integer(c_int), value, intent(in) :: l
      real(c_double), value, intent(in) :: x
      real(c_double) :: legendre_polynomial
    end function legendre_polynomial

    ! associated Laguerre polynomial L^{(a)}_{n}(x)
    function laguerre(n,a,x) bind(c,name='gsl_sf_laguerre_n')
      import c_int, c_double
      integer(c_int), value, intent(in) :: n
      real(c_double), value, intent(in) :: a, x
      real(c_double) :: laguerre
    end function laguerre

    ! Gauss-Legendre quadrature
    function gauss_legendre_allocate(n) &
          & bind(c,name='gsl_integration_glfixed_table_alloc')
      import c_int, c_ptr
      integer(c_int), value, intent(in) :: n
      type(c_ptr) :: gauss_legendre_allocate
    end function gauss_legendre_allocate
    function gauss_legendre_ith_point_weight(a,b,i,xi,wi,t) &
          & bind(c,name='gsl_integration_glfixed_point')
      import c_int, c_double, c_ptr
      real(c_double), value, intent(in) :: a, b
      integer(c_int), value, intent(in) :: i
      real(c_double) :: xi, wi
      type(c_ptr), value, intent(in) :: t
      integer(c_int) :: gauss_legendre_ith_point_weight
    end function gauss_legendre_ith_point_weight
    subroutine gauss_legendre_release(t) &
          & bind(c,name='gsl_integration_glfixed_table_free')
      import c_ptr
      type(c_ptr), value :: t
    end subroutine gauss_legendre_release

    ! open, read, write, and close gzip file (additional interface gzip_open below)
    function gz_open(filename, mode) bind(c, name='gzopen')
      import c_char, c_ptr
      character(c_char) :: filename(*), mode(*)
      type(c_ptr) :: gz_open
    end function gz_open
    function gzip_read(f, buf, len) bind(c, name='gzread')
      import c_int, c_char, c_ptr
      type(c_ptr), value :: f
      character(c_char) :: buf(*)
      integer(c_int), value, intent(in) :: len
      type(c_ptr) :: gzip_read
    end function gzip_read
    function gzip_readline(f, buf, len) bind(c, name='gzgets')
      import c_int, c_char, c_ptr
      type(c_ptr), value :: f
      character(c_char) :: buf(*)
      integer(c_int), value, intent(in) :: len
      type(c_ptr) :: gzip_readline
    end function gzip_readline
    function gzip_write(f, buf, len) bind(c, name='gzwrite')
      import c_int, c_char, c_ptr
      type(c_ptr), value :: f
      character(c_char) :: buf(*)
      integer(c_int), value, intent(in) :: len
      integer(c_int) :: gzip_write
    end function gzip_write
    function gzip_close( f ) bind(c, name='gzclose')
      import c_ptr
      type(c_ptr), value :: f
      type(c_ptr) :: gzip_close
    end function gzip_close
  end interface
contains

  subroutine skip_comment(nfile, comment)
    implicit none
    integer,intent(in)::nfile
    character(*), intent(in) :: comment
    type(sys) :: s
    character(20) :: line
    read(nfile,'(a)') line
    do while  (s%find(line, comment))
      read(nfile,'(a)') line
    end do
    backspace(nfile)
  end subroutine skip_comment

  real(8) function delta(i1, i2)
    integer, intent(in) :: i1, i2
    delta = 0.d0
    if(i1 == i2) delta = 1.d0
  end function delta

  function hat(j) result(r)
    real(8) :: r
    integer, intent(in) :: j
    r = dsqrt(dble(j + 1))
  end function hat

  ! 3-j symbol for Wigner-Eckert theorem
  real(8) function geometry_part(jbra, jop, jket, mbra, mop, mket) result(r)
    integer, intent(in) :: jbra, mbra, jop, mop, jket, mket
    r = (-1.d0) ** ((jbra - mbra)/2) * &
      & tjs(jbra, jop, jket, -mbra, mop, mket)
  end function geometry_part

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

  logical function triag(i,j,k)
      implicit none
      integer,intent(in)::i,j,k
      triag = ((i-(j+k))*(i-abs(j-k)) > 0)
  end function triag

  subroutine init_dtrinomial()
    ! cache for Talmi-Moshinksy bracket
    integer :: i, j, k
    allocate(dtrinomial(0:n_trinomial, 0:n_trinomial, 0:n_trinomial))
    !$omp parallel do private( i, j, k ) schedule (dynamic)
    do k = 0, n_trinomial
      do j = 0, n_trinomial
        do i = 0, n_trinomial
          dtrinomial(i, j, k) = dtrinomial_func(i, j, k)
        end do
      end do
    end do
  end subroutine init_dtrinomial

  subroutine fin_dtrinomial()
    deallocate(dtrinomial)
  end subroutine fin_dtrinomial

  function dcg(j1, m1, j2, m2, j3, m3) result(s)
    !
    !  Clebsch-Gordan coefficient
    !  dcg(j1, m1, j2, m2, j3, m3)
    !  = ((j1)/2, (m1)/2, (j2)/2, (m2)/2 | (j3)/2, (m3)/2)
    integer, intent(in) :: j1, j2, j3, m1, m2, m3
    real(8) :: s
    s = coupling_3j(j1,j2,j3,m1,m2,-m3) * hat(j3) * (-1.d0) ** ((j1-j2+m3)/2)
  end function dcg

  function tjs(j1, j2, j3, m1, m2, m3) result(r)
    real(8) :: r
    integer, intent(in) :: j1, j2, j3, m1, m2, m3
    r = coupling_3j(j1,j2,j3,m1,m2,m3)
  end function tjs

  function sjs(j1, j2, j3, l1, l2, l3) result(s)
    !
    !  6j coefficient
    !  d6j(j1, j2, j3, l1, l2, l3) = {(j1)/2 (j2)/2 (j3)/2}
    !                                {(l1)/2 (l2)/3 (l3)/2}
    integer, intent(in) :: j1, j2, j3, l1, l2, l3
    real(8) :: s
    s = coupling_6j(j1,j2,j3,l1,l2,l3)
  end function sjs

  function snj(j11, j12, j13, j21, j22, j23, j31, j32, j33) result(s)
    !
    !  9j coefficient
    !  d9j(j11, j12, j13, j21, j22, j23, j31, j32, j33)
    !    {(j11)/2 (j12)/2 (j13)/2}
    !  = {(j21)/2 (j22)/2 (j23)/2}
    !    {(j31)/2 (j32)/2 (j33)/2}
    integer, intent(in) :: j11, j12, j13, j21, j22, j23, j31, j32, j33
    real(8) :: s
    s = coupling_9j(j11,j12,j13,j21,j22,j23,j31,j32,j33)
  end function snj

  function gmosh(nl, ll, nr, lr, n1, l1, n2, l2, lm, d) result(r)
    !
    ! reference G.P.Kamuntavicius, R.K.Kalinauskas, B.R.Barrett, S.Mickevicius, and D. Germanas, Nucl. Phys. A 695, 191 (2001).
    !
    real(8) :: r
    integer, intent(in) :: nl, ll, nr, lr, n1, l1, n2, l2, lm
    real(8), intent(in) :: d
    integer :: ee, er, e1, e2, m, ed, eb, ec, ea, ld, lb, lc, la
    real(8) :: t, s

    r = 0.d0
    if(.not. allocated(dtrinomial) ) then
      write(*,'(a)') "you need to call init_dtrinomial first!"
      return
    end if
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

  function dtrinomial_func(i, j, k) result(s)
    !
    !  trinomial coefficient: i!! / j!! / k!!
    !
    integer, intent(in) :: i, j, k
    real(8) :: s
    integer :: m
    s = 1.d0
    m = max(i, j, k)
    if(m == 0) return
    if(m > n_trinomial) then

      write(*,'(a)') 'in trinomial_func, index is too large'
      return

    end if
    s = double_factorial(i) / (double_factorial(j) * double_factorial(k))
  end function dtrinomial_func

  function ho_radial_wf_norm(n,l,anu,r) result(s)
    ! R = sqrt( 2 * nu * Gamma(n+1) / Gamma(n+l+1.5) ) x^{l+1} e^{ -x^2/2 } L^{n}_{l+0.5}( x^2 )
    ! x = r * nu or p * nu
    ! nu = sqrt(m w / h) or sqrt( h / mw )
    ! Note that anu = nu^2
    integer,intent(in) :: n,l
    real(8),intent(in) :: anu,r
    real(8) :: s
    real(8) :: prefact, exp_component, nu, x, a
    nu = sqrt(anu)
    x = nu * r
    a = dble(l)+0.5d0
    prefact = sqrt( 2.d0 * nu )
    exp_component = 0.5d0*ln_gamma(dble(n+1)) - 0.5d0*ln_gamma(dble(n+l)+1.5d0) - 0.5d0*x**2 &
        & + dble(l+1) * log(x)
    s = prefact * exp(exp_component) * laguerre(n,a,x**2)
  end function ho_radial_wf_norm

  subroutine gauss_legendre(x1,x2,x,w,n)
    ! input:
    ! x1   : lower limit of the integration interval
    ! x2   : upper limit ---------- "" -------------
    ! n    : the desired number of mesh points
    ! output :
    ! x     : gauss-legendre mesh points on the interval (x1,x2)
    ! w     : the corresponding weights
    integer, intent(in) :: n
    real(8), intent(in) :: x1, x2
    real(8), intent(inout) :: x(:), w(:)
    real(8) :: xi, wi
    integer :: info, i
    type(c_ptr) :: t
    t = gauss_legendre_allocate(n)
    do i = 1, n
      info = gauss_legendre_ith_point_weight(x1,x2,i-1,xi,wi,t)
      x(i) = xi
      w(i) = wi
    end do
    call gauss_legendre_release(t)
  end subroutine gauss_legendre

  function gzip_open( filename, mode ) result(p)
    character(*), intent(in) :: filename, mode
    type(c_ptr) :: p
    p = gz_open(trim(filename)//achar(0), trim(mode)//achar(0))
  end function gzip_open
end module MyLibrary

