module common_library
  use class_sys, only: sy
  implicit none
  type(sy), private :: sys
  integer, private, parameter :: n_dbinomial = 200, n_triag = 100
  real(8), private, allocatable  :: dbinomial(:,:), triangle_c(:,:,:)
  real(8), parameter :: pi = 3.141592741012573d0 ! \pi
  real(8), parameter:: hc = 197.32705d0         ! \hbar c [MeV fm]
  real(8), parameter :: alpha = 137.035999d0     ! electric fine structure constant
  real(8) :: amp      ! proton mass [MeV] can be changed in LQCD calc.
  real(8) :: amn      ! neutron mass [MeV] can be changed in LQCD calc.
  real(8) :: rmass    ! reduced mass
  real(8) :: amnucl   ! averaged nucleon mass

contains
  subroutine set_physics_constant()
    amp = 938.27231d0        ! proton mass [MeV] can be changed in LQCD calc.
    amn = 939.56563d0        ! neutron mass [MeV] can be changed in LQCD calc.
    rmass = (amp * amn) / (amp + amn)
    amnucl = (amp + amn) * 0.5d0
  end subroutine set_physics_constant

  real(8) function hat(i)
    integer, intent(in) :: i
    hat = dsqrt(dble(i + 1))
  end function hat

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

    allocate( dbinomial(0:n_dbinomial, 0:n_dbinomial) )

    dbinomial(:,:) = 0.d0

    !$omp parallel do private(n, m) schedule(dynamic)
    do n = 0, n_dbinomial
      do m = 0, n
        dbinomial(n, m) = dbinomial_func(n, m)
      end do
    end do

    allocate(triangle_c(0:n_triag, 0:n_triag, 0:n_triag))
    triangle_c(:,:,:) = 0.d0
    !$omp parallel do private( i, j, k, d, info ) schedule (dynamic)
    do k = 0, n_triag
      do j = 0, n_triag
        do i = 0, n_triag
          call triangle(i, j, k, d, info)
          if (info == 1) then
            triangle_c(i, j, k) = 0.d0
            cycle
          end if
          triangle_c(i, j, k) = d
        end do
      end do
    end do

  end subroutine init_dbinomial_triangle

  subroutine fin_dbinomial_triangle()
    deallocate(dbinomial, triangle_c)
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
    !
    !  using the formula by Racah (1942)
    !
    implicit none
    integer, intent(in) :: j1, j2, j3, m1, m2, m3
    integer :: js, jm1, jm2, jm3, k1, k2, k3
    integer :: iz, izmin, izmax, isign
    double precision :: s
    double precision :: tmp, delta

    s = 0.0d0
    jm1 = j1 - m1
    jm2 = j2 - m2
    jm3 = j3 - m3

    ! error and trivial-value check
    if (abs(m1) > j1 .or. abs(m2) > j2 .or. abs(m3) > j3 .or. &
        & mod(jm1, 2) /= 0 .or. mod(jm2, 2) /= 0 .or. mod(jm3, 2) /= 0) then
      write(*,*) 'error [dcg]: invalid j or m'
      write(*,'(1a, 1i4, 2x, 1a, 1i4)') 'j1 =', j1, 'm1 =', m1
      write(*,'(1a, 1i4, 2x, 1a, 1i4)') 'j2 =', j2, 'm2 =', m2
      write(*,'(1a, 1i4, 2x, 1a, 1i4)') 'j3 =', j3, 'm3 =', m3
      stop
    end if

    ! call triangle(j1, j2, j3, delta, info)
    ! if (info == 1) return
    if (max(j1, j2, j3) > n_triag) stop 'increase n_triag'
    delta = triangle_c( j1, j2, j3 )
    if (delta == 0.d0) return

    if (m3 /= m1+m2) return

    jm1 = jm1 / 2
    jm2 = jm2 / 2
    jm3 = jm3 / 2
    js = (j1 + j2 + j3)/2
    k1 = (j2 + j3 - j1)/2
    k2 = (j3 + j1 - j2)/2
    k3 = (j1 + j2 - j3)/2

    if (max(j1, j2, j3) > n_dbinomial) stop 'increase n_dbinomial'

    tmp = sqrt(dbinomial(j1, k2)/dbinomial(j1, jm1)) &
        & * sqrt(dbinomial(j2, k3)/dbinomial(j2, jm2)) &
        & * sqrt(dbinomial(j3, k1)/dbinomial(j3, jm3)) &
        & * sqrt((j3+1.0d0)) * delta

    izmin = max(0, jm1-k2, k3-jm2)
    izmax = min(k3, jm1, j2-jm2)

    if (izmax > n_dbinomial) stop 'increase n_dbinomial'

    isign = (-1)**izmin
    do iz = izmin, izmax
      s = s + isign * dbinomial(k3,iz) * dbinomial(k2,jm1-iz) &
          & * dbinomial(k1,j2-jm2-iz)
      isign = isign * (-1)
    end do

    s = s * tmp

  end function dcg

  function sjs(j1, j2, j3, l1, l2, l3) result(s)
    !
    !  6j coefficient
    !
    !  d6j(j1, j2, j3, l1, l2, l3) = {(j1)/2 (j2)/2 (j3)/2}
    !                                {(l1)/2 (l2)/3 (l3)/2}
    !
    !  see I. Talmi, Simple Models of Complex Nuclei, p. 158
    !
    implicit none
    integer, intent(in) :: j1, j2, j3, l1, l2, l3
    double precision :: s
    double precision :: d
    integer :: izmin, izmax, iz, isign
    integer :: js, k1, k2, k3, jl1, jl2, jl3

    s = 0.0d0

    if (max(j1, j2, j3, l1, l2, l3) > n_triag) stop 'increase n_triag'

    d = triangle_c(j1, j2, j3)
    if (d == 0.d0) return

    d =    triangle_c(j1, l2, l3) * triangle_c(l1, j2, l3) &
        * triangle_c(l1, l2, j3) / d
    if ( d == 0.d0 ) return

    ! call triangle(j1, j2, j3, deltas, infos)
    ! call triangle(j1, l2, l3, delta1, info1)
    ! call triangle(l1, j2, l3, delta2, info2)
    ! call triangle(l1, l2, j3, delta3, info3)
    ! if (infos == 1 .or. info1 == 1 .or. info2 == 1 .or. info3 == 1) return

    js = (j1 + j2 + j3)/2
    k1 = (j2 + j3 - j1)/2
    k2 = (j3 + j1 - j2)/2
    k3 = (j1 + j2 - j3)/2
    jl1 = (j1 + l2 + l3)/2
    jl2 = (l1 + j2 + l3)/2
    jl3 = (l1 + l2 + j3)/2

    izmin = max(0, js, jl1, jl2, jl3)
    izmax = min(k1+jl1, k2+jl2, k3+jl3)

    if (izmax+1 > n_dbinomial) stop 'increase n_dbinomial'

    isign = (-1)**izmin
    do iz = izmin, izmax
      s = s + isign * dbinomial(iz+1, iz-js) &
          & * dbinomial(k1, iz-jl1) * dbinomial(k2, iz-jl2) &
          & * dbinomial(k3, iz-jl3)
      isign = isign * (-1)
    end do
    ! s = s * delta1 * delta2 * delta3 / deltas
    s = s * d

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
    !  see I. Talmi, Simple Models of Complex Nuclei, p. 968
    !
    implicit none
    integer, intent(in) :: j11, j12, j13, j21, j22, j23, j31, j32, j33
    double precision :: s
    integer :: k, kmin, kmax

    kmin = max(abs(j11-j33), abs(j12-j23), abs(j21-j32))
    kmax = min(j11+j33, j12+j23, j21+j32)

    s = 0.0d0
    do k = kmin, kmax, 2
      s = s + (k+1.0d0) &
          & * sjs(j11, j12, j13, j23, j33, k) &
          & * sjs(j21, j22, j23, j12, k, j32) &
          & * sjs(j31, j32, j33, k, j11, j21)
    end do
    s = s * (-1)**kmin

  end function snj

  function dbinomial_func(n, m) result(s)
    !
    !  binomial coefficient: n_C_m
    !  s: double precision
    !
    integer, intent(in) :: n, m
    double precision :: s, s1, s2
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

  end function dbinomial_func

  subroutine triangle(j1, j2, j3, delta, info)
    !
    !  triangle often used in calculation of 3j, 6j etc.
    !  delta
    !  = sqrt(((j1+j2-j3)/2)!((j1-j2+j3)/2)!((-j1+j2+j3)/2)!/(1+(j1+j2+j3)/2)!)
    !
    implicit none
    integer, intent(in) :: j1, j2, j3
    double precision, intent(out) :: delta
    integer, intent(out) :: info
    integer :: js, k1, k2, k3

    info = 0
    js = j1 + j2 + j3
    k1 = j2 + j3 - j1
    k2 = j3 + j1 - j2
    k3 = j1 + j2 - j3

    if (j1 < 0 .or. j2 < 0 .or. j3 < 0) then
      write(*,*) 'error [triangle]: invalid j'
      write(*,'(1a, 1i4, 2x, 1a, 1i4, 2x, 1a, 1i4)') &
          & 'j1 =', j1, 'j2 =', j2, 'j3 =', j3
      stop
    end if
    if (k1 < 0 .or. k2 < 0 .or. k3 <0 .or. mod(js, 2) /=0) then
      info = 1
      return
    endif

    ! exclude too large arguments to prevent from round-off error
    if (js > 300) then
      write(*,'(1a, 1i5, 1a)') '[triangle]: j1+j2+j3 =', js, ' is too large'
      stop
    end if

    js = js / 2
    k1 = k1 / 2

    if (js > n_dbinomial) stop 'increase n_dbinomial'

    delta = 1.0d0 / &
        & (sqrt(dbinomial(js, j3)) * sqrt(dbinomial(j3, k1)) * sqrt(js+1.0d0))

  end subroutine triangle

end module common_library
