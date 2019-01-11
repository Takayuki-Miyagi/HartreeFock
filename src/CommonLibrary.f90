module CommonLibrary
  use ClassSys, only: sys
  implicit none

  integer, private, parameter :: n_dbinomial = 1000, n_triag = 500, n_trinomial = 100
  real(8), private, allocatable  :: dbinomial(:,:), triangle_c(:,:,:), dtrinomial(:,:,:)
  real(8), parameter, public :: pi = 3.141592741012573d0 ! \pi
  real(8), parameter, public :: hc = 197.32705d0         ! \hbar c [MeV fm]
  real(8), public :: amp = 938.27231d0        ! proton mass [MeV] can be changed in LQCD calc.
  real(8), public :: amn = 939.56563d0        ! neutron mass [MeV] can be changed in LQCD calc.
  real(8), parameter, public :: alpha = 137.035999d0     ! electric fine structure constant
  real(8) :: rmass
  real(8) :: amnucl

  real(8), private, allocatable :: gamal(:), dsq(:,:)
  real(8) :: f_mb(200), g_mb(200), w_mb(200)

contains
  subroutine set_physics_constant()
    amp = 938.27231d0        ! proton mass [MeV] can be changed in LQCD calc.
    amn = 939.56563d0        ! neutron mass [MeV] can be changed in LQCD calc.
    rmass = (amp * amn) / (amp + amn)
    amnucl = (amp + amn) * 0.5d0
  end subroutine set_physics_constant

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

  real(8) function Del(i1, i2)
    integer, intent(in) :: i1, i2
    Del = dsqrt(1.d0 + delta(i1,i2))
  end function Del

  real(8) function delta(i1, i2)
    integer, intent(in) :: i1, i2
    delta = 0.d0
    if(i1 == i2) delta = 1.d0
  end function delta

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

  !c========================================
  !
  !
  !      this routine calculates gauss-legendre mesh points and weights
  !      input:
  !      x1   : lower limit of the integration interval
  !      x2   : upper limit ---------- "" -------------
  !      n    : the desired number of mesh points
  !      output :
  !      x     : gauss-legendre mesh points on the interval (x1,x2)
  !      w     : the corresponding weights
  !      from  : numerical recipes
  !      f90 version : m. hjorth-jensen
  !
  subroutine gauss_legendre(x1,x2,x,w,n)
    implicit none
    integer, intent(in) :: n
    integer :: i, j, m
    real(8), intent(in) :: x1, x2
    real(8), intent(inout) :: x, w
    real(8) :: eps
    dimension :: x(n), w(n)
    parameter (eps=3.d-14)
    real(8) :: p1,p2,p3,pp,xl,xm,z,z1

    m=(n+1)/2
    xm=0.5d0*(x2+x1)
    xl=0.5d0*(x2-x1)
    pp = 0.d0
    do i=1,m
      z1=0.
      z=cos(pi*(i-0.25d0)/(n+0.5d0))
      do while ( abs(z-z1) > eps)
        p1=1.d0
        p2=0.d0
        do j=1,n
          p3=p2
          p2=p1
          p1=((2.0d0*j-1.d0)*z*p2-(j-1.0d0)*p3)/j
        enddo
        pp=n*(z*p1-p2)/(z*z-1.d0)
        z1=z
        z=z-p1/pp
      enddo
      x(i)=xm-xl*z
      x(n+1-i)=xm+xl*z
      w(i)=2.d0*xl/((1.0d0-z*z)*pp*pp)
      w(n+1-i)=w(i)
    enddo

  end subroutine gauss_legendre

  subroutine waver(n,l,anu,r,wave)
    !     the radial part wave function (wave=rr(r)) int wave(r)**2 dr =1)
    !     anu is the size parameter of nuclear well
    !     anu=mass times omega over h bar
    integer,intent(in) :: n,l
    integer :: nnn
    real(8),intent(in) :: anu,r
    real(8),intent(out) :: wave
    real(8) :: aq,bq,guerpq,sig,zzq,scale, dlanu, zz
    real(8) :: wavel, guerp, a, b

    dlanu=dlog(anu)
    zz=anu*r*r
    wavel=0.25d0*dlanu-zz/2.d0 &
        &     +dble(l+1)*(0.5d0*dlanu+dlog(r))
    if (n==0) then
      guerp=dexp(0.5d0*(dlog(2.d0)-gamal(2*l+3)))
      wave=dexp(wavel)*guerp
    elseif (n==1) then
      guerp=dexp(0.5d0*(dlog(2.d0)-gamal(2*l+5))) &
          &                   *(dble(l)+1.5d0-zz)
      wave=dexp(wavel)*guerp
    else
      if (n>220.or.r>21.d0) then
        zzq=zz
        guerpq=exp(0.5d0*(log(2.0)-gamal(2*l+5))) &
            &           *(dble(l)+1.5-zzq)
        aq=exp(0.5d0*(log(2.0d0)-gamal(2*l+3)))
        scale=1.d0
        do nnn=2,n
          if (abs(aq)>1.0d150.or.abs(guerpq)>1.0d+150) then
            scale=scale*1.0d+150
            aq=aq*1.0d-150
            guerpq=guerpq*1.0d-150
          endif
          bq=(dble(l+2*nnn)-0.5d0-zzq)/dsq(nnn,l)*guerpq &
              &              -dsq(nnn-1,l)/dsq(nnn,l)*aq
          aq=guerpq
          guerpq=bq
        end do
        if (guerpq>=0.0d0) then
          sig=1.0d0
        else
          sig=-1.0d0
        endif
        guerpq=log(abs(guerpq))+log(scale)
        wavel=wavel+guerpq
        wave=sig*exp(wavel)
      else
        guerp=dexp(0.5d0*(dlog(2.d0)-gamal(2*l+5))) &
            &           *(dble(l)+1.5d0-zz)
        a=dexp(0.5d0*(dlog(2.d0)-gamal(2*l+3)))
        do nnn=2,n
          b=(dble(l+2*nnn)-0.5d0-zz)/dsq(nnn,l)*guerp &
              &              -dsq(nnn-1,l)/dsq(nnn,l)*a
          a=guerp
          guerp=b
        end do
        wave=dexp(wavel)*guerp
      endif
    endif
  end subroutine waver

  subroutine init_ho_function(nmax)
    implicit none
    integer,intent(in) :: nmax
    integer ::i1,il,in, maxgam

    maxgam = 2 * 400 + 3
    ! maxgam = 2 * nrmax + 3
    if (allocated(dsq)) deallocate(dsq)
    if (allocated(gamal)) deallocate(gamal)
    allocate(dsq(nmax,0:nmax))
    allocate(gamal(maxgam))
    gamal(2)=0.d0
    gamal(1)=0.5d0*dlog(3.14159265358979312d0)
    do i1=3,maxgam
      gamal(i1)=dlog(dble(i1)/2.d0-1.d0)+gamal(i1-2)
    end do
    do il=0,nmax
      do in=1,nmax
        dsq(in,il)=dsqrt(dble(in)*(dble(il+in)+0.5d0))
      end do
    end do
  end subroutine init_ho_function

  subroutine fin_ho_function()
    deallocate(dsq)
    deallocate(gamal)
  end subroutine fin_ho_function

  subroutine init_dbinomial_triangle()
    ! cache for tuning 6-j coefficient
    integer :: n, m, i, j, k, info
    real(8) :: d, a

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
    !!$omp parallel do private( i, j, k, d, info ) schedule (dynamic)
    !do k = 0, n_trinomial
    !  do j = 0, n_trinomial
    !    do i = 0, n_trinomial
    !      dtrinomial(i, j, k) = dtrinomial_func(i, j, k)
    !    end do
    !  end do
    !end do

    f_mb(1)=0.
    g_mb(1)=LOG(0.5D0)
    w_mb(1)=0.
    DO i=2,200
       a=dble(i-1)
       f_mb(i)=f_mb(i-1)+LOG(a)
       g_mb(i)=g_mb(i-1)+LOG(a+0.5D0)
       w_mb(i)=LOG(a+a+1.)
    ENDDO

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

  function tjs(j1, j2, j3, m1, m2, m3) result(r)
    real(8) :: r
    integer, intent(in) :: j1, j2, j3, m1, m2, m3

    r = (-1.d0) ** ((j1 - j2 - m3)/2) * &
      & dcg(j1, m1, j2, m2, j3, -m3) / hat(j3)
  end function tjs

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

  real(8) function gmosh_old (n,l,nc,lc,n1,l1,n2,l2,lr,dd) result(gmosh)
    implicit none
    integer, intent(in) :: n,l,nc,lc,n1,l1,n2,l2,lr
    real(8), intent(in) :: dd
    integer :: ip,ixf,ix, iyi, iyf, j1f,j2,k1i,k1f,m1f,iy,m2f,k2, &
         m2,m2i,m1,j1,k2f,k2i,k1
    real(8) :: dl, d1l, bb, ba, anorm, y, p, bc, cfac, bm , &
         sm, s, sxy, bxy, d

    d = 1.d0/dd
    gmosh=0.
    if(n+n+nc+nc+l+lc-n1-n1-n2-n2-l1-l2 /= 0 ) return
    if(l+lc-lr < 0 ) return
    if(l1+l2-lr < 0 ) return
    if(abs(l-lc)-lr > 0 ) return
    if(abs(l1-l2)-lr > 0 ) return
    dl=log(d)
    d1l=log(d+1.)
    bb=f_mb(n1+1)+f_mb(n2+1)+f_mb(n+1)-f_mb(nc+1)+ &
         g_mb(n1+l1+1)+g_mb(n2+l2+1) &
         -g_mb(n+l+1)-g_mb(nc+lc+1)
    ba=w_mb(l1+1)+w_mb(l2+1)+w_mb(lc+1)+w_mb(l+1)+ &
         f_mb(l1+l2-lr+1)+f_mb(l+lc+lr+2) &
         +f_mb(l+lc-lr+1)+f_mb(lc+lr-l+1)+ &
         f_mb(lr+l-lc+1)-f_mb(l1+l2+lr+2) &
         -f_mb(l1+lr-l2+1)-f_mb(l2+lr-l1+1)-dble(l)*d1l
    ip=lr+n+n1+n2
    p=1+2*(ip/2*2-ip)
    anorm=p*exp(0.5d0*(bb+ba))
    y=0.
    j1f=l+1
    do j1=1,j1f
       j2=l+2-j1
       k1i=abs(l1-j1+1)+1
       k1f=l1+j1
       do k1=k1i,k1f,2
          m1f=n1-(j1+k1-l1)/2+2
          if(m1f-1 < 0 ) cycle
          k2i=max(abs(l2-j2+1),abs(lc-k1+1))+1
          k2f=min(l2+j2,lc+k1)
          if(k2i-k2f > 0 ) cycle
          do k2=k2i,k2f,2
             m2f=n2-(j2+k2-l2)/2+2
             if(m2f-1 < 0 ) cycle
             ip=j2-1+(l1+k1+j1+l2+k2+j2)/2
             p=1+2*(ip/2*2-ip)
             bc=0.5d0*(dble(k1+j2-2)*dl-dble(k1+k2-2)*d1l) &
                  +f_mb(k1+l1-j1+1)+f_mb(k1+k2-lc-1)+ &
                  f_mb(k2+l2-j2+1)-f_mb(k1+l1+j1)-f_mb(k1+k2+lc)- &
                  f_mb(k2+l2+j2)+w_mb(k1)+w_mb(k2)+f_mb((k1+l1+j1)/2)+ &
                  f_mb((k1+k2+lc)/2)+f_mb((k2+l2+j2)/2)- &
                  f_mb((k1+l1-j1)/2+1)-f_mb((l1+j1-k1)/2+1)- &
                  f_mb((j1+k1-l1)/2)-f_mb((k1+k2-lc)/2)- &
                  f_mb((k2+lc-k1)/2+1)-f_mb((lc+k1-k2)/2+1) &
                  -f_mb((k2+l2-j2)/2+1)-f_mb((l2+j2-k2)/2+1)- &
                  f_mb((j2+k2-l2)/2)
             cfac=p*exp(bc)
             sxy=0.
             ixf=min(k1+k1,k1+k2-lc)-1
             do ix=1,ixf
                iyi=max(1,ix+j1+l2-k1-lr)
                iyf=min(l2+l2+1,l1+l2-lr+1,l2+lc+ix-k1-j2+2)
                if(iyi-iyf > 0 ) cycle
                do iy=iyi,iyf
                   ip=ix+iy
                   p=1+2*(ip/2*2-ip)
                   bxy=f_mb(k1+k1-ix)+f_mb(l2+l2-iy+2)+ &
                        f_mb(k2+lc-k1+ix)+f_mb(l1+lr-l2+iy) &
                        -f_mb(ix)-f_mb(iy)-f_mb(k1+k2-lc-ix)- &
                        f_mb(l1+l2-lr-iy+2)-f_mb(k1-l2+lr-j1+iy-ix+1)- &
                        f_mb(l2-k1+lc-j2+ix-iy+3)
                   sxy=sxy+p*exp(bxy)
                enddo
             enddo
             s=cfac*sxy
             sm=0.
             do m1=1,m1f
                m2i=max(1,nc-m1-(k1+k2-lc)/2+3)
                if(m2i-m2f > 0 ) cycle
                do m2=m2i,m2f
                   ip=m1+m2
                   p=1+2*(ip/2*2-ip)
                   bm=dble(m1-1)*dl-dble(m1+m2-2)*d1l+g_mb(1) &
                        +g_mb(m1+m2+(k1+k2+lc)/2-2)-g_mb(k1+m1-1)- &
                        g_mb(k2+m2-1)+f_mb(m1+m2+(k1+k2-lc)/2-2)- &
                        f_mb(m1)-f_mb(m2)-f_mb(n1-m1-(j1+k1-l1)/2+3)- &
                        f_mb(n2-m2-(j2+k2-l2)/2+3) &
                        -f_mb(m1+m2-nc+(k1+k2-lc)/2-2)
                   sm=sm+p*exp(bm)
                enddo
             enddo
             y=y+s*sm
          enddo
       enddo
    enddo
    gmosh=anorm*y
    !gmosh = gmosh * (-1.d0) ** (lc-lr)
    !gmosh = gmosh * (-1.d0) ** (l1-lr)
    gmosh = gmosh * (-1.d0) ** (l+lc+l1)
    !gmosh = gmosh * (-1.d0) ** (l+lc+l2)

  end function gmosh_old

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

  function dtrinomial_func(i, j, k) result(s)
    !
    !  trinomial coefficient: i!! / j!! / k!!
    !
    integer, intent(in) :: i, j, k
    real(8) :: s, di, dj, dk
    integer :: m, l

    s = 1.d0
    m = max(i, j, k)
    if(m == 0) return
    if(m > n_trinomial) then

      write(*,'(a)') 'in trinomial_func, index is too large'

    end if

    di = 1.d0
    do l = i, 0, -2
      if(l == 0) exit
      di = di * dble(l)
    end do

    dj = 1.d0
    do l = j, 0, -2
      if(l == 0) exit
      dj = dj * dble(l)
    end do

    dk = 1.d0
    do l = k, 0, -2
      if(l == 0) exit
      dk = dk * dble(l)
    end do
    s = di / (dj * dk)
  end function dtrinomial_func

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
    if (js > n_triag*3) then
      write(*,'(1a, 1i5, 1a)') '[triangle]: j1+j2+j3 =', js, ' is too large'
      stop
    end if

    js = js / 2
    k1 = k1 / 2

    if (js > n_dbinomial) stop 'increase n_dbinomial'

    delta = 1.0d0 / &
        & (sqrt(dbinomial(js, j3)) * sqrt(dbinomial(j3, k1)) * sqrt(js+1.0d0))

  end subroutine triangle

  function hat(j) result(r)
    real(8) :: r
    integer, intent(in) :: j
    r = dsqrt(dble(j + 1))
  end function hat

end module CommonLibrary

