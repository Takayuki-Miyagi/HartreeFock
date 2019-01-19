module MBPT
  use omp_lib
  use ModelSpace
  use Operators
  implicit none

  public :: MBPTEnergy

  private :: CalcEnergyCorr
  private :: energy_second
  private :: denom
  private :: Eket
  private :: Pari
  private :: Tz
  private :: energy_third_pp
  private :: energy_third_hh
  private :: energy_third_ph
  private :: cross_couple

  type :: MBPTEnergy
    real(8) :: e_0 = 0.d0
    real(8) :: e_2 = 0.d0
    real(8) :: e_3 = 0.d0
    real(8) :: e_3_pp = 0.d0
    real(8) :: e_3_hh = 0.d0
    real(8) :: e_3_ph = 0.d0
  contains
    procedure :: calc => CalcEnergyCorr
    procedure :: energy_second
    procedure :: energy_third
    procedure :: energy_third_pp
    procedure :: energy_third_hh
    procedure :: energy_third_ph
  end type MBPTEnergy
contains
  subroutine CalcEnergyCorr(this,ms,hamil)
    class(MBPTEnergy), intent(inout) :: this
    type(MSpace), intent(in) :: ms
    type(Op), intent(in) :: hamil

    write(*,*)
    write(*,'(a)') " Many-body perturbation calculation up to 3rd order"
    write(*,*)

    this%e_0 = hamil%zero
    call this%energy_second(ms,hamil)
    write(*,'(a,f18.8)') "Second order correction: ", this%e_2

    call this%energy_third(ms,hamil)
    write(*,'(a,4f18.8)') "Third order corrections pp, hh, ph, and total: ", &
        & this%e_3_pp, this%e_3_hh, this%e_3_ph, this%e_3
    write(*,'(a,f18.8)') "Third order MBPT energy: ", this%e_0 + this%e_2 + this%e_3

  end subroutine CalcEnergyCorr

  subroutine energy_second(this,ms,h)
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
    ! <ij||ab> <ab||ij> / 4 denominator
    use Profiler, only: timer
    class(MBPTEnergy), intent(inout) :: this
    type(MSPace), intent(in) :: ms
    type(Op), intent(in) :: h
    integer, allocatable :: p1(:), p2(:), h1(:), h2(:)
    integer :: a, b, i, j
    integer :: ja, jb, ji, jj, J2
    integer :: nidx, idx
    integer :: cnt
    real(8) :: vsum, v, norm, ti

    ti = omp_get_wtime()
    cnt = 0
    do a = 1, size(ms%particles)
      do b = 1, size(ms%particles)
        do i = 1, size(ms%holes)
          do j= 1, size(ms%holes)
            if(Pari(ms%sps,ms%particles(a),ms%particles(b)) /= &
                & Pari(ms%sps,ms%holes(i),ms%holes(j))) cycle
            if(Tz(ms%sps,  ms%particles(a),ms%particles(b)) /= &
                & Tz(ms%sps,ms%holes(i),ms%holes(j))) cycle
            if(Eket(ms%sps,ms%particles(a),ms%particles(b)) > ms%e2max) cycle
            if(Eket(ms%sps,ms%holes(    i),ms%holes(    j)) > ms%e2max) cycle
            cnt = cnt+1
          end do
        end do
      end do
    end do
    nidx = cnt
    allocate(p1(nidx))
    allocate(p2(nidx))
    allocate(h1(nidx))
    allocate(h2(nidx))
    cnt = 0
    do a = 1, size(ms%particles)
      do b = 1, size(ms%particles)
        do i = 1, size(ms%holes)
          do j= 1, size(ms%holes)
            if(Pari(ms%sps,ms%particles(a),ms%particles(b)) /= &
                & Pari(ms%sps,ms%holes(i),ms%holes(j))) cycle
            if(Tz(ms%sps,  ms%particles(a),ms%particles(b)) /= &
                & Tz(ms%sps,ms%holes(i),ms%holes(j))) cycle
            if(Eket(ms%sps,ms%particles(a),ms%particles(b)) > ms%e2max) cycle
            if(Eket(ms%sps,ms%holes(    i),ms%holes(    j)) > ms%e2max) cycle
            cnt = cnt+1
            p1(cnt) = ms%particles(a)
            p2(cnt) = ms%particles(b)
            h1(cnt) = ms%holes(    i)
            h2(cnt) = ms%holes(    j)
          end do
        end do
      end do
    end do




    vsum = 0.d0
    !$omp parallel
    !$omp do private(idx, a, b, i, j, ja, jb, ji, jj, norm, &
    !$omp &  v, J2) reduction(+:vsum)
    do idx = 1, nidx
      a = p1(idx)
      b = p2(idx)
      i = h1(idx)
      j = h2(idx)
      ja = ms%sps%orb(a)%j
      jb = ms%sps%orb(b)%j
      ji = ms%sps%orb(i)%j
      jj = ms%sps%orb(j)%j

      norm = 1.d0
      if(a==b) norm = norm * 2.d0
      if(i==j) norm = norm * 2.d0

      v = 0.d0
      do J2 = max(abs(ja-jb),abs(ji-jj))/2, min(ja+jb,ji+jj)/2
        if(a==b .and. mod(J2,2)==1) cycle
        if(i==j .and. mod(J2,2)==1) cycle
        v = v + dble(2*J2+1) * h%two%GetTwBME(ms%sps,ms%two,a,b,i,j,J2) **2
      end do
      vsum = vsum + v * norm * 0.25d0 / denom(ms%sps,ms%one,h%one,i,j,a,b)
    end do
    !$omp end do
    !$omp end parallel
    this%e_2 = vsum

    deallocate(p1,p2,h1,h2)
    call timer%Add("Second order MBPT",omp_get_wtime()-ti)
  end subroutine energy_second

  function denom(sps,one,f,h1,h2,p1,p2) result(r)
    type(Orbits), intent(in) :: sps
    type(OneBodySpace), intent(in) :: one
    type(NBodyPart), intent(in) :: f
    integer, intent(in) :: h1, h2, p1, p2
    real(8) :: r

    r = f%GetOBME(sps,one,h1,h1) + f%GetOBME(sps,one,h2,h2) - &
        & f%GetOBME(sps,one,p1,p1) - f%GetOBME(sps,one,p2,p2)
  end function denom

  function Eket(sps,a,b)
    type(Orbits), intent(in) :: sps
    integer, intent(in) :: a, b
    integer :: Eket
    Eket = sps%orb(a)%e + sps%orb(b)%e
  end function Eket

  function Pari(sps,a,b) result(p)
    type(Orbits), intent(in) :: sps
    integer, intent(in) :: a, b
    integer :: p
    p = (-1) ** (sps%orb(a)%l + sps%orb(b)%l)
  end function Pari

  function Tz(sps,a,b)
    type(Orbits), intent(in) :: sps
    integer, intent(in) :: a, b
    integer :: Tz
    Tz = sps%orb(a)%z + sps%orb(b)%z
  end function Tz

  subroutine energy_third(this,ms,h)
    class(MBPTEnergy), intent(inout) :: this
    type(MSPace), intent(in) :: ms
    type(Op), intent(in) :: h

    call this%energy_third_hh(ms,h)
    call this%energy_third_pp(ms,h)
    call this%energy_third_ph(ms,h)
    this%e_3 = this%e_3_hh + this%e_3_pp + this%e_3_ph
  end subroutine energy_third

  subroutine energy_third_pp(this,ms,h)
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
    use Profiler, only: timer
    class(MBPTEnergy), intent(inout) :: this
    type(MSPace), intent(in) :: ms
    type(Op), intent(in) :: h
    integer, allocatable :: p1(:), p2(:), p3(:), p4(:)
    integer, allocatable :: h1(:), h2(:)
    integer :: np, nh, nidx, idx, cnt
    integer :: a, b, c, d, i, j
    integer :: ja, jb, jc, jd, ji, jj
    integer :: J2
    real(8) :: v, vsum, norm, ti

    ti = omp_get_wtime()
    call timer%cmemory()

    nh = size(ms%holes)
    np = size(ms%particles)

    nidx = nh**2 * np**4
    cnt = 0
    do a = 1, np
      do b = 1, np
        do c = 1, np
          do d = 1, np
            do i = 1, nh
              do j = 1, nh
                if(Pari(ms%sps,ms%holes(i),ms%holes(j)) /= &
                    & Pari(ms%sps,ms%particles(a),ms%particles(b))) cycle
                if(Pari(ms%sps,ms%particles(a),ms%particles(b)) /= &
                    & Pari(ms%sps,ms%particles(c),ms%particles(d))) cycle
                if(Pari(ms%sps,ms%particles(c),ms%particles(d)) /= &
                    & Pari(ms%sps,ms%holes(i),ms%holes(j))) cycle
                if(Tz(ms%sps,ms%holes(i),ms%holes(j)) /= &
                    & Tz(ms%sps,ms%particles(a),ms%particles(b))) cycle
                if(Tz(ms%sps,ms%particles(a),ms%particles(b)) /= &
                    & Tz(ms%sps,ms%particles(c),ms%particles(d))) cycle
                if(Tz(ms%sps,ms%particles(c),ms%particles(d)) /= &
                    & Tz(ms%sps,ms%holes(i),ms%holes(j))) cycle
                if(Eket(ms%sps,ms%holes(i),ms%holes(j)) > ms%e2max) cycle
                if(Eket(ms%sps,ms%particles(a),ms%particles(b)) > ms%e2max) cycle
                if(Eket(ms%sps,ms%particles(c),ms%particles(d)) > ms%e2max) cycle
                cnt = cnt + 1
              end do
            end do
          end do
        end do
      end do
    end do
    nidx = cnt
    allocate(p1(nidx))
    allocate(p2(nidx))
    allocate(p3(nidx))
    allocate(p4(nidx))
    allocate(h1(nidx))
    allocate(h2(nidx))
    cnt = 0
    do a = 1, np
      do b = 1, np
        do c = 1, np
          do d = 1, np
            do i = 1, nh
              do j = 1, nh
                if(Pari(ms%sps,ms%holes(i),ms%holes(j)) /= &
                    & Pari(ms%sps,ms%particles(a),ms%particles(b))) cycle
                if(Pari(ms%sps,ms%particles(a),ms%particles(b)) /= &
                    & Pari(ms%sps,ms%particles(c),ms%particles(d))) cycle
                if(Pari(ms%sps,ms%particles(c),ms%particles(d)) /= &
                    & Pari(ms%sps,ms%holes(i),ms%holes(j))) cycle
                if(Tz(ms%sps,ms%holes(i),ms%holes(j)) /= &
                    & Tz(ms%sps,ms%particles(a),ms%particles(b))) cycle
                if(Tz(ms%sps,ms%particles(a),ms%particles(b)) /= &
                    & Tz(ms%sps,ms%particles(c),ms%particles(d))) cycle
                if(Tz(ms%sps,ms%particles(c),ms%particles(d)) /= &
                    & Tz(ms%sps,ms%holes(i),ms%holes(j))) cycle
                if(Eket(ms%sps,ms%holes(i),ms%holes(j)) > ms%e2max) cycle
                if(Eket(ms%sps,ms%particles(a),ms%particles(b)) > ms%e2max) cycle
                if(Eket(ms%sps,ms%particles(c),ms%particles(d)) > ms%e2max) cycle
                cnt = cnt + 1
                p1(cnt) = ms%particles(a)
                p2(cnt) = ms%particles(b)
                p3(cnt) = ms%particles(c)
                p4(cnt) = ms%particles(d)
                h1(cnt) = ms%holes(    i)
                h2(cnt) = ms%holes(    j)
              end do
            end do
          end do
        end do
      end do
    end do

    !$omp parallel
    !$omp do private(idx, a, b, c, d, i, j, &
    !$omp &  ja, jb, jc, jd, ji, jj, norm, v, J2) reduction(+:vsum)
    do idx = 1, nidx
      a = p1(idx)
      b = p2(idx)
      c = p3(idx)
      d = p4(idx)
      i = h1(idx)
      j = h2(idx)

      ja = ms%sps%orb(a)%j
      jb = ms%sps%orb(b)%j
      jc = ms%sps%orb(c)%j
      jd = ms%sps%orb(d)%j
      ji = ms%sps%orb(i)%j
      jj = ms%sps%orb(j)%j
      norm = 1.d0
      if(a==b) norm = norm * 2.d0
      if(c==d) norm = norm * 2.d0
      if(i==j) norm = norm * 2.d0

      v = 0.d0
      do J2 = max(abs(ja-jb),abs(jc-jd),abs(ji-jj))/2, &
            & min(   (ja+jb),   (jc+jd),   (ji+jj))/2
        if(a==b .and. mod(J2,2)==1) cycle
        if(c==d .and. mod(J2,2)==1) cycle
        if(i==j .and. mod(J2,2)==1) cycle
        v = v + dble(2*J2+1) * h%two%GetTwBME(ms%sps,ms%two,i,j,a,b,J2) * &
            & h%two%GetTwBME(ms%sps,ms%two,a,b,c,d,J2) * &
            & h%two%GetTwBME(ms%sps,ms%two,c,d,i,j,J2)
      end do
      vsum = vsum + v * norm * 0.125d0 / &
          & ( denom(ms%sps,ms%one,h%one,i,j,a,b) * denom(ms%sps,ms%one,h%one,i,j,c,d) )
    end do
    !$omp end do
    !$omp end parallel

    this%e_3_pp = vsum

    call timer%countup_memory('temporary array for MBPT pp ladder')
    deallocate(p1,p2,p3,p4,h1,h2)
    call timer%Add("Third order MBPT pp ladder",omp_get_wtime()-ti)
  end subroutine energy_third_pp

  subroutine energy_third_hh(this,ms,h)
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
    use Profiler, only: timer
    class(MBPTEnergy), intent(inout) :: this
    type(MSPace), intent(in) :: ms
    type(Op), intent(in) :: h
    integer, allocatable :: p1(:), p2(:)
    integer, allocatable :: h1(:), h2(:), h3(:), h4(:)
    integer :: np, nh, nidx, idx, cnt
    integer :: a, b, i, j, k, l
    integer :: ja, jb, ji, jj, jk, jl
    integer :: J2
    real(8) :: v, vsum, norm, ti

    ti = omp_get_wtime()
    call timer%cmemory()

    nh = size(ms%holes)
    np = size(ms%particles)

    cnt = 0
    do a = 1, np
      do b = 1, np
        do i = 1, nh
          do j = 1, nh
            do k = 1, nh
              do l = 1, nh
                if(Pari(ms%sps,ms%holes(i),ms%holes(j)) /= &
                    & Pari(ms%sps,ms%particles(a),ms%particles(b))) cycle
                if(Pari(ms%sps,ms%holes(k),ms%holes(l)) /= &
                    & Pari(ms%sps,ms%holes(i),ms%holes(j))) cycle
                if(Pari(ms%sps,ms%particles(a),ms%particles(b)) /= &
                    & Pari(ms%sps,ms%holes(k),ms%holes(l))) cycle
                if(Tz(ms%sps,ms%holes(i),ms%holes(j)) /= &
                    & Tz(ms%sps,ms%particles(a),ms%particles(b))) cycle
                if(Tz(ms%sps,ms%holes(k),ms%holes(l)) /= &
                    & Tz(ms%sps,ms%holes(i),ms%holes(j))) cycle
                if(Tz(ms%sps,ms%particles(a),ms%particles(b)) /= &
                    & Tz(ms%sps,ms%holes(k),ms%holes(l))) cycle
                if(Eket(ms%sps,a,b) > ms%e2max) cycle
                if(Eket(ms%sps,i,j) > ms%e2max) cycle
                if(Eket(ms%sps,k,l) > ms%e2max) cycle
                cnt = cnt + 1
              end do
            end do
          end do
        end do
      end do
    end do

    nidx = cnt
    allocate(p1(nidx))
    allocate(p2(nidx))
    allocate(h1(nidx))
    allocate(h2(nidx))
    allocate(h3(nidx))
    allocate(h4(nidx))
    cnt = 0
    do a = 1, np
      do b = 1, np
        do i = 1, nh
          do j = 1, nh
            do k = 1, nh
              do l = 1, nh
                if(Pari(ms%sps,ms%holes(i),ms%holes(j)) /= &
                    & Pari(ms%sps,ms%particles(a),ms%particles(b))) cycle
                if(Pari(ms%sps,ms%holes(k),ms%holes(l)) /= &
                    & Pari(ms%sps,ms%holes(i),ms%holes(j))) cycle
                if(Pari(ms%sps,ms%particles(a),ms%particles(b)) /= &
                    & Pari(ms%sps,ms%holes(k),ms%holes(l))) cycle
                if(Tz(ms%sps,ms%holes(i),ms%holes(j)) /= &
                    & Tz(ms%sps,ms%particles(a),ms%particles(b))) cycle
                if(Tz(ms%sps,ms%holes(k),ms%holes(l)) /= &
                    & Tz(ms%sps,ms%holes(i),ms%holes(j))) cycle
                if(Tz(ms%sps,ms%particles(a),ms%particles(b)) /= &
                    & Tz(ms%sps,ms%holes(k),ms%holes(l))) cycle
                if(Eket(ms%sps,ms%particles(a),ms%particles(b)) > ms%e2max) cycle
                if(Eket(ms%sps,ms%holes(i),ms%holes(j)) > ms%e2max) cycle
                if(Eket(ms%sps,ms%holes(k),ms%holes(l)) > ms%e2max) cycle
                cnt = cnt + 1
                p1(cnt) = ms%particles(a)
                p2(cnt) = ms%particles(b)
                h1(cnt) = ms%holes(    i)
                h2(cnt) = ms%holes(    j)
                h3(cnt) = ms%holes(    k)
                h4(cnt) = ms%holes(    l)
              end do
            end do
          end do
        end do
      end do
    end do

    !$omp parallel
    !$omp do private(idx, a, b, i, j, k, l, &
    !$omp &  ja, jb, ji, jj, jk, jl, norm, v, J2) reduction(+:vsum)
    do idx = 1, nidx
      a = p1(idx)
      b = p2(idx)
      i = h1(idx)
      j = h2(idx)
      k = h3(idx)
      l = h4(idx)

      ja = ms%sps%orb(a)%j
      jb = ms%sps%orb(b)%j
      ji = ms%sps%orb(i)%j
      jj = ms%sps%orb(j)%j
      jk = ms%sps%orb(k)%j
      jl = ms%sps%orb(l)%j
      norm = 1.d0
      if(a==b) norm = norm * 2.d0
      if(i==j) norm = norm * 2.d0
      if(k==l) norm = norm * 2.d0

      v = 0.d0
      do J2 = max(abs(ja-jb),abs(ji-jj),abs(jk-jl))/2, &
            & min(   (ja+jb),   (ji+jj),   (jk+jl))/2
        if(a==b .and. mod(J2,2)==1) cycle
        if(i==j .and. mod(J2,2)==1) cycle
        if(k==l .and. mod(J2,2)==1) cycle
        v = v + dble(2*J2+1) * h%two%GetTwBME(ms%sps,ms%two,i,j,a,b,J2) * &
            & h%two%GetTwBME(ms%sps,ms%two,k,l,i,j,J2) * &
            & h%two%GetTwBME(ms%sps,ms%two,a,b,k,l,J2)
      end do
      vsum = vsum + v * norm * 0.125d0 / &
          & ( denom(ms%sps,ms%one,h%one,i,j,a,b) * denom(ms%sps,ms%one,h%one,k,l,a,b) )
    end do
    !$omp end do
    !$omp end parallel

    this%e_3_hh = vsum

    call timer%countup_memory('temporary array for MBPT hh ladder')
    deallocate(p1,p2,h1,h2,h3,h4)
    call timer%Add("Third order MBPT hh ladder",omp_get_wtime()-ti)
  end subroutine energy_third_hh

  subroutine energy_third_ph(this,ms,h)
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
    use Profiler, only: timer
    class(MBPTEnergy), intent(inout) :: this
    type(MSPace), intent(in) :: ms
    type(Op), intent(in) :: h
    integer, allocatable :: p1(:), p2(:), p3(:)
    integer, allocatable :: h1(:), h2(:), h3(:)
    integer :: np, nh, nidx, idx, cnt
    integer :: a, b, c, i, j, k
    integer :: ja, jb, jc, ji, jj, jk
    integer :: J2
    real(8) :: v, vsum, norm, ti

    ti = omp_get_wtime()
    call timer%cmemory()

    nh = size(ms%holes)
    np = size(ms%particles)

    cnt = 0
    do a = 1, np
      do b = 1, np
        do c = 1, np
          do i = 1, nh
            do j = 1, nh
              do k = 1, nh
                if(Pari(ms%sps,ms%holes(i),ms%holes(j)) /= &
                    & Pari(ms%sps,ms%particles(a),ms%particles(b))) cycle
                if(Pari(ms%sps,ms%holes(k),ms%particles(b)) /= &
                    & Pari(ms%sps,ms%holes(i),ms%particles(c))) cycle
                if(Pari(ms%sps,ms%particles(a),ms%particles(c)) /= &
                    & Pari(ms%sps,ms%holes(k),ms%holes(j))) cycle
                if(Tz(ms%sps,ms%holes(i),ms%holes(j)) /= &
                    & Tz(ms%sps,ms%particles(a),ms%particles(b))) cycle
                if(Tz(ms%sps,ms%holes(k),ms%particles(b)) /= &
                    & Tz(ms%sps,ms%holes(i),ms%particles(c))) cycle
                if(Tz(ms%sps,ms%particles(a),ms%particles(c)) /= &
                    & Tz(ms%sps,ms%holes(k),ms%holes(j))) cycle
                if(Eket(ms%sps,ms%particles(a),ms%particles(b)) > ms%e2max) cycle
                if(Eket(ms%sps,ms%holes(    i),ms%holes(    j)) > ms%e2max) cycle
                if(Eket(ms%sps,ms%holes(    k),ms%particles(b)) > ms%e2max) cycle
                if(Eket(ms%sps,ms%holes(    i),ms%particles(c)) > ms%e2max) cycle
                if(Eket(ms%sps,ms%particles(a),ms%particles(c)) > ms%e2max) cycle
                if(Eket(ms%sps,ms%holes(    k),ms%holes(    j)) > ms%e2max) cycle
                cnt = cnt + 1
              end do
            end do
          end do
        end do
      end do
    end do

    nidx = cnt
    allocate(p1(nidx))
    allocate(p2(nidx))
    allocate(p3(nidx))
    allocate(h1(nidx))
    allocate(h2(nidx))
    allocate(h3(nidx))

    cnt = 0
    do a = 1, np
      do b = 1, np
        do c = 1, np
          do i = 1, nh
            do j = 1, nh
              do k = 1, nh
                if(Pari(ms%sps,ms%holes(i),ms%holes(j)) /= &
                    & Pari(ms%sps,ms%particles(a),ms%particles(b))) cycle
                if(Pari(ms%sps,ms%holes(k),ms%particles(b)) /= &
                    & Pari(ms%sps,ms%holes(i),ms%particles(c))) cycle
                if(Pari(ms%sps,ms%particles(a),ms%particles(c)) /= &
                    & Pari(ms%sps,ms%holes(k),ms%holes(j))) cycle
                if(Tz(ms%sps,ms%holes(i),ms%holes(j)) /= &
                    & Tz(ms%sps,ms%particles(a),ms%particles(b))) cycle
                if(Tz(ms%sps,ms%holes(k),ms%particles(b)) /= &
                    & Tz(ms%sps,ms%holes(i),ms%particles(c))) cycle
                if(Tz(ms%sps,ms%particles(a),ms%particles(c)) /= &
                    & Tz(ms%sps,ms%holes(k),ms%holes(j))) cycle
                if(Eket(ms%sps,ms%particles(a),ms%particles(b)) > ms%e2max) cycle
                if(Eket(ms%sps,ms%holes(    i),ms%holes(    j)) > ms%e2max) cycle
                if(Eket(ms%sps,ms%holes(    k),ms%particles(b)) > ms%e2max) cycle
                if(Eket(ms%sps,ms%holes(    i),ms%particles(c)) > ms%e2max) cycle
                if(Eket(ms%sps,ms%particles(a),ms%particles(c)) > ms%e2max) cycle
                if(Eket(ms%sps,ms%holes(    k),ms%holes(    j)) > ms%e2max) cycle
                cnt = cnt + 1
                p1(cnt) = ms%particles(a)
                p2(cnt) = ms%particles(b)
                p3(cnt) = ms%particles(c)
                h1(cnt) = ms%holes(    i)
                h2(cnt) = ms%holes(    j)
                h3(cnt) = ms%holes(    k)
              end do
            end do
          end do
        end do
      end do
    end do

    !$omp parallel
    !$omp do private(idx, a, b, c, i, j, k, &
    !$omp &  ja, jb, jc, ji, jj, jk, norm, v, J2) reduction(+:vsum)
    do idx = 1, nidx
      a = p1(idx)
      b = p2(idx)
      c = p3(idx)
      i = h1(idx)
      j = h2(idx)
      k = h3(idx)

      ja = ms%sps%orb(a)%j
      jb = ms%sps%orb(b)%j
      jc = ms%sps%orb(c)%j
      ji = ms%sps%orb(i)%j
      jj = ms%sps%orb(j)%j
      jk = ms%sps%orb(k)%j

      norm = 1.d0
      if(i==j) norm = norm * sqrt(2.d0)
      if(j==k) norm = norm * sqrt(2.d0)
      if(a==b) norm = norm * sqrt(2.d0)
      if(a==c) norm = norm * sqrt(2.d0)

      v = 0.d0
      do J2 = max(abs(ji-jb),abs(ja-jj),abs(jk-jc))/2, &
            & min(   (ji+jb),   (ja+jj),   (jk+jc))/2
        v = v + dble(2*J2+1) * &
            & cross_couple(ms%sps,ms%two,h%two,i,j,a,b,J2) * &
            & cross_couple(ms%sps,ms%two,h%two,k,b,i,c,J2) * &
            & cross_couple(ms%sps,ms%two,h%two,a,c,k,j,J2)
      end do
      vsum = vsum - v *  norm / &
          & ( denom(ms%sps,ms%one,h%one,k,j,a,c) * denom(ms%sps,ms%one,h%one,i,j,a,b) )
    end do
    !$omp end do
    !$omp end parallel

    this%e_3_ph = vsum
    call timer%countup_memory('temporary array for MBPT ph ladder')
    deallocate(p1,p2,p3,h1,h2,h3)
    call timer%Add("Third order MBPT ph ladder",omp_get_wtime()-ti)
  end subroutine energy_third_ph

  function cross_couple(sps,two,v,i,j,a,b,L) result(r)
    use CommonLibrary, only: sjs
    type(Orbits), intent(in) :: sps
    type(TwoBodySpace), intent(in) :: two
    type(NBodyPart), intent(in) :: v
    integer, intent(in) :: i, j, a, b, L
    integer :: ji, jj, ja, jb
    integer :: J2
    real(8) :: r

    r = 0.d0

    ji = sps%orb(i)%j
    jj = sps%orb(j)%j
    ja = sps%orb(a)%j
    jb = sps%orb(b)%j

    do J2 = max(abs(ji-jb),abs(ja-jj))/2, min((ji+jb),(ja+jj))/2
      r = r + dble(2*J2+1) * &
          & sjs(ji, jj, 2*J2, ja, jb, 2*L) * &
          & v%GetTwBME(sps,two,i,j,a,b,J2)
    end do
  end function cross_couple

end module MBPT
