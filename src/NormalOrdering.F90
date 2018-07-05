module NormalOrdering
  use InputParameters, only: parameters
  use ScalarOperator, only: ScalarOperators, NBodyScalars
  use ModelSpace, only: spo_pn, MSpace, OneBodySpace, TwoBodySpace, ThreeBodySpace, &
      & OneBodyChannel, TwoBodyChannel, ThreeBodyChannel
  use MatrixDouble, only: DMat
  use VectorDouble, only: DVec
  use LinAlgLib
  implicit none
  private :: NormOrd_idx, NormOrd_ket, InitNOThree2Two, FinNOThree2Two, &
      & HFNormOrd_ch, hole, braket, store_idx, InitNOThree2Two_HFbasis, &
      & FinNOThree2Two_HFbasis, Del, ScalarNOThree2Two, &
      & ScalarNOThree2Two_HFBasis, ScalarNOTwo2One, ScalarNOOne2Zero
  public :: NOThree2Two, NOThree2Two_HFbasis, NormOrd

  interface NormOrd
    module procedure :: ScalarNOThree2Two, ScalarNOThree2Two_HFBasis, &
          & ScalarNOTwo2One, ScalarNOOne2Zero
  end interface NormOrd

  !------------------------------
  type :: NormOrd_idx
    integer :: n
    integer, allocatable :: num(:,:)
  end type NormOrd_idx

  type :: NormOrd_ket
    type(NormOrd_idx), allocatable :: ket(:)
  end type NormOrd_ket

  type :: NOThree2Two
    type(NormOrd_ket), allocatable :: jptz(:)
    real(8) :: usedmem = 0.d0
  contains
    procedure :: init => InitNOThree2Two
    procedure :: fin => FinNOThree2Two
  end type NOThree2Two
  !------------------------------


  !------------------------------
  type :: braket
    integer :: n
    integer, allocatable :: num(:,:)
  contains
    procedure :: store_idx
  end type braket

  type :: hole
    integer :: n
    integer, allocatable :: h2label(:), label2h(:)
    type(braket), allocatable :: h(:)
  end type  hole

  type :: HFNormOrd_ch
    type(hole), allocatable :: idx(:)
  end type HFNormOrd_ch

  type :: NOThree2Two_HFbasis
    type(HFNormOrd_ch), allocatable :: jptz(:)
    real(8) :: usedmem
  contains
    procedure :: init => InitNOThree2Two_HFbasis
    procedure :: fin => FinNOThree2Two_HFbasis
  end type NOThree2Two_HFbasis
  !------------------------------

contains

  subroutine InitNOThree2Two(this, params, sps, two, thr)
    class(NOThree2Two), intent(inout) :: this
    type(parameters), intent(in) :: params
    type(spo_pn), intent(in) :: sps
    type(TwoBodySpace), intent(in) :: two
    type(ThreeBodySpace), intent(in) :: thr
    integer :: ich2, j2, p2, tz2, n
    integer :: ich3, j3, p3, tz3, ket
    integer :: e3max, st, a, b, c, loop, num
    integer :: i, jlab, idx
    real(8) :: cnt

    e3max = params%e3max
    allocate(this%jptz(two%n))
    cnt = 0.d0
    do ich2 = 1, two%n
      j2 = two%j(ich2)
      p2 = two%p(ich2)
      tz2 = two%tz(ich2)
      n = two%jptz(ich2)%n
      allocate(this%jptz(ich2)%ket(n))
      !$omp parallel
      !$omp do private(st, a, b, loop, num, c, p3, tz3, j3, ich3, jlab, idx, &
      !$omp &  i, ket, cnt)
      do st = 1, n
        a = two%jptz(ich2)%n2label1(st)
        b = two%jptz(ich2)%n2label2(st)
        this%jptz(ich2)%ket(st)%n = 0
        do loop = 1, 2
          num = 0
          do c = 1, sps%n
            if(abs(sps%NOCoef(c)) < 1.d-4) cycle
            if(sps%nshell(a) + sps%nshell(b) + sps%nshell(c) > e3max) cycle
            p3 = p2 * (-1) ** sps%ll(c)
            tz3 = 2*tz2 + sps%itz(c)
            do j3 = iabs(2*j2 - sps%jj(c)), (2*j2 + sps%jj(c)), 2
              ich3 = thr%jptz2n(j3,p3,tz3)
              if(ich3 == 0) cycle
              call find_cfp(a,b,c,j2,thr%jptz(ich3),jlab,idx)
              if(idx == 0) cycle
              do i = 1, thr%jptz(ich3)%idx(idx)%nphys
                ket = thr%jptz(ich3)%idx(idx)%labels2n(i)
                num = num + 1
                if(loop == 1) cycle

                this%jptz(ich2)%ket(st)%num(1, num) = ich3
                this%jptz(ich2)%ket(st)%num(2, num) = j3
                this%jptz(ich2)%ket(st)%num(3, num) = c
                this%jptz(ich2)%ket(st)%num(4, num) = idx
                this%jptz(ich2)%ket(st)%num(5, num) = jlab
                this%jptz(ich2)%ket(st)%num(6, num) = i
                this%jptz(ich2)%ket(st)%num(7, num) = ket

              end do
            end do
          end do
          if(loop == 1) then
            this%jptz(ich2)%ket(st)%n = num
            allocate(this%jptz(ich2)%ket(st)%num(7, num))
            cnt = cnt + dble(this%jptz(ich2)%ket(st)%n * 7)
          end if
        end do
      end do
      !$omp end do
      !$omp end parallel
    end do
    this%usedmem = cnt * 4.d0 / (1024.d0 ** 3)
  end subroutine InitNOThree2Two

  subroutine FinNOThree2Two(this)
    class(NOThree2Two), intent(inout) :: this
    integer :: ich, st, n
    do ich = 1, size(this%jptz)
      n = size(this%jptz(ich)%ket)
      do st = 1, n
        if(this%jptz(ich)%ket(n)%n < 1) cycle
        deallocate(this%jptz(ich)%ket(st)%num)
      end do
      deallocate(this%jptz(ich)%ket)
    end do
    deallocate(this%jptz)
  end subroutine FinNOThree2Two

  subroutine InitNOThree2Two_HFbasis(this, params, sps, one, two)
    class(NOThree2Two_HFbasis), intent(inout) :: this
    type(parameters), intent(in) :: params
    type(spo_pn), intent(in) :: sps
    type(OneBodySpace), intent(in) :: one
    type(TwoBodySpace), intent(in) :: two
    integer :: ich2, j2, p2, tz2, n, bra2, num
    integer :: a, b, c, e3max, loop
    real(8) :: cnt
    e3max = params%e3max_3nf
    cnt = 0.d0
    allocate(this%jptz(two%n))
    do ich2 = 1, two%n
      j2 = two%j(ich2)
      p2 = two%p(ich2)
      tz2 = two%tz(ich2)
      n = two%jptz(ich2)%n
      allocate(this%jptz(ich2)%idx(n))
      do bra2 = 1, n
        a = two%jptz(ich2)%n2label1(bra2)
        b = two%jptz(ich2)%n2label2(bra2)

        do loop = 1, 2
          num = 0
          do c = 1, sps%n
            if(abs(sps%NOCoef(c)) < 1.d-4) cycle
            if(sps%nshell(a) + sps%nshell(b) + sps%nshell(c) > params%e3max_3nf) cycle
            num = num + 1
            if(loop == 1) cycle
            this%jptz(ich2)%idx(bra2)%label2h(num) = c
            this%jptz(ich2)%idx(bra2)%h2label(c) = num
            call this%jptz(ich2)%idx(bra2)%h(num)%store_idx(params,sps,one,j2,p2,tz2,a,b,c)
            cnt = cnt + dble(this%jptz(ich2)%idx(bra2)%h(num)%n * 9) * 4.d0
          end do
          if(loop == 1) then
            this%jptz(ich2)%idx(bra2)%n = num
            allocate(this%jptz(ich2)%idx(bra2)%h(num))
            allocate(this%jptz(ich2)%idx(bra2)%label2h(num))
            allocate(this%jptz(ich2)%idx(bra2)%h2label(sps%n))
            this%jptz(ich2)%idx(bra2)%label2h = 0
            this%jptz(ich2)%idx(bra2)%h2label = 0
          end if
        end do

      end do
    end do
    this%usedmem = cnt / (1024.d0 ** 3)
  end subroutine InitNOThree2Two_HFbasis

  subroutine FinNOThree2Two_HFbasis(this)
    class(NOThree2Two_HFbasis), intent(inout) :: this
    integer :: ich, bra2, num
    do ich = 1, size(this%jptz)
      do bra2 = 1, size(this%jptz(ich)%idx)
        do num = 1, this%jptz(ich)%idx(bra2)%n
          deallocate(this%jptz(ich)%idx(bra2)%h(num)%num)
        end do
        deallocate(this%jptz(ich)%idx(bra2)%h)
      end do
      deallocate(this%jptz(ich)%idx)
    end do
    deallocate(this%jptz)
  end subroutine FinNOThree2Two_HFbasis

  function ScalarNOThree2Two(params, two, thr, sc3, no) result(sc2)
    type(NBodyScalars) :: sc2
    type(parameters), intent(in) :: params
    type(TwoBodySpace), intent(in) :: two
    type(ThreeBodySpace), intent(in) :: thr
    type(NBodyScalars), intent(in) :: sc3
    type(NOThree2Two), intent(in) :: no
    integer :: e3max, ich2, j2, p2, itz2, n2, j3
    integer :: a, b, c, d, eb, ek, ii, jj
    integer :: bra2, ket2, bra3, ket3
    integer :: ichb, ichk, idxb, idxk
    integer :: ljab, ljcd, i123, i456
    real(8), allocatable :: w2(:,:)
    real(8) :: w12
    call sc2%init(two%SpinParityTz)
    e3max = params%e3max_3nf
    do ich2 = 1, two%n
      j2 = two%j(ich2)
      p2 = two%p(ich2)
      itz2 = two%tz(ich2)
      n2 = two%jptz(ich2)%n
      allocate(w2(n2,n2))
      w2(:,:) = 0.d0
      do bra2 = 1, n2
        a = two%jptz(ich2)%n2label1(bra2)
        b = two%jptz(ich2)%n2label2(bra2)
        do ket2 = 1, bra2
          c = two%jptz(ich2)%n2label1(ket2)
          d = two%jptz(ich2)%n2label2(ket2)
          w12 = 0.d0
          do ii = 1, no%Jptz(ich2)%ket(bra2)%n
            ichb = no%jptz(ich2)%ket(bra2)%num(1,ii)
            j3   = no%jptz(ich2)%ket(bra2)%num(2,ii)
            eb   = no%jptz(ich2)%ket(bra2)%num(3,ii)
            idxb = no%jptz(ich2)%ket(bra2)%num(4,ii)
            ljab = no%jptz(ich2)%ket(bra2)%num(5,ii)
            i123 = no%jptz(ich2)%ket(bra2)%num(6,ii)
            bra3 = no%jptz(ich2)%ket(bra2)%num(7,ii)
            do jj = 1, no%Jptz(ich2)%ket(ket2)%n
              ichk = no%jptz(ich2)%ket(ket2)%num(1,jj)
              j3   = no%jptz(ich2)%ket(ket2)%num(2,jj)
              ek   = no%jptz(ich2)%ket(ket2)%num(3,jj)
              idxk = no%jptz(ich2)%ket(ket2)%num(4,jj)
              ljcd = no%jptz(ich2)%ket(ket2)%num(5,jj)
              i456 = no%jptz(ich2)%ket(ket2)%num(6,jj)
              ket3 = no%jptz(ich2)%ket(ket2)%num(7,jj)
              if(ichb /= ichk) cycle
              if(eb /= ek) cycle
              w12 = w12 + dble(j3 + 1) * sc3%jptz(ichb)%m(bra3, ket3) * &
                  & thr%jptz(ichb)%idx(idxb)%cfp(ljab, i123) * &
                  & thr%jptz(ichb)%idx(idxk)%cfp(ljcd, i456)
            end do
          end do

          w2(bra2, ket2) = 6.d0 * w12 / (dble(2*j2 + 1) * Del(a,b) * Del(c,d))
          w2(ket2, bra2) = w2(bra2, ket2)
        end do
      end do
      sc2%jptz(ich2)%m = w2
      deallocate(w2)
    end do
  end function ScalarNOThree2Two

  function ScalarNOThree2Two_HFBasis(params, sps, HF, two, thbme, no) result(sc2)
    use read_3bme, only: iThreeBodyScalar
#ifdef MPI
    use mpi
    use MPIFunction, only: master_slave, ierr
#else
    use MPIFunction, only: master_slave
#endif
    type(NBodyScalars) :: sc2
    type(parameters), intent(in) :: params
    type(spo_pn), intent(in) :: sps
    type(NBodyScalars), intent(in) :: HF
    type(iThreeBodyScalar), intent(in) :: thbme
    type(TwoBodySpace), intent(in) :: two
    type(NOThree2Two_HFbasis), intent(in) :: no
    integer, allocatable :: slranks(:)
    integer :: e3max
#ifdef MPI
    integer :: ich
#endif
    call sc2%init(two%SpinParityTz)
    e3max = params%e3max_3nf
    call master_slave(insideNO2, two%n, slranks)
#ifdef MPI
    call mpi_bcast(slranks(1), two%ichmax, mpi_integer, 0, mpi_comm_world, ierr)
    do ich = 1, two%ichmax
      n = two%jptz(ich)%n
      call mpi_bcast(twbme%jptz(ich)%m(1,1),n**2,mpi_real8,slranks(ich),mpi_comm_world,ierr)
    end do
#endif
  contains
    subroutine insideNO2(ich2)
      use ScalarOperator, only: get3BMEpn
      integer, intent(in) :: ich2
      integer :: n, m
      integer :: j2, p2, itz2, n2
      integer :: p3, itz3, j3
      integer :: a, b, c, d, e, f
      integer :: bra2, ket2
      integer :: ich1, ehf, eho, fho
      real(8), allocatable :: w2(:,:)
      real(8) :: w12, v
      j2 = two%j(ich2)
      p2 = two%p(ich2)
      itz2 = two%tz(ich2)
      n2 = two%jptz(ich2)%n
      allocate(w2(n2,n2))
      w2(:,:) = 0.d0
      !$omp parallel
      !$omp do private(bra2, a, b, ket2, c, d, w12, n, &
      !$omp & ich1, j3, p3, itz3, e, f, eho, fho, ehf, v) schedule(dynamic)
      do bra2 = 1, n2
        a = two%jptz(ich2)%n2label1(bra2)
        b = two%jptz(ich2)%n2label2(bra2)
        do ket2 = 1, bra2
          c = two%jptz(ich2)%n2label1(ket2)
          d = two%jptz(ich2)%n2label2(ket2)

          w12 = 0.d0
          if(no%jptz(ich2)%idx(bra2)%n < 1) cycle
          do n = 1, no%jptz(ich2)%idx(bra2)%n
            if(no%jptz(ich2)%idx(bra2)%h(n)%n < 1) cycle
            do m = 1, no%jptz(ich2)%idx(bra2)%h(n)%n
              ich1 = no%jptz(ich2)%idx(bra2)%h(n)%num(1, m)
              j3   = no%jptz(ich2)%idx(bra2)%h(n)%num(2, m)
              p3   = no%jptz(ich2)%idx(bra2)%h(n)%num(3, m)
              itz3 = no%jptz(ich2)%idx(bra2)%h(n)%num(4, m)
              e    = no%jptz(ich2)%idx(bra2)%h(n)%num(5, m)
              f    = no%jptz(ich2)%idx(bra2)%h(n)%num(6, m)
              eho  = no%jptz(ich2)%idx(bra2)%h(n)%num(7, m)
              fho  = no%jptz(ich2)%idx(bra2)%h(n)%num(8, m)
              ehf  = no%jptz(ich2)%idx(bra2)%h(n)%num(9, m)
              v = Get3BMEpn(sps,thbme,j3,p3,itz3,a,b,e,j2,c,d,f,j2)
              w12 = w12 + dble(j3 + 1) * v * &
                  & HF%jptz(ich1)%m(eho, ehf) * &
                  & HF%jptz(ich1)%m(fho, ehf)
            end do
          end do
          w2(bra2, ket2) = w12 / (dble(2*j2 + 1) * Del(a,b) * Del(c,d))
          w2(ket2, bra2) = w2(bra2, ket2)
        end do
      end do
      !$omp end do
      !$omp end parallel
      sc2%jptz(ich2)%m = w2
      deallocate(w2)
    end subroutine insideNO2
  end function ScalarNOThree2Two_HFBasis

  function ScalarNOTwo2One(params, sps, one, two, sc2) result(sc1)
    type(NBodyScalars) :: sc1
    type(parameters), intent(in) :: params
    type(spo_pn), intent(in) :: sps
    type(OneBodySpace), intent(in) :: one
    type(TwoBodySpace), intent(in) :: two
    type(NBodyScalars), intent(in) :: sc2
    integer :: ich1, ich2
    integer :: a, b, c, n, bra1, ket1, bra2, ket2
    integer :: iphase
    integer :: j1, p1, itz1
    integer :: j2, p2, itz2
    real(8) :: w1, v2, vv
    call sc1%init(one%SpinParityTz)
    do ich1 = 1, one%n
      j1 = one%j(ich1)
      p1 = one%p(ich1)
      itz1 = one%tz(ich1)
      n = one%jptz(ich1)%n
      do bra1 = 1, n
        a = one%jptz(ich1)%n2label(bra1)
        do ket1 = 1, bra1
          b = one%jptz(ich1)%n2label(ket1)
          w1 = 0.d0
          do c = 1, sps%n
            if(abs(sps%NoCoef(c)) < 1.d-4) cycle
            if(sps%nshell(a) + sps%nshell(c) > params%e2max) cycle
            if(sps%nshell(b) + sps%nshell(c) > params%e2max) cycle
            p2 = p1 * (-1) ** (sps%ll(c))
            itz2 = (itz1 + sps%itz(c)) / 2

            vv = 0.d0
            do j2 = iabs(j1 - sps%jj(c))/2, (j1 + sps%jj(c))/2
              if(a == c .and. mod(j2, 2) == 1) cycle
              if(b == c .and. mod(j2, 2) == 1) cycle
              ich2 = two%jptz2n(j2, p2, itz2)
              bra2 = two%jptz(ich2)%labels2n(a,c)
              ket2 = two%jptz(ich2)%labels2n(b,c)
              if(bra2 * ket2 == 0) cycle
              iphase = two%jptz(ich2)%iphase(a,c) * &
                  &    two%jptz(ich2)%iphase(b,c)
              v2 = sc2%jptz(ich2)%m(bra2, ket2) * dble(iphase)
              vv = vv + dble(2*j2 + 1) * v2
            end do
            w1 = w1 + vv * Del(a,c) * Del(b,c) * sps%NoCoef(c)
          end do
          sc1%jptz(ich1)%m(bra1, ket1) = w1 / dble(j1 + 1)
          sc1%jptz(ich1)%m(ket1, bra1) = w1 / dble(j1 + 1)
        end do
      end do
    end do
  end function ScalarNOTwo2One

  function ScalarNOOne2Zero(sps, one, sc1) result(zero)
    type(spo_pn), intent(in) :: sps
    real(8) :: zero
    type(OneBodySpace), intent(in) :: one
    type(NBodyScalars), intent(in) :: sc1
    integer :: i1, n, ich1
    zero = 0.d0
    do i1 = 1, sps%n
      if(abs(sps%NoCoef(i1)) < 1.d-4) cycle
      ich1 = one%label2jptz(i1)
      n = one%jptz(ich1)%label2n(i1)
      zero = zero + sc1%jptz(ich1)%m(n,n) * sps%NoCoef(i1) * dble(sps%jj(i1) + 1)
    end do
  end function ScalarNOOne2Zero

  subroutine find_cfp(a,b,c,jab,thr, ii, idx)
    integer, intent(in) :: a, b, c, jab
    type(ThreeBodyChannel), intent(in) :: thr
    integer, intent(out) :: ii, idx
    integer :: abc(3)
    integer :: k
    integer :: i1, i2, i3, j12
    abc = (/a, b, c/)
    call sort_descending(abc)
    idx = 0
    if(.not.allocated(thr%labels2nsub)) return
    idx = thr%labels2nsub(abc(1), abc(2), abc(3))
    if(idx == 0) return
    do k = 1, thr%idx(idx)%n
      i1 = thr%idx(idx)%n2label1(k)
      i2 = thr%idx(idx)%n2label2(k)
      i3 = thr%idx(idx)%n2label3(k)
      j12= thr%idx(idx)%n2label4(k)
      if(a == i1 .and. b == i2 .and. c == i3 .and. jab == j12) then
        ii = k
        return
      end if
    end do
  end subroutine find_cfp

  subroutine sort_descending(ix)
    integer, intent(inout) :: ix(:)
    integer :: i, j, n, k
    n = size(ix)
    do i = 1, n - 1
      do j = i + 1, N
        if(ix(i) < ix(j)) then
          k = ix(i)
          ix(i) = ix(j)
          ix(j) = k
        end if
      end do
    end do
  end subroutine sort_descending

  subroutine store_idx(this, params, sps, one, j2, p2, tz2, a,b,c)
    class(braket), intent(inout) :: this
    type(spo_pn), intent(in) :: sps
    type(OneBodySpace), intent(in) :: one
    type(parameters), intent(in) :: params
    integer, intent(in) :: j2, p2, tz2, a, b, c
    integer :: j1, p1, itz1, ich1
    integer :: j3, p3, itz3, num
    integer :: loop, eho1, eho2
    integer :: e3max, e3cut

    e3max = params%e3max_3nf
    e3cut = params%e3cut
    j1 = sps%jj(c)
    p1 = (-1) ** sps%ll(c)
    itz1 = sps%itz(c)
    ich1 = one%jptz2n(j1, p1, itz1)
    do loop = 1, 2
      num = 0
      do eho1 = 1, sps%n
        if(j1 /= sps%jj(eho1)) cycle
        if(p1 /= (-1) ** sps%ll(eho1)) cycle
        if(itz1 /= sps%itz(eho1)) cycle
        if(sps%nshell(a) + sps%nshell(b) + sps%nshell(eho1) > e3max) cycle
        if(sps%nshell(a) + sps%nshell(b) + sps%nshell(eho1) > e3cut) cycle
        do eho2 = 1, sps%n
          if(j1 /= sps%jj(eho2)) cycle
          if(p1 /= (-1) ** sps%ll(eho2)) cycle
          if(itz1 /= sps%itz(eho2)) cycle
          if(sps%nshell(a) + sps%nshell(b) + sps%nshell(eho2) > e3max) cycle
          if(sps%nshell(a) + sps%nshell(b) + sps%nshell(eho2) > e3cut) cycle


          p3 = p1 * p2
          itz3 = 2*tz2 + itz1
          do j3 = iabs(2*j2 - j1), (2*j2 + j1), 2
            num = num + 1
            if(loop == 1) cycle
            this%num(1, num) = ich1
            this%num(2, num) = j3
            this%num(3, num) = p3
            this%num(4, num) = itz3
            this%num(5, num) = eho1
            this%num(6, num) = eho2
            this%num(7, num) = one%jptz(ich1)%label2n(eho1)
            this%num(8, num) = one%jptz(ich1)%label2n(eho2)
            this%num(9, num) = one%jptz(ich1)%label2n(c)
          end do
        end do
      end do
      if(loop == 1) then
        this%n = num
        allocate(this%num(9,this%n))
      end if
    end do
  end subroutine store_idx

  real(8) function Del(i1, i2)
    integer, intent(in) :: i1, i2
    Del = 1.d0
    if(i1 == i2) Del = dsqrt(2.d0)
  end function Del
end module NormalOrdering
