module HartreeFock
  use LinAlgLib
  use Operators
  implicit none

  private :: InitMonopoleVec2
  private :: InitMonopoleVec3
  private :: FinMonopoleVec
  private :: InitMonopole2
  private :: InitMonopole3
  private :: FinMonopole
  private :: GetIndex2
  private :: GetSpLabels2
  private :: GetIndex3
  private :: GetSpLabels3

  type :: HFSolver
    logical :: is_three_body
    type(DMat), allocatable :: rho(:)
  end type HFSolver

  type, private :: MonopoleVec
    integer :: n_dim = 0
    integer :: J, P, Z
    real(8), allocatable :: v(:)
    integer(8), allocatable :: idx(:)
  contains
    procedure :: InitMonopoleVec2
    procedure :: InitMonopoleVec3
    procedure :: FinMonopoleVec
  end type MonopoleVec

  type, private :: Monopole
    type(MonopoleVec), allocatable :: jpz(:)
    integer, allocatable :: jpz2ch(:,:,:)
    integer :: NChan
  contains
    procedure :: InitMonopole2
    procedure :: InitMonopole3
    procedure :: FinMonopole
  end type Monopole

contains

  subroutine FinMonopole(this)
    class(Monopole), intent(inout) :: this
    integer :: ch
    do ch = 1, this%NChan
      call this%jpz(ch)%FinMonopoleVec()
    end do
    deallocate(this%jpz)
    deallocate(this%jpz2ch)
  end subroutine FinMonopole

  subroutine InitMonopole2(this, ms, vnn)
    class(Monopole), intent(inout) :: this
    type(MSpace), intent(in) :: ms
    type(NBodyPart), intent(in) :: vnn
    integer :: p, j, z, ch, n
    integer :: i1, l1, j1, z1, e1
    integer :: i2, l2, j2, z2, e2
    integer :: i3, l3, j3, z3, e3
    integer :: i4, l4, j4, z4, e4
    integer, allocatable :: jj(:), pp(:), zz(:), nn(:)

    allocate(this%jpz2ch(1:2*min(ms%lmax, ms%emax)+1,-1:1,-1:1))
    this%jpz2ch(:,:,:) = 0

    ch = 0
    do j = 1, 2*min(ms%lmax, ms%emax)+1, 2
      do p = 1, -1, -2
        do z = -1, 1, 2

          n = 0
          do i1 = 1, ms%sps%norbs
            l1 = ms%sps%orb(i1)%l
            j1 = ms%sps%orb(i1)%j
            z1 = ms%sps%orb(i1)%z
            e1 = ms%sps%orb(i1)%e
            if(j1 /= j) cycle
            if((-1)**l1 /= p) cycle
            if(z1 /= z) cycle

            do i3 = 1, ms%sps%norbs
              l3 = ms%sps%orb(i3)%l
              j3 = ms%sps%orb(i3)%j
              z3 = ms%sps%orb(i3)%z
              e3 = ms%sps%orb(i3)%e
              if(j3 /= j) cycle
              if((-1)**l3 /= p) cycle
              if(z3 /= z) cycle

              do i2 = 1, ms%sps%norbs
                l2 = ms%sps%orb(i2)%l
                j2 = ms%sps%orb(i2)%j
                z2 = ms%sps%orb(i2)%z
                e2 = ms%sps%orb(i2)%e

                do i4 = 1, ms%sps%norbs
                  l4 = ms%sps%orb(i4)%l
                  j4 = ms%sps%orb(i4)%j
                  z4 = ms%sps%orb(i4)%z
                  e4 = ms%sps%orb(i4)%e
                  if(j2 /= j4) cycle
                  if((-1)**l2 /= (-1)**l4) cycle
                  if(z2 /= z4) cycle

                  if(e1 + e2 > ms%e2max) cycle
                  if(e3 + e4 > ms%e2max) cycle

                  n = n + 1

                end do
              end do
            end do
          end do
          if(n /= 0) then
            ch = ch + 1
            this%jpz2ch(J,P,Z) = ch
          end if

        end do
      end do
    end do

    this%NChan = ch
    allocate(this%jpz(ch))
    allocate(jj(ch))
    allocate(pp(ch))
    allocate(zz(ch))
    allocate(nn(ch))

    ch = 0
    do j = 1, 2*min(ms%lmax, ms%emax)+1, 2
      do p = 1, -1, -2
        do z = -1, 1, 2

          n = 0
          do i1 = 1, ms%sps%norbs
            l1 = ms%sps%orb(i1)%l
            j1 = ms%sps%orb(i1)%j
            z1 = ms%sps%orb(i1)%z
            e1 = ms%sps%orb(i1)%e
            if(j1 /= j) cycle
            if((-1)**l1 /= p) cycle
            if(z1 /= z) cycle

            do i3 = 1, ms%sps%norbs
              l3 = ms%sps%orb(i3)%l
              j3 = ms%sps%orb(i3)%j
              z3 = ms%sps%orb(i3)%z
              e3 = ms%sps%orb(i3)%e
              if(j3 /= j) cycle
              if((-1)**l3 /= p) cycle
              if(z3 /= z) cycle

              do i2 = 1, ms%sps%norbs
                l2 = ms%sps%orb(i2)%l
                j2 = ms%sps%orb(i2)%j
                z2 = ms%sps%orb(i2)%z
                e2 = ms%sps%orb(i2)%e

                do i4 = 1, ms%sps%norbs
                  l4 = ms%sps%orb(i4)%l
                  j4 = ms%sps%orb(i4)%j
                  z4 = ms%sps%orb(i4)%z
                  e4 = ms%sps%orb(i4)%e
                  if(j2 /= j4) cycle
                  if((-1)**l2 /= (-1)**l4) cycle
                  if(z2 /= z4) cycle

                  if(e1 + e2 > ms%e2max) cycle
                  if(e3 + e4 > ms%e2max) cycle

                  n = n + 1

                end do
              end do
            end do
          end do
          if(n /= 0) then
            ch = ch + 1
            jj(ch) = j
            pp(ch) = p
            zz(ch) = z
            nn(ch) = n
          end if

        end do
      end do
    end do

    do ch = 1, this%NChan
      call this%jpz(ch)%InitMonopoleVec2(ms, vnn, jj(ch),pp(ch),zz(ch),nn(ch))
    end do

    deallocate(jj)
    deallocate(pp)
    deallocate(zz)
    deallocate(nn)

  end subroutine InitMonopole2

  subroutine InitMonopole3(this, ms, v3n)
    use CommonLibrary, only: triag
    class(Monopole), intent(inout) :: this
    type(MSpace), intent(in) :: ms
    type(NBodyPartSp), intent(in) :: v3n
    integer :: p, j, z, ch, n
    integer :: i1, l1, j1, z1, e1
    integer :: i2, l2, j2, z2, e2
    integer :: i3, l3, j3, z3, e3
    integer :: i4, l4, j4, z4, e4
    integer :: i5, l5, j5, z5, e5
    integer :: i6, l6, j6, z6, e6
    integer, allocatable :: jj(:), pp(:), zz(:), nn(:)

    allocate(this%jpz2ch(1:2*min(ms%lmax, ms%emax)+1,-1:1,-1:1))
    this%jpz2ch(:,:,:) = 0

    ch = 0
    do j = 1, 2*min(ms%lmax, ms%emax)+1, 2
      do p = 1, -1, -2
        do z = -1, 1, 2

          n = 0
          do i1 = 1, ms%sps%norbs
            l1 = ms%sps%orb(i1)%l
            j1 = ms%sps%orb(i1)%j
            z1 = ms%sps%orb(i1)%z
            e1 = ms%sps%orb(i1)%e

            do i4 = 1, ms%sps%norbs
              l4 = ms%sps%orb(i4)%l
              j4 = ms%sps%orb(i4)%j
              z4 = ms%sps%orb(i4)%z
              e4 = ms%sps%orb(i4)%e

              if(j4 /= j1) cycle
              if((-1)**l4 /= (-1)**l1) cycle
              if(z4 /= z1) cycle

              if(j4 /= j) cycle
              if((-1)**l4 /= P) cycle
              if(z4 /= z) cycle

              do i2 = 1, ms%sps%norbs
                l2 = ms%sps%orb(i2)%l
                j2 = ms%sps%orb(i2)%j
                z2 = ms%sps%orb(i2)%z
                e2 = ms%sps%orb(i2)%e

                do i5 = 1, ms%sps%norbs
                  l5 = ms%sps%orb(i5)%l
                  j5 = ms%sps%orb(i5)%j
                  z5 = ms%sps%orb(i5)%z
                  e5 = ms%sps%orb(i5)%e
                  if(j2 /= j5) cycle
                  if((-1)**l2 /= (-1)**l5) cycle
                  if(z2 /= z5) cycle

                  if(e1 + e2 > ms%e2max) cycle
                  if(e4 + e5 > ms%e2max) cycle

                  do i3 = 1, ms%sps%norbs
                    l3 = ms%sps%orb(i3)%l
                    j3 = ms%sps%orb(i3)%j
                    z3 = ms%sps%orb(i3)%z
                    e3 = ms%sps%orb(i3)%e
                    do i6 = 1, ms%sps%norbs
                      l6 = ms%sps%orb(i6)%l
                      j6 = ms%sps%orb(i6)%j
                      z6 = ms%sps%orb(i6)%z
                      e6 = ms%sps%orb(i6)%e
                      if(j3 /= j6) cycle
                      if((-1)**l3 /= (-1)**l6) cycle
                      if(z3 /= z6) cycle

                      if(e2+e3 > ms%e2max) cycle
                      if(e5+e6 > ms%e2max) cycle
                      if(e1+e3 > ms%e2max) cycle
                      if(e4+e6 > ms%e2max) cycle
                      if(e1+e2+e3 > ms%e3max) cycle
                      if(e4+e5+e6 > ms%e3max) cycle

                      n = n + 1
                    end do
                  end do

                end do
              end do
            end do
          end do
          if(n /= 0) then
            ch = ch + 1
            this%jpz2ch(J,P,Z) = ch
          end if

        end do
      end do
    end do

    this%NChan = ch
    allocate(this%jpz(ch))
    allocate(jj(ch))
    allocate(pp(ch))
    allocate(zz(ch))
    allocate(nn(ch))

    ch = 0
    do j = 1, 2*min(ms%lmax, ms%emax)+1, 2
      do p = 1, -1, -2
        do z = -1, 1, 2

          n = 0
          do i1 = 1, ms%sps%norbs
            l1 = ms%sps%orb(i1)%l
            j1 = ms%sps%orb(i1)%j
            z1 = ms%sps%orb(i1)%z
            e1 = ms%sps%orb(i1)%e

            do i4 = 1, ms%sps%norbs
              l4 = ms%sps%orb(i4)%l
              j4 = ms%sps%orb(i4)%j
              z4 = ms%sps%orb(i4)%z
              e4 = ms%sps%orb(i4)%e

              if(j4 /= j1) cycle
              if((-1)**l4 /= (-1)**l1) cycle
              if(z4 /= z1) cycle

              if(j4 /= j) cycle
              if((-1)**l4 /= P) cycle
              if(z4 /= z) cycle

              do i2 = 1, ms%sps%norbs
                l2 = ms%sps%orb(i2)%l
                j2 = ms%sps%orb(i2)%j
                z2 = ms%sps%orb(i2)%z
                e2 = ms%sps%orb(i2)%e

                do i5 = 1, ms%sps%norbs
                  l5 = ms%sps%orb(i5)%l
                  j5 = ms%sps%orb(i5)%j
                  z5 = ms%sps%orb(i5)%z
                  e5 = ms%sps%orb(i5)%e
                  if(j2 /= j5) cycle
                  if((-1)**l2 /= (-1)**l5) cycle
                  if(z2 /= z5) cycle

                  if(e1 + e2 > ms%e2max) cycle
                  if(e4 + e5 > ms%e2max) cycle

                  do i3 = 1, ms%sps%norbs
                    l3 = ms%sps%orb(i3)%l
                    j3 = ms%sps%orb(i3)%j
                    z3 = ms%sps%orb(i3)%z
                    e3 = ms%sps%orb(i3)%e
                    do i6 = 1, ms%sps%norbs
                      l6 = ms%sps%orb(i6)%l
                      j6 = ms%sps%orb(i6)%j
                      z6 = ms%sps%orb(i6)%z
                      e6 = ms%sps%orb(i6)%e
                      if(j3 /= j6) cycle
                      if((-1)**l3 /= (-1)**l6) cycle
                      if(z3 /= z6) cycle

                      if(e2+e3 > ms%e2max) cycle
                      if(e5+e6 > ms%e2max) cycle
                      if(e1+e3 > ms%e2max) cycle
                      if(e4+e6 > ms%e2max) cycle
                      if(e1+e2+e3 > ms%e3max) cycle
                      if(e4+e5+e6 > ms%e3max) cycle

                      n = n + 1
                    end do
                  end do

                end do
              end do
            end do
          end do
          if(n /= 0) then
            ch = ch + 1
            jj(ch) = j
            pp(ch) = p
            zz(ch) = z
            nn(ch) = n
          end if

        end do
      end do
    end do

    do ch = 1, this%NChan
      call this%jpz(ch)%InitMonopoleVec3(ms,v3n,jj(ch),pp(ch),zz(ch),nn(ch))
    end do

    deallocate(jj)
    deallocate(pp)
    deallocate(zz)
    deallocate(nn)

  end subroutine InitMonopole3

  subroutine FinMonopoleVec(this)
    class(MonopoleVec), intent(inout) :: this
    deallocate(this%v)
    deallocate(this%idx)
  end subroutine FinMonopoleVec

  subroutine InitMonopoleVec2(this, ms, vnn, J, P, Z, n)
    class(MonopoleVec), intent(inout) :: this
    type(MSpace), intent(in) :: ms
    type(NBodyPart), intent(in) :: vnn
    integer, intent(in) :: J, P, Z, N
    integer :: cnt, idx
    integer(8) :: i1, i2, i3, i4, num
    integer :: l1, j1, z1, e1
    integer :: l2, j2, z2, e2
    integer :: l3, j3, z3, e3
    integer :: l4, j4, z4, e4
    integer :: JJ
    real(8) :: norm, v
    this%n_dim = N
    this%J = J
    this%P = P
    this%Z = Z

    allocate(this%v(N))
    allocate(this%idx(N))

    cnt = 0
    do i1 = 1, ms%sps%norbs
      l1 = ms%sps%orb(i1)%l
      j1 = ms%sps%orb(i1)%j
      z1 = ms%sps%orb(i1)%z
      e1 = ms%sps%orb(i1)%e
      if(j1 /= j) cycle
      if((-1)**l1 /= p) cycle
      if(z1 /= z) cycle

      do i3 = 1, ms%sps%norbs
        l3 = ms%sps%orb(i3)%l
        j3 = ms%sps%orb(i3)%j
        z3 = ms%sps%orb(i3)%z
        e3 = ms%sps%orb(i3)%e
        if(j3 /= j) cycle
        if((-1)**l3 /= p) cycle
        if(z3 /= z) cycle

        do i2 = 1, ms%sps%norbs
          l2 = ms%sps%orb(i2)%l
          j2 = ms%sps%orb(i2)%j
          z2 = ms%sps%orb(i2)%z
          e2 = ms%sps%orb(i2)%e

          do i4 = 1, ms%sps%norbs
            l4 = ms%sps%orb(i4)%l
            j4 = ms%sps%orb(i4)%j
            z4 = ms%sps%orb(i4)%z
            e4 = ms%sps%orb(i4)%e
            if(j2 /= j4) cycle
            if((-1)**l2 /= (-1)**l4) cycle
            if(z2 /= z4) cycle

            if(e1 + e2 > ms%e2max) cycle
            if(e3 + e4 > ms%e2max) cycle

            cnt = cnt + 1
            this%idx(cnt) = GetIndex2(i1,i2,i3,i4)

          end do
        end do
      end do
    end do

    !$omp parallel
    !$omp do private(idx,num,i1,i2,i3,i4,j1,j2,v,JJ)
    do idx = 1, this%n_dim
      num = this%idx(cnt)
      call GetSpLabels2(num,i1,i2,i3,i4)
      j1 = ms%sps%orb(i1)%j
      j2 = ms%sps%orb(i2)%j
      v = 0.d0
      norm = 1.d0
      if(i1 == i2) norm = norm * dsqrt(2.d0)
      if(i3 == i4) norm = norm * dsqrt(2.d0)
      do JJ = abs(j1-j2)/2, (j1+j2)/2
        v = v + dble(2*JJ+1) * &
            & vnn%GetTwBME(ms%sps,ms%two,i1,i2,i3,i4,JJ)
      end do
      this%v(idx) = v * norm
    end do
    !$omp end do
    !$omp end parallel
  end subroutine InitMonopoleVec2

  subroutine InitMonopoleVec3(this, ms, v3n, J, P, Z, n)
    use CommonLibrary, only: triag
    class(MonopoleVec), intent(inout) :: this
    type(MSpace), intent(in) :: ms
    type(NBodyPartSp), intent(in) :: v3n
    integer, intent(in) :: J, P, Z, N
    integer :: cnt, idx
    integer(8) :: i1, i2, i3, i4, i5, i6, num
    integer :: l1, j1, z1, e1
    integer :: l2, j2, z2, e2
    integer :: l3, j3, z3, e3
    integer :: l4, j4, z4, e4
    integer :: l5, j5, z5, e5
    integer :: l6, j6, z6, e6
    integer :: JJ, JJJ
    real(8) :: norm, v
    this%n_dim = N
    this%J = J
    this%P = P
    this%Z = Z

    allocate(this%v(N))
    allocate(this%idx(N))

    cnt = 0
    do i1 = 1, ms%sps%norbs
      l1 = ms%sps%orb(i1)%l
      j1 = ms%sps%orb(i1)%j
      z1 = ms%sps%orb(i1)%z
      e1 = ms%sps%orb(i1)%e

      do i4 = 1, ms%sps%norbs
        l4 = ms%sps%orb(i4)%l
        j4 = ms%sps%orb(i4)%j
        z4 = ms%sps%orb(i4)%z
        e4 = ms%sps%orb(i4)%e
        if(j4 /= j1) cycle
        if((-1)**l4 /= (-1)**l1) cycle
        if(z4 /= z1) cycle

        if(j4 /= j) cycle
        if((-1)**l4 /= p) cycle
        if(z4 /= z) cycle

        do i2 = 1, ms%sps%norbs
          l2 = ms%sps%orb(i2)%l
          j2 = ms%sps%orb(i2)%j
          z2 = ms%sps%orb(i2)%z
          e2 = ms%sps%orb(i2)%e

          do i5 = 1, ms%sps%norbs
            l5 = ms%sps%orb(i5)%l
            j5 = ms%sps%orb(i5)%j
            z5 = ms%sps%orb(i5)%z
            e5 = ms%sps%orb(i5)%e
            if(j2 /= j5) cycle
            if((-1)**l2 /= (-1)**l5) cycle
            if(z2 /= z5) cycle

            if(e1 + e2 > ms%e2max) cycle
            if(e3 + e4 > ms%e2max) cycle

            do i3 = 1, ms%sps%norbs
              l3 = ms%sps%orb(i3)%l
              j3 = ms%sps%orb(i3)%j
              z3 = ms%sps%orb(i3)%z
              e3 = ms%sps%orb(i3)%e
              do i6 = 1, ms%sps%norbs
                l6 = ms%sps%orb(i6)%l
                j6 = ms%sps%orb(i6)%j
                z6 = ms%sps%orb(i6)%z
                e6 = ms%sps%orb(i6)%e
                if(j3 /= j6) cycle
                if((-1)**l3 /= (-1)**l6) cycle
                if(z3 /= z6) cycle

                if(e1+e2+e3 > ms%e3max) cycle
                if(e4+e5+e6 > ms%e3max) cycle

                cnt = cnt + 1
                this%idx(cnt) = GetIndex3(i1,i2,i3,i4,i5,i6)

              end do
            end do
          end do

        end do
      end do
    end do

    !$omp parallel
    !$omp do private(idx,num,i1,i2,i3,i4,i5,i6,j1,j2,j3,v,JJ,JJJ)
    do idx = 1, this%n_dim
      num = this%idx(cnt)
      call GetSpLabels3(num,i1,i2,i3,i4,i5,i6)
      j1 = ms%sps%orb(i1)%j
      j2 = ms%sps%orb(i2)%j
      j3 = ms%sps%orb(i3)%j
      v = 0.d0
      do JJ = abs(j1-j2)/2, (j1+j2)/2
        do JJJ = abs(2*JJ-j3), (2*JJ+j3), 2
          v = v + dble(JJJ+1) * v3n%GetThBME(ms,i1,i2,i3,JJ,i4,i5,i6,JJ,JJJ) * norm
        end do
      end do
      this%v(idx) = v / dble(j1+1)
    end do
    !$omp end do
    !$omp end parallel
  end subroutine InitMonopoleVec3

  function GetIndex2(i1,i2,i3,i4) result(r)
    integer(8), intent(in) :: i1, i2, i3, i4
    integer(8) :: r
    r = i1 + rshift(i2,10) + rshift(i3,20) + rshift(i4,30)
  end function GetIndex2

  subroutine GetSpLabels2(idx,i1,i2,i3,i4)
    integer(8), intent(in)  :: idx
    integer(8), intent(out) :: i1, i2, i3, i4
    i1 = mod(idx,1024)
    i2 = mod(rshift(idx,10),1024)
    i3 = mod(rshift(idx,20),1024)
    i4 = mod(rshift(idx,30),1024)
  end subroutine GetSpLabels2

  function GetIndex3(i1,i2,i3,i4,i5,i6) result(r)
    integer(8), intent(in) :: i1, i2, i3, i4, i5, i6
    integer(8) :: r
    r = i1 + rshift(i2,10) + rshift(i3,20) + rshift(i4,30) + &
        &    rshift(i5,40) + rshift(i6,50)
  end function GetIndex3

  subroutine GetSpLabels3(idx,i1,i2,i3,i4,i5,i6)
    integer(8), intent(in)  :: idx
    integer(8), intent(out) :: i1, i2, i3, i4, i5, i6
    i1 = mod(idx,1024)
    i2 = mod(rshift(idx,10),1024)
    i3 = mod(rshift(idx,20),1024)
    i4 = mod(rshift(idx,30),1024)
    i5 = mod(rshift(idx,40),1024)
    i6 = mod(rshift(idx,50),1024)
  end subroutine GetSpLabels3
end module HartreeFock
