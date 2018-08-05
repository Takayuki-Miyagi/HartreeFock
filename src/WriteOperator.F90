module WriteOperator
  use class_sys, only: sy
  use class_stopwatch
  use MPIFunction, only: myrank
  use InputParameters, only: parameters
  use ModelSpace, only: spo_pn, MSpace, OneBodySpace, TwoBodySpace, ThreeBodySpace, &
      & OneBodyChannel, TwoBodyChannel, ThreeBodyChannel
  use read_3bme, only: iThreeBodyScalar
  use ScalarOperator
  implicit none
  type(sy) :: sys
  private :: write_snt_bin, write_snt_txt
  public :: write_hamil
contains
  subroutine write_hamil(params, sps, ms, hamil)
    type(parameters), intent(in) :: params
    type(spo_pn), intent(in) :: sps
    type(MSpace), intent(in) :: ms
    type(ScalarOperators), intent(in) :: hamil
    if(sys%find(params%no2bhfile, '.bin')) call write_snt_bin(params, sps, ms, hamil)
    if(sys%find(params%no2bhfile, '.txt')) call write_snt_txt(params, sps, ms, hamil)
  end subroutine write_hamil

  subroutine write_snt_bin(params, sps, ms, hamil)
    type(parameters), intent(in) :: params
    type(spo_pn), intent(in) :: sps
    type(MSpace), intent(in) :: ms
    type(ScalarOperators), intent(in) :: hamil
    integer :: iunit = 19
    integer :: i, num
    integer :: ich
    integer :: n
    integer :: bra, ket
    integer, allocatable :: aa(:), bb(:), cc(:), dd(:), jj(:)
    real(8), allocatable :: v1save(:)
    integer, allocatable :: nnum(:)

    call start_stopwatch(time_io_write)
    open(iunit, form = "unformatted", file = params%no2bhfile, status = "replace", access='stream')
    write(iunit) sps%n/2, sps%n/2, 0, 0
    allocate(nnum(sps%n))
    do i = 1, sps%n
      nnum(i) = i
    end do
    write(iunit) nnum
    write(iunit) sps%nn
    write(iunit) sps%ll
    write(iunit) sps%jj
    write(iunit) sps%itz
    deallocate(nnum)

    write(iunit) hamil%zero
    num = 0
    do ich = 1, ms%one%n
      num = num + ms%one%jptz(ich)%n * (ms%one%jptz(ich)%n + 1) / 2
    end do
    write(iunit) num, 10, params%hw
    allocate(aa(num), bb(num), v1save(num))
    num = 0
    do ich = 1, ms%one%n
      n = ms%one%jptz(ich)%n
      do bra = 1, n
        do ket = 1, bra
          num = num + 1
          aa(num) = ms%one%jptz(ich)%n2label(bra)
          bb(num) = ms%one%jptz(ich)%n2label(ket)
          v1save(num) = hamil%one%jptz(ich)%m(bra, ket)
        end do
      end do
    end do
    write(iunit) aa
    write(iunit) bb
    write(iunit) v1save
    deallocate(aa, bb, v1save)

    num = 0
    do ich = 1, ms%two%n
      n = ms%two%jptz(ich)%n
      num = num + n * (n + 1) / 2
    end do
    write(iunit) num
    allocate(aa(num), bb(num), cc(num), dd(num), jj(num))
    allocate(v1save(num))
    num = 0
    do ich = 1, ms%two%n
      n = ms%two%jptz(ich)%n
      do bra = 1, n
        do ket = 1, bra
          num = num + 1
          aa(num) = ms%two%jptz(ich)%n2label1(bra)
          bb(num) = ms%two%jptz(ich)%n2label2(bra)
          cc(num) = ms%two%jptz(ich)%n2label1(ket)
          dd(num) = ms%two%jptz(ich)%n2label2(ket)
          jj(num) = ms%two%j(ich)
          v1save(num) = hamil%two%jptz(ich)%m(bra, ket)
        end do
      end do
    end do
    write(iunit) aa
    write(iunit) bb
    write(iunit) cc
    write(iunit) dd
    write(iunit) jj
    write(iunit) v1save
    close(iunit)
    call stop_stopwatch(time_io_write)
    deallocate(aa, bb, cc, dd, jj)
    deallocate(v1save)
  end subroutine write_snt_bin

  subroutine write_snt_txt(params, sps, ms, hamil)
    type(parameters), intent(in) :: params
    type(spo_pn), intent(in) :: sps
    type(MSpace), intent(in) :: ms
    type(ScalarOperators), intent(in) :: hamil
    integer :: iunit = 19
    integer :: i, num
    integer :: ich
    integer :: n
    integer :: bra, ket
    integer, allocatable :: nnum(:)

    call start_stopwatch(time_io_write)
    open(iunit, file = params%no2bhfile, status = "replace")
    call params%PrtParams(iunit)
    write(iunit, '(a)') '! n_proton, n_neutron, core_proton, core_neutron'
    write(iunit, '(a)') '! num, n, l, j, tz'
    write(iunit, '(4i4)') sps%n/2, sps%n/2, 0, 0
    allocate(nnum(sps%n))
    do i = 1, sps%n
      nnum(i) = i
    end do
    do i = 1, sps%n
      write(iunit, '(5i4)') nnum(i), sps%nn(i), sps%ll(i), &
          & sps%jj(i), sps%itz(i)
    end do
    deallocate(nnum)

    write(iunit, '(a)') '! Zero-body term (MeV)'
    write(iunit, '(f15.8)') hamil%zero
    num = 0
    do ich = 1, ms%one%n
      num = num + ms%one%jptz(ich)%n * (ms%one%jptz(ich)%n + 1) / 2
    end do
    write(iunit, '(a)') '! One-body term'
    write(iunit, '(a)') '! num, method1 = 0, hw'
    write(iunit, '(a)') '! i, j, <i|f|j>'
    write(iunit, '(2i5, f8.4)') num, 0, params%hw
    do ich = 1, ms%one%n
      n = ms%one%jptz(ich)%n
      do bra = 1, n
        do ket = 1, n
          write(iunit, '(2i4, f15.8)') ms%one%jptz(ich)%n2label(bra), &
              & ms%one%jptz(ich)%n2label(ket), hamil%one%jptz(ich)%m(bra,ket)
        end do
      end do
    end do

    num = 0
    do ich = 1, ms%two%n
      n = ms%two%jptz(ich)%n
      num = num + n * (n + 1) / 2
    end do
    write(iunit, '(a)') '! Two-body term'
    write(iunit, '(a)') '! num, method2 = 0, hw'
    write(iunit, '(a)') '! i, j, k, l, JJ, <ij:JJ|v|kl:JJ>'
    write(iunit, '(1i10, i5, f8.4)') num, 0, params%hw
    num = 0
    do ich = 1, ms%two%n
      n = ms%two%jptz(ich)%n
      do bra = 1, n
        do ket = 1, bra
          num = num + 1
          write(iunit, '(5i4, f15.8)') ms%two%jptz(ich)%n2label1(bra), &
              & ms%two%jptz(ich)%n2label2(bra), &
              & ms%two%jptz(ich)%n2label1(ket), &
              & ms%two%jptz(ich)%n2label2(ket), &
              & ms%two%j(ich), &
              & hamil%two%jptz(ich)%m(bra,ket)
        end do
      end do
    end do
    close(iunit)
    call stop_stopwatch(time_io_write)

  end subroutine write_snt_txt
end module WriteOperator
