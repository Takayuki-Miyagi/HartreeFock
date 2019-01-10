! Define one-, two-, and three-body model space
! Note that the three-body space is orthonormalized
! with coefficient of fractional parentage (cfp), which will not be used
! in Hartree-Fock calculation. However, this definition may be useful for
! beyond Hartree-Fock calculation methods.
module ModelSpace
  use SingleParticleState
  implicit none

  type :: jpz
    integer :: j = -1
    integer :: p = 0
    integer :: z = 100
  end type jpz

  type, extends(jpz) :: OneBodyChannel
    integer :: n
    integer, allocatable :: n2spi(:)
    integer, allocatable :: spi2n(:)
  end type OneBodyChannel

  type :: OneBodySpace
    type(OneBodyChannel), allocatable :: jpz(:)
    integer, allocatable :: jpz2ch(:,:,:)
  end type OneBodySpace

  type, extends(jpz) :: TwoBodyChannel
    integer :: n
    integer, allocatable :: n2spi1(:)
    integer, allocatable :: n2spi2(:)
    integer, allocatable :: spis2n(:,:)
    integer, allocatable :: iphase(:,:)
  end type TwoBodyChannel

  type :: TwoBodySpace
    type(TwoBodyChannel), allocatable :: jpz(:)
    integer, allocatable :: jpz2ch(:,:,:)
  end type TwoBodySpace

  type :: AdditionalQN
    integer :: north, nphys
    integer, allocatable :: Jab2n(:)
    integer, allocatable :: n2spi1(:)
    integer, allocatable :: n2spi2(:)
    integer, allocatable :: n2spi3(:)
    integer, allocatable :: n2Jab(:)
    real(8), allocatable :: cfp(:,:)
  end type AdditionalQN

  type, extends(jpz) :: ThreeBodyChannel
    integer :: n, n_idx
    type(AdditionalQN), allocatable :: idx(:)
    integer, allocatable :: n2spi1(:)
    integer, allocatable :: n2spi2(:)
    integer, allocatable :: n2spi3(:)
    integer, allocatable :: n2labl(:)
    integer, allocatable :: spis2idx(:,:,:)
  end type ThreeBodyChannel

  type :: ThreeBodySpace
    type(TwoBodyChannel), allocatable :: jpz(:)
    integer, allocatable :: jpz2ch(:,:,:)
  end type ThreeBodySpace

  type :: MSpace
    type(Orbits) :: sps
    type(OneBodySpace) :: one
    type(TwoBodySpace) :: two
    type(ThreeBodySpace) :: thr
    logical :: is_constructed=.false.
    real(8), allocatable :: NOCoef(:)
    integer, allocatable :: hole_orbts(:)
    integer, allocatable :: particle_orbts(:)
  end type MSpace

  type, private :: iter3
    integer, allocatable :: ii1(:), ii2(:), ii3(:)
  contains
    procedure :: GenIter
  end type iter3
contains
  subroutine GenIter(ite, i1, i2, i3, num)
    class(iter3), intent(inout) :: ite
    integer, intent(in) :: i1, i2, i3
    integer, intent(out) :: num
    integer :: icase
    if(allocated(ite%ii1)) deallocate(ite%ii1)
    if(allocated(ite%ii2)) deallocate(ite%ii2)
    if(allocated(ite%ii3)) deallocate(ite%ii3)

    icase = -100
    if(i2 == i3 .and. i1 == i2) icase = 1
    if(i2 == i3 .and. i1 /= i2) icase = 2
    if(i2 /= i3 .and. i1 == i2) icase = 3
    if(i2 /= i3 .and. i1 == i3) icase = 4
    if(i2 /= i3 .and. i1 /= i2 .and. i1 /= i3) icase = 5

    select case(icase)
    case(1)
      num = 1
      allocate(ite%ii1(num), ite%ii2(num), ite%ii3(num))
      ite%ii1(1) = i1; ite%ii2(1) = i1; ite%ii3(1) = i1
    case(2)
      num = 3
      allocate(ite%ii1(num), ite%ii2(num), ite%ii3(num))
      ite%ii1(1) = i1; ite%ii2(1) = i2; ite%ii3(1) = i2
      ite%ii1(2) = i2; ite%ii2(2) = i1; ite%ii3(2) = i2
      ite%ii1(3) = i2; ite%ii2(3) = i2; ite%ii3(3) = i1
    case(3)
      num = 3
      allocate(ite%ii1(num), ite%ii2(num), ite%ii3(num))
      ite%ii1(1) = i1; ite%ii2(1) = i1; ite%ii3(1) = i3
      ite%ii1(2) = i3; ite%ii2(2) = i1; ite%ii3(2) = i1
      ite%ii1(3) = i1; ite%ii2(3) = i3; ite%ii3(3) = i1
    case(4)
      num = 3
      allocate(ite%ii1(num), ite%ii2(num), ite%ii3(num))
      ite%ii1(1) = i1; ite%ii2(1) = i2; ite%ii3(1) = i1
      ite%ii1(2) = i1; ite%ii2(2) = i1; ite%ii3(2) = i2
      ite%ii1(3) = i2; ite%ii2(3) = i1; ite%ii3(3) = i1
    case(5)
      num = 6
      allocate(ite%ii1(num), ite%ii2(num), ite%ii3(num))
      ite%ii1(1) = i1; ite%ii2(1) = i2; ite%ii3(1) = i3
      ite%ii1(2) = i3; ite%ii2(2) = i1; ite%ii3(2) = i2
      ite%ii1(3) = i2; ite%ii2(3) = i3; ite%ii3(3) = i1
      ite%ii1(4) = i2; ite%ii2(4) = i1; ite%ii3(4) = i3
      ite%ii1(5) = i3; ite%ii2(5) = i2; ite%ii3(5) = i1
      ite%ii1(6) = i1; ite%ii2(6) = i3; ite%ii3(6) = i2
    case default
      write(*,*) 'Error in GenIter:'
      return
    end select
  end subroutine GenIter
end module ModelSpace
