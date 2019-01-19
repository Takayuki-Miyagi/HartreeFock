module MBPT
  use ModelSpace
  use Operators
  implicit none

  type :: MBPTEnergy
    real(8) :: e_0 = 0.d0
    real(8) :: e_2 = 0.d0
    real(8) :: e_3 = 0.d0
    real(8) :: e_3_pp = 0.d0
    real(8) :: e_3_hh = 0.d0
    real(8) :: e_3_ph = 0.d0
  contains
    procedure :: calc => CalcEnergyCorr
  end type MBPTEnergy
contains
  subroutine CalcEnergyCorr(this,ms,hamil)
    class(MBPTEnergy), intent(inout) :: this
    type(MSpace), intent(in) :: ms
    type(Op), intent(in) :: hamil

    this%e_0 = hamil%zero
    call energy_second()
    call energy_third()

  end subroutine CalcEnergyCorr

end module MBPT
