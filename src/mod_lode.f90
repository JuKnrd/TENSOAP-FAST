module lode

 type LODE_Model
  ! Variables for running LODE calculations
  logical nonorm,fixed_cell
  real*8 sigewald
  integer radsize,lebsize
 end type LODE_Model

 contains

!***************************************************************************************************

end module
