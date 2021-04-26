!> @brief Subcolumn generator for calculating stochastic clouds.
module stochastic_clouds
use mo_rte_kind, only: wp

use incomplete_beta, only: IncompleteBeta
implicit none
private


!> @brief Total water probability density function.
type, public :: TotalWaterPDF
  type(IncompleteBeta), pointer :: beta !< Incomplete beta distribution.
  integer :: p !< Incomplete beta distribution shape parameter.
  integer :: q !< Incomplete beta distribution shape parameter.
  contains
  procedure, public :: construct
  procedure, public :: destruct
  procedure, public :: sample_condensate
  procedure, private :: specific_saturation_humidity
  procedure, private :: width
endtype TotalWaterPDF


public :: overlap_parameter


contains


!> @brief Calculates random array for cloudiness from equation 1 from
!! doi: 10.1256/qj.03.99.
subroutine cloudiness(overlap_parameter, x)

  real(kind=wp), dimension(:), intent(in) :: overlap_parameter !< Overlap parameter between adjacent layers.
  real(kind=wp), dimension(:), intent(out) :: x !< A random number for each layer.

  integer :: i
  real(kind=wp), dimension(size(x)-1) :: r

  call random_number(x)
  call random_number(r)
  do i = 1, size(r)
    if (r(i) .le. overlap_parameter(i)) then
      x(i+1) = x(i)
    endif
  enddo
end subroutine cloudiness


!> @brief Constructs a TotalWaterPDF object.
subroutine construct(self, p, q, beta)

  class(TotalWaterPDF), intent(inout) :: self
  integer, intent(in) :: p !< Incomplete beta distribution shape parameter.
  integer, intent(in) :: q !< Incomplete beta distribution shape parameter.
  type(IncompleteBeta), intent(in), target :: beta !< IncompleteBeta object.

  self%p = p
  self%q = q
  self%beta => beta
end subroutine construct


!> @brief Destructs a TotalWaterPDF object.
subroutine destruct(self)

  class(TotalWaterPDF), intent(inout) :: self

  self%beta => null()
end subroutine destruct


!> @brief Calculates the overlap parameter defined in equation 2 from
!!        doi: 10.1029/2004JD005100.
pure subroutine overlap_parameter(altitude, scale_length, alpha)

  real(kind=wp), dimension(:), intent(in) :: altitude !< Layer altitude.
  real(kind=wp), intent(in) :: scale_length !< Scale length.
  real(kind=wp), dimension(:), intent(out) :: alpha !< Overlap parameter between adjacent layers.

  integer :: n

  n = size(altitude)
  alpha(:) = exp(-1.*abs(altitude(:n-1) - altitude(2:))/scale_length)
end subroutine overlap_parameter


!> @brief Draws samples of liquid and ice condensate mixing ratios from the total water
!!        (vapor + cloud liquid + ice) mixing ratio probability distribution function
!!        whose mean total condensate amount equals the sum of the input cloud liquid
!!        and ice condensate amounts and mean saturation humidity equals one minus
!!        the input cloud fraction for each layer.  This method is detailed in the
!!        appendix of doi: 10.1175/MWR3257.1.
subroutine sample_condensate(self, cloud_fraction, lwc, iwc, overlap, ql, qi)

  class(TotalWaterPDF), intent(in) :: self
  real(kind=wp), dimension(:), intent(in) :: cloud_fraction !< Saturated volume fraction.
  real(kind=wp), dimension(:), intent(in) :: lwc !< Cloud liquid water condensate mixing ratio in dry air [kg/kg].
  real(kind=wp), dimension(:), intent(in) :: iwc !< Cloud ice water condensate mixing ratio in dry air [kg/kg].
  real(kind=wp), dimension(:), intent(in) :: overlap !< Overlap parameter between adjacent layers.
  real(kind=wp), dimension(:), intent(out) :: ql !< Cloud liquid water condensate mixing ratio in dry air [kg/kg].
  real(kind=wp), dimension(:), intent(out) :: qi !< Cloud ice water condensate mixing ratio in dry air [kg/kg].

  real(kind=wp), dimension(size(cloud_fraction)) :: liquid_fraction
  real(kind=wp), dimension(size(cloud_fraction)) :: qs
  real(kind=wp), dimension(size(cloud_fraction)) :: total_condensate
  real(kind=wp), dimension(size(cloud_fraction)) :: width
  real(kind=wp), dimension(size(cloud_fraction)) :: x

  call cloudiness(overlap, x)
  where (x .gt. (1. - cloud_fraction))
    qs(:) = self%specific_saturation_humidity(cloud_fraction(:))
    width(:) = self%width(cloud_fraction(:), lwc(:), iwc(:), qs(:))
    total_condensate(:) = width(:)*(self%beta%inverse(self%p, self%q, x(:)) - qs(:))
    liquid_fraction = lwc(:)/(lwc(:) + iwc(:))
    ql(:) = total_condensate(:)*liquid_fraction(:)
    qi(:) = total_condensate(:)*(1. - liquid_fraction(:))
  elsewhere
    ql(:) = 0.
    qi(:) = 0.
  endwhere
end subroutine sample_condensate


!> @brief Calculates normalized saturation specific humidity from equation A1 from
!!        doi: 10.1175/MWR3257.1$ 
elemental pure function specific_saturation_humidity(self, cloud_fraction) &
  result(qs)

  class(TotalWaterPDF), intent(in) :: self
  real(kind=wp), intent(in) :: cloud_fraction !< saturated volume fraction.
  real(kind=wp) :: qs !< Normalized saturation specific humidity in dry air [kg/kg].

  qs = self%beta%inverse(self%p, self%q, 1. - cloud_fraction)
end function specific_saturation_humidity


!> @brief Calculates the width of the total water probability distribution function (b - a)
!!        from equation A2 from doi: 10.1175/MWR3257.1, ignoring the parameter alpha.
elemental pure function width(self, cloud_fraction, lwc, iwc, qs) &
  result(w)

  class(TotalWaterPDF), intent(in) :: self
  real(kind=wp), intent(in) :: cloud_fraction !< Saturated volume fraction.
  real(kind=wp), intent(in) :: lwc !< Cloud liquid water condensate mixing ratio in dry air [kg/kg].
  real(kind=wp), intent(in) :: iwc !< Cloud ice water condensate mixing ratio in dry air [kg/kg].
  real(kind=wp), intent(in) :: qs !< Normalized saturation specific humidity in dry air [kg/kg].
  real(kind=wp) :: w !< Distribution width (b - a).

  w = (lwc + iwc)/((real(self%p)/real(self%p + self%q))* &
      (1. - self%beta%value(self%p + 1, self%q, qs)) - qs*cloud_fraction)
end function width


end module stochastic_clouds
