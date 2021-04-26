!> @brief Slingo liquid water cloud optics parameterization from
!! doi: 10.1175/1520-0469(1989)046<1419:AGPFTS>2.0.CO;2.
module slingo_liquid_cloud
use mo_rte_kind, only: wp

use netcdf_utils, only: close_dataset, open_dataset, read_attribute, read_variable
use optics_utils, only: OpticalProperties
implicit none
private


!> @brief Slingo liquid water cloud optics parameterization.
type, public :: Slingo
  real(kind=wp), dimension(:), allocatable :: a !< Optical depth coefficients [g-1 m2] (band).
  real(kind=wp), dimension(:), allocatable :: b !< Optical depth coefficients [g-1 m2 micron] (band).
  real(kind=wp), dimension(:,:), allocatable :: band_limits !< Parameterization band limits [cm-1] (2, band).
  real(kind=wp), dimension(:), allocatable :: bands !< Parameterization band centers [cm-1].
  real(kind=wp), dimension(:), allocatable :: c !< Single-scatter albedo coefficients (band).
  real(kind=wp), dimension(:), allocatable :: d !< Single-scatter albedo coefficients [micron-1] (band).
  real(kind=wp), dimension(:), allocatable :: e !< Asymmetry factor coefficients (band).
  real(kind=wp), dimension(:), allocatable :: f !< Asymmetry factor coefficients [micron-1] (band).
  real(kind=wp) :: max_radius !< Maximum radius defined in parameterization [micron].
  real(kind=wp) :: min_radius !< Minimum radius defined in parameterization [micron].
  contains
  procedure, public :: construct
  procedure, public :: destruct
  procedure, public :: optics
  procedure, private :: optics_
endtype Slingo


contains


!> @brief Constructs a Slingo object.
subroutine construct(self, path)

  class(Slingo), intent(inout) :: self
  character(len=*), intent(in) :: path !< Path to input dataset.

  integer :: ncid

  ncid = open_dataset(path)
  call read_attribute(ncid, "max_radius", self%max_radius)
  call read_attribute(ncid, "min_radius", self%min_radius)
  call read_variable(ncid, "band_bnds", self%band_limits)
  allocate(self%bands(size(self%band_limits, 2)))
  self%bands(:) = 0.5*(self%band_limits(1,:) + self%band_limits(2,:))
  call read_variable(ncid, "a", self%a)
  call read_variable(ncid, "b", self%b)
  call read_variable(ncid, "c", self%c)
  call read_variable(ncid, "d", self%d)
  call read_variable(ncid, "e", self%e)
  call read_variable(ncid, "f", self%f)
  call close_dataset(ncid)
end subroutine construct


!> @brief Destructs a Slingo object.
subroutine destruct(self)

  class(Slingo), intent(inout) :: self

  if (allocated(self%a)) deallocate(self%a)
  if (allocated(self%band_limits)) deallocate(self%band_limits)
  if (allocated(self%bands)) deallocate(self%bands)
  if (allocated(self%b)) deallocate(self%b)
  if (allocated(self%c)) deallocate(self%c)
  if (allocated(self%d)) deallocate(self%d)
  if (allocated(self%e)) deallocate(self%e)
  if (allocated(self%f)) deallocate(self%f)
end subroutine destruct


!> @brief Calculates cloud optics.
elemental pure subroutine optics(self, water_concentration, equivalent_radius, &
                                 optical_properties)

  class(Slingo), intent(in) :: self
  real(kind=wp), intent(in) :: water_concentration !< Water concentration [g m-3].
  real(kind=wp), intent(in) :: equivalent_radius !< Particle equivalent radius [micron].
  type(OpticalProperties), intent(inout) :: optical_properties !< Optical properties.

  integer :: i

  do i = 1, size(optical_properties%bands)
    call self%optics_(water_concentration, equivalent_radius, i, &
                      optical_properties%extinction_coefficient(i), &
                      optical_properties%single_scatter_albedo(i), &
                      optical_properties%asymmetry_factor(i))
  enddo
end subroutine optics


!> @brief Calculates cloud optics.
elemental pure subroutine optics_(self, water_concentration, equivalent_radius, band, &
                                  extinction_coefficient, single_scatter_albedo, &
                                  asymmetry_factor)

  class(Slingo), intent(in) :: self
  real(kind=wp), intent(in) :: water_concentration !< Water concentration [g m-3].
  real(kind=wp), intent(in) :: equivalent_radius !< Droplet equivalent radius [micron].
  integer, intent(in) :: band !< Band index.
  real(kind=wp), intent(out) :: extinction_coefficient !< Extinction coefficient.
  real(kind=wp), intent(out) :: single_scatter_albedo !< Single-scatter albedo.
  real(kind=wp), intent(out) :: asymmetry_factor !< Asymmetry factor.

  real(kind=wp) :: r

  r = min(self%max_radius, max(self%min_radius, equivalent_radius))
  extinction_coefficient = water_concentration*(self%a(band) + self%b(band)/r) !Equation 1.
  single_scatter_albedo = 1. - (self%c(band) + self%d(band)*r) !Equation 2.
  asymmetry_factor = self%e(band) + self%f(band)*r !Equation 3.
end subroutine optics_


end module slingo_liquid_cloud
