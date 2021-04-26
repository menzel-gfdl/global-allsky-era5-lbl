!> @brief Hu and Stamnes cloud liquid water parameterization from
!! doi: 10.1175/1520-0442(1993)006<0728:AAPOTR>2.0.CO;2
module hu_stamnes
use mo_rte_kind, only: wp

use netcdf_utils, only: close_dataset, open_dataset, read_attribute, read_variable
use optics_utils, only: OpticalProperties
implicit none
private


!> @brief Hu and Stamnes cloud liquid water parameterization.
type, public :: HuStamnes
  real(kind=wp), dimension(:,:), allocatable :: a1 !< Exctinction coefficient parameter (band, radius).
  real(kind=wp), dimension(:,:), allocatable :: a2 !< Single-scatter albedo parameter (band, radius).
  real(kind=wp), dimension(:,:), allocatable :: a3 !< Asymmetry factor parameter (band, radius).
  real(kind=wp), dimension(:,:), allocatable :: band_limits !< Lower/upper bounds of parameterization [cm-1] (2, band).
  real(kind=wp), dimension(:), allocatable :: bands !< Parameterization band centers [cm-1] (band).
  real(kind=wp), dimension(:,:), allocatable :: b1 !< Exctinction coefficient parameter (band, radius).
  real(kind=wp), dimension(:,:), allocatable :: b2 !< Single-scatter albedo parameter (band, radius).
  real(kind=wp), dimension(:,:), allocatable :: b3 !< Asymmetry factor parameter (band, radius).
  real(kind=wp), dimension(:,:), allocatable :: c1 !< Exctinction coefficient parameter (band, radius).
  real(kind=wp), dimension(:,:), allocatable :: c2 !< Single-scatter albedo parameter (band, radius).
  real(kind=wp), dimension(:,:), allocatable :: c3 !< Asymmetry factor parameter (band, radius).
  integer :: last_ir_band !< Index of last infrared (longwave) band.
  real(kind=wp) :: max_radius !< Maximum radius defined in parameterization [micron].
  real(kind=wp) :: min_radius !< Minimum radius defined in parameterization [micron].
  real(kind=wp), dimension(:), allocatable :: radii !< Radius bins [micron] for parameterization (radius).
  contains
  procedure, public :: construct
  procedure, public :: destruct
  procedure, public :: optics
  procedure, private :: optics_
endtype HuStamnes


contains


!> @brief Constructs a HuStamnes object.
subroutine construct(self, path)

  class(HuStamnes), intent(inout) :: self
  character(len=*), intent(in) :: path !< Path to input dataset.

  integer :: n
  character(len=:), allocatable :: name
  integer :: ncid
  real(kind=wp), dimension(2) :: radii
  real(kind=wp), dimension(:,:), allocatable :: radius_bounds

  ncid = open_dataset(path)
  call read_attribute(ncid, "radius", "bounds", name)
  call read_attribute(ncid, "radius", "valid_range", radii)
  self%min_radius = radii(1)
  self%max_radius = radii(2)
  call read_variable(ncid, name, radius_bounds)
  deallocate(name)
  n = size(radius_bounds, 2)
  allocate(self%radii(n + 1))
  self%radii(:n) = radius_bounds(1,:)
  self%radii(n+1) = radius_bounds(2,n)
  deallocate(radius_bounds)
  call read_variable(ncid, "band_bnds", self%band_limits)
  call read_attribute(ncid, "band_bnds", "last_IR_band", self%last_ir_band)
  call read_variable(ncid, "band", self%bands)
  call read_variable(ncid, "a1", self%a1)
  call read_variable(ncid, "a2", self%a2)
  call read_variable(ncid, "a3", self%a3)
  call read_variable(ncid, "b1", self%b1)
  call read_variable(ncid, "b2", self%b2)
  call read_variable(ncid, "b3", self%b3)
  call read_variable(ncid, "c1", self%c1)
  call read_variable(ncid, "c2", self%c2)
  call read_variable(ncid, "c3", self%c3)
  call close_dataset(ncid)
end subroutine construct


!> @brief Destructs a HuStamnes object.
subroutine destruct(self)

  class(HuStamnes), intent(inout) :: self

  if (allocated(self%a1)) deallocate(self%a1)
  if (allocated(self%a2)) deallocate(self%a2)
  if (allocated(self%a3)) deallocate(self%a3)
  if (allocated(self%band_limits)) deallocate(self%band_limits)
  if (allocated(self%bands)) deallocate(self%bands)
  if (allocated(self%b1)) deallocate(self%b1)
  if (allocated(self%b2)) deallocate(self%b2)
  if (allocated(self%b3)) deallocate(self%b3)
  if (allocated(self%c1)) deallocate(self%c1)
  if (allocated(self%c2)) deallocate(self%c2)
  if (allocated(self%c3)) deallocate(self%c3)
  if (allocated(self%radii)) deallocate(self%radii)
end subroutine destruct


!> @brief Calculates cloud optics.
elemental pure subroutine optics(self, water_concentration, equivalent_radius, &
                                 optical_properties)

  class(HuStamnes), intent(in) :: self
  real(kind=wp), intent(in) :: water_concentration !< Water concentration [g m-3].
  real(kind=wp), intent(in) :: equivalent_radius !< Particle equivalent radius [micron].
  type(OpticalProperties), intent(inout) :: optical_properties

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

  class(HuStamnes), intent(in) :: self
  real(kind=wp), intent(in) :: water_concentration !< Water concentration [g m-3].
  real(kind=wp), intent(in) :: equivalent_radius !< Droplet equivalent radius [micron].
  integer, intent(in) :: band !< Band index.
  real(kind=wp), intent(out) :: extinction_coefficient !< Extinction coefficient [m-1].
  real(kind=wp), intent(out) :: single_scatter_albedo !< Single-scatter albedo.
  real(kind=wp), intent(out) :: asymmetry_factor !< Asymmetry factor.

  real(kind=wp), parameter :: m_to_km = 1.e-3
  integer :: i
  real(kind=wp) :: r

  r = min(self%max_radius, max(self%min_radius, equivalent_radius))
  do i = 2, size(self%radii)
    if (self%radii(i) .gt. r) then
      exit
    endif
  enddo
  i = i - 1
  extinction_coefficient = water_concentration*m_to_km*(self%a1(band,i)* &
                           (r**self%b1(band,i)) + self%c1(band,i)) !Equation 13.
  single_scatter_albedo = 1. - (self%a2(band,i)*(r**self%b2(band,i)) + self%c2(band,i)) !Equation 14.
  asymmetry_factor = self%a3(band,i)*(r**self%b3(band,i)) + self%c3(band,i) !Equation 15.
end subroutine optics_


end module hu_stamnes
