module optics_utils
use mo_rte_kind, only: wp
implicit none
private


!> @brief Optical properties.
type, public :: OpticalProperties
  real(kind=wp), dimension(:), allocatable :: bands !< Band center wavenumbers [cm-1] (band).
  real(kind=wp), dimension(:,:), allocatable :: band_limits !< Band edge wavenumbers [cm-1] (2, band).
  real(kind=wp), dimension(:), allocatable :: extinction_coefficient !< Extinction coefficient (band).
  real(kind=wp), dimension(:), allocatable :: single_scatter_albedo !< Single-scatter albedo (band).
  real(kind=wp), dimension(:), allocatable :: asymmetry_factor !< Asymmetry factor (band).
  contains
  procedure, public :: construct
  procedure, public :: destruct
  procedure, public :: thick_average
endtype OpticalProperties


contains


!> @brief Constructs an OpticalProperties object.
pure subroutine construct(self, bands, band_limits)

  class(OpticalProperties), intent(inout) :: self
  real(kind=wp), dimension(:), intent(in) :: bands !< Band center wavenumbers [cm-1]
  real(kind=wp), dimension(2,size(bands)), intent(in) :: band_limits !< Band edge wavenumbers [cm-1].

  integer :: n

  n = size(bands)
  allocate(self%bands(n))
  self%bands(:) = bands(:)
  allocate(self%band_limits(2, n))
  self%band_limits(:,:) = band_limits(:,:)
  allocate(self%extinction_coefficient(n))
  allocate(self%single_scatter_albedo(n))
  allocate(self%asymmetry_factor(n))
end subroutine construct


!> @brief Destructs an OpticalProperties object.
subroutine destruct(self)

  class(OpticalProperties), intent(inout) :: self

  if (allocated(self%bands)) deallocate(self%bands)
  if (allocated(self%band_limits)) deallocate(self%band_limits)
  if (allocated(self%extinction_coefficient)) deallocate(self%extinction_coefficient)
  if (allocated(self%single_scatter_albedo)) deallocate(self%single_scatter_albedo)
  if (allocated(self%asymmetry_factor)) deallocate(self%asymmetry_factor)
end subroutine destruct


!> @brief Linearly interpolates.
pure function interp(x, y, newx) result(newy)

  real(kind=wp), dimension(2), intent(in) :: x !< Domain values.
  real(kind=wp), dimension(2), intent(in) :: y !< Function values
  real(kind=wp), intent(in) :: newx !< New domain value.
  real(kind=wp) :: newy !< Interpolated function value.

  real(kind=wp) :: m
  real(kind=wp) :: b

  m = (y(2) - y(1))/(x(2) - x(1))
  b = y(1) - m*x(1)
  newy = m*newx + b
end function interp


!> @brief Creates a new OpticalProperties object on a new set of bands from the current
!!        OpticalProperties object.
elemental pure subroutine thick_average(self, optics, starting_band, ending_band)

  class(OpticalProperties), intent(in) :: self
  type(OpticalProperties), intent(inout) :: optics !< Optical properties on input bands.
  integer, intent(in), optional :: starting_band !< Index of lowest band to interpolate in.
  integer, intent(in), optional :: ending_band !< Index of highest band to interpolate in.

  integer :: a
  integer :: b
  integer :: c
  integer :: i
  integer :: j
  integer :: k
  integer :: n

  if (present(starting_band)) then
    a = starting_band
  else
    a = 1
  endif
  if (present(ending_band)) then
    b = ending_band
  else
    b = size(self%bands)
  endif

  !Linearly interpolate for now.
  n = size(optics%bands)
  do i = 1, n
    if (self%band_limits(1,a) .lt. optics%bands(i)) exit
  enddo
  if (i .gt. 1) then
    !Use value from the first band in all new bands that are less than the first band's
    !lower limit.
    optics%extinction_coefficient(:i-1) = self%extinction_coefficient(a)
    optics%single_scatter_albedo(:i-1) = self%single_scatter_albedo(a)
    optics%asymmetry_factor(:i-1) = self%asymmetry_factor(a)
  endif
  c = a
  do j = i, n
    if (optics%bands(j) .gt. self%band_limits(2,b)) exit
    do k = c, b
      if (self%bands(k) .gt. optics%bands(j)) exit
    enddo
    c = k
    if (k .eq. a) then
      if (optics%bands(j) .le. self%bands(k)) then
        !Inside the first band, but left of the band center.
        optics%extinction_coefficient(j) = self%extinction_coefficient(k)
        optics%single_scatter_albedo(j) = self%single_scatter_albedo(k)
        optics%asymmetry_factor(j) = self%asymmetry_factor(k)
        cycle
      endif
      k = a + 1
    endif
    if (k .eq. b + 1) then
      k = b
      if (optics%bands(j) .ge. self%bands(k)) then
        !Inside the last band, but right of the band center.
        optics%extinction_coefficient(j) = self%extinction_coefficient(k)
        optics%single_scatter_albedo(j) = self%single_scatter_albedo(k)
        optics%asymmetry_factor(j) = self%asymmetry_factor(k)
        cycle
      endif
    endif
    optics%extinction_coefficient(j) = interp(self%bands(k-1:k), &
                                              self%extinction_coefficient(k-1:k), &
                                              optics%bands(j))
    optics%single_scatter_albedo(j) = interp(self%bands(k-1:k), &
                                             self%single_scatter_albedo(k-1:k), &
                                             optics%bands(j))
    optics%asymmetry_factor(j) = interp(self%bands(k-1:k), &
                                        self%asymmetry_factor(k-1:k), &
                                        optics%bands(j))
  enddo
  if (j .le. n) then
    !Use value from the last band in all new bands that are greater than the last band's
    !upper limit.
    optics%extinction_coefficient(j:) = self%extinction_coefficient(b)
    optics%single_scatter_albedo(j:) = self%single_scatter_albedo(b)
    optics%asymmetry_factor(j:) = self%asymmetry_factor(b)
  endif
end subroutine thick_average


end module optics_utils
