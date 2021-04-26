!> @brief NetCDF utilities
module netcdf_utils
use, intrinsic :: iso_fortran_env, only: error_unit
use netcdf

use mo_rte_kind, only: wl, wp
implicit none
private


public :: add_dimension
public :: add_variable
public :: close_dataset
public :: create_dataset
public :: dimension_length
public :: open_dataset
public :: read_attribute
public :: read_variable
public :: variable_exists
public :: write_variable


interface read_attribute
  module procedure read_attribute_int
  module procedure read_attribute_string
  module procedure read_attribute_1d_real
  module procedure read_global_attribute_real
end interface read_attribute


interface read_variable
  module procedure read_variable_1d_string
  module procedure read_variable_0d_int
  module procedure read_variable_1d_int
  module procedure read_variable_2d_int
  module procedure read_variable_3d_int
  module procedure read_variable_0d_real
  module procedure read_variable_1d_real
  module procedure read_variable_2d_real
  module procedure read_variable_3d_real
  module procedure read_variable_4d_real
  module procedure read_variable_1d_int_to_bool
end interface read_variable


interface write_variable
  module procedure write_variable_1d_real
  module procedure write_variable_2d_real
  module procedure write_variable_3d_real
  module procedure write_variable_4d_real
end interface write_variable


contains


!> @brief Add a dimension to a netCDF dataset.
function add_dimension(ncid, name, length) result(dimid)

  integer, intent(in) :: ncid !< NetCDF id.
  character(len=*), intent(in) :: name
  integer, intent(in), optional :: length
  integer :: dimid

  integer :: err
  integer :: s

  if (present(length)) then
    s = length
  else
    s = nf90_unlimited
  endif
  err = nf90_def_dim(ncid, trim(name), s, dimid)
  call netcdf_catch(err)
end function add_dimension


!> @brief Add a variable to a netCDF dataset.
function add_variable(ncid, dimid, name, standard_name, units, fill_value, &
                      positive, axis, calendar) result(varid)

  integer, intent(in) :: ncid !< NetCDF id.
  integer, dimension(:), intent(in) :: dimid !< Dimension ids.
  character(len=*), intent(in) :: name !< Variable name.
  character(len=*), intent(in) :: standard_name !< Variable standard name.
  character(len=*), intent(in), optional :: units !< Variable units.
  real(kind=wp), intent(in), optional :: fill_value !< Fill value.
  character(len=*), intent(in), optional :: positive !< Vertical sense.
  character(len=1), intent(in), optional :: axis !< Axis id.
  character(len=*), intent(in), optional :: calendar !< Calendar.
  integer :: varid

  integer :: err

  err = nf90_def_var(ncid, trim(name), nf90_double, dimid, varid)
  call netcdf_catch(err)
  err = nf90_put_att(ncid, varid, "standard_name", trim(standard_name))
  call netcdf_catch(err)
  if (present(units)) then
    err = nf90_put_att(ncid, varid, "units", trim(units))
    call netcdf_catch(err)
  endif
  if (present(fill_value)) then
    err = nf90_put_att(ncid, varid, "_FillValue", fill_value)
    call netcdf_catch(err)
  endif
  if (present(positive)) then
    err = nf90_put_att(ncid, varid, "positive", trim(positive))
    call netcdf_catch(err)
  endif
  if (present(axis)) then
    err = nf90_put_att(ncid, varid, "axis", trim(axis))
    call netcdf_catch(err)
  endif
  if (present(calendar)) then
    err = nf90_put_att(ncid, varid, "calendar", trim(calendar))
    call netcdf_catch(err)
  endif
end function add_variable


!> @brief Closes netCDF dataset.
subroutine close_dataset(ncid)

  integer, intent(in) :: ncid !< NetCDF id.

  integer :: err

  err = nf90_close(ncid)
  call netcdf_catch(err)
end subroutine close_dataset


!> @brief Creates netCDF dataset.
function create_dataset(path) result(ncid)

  character(len=*), intent(in) :: path !< File path.
  integer :: ncid !< NetCDF id.

  integer :: err

  err = nf90_create(trim(path), nf90_netcdf4, ncid)
  call netcdf_catch(err)
end function create_dataset


!> @brief Gets the length of a dimension.
function dimension_length(ncid, name) result(length)

  integer, intent(in) :: ncid
  character(len=*), intent(in) :: name
  integer :: length

  integer :: dimid
  integer :: err

  err = nf90_inq_dimid(ncid, trim(name), dimid)
  call netcdf_catch(err)
  err = nf90_inquire_dimension(ncid, dimid, len=length)
  call netcdf_catch(err)
  return
end function dimension_length


!> @brief Crashes if any netCDF errors are detected.
subroutine netcdf_catch(err)

  integer, intent(in) :: err !< Code returned from netCDF functions.

  if (err .ne. nf90_noerr) then
    write(error_unit, *) nf90_strerror(err)
    stop 1
  endif
end subroutine netcdf_catch


!> @brief Opens netCDF dataset.
function open_dataset(path) result(ncid)

  character(len=*), intent(in) :: path !< Path to dataset.
  integer :: ncid !< NetCDF id.

  integer :: err

  err = nf90_open(trim(path), nf90_nowrite, ncid)
  call netcdf_catch(err)
end function open_dataset


!> @brief Reads an attribute from netCDF dataset.
subroutine read_attribute_int(ncid, variable, name, buffer)

  integer, intent(in) :: ncid !< NetCDF id.
  character(len=*), intent(in) :: variable !< Variable name.
  character(len=*), intent(in) :: name !< Attribute name.
  integer, intent(out) :: buffer

  integer :: err
  integer :: varid

  err = nf90_inq_varid(ncid, trim(variable), varid)
  call netcdf_catch(err)
  err = nf90_get_att(ncid, varid, trim(name), buffer)
  call netcdf_catch(err)
end subroutine read_attribute_int


!> @brief Reads an attribute from netCDF dataset.
subroutine read_attribute_string(ncid, variable, name, buffer)

  integer, intent(in) :: ncid !< NetCDF id.
  character(len=*), intent(in) :: variable !< Variable name.
  character(len=*), intent(in) :: name !< Attribute name.
  character(len=:), allocatable, intent(inout) :: buffer

  integer :: err
  integer :: n
  integer :: varid

  err = nf90_inq_varid(ncid, trim(variable), varid)
  call netcdf_catch(err)
  err = nf90_inquire_attribute(ncid, varid, trim(name), len=n)
  call netcdf_catch(err)
  if (allocated(buffer)) deallocate(buffer)
  allocate(character(len=n) :: buffer)
  err = nf90_get_att(ncid, varid, trim(name), buffer)
  call netcdf_catch(err)
end subroutine read_attribute_string


!> @brief Reads an attribute from netCDF dataset.
subroutine read_attribute_1d_real(ncid, variable, name, buffer)

  integer, intent(in) :: ncid !< NetCDF id.
  character(len=*), intent(in) :: variable !< Variable name.
  character(len=*), intent(in) :: name !< Attribute name.
  real(kind=wp), dimension(:), intent(out) :: buffer

  integer :: err
  integer :: varid

  err = nf90_inq_varid(ncid, trim(variable), varid)
  call netcdf_catch(err)
  err = nf90_get_att(ncid, varid, trim(name), buffer)
  call netcdf_catch(err)
end subroutine read_attribute_1d_real


!> @brief Reads an attribute from netCDF dataset.
subroutine read_global_attribute_real(ncid, name, buffer)

  integer, intent(in) :: ncid !< NetCDF id.
  character(len=*), intent(in) :: name !< Attribute name.
  real(kind=wp), intent(out) :: buffer

  integer :: err

  err = nf90_get_att(ncid, nf90_global, trim(name), buffer)
  call netcdf_catch(err)
end subroutine read_global_attribute_real


!> @brief Reads variable from netCDF dataset.
subroutine read_variable_1d_string(ncid, name, buffer, start, counts)

  integer, intent(in) :: ncid !< NetCDF id.
  character(len=*), intent(in) :: name !< Variable name.
  character(len=*), dimension(:), allocatable, intent(inout) :: buffer
  integer, dimension(:), intent(in), optional :: start !< Corner indices
  integer, dimension(:), intent(in), optional :: counts !< Edge lengths.

  integer, dimension(nf90_max_dims) :: dimids
  integer :: err
  integer :: i
  integer :: n
  integer, dimension(nf90_max_dims) :: sizes
  integer, dimension(nf90_max_dims) :: starts
  integer :: varid

  err = nf90_inq_varid(ncid, trim(name), varid)
  call netcdf_catch(err)
  err = nf90_inquire_variable(ncid, varid, ndims=n)
  call netcdf_catch(err)
  if (present(counts)) then
    sizes(:n) = counts(:n)
  else
    err = nf90_inquire_variable(ncid, varid, dimids=dimids)
    call netcdf_catch(err)
    do i = 1, n
      err = nf90_inquire_dimension(ncid, dimids(i), len=sizes(i))
      call netcdf_catch(err)
    enddo
  endif
  if (allocated(buffer)) deallocate(buffer)
  allocate(buffer(sizes(2)))
  err = nf90_get_var(ncid, varid, buffer, start, sizes)
  call netcdf_catch(err)
end subroutine read_variable_1d_string


!> @brief Reads variable from netCDF dataset.
subroutine read_variable_0d_int(ncid, name, buffer, start)

  integer, intent(in) :: ncid !< NetCDF id.
  character(len=*), intent(in) :: name !< Variable name.
  integer, intent(inout) :: buffer
  integer, dimension(:), intent(in), optional :: start !< Corner indices

  integer :: add
  integer :: err
  integer :: n
  integer :: scales
  integer :: varid

  err = nf90_inq_varid(ncid, trim(name), varid)
  call netcdf_catch(err)
  err = nf90_get_var(ncid, varid, buffer, start)
  call netcdf_catch(err)
  err = nf90_get_att(ncid, varid, "scale_factor", scales)
  if (err .eq. nf90_enotatt) then
    scales = 1
  else
    call netcdf_catch(err)
  endif
  err = nf90_get_att(ncid, varid, "add_offset", add)
  if (err .eq. nf90_enotatt) then
    add = 0
  else
    call netcdf_catch(err)
  endif
  buffer = buffer*scales + add
end subroutine read_variable_0d_int


!> @brief Reads variable from netCDF dataset.
subroutine read_variable_1d_int(ncid, name, buffer, start, counts)

  integer, intent(in) :: ncid !< NetCDF id.
  character(len=*), intent(in) :: name !< Variable name.
  integer, dimension(:), allocatable, intent(inout) :: buffer
  integer, dimension(:), intent(in), optional :: start !< Corner indices
  integer, dimension(:), intent(in), optional :: counts !< Edge lengths.

  integer :: add
  integer, dimension(nf90_max_dims) :: dimids
  integer :: err
  integer :: i
  integer :: n
  integer :: scales
  integer, dimension(nf90_max_dims) :: sizes
  integer :: varid

  err = nf90_inq_varid(ncid, trim(name), varid)
  call netcdf_catch(err)
  err = nf90_inquire_variable(ncid, varid, ndims=n)
  call netcdf_catch(err)
  if (present(counts)) then
    sizes(:n) = counts(:n)
  else
    err = nf90_inquire_variable(ncid, varid, dimids=dimids)
    call netcdf_catch(err)
    do i = 1, n
      err = nf90_inquire_dimension(ncid, dimids(i), len=sizes(i))
      call netcdf_catch(err)
    enddo
  endif
  if (allocated(buffer)) deallocate(buffer)
  allocate(buffer(sizes(1)))
  err = nf90_get_var(ncid, varid, buffer, start, sizes)
  call netcdf_catch(err)
  err = nf90_get_att(ncid, varid, "scale_factor", scales)
  if (err .eq. nf90_enotatt) then
    scales = 1
  else
    call netcdf_catch(err)
  endif
  err = nf90_get_att(ncid, varid, "add_offset", add)
  if (err .eq. nf90_enotatt) then
    add = 0
  else
    call netcdf_catch(err)
  endif
  buffer(:) = buffer(:)*scales + add
end subroutine read_variable_1d_int


!> @brief Reads variable from netCDF dataset.
subroutine read_variable_2d_int(ncid, name, buffer, start, counts)

  integer, intent(in) :: ncid !< NetCDF id.
  character(len=*), intent(in) :: name !< Variable name.
  integer, dimension(:,:), allocatable, intent(inout) :: buffer
  integer, dimension(:), intent(in), optional :: start !< Corner indices
  integer, dimension(:), intent(in), optional :: counts !< Edge lengths.

  integer :: add
  integer, dimension(nf90_max_dims) :: dimids
  integer :: err
  integer :: i
  integer :: n
  integer :: scales
  integer, dimension(nf90_max_dims) :: sizes
  integer :: varid

  err = nf90_inq_varid(ncid, trim(name), varid)
  call netcdf_catch(err)
  err = nf90_inquire_variable(ncid, varid, ndims=n)
  call netcdf_catch(err)
  if (present(counts)) then
    sizes(:n) = counts(:n)
  else
    err = nf90_inquire_variable(ncid, varid, dimids=dimids)
    call netcdf_catch(err)
    do i = 1, n
      err = nf90_inquire_dimension(ncid, dimids(i), len=sizes(i))
      call netcdf_catch(err)
    enddo
  endif
  if (allocated(buffer)) deallocate(buffer)
  allocate(buffer(sizes(1), sizes(2)))
  err = nf90_get_var(ncid, varid, buffer, start, sizes)
  call netcdf_catch(err)
  err = nf90_get_att(ncid, varid, "scale_factor", scales)
  if (err .eq. nf90_enotatt) then
    scales = 1
  else
    call netcdf_catch(err)
  endif
  err = nf90_get_att(ncid, varid, "add_offset", add)
  if (err .eq. nf90_enotatt) then
    add = 0
  else
    call netcdf_catch(err)
  endif
  buffer(:,:) = buffer(:,:)*scales + add
end subroutine read_variable_2d_int


!> @brief Reads variable from netCDF dataset.
subroutine read_variable_3d_int(ncid, name, buffer, start, counts)

  integer, intent(in) :: ncid !< NetCDF id.
  character(len=*), intent(in) :: name !< Variable name.
  integer, dimension(:,:,:), allocatable, intent(inout) :: buffer
  integer, dimension(:), intent(in), optional :: start !< Corner indices
  integer, dimension(:), intent(in), optional :: counts !< Edge lengths.

  integer :: add
  integer, dimension(nf90_max_dims) :: dimids
  integer :: err
  integer :: i
  integer :: n
  integer :: scales
  integer, dimension(nf90_max_dims) :: sizes
  integer :: varid

  err = nf90_inq_varid(ncid, trim(name), varid)
  call netcdf_catch(err)
  err = nf90_inquire_variable(ncid, varid, ndims=n)
  call netcdf_catch(err)
  if (present(counts)) then
    sizes(:n) = counts(:n)
  else
    err = nf90_inquire_variable(ncid, varid, dimids=dimids)
    call netcdf_catch(err)
    do i = 1, n
      err = nf90_inquire_dimension(ncid, dimids(i), len=sizes(i))
      call netcdf_catch(err)
    enddo
  endif
  if (allocated(buffer)) deallocate(buffer)
  allocate(buffer(sizes(1), sizes(2), sizes(3)))
  err = nf90_get_var(ncid, varid, buffer, start, sizes)
  call netcdf_catch(err)
  err = nf90_get_att(ncid, varid, "scale_factor", scales)
  if (err .eq. nf90_enotatt) then
    scales = 1
  else
    call netcdf_catch(err)
  endif
  err = nf90_get_att(ncid, varid, "add_offset", add)
  if (err .eq. nf90_enotatt) then
    add = 0
  else
    call netcdf_catch(err)
  endif
  buffer(:,:,:) = buffer(:,:,:)*scales + add
end subroutine read_variable_3d_int


!> @brief Reads variable from netCDF dataset.
subroutine read_variable_0d_real(ncid, name, buffer, start)

  integer, intent(in) :: ncid !< NetCDF id.
  character(len=*), intent(in) :: name !< Variable name.
  real(kind=wp), intent(inout) :: buffer
  integer, dimension(:), intent(in), optional :: start !< Corner indices

  real(kind=wp) :: add
  integer :: err
  real(kind=wp) :: scales
  integer :: varid

  err = nf90_inq_varid(ncid, trim(name), varid)
  call netcdf_catch(err)
  err = nf90_get_var(ncid, varid, buffer, start)
  call netcdf_catch(err)
  err = nf90_get_att(ncid, varid, "scale_factor", scales)
  if (err .eq. nf90_enotatt) then
    scales = 1._wp
  else
    call netcdf_catch(err)
  endif
  err = nf90_get_att(ncid, varid, "add_offset", add)
  if (err .eq. nf90_enotatt) then
    add = 0._wp
  else
    call netcdf_catch(err)
  endif
  buffer = buffer*scales + add
end subroutine read_variable_0d_real


!> @brief Reads variable from netCDF dataset.
subroutine read_variable_1d_real(ncid, name, buffer, start, counts)

  integer, intent(in) :: ncid !< NetCDF id.
  character(len=*), intent(in) :: name !< Variable name.
  real(kind=wp), dimension(:), allocatable, intent(inout) :: buffer
  integer, dimension(:), intent(in), optional :: start !< Corner indices
  integer, dimension(:), intent(in), optional :: counts !< Edge lengths.

  real(kind=wp) :: add
  integer, dimension(nf90_max_dims) :: dimids
  integer :: err
  integer :: i
  integer :: n
  real(kind=wp) :: scales
  integer, dimension(nf90_max_dims) :: sizes
  integer :: varid

  err = nf90_inq_varid(ncid, trim(name), varid)
  call netcdf_catch(err)
  err = nf90_inquire_variable(ncid, varid, ndims=n)
  call netcdf_catch(err)
  if (present(counts)) then
    sizes(:n) = counts(:n)
  else
    err = nf90_inquire_variable(ncid, varid, dimids=dimids)
    call netcdf_catch(err)
    do i = 1, n
      err = nf90_inquire_dimension(ncid, dimids(i), len=sizes(i))
      call netcdf_catch(err)
    enddo
  endif
  if (allocated(buffer)) deallocate(buffer)
  allocate(buffer(sizes(1)))
  err = nf90_get_var(ncid, varid, buffer, start, sizes)
  call netcdf_catch(err)
  err = nf90_get_att(ncid, varid, "scale_factor", scales)
  if (err .eq. nf90_enotatt) then
    scales = 1._wp
  else
    call netcdf_catch(err)
  endif
  err = nf90_get_att(ncid, varid, "add_offset", add)
  if (err .eq. nf90_enotatt) then
    add = 0._wp
  else
    call netcdf_catch(err)
  endif
  buffer(:) = buffer(:)*scales + add
end subroutine read_variable_1d_real


!> @brief Reads variable from netCDF dataset.
subroutine read_variable_2d_real(ncid, name, buffer, start, counts)

  integer, intent(in) :: ncid !< NetCDF id.
  character(len=*), intent(in) :: name !< Variable name.
  real(kind=wp), dimension(:,:), allocatable, intent(inout) :: buffer
  integer, dimension(:), intent(in), optional :: start !< Corner indices
  integer, dimension(:), intent(in), optional :: counts !< Edge lengths.

  real(kind=wp) :: add
  integer, dimension(nf90_max_dims) :: dimids
  integer :: err
  integer :: i
  integer :: n
  real(kind=wp) :: scales
  integer, dimension(nf90_max_dims) :: sizes
  integer :: varid

  err = nf90_inq_varid(ncid, trim(name), varid)
  call netcdf_catch(err)
  err = nf90_inquire_variable(ncid, varid, ndims=n)
  call netcdf_catch(err)
  if (present(counts)) then
    sizes(:n) = counts(:n)
  else
    err = nf90_inquire_variable(ncid, varid, dimids=dimids)
    call netcdf_catch(err)
    do i = 1, n
      err = nf90_inquire_dimension(ncid, dimids(i), len=sizes(i))
      call netcdf_catch(err)
    enddo
  endif
  if (allocated(buffer)) deallocate(buffer)
  allocate(buffer(sizes(1), sizes(2)))
  err = nf90_get_var(ncid, varid, buffer, start, sizes)
  call netcdf_catch(err)
  err = nf90_get_att(ncid, varid, "scale_factor", scales)
  if (err .eq. nf90_enotatt) then
    scales = 1._wp
  else
    call netcdf_catch(err)
  endif
  err = nf90_get_att(ncid, varid, "add_offset", add)
  if (err .eq. nf90_enotatt) then
    add = 0._wp
  else
    call netcdf_catch(err)
  endif
  buffer(:,:) = buffer(:,:)*scales + add
end subroutine read_variable_2d_real


!> @brief Reads variable from netCDF dataset.
subroutine read_variable_3d_real(ncid, name, buffer, start, counts)

  integer, intent(in) :: ncid !< NetCDF id.
  character(len=*), intent(in) :: name !< Variable name.
  real(kind=wp), dimension(:,:,:), allocatable, intent(inout) :: buffer
  integer, dimension(:), intent(in), optional :: start !< Corner indices
  integer, dimension(:), intent(in), optional :: counts !< Edge lengths.

  real(kind=wp) :: add
  integer, dimension(nf90_max_dims) :: dimids
  integer :: err
  integer :: i
  integer :: n
  real(kind=wp) :: scales
  integer, dimension(nf90_max_dims) :: sizes
  integer :: varid

  err = nf90_inq_varid(ncid, trim(name), varid)
  call netcdf_catch(err)
  err = nf90_inquire_variable(ncid, varid, ndims=n)
  call netcdf_catch(err)
  if (present(counts)) then
    sizes(:n) = counts(:n)
  else
    err = nf90_inquire_variable(ncid, varid, dimids=dimids)
    call netcdf_catch(err)
    do i = 1, n
      err = nf90_inquire_dimension(ncid, dimids(i), len=sizes(i))
      call netcdf_catch(err)
    enddo
  endif
  if (allocated(buffer)) deallocate(buffer)
  allocate(buffer(sizes(1), sizes(2), sizes(3)))
  err = nf90_get_var(ncid, varid, buffer, start, sizes)
  call netcdf_catch(err)
  err = nf90_get_att(ncid, varid, "scale_factor", scales)
  if (err .eq. nf90_enotatt) then
    scales = 1._wp
  else
    call netcdf_catch(err)
  endif
  err = nf90_get_att(ncid, varid, "add_offset", add)
  if (err .eq. nf90_enotatt) then
    add = 0._wp
  else
    call netcdf_catch(err)
  endif
  buffer(:,:,:) = buffer(:,:,:)*scales + add
end subroutine read_variable_3d_real


!> @brief Reads variable from netCDF dataset.
subroutine read_variable_4d_real(ncid, name, buffer, start, counts)

  integer, intent(in) :: ncid !< NetCDF id.
  character(len=*), intent(in) :: name !< Variable name.
  real(kind=wp), dimension(:,:,:,:), allocatable, intent(inout) :: buffer
  integer, dimension(:), intent(in), optional :: start !< Corner indices
  integer, dimension(:), intent(in), optional :: counts !< Edge lengths.

  real(kind=wp) :: add
  integer, dimension(nf90_max_dims) :: dimids
  integer :: err
  integer :: i
  integer :: n
  real(kind=wp) :: scales
  integer, dimension(nf90_max_dims) :: sizes
  integer :: varid

  err = nf90_inq_varid(ncid, trim(name), varid)
  call netcdf_catch(err)
  err = nf90_inquire_variable(ncid, varid, ndims=n)
  call netcdf_catch(err)
  if (present(counts)) then
    sizes(:n) = counts(:n)
  else
    err = nf90_inquire_variable(ncid, varid, dimids=dimids)
    call netcdf_catch(err)
    do i = 1, n
      err = nf90_inquire_dimension(ncid, dimids(i), len=sizes(i))
      call netcdf_catch(err)
    enddo
  endif
  if (allocated(buffer)) deallocate(buffer)
  allocate(buffer(sizes(1), sizes(2), sizes(3), sizes(4)))
  err = nf90_get_var(ncid, varid, buffer, start, sizes)
  call netcdf_catch(err)
  err = nf90_get_att(ncid, varid, "scale_factor", scales)
  if (err .eq. nf90_enotatt) then
    scales = 1._wp
  else
    call netcdf_catch(err)
  endif
  err = nf90_get_att(ncid, varid, "add_offset", add)
  if (err .eq. nf90_enotatt) then
    add = 0._wp
  else
    call netcdf_catch(err)
  endif
  buffer(:,:,:,:) = buffer(:,:,:,:)*scales + add
end subroutine read_variable_4d_real


!> @brief Writes variable to netCDF dataset.
subroutine write_variable_1d_real(ncid, varid, buffer, start, counts)

  integer, intent(in) :: ncid !< NetCDF id.
  integer, intent(in) :: varid !< Variable id.
  real(kind=wp), dimension(:), intent(in) :: buffer
  integer, dimension(:), intent(in), optional :: start !< Corner indices
  integer, dimension(:), intent(in), optional :: counts !< Edge lengths.

  integer :: err

  err = nf90_put_var(ncid, varid, buffer, start, counts)
  call netcdf_catch(err)
end subroutine write_variable_1d_real


!> @brief Writes variable to netCDF dataset.
subroutine write_variable_2d_real(ncid, varid, buffer, start, counts)

  integer, intent(in) :: ncid !< NetCDF id.
  integer, intent(in) :: varid !< Variable id.
  real(kind=wp), dimension(:,:), intent(in) :: buffer
  integer, dimension(:), intent(in), optional :: start !< Corner indices
  integer, dimension(:), intent(in), optional :: counts !< Edge lengths.

  integer :: err

  err = nf90_put_var(ncid, varid, buffer, start, counts)
  call netcdf_catch(err)
end subroutine write_variable_2d_real


!> @brief Writes variable to netCDF dataset.
subroutine write_variable_3d_real(ncid, varid, buffer, start, counts)

  integer, intent(in) :: ncid !< NetCDF id.
  integer, intent(in) :: varid !< Variable id.
  real(kind=wp), dimension(:,:,:), intent(in) :: buffer
  integer, dimension(:), intent(in), optional :: start !< Corner indices
  integer, dimension(:), intent(in), optional :: counts !< Edge lengths.

  integer :: err

  err = nf90_put_var(ncid, varid, buffer, start, counts)
  call netcdf_catch(err)
end subroutine write_variable_3d_real


!> @brief Writes variable to netCDF dataset.
subroutine write_variable_4d_real(ncid, varid, buffer, start, counts)

  integer, intent(in) :: ncid !< NetCDF id.
  integer, intent(in) :: varid !< Variable id.
  real(kind=wp), dimension(:,:,:,:), intent(in) :: buffer
  integer, dimension(:), intent(in), optional :: start !< Corner indices
  integer, dimension(:), intent(in), optional :: counts !< Edge lengths.

  integer :: err

  err = nf90_put_var(ncid, varid, buffer, start, counts)
  call netcdf_catch(err)
end subroutine write_variable_4d_real


!> @brief Reads variable from netCDF dataset.
subroutine read_variable_1d_int_to_bool(ncid, name, buffer, start, counts)

  integer, intent(in) :: ncid !< NetCDF id.
  character(len=*), intent(in) :: name !< Variable name.
  logical(wl), dimension(:), allocatable, intent(inout) :: buffer
  integer, dimension(:), intent(in), optional :: start !< Corner indices
  integer, dimension(:), intent(in), optional :: counts !< Edge lengths.

  integer, dimension(:), allocatable :: buffer1d
  integer :: i

  call read_variable(ncid, name, buffer1d, start, counts)
  if (allocated(buffer)) deallocate(buffer)
  allocate(buffer(size(buffer1d)))
  do i = 1, size(buffer1d)
    if (buffer1d(i) .eq. 0) then
      buffer(i) = .false.
    else
      buffer(i) = .true.
    endif
  enddo
  deallocate(buffer1d)
end subroutine read_variable_1d_int_to_bool


!> @brief Determine if a variable is in a netCDF dataset.
function variable_exists(ncid, name) result(has_variable)

  integer, intent(in) :: ncid !< NetCDF id.
  character(len=*), intent(in) :: name !< Variable name.
  logical :: has_variable

  integer :: err
  integer :: varid

  err = nf90_inq_varid(ncid, trim(name), varid)
  if (err .eq. nf90_enotvar) then
    has_variable = .false.
  else
    call netcdf_catch(err)
    has_variable = .true.
  endif
end function variable_exists


end module netcdf_utils
