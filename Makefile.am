lib_LTLIBRARIES = libclouds.la
bin_PROGRAMS = global_allsky_era5_lbl

AM_FCFLAGS = -I$(top_srcdir)/dependencies/GRTCODE/fortran-bindings \
             -I$(top_srcdir)/dependencies/rte-rrtmgp/build

#Cloud optics library.
clouds = $(top_srcdir)/src/clouds
libclouds_la_SOURCES = $(clouds)/hu_stamnes.f90 \
                       $(clouds)/ice_cloud_optics.f90 \
                       $(clouds)/incomplete_beta.f90 \
                       $(clouds)/netcdf_utils.f90 \
                       $(clouds)/optics_utils.f90 \
                       $(clouds)/slingo.f90 \
                       $(clouds)/stochastic_clouds.f90
#libclouds_la_LIBADD = -lrte -lnetcdff
$(clouds)/hu_stamnes.lo: $(clouds)/optics_utils.lo \
                         $(clouds)/netcdf_utils.lo
$(clouds)/ice_cloud_optics.lo: $(clouds)/optics_utils.lo \
                               $(clouds)/netcdf_utils.lo
$(clouds)/incomplete_beta.lo: $(clouds)/netcdf_utils.lo
$(clouds)/slingo.lo: $(clouds)/optics_utils.lo \
                     $(clouds)/netcdf_utils.lo
$(clouds)/stochastic_clouds.lo: $(clouds)/incomplete_beta.lo


#Global allsky ERA5 line-by-line executable.
era5 = $(top_srcdir)/src
global_allsky_era5_lbl_SOURCES = $(era5)/argparse.f90 \
                                 $(era5)/era5.f90 \
                                 $(era5)/main.f90 \
                                 $(top_srcdir)/dependencies/rte-rrtmgp/examples/mo_simple_netcdf.F90 \
                                 $(top_srcdir)/dependencies/rte-rrtmgp/examples/mo_load_coefficients.F90
global_allsky_era5_lbl_LDADD = $(top_builddir)/libclouds.la \
                               $(top_srcdir)/dependencies/GRTCODE/fortran-bindings/libgrtcode_fortran.a \
                               $(top_srcdir)/dependencies/rte-rrtmgp/build/librrtmgp.a \
                               $(top_srcdir)/dependencies/rte-rrtmgp/build/librte.a \
                               -lnetcdff
$(era5)/era5.$(OBJEXT): $(era5)/argparse.$(OBJEXT)
$(top_srcdir)/dependencies/rte-rrtmgp/examples/mo_load_coefficients.$(OBJEXT): $(top_srcdir)/dependencies/rte-rrtmgp/examples/mo_simple_netcdf.$(OBJEXT)
$(era5)/main.$(OBJEXT): $(era5)/era5.$(OBJEXT) \
                        $(era5)/argparse.$(OBJEXT) \
                        $(top_srcdir)/dependencies/rte-rrtmgp/examples/mo_load_coefficients.$(OBJEXT)


CLEANFILES = *.mod *__genmod.f90
