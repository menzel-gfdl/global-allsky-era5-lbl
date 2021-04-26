AC_PREREQ([2.59])
AC_INIT([global_allsky_era5_lbl], [alpha])
AC_CANONICAL_HOST
AC_CANONICAL_BUILD
AC_CANONICAL_TARGET
AM_INIT_AUTOMAKE([foreign dist-zip subdir-objects])
AC_CONFIG_MACRO_DIR([m4])
LT_PREREQ([2.4])
LT_INIT()
AC_PROG_FC
AC_PROG_INSTALL
AC_SYS_LARGEFILE
AC_LANG_PUSH(Fortran)
AC_FC_FREEFORM()
AC_FC_LINE_LENGTH([unlimited])
AC_CHECK_LIB([netcdff], [nf_create], [], [AC_MSG_ERROR([Cannot find netCDF4 fortran library. Set LDFLAGS.])])
AC_CONFIG_FILES([Makefile])
AC_OUTPUT()
