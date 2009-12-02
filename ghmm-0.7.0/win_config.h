/* win_config.h   Windows specific settings.  */

/* Define to empty if the keyword does not work.  */
#undef const

/* Define to `unsigned' if <sys/types.h> doesn't define.  */
#undef size_t

/* Define if you have the ANSI C header files.  */
#define STDC_HEADERS

/* Define if you can safely include both <sys/time.h> and <time.h>.  */
#define TIME_WITH_SYS_TIME

/* Define if your <sys/time.h> declares struct tm.  */
#undef TM_IN_SYS_TIME

/* Define where phi data file lies */
#define PHI_DATA_FILE "ghmm/PHI_001_20.dat"

/* Define if you have the <cerrno> header file.  */
#define HAVE_CERRNO

/* Define if you have the <cmath> header file.  */
#define HAVE_CMATH

/* Define if you have the <cstdarg> header file.  */
#define HAVE_CSTDARG

/* Define if you have the <cstdio> header file.  */
#define HAVE_CSTDIO

/* Define if you have the <cstdlib> header file.  */
#define HAVE_CSTDLIB

/* Define if you have the <cstring> header file.  */
#define HAVE_CSTRING

/* Define if you have the <dlfcn.h> header file.  */
#define HAVE_DLFCN_H

/* Define if you have the pthread library (-lpthread).  */
#undef HAVE_LIBPTHREAD

/* Name of package */
#define PACKAGE "ghmm"

/* Version number of package */
#define VERSION 0.3.0

/* include experimental features */
#undef __EXPERIMENTAL__

/* define if the compiler implements namespaces */
#define HAVE_NAMESPACES

/* define if the compiler supports Standard Template Library */
#define HAVE_STL

/* struct gsl_interval exists */
#undef HAVE_GSL_INTERVAL

/* root solver allocation takes only one argument */
#undef GSL_ROOT_FSLOVER_ALLOC_WITH_ONE_ARG

/* gsl_histogram_set_ranges_uniform is defined */
#undef GSL_HISTOGRAM_SET_RANGES_UNIFORM

/* gsl_sf_erf exists */
#define HAVE_GSL_SF_ERF

/* gsl_sf_erfc exists */
#define HAVE_GSL_SF_ERFC

/* use GSL functions instead of ghmm interpolation algorithms */
#define DO_WITH_GSL

