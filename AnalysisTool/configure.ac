## LibExample:
## an example of using Automake and Libtool to build a shared library

AC_INIT([LibAnalysisTool],[1.7.1],[hindrichs@physik.rwth-aachen.de])
#AM_INIT_AUTOMAKE([foreign -Wall -Werror])
AM_INIT_AUTOMAKE
#AC_PREREQ([2.67])

AC_ARG_ENABLE(mpi, [AS_HELP_STRING([--enable-mpi], [use mpi])],
		[MPI_CFLAGS="-DUSE_MPI -I$MPI_INCLUDE"
		 MPI_LIBS="-L$MPI_LIBDIR -lmpi -lmpi_cxx"
],
[])
AC_SUBST([MPI_CFLAGS])
AC_SUBST([MPI_LIBS])

AC_ARG_WITH([root-include-path],
	      [AS_HELP_STRING([--with-root-include-path], [location of the root headers])],
	      [ROOT_CFLAGS="-I$withval"],
	      [ROOT_CFLAGS='-I $(shell root-config --incdir)'])
AC_SUBST([ROOT_CFLAGS])

AC_ARG_WITH([root-lib-path],
	      [AS_HELP_STRING([--with-root-lib-path], [location of the root libraries])],
	      [ROOT_LIBS="-L$withval $(shell root-config --libs) -lMinuit -lMathCore -lrt"],
              [ROOT_LIBS='$(shell root-config --libs) -lMinuit -lMathCore -lrt'])
AC_SUBST([ROOT_LIBS])

AC_ARG_WITH([boost-include-path],
  [AS_HELP_STRING([--with-boost-include-path],
    [location of the boost headers])],
  [BOOST_CFLAGS="-I$withval"],
  [BOOST_CFLAGS=''])
AC_SUBST([BOOST_CFLAGS])

AC_ARG_WITH([boost-lib-path],
  [AS_HELP_STRING([--with-boost-lib-path], [location of the boost libraries])],
  [BOOST_LIBS="-L$withval -lboost_regex"],
  [BOOST_LIBS='-lboost_regex'])
AC_SUBST([BOOST_LIBS])

AC_CONFIG_MACRO_DIR([m4])

AC_PROG_CXX
AC_PROG_CC
AC_PROG_CPP
AC_PROG_INSTALL
AC_PROG_LN_S
AC_PROG_MAKE_SET
AC_PROG_LIBTOOL
#AM_INIT_AUTOMAKE([1.10 -Wall no-define])

#AC_SUBST(INTI_CFLAGS)
#AC_SUBST(INTI_LIBS)

#LT_INIT

#AC_SUBST([EXAMPLE_SO_VERSION], [1:1:1])
#AC_SUBST([EXAMPLE_API_VERSION], [1.0])

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
                  
