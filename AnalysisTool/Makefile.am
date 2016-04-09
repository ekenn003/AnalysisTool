ACLOCAL_AMFLAGS = -I m4
#AM_CXXFLAGS = $(INTI_CFLAGS)
AM_CXXFLAGS = -Werror
lib_LTLIBRARIES = libAnalysisTool.la
libAnalysisTool_la_SOURCES = Analyse.cc AnalysisBox.cc DataForm.cc Analyse.h AnalysisBox.h DataForm.h
libAnalysisTool_la_CPPFLAGS = $(ROOT_CFLAGS) $(BOOST_CFLAGS) $(MPI_CFLAGS)
libAnalysisTool_la_LIBADD = $(ROOT_LIBS) $(BOOST_LIBS) $(MPI_LIBS)
AnalysisTool_includedir = $(includedir)/AnalysisTool
nobase_AnalysisTool_include_HEADERS = Analyse.h AnalysisBox.h DataForm.h


bin_SCRIPTS = AnalysisToolUseThis
CLEANFILES = $(bin_SCRIPTS)
EXTRA_DIST = AnalysisToolUseThis.in

AnalysisToolUseThis: AnalysisToolUseThis.in
	sed -e 's,[@]prefix[@],$(prefix),g' -e 's,[@]VERSION[@],$(VERSION),g' < AnalysisToolUseThis.in > AnalysisToolUseThis

bin_PROGRAMS = merger
merger_SOURCES = src/merger.cc
AM_CPPFLAGS = $(ROOT_CFLAGS) $(BOOST_CFLAGS) $(MPI_CFLAGS)
merger_LDADD = $(ROOT_LIBS) $(BOOST_LIBS) $(MPI_LIBS) -lAnalysisTool

