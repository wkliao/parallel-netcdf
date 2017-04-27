#
# Copyright (C) 2012, Northwestern University and Argonne National Laboratory
# See COPYRIGHT notice in top-level directory.
#
# $Id$
#
# @configure_input@

bin_PROGRAMS = ncmpidiff
ncmpidiff_SOURCES = ncmpidiff.c
ncmpidiff_LDADD = $(top_builddir)/src/libs/libpnetcdf.la

# According Automake Manual, Section 8.3.5 Libtool Convenience Libraries,
# using nodist_EXTRA to trick automake to use C++ compiler as the linker
if HAS_MPICXX
nodist_EXTRA_ncmpidiff_SOURCES = dummycxx.cxx
endif

$(top_builddir)/src/libs/libpnetcdf.la:
	set -e; cd $(top_builddir)/src/libs && $(MAKE) $(MFLAGS)

dist_man_MANS = ncmpidiff.1

CLEANFILES = *.gcno

dist-hook:
	$(SED) -i "s/PNETCDF_RELEASE_VERSION/$(PNETCDF_VERSION)/" $(distdir)/ncmpidiff.1
	$(SED) -i "s/PNETCDF_RELEASE_DATE/`date '+%e %b %Y'`/"    $(distdir)/ncmpidiff.1
	$(SED) -i "s/PNETCDF_RELEASE_DATE2/`date '+%Y-%m-%d'`/"   $(distdir)/ncmpidiff.1

tests-local: all

