#
# Copyright (C) 2012, Northwestern University and Argonne National Laboratory
# See COPYRIGHT notice in top-level directory.
#
# $Id$
#
# @configure_input@

AM_CPPFLAGS = -I$(top_builddir)/src/drivers/ncmpi -I$(top_srcdir)/src/drivers/ncmpi -I${top_srcdir}/src/include

bin_PROGRAMS = ncmpivalid
ncmpivalid_SOURCES = ncmpivalid.c
ncmpivalid_LDADD = $(top_builddir)/src/libs/libpnetcdf.la

# According Automake Manual, Section 8.3.5 Libtool Convenience Libraries,
# using nodist_EXTRA to trick automake to use C++ compiler as the linker
if HAS_MPICXX
nodist_EXTRA_ncmpivalid_SOURCES = dummycxx.cxx
endif

$(top_builddir)/src/libs/libpnetcdf.la:
	set -e; cd $(top_builddir)/src/libs && $(MAKE) $(MFLAGS)

dist_man_MANS = ncmpivalid.1

CLEANFILES = *.gcno

tests-local: all

dist-hook:
	$(SED) -i "s/PNETCDF_RELEASE_VERSION/$(PNETCDF_VERSION)/" $(distdir)/ncmpivalid.1
	$(SED) -i "s/PNETCDF_RELEASE_DATE/`date '+%e %b %Y'`/"    $(distdir)/ncmpivalid.1
	$(SED) -i "s/PNETCDF_RELEASE_DATE2/`date '+%Y-%m-%d'`/"   $(distdir)/ncmpivalid.1

