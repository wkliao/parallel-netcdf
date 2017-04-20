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

dist_man_MANS = ncmpidiff.1

CLEANFILES = *.gcno

dist-hook:
	sed -i "s/PNETCDF_RELEASE_VERSION/$(PNETCDF_VERSION)/" $(distdir)/ncmpidiff.1
	sed -i "s/PNETCDF_RELEASE_DATE/`date '+%e %b %Y'`/"    $(distdir)/ncmpidiff.1
	sed -i "s/PNETCDF_RELEASE_DATE2/`date '+%Y-%m-%d'`/"   $(distdir)/ncmpidiff.1

tests-local: all

.PHONY: $(top_builddir)/src/libs/libpnetcdf.la
