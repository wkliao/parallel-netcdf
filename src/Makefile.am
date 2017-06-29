#
# Copyright (C) 2003, Northwestern University and Argonne National Laboratory
# See COPYRIGHT notice in top-level directory.
#
# $Id$
#
# @configure_input@

# adding . to the list of SUBDIRS causes the directories built before .
# This is useful for "make -j" option, as libs should be built after all
# sub-libraries
SUBDIRS = include dispatchers drivers binding packaging . libs utils

# For VPATH build (parallel build), try delete all sub-directories
RMDIRS = include dispatchers drivers binding packaging
distclean-local:
	if [ "$(abs_builddir)" != "$(abs_srcdir)" ] ; then \
	   for d in $(RMDIRS) ; do \
	       if [ "$$d" = "." ] ; then \
	          continue ; \
	       fi ; \
	       rmdir $$d || true ; \
	   done ; \
	fi
