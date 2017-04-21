#
# Copyright (C) 2003, Northwestern University and Argonne National Laboratory
# See COPYRIGHT notice in top-level directory.
#
# $Id$
#
# @configure_input@

SUBDIRS = ncmpigen ncmpidump ncmpidiff ncmpivalid pnetcdf_version ncoffsets

# The script shows the end users how pnetcdf is built
bin_SCRIPTS = pnetcdf-config

EXTRA_DIST = pnetcdf-config.in \
             check_install.in

dist-hook:
	$(SED) -i "s/PNETCDF_RELEASE_DATE/`date '+%e %b %Y'`/" $(distdir)/pnetcdf-config.in

