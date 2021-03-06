#
# Copyright (C) 2003, Northwestern University and Argonne National Laboratory
# See COPYRIGHT notice in top-level directory.
#
# $Id$
#
# @configure_input@

SUFFIXES = .o .c .f

AM_DEFAULT_SOURCE_EXT = .c

AM_CPPFLAGS  = -I$(srcdir)/../common
AM_CPPFLAGS += -I$(top_srcdir)/src/include
AM_CPPFLAGS += -I$(top_builddir)/src/include
LDADD  = $(top_builddir)/src/libs/libpnetcdf.la
LDADD += ../common/libtestutils.la

AM_CFLAGS =
if IS_BIGENDIAN
   # Do not add to AM_CPPFLAGS, as it will also be used by Fortran programs
   AM_CFLAGS += -DWORDS_BIGENDIAN
   # AM_FFLAGS  += $(FC_DEFINE)WORDS_BIGENDIAN
   # AM_FCFLAGS += $(FC_DEFINE)WORDS_BIGENDIAN
endif

if ENABLE_LARGE_REQ
   AM_CFLAGS += -DENABLE_LARGE_REQ
endif

check_PROGRAMS = large_files \
                 large_var \
                 large_dims_vars_attrs \
                 high_dim_var \
                 tst_cdf5_begin \
                 large_coalesce

if HAS_FORTRAN
   check_PROGRAMS += bigrecords
   AM_FFLAGS = -I$(top_builddir)/src/binding/f77
   bigrecords_SOURCES = bigrecords.f
endif

EXTRA_DIST = wrap_runs.sh

NC_FILES = $(check_PROGRAMS:%=$(TESTOUTDIR)/%.nc)

CLEANFILES = $(NC_FILES) core core.* *.gcda *.gcno *.gcov gmon.out

../common/libtestutils.la:
	set -e; cd ../common && $(MAKE) $(MFLAGS) tests

if RUN_LARGE_FILE_TEST
TESTPROGRAMS = $(check_PROGRAMS)
else
TESTPROGRAMS =
endif

# autimake 1.11.3 has not yet implemented AM_TESTS_ENVIRONMENT
# For newer versions, we can use AM_TESTS_ENVIRONMENT instead
# AM_TESTS_ENVIRONMENT  = TESTPROGRAMS="$(TESTPROGRAMS)" ; export TESTPROGRAMS;
# AM_TESTS_ENVIRONMENT += TESTSEQRUN="$(TESTSEQRUN)" ; export TESTSEQRUN;
# AM_TESTS_ENVIRONMENT += TESTOUTDIR="$(TESTOUTDIR)" ; export TESTOUTDIR;
TESTS_ENVIRONMENT  = TESTPROGRAMS="$(TESTPROGRAMS)" ; export TESTPROGRAMS;
TESTS_ENVIRONMENT += TESTSEQRUN="$(TESTSEQRUN)" ; export TESTSEQRUN;
TESTS_ENVIRONMENT += TESTOUTDIR="$(TESTOUTDIR)" ; export TESTOUTDIR;

TEST_EXTENSIONS = .sh
LOG_COMPILER = $(srcdir)/wrap_runs.sh
SH_LOG_COMPILER =

TESTS = $(TESTPROGRAMS)

# Some of these tests are designed to run on one processes,
# Run them on 4 processes to see if they can handle well

TESTMPIRUN2  = `echo $(TESTMPIRUN) | $(SED) -e 's/NP/2/g'`
TESTMPIRUN4  = `echo $(TESTMPIRUN) | $(SED) -e 's/NP/4/g'`
TESTMPIRUN6  = `echo $(TESTMPIRUN) | $(SED) -e 's/NP/6/g'`

ptest ptest4: $(TESTPROGRAMS)
	for j in 0 ; do { \
	export PNETCDF_SAFE_MODE=$$j ; \
	set -e ; for i in $(TESTPROGRAMS) ; do ( \
	$(RM) -f $(TESTOUTDIR)/$$i.nc ; \
	$(TESTMPIRUN4) ./$$i $(TESTOUTDIR)/$$i.nc ; \
	$(RM) -f $(TESTOUTDIR)/$$i.nc ; \
	) ; done ; } ; done

ptest2: $(TESTPROGRAMS)
	for j in 0 ; do { \
	export PNETCDF_SAFE_MODE=$$j ; \
	set -e ; for i in $(TESTPROGRAMS); do ( \
	$(RM) -f $(TESTOUTDIR)/$$i.nc ; \
	$(TESTMPIRUN2) ./$$i $(TESTOUTDIR)/$$i.nc ; \
	$(RM) -f $(TESTOUTDIR)/$$i.nc ; \
	) ; done ; } ; done

ptest6: $(TESTPROGRAMS)
	for j in 0 ; do { \
	export PNETCDF_SAFE_MODE=$$j ; \
	set -e ; for i in $(TESTPROGRAMS); do ( \
	$(RM) -f $(TESTOUTDIR)/$$i.nc ; \
	$(TESTMPIRUN6) ./$$i $(TESTOUTDIR)/$$i.nc ; \
	$(RM) -f $(TESTOUTDIR)/$$i.nc ; \
	) ; done ; } ; done

ptests: ptest2 ptest4 ptest6
ptest8 ptest10:

# build check targets but not invoke
tests-local: all $(check_PROGRAMS)

.PHONY: ptest ptests ptest2 ptest4 ptest6 ptest8 ptest10

