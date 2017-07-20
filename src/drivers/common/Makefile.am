#
# Copyright (C) 2012, Northwestern University and Argonne National Laboratory
# See COPYRIGHT notice in top-level directory.
#
# $Id$
#
# @configure_input@

SUFFIXES = .a .o .c .m4 .h

AM_CPPFLAGS  = -I${top_srcdir}/src/include
AM_CPPFLAGS += -I${top_srcdir}/src/include
AM_CPPFLAGS += -I${top_srcdir}/src/drivers/include

if PNETCDF_DEBUG
   AM_CPPFLAGS += -DPNETCDF_DEBUG
endif

noinst_LTLIBRARIES = libcommon.la

H_SRCS = utf8proc_data.h

C_SRCS = utf8proc.c \
         mem_alloc.c \
         dtype_decode.c \
         create_imaptype.c \
         error_mpi2nc.c \
         check_start_count_stride.c \
         check_name.c

libcommon_la_SOURCES = $(C_SRCS) $(H_SRCS)

CLEANFILES = *.gcno

tests-local: all

