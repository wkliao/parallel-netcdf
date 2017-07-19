/*
 *  Copyright (C) 2017, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* $Id$ */

#ifndef _COMMON_H
#define _COMMON_H

#include <mpi.h>

/*
 * Macros for dealing with flag bits.
 */
#define fSet(t, f)      ((t) |=  (f))
#define fClr(t, f)      ((t) &= ~(f))
#define fIsSet(t, f)    ((t) &   (f))
#define fMask(t, f)     ((t) & ~ (f))

#ifndef MAX
#define MAX(mm,nn) (((mm) > (nn)) ? (mm) : (nn))
#endif
#ifndef MIN
#define MIN(mm,nn) (((mm) < (nn)) ? (mm) : (nn))
#endif

extern void *
NCI_Malloc_fn(size_t size, const int lineno, const char *func,
              const char *filename);

extern void *
NCI_Calloc_fn(size_t nelem, size_t elsize, const int lineno, const char *func,
              const char *filename);

extern void *
NCI_Realloc_fn(void *ptr, size_t size, const int lineno, const char *func,
               const char *filename);

extern void
NCI_Free_fn(void *ptr, const int lineno, const char *func,
            const char *filename);

#define NCI_Malloc(a)    NCI_Malloc_fn(a,__LINE__,__func__,__FILE__)
#define NCI_Calloc(a,b)  NCI_Calloc_fn(a,b,__LINE__,__func__,__FILE__)
#define NCI_Realloc(a,b) NCI_Realloc_fn(a,b,__LINE__,__func__,__FILE__)
#define NCI_Free(a)      NCI_Free_fn(a,__LINE__,__func__,__FILE__)

extern int
ncmpii_inq_malloc_size(size_t *size);

extern int
ncmpii_inq_malloc_max_size(size_t *size);

extern int
ncmpii_inq_malloc_list(void);

extern int
ncmpii_dtype_decode(MPI_Datatype dtype, MPI_Datatype *ptype, int *el_size,
                    MPI_Offset *nelems, int *isderived,
                    int *iscontig_of_ptypes);

extern int
ncmpii_create_imaptype(int ndims, const MPI_Offset *count,
                       const MPI_Offset *imap, MPI_Datatype ptype,
                       MPI_Datatype *imaptype);

extern int
ncmpii_error_mpi2nc(int mpi_errorcode, char *msg);

extern int
ncmpii_start_count_stride_check(int format, int api, int ndims, int numrecs,
                const MPI_Offset *shape, const MPI_Offset *start,
                const MPI_Offset *count, const MPI_Offset *stride,
                const int rw_flag);

#endif
