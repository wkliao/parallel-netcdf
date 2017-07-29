/*
 *  Copyright (C) 2017, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* $Id$ */

#ifndef _COMMON_H
#define _COMMON_H

#include <mpi.h>
#include <pnetcdf.h>

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

/* useful for aligning memory */
#define _RNDUP(x, unit)      ((((x) + (unit) - 1) / (unit)) * (unit))
#define _RNDDOWN(x, unit)    ((x) - ((x)%(unit)))

/* #define M_RND_UNIT   (sizeof(double))
 * SIZEOF_DOUBLE is defined in config.h
 */
#define M_RND_UNIT        SIZEOF_DOUBLE
#define M_RNDUP(x)        _RNDUP(x, M_RND_UNIT)
#define M_RNDDOWN(x)      _RNDDOWN(x, M_RND_UNIT)
#define D_RNDUP(x, align) _RNDUP(x, (off_t)(align))

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

extern int
ncmpii_check_name(const char *name, int file_ver);

extern MPI_Datatype
ncmpii_nc2mpitype(nc_type xtype);

extern int
ncmpii_xlen_nc_type(nc_type xtype, int *size);

extern int
ncmpii_buftype_decode(int ndims, nc_type xtype, const MPI_Offset *count,
                      MPI_Datatype buftype, MPI_Datatype *ptype,
                      MPI_Offset *bufcount, MPI_Offset *bnelems,
                      MPI_Offset *nbytes, int *el_size,
                      int *buftype_is_contig);

extern int
ncmpii_pack(int ndims, const MPI_Offset *count, const MPI_Offset *imap,
            void *buf, MPI_Offset bufcount, MPI_Datatype buftype,
            MPI_Offset *bnelems, MPI_Datatype *ptype, void **cbuf);
#endif
