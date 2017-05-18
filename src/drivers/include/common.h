/*
 *  Copyright (C) 2017, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* $Id$ */

#ifndef _COMMON_H
#define _COMMON_H

#ifndef MAX
#define MAX(mm,nn) (((mm) > (nn)) ? (mm) : (nn))
#endif
#ifndef MIN
#define MIN(mm,nn) (((mm) < (nn)) ? (mm) : (nn))
#endif

#ifdef PNETCDF_DEBUG
#define DEBUG_RETURN_ERROR(err) {                                       \
    char *_env_str = getenv("PNETCDF_VERBOSE_DEBUG_MODE");              \
    if (_env_str != NULL && *_env_str != '0') {                         \
        int _rank;                                                      \
        MPI_Comm_rank(MPI_COMM_WORLD, &_rank);                          \
        fprintf(stderr, "Rank %d: %s error at line %d of %s in %s\n",   \
        _rank,ncmpi_strerrno(err),__LINE__,__func__,__FILE__);          \
    }                                                                   \
    return err;                                                         \
}
#define DEBUG_ASSIGN_ERROR(status, err) {                               \
    char *_env_str = getenv("PNETCDF_VERBOSE_DEBUG_MODE");              \
    if (_env_str != NULL && *_env_str != '0') {                         \
        int _rank;                                                      \
        MPI_Comm_rank(MPI_COMM_WORLD, &_rank);                          \
        fprintf(stderr, "Rank %d: %s error at line %d of %s in %s\n",   \
        _rank,ncmpi_strerrno(err),__LINE__,__func__,__FILE__);          \
    }                                                                   \
    status = err;                                                       \
}
#define DEBUG_TRACE_ERROR {                                             \
    char *_env_str = getenv("PNETCDF_VERBOSE_DEBUG_MODE");              \
    if (_env_str != NULL && *_env_str != '0') {                         \
        int _rank;                                                      \
        MPI_Comm_rank(MPI_COMM_WORLD, &_rank);                          \
        fprintf(stderr, "Rank %d: %s error at line %d of %s in %s\n",   \
        _rank,ncmpi_strerrno(err),__LINE__,__func__,__FILE__);          \
    }                                                                   \
}
#else
#define DEBUG_RETURN_ERROR(err) return err;
#define DEBUG_ASSIGN_ERROR(status, err) status = err;
#define DEBUG_TRACE_ERROR
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

#endif
