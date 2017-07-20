/*
 *  Copyright (C) 2003, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* $Id$ */

/*
 * This file implements the corresponding APIs defined in
 * src/dispatchers/dimension.c
 *
 * ncmpi_def_dim()    : dispatcher->def_dim()
 * ncmpi_inq_dimid()  : dispatcher->inq_dimid()
 * ncmpi_inq_dim()    : dispatcher->inq_dim()
 * ncmpi_rename_dim() : dispatcher->rename_dim()
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#ifdef HAVE_STDLIB_H
#include <stdlib.h>
#endif
#include <stdio.h>
#include <string.h>
#include <assert.h>

#include <mpi.h>

#include <pnc_debug.h>
#include <common.h>
#include "nc.h"
#include "ncx.h"
#include <utf8proc.h>

/*----< dup_NC_dim() >-------------------------------------------------------*/
static int
dup_NC_dim(const NC_dim *rdimp, NC_dim **dimp)
{
    *dimp = (NC_dim*) NCI_Malloc(sizeof(NC_dim));
    if (*dimp == NULL) DEBUG_RETURN_ERROR(NC_ENOMEM)

    (*dimp)->size     = rdimp->size;
    (*dimp)->name_len = strlen(rdimp->name)+1;
    (*dimp)->name     = (char*) NCI_Malloc((*dimp)->name_len);
    if ((*dimp)->name == NULL) DEBUG_RETURN_ERROR(NC_ENOMEM)
    strcpy((*dimp)->name, rdimp->name);

    return NC_NOERR;
}

#ifdef SEARCH_NAME_LINEARLY
/*----< NC_finddim() >-------------------------------------------------------*/
/*
 * Step thru NC_DIMENSION array, seeking match on name.
 * If found, set the dim ID pointed by dimidp, otherwise return NC_EBADDIM
 */
static int
NC_finddim(const NC_dimarray *ncap,
           const char        *name,  /* normalized dim name */
           int               *dimidp)
{
    int dimid;
    size_t nchars;

    if (ncap->ndefined == 0) return NC_EBADDIM;

    /* note that the number of dimensions allowed is < 2^32 */
    nchars = strlen(name);
    for (dimid=0; dimid<ncap->ndefined; dimid++) {
        if (ncap->value[dimid]->name_len == nchars &&
            strcmp(ncap->value[dimid]->name, name) == 0) {
            /* found the matched name */
            if (dimidp != NULL) *dimidp = dimid;
            return NC_NOERR; /* found it */
        }
    }
    return NC_EBADDIM; /* the name is not found */
}
#else
/*----< NC_finddim() >-------------------------------------------------------*/
/*
 * Search name from hash table ncap->nameT.
 * If found, set the dim ID pointed by dimidp, otherwise return NC_EBADDIM
 */
static int
NC_finddim(const NC_dimarray *ncap,
           const char        *name,  /* normalized dim name */
           int               *dimidp)
{
    int i, key, dimid;
    size_t nchars;

    if (ncap->ndefined == 0) return NC_EBADDIM;

    /* hash the dim name into a key for name lookup */
    key = HASH_FUNC(name);

    /* check the list using linear search */
    nchars = strlen(name);
    for (i=0; i<ncap->nameT[key].num; i++) {
        dimid = ncap->nameT[key].list[i];
        if (ncap->value[dimid]->name_len == nchars &&
            strcmp(name, ncap->value[dimid]->name) == 0) {
            if (dimidp != NULL) *dimidp = dimid;
            return NC_NOERR; /* the name already exists */
        }
    }
    return NC_EBADDIM; /* the name has never been used */
}
#endif

/* dimarray */

/*----< ncmpio_free_NC_dimarray() >------------------------------------------*/
inline void
ncmpio_free_NC_dimarray(NC_dimarray *ncap)
{
    int i;

    assert(ncap != NULL);
    if (ncap->nalloc == 0) return;

    assert(ncap->value != NULL);
    for (i=0; i<ncap->ndefined; i++) {
        NCI_Free(ncap->value[i]->name);
        NCI_Free(ncap->value[i]);
    }

    NCI_Free(ncap->value);
    ncap->value    = NULL;
    ncap->nalloc   = 0;
    ncap->ndefined = 0;

    /* free space allocated for dim name lookup table */
    ncmpio_hash_table_free(ncap->nameT);
}

/*----< ncmpio_dup_NC_dimarray() >-------------------------------------------*/
int
ncmpio_dup_NC_dimarray(NC_dimarray *ncap, const NC_dimarray *ref)
{
    int i, status=NC_NOERR;

    assert(ref != NULL);
    assert(ncap != NULL);

    if (ref->nalloc == 0) {
        ncap->nalloc   = 0;
        ncap->ndefined = 0;
        ncap->value    = NULL;
        return NC_NOERR;
    }

    /* allocate array of NC_dim objects */
    if (ref->nalloc > 0) {
        ncap->value = (NC_dim**) NCI_Calloc(ref->nalloc, sizeof(NC_dim*));
        if (ncap->value == NULL) DEBUG_RETURN_ERROR(NC_ENOMEM)
        ncap->nalloc = ref->nalloc;
    }

    /* duplicate each NC_dim objects */
    ncap->ndefined = 0;
    for (i=0; i<ref->ndefined; i++) {
        status = dup_NC_dim(ref->value[i], &ncap->value[i]);
        if (status != NC_NOERR) {
            ncmpio_free_NC_dimarray(ncap);
            return status;
        }
        ncap->ndefined++;
    }
    assert(ncap->ndefined == ref->ndefined);

    /* duplicate dim name lookup table */
    ncmpio_hash_table_copy(ncap->nameT, ref->nameT);

    return NC_NOERR;
}

/*----< incr_NC_dimarray() >-------------------------------------------------*/
/*
 * Add a new NC_dim handle to the end of an array of handles
 */
static int
incr_NC_dimarray(NC_dimarray *ncap,
                 NC_dim      *newdimp)
{
    NC_dim **vp;

    assert(ncap != NULL);

    if (ncap->nalloc == 0) {
        assert(ncap->ndefined == 0);
        vp = (NC_dim **) NCI_Malloc(NC_ARRAY_GROWBY * sizeof(NC_dim*));
        if (vp == NULL) DEBUG_RETURN_ERROR(NC_ENOMEM)

        ncap->value = vp;
        ncap->nalloc = NC_ARRAY_GROWBY;
    }
    else if (ncap->ndefined + 1 > ncap->nalloc) {
        vp = (NC_dim **) NCI_Realloc(ncap->value,
                         (ncap->nalloc + NC_ARRAY_GROWBY) * sizeof(NC_dim*));
        if (vp == NULL) DEBUG_RETURN_ERROR(NC_ENOMEM)

        ncap->value = vp;
        ncap->nalloc += NC_ARRAY_GROWBY;
    }
    /* else here means some space still available */

    if (newdimp != NULL) {
        ncap->value[ncap->ndefined] = newdimp;
        ncap->ndefined++;
    }
    return NC_NOERR;
}

/*----< inq_NC_dim() >-------------------------------------------------------*/
static NC_dim*
inq_NC_dim(const NC_dimarray *ncap, int dimid)
{
    /* returns the dimension ID defined earlier */
    assert(ncap != NULL);

    if (dimid < 0 || ncap->ndefined == 0 || dimid >= ncap->ndefined)
        return NULL;

    assert(ncap->value != NULL);

    return ncap->value[dimid];
}

/*----< ncmpio_def_dim() >---------------------------------------------------*/
int
ncmpio_def_dim(void       *ncdp,    /* IN:  NC object */
               const char *name,    /* IN:  name of dimension */
               MPI_Offset  size,    /* IN:  dimension size */
               int        *dimidp)  /* OUT: dimension ID */
{
    int dimid, err=NC_NOERR;
    char *nname=NULL;  /* normalized name */
    NC *ncp=(NC*)ncdp;
    NC_dim *dimp=NULL;

#if 0
    /* must be called in define mode */
    if (!NC_indef(ncp)) {
        DEBUG_ASSIGN_ERROR(err, NC_ENOTINDEFINE)
        goto err_check;
    }

    if (name == NULL || *name == 0) {
        DEBUG_ASSIGN_ERROR(err, NC_EBADNAME)
        goto err_check;
    }

    if (strlen(name) > NC_MAX_NAME) {
        DEBUG_ASSIGN_ERROR(err, NC_EMAXNAME)
        goto err_check;
    }

    /* check if the name string is legal for the netcdf format */
    err = ncmpii_check_name(name, ncp->format);
    if (err != NC_NOERR) {
        DEBUG_TRACE_ERROR
        goto err_check;
    }

    /* MPI_Offset is usually a signed value, but serial netcdf uses
     * size_t -- normally unsigned
     * In 1999 ISO C standard, size_t is an unsigned integer type of at least
     * 16 bit. */
    if (ncp->format == 2) { /* CDF-2 format, max is 2^32-4 */
        if (size > X_UINT_MAX - 3 || (size < 0))
            /* "-3" handles rounded-up size */
            err = NC_EDIMSIZE;
    } else if (ncp->format == 5) { /* CDF-5 format*/
        if (size < 0)
            err = NC_EDIMSIZE;
    } else { /* CDF-1 format, max is 2^31-4 */
        if (size > X_INT_MAX - 3 || (size < 0))
            /* "-3" handles rounded-up size */
            err = NC_EDIMSIZE;
    }
    if (err != NC_NOERR) {
        DEBUG_TRACE_ERROR
        goto err_check;
    }

    if (size == NC_UNLIMITED) {
#if 0
        /* check for any existing unlimited dimension, netcdf allows
         * one per file
         */
        dimid = ncmpio_find_NC_Udim(&ncp->dims, &dimp);
        if (dimid != -1) {
            DEBUG_ASSIGN_ERROR(err, NC_EUNLIMIT) /* found an existing one */
            goto err_check;
        }
#endif
        if (ncp->dims.unlimited_id != -1) {
            DEBUG_ASSIGN_ERROR(err, NC_EUNLIMIT) /* already defined */
            goto err_check;
        }
    }

    /* Note we no longer limit the number of dimensions, as CDF file formats
     * impose no such limit. Thus, the value of NC_MAX_DIMS has been changed
     * to NC_MAX_INT, as NC_dimarray.ndefined is of type signed int and so is
     * ndims argument in ncmpi_inq_varndims()
     */
    if (ncp->dims.ndefined == NC_MAX_DIMS) {
        DEBUG_ASSIGN_ERROR(err, NC_EMAXDIMS)
        goto err_check;
    }

    /* create a normalized character string */
    nname = (char *)ncmpii_utf8proc_NFC((const unsigned char *)name);
    if (nname == NULL) {
        DEBUG_ASSIGN_ERROR(err, NC_ENOMEM)
        goto err_check;
    }

    /* check if the name string is previously used */
    err = NC_finddim(&ncp->dims, nname, NULL);
    if (err != NC_EBADDIM) {
        DEBUG_ASSIGN_ERROR(err, NC_ENAMEINUSE)
        goto err_check;
    }
    else
        err = NC_NOERR;

err_check:
    if (ncp->safe_mode) {
        int root_name_len, rank, status, mpireturn;
        char *root_name=NULL;
        MPI_Offset root_size;

        /* check the error so far across processes */
        TRACE_COMM(MPI_Allreduce)(&err, &status, 1, MPI_INT, MPI_MIN,ncp->comm);
        if (mpireturn != MPI_SUCCESS) {
            NCI_Free(nname);
            return ncmpii_error_mpi2nc(mpireturn, "MPI_Allreduce");
        }
        if (status != NC_NOERR) {
            NCI_Free(nname);
            return status;
        }
        assert(nname != NULL);
        assert(name != NULL);

        MPI_Comm_rank(ncp->comm, &rank);

        /* check if name is consistent among all processes */
        root_name_len = 1;
        if (rank == 0 && name != NULL) root_name_len += strlen(name);
        TRACE_COMM(MPI_Bcast)(&root_name_len, 1, MPI_INT, 0, ncp->comm);
        if (mpireturn != MPI_SUCCESS) {
            if (nname != NULL) NCI_Free(nname);
            return ncmpii_error_mpi2nc(mpireturn, "MPI_Bcast root_name_len");
        }

        if (rank == 0)
            root_name = (char*)name;
        else
            root_name = (char*) NCI_Malloc((size_t)root_name_len);
        TRACE_COMM(MPI_Bcast)(root_name, root_name_len, MPI_CHAR, 0, ncp->comm);
        if (mpireturn != MPI_SUCCESS) {
            if (nname != NULL) NCI_Free(nname);
            if (rank > 0) NCI_Free(root_name);
            return ncmpii_error_mpi2nc(mpireturn, "MPI_Bcast");
        }
        if (err == NC_NOERR && rank > 0 && strcmp(root_name, name))
            DEBUG_ASSIGN_ERROR(err, NC_EMULTIDEFINE_DIM_NAME)
        if (rank > 0) NCI_Free(root_name);

        /* check if sizes are consistent across all processes */
        root_size = size;
        TRACE_COMM(MPI_Bcast)(&root_size, 1, MPI_OFFSET, 0, ncp->comm);
        if (mpireturn != MPI_SUCCESS) {
            if (nname != NULL) NCI_Free(nname);
            return ncmpii_error_mpi2nc(mpireturn, "MPI_Bcast");
        }
        if (err == NC_NOERR && root_size != size)
            DEBUG_ASSIGN_ERROR(err, NC_EMULTIDEFINE_DIM_SIZE)

        /* find min error code across processes */
        TRACE_COMM(MPI_Allreduce)(&err, &status, 1, MPI_INT, MPI_MIN,ncp->comm);
        if (mpireturn != MPI_SUCCESS) {
            NCI_Free(nname);
            return ncmpii_error_mpi2nc(mpireturn, "MPI_Allreduce");
        }
        if (status != NC_NOERR) {
            NCI_Free(nname);
            return status;
        }
    }

    if (err != NC_NOERR) {
        if (nname != NULL) NCI_Free(nname);
        return err;
    }
#else
    /* create a normalized character string */
    nname = (char *)ncmpii_utf8proc_NFC((const unsigned char *)name);
    if (nname == NULL) DEBUG_RETURN_ERROR(NC_ENOMEM)
#endif

    /* create a new dimension object (dimp->name points to nname) */
    dimp = (NC_dim*) NCI_Malloc(sizeof(NC_dim));
    if (dimp == NULL) {
        NCI_Free(nname);
        DEBUG_RETURN_ERROR(NC_ENOMEM)
    }
    dimp->size     = size;
    dimp->name     = nname;
    dimp->name_len = strlen(nname);

#ifndef SEARCH_NAME_LINEARLY
    ncmpio_hash_insert(ncp->dims.nameT, nname, ncp->dims.ndefined);
#endif

    /* Add a new dim handle to the end of handle array */
    err = incr_NC_dimarray(&ncp->dims, dimp);
    if (err != NC_NOERR) {
        NCI_Free(dimp->name);
        NCI_Free(dimp);
        DEBUG_RETURN_ERROR(err)
    }

    /* ncp->dims.ndefined has been increased in incr_NC_dimarray() */
    dimid = (int)ncp->dims.ndefined - 1;

    if (size == NC_UNLIMITED) ncp->dims.unlimited_id = dimid;

    if (dimidp != NULL) *dimidp = dimid;

    return err;
}

/*----< ncmpio_inq_dimid() >-------------------------------------------------*/
int
ncmpio_inq_dimid(void       *ncdp,
                 const char *name,
                 int        *dimid)
{
    int err=NC_NOERR;
    char *nname=NULL; /* normalized name */
    NC *ncp=(NC*)ncdp;

    if (name == NULL || *name == 0 || strlen(name) > NC_MAX_NAME)
        DEBUG_RETURN_ERROR(NC_EBADNAME)

    /* create a normalized character string */
    nname = (char *)ncmpii_utf8proc_NFC((const unsigned char *)name);
    if (nname == NULL) DEBUG_RETURN_ERROR(NC_ENOMEM)

    err = NC_finddim(&ncp->dims, nname, dimid);
    NCI_Free(nname);

    return err;
}

/*----< ncmpio_inq_dim() >---------------------------------------------------*/
int
ncmpio_inq_dim(void       *ncdp,
               int         dimid,
               char       *name,
               MPI_Offset *sizep)
{
    NC *ncp=(NC*)ncdp;
    NC_dim *dimp=NULL;

    dimp = inq_NC_dim(&ncp->dims, dimid);
    if (dimp == NULL) DEBUG_RETURN_ERROR(NC_EBADDIM)

    if (name != NULL)
        /* in PnetCDF, name is always NULL character terminated */
        strcpy(name, dimp->name);

    if (sizep != NULL) {
        if (dimp->size == NC_UNLIMITED)
            *sizep = ncp->numrecs;
        else
            *sizep = dimp->size;
    }
    return NC_NOERR;
}

/*----< ncmpio_rename_dim() >-------------------------------------------------*/
/* This API is collective and can be called in either define or data mode..
 * If the new name is longer than the old name, the netCDF dataset must be in
 * the define mode.
 */
int
ncmpio_rename_dim(void       *ncdp,
                  int         dimid,
                  const char *newname)
{
    int err=NC_NOERR;
    char *nnewname=NULL; /* normalized newname */
    size_t nnewname_len=0;
    NC *ncp=(NC*)ncdp;
    NC_dim *dimp=NULL;

#if 0
    /* check file's write permission */
    if (NC_readonly(ncp)) {
        DEBUG_ASSIGN_ERROR(err, NC_EPERM)
        goto err_check;
    }

    if (newname == NULL || *newname == 0) {
        DEBUG_ASSIGN_ERROR(err, NC_EBADNAME)
        goto err_check;
    }

    if (strlen(newname) > NC_MAX_NAME) {
        DEBUG_ASSIGN_ERROR(err, NC_EMAXNAME)
        goto err_check;
    }

    /* check whether newname is legal */
    err = ncmpii_check_name(newname, ncp->format);
    if (err != NC_NOERR) {
        DEBUG_TRACE_ERROR
        goto err_check;
    }
#endif
    /* create a normalized character string */
    nnewname = (char *)ncmpii_utf8proc_NFC((const unsigned char *)newname);
    if (nnewname == NULL) {
        DEBUG_ASSIGN_ERROR(err, NC_ENOMEM)
        goto err_check;
    }
    nnewname_len = strlen(nnewname);
#if 0
    /* check whether newname is already in use */
    err = NC_finddim(&ncp->dims, nnewname, NULL);
    if (err != NC_EBADDIM) { /* expecting NC_EBADDIM */
        DEBUG_ASSIGN_ERROR(err, NC_ENAMEINUSE)
        goto err_check;
    }
    else err = NC_NOERR;  /* reset err */
#endif
    /* retrieve dim object */
    dimp = inq_NC_dim(&ncp->dims, dimid);
    if (dimp == NULL) {
        DEBUG_ASSIGN_ERROR(err, NC_EBADDIM)
        goto err_check;
    }

    if (! NC_indef(ncp) && dimp->name_len < nnewname_len) {
        /* when in data mode, newname cannot be longer than the old one */
        DEBUG_ASSIGN_ERROR(err, NC_ENOTINDEFINE)
        goto err_check;
    }

#ifndef SEARCH_NAME_LINEARLY
    /* update dim name lookup table, by removing the old name and add
     * the new name */
    err = ncmpio_update_name_lookup_table(ncp->dims.nameT, dimid,
                             ncp->dims.value[dimid]->name, nnewname);
    if (err != NC_NOERR) {
        DEBUG_TRACE_ERROR
        goto err_check;
    }
#endif

err_check:
    if (ncp->safe_mode) {
        int status, mpireturn;

        /* check the error so far across processes */
        TRACE_COMM(MPI_Allreduce)(&err, &status, 1, MPI_INT, MPI_MIN,ncp->comm);
        if (mpireturn != MPI_SUCCESS) {
            NCI_Free(nnewname);
            return ncmpii_error_mpi2nc(mpireturn, "MPI_Allreduce");
        }
        if (status != NC_NOERR) {
            NCI_Free(nnewname);
            return status;
        }
#if 0
        int rank, root_name_len, root_dimid;
        char *root_name=NULL;
        assert(nnewname != NULL);

        MPI_Comm_rank(ncp->comm, &rank);

        /* check if newname is consistent among all processes */
        root_name_len = 1;
        if (rank == 0 && newname != NULL) root_name_len += strlen(newname);
        TRACE_COMM(MPI_Bcast)(&root_name_len, 1, MPI_INT, 0, ncp->comm);
        if (mpireturn != MPI_SUCCESS)
            return ncmpii_error_mpi2nc(mpireturn, "MPI_Bcast root_name_len");

        if (rank == 0)
            root_name = (char*)newname;
        else
            root_name = (char*) NCI_Malloc((size_t)root_name_len);
        TRACE_COMM(MPI_Bcast)(root_name, root_name_len, MPI_CHAR, 0, ncp->comm);
        if (mpireturn != MPI_SUCCESS) {
            if (rank > 0) NCI_Free(root_name);
            return ncmpii_error_mpi2nc(mpireturn, "MPI_Bcast");
        }
        if (err == NC_NOERR && rank > 0 && strcmp(root_name, newname))
            DEBUG_ASSIGN_ERROR(err, NC_EMULTIDEFINE_DIM_NAME)
        if (rank > 0) NCI_Free(root_name);

        /* check if dimid is consistent across all processes */
        root_dimid = dimid;
        TRACE_COMM(MPI_Bcast)(&root_dimid, 1, MPI_INT, 0, ncp->comm);
        if (mpireturn != MPI_SUCCESS)
            return ncmpii_error_mpi2nc(mpireturn, "MPI_Bcast");
        if (err == NC_NOERR && root_dimid != dimid)
            DEBUG_ASSIGN_ERROR(err, NC_EMULTIDEFINE_FNC_ARGS)

        /* find min error code across processes */
        TRACE_COMM(MPI_Allreduce)(&err, &status, 1, MPI_INT, MPI_MIN,ncp->comm);
        if (mpireturn != MPI_SUCCESS) {
            NCI_Free(nnewname);
            return ncmpii_error_mpi2nc(mpireturn, "MPI_Allreduce");
        }
        if (status != NC_NOERR) {
            NCI_Free(nnewname);
            return status;
        }
#endif
    }

    if (err != NC_NOERR) {
        if (nnewname != NULL) NCI_Free(nnewname);
        return err;
    }

    /* replace the old name with new name */
    assert(dimp != NULL);
    NCI_Free(dimp->name);
    dimp->name     = nnewname;
    dimp->name_len = nnewname_len;

    if (! NC_indef(ncp)) { /* when file is in data mode */
        /* Let root write the entire header to the file. Note that we cannot
         * just update the variable name in its space occupied in the file
         * header, because if the file space occupied by the name shrinks, all
         * the metadata following it must be moved ahead.
         */
        err = ncmpio_write_header(ncp);
        if (err != NC_NOERR) DEBUG_RETURN_ERROR(err)
    }

    return err;
}
