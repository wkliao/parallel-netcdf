/*
 *  Copyright (C) 2015, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* $Id$ */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>  /* strcasecmp() */
#include <assert.h>
#include <errno.h>
#include <mpi.h>

#include <pnc_debug.h>
#include <common.h>
#include "nc.h"

/*----< ncmpio_set_pnetcdf_hints() >-----------------------------------------*/
/* this is where the I/O hints designated to pnetcdf are extracted */
void ncmpio_set_pnetcdf_hints(NC *ncp, MPI_Info info)
{
    char value[MPI_MAX_INFO_VAL];
    int  flag;

    if (info == MPI_INFO_NULL) return;

    /* nc_header_align_size, nc_var_align_size, and r_align * take effect when
     * a file is created or opened and later adding more header or variable
     * data */

    /* extract PnetCDF hints from user info object */
    MPI_Info_get(info, "nc_header_align_size", MPI_MAX_INFO_VAL-1, value,
                 &flag);
    if (flag) {
        errno = 0;  /* errno must set to zero before calling strtoll */
        ncp->h_align = strtoll(value,NULL,10);
        if (errno != 0) ncp->h_align = 0;
        else if (ncp->h_align < 0) ncp->h_align = 0;
    }

    MPI_Info_get(info, "nc_var_align_size", MPI_MAX_INFO_VAL-1, value, &flag);
    if (flag) {
        errno = 0;  /* errno must set to zero before calling strtoll */
        ncp->v_align = strtoll(value,NULL,10);
        if (errno != 0) ncp->v_align = 0;
        else if (ncp->v_align < 0) ncp->v_align = 0;
    }

    MPI_Info_get(info, "nc_record_align_size", MPI_MAX_INFO_VAL-1, value,
                 &flag);
    if (flag) {
        errno = 0;  /* errno must set to zero before calling strtoll */
        ncp->r_align = strtoll(value,NULL,10);
        if (errno != 0) ncp->r_align = 0;
        else if (ncp->r_align < 0) ncp->r_align = 0;
    }

    /* get header reading chunk size from info */
    MPI_Info_get(info, "nc_header_read_chunk_size", MPI_MAX_INFO_VAL-1, value,
                 &flag);
    if (flag) {
        errno = 0;  /* errno must set to zero before calling strtoll */
        ncp->chunk = strtoll(value,NULL,10);
        if (errno != 0) ncp->chunk = 0;
        else if (ncp->chunk < 0) ncp->chunk = 0;
    }

#ifdef ENABLE_SUBFILING
    MPI_Info_get(info, "pnetcdf_subfiling", MPI_MAX_INFO_VAL-1, value, &flag);
    if (flag && strcasecmp(value, "enable") == 0)
        ncp->subfile_mode = 1;

    MPI_Info_get(info, "nc_num_subfiles", MPI_MAX_INFO_VAL-1, value, &flag);
    if (flag) {
        errno = 0;
        ncp->num_subfiles = strtoll(value,NULL,10);
        if (errno != 0) ncp->num_subfiles = 0;
        else if (ncp->num_subfiles < 0) ncp->num_subfiles = 0;
    }
    if (ncp->subfile_mode == 0) ncp->num_subfiles = 0;
#endif
}

/*----< ncmpio_last_offset() >-----------------------------------------------*/
/* Returns the file offset of the last byte + 1 accessed by this request.
 * If count is NULL, this is equivalent to the starting offset of this
 * request. Note zero-length request should never call this subroutine.
 */
int
ncmpio_last_offset(const NC         *ncp,
                   const NC_var     *varp,
                   const MPI_Offset  start[],   /* [varp->ndims] */
                   const MPI_Offset  count[],   /* [varp->ndims] */
                   const MPI_Offset  stride[],  /* [varp->ndims] */
                   const int         reqMode,
                   MPI_Offset       *offset_ptr) /* OUT: file offset */
{
    MPI_Offset offset, *last_indx=NULL;
    int i, ndims, firstDim = 0;

    offset = varp->begin; /* beginning file offset of this variable */
    ndims  = varp->ndims; /* number of dimensions of this variable */

    if (ndims == 0) {
        *offset_ptr = varp->begin + varp->xsz;
        return NC_NOERR;
    }

    if (count != NULL) {
        last_indx = (MPI_Offset*) NCI_Malloc((size_t)ndims * SIZEOF_MPI_OFFSET);

        if (stride != NULL) {
            for (i=0; i<ndims; i++) {
                assert(count[i] > 0);
                last_indx[i] = start[i] + (count[i] - 1) * stride[i];
            }
        }
        else { /* stride == NULL */
            for (i=0; i<ndims; i++) {
                assert(count[i] > 0);
                last_indx[i] = start[i] + count[i] - 1;
            }
        }
    }
    else { /* when count == NULL stride is of no use */
        last_indx = (MPI_Offset*) start;
    }

    /* check whether last_indx is valid */

    firstDim = 0;
    /* check NC_EINVALCOORDS for record dimension */
    if (varp->shape[0] == NC_UNLIMITED) {
        if (ncp->format < 5 && last_indx[0] > NC_MAX_UINT) { /* CDF-1 and 2 */
            if (count != NULL) NCI_Free(last_indx);
            DEBUG_RETURN_ERROR(NC_EINVALCOORDS)
        }
        /* for record variable, [0] is the NC_UNLIMITED dimension */
        if (fIsSet(reqMode, NC_REQ_RD) && last_indx[0] >= ncp->numrecs) {
            /* read cannot go beyond current numrecs */
            if (count != NULL) NCI_Free(last_indx);
            DEBUG_RETURN_ERROR(NC_EINVALCOORDS)
        }
        firstDim = 1; /* done for checking the record dimension */
    }
    /* continue to check NC_EINVALCOORDS for the rest dimensions */
    for (i=firstDim; i<ndims; i++) {
        if (last_indx[i] < 0 || last_indx[i] >= varp->shape[i]) {
            if (count != NULL) NCI_Free(last_indx);
            DEBUG_RETURN_ERROR(NC_EINVALCOORDS)
        }
    }

    if (varp->shape[0] == NC_UNLIMITED)
        offset += last_indx[0] * ncp->recsize;
    else
        offset += last_indx[ndims-1] * varp->xsz;

    if (ndims > 1) {
        if (IS_RECVAR(varp))
            offset += last_indx[ndims - 1] * varp->xsz;
        else
            offset += last_indx[0] * varp->dsizes[1] * varp->xsz;

        for (i=1; i<ndims-1; i++)
            offset += last_indx[i] * varp->dsizes[i+1] * varp->xsz;
    }

    if (count != NULL) NCI_Free(last_indx);

    *offset_ptr = offset;
    return NC_NOERR;
}

