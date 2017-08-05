/*
 *  Copyright (C) 2015, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* $Id$ */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <stdio.h>
#ifdef HAVE_STDLIB_H
#include <stdlib.h>
#endif
#include <string.h>
#include <strings.h>  /* strcasecmp() */
#include <assert.h>
#include <mpi.h>

#include <pnc_debug.h>
#include <common.h>
#include "nc.h"
#include "ncx.h"

#if 0
/*----< ncmpio_sanity_check() >----------------------------------------------*/
/* check the following errors and in that precedence.
 * NC_EBADID, NC_EPERM, NC_EINDEFINE, NC_EINDEP/NC_ENOTINDEP, NC_ENOTVAR,
 * NC_ECHAR, NC_EINVALCOORDS, NC_EEDGE, NC_ESTRIDE, NC_EINVAL.
 */
int ncmpio_sanity_check(NC                *ncp,
                        int               varid,
                        const MPI_Offset  bufcount,
                        MPI_Datatype      buftype,  /* internal datatype */
                        int               reqMode,
                        NC_var          **varp)  /* OUT */
{
    /* all errors detected here are fatal, must return immediately */
    int err=NC_NOERR;

    /* check file write permission if this is write request */
    if (fIsSet(reqMode, NC_REQ_WR) && NC_readonly(ncp))
        DEBUG_RETURN_ERROR(NC_EPERM)

    if (fIsSet(reqMode, NC_REQ_BLK)) {
        /* blocking APIs must be called in data mode */
        if (NC_indef(ncp))
            DEBUG_RETURN_ERROR(NC_EINDEFINE)

        /* for blocking APIs, check if in the right collective or independent
         * mode, nonblocking APIs can be called in either mode */
        if (fIsSet(reqMode, NC_REQ_INDEP) && !NC_indep(ncp))
            DEBUG_RETURN_ERROR(NC_ENOTINDEP)
        else if (fIsSet(reqMode, NC_REQ_COLL) && NC_indep(ncp))
            DEBUG_RETURN_ERROR(NC_EINDEP)
    }

    if (fIsSet(reqMode, NC_REQ_ZERO)) return NC_NOERR;

    /* check if varid is valid (check NC_ENOTVAR) */
    err = ncmpio_NC_lookupvar(ncp, varid, varp);
    if (err != NC_NOERR) return err;

    /* check NC_ECHAR */
    if (fIsSet(reqMode, NC_REQ_FLEX)) {
        /* when buftype == MPI_DATATYPE_NULL, bufcount is ignored and this API
         * assumes argument buf's data type matches the data type of variable
         * defined in the file - no data conversion will be done.
         */
        if (buftype != MPI_DATATYPE_NULL) {
            int isderived, el_size, buftype_is_contig;
            MPI_Datatype ptype;
            MPI_Offset   bnelems=0;

            err = ncmpii_dtype_decode(buftype, &ptype, &el_size, &bnelems,
                                      &isderived, &buftype_is_contig);
            if (err != NC_NOERR) return err;

            err = NCMPII_ECHAR((*varp)->type, ptype);
            if (err != NC_NOERR) return err;
        }
        /* else case: itype matches xtype */

        /* for flexible APIs, bufcount cannot be negative */
        if (bufcount < 0) DEBUG_RETURN_ERROR(NC_EINVAL)
    }
    else { /* called from a high-level API */
        err = NCMPII_ECHAR((*varp)->type, buftype);
        if (err != NC_NOERR) return err;
    }
    return NC_NOERR;
}
#endif

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

