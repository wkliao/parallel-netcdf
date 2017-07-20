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
    if (flag && strcasecmp(value, "disable") == 0)
        ncp->subfile_mode = 0;

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

/*----< ncmpio_xlen_nc_type() >----------------------------------------------*/
/* return the length of external NC data type */
int
ncmpio_xlen_nc_type(nc_type type) {
    switch(type) {
        case NC_BYTE:
        case NC_CHAR:
        case NC_UBYTE:  return X_SIZEOF_CHAR;
        case NC_SHORT:  return X_SIZEOF_SHORT;
        case NC_USHORT: return X_SIZEOF_USHORT;
        case NC_INT:    return X_SIZEOF_INT;
        case NC_UINT:   return X_SIZEOF_UINT;
        case NC_FLOAT:  return X_SIZEOF_FLOAT;
        case NC_DOUBLE: return X_SIZEOF_DOUBLE;
        case NC_INT64:  return X_SIZEOF_INT64;
        case NC_UINT64: return X_SIZEOF_UINT64;
        default: DEBUG_RETURN_ERROR(NC_EBADTYPE);
    }
    return NC_NOERR;
}

/*----< ncmpio_last_offset() >-----------------------------------------------*/
/* returns the file offset of the last byte accessed of this request
 * If counts is NULL, this is equivalent to the starting offset of this
 * request
 */
int
ncmpio_last_offset(NC               *ncp,
                   NC_var           *varp,
                   const MPI_Offset  starts[],   /* [varp->ndims] */
                   const MPI_Offset  counts[],   /* [varp->ndims] */
                   const MPI_Offset  strides[],  /* [varp->ndims] */
                   const int         reqMode,
                   MPI_Offset       *offset_ptr) /* return file offset */
{
    MPI_Offset offset, *end_off=NULL;
    int status, i, ndims;

    offset = varp->begin; /* beginning file offset of this variable */
    ndims  = varp->ndims; /* number of dimensions of this variable */

    if (counts != NULL) {
        end_off = (MPI_Offset*) NCI_Malloc((size_t)ndims * SIZEOF_MPI_OFFSET);

        if (strides != NULL) {
            for (i=0; i<ndims; i++)
                end_off[i] = starts[i] + (counts[i] - 1) * strides[i];
        }
        else { /* strides == NULL */
            for (i=0; i<ndims; i++)
                end_off[i] = starts[i] + counts[i] - 1;
        }
    }
    else { /* when counts == NULL strides is of no use */
        end_off = (MPI_Offset*) starts;
    }

    /* check whether end_off is valid */
    status = ncmpii_start_count_stride_check(ncp->format, API_VAR1, varp->ndims,
                                             ncp->numrecs, varp->shape,
                                             end_off, NULL, NULL, reqMode);
    if (status != NC_NOERR) {
#ifdef CDEBUG
        printf("%s(): ncmpii_start_count_stride_check() fails\n",__func__);
#endif
        if (end_off != NULL && end_off != starts) NCI_Free(end_off);
        return status;
    }

    if (ndims > 0) {
        if (IS_RECVAR(varp))
            /* no need to check recsize here: if MPI_Offset is only 32 bits we
               will have had problems long before here */
            offset += end_off[0] * ncp->recsize;
        else
            offset += end_off[ndims-1] * varp->xsz;

        if (ndims > 1) {
            if (IS_RECVAR(varp))
                offset += end_off[ndims - 1] * varp->xsz;
            else
                offset += end_off[0] * varp->dsizes[1] * varp->xsz;

            for (i=1; i<ndims-1; i++)
                offset += end_off[i] * varp->dsizes[i+1] * varp->xsz;
        }
    }
    if (counts != NULL && end_off != NULL)
        NCI_Free(end_off);

    *offset_ptr = offset;
    return NC_NOERR;
}

