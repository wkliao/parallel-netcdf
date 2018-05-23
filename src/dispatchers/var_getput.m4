dnl Process this m4 file to produce 'C' language file.
dnl
dnl If you see this line, you can ignore the next one.
/* Do not edit this file. It is produced from the corresponding .m4 source */
dnl
/*
 *  Copyright (C) 2017, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* $Id$ */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <stdlib.h>

#include <pnetcdf.h>
#include <dispatch.h>
#include <pnc_debug.h>
#include <common.h>

include(`foreach.m4')dnl
include(`utils.m4')dnl
dnl
define(`GOTO_CHECK',`{ DEBUG_ASSIGN_ERROR(err, $1) goto err_check; }')dnl

#define GET_ONE_COUNT(ndims, count) {                                    \
    int _i;                                                              \
    count = (MPI_Offset*) NCI_Malloc((size_t)ndims * SIZEOF_MPI_OFFSET); \
    for (_i=0; _i<ndims; _i++)                                           \
        count[_i] = 1;                                                   \
}

#define GET_FULL_DIMENSIONS(pncp, varp, start, count) {                       \
    int _i;                                                                   \
    start = (MPI_Offset*) NCI_Malloc((size_t)varp.ndims*2*SIZEOF_MPI_OFFSET); \
    count = start + varp.ndims;                                               \
                                                                              \
    for (_i=0; _i<varp.ndims; _i++) {                                         \
        count[_i] = varp.shape[_i];                                           \
        start[_i] = 0;                                                        \
    }                                                                         \
    if (varp.recdim >= 0) { /* find current numrec if varp is record var */   \
        MPI_Offset numrecs;                                                   \
        err = pncp->driver->inq_dim(pncp->ncp, varp.recdim, NULL, &numrecs);  \
        if (err != NC_NOERR) {                                                \
            reqMode |= NC_REQ_ZERO;                                           \
            NCI_Free(start);                                                  \
        }                                                                     \
        else                                                                  \
            count[0] = numrecs;                                               \
    }                                                                         \
}


/*----< check_EINVALCOORDS() >-----------------------------------------------*/
static
int check_EINVALCOORDS(MPI_Offset start,
                       MPI_Offset count,
                       MPI_Offset shape)
{
#ifdef RELAX_COORD_BOUND
    if (start < 0 || start > shape)
        DEBUG_RETURN_ERROR(NC_EINVALCOORDS)
    if (start == shape && count > 0)
        DEBUG_RETURN_ERROR(NC_EINVALCOORDS)
#else
    if (start < 0 || start >= shape)
        DEBUG_RETURN_ERROR(NC_EINVALCOORDS)
#endif
    return NC_NOERR;
}


typedef enum {
    API_GET,
    API_PUT,
    API_IGET,
    API_IPUT,
    API_BPUT
} IO_type;

/*----< check_EEDGE() >------------------------------------------------------*/
static
int check_EEDGE(const MPI_Offset *start,
                const MPI_Offset *count,
                const MPI_Offset *stride,
                const MPI_Offset *shape)
{
    if (*count > *shape || *start + *count > *shape)
        DEBUG_RETURN_ERROR(NC_EEDGE);
    if (stride == NULL) { /* vars APIs but stride is NULL */
        if (*count > *shape || *start + *count > *shape)
            DEBUG_RETURN_ERROR(NC_EEDGE)
    }
    else { /* for vars/varm APIs */
        if (*count > 0 && *start + (*count - 1) * (*stride) >= *shape)
            DEBUG_RETURN_ERROR(NC_EEDGE)
    }
    return NC_NOERR;
}

/*----< check_start_count_stride() >-----------------------------------------*/
static
int check_start_count_stride(PNC              *pncp,
                             int               varid,
                             int               isRead,
                             NC_api            api_kind, /* var1/vara/vars */
                             const MPI_Offset *start,
                             const MPI_Offset *count,
                             const MPI_Offset *stride)
{
    /* only var1, vara, vars, and varm APIs will reach here */
    int i, err, ndims, firstDim;
    MPI_Offset *shape=NULL;

    shape = pncp->vars[varid].shape;
    /* if record variable, obtain the current size of record dimension */
    if (pncp->vars[varid].recdim >= 0) {
        err = pncp->driver->inq_dim(pncp->ncp, pncp->vars[varid].recdim, NULL,
                                    &shape[0]);
        if (err != NC_NOERR) return err;
    }

    /* Check NC_EINVALCOORDS error for argument start[]
     * for API var1/vara/vars/varm, start cannot be NULL, except for scalars
     * and negative start[] is illegal */
    if (start == NULL || start[0] < 0) DEBUG_RETURN_ERROR(NC_EINVALCOORDS)

    firstDim = 0;
    /* check NC_EINVALCOORDS for record dimension */
    if (pncp->vars[varid].recdim >= 0) {
        if (pncp->format < NC_FORMAT_CDF5 && start[0] > NC_MAX_UINT)
            DEBUG_RETURN_ERROR(NC_EINVALCOORDS) /* CDF-1 and 2 */

        /* for record variable, [0] is the NC_UNLIMITED dimension */
        /* read cannot go beyond current numrecs */
        if (isRead) {
            MPI_Offset len = (count == NULL) ? 1 : count[0];
            err = check_EINVALCOORDS(start[0], len, shape[0]);
            if (err != NC_NOERR) return err;
        }
        firstDim = 1; /* done for checking the record dimension */
    }

    /* continue to check NC_EINVALCOORDS for the rest dimensions */
    ndims = pncp->vars[varid].ndims;
    for (i=firstDim; i<ndims; i++) {
        MPI_Offset len = (count == NULL) ? 1 : count[i];
        err = check_EINVALCOORDS(start[i], len, shape[i]);
        if (err != NC_NOERR) return err;
    }

    /* check NC_EEDGE error for argument count[] */

    if (count == NULL) {
        if (api_kind == API_VARA || api_kind == API_VARS ||
            api_kind == API_VARM)
            /* vara/vars/varm, count cannot be NULL */
            DEBUG_RETURN_ERROR(NC_EEDGE)
    }
    else {
        firstDim = 0;
        /* check record dimension */
        if (pncp->vars[varid].recdim >= 0) {
            if (count[0] < 0)  /* no negative count[] */
                DEBUG_RETURN_ERROR(NC_ENEGATIVECNT)

            /* for record variable, [0] is the NC_UNLIMITED dimension */
            /* read cannot go beyond current numrecs */
            if (isRead) {
                err = check_EEDGE(start, count, stride, shape);
                if (err != NC_NOERR) return err;
            }
            firstDim = 1; /* skip checking the record dimension */
        }

        /* continue to check NC_EEDGE for the rest dimensions */
        for (i=firstDim; i<ndims; i++) {
            if (shape[i] < 0) DEBUG_RETURN_ERROR(NC_EEDGE)
            if (count[i] < 0) /* no negative count[] */
                DEBUG_RETURN_ERROR(NC_ENEGATIVECNT)
            if (stride == NULL)
                err = check_EEDGE(start+i, count+i, NULL, shape+i);
            else
                err = check_EEDGE(start+i, count+i, stride+i, shape+i);
            if (err != NC_NOERR) return err;
        }

        /* Check NC_ESTRIDE for non-positive values. We did not check
         * stride[i] >= shape[i], as it is caught as NC_EEDGE error above */
        if (stride != NULL) {
            for (i=0; i<ndims; i++) {
                if (stride[i] <= 0) DEBUG_RETURN_ERROR(NC_ESTRIDE)
            }
        }
    }
    return NC_NOERR;
}

/*----< sanity_check() >-----------------------------------------------------*/
static
int sanity_check(PNC          *pncp,
                 int           varid,
                 IO_type       io,       /* get/put/iget/iput/bput */
                 MPI_Datatype  itype,    /* internal data type */
                 int           isColl)   /* collective or indepdnent API */
{
    /* check file write permission for put APIs */
    if (io == API_PUT || io == API_IPUT || io == API_BPUT)
        if (pncp->flag & NC_MODE_RDONLY) DEBUG_RETURN_ERROR(NC_EPERM)

    /* blocking get/put APIs must be called in data mode */
    if (io == API_PUT || io == API_GET)
        if (pncp->flag & NC_MODE_DEF) DEBUG_RETURN_ERROR(NC_EINDEFINE)

    /* for blocking APIs, check if in collective or independent mode */
    if (io == API_PUT || io == API_GET) {
        if (isColl) { /* check if file is currently in collective data mode */
            if (pncp->flag & NC_MODE_INDEP) DEBUG_RETURN_ERROR(NC_EINDEP)
        }
        else { /* check if file is currently in independent data mode */
            if (!(pncp->flag & NC_MODE_INDEP)) DEBUG_RETURN_ERROR(NC_ENOTINDEP)
        }
    }

    /* variable NC_GLOBAL is illegal in get/put APIs */
    if (varid == NC_GLOBAL) DEBUG_RETURN_ERROR(NC_EGLOBAL)

    /* check whether variable ID is valid */
    if (varid < 0 || varid >= pncp->nvars) DEBUG_RETURN_ERROR(NC_ENOTVAR)

    /* MPI_DATATYPE_NULL in this case represent a flexible API */
    if (itype == MPI_DATATYPE_NULL) return NC_NOERR;

    /* check itype against xtype for NC_ECHAR */
    if (itype == MPI_CHAR) {
        if (pncp->vars[varid].xtype != NC_CHAR) DEBUG_RETURN_ERROR(NC_ECHAR)
    }
    else {
        if (pncp->vars[varid].xtype == NC_CHAR) DEBUG_RETURN_ERROR(NC_ECHAR)
    }
    return NC_NOERR;
}

/*----< allreduce_error() >--------------------------------------------------*/
/* This subroutine is for safe mode to check errors across all processes */
static
int allreduce_error(PNC *pncp, int err)
{
    int minE, mpireturn;
    TRACE_COMM(MPI_Allreduce)(&err, &minE, 1, MPI_INT, MPI_MIN, pncp->comm);
    if (mpireturn != MPI_SUCCESS)
        return ncmpii_error_mpi2nc(mpireturn, "MPI_Allreduce");
    return minE;
}

define(`IO_TYPE', `ifelse(`$1',  `get', `API_GET',
                          `$1',  `put', `API_PUT',
                          `$1', `iget', `API_IGET',
                          `$1', `iput', `API_IPUT',
                          `$1', `bput', `API_BPUT')')dnl
dnl
define(`IS_COLL', `ifelse(`$1',`_all',`1',`0')')dnl
define(`IS_READ', `ifelse(`$1',`get',`1',`$1',`iget',`1',`0')')dnl
define(`IndexArgs', `ifelse(`$1', `',  `NULL, NULL, NULL',
                            `$1', `1', `start, NULL, NULL',
                            `$1', `a', `start, count, NULL',
                            `$1', `s', `start, count, stride',
                            `$1', `m', `start, count, stride')')dnl
dnl
define(`FLEX_ARG',`ifelse(`$1',`',`bufcount, buftype',`-1, ITYPE2MPI($1)')')dnl
dnl
define(`IO_MODE',`ifelse(`$1', `get',`NC_REQ_RD',`$1', `put',`NC_REQ_WR',
                         `$1',`iget',`NC_REQ_RD',`$1',`iput',`NC_REQ_WR',
                         `$1',`bput',`NC_REQ_WR')')dnl
dnl
define(`NB_MODE',`ifelse(`$1', `get',`NC_REQ_BLK',`$1', `put',`NC_REQ_BLK',
                         `$1',`iget',`NC_REQ_NBI',`$1',`iput',`NC_REQ_NBI',
                         `$1',`bput',`NC_REQ_NBB')')dnl
dnl
define(`FLEX_MODE',`ifelse(`$1',`',`NC_REQ_FLEX',`NC_REQ_HL')')dnl
define(`COLL_MODE',`ifelse(`$1',`',`NC_REQ_INDEP',`NC_REQ_COLL')')dnl
dnl
dnl
define(`APINAME',`ifelse(`$3',`',`ncmpi_$1_var$2$4',`ncmpi_$1_var$2_$3$4')')dnl
dnl
dnl GETPUT_API(get/put, `'/1/a/s/m, `'/itype, `'/_all)
dnl
define(`GETPUT_API',dnl
`dnl
/*----< APINAME($1,$2,$3,$4)() >---------------------------------------------*/
/* This API is ifelse(`$4',`',`an independent',`a collective') subroutine. */
int
APINAME($1,$2,$3,$4)(int ncid,
                     int varid,
                     ArgKind($2)
                     BufArgs($1,$3))
{
    int status, err, reqMode=0;
    PNC *pncp;
    ifelse(`$2',`',`',`NC_api api_kind=API_KIND($2);')
    ifelse(`$2',`',`MPI_Offset *start=NULL, *count=NULL;',
           `$2',`1',`MPI_Offset *count=NULL;')

    /* check if ncid is valid.
     * For invalid ncid, we must return error now, as there is no way to
     * continue with invalid ncp. However, collective APIs might hang if this
     * error occurs only on a subset of processes
     */
    err = PNC_check_id(ncid, &pncp);
    if (err != NC_NOERR) return err;

    err = sanity_check(pncp, varid, IO_TYPE($1), ITYPE2MPI($3), IS_COLL($4));

    ifelse(`$2',`m',`if (imap == NULL && stride != NULL) api_kind = API_VARS;
    else if (imap == NULL && stride == NULL) api_kind = API_VARA;',
           `$2',`s',`if (stride == NULL) api_kind = API_VARA;')

    ifelse(`$2',`',`',`/* not-scalar variable checks start, count, stride */
    if (err == NC_NOERR && pncp->vars[varid].ndims > 0)
        err = check_start_count_stride(pncp, varid, IS_READ($1), api_kind,
                                       IndexArgs($2));')

    ifelse(`$4',`',`/* for independent API, return now if error encountered */
    if (err != NC_NOERR) return err;
    ifelse(`$3',`',`
    /* independent flexible API, return now if zero-length request */
    if (bufcount == 0) return NC_NOERR;')',`
    /* collective APIs and safe mode enabled, check errors across all procs */
    if (pncp->flag & NC_MODE_SAFE) {
        err = allreduce_error(pncp, err);
        if (err != NC_NOERR) return err;
    }
    else if (err == NC_EPERM || err == NC_EINDEFINE || err == NC_EINDEP ||
             err == NC_ENOTINDEP) /* cannot continue if fatal errors */
        return err;
    else if (err != NC_NOERR) /* other errors, participate collective call */
        reqMode |= NC_REQ_ZERO;')

    reqMode |= IO_MODE($1) | NB_MODE($1) | FLEX_MODE($3) | COLL_MODE($4);

    ifelse(`$2',`',`if (err == NC_NOERR)
        GET_FULL_DIMENSIONS(pncp, pncp->vars[varid], start, count)',
           `$2',`1',`if (err == NC_NOERR)
        GET_ONE_COUNT(pncp->vars[varid].ndims, count)')

    /* call the subroutine that implements APINAME($1,$2,$3,$4)() */
    status = pncp->driver->`$1'_var(pncp->ncp, varid, start, count,
                                    ArgStrideMap($2), buf,
                                    FLEX_ARG($3), reqMode);

    ifelse(`$2',`',`if (err == NC_NOERR) NCI_Free(start);',
           `$2',`1',`if (err == NC_NOERR) NCI_Free(count);')

    return ifelse(`$4',`',`status;',`(err != NC_NOERR) ? err : status; /* first error encountered */')
}
')dnl
dnl
dnl
foreach(`kind', (, 1, a, s, m),
        `foreach(`putget', (put, get),
                 `foreach(`collindep', (, _all),
                          `foreach(`iType', (`',ITYPE_LIST),
                                   `GETPUT_API(putget,kind,iType,collindep)'
)')')')dnl
dnl
/* ncmpi_get/put_varn_<type>_<mode> API:
 *    type:   data type of I/O buffer, buf
 *    mode:   independent (<nond>) or collective (_all)
 *
 * arguments:
 *    num:    number of start and count pairs
 *    starts: an 2D array of size [num][ndims]. Each starts[i][*] indicates
 *            the starting array indices for a subarray request. ndims is
 *            the number of dimensions of the defined netCDF variable.
 *    counts: an 2D array of size [num][ndims]. Each counts[i][*] indicates
 *            the number of array elements to be accessed. This argument
 *            can be NULL, equivalent to counts with all 1s.
 *    bufcount and buftype: these 2 arguments are only available for flexible
 *            APIs, indicating the I/O buffer memory layout. When buftype is
 *            MPI_DATATYPE_NULL, bufcount is ignored and the data type of buf
 *            is considered matched the variable data type defined in the file.
 */
dnl
define(`NAPINAME',`ifelse(`$2',`',`ncmpi_$1_varn$3',`ncmpi_$1_varn_$2$3')')dnl
dnl
dnl VARN(get/put, `'/iType, `'/_all)
dnl
define(`VARN',dnl
`dnl
/*----< NAPINAME($1,$2,$3)() >-----------------------------------------------*/
/* This API is ifelse(`$3',`',`an independent',`a collective') subroutine. */
int
NAPINAME($1,$2,$3)(int                ncid,
                   int                varid,
                   int                num,
                   MPI_Offset* const *starts,
                   MPI_Offset* const *counts,
                   BufArgs($1,$2))
{
    int i, err, status, reqMode=0;
    PNC *pncp;

    /* check if ncid is valid.
     * For invalid ncid, we must return error now, as there is no way to
     * continue with invalid ncp. However, collective APIs might hang if this
     * error occurs only on a subset of processes
     */
    err = PNC_check_id(ncid, &pncp);
    if (err != NC_NOERR) return err;

    err = sanity_check(pncp, varid, IO_TYPE($1), ITYPE2MPI($2), IS_COLL($3));

    if (num > 0 && starts == NULL) DEBUG_ASSIGN_ERROR(err, NC_ENULLSTART)

    if (err == NC_NOERR && pncp->vars[varid].ndims > 0) {
        NC_api api = (counts == NULL) ? API_VAR1 : API_VARA;
        for (i=0; i<num; i++) {
            MPI_Offset *len = (counts == NULL) ? NULL : counts[i];
            err = check_start_count_stride(pncp, varid, IS_READ($1),
                                           api, starts[i], len, NULL);
            if (err != NC_NOERR) break;
        }
    }

    ifelse(`$3',`',`/* for independent API, return now if error encountered */
    if (err != NC_NOERR) return err;
    /* for independent API, return now if zero-length request */
    if (num == 0) return NC_NOERR;',`
    /* In safe mode, check errors across all processes */
    if (pncp->flag & NC_MODE_SAFE) {
        err = allreduce_error(pncp, err);
        if (err != NC_NOERR) return err;
    }
    else if (err == NC_EPERM || err == NC_EINDEFINE || err == NC_EINDEP ||
             err == NC_ENOTINDEP) /* cannot continue if fatal errors */
        return err;
    else if (err != NC_NOERR) /* other errors, participate collective call */
        reqMode |= NC_REQ_ZERO;')

    reqMode |= IO_MODE($1) | NB_MODE($1) | FLEX_MODE($2) | COLL_MODE($3);

    /* calling the subroutine that implements NAPINAME($1,$2,$3)() */
    status = pncp->driver->`$1'_varn(pncp->ncp, varid, num, starts, counts,
                                     buf, FLEX_ARG($2), reqMode);

    return ifelse(`$3',`',`status;',`(err != NC_NOERR) ? err : status; /* first error encountered */')
}
')dnl
dnl
foreach(`putget', (put, get),
        `foreach(`iType', (`',ITYPE_LIST),
                 `foreach(`collindep', (, _all),
                          `VARN(putget,iType,collindep)'
)')')dnl
dnl
define(`MStartCount',`ifelse(`$1', `',  `NULL, NULL',
                             `$1', `1', `starts[i], NULL',
                             `$1', `a', `starts[i], counts[i]',
                             `$1', `s', `starts[i], counts[i]',
                             `$1', `m', `starts[i], counts[i]')')dnl
dnl
define(`MAPINAME',`ifelse(`$3',`',`ncmpi_m$1_var$2$4',`ncmpi_m$1_var$2_$3$4')')dnl
dnl MVAR(put/get, `'/1/a/s/m, `'/iType, `'/_all)
dnl
define(`MVAR',dnl
`dnl
/*----< MAPINAME($1,$2,$3,$4)() >--------------------------------------------*/
/* This API is ifelse(`$4',`',`an independent',`a collective') subroutine. */
int
MAPINAME($1,$2,$3,$4)(int                ncid,
                      int                nvars,
                      int               *varids,
   ifelse(`$2', `1', `MPI_Offset* const *starts,',
          `$2', `a', `MPI_Offset* const *starts,
                      MPI_Offset* const *counts,',
          `$2', `s', `MPI_Offset* const *starts,
                      MPI_Offset* const *counts,
                      MPI_Offset* const *strides,',
          `$2', `m', `MPI_Offset* const *starts,
                      MPI_Offset* const *counts,
                      MPI_Offset* const *strides,
                      MPI_Offset* const *imaps,')
   ifelse(`$3', `',
    `ifelse($1,`get',`void **bufs,',`void* const *bufs,')
                      const MPI_Offset *bufcounts,
                      const MPI_Datatype *buftypes',
    `ifelse($1,`get',`NC2ITYPE($3) **bufs',
                     `NC2ITYPE($3)* const *bufs')'))
{
    int i, reqMode=0, status=NC_NOERR, err, *reqs;
    PNC *pncp;
    ifelse(`$2',`',`',`NC_api api_kind=API_KIND($2);')

    /* check if ncid is valid.
     * For invalid ncid, we must return error now, as there is no way to
     * continue with invalid ncp. However, collective APIs might hang if this
     * error occurs only on a subset of processes
     */
    err = PNC_check_id(ncid, &pncp);
    if (err != NC_NOERR) return err;

    ifelse(`$4',`',`/* for independent API, return now if zero-length request */
    if (nvars == 0) return NC_NOERR;')

    ifelse(`$2',`m',`if (imaps == NULL && strides != NULL) api_kind = API_VARS;
    else if (imaps == NULL && strides == NULL) api_kind = API_VARA;',
           `$2',`s',`if (strides == NULL) api_kind = API_VARA;')

    for (i=0; i<nvars; i++) {
        err = sanity_check(pncp, varids[i], IO_TYPE($1), ITYPE2MPI($3), IS_COLL($4));
        if (err != NC_NOERR) break;

        ifelse(`$2',`',`',`/* checks start, count, stride for non-scalars */
        if (pncp->vars[varids[i]].ndims > 0) {
            MPI_Offset *stride=NULL;
            ifelse(`$2',`m',`if (strides != NULL) stride = strides[i];',
                   `$2',`s',`if (strides != NULL) stride = strides[i];')
            err = check_start_count_stride(pncp, varids[i], IS_READ($1),
                                           api_kind, MStartCount($2), stride);
            if (err != NC_NOERR) break;
        }')
    }

    reqMode |= IO_MODE($1) | NC_REQ_NBI | FLEX_MODE($3) | COLL_MODE($4);

    ifelse(`$4',`',`/* for independent API, return now if error encountered */
    if (err != NC_NOERR) return err;',`
    /* In safe mode, check errors across all processes */
    if (pncp->flag & NC_MODE_SAFE) {
        err = allreduce_error(pncp, err);
        if (err != NC_NOERR) return err;
    }
    else if (err == NC_EPERM || err == NC_EINDEFINE || err == NC_EINDEP ||
             err == NC_ENOTINDEP) /* cannot continue if fatal errors */
        return err;
    else if (err != NC_NOERR) { /* other errors, participate collective call */
        status = pncp->driver->wait(pncp->ncp, 0, NULL, NULL, reqMode);
        return err;
    }')

    reqs = (int*) NCI_Malloc((size_t)nvars * SIZEOF_INT);
    for (i=0; i<nvars; i++) {
        MPI_Offset *start, *count, *stride=NULL, *imap=NULL;

        /* call the nonblocking subroutines */
        ifelse(`$2',`',`GET_FULL_DIMENSIONS(pncp, pncp->vars[varids[i]], start, count)
        if (err != NC_NOERR) break;',
               `$2',`1',`GET_ONE_COUNT(pncp->vars[varids[i]].ndims, count)
        start = starts[i];',`start = starts[i]; count = counts[i];')
        ifelse(`$2',`s',`if (strides != NULL) stride = strides[i];',
               `$2',`m',`if (strides != NULL) stride = strides[i];
        if (imaps != NULL) imap = imaps[i];')

        err = pncp->driver->i`$1'_var(pncp->ncp, varids[i], start, count,
                                      stride, imap, bufs[i],
                                      ifelse(`$3',`',`bufcounts[i],buftypes[i]',
                                                     `-1, ITYPE2MPI($3)'),
                                      &reqs[i], reqMode);
        ifelse(`$2',`',`NCI_Free(start);',`$2',`1',`NCI_Free(count);')
        if (err != NC_NOERR) break;
    }
    status = pncp->driver->wait(pncp->ncp, i, reqs, NULL, reqMode);
    NCI_Free(reqs);

    return (err != NC_NOERR) ? err : status;
}
')dnl
dnl
foreach(`kind', (, 1, a, s, m),
        `foreach(`putget', (put, get),
                 `foreach(`collindep', (, _all),
                          `foreach(`iType', (`',ITYPE_LIST),
                                   `MVAR(putget,kind,iType,collindep)'
)')')')dnl
dnl
define(`IAPINAME',`ifelse(`$3',`',`ncmpi_$1_var$2',`ncmpi_$1_var$2_$3')')dnl
dnl
dnl IGETPUT_API(iget/iput/bput, `'/1/a/s/m, `'/iType)
dnl
define(`IGETPUT_API',dnl
`dnl
/*----< IAPINAME($1,$2,$3)() >-----------------------------------------------*/
/* This API is an independent subroutine, which can be called in either
 * collective or independent data mode or even in define mode.
 */
int
IAPINAME($1,$2,$3)(int ncid,
                   int varid,
                   ArgKind($2)
                   BufArgs(substr($1,1),$3),
                   int *reqid)
{
    int err, reqMode;
    PNC *pncp;
    ifelse(`$2',`',`',`NC_api api_kind=API_KIND($2);')
    ifelse(`$2',`',`MPI_Offset *start, *count;',`$2',`1',`MPI_Offset *count;')

    /* check if ncid is valid.
     * For invalid ncid, we must return error now, as there is no way to
     * continue with invalid ncp. However, collective APIs might hang if this
     * error occurs only on a subset of processes
     */
    err = PNC_check_id(ncid, &pncp);
    if (err != NC_NOERR) return err;

    if (reqid != NULL) *reqid = NC_REQ_NULL;

    err = sanity_check(pncp, varid, IO_TYPE($1), ITYPE2MPI($3), 0);
    if (err != NC_NOERR) return err;

    ifelse(`$1',`bput',`/* check if buffer has been attached */
    MPI_Offset buf_size;
    err = pncp->driver->inq_misc(pncp->ncp, NULL, NULL, NULL, NULL,
                                 NULL, NULL, NULL, NULL, NULL, NULL,
                                 NULL, NULL, NULL, NULL, &buf_size);
    if (err != NC_NOERR) return err;')

    ifelse(`$2',`m',`if (imap == NULL && stride != NULL) api_kind = API_VARS;
    else if (imap == NULL && stride == NULL) api_kind = API_VARA;',
           `$2',`s',`if (stride == NULL) api_kind = API_VARA;')

    ifelse(`$2',`',`',`/* not-scalar variable checks start, count, stride */
    if (pncp->vars[varid].ndims > 0) {
        err = check_start_count_stride(pncp, varid, IS_READ($1), api_kind,
                                       IndexArgs($2));
        if (err != NC_NOERR) return err;
    }')

    ifelse(`$3',`',`if (bufcount == 0) return NC_NOERR;')

    reqMode = IO_MODE($1) | NB_MODE($1) | FLEX_MODE($3);

    ifelse(`$2',`',`GET_FULL_DIMENSIONS(pncp, pncp->vars[varid], start, count)
                    if (err != NC_NOERR) return err;',
           `$2',`1',`GET_ONE_COUNT(pncp->vars[varid].ndims, count)')

    /* calling the subroutine that implements IAPINAME($1,$2,$3)() */
    err = pncp->driver->`$1'_var(pncp->ncp, varid, start, count,
                                 ArgStrideMap($2), buf,
                                 FLEX_ARG($3), reqid, reqMode);

    ifelse(`$2',`',`NCI_Free(start);',`$2',`1',`NCI_Free(count);')
    return err;
}
')dnl
dnl
foreach(`kind', (, 1, a, s, m),
        `foreach(`putget', (iput, iget, bput),
                 `foreach(`iType', (`',ITYPE_LIST),
                          `IGETPUT_API(putget,kind,iType)'
)')')dnl
dnl
/* ncmpi_iget/iput_varn_<type>_<mode> API:
 *    type:   data type of I/O buffer, buf
 *    mode:   indpendent (<nond>) or collective (_all)
 *
 * arguments:
 *    num:    number of start and count pairs
 *    starts: an 2D array of size [num][ndims]. Each starts[i][*] indicates
 *            the starting array indices for a subarray request. ndims is
 *            the number of dimensions of the defined netCDF variable.
 *    counts: an 2D array of size [num][ndims]. Each counts[i][*] indicates
 *            the number of array elements to be accessed. This argument
 *            can be NULL, equivalent to counts with all 1s.
 *    bufcount and buftype: these 2 arguments are only available for flexible
 *            APIs, indicating the I/O buffer memory layout. When buftype is
 *            MPI_DATATYPE_NULL, bufcount is ignored and the data type of buf
 *            is considered matched the variable data type defined in the file.
 *    reqid:  request ID returned to user
 */
dnl
define(`INAPINAME',`ifelse(`$2',`',`ncmpi_$1_varn',`ncmpi_$1_varn_$2')')dnl
dnl
dnl IVARN(iget/iput/bput, `'/iType)
dnl
define(`IVARN',dnl
`dnl
/*----< INAPINAME($1,$2)() >--------------------------------------------------*/
/* This API is an independent subroutine, which can be called in either
 * collective or independent data mode or even in define mode.
 */
int
INAPINAME($1,$2)(int                ncid,
                 int                varid,
                 int                num,
                 MPI_Offset* const *starts,
                 MPI_Offset* const *counts,
                 BufArgs(substr($1,1),$2),
                 int               *reqid)
{
    int i, err, reqMode;
    PNC *pncp;

    /* check if ncid is valid.
     * For invalid ncid, we must return error now, as there is no way to
     * continue with invalid ncp. However, collective APIs might hang if this
     * error occurs only on a subset of processes
     */
    err = PNC_check_id(ncid, &pncp);
    if (err != NC_NOERR) return err;

    if (reqid != NULL) *reqid = NC_REQ_NULL;

    err = sanity_check(pncp, varid, IO_TYPE($1), ITYPE2MPI($2), 0);
    if (err != NC_NOERR) return err;

    if (num == 0) return NC_NOERR;

    ifelse(`$1',`bput',`/* check if buffer has been attached */
    MPI_Offset buf_size;
    err = pncp->driver->inq_misc(pncp->ncp, NULL, NULL, NULL, NULL,
                                 NULL, NULL, NULL, NULL, NULL, NULL,
                                 NULL, NULL, NULL, NULL, &buf_size);
    if (err != NC_NOERR) return err;')

    if (num > 0 && starts == NULL) DEBUG_RETURN_ERROR(NC_ENULLSTART)

    if (pncp->vars[varid].ndims > 0) {
        NC_api api = (counts == NULL) ? API_VAR1 : API_VARA;
        for (i=0; i<num; i++) {
            MPI_Offset *len = (counts == NULL) ? NULL : counts[i];
            err = check_start_count_stride(pncp, varid, IS_READ($1),
                                           api, starts[i], len, NULL);
            if (err != NC_NOERR) return err;
        }
    }

    ifelse(`$2',`',`if (bufcount == 0) return NC_NOERR;')

    reqMode = IO_MODE($1) | NB_MODE($1) | FLEX_MODE($2);

    /* calling the subroutine that implements INAPINAME($1,$2)() */
    return pncp->driver->`$1'_varn(pncp->ncp, varid, num, starts, counts,
                                   buf, FLEX_ARG($2), reqid, reqMode);
}
')dnl
dnl
foreach(`putget', (iget, iput, bput),
        `foreach(`iType', (`',ITYPE_LIST),
                 `IVARN(putget,iType)'
)')dnl
dnl
dnl VARD(get/put, `'/_all)
dnl
define(`VARD',dnl
`dnl
/*----< ncmpi_$1_vard$2() >--------------------------------------------------*/
/* This API is ifelse(`$4',`',`an independent',`a collective') subroutine. */
int
ncmpi_$1_vard$2(int           ncid,
                int           varid,
                MPI_Datatype  filetype,  /* access layout to the variable in file */
                ifelse($1, `get', `void *buf', `const void *buf'),
                MPI_Offset    bufcount,
                MPI_Datatype  buftype)   /* data type of the buffer */
{
    int err, status, reqMode=0;
    PNC *pncp;

    /* check if ncid is valid.
     * For invalid ncid, we must return error now, as there is no way to
     * continue with invalid ncp. However, collective APIs might hang if this
     * error occurs only on a subset of processes
     */
    err = PNC_check_id(ncid, &pncp);
    if (err != NC_NOERR) return err;

    err = sanity_check(pncp, varid, IO_TYPE($1), MPI_DATATYPE_NULL, IS_COLL($2));

    ifelse(`$2',`',
    `/* for independent API, return now if error encountered or zero request */
    if (err != NC_NOERR) return err;
    if (bufcount == 0) return NC_NOERR;',
    `/* In safe mode, check errors across all processes */
    if (pncp->flag & NC_MODE_SAFE) {
        err = allreduce_error(pncp, err);
        if (err != NC_NOERR) return err;
    }
    else if (err == NC_EPERM || err == NC_EINDEFINE || err == NC_EINDEP ||
             err == NC_ENOTINDEP) /* cannot continue if fatal errors */
        return err;
    else if (err != NC_NOERR) /* other errors, participate collective call */
        reqMode |= NC_REQ_ZERO;')

    reqMode |= IO_MODE($1) | NC_REQ_BLK | NC_REQ_FLEX | COLL_MODE($2);

    /* calling the subroutine that implements ncmpi_$1_vard$2() */
    status = pncp->driver->$1_vard(pncp->ncp, varid, filetype, buf,
                                   bufcount, buftype, reqMode);

    return ifelse(`$2',`',`status;',`(err != NC_NOERR) ? err : status; /* first error encountered */')
}
')
dnl
foreach(`putget', (put, get),
        `foreach(`collindep', (, _all),
                 `VARD(putget,collindep)'
)')dnl
dnl
