dnl Process this m4 file to produce 'C' language file.
dnl
dnl If you see this line, you can ignore the next one.
/* Do not edit this file. It is produced from the corresponding .m4 source */
dnl
/*
 *  Copyright (C) 2014, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* $Id$ */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <stdio.h>
#include <unistd.h>
#ifdef HAVE_STDLIB_H
#include <stdlib.h>
#endif
#include <assert.h>

#include <mpi.h>

#include "nc.h"
#include "ncx.h"
#include "ncmpidtype.h"
#include "macro.h"

/*----< ncmpii_getput_varn() >------------------------------------------------*/
static int
ncmpii_getput_varn(NC                *ncp,
                   NC_var            *varp,
                   int                num,
                   MPI_Offset* const *starts,  /* [num][varp->ndims] */
                   MPI_Offset* const *counts,  /* [num][varp->ndims] */
                   void              *buf,
                   MPI_Offset         bufcount,
                   MPI_Datatype       buftype,   /* data type of the buffer */
                   int                rw_flag,   /* WRITE_REQ or READ_REQ */
                   int                io_method) /* COLL_IO or INDEP_IO */
{
    int i, j, el_size, status=NC_NOERR, min_st, err, free_cbuf=0;
    int req_id=NC_REQ_NULL, st, isSameGroup, position;
    void *cbuf=NULL;
    char *bufp;
    MPI_Offset packsize=0, **_counts=NULL;
    MPI_Datatype ptype;

    /* check for zero-size request */
    if (num == 0 || bufcount == 0) goto err_check;

    /* it is illegal for starts to be NULL */
    if (starts == NULL) {
        DEBUG_ASSIGN_ERROR(status, NC_ENULLSTART)
        goto err_check;
    }
    else { /* it is illegal for any starts[i] to be NULL */
        for (i=0; i<num; i++) {
            if (starts[i] == NULL) {
                DEBUG_ASSIGN_ERROR(status, NC_ENULLSTART)
                goto err_check;
            }
        }
    }

    if (buftype == MPI_DATATYPE_NULL) {
        /* In this case, bufcount is ignored and will be recalculated to match
         * counts[]. Note buf's data type must match the data type of
         * variable defined in the file - no data conversion will be done.
         */
        if (counts == NULL)
            bufcount = 1;
        else {
            bufcount = 0;
            for (j=0; j<num; j++) {
                MPI_Offset bufcount_j = 1;
                if (counts[i] == NULL) {
                    DEBUG_ASSIGN_ERROR(status, NC_ENULLCOUNT)
                    goto err_check;
                }
                for (i=0; i<varp->ndims; i++) {
                    if (counts[j][i] < 0) { /* no negative counts[][] */
                        DEBUG_ASSIGN_ERROR(status, NC_ENEGATIVECNT)
                        goto err_check;
                    }
                    bufcount_j *= counts[j][i];
                }
                bufcount += bufcount_j;
            }
        }
        /* assign buftype match with the variable's data type */
        buftype = ncmpii_nc2mpitype(varp->type);
    }

    cbuf = buf;
    if (bufcount > 0) { /* flexible API is used */
        /* pack buf into cbuf, a contiguous buffer */
        int isderived, iscontig_of_ptypes;
        MPI_Offset bnelems=0;

        /* ptype (primitive MPI data type) from buftype
         * el_size is the element size of ptype
         * bnelems is the total number of ptype elements in buftype
         */
        status = ncmpii_dtype_decode(buftype, &ptype, &el_size, &bnelems,
                                     &isderived, &iscontig_of_ptypes);

        if (status != NC_NOERR) goto err_check;

        if (bufcount != (int)bufcount) {
            DEBUG_ASSIGN_ERROR(status, NC_EINTOVERFLOW)
            goto err_check;
        }

        /* check if buftype is contiguous, if not, pack to one, cbuf */
        if (! iscontig_of_ptypes && bnelems > 0) {
            position = 0;
            packsize  = bnelems*el_size;
            if (packsize != (int)packsize) {
                DEBUG_ASSIGN_ERROR(status, NC_EINTOVERFLOW)
                goto err_check;
            }
            cbuf = NCI_Malloc((size_t)packsize);
            free_cbuf = 1;
            if (rw_flag == WRITE_REQ)
                MPI_Pack(buf, (int)bufcount, buftype, cbuf, (int)packsize,
                         &position, MPI_COMM_SELF);
        }
    }
    else {
        /* this subroutine is called from a high-level API */
        status = NCMPII_ECHAR(varp->type, buftype);
        if (status != NC_NOERR) goto err_check;

        ptype = buftype;
        el_size = ncmpix_len_nctype(varp->type);
    }

    /* We allow counts == NULL and treat this the same as all 1s */
    if (counts == NULL) {
        _counts    = (MPI_Offset**) NCI_Malloc((size_t)num * sizeof(MPI_Offset*));
        _counts[0] = (MPI_Offset*)  NCI_Malloc((size_t)(num * varp->ndims *
                                                        SIZEOF_MPI_OFFSET));
        for (i=1; i<num; i++)
            _counts[i] = _counts[i-1] + varp->ndims;
        for (i=0; i<num; i++)
            for (j=0; j<varp->ndims; j++)
                _counts[i][j] = 1;
    }
    else
        _counts = (MPI_Offset**) counts;

    /* break buf into num pieces */
    isSameGroup=0;
    bufp = (char*)cbuf;
    for (i=0; i<num; i++) {
        MPI_Offset buflen;
        for (buflen=1, j=0; j<varp->ndims; j++) {
            if (_counts[i][j] < 0) { /* any negative counts[][] is illegal */
                DEBUG_ASSIGN_ERROR(status, NC_ENEGATIVECNT)
                goto err_check;
            }
            buflen *= _counts[i][j];
        }
        if (buflen == 0) continue;
        status = ncmpii_igetput_varm(ncp, varp, starts[i], _counts[i], NULL,
                                     NULL, bufp, buflen, ptype, &req_id,
                                     rw_flag, 0, isSameGroup);
        if (status != NC_NOERR) goto err_check;

        /* use isSamegroup so we end up with one nonblocking request (only the
         * first request gets a request ID back, the rest reuse the same ID.
         * This single ID represents num nonblocking requests */
        isSameGroup=1;
        bufp += buflen * el_size;
    }

err_check:
    if (_counts != NULL && _counts != counts) {
        NCI_Free(_counts[0]);
        NCI_Free(_counts);
    }

    if (ncp->safe_mode == 1 && io_method == COLL_IO) {
        int mpireturn;
        TRACE_COMM(MPI_Allreduce)(&status, &min_st, 1, MPI_INT, MPI_MIN,
                                  ncp->nciop->comm);
        if (mpireturn != MPI_SUCCESS)
            return ncmpii_handle_error(mpireturn, "MPI_Allreduce"); 

        if (min_st != NC_NOERR) {
            if (req_id != NC_REQ_NULL) /* cancel pending nonblocking request */
                ncmpii_cancel(ncp, 1, &req_id, &st);
            if (free_cbuf) NCI_Free(cbuf);
            return status;
        }
    }

    if (io_method == INDEP_IO && status != NC_NOERR) {
        if (req_id != NC_REQ_NULL) /* cancel pending nonblocking request */
            ncmpii_cancel(ncp, 1, &req_id, &st);
        if (free_cbuf) NCI_Free(cbuf);
        return status;
    }

    num = 1;
    if (status != NC_NOERR)
        /* This can only be reached for COLL_IO and safe_mode == 0.
           Set num=0 just so this process can participate the collective
           calls in wait_all */
        num = 0;

    err = ncmpiio_wait(ncp, io_method, num, &req_id, &st);

    /* unpack to user buf, if buftype is noncontiguous */
    if (status == NC_NOERR && rw_flag == READ_REQ && free_cbuf) {
        position = 0;
        MPI_Unpack(cbuf, (int)packsize, &position, buf, (int)bufcount, buftype,
                   MPI_COMM_SELF);
    }

    /* return the first error, if there is one */
    if (status == NC_NOERR) status = err;
    if (status == NC_NOERR) status = st;

    if (free_cbuf) NCI_Free(cbuf);

    return status;
}

include(`utils.m4')

dnl
define(`VARN',dnl
`dnl
/*----< ncmpii_$1_varn() >--------------------------------------------------*/
int
ncmpii_$1_varn(void              *ncdp,
               int                varid,
               int                num,
               MPI_Offset* const *starts,
               MPI_Offset* const *counts,
               ifelse(`$1',`put',`const') void *buf,
               MPI_Offset         bufcount,
               MPI_Datatype       buftype,
               nc_type            itype,
               int                io_method)
{
    int     err, status;
    NC     *ncp=(NC*)ncdp;
    NC_var *varp=NULL;

    status = ncmpii_sanity_check(ncp, varid, NULL, NULL, NULL, bufcount,
                                 buftype, API_VARN, (itype==NC_NAT), 1,
                                 ReadWrite($1), io_method, &varp);
    if (status != NC_NOERR) {
        if (io_method == INDEP_IO ||
            status == NC_EBADID ||
            status == NC_EPERM ||
            status == NC_EINDEFINE ||
            status == NC_EINDEP ||
            status == NC_ENOTINDEP)
            return status;  /* fatal error, cannot continue */

        /* for collective API, participate the collective I/O with zero-length
         * request for this process */
        err = ncmpii_getput_zero_req(ncp, ReadWrite($1));
        assert(err == NC_NOERR);

        /* return the error code from sanity check */
        return status;
    }

    return ncmpii_getput_varn(ncp, varp, num, starts, counts, (void*)buf,
                              bufcount, buftype, ReadWrite($1), io_method);
}
')dnl

VARN(put)
VARN(get)

