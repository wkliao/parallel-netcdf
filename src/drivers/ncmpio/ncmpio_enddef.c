/*
 *  Copyright (C) 2003, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* $Id$ */

/*
 * This file implements the corresponding APIs defined in src/dispatchers/file.c
 *
 * ncmpi_enddef()           : dispatcher->enddef()
 * ncmpi__enddef()          : dispatcher->_enddef()
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>  /* strtol() */
#include <string.h>  /* memset() */
#include <assert.h>
#include <errno.h>

#include <mpi.h>

#include <pnc_debug.h>
#include <common.h>
#include "nc.h"
#include "ncx.h"
#ifdef ENABLE_SUBFILING
#include "subfile.h"
#endif

#include "ncio.h"
#include "fbits.h"
#include "rnd.h"    /* D_RNDUP(), M_RNDUP() */


/*----< ncmpiio_move() >-----------------------------------------------------*/
static int
ncmpiio_move(ncio *const nciop,
             MPI_Offset  to,
             MPI_Offset  from,
             MPI_Offset  nbytes)
{
    int rank, nprocs, bufcount, mpireturn, err, status=NC_NOERR, min_st;
    void *buf;
    int chunk_size=1048576; /* move 1 MB per process at a time */
    MPI_Status mpistatus;

    MPI_Comm_size(nciop->comm, &nprocs);
    MPI_Comm_rank(nciop->comm, &rank);

    /* if the file striping unit size is known (obtained from MPI-IO), then
     * we use that instead of 1 MB */
    if (nciop->striping_unit > 0) chunk_size = nciop->striping_unit;

    /* buf will be used as a temporal buffer to move data in chunks, i.e.
     * read a chunk and later write to the new location */
    buf = NCI_Malloc((size_t)chunk_size);
    if (buf == NULL) DEBUG_RETURN_ERROR(NC_ENOMEM)

    /* make fileview entire file visible */
    TRACE_IO(MPI_File_set_view)(nciop->collective_fh, 0, MPI_BYTE, MPI_BYTE,
                                "native", MPI_INFO_NULL);

    /* move the variable starting from its tail toward its beginning */
    while (nbytes > 0) {
        int get_size=0;

        /* calculate how much to move at each time */
        bufcount = chunk_size;
        if (nbytes < (MPI_Offset)nprocs * chunk_size) {
            /* handle the last group of chunks */
            MPI_Offset rem_chunks = nbytes / chunk_size;
            if (rank > rem_chunks) /* these processes do not read/write */
                bufcount = 0;
            else if (rank == rem_chunks) /* this process reads/writes less */
                bufcount = (int)(nbytes % chunk_size);
            nbytes = 0;
        }
        else {
            nbytes -= chunk_size*nprocs;
        }

        /* explicitly initialize mpistatus object to 0, see comments below */
        memset(&mpistatus, 0, sizeof(MPI_Status));

        /* read the original data @ from+nbytes+rank*chunk_size */
        TRACE_IO(MPI_File_read_at_all)(nciop->collective_fh,
                                       from+nbytes+rank*chunk_size,
                                       buf, bufcount, MPI_BYTE, &mpistatus);
        if (mpireturn != MPI_SUCCESS) {
            err = ncmpii_handle_error(mpireturn, "MPI_File_read_at_all");
            if (err == NC_EFILE) DEBUG_ASSIGN_ERROR(status, NC_EREAD)
        }
        else {
            /* for zero-length read, MPI_Get_count may report incorrect result
             * for some MPICH version, due to the uninitialized MPI_Status
             * object passed to MPI-IO calls. Thus we initialize it above to
             * work around. Otherwise we can just use:
            nciop->get_size += bufcount;
             */
            MPI_Get_count(&mpistatus, MPI_BYTE, &get_size);
            nciop->get_size += get_size;
        }

        /* MPI_Barrier(nciop->comm); */
        /* important, in case new region overlaps old region */
        TRACE_COMM(MPI_Allreduce)(&status, &min_st, 1, MPI_INT, MPI_MIN, nciop->comm);
        status = min_st;
        if (status != NC_NOERR) break;

        /* write to new location @ to+nbytes+rank*chunk_size
         *
         * Ideally, we should write the amount of get_size returned from a call
         * to MPI_Get_count in the below MPI write. This is in case some
         * variables are defined but never been written. The value returned by
         * MPI_Get_count is supposed to be the actual amount read by the MPI
         * read call. If partial data (or none) is available for read, then we
         * should just write that amount. Note this MPI write is collective,
         * and thus all processes must participate the call even if get_size
         * is 0. However, in some MPICH versions MPI_Get_count fails to report
         * the correct value due to an internal error that fails to initialize
         * the MPI_Status object. Therefore, the solution can be either to
         * explicitly initialize the status object to zeros, or to just use
         * bufcount for write. Note that the latter will write the variables
         * that have not been written before. Below uses the former option.
         */
        TRACE_IO(MPI_File_write_at_all)(nciop->collective_fh,
                                        to+nbytes+rank*chunk_size,
                                        buf, get_size /* bufcount */,
                                        MPI_BYTE, &mpistatus);
        if (mpireturn != MPI_SUCCESS) {
            err = ncmpii_handle_error(mpireturn, "MPI_File_write_at_all");
            if (err == NC_EFILE) DEBUG_ASSIGN_ERROR(status, NC_EWRITE)
        }
        else {
            /* for zero-length read, MPI_Get_count may report incorrect result
             * for some MPICH version, due to the uninitialized MPI_Status
             * object passed to MPI-IO calls. Thus we initialize it above to
             * work around. Otherwise we can just use:
            nciop->put_size += bufcount;
             */
            int put_size;
            MPI_Get_count(&mpistatus, MPI_BYTE, &put_size);
            nciop->put_size += put_size;
        }
        TRACE_COMM(MPI_Allreduce)(&status, &min_st, 1, MPI_INT, MPI_MIN, nciop->comm);
        status = min_st;
        if (status != NC_NOERR) break;
    }
    NCI_Free(buf);
    return status;
}

/*----< move_fixed_vars() >--------------------------------------------------*/
/* move one fixed variable at a time, only when the new begin > old begin */
static int
move_fixed_vars(NC *ncp, NC *old)
{
    int i, err, status=NC_NOERR;

    /* move starting from the last fixed variable */
    for (i=old->vars.ndefined-1; i>=0; i--) {
        if (IS_RECVAR(old->vars.value[i])) continue;

        MPI_Offset from = old->vars.value[i]->begin;
        MPI_Offset to   = ncp->vars.value[i]->begin;
        if (to > from) {
            err = ncmpiio_move(ncp->nciop, to, from, ncp->vars.value[i]->len);
            if (status == NC_NOERR) status = err;
        }
    }
    return status;
}

/*----< move_recs_r() >------------------------------------------------------*/
/*
 * Move the record variables down,
 * re-arrange records as needed
 * Fill as needed.
 */
static int
move_recs_r(NC *ncp, NC *old) {
    int status;
    MPI_Offset recno;
    const MPI_Offset nrecs = ncp->numrecs;
    const MPI_Offset ncp_recsize = ncp->recsize;
    const MPI_Offset old_recsize = old->recsize;
    const off_t ncp_off = ncp->begin_rec;
    const off_t old_off = old->begin_rec;

    assert(ncp_recsize >= old_recsize);

    if (ncp_recsize == old_recsize) {
        if (ncp_recsize == 0) /* no record variable defined yet */
            return NC_NOERR;

        /* No new record variable inserted, move all record variables as a whole */
        status = ncmpiio_move(ncp->nciop, ncp_off, old_off, ncp_recsize * nrecs);
        if (status != NC_NOERR)
            return status;
    } else {
        /* new record variables inserted, move one whole record at a time */
        for (recno = nrecs-1; recno >= 0; recno--) {
            status = ncmpiio_move(ncp->nciop,
                                  ncp_off+recno*ncp_recsize,
                                  old_off+recno*old_recsize,
                                  old_recsize);
            if (status != NC_NOERR)
                return status;
        }
    }

    return NC_NOERR;
}

/*----< NC_begins() >--------------------------------------------------------*/
/*
 * This function is only called at enddef().
 * It computes each variable's 'begin' offset, and sets/updates the followings:
 *    ncp->xsz                   ---- header size
 *    ncp->vars.value[*]->begin  ---- each variable's 'begin' offset
 *    ncp->begin_var             ---- offset of first non-record variable
 *    ncp->begin_rec             ---- offset of first     record variable
 *    ncp->recsize               ---- sum of single records
 *    ncp->numrecs               ---- number of records (set only if new file)
 */
static int
NC_begins(NC *ncp)
{
    int i, j, rank, mpireturn;
    MPI_Offset end_var=0;
    NC_var *last = NULL;
    NC_var *first_var = NULL;       /* first "non-record" var */

    /* CDF file format determines the size of variable's "begin" in the header */

    /* get the true header size (not header extent) */
    MPI_Comm_rank(ncp->nciop->comm, &rank);
    ncp->xsz = ncmpii_hdr_len_NC(ncp);

    if (ncp->safe_mode) { /* this consistency check is redundant as metadata is
                             kept consistent at all time when safe mode is on */
        int err, status;
        MPI_Offset root_xsz = ncp->xsz;

        /* only root's header size matters */
        TRACE_COMM(MPI_Bcast)(&root_xsz, 1, MPI_OFFSET, 0, ncp->nciop->comm);
        if (mpireturn != MPI_SUCCESS)
            return ncmpii_handle_error(mpireturn, "MPI_Bcast"); 

        err = NC_NOERR;
        if (root_xsz != ncp->xsz) err = NC_EMULTIDEFINE;

        /* find min error code across processes */
        TRACE_COMM(MPI_Allreduce)(&err, &status, 1, MPI_INT, MPI_MIN, ncp->nciop->comm);
        if (mpireturn != MPI_SUCCESS)
            return ncmpii_handle_error(mpireturn, "MPI_Allreduce");
        if (status != NC_NOERR) DEBUG_RETURN_ERROR(status)
    }

    /* This function is called in ncmpi_enddef(), which can happen either when
     * creating a new file or opening an existing file with metadata modified.
     * For the former case, ncp->begin_var == 0 here.
     * For the latter case, we set begin_var a new value only if the new header
     * grows out of its extent or the start of non-record variables is not
     * aligned as requested by ncp->h_align.
     * Note ncp->xsz is header size and ncp->begin_var is header extent.
     * Add the minimum header free space requested by user.
     */
    ncp->begin_var = D_RNDUP(ncp->xsz + ncp->h_minfree, ncp->h_align);

    if (ncp->old != NULL) {
        /* If this define mode was entered from a redef(), we check whether
         * the new begin_var against the old begin_var. We do not shrink
         * the header extent.
         */
        if (ncp->begin_var < ncp->old->begin_var)
            ncp->begin_var = ncp->old->begin_var;
    }

    /* ncp->begin_var is the aligned starting file offset of the first
       variable, also the extent of file header */

    /* Now calculate the starting file offsets for all variables.
       loop thru vars, first pass is for the 'non-record' vars */
    end_var = ncp->begin_var;
    for (j=0, i=0; i<ncp->vars.ndefined; i++) {
        if (IS_RECVAR(ncp->vars.value[i]))
            /* skip record variables on this pass */
            continue;
        if (first_var == NULL) first_var = ncp->vars.value[i];

        /* for CDF-1 check if over the file size limit 32-bit integer */
        if (ncp->format == 1 && end_var > X_OFF_MAX)
            DEBUG_RETURN_ERROR(NC_EVARSIZE)

        /* this will pad out non-record variables with zero to the
         * requested alignment.  record variables are a bit trickier.
         * we don't do anything special with them */
        ncp->vars.value[i]->begin = D_RNDUP(end_var, ncp->v_align);

        if (ncp->old != NULL) {
            /* move to the next fixed variable */
            for (; j<ncp->old->vars.ndefined; j++)
                if (!IS_RECVAR(ncp->old->vars.value[j]))
                    break;
            if (j < ncp->old->vars.ndefined) {
                if (ncp->vars.value[i]->begin < ncp->old->vars.value[j]->begin)
                    /* the first ncp->vars.ndefined non-record variables should
                       be the same. If the new begin is smaller, reuse the old
                       begin */
                    ncp->vars.value[i]->begin = ncp->old->vars.value[j]->begin;
                j++;
            }
        }
        /* end_var is the end offset of variable i */
        end_var = ncp->vars.value[i]->begin + ncp->vars.value[i]->len;
    }

    /* end_var now is pointing to the end of last non-record variable */

    /* only (re)calculate begin_rec if there is not sufficient
     * space at end of non-record variables or if start of record
     * variables is not aligned as requested by ncp->r_align.
     * If the existing begin_rec is already >= index, then leave the
     * begin_rec as is (in case some non-record variables are deleted)
     */
    if (ncp->begin_rec < end_var ||
        ncp->begin_rec != D_RNDUP(ncp->begin_rec, ncp->v_align))
        ncp->begin_rec = D_RNDUP(end_var, ncp->v_align);

    /* expand free space for fixed variable section */
    if (ncp->begin_rec < end_var + ncp->v_minfree)
        ncp->begin_rec = D_RNDUP(end_var + ncp->v_minfree, ncp->v_align);

    /* align the starting offset for record variable section */
    if (ncp->r_align > 1)
        ncp->begin_rec = D_RNDUP(ncp->begin_rec, ncp->r_align);

    if (ncp->old != NULL) {
        /* check whether the new begin_rec is smaller */
        if (ncp->begin_rec < ncp->old->begin_rec)
            ncp->begin_rec = ncp->old->begin_rec;
    }

    if (first_var != NULL)
        ncp->begin_var = first_var->begin;
    else
        ncp->begin_var = ncp->begin_rec;

    end_var = ncp->begin_rec;
    /* end_var now is pointing to the beginning of record variables
     * note that this can be larger than the end of last non-record variable
     */

    ncp->recsize = 0;

    /* TODO: alignment for record variables (maybe using a new hint) */

    /* loop thru vars, second pass is for the 'record' vars,
     * re-calculate the starting offset for each record variable */
    for (j=0, i=0; i<ncp->vars.ndefined; i++) {
        if (!IS_RECVAR(ncp->vars.value[i]))
            /* skip non-record variables on this pass */
            continue;

        /* X_OFF_MAX is the max of 32-bit integer */
        if (ncp->format == 1 && end_var > X_OFF_MAX)
            DEBUG_RETURN_ERROR(NC_EVARSIZE)

        /* A few attempts at aligning record variables have failed
         * (either with range error or 'value read not that expected',
         * or with an error in ncmpi_redef()).  Not sufficient to align
         * 'begin', but haven't figured out what else to adjust */
        ncp->vars.value[i]->begin = end_var;

        if (ncp->old != NULL) {
            /* move to the next record variable */
            for (; j<ncp->old->vars.ndefined; j++)
                if (IS_RECVAR(ncp->old->vars.value[j]))
                    break;
            if (j < ncp->old->vars.ndefined) {
                if (ncp->vars.value[i]->begin < ncp->old->vars.value[j]->begin)
                    /* if the new begin is smaller, use the old begin */
                    ncp->vars.value[i]->begin = ncp->old->vars.value[j]->begin;
                j++;
            }
        }
        end_var += ncp->vars.value[i]->len;
        /* end_var is the end offset of record variable i */

        /* check if record size must fit in 32-bits (for CDF-1) */
#if SIZEOF_OFF_T == SIZEOF_SIZE_T && SIZEOF_SIZE_T == 4
        if (ncp->recsize > X_UINT_MAX - ncp->vars.value[i]->len)
            DEBUG_RETURN_ERROR(NC_EVARSIZE)
#endif
        ncp->recsize += ncp->vars.value[i]->len;
        last = ncp->vars.value[i];
    }

    /*
     * for special case (Check CDF-1 and CDF-2 file format specifications.)
     * "A special case: Where there is exactly one record variable, we drop the
     * requirement that each record be four-byte aligned, so in this case there
     * is no record padding."
     */
    if (last != NULL) {
        if (ncp->recsize == last->len) {
            /* exactly one record variable, pack value */
            ncp->recsize = *last->dsizes * last->xsz;
        }
#if 0
        else if (last->len == UINT32_MAX) { /* huge last record variable */
            ncp->recsize += *last->dsizes * last->xsz;
        }
#endif
    }

/* below is only needed if alignment is performed on record variables */
#if 0
    /*
     * for special case of exactly one record variable, pack value
     */
    /* if there is exactly one record variable, then there is no need to
     * pad for alignment -- there's nothing after it */
    if (last != NULL && ncp->recsize == last->len)
        ncp->recsize = *last->dsizes * last->xsz;
#endif

    if (NC_IsNew(ncp))
        NC_set_numrecs(ncp, 0);

    return NC_NOERR;
}

/*----< write_NC() >---------------------------------------------------------*/
/*
 * This function is collective and only called by enddef().
 * Write out the header
 * 1. Call ncmpii_hdr_put_NC() to copy the header object, ncp, to a buffer.
 * 2. Call NC_check_header() to check if header is consistent across all
 *    processes.
 * 3. Process rank 0 writes the header to file.
 * This is a collective call.
 */
static int
write_NC(NC *ncp)
{
    void *buf;
    int status, mpireturn, err, rank;
    MPI_Offset local_xsz;

    assert(!NC_readonly(ncp));

    /* In NC_begins(), root's ncp->xsz, root's header size, has been
     * broadcasted, so ncp->xsz is now root's header size. To check any
     * inconsistency in file header, we need to calculate local's header
     * size by calling ncmpii_hdr_len_NC()./
     */
    local_xsz = ncmpii_hdr_len_NC(ncp);

    /* Note valgrind will complain about uninitialized buf below, but buf will
     * be first filled with header of size ncp->xsz and later write to file.
     * So, no need to change to NCI_Calloc for the sake of valgrind.
     */
    buf = NCI_Malloc((size_t)local_xsz); /* buffer for local header object */

    /* copy the entire local header object to buffer */
    status = ncmpii_hdr_put_NC(ncp, buf);
    if (status != NC_NOERR) { /* a fatal error */
        NCI_Free(buf);
        return status;
    }

    MPI_Comm_rank(ncp->nciop->comm, &rank);

#ifdef _DIFF_HEADER
    if (ncp->safe_mode == 0) {
        int h_size=(int)ncp->xsz;
        void *root_header;

        /* check header against root's */
        if (rank == 0)
            root_header = buf;
        else
            root_header = (void*) NCI_Malloc((size_t)h_size);

        TRACE_COMM(MPI_Bcast)(root_header, h_size, MPI_BYTE, 0, ncp->nciop->comm);
        if (mpireturn != MPI_SUCCESS)
            return ncmpii_handle_error(mpireturn, "MPI_Bcast"); 

        if (rank > 0) {
            if (h_size != local_xsz || memcmp(buf, root_buf, h_size))
                status = NC_EMULTIDEFINE;
            NCI_Free(root_buf);
        }

        /* report error if header is inconsistency */
        TRACE_COMM(MPI_Allreduce)(&status, &err, 1, MPI_INT, MPI_MIN, ncp->nciop->comm);
        if (mpireturn != MPI_SUCCESS) {
            ncmpii_handle_error(mpireturn,"MPI_Allreduce");
            DEBUG_RETURN_ERROR(NC_EMPI)
        }
        if (err != NC_NOERR) {
            /* TODO: this error return is harsh. Maybe relax for inconsistent
             * attribute contents? */
            if (status == NC_NOERR) status = err;
            NCI_Free(buf);
            return status;
        }
    }
#endif

#ifdef _CHECK_HEADER_IN_DETAIL
    /* check the header consistency across all processes and sync header.
     * When safe_mode is on:
     *   The returned status on root can be either NC_NOERR (all headers are
     *   consistent) or NC_EMULTIDEFINE (some headers are inconsistent).
     *   The returned status on non-root processes can be NC_NOERR, fatal
     *   error (>-250), or inconsistency error (-250 to -269).
     * When safe_mode is off:
     *   The returned status on root is always NC_NOERR
     *   The returned status on non-root processes can be NC_NOERR, fatal
     *   error (>-250), or inconsistency error (-250 to -269).
     * For fatal error, we should stop. For others, we can continue.
     */
    int max_err;
    status = NC_check_header(ncp, buf, local_xsz);

    /* check for fatal error */
    err =  (status != NC_NOERR && !ErrIsHeaderDiff(status)) ? 1 : 0;
    max_err = err;

    if (ncp->safe_mode == 1) {
        TRACE_COMM(MPI_Allreduce)(&err, &max_err, 1, MPI_INT, MPI_MAX,
                                  ncp->nciop->comm);
        if (mpireturn != MPI_SUCCESS) {
            ncmpii_handle_error(mpireturn,"MPI_Allreduce");
            DEBUG_RETURN_ERROR(NC_EMPI)
        }
    }
    if (max_err == 1) { /* some processes encounter a fatal error */
        NCI_Free(buf);
        return status;
    }
#endif
    /* For non-fatal error, we continue to write header to the file, as now the
     * header object in memory has been sync-ed across all processes. */

    /* only rank 0's header gets written to the file */
    if (rank == 0) {
        /* rank 0's fileview already includes the file header */
        MPI_Status mpistatus;
        if (ncp->xsz != (int)ncp->xsz)
            DEBUG_ASSIGN_ERROR(status, NC_EINTOVERFLOW)

        TRACE_IO(MPI_File_write_at)(ncp->nciop->collective_fh, 0, buf,
                                    (int)ncp->xsz, MPI_BYTE, &mpistatus);
        if (mpireturn != MPI_SUCCESS) {
            err = ncmpii_handle_error(mpireturn, "MPI_File_write_at");
            /* write has failed, which is more serious than inconsistency */
            if (err == NC_EFILE) DEBUG_ASSIGN_ERROR(status, NC_EWRITE)
        }
        else {
            ncp->nciop->put_size += ncp->xsz;
        }
    }

    if (ncp->safe_mode == 1) {
        /* broadcast root's status, because only root writes to the file */
        int root_status = status;
        TRACE_COMM(MPI_Bcast)(&root_status, 1, MPI_INT, 0, ncp->nciop->comm);
        /* root's write has failed, which is more serious than inconsistency */
        if (root_status == NC_EWRITE) DEBUG_ASSIGN_ERROR(status, NC_EWRITE)
    }

    fClr(ncp->flags, NC_NDIRTY);
    NCI_Free(buf);

    return status;
}

/* Many subroutines called in ncmpii__enddef() are collective. We check the
 * error codes of all processes only in safe mode, so the program can stop
 * collectively, if any one process got an error. However, when safe mode is
 * off, we simply return the error and program may hang if some processes
 * do not get error and proceed to the next subroutine call.
 */ 
#define CHECK_ERROR(err) {                                              \
    if (ncp->safe_mode == 1) {                                          \
        int status;                                                     \
        TRACE_COMM(MPI_Allreduce)(&err, &status, 1, MPI_INT, MPI_MIN,   \
                                  ncp->nciop->comm);                    \
        if (mpireturn != MPI_SUCCESS)                                   \
            return ncmpii_handle_error(mpireturn, "MPI_Allreduce");     \
        if (status != NC_NOERR) return status;                          \
    }                                                                   \
    else if (err != NC_NOERR)                                           \
        return err;                                                     \
}

/*----< ncmpii__enddef() >---------------------------------------------------*/
/* This is a collective subroutine. */
int
ncmpii__enddef(void       *ncdp,
               MPI_Offset  h_minfree,
               MPI_Offset  v_align,
               MPI_Offset  v_minfree,
               MPI_Offset  r_align)
{
    int i, flag, striping_unit, mpireturn, err=NC_NOERR, status=NC_NOERR;
    char value[MPI_MAX_INFO_VAL];
    MPI_Offset h_align, all_var_size;
    NC *ncp = (NC*)ncdp;
#ifdef ENABLE_SUBFILING
    NC *ncp_sf=NULL;
#endif

    if (!NC_indef(ncp)) /* must currently in define mode */
        DEBUG_RETURN_ERROR(NC_ENOTINDEFINE)

    if (NC_readonly(ncp)) /* must have write permission */
        DEBUG_RETURN_ERROR(NC_EPERM)

    if (ncp->safe_mode) {
        /* check if h_minfree, v_align, v_minfree, and r_align are consistent
         * among all processes */
        MPI_Offset root_args[4];

        root_args[0] = h_minfree;
        root_args[1] = v_align;
        root_args[2] = v_minfree;
        root_args[3] = r_align;
        TRACE_COMM(MPI_Bcast)(&root_args, 4, MPI_OFFSET, 0, ncp->nciop->comm);
        if (mpireturn != MPI_SUCCESS)
            return ncmpii_handle_error(mpireturn, "MPI_Bcast");

        if (root_args[0] != h_minfree ||
            root_args[1] != v_align   ||
            root_args[2] != v_minfree ||
            root_args[3] != r_align)
            DEBUG_ASSIGN_ERROR(err, NC_EMULTIDEFINE_FNC_ARGS)
        else
            err = NC_NOERR;

        /* find min error code across processes */
        TRACE_COMM(MPI_Allreduce)(&err, &status, 1, MPI_INT, MPI_MIN, ncp->nciop->comm);
        if (mpireturn != MPI_SUCCESS)
            return ncmpii_handle_error(mpireturn, "MPI_Allreduce");

        if (status != NC_NOERR) return status;
    }

    if (h_minfree > 0) ncp->h_minfree = h_minfree;
    if (v_minfree > 0) ncp->v_minfree = v_minfree;

    /* calculate a good align size for PnetCDF level hints:
     * header_align_size and var_align_size based on the MPI-IO hint
     * striping_unit. This hint can be either supplied by the user or obtained
     * from MPI-IO (for example, ROMIO's Lustre driver makes a system call to
     * get the striping parameters of a file).
     */
    MPI_Info_get(ncp->nciop->mpiinfo, "striping_unit", MPI_MAX_INFO_VAL-1,
                 value, &flag);
    striping_unit = 0;
    if (flag) {
        errno = 0;
        striping_unit = (int)strtol(value,NULL,10);
        if (errno != 0) striping_unit = 0;
    }
    ncp->nciop->striping_unit = striping_unit;

    all_var_size = 0;  /* sum of all defined variables */
    for (i=0; i<ncp->vars.ndefined; i++)
        all_var_size += ncp->vars.value[i]->len;

    /* ncp->h_align, ncp->v_align, ncp->r_align, and ncp->chunk have been
     * set during file create/open */

    if (ncp->h_align == 0) { /* user info does not set nc_header_align_size */
        if (striping_unit &&
            all_var_size > HEADER_ALIGNMENT_LB * striping_unit)
            /* if striping_unit is available and file size sufficiently large */
            ncp->h_align = striping_unit;
        else
            ncp->h_align = FILE_ALIGNMENT_DEFAULT;
    }
    /* else respect user hint */

    if (ncp->v_align == 0) { /* user info does not set nc_var_align_size */
        if (v_align > 0) /* else respect user hint */
            ncp->v_align = v_align;
        else if (striping_unit &&
                 all_var_size > HEADER_ALIGNMENT_LB * striping_unit)
            /* if striping_unit is available and file size sufficiently large */
            ncp->v_align = striping_unit;
        else
            ncp->v_align = FILE_ALIGNMENT_DEFAULT;
    }

    if (ncp->r_align == 0) { /* user info does not set nc_record_align_size */
        if (r_align > 0) /* else respect user hint */
            ncp->r_align = r_align;
        if (striping_unit)
            ncp->r_align = striping_unit;
        else
            ncp->r_align = FILE_ALIGNMENT_DEFAULT;
    }

    /* all CDF formats require 4-bytes alignment */
    if (ncp->h_align == 0) ncp->h_align = 4;
    else                   ncp->h_align = D_RNDUP(ncp->h_align, 4);
    if (ncp->v_align == 0) ncp->v_align = 4;
    else                   ncp->v_align = D_RNDUP(ncp->v_align, 4);
    if (ncp->r_align == 0) ncp->r_align = 4;
    else                   ncp->r_align = D_RNDUP(ncp->r_align, 4);

    /* reflect the hint changes to the MPI info object, so the user can inquire
     * what the true hint values are being used
     */
    sprintf(value, "%lld", ncp->h_align);
    MPI_Info_set(ncp->nciop->mpiinfo, "nc_header_align_size", value);
    sprintf(value, "%lld", ncp->v_align);
    MPI_Info_set(ncp->nciop->mpiinfo, "nc_var_align_size", value);
    sprintf(value, "%lld", ncp->r_align);
    MPI_Info_set(ncp->nciop->mpiinfo, "nc_record_align_size", value);

#ifdef ENABLE_SUBFILING
    /* num of subfiles has been determined already */
    if (ncp->nc_num_subfiles > 1) {
        /* TODO: should return subfile-related msg when there's an error */
        err = ncmpii_subfile_partition(ncp, &ncp->ncid_sf);
        CHECK_ERROR(err)
    }
#endif

    /* check whether sizes of all variables are legal */
    err = ncmpii_NC_check_vlens(ncp);
    CHECK_ERROR(err)

    /* When ncp->old == NULL, this enddef is called the first time after file
     * create call. In this case, we compute each variable's 'begin', starting
     * file offset as well as the offsets of record variables.
     * When ncp->old != NULL, this enddef is called after a redef. In this
     * case, we re-used all variable offsets as many as possible.
     *
     * Note in NC_begins, root broadcasts ncp->xsz, the file header size, to
     * all processes.
     */
    err = NC_begins(ncp);
    CHECK_ERROR(err)

#ifdef ENABLE_SUBFILING
    if (ncp->nc_num_subfiles > 1) {
        /* get ncp info for the subfile */
        PNC *pncp;
        err = PNC_check_id(ncp->ncid_sf, &pncp);
        if (err != NC_NOERR) DEBUG_RETURN_ERROR(err)

        ncp_sf = (NC*)pncp->ncp;
        err = NC_begins(ncp_sf);
        CHECK_ERROR(err)
    }
#endif

    if (ncp->old != NULL) {
        /* The current define mode was entered from ncmpi_redef, not from
         * ncmpi_create. We must check if header has been expanded.
         */

        assert(!NC_IsNew(ncp));
        assert(fIsSet(ncp->flags, NC_INDEF));
        assert(ncp->begin_rec >= ncp->old->begin_rec);
        assert(ncp->begin_var >= ncp->old->begin_var);
        assert(ncp->vars.ndefined >= ncp->old->vars.ndefined);
        /* ncp->numrecs has already sync-ed in ncmpi_redef */

        if (ncp->vars.ndefined > 0) { /* no. record and non-record variables */
            if (ncp->begin_var > ncp->old->begin_var) {
                /* header size increases, shift the entire data part down */
                /* shift record variables first */
                err = move_recs_r(ncp, ncp->old);
                CHECK_ERROR(err)

                /* shift non-record variables */
                /* err = move_vars_r(ncp, ncp->old); */
                err = move_fixed_vars(ncp, ncp->old);
                CHECK_ERROR(err)
            }
            else if (ncp->begin_rec > ncp->old->begin_rec ||
                     ncp->recsize   > ncp->old->recsize) {
                /* number of non-record variables increases, or
                   number of records of record variables increases,
                   shift and move all record variables down */
                err = move_recs_r(ncp, ncp->old);
                CHECK_ERROR(err)
            }
        }
    } /* ... ncp->old != NULL */

    /* first sync header objects in memory across all processes, and then root
     * writes the header to file. Note safe_mode error check will be done in
     * write_NC() */
    status = write_NC(ncp);

    /* we should continue to exit define mode, even if header is inconsistent
     * among processes, so the program can proceed, say to close file properly.
     * However, if ErrIsHeaderDiff(status) is true, this error should
     * be considered fatal, as inconsistency is about the data structure,
     * rather then contents (such as attribute values) */

#ifdef ENABLE_SUBFILING
    /* write header to subfile */
    if (ncp->nc_num_subfiles > 1) {
        err = write_NC(ncp_sf);
        if (status == NC_NOERR) status = err;
    }
#endif

    /* update the total number of record variables */
    ncp->vars.num_rec_vars = 0;
    for (i=0; i<ncp->vars.ndefined; i++)
        ncp->vars.num_rec_vars += IS_RECVAR(ncp->vars.value[i]);

    /* fill variables according to their fill mode settings */
    if (ncp->vars.ndefined > 0) {
        err = ncmpii_fill_vars(ncp);
        if (status == NC_NOERR) status = err;
    }

    if (ncp->old != NULL) {
        ncmpii_free_NC(ncp->old);
        ncp->old = NULL;
    }
    fClr(ncp->flags, NC_CREAT | NC_INDEF);

#ifdef ENABLE_SUBFILING
    if (ncp->nc_num_subfiles > 1)
        fClr(ncp_sf->flags, NC_CREAT | NC_INDEF);
#endif

    /* If the user sets NC_SHARE, we enforce a stronger data consistency */
    if (NC_doFsync(ncp))
        ncmpiio_sync(ncp->nciop); /* calling MPI_File_sync() */

    return status;
}

/*----< ncmpii_enddef() >----------------------------------------------------*/
/* This is a collective subroutine. */
int
ncmpii_enddef(void *ncdp)
{
    return ncmpii__enddef(ncdp, 0, 0, 0, 0);
}

