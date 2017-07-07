/*
 *  Copyright (C) 2003, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* $Id$ */

/*
 * This file implements the corresponding APIs defined in src/dispatchers/file.c
 *
 * ncmpi_redef()            : dispatcher->redef()
 * ncmpi_close()            : dispatcher->close()
 * ncmpi_sync()             : dispatcher->sync()
 * ncmpi_abort()            : dispatcher->abort()
 * ncmpi_inq()              : dispatcher->inq()
 * ncmpi_begin_indep_data() : dispatcher->begin_indep_data()
 * ncmpi_end_indep_data()   : dispatcher->end_indep_data()
 * ncmpi_inq()              : dispatcher->inq()
 * ncmpi_inq_xxx()          : dispatcher->inq_misc()
 * ncmpi_sync_numrecs()     : dispatcher->sync_numrecs()
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>  /* strcpy() */
#include <assert.h>
#include <errno.h>

#include <mpi.h>

#include <pnc_debug.h>
#include <common.h>
#include "nc.h"
#include "ncio.h"
#include "ncx.h"
#include "fbits.h"
#ifdef ENABLE_SUBFILING
#include "subfile.h"
#endif

/*----< ncmpii_redef() >-----------------------------------------------------*/
/* This is a collective subroutine. */
int
ncmpii_redef(void *ncdp)
{
    int err=NC_NOERR;
    NC *ncp = (NC*)ncdp;

    if (NC_readonly(ncp)) DEBUG_RETURN_ERROR(NC_EPERM) /* read-only */
    /* if open mode is inconsistent, then this return might cause parallel
     * program to hang */

    /* cannot be in define mode, must enter from data mode */
    if (NC_indef(ncp)) DEBUG_RETURN_ERROR(NC_EINDEFINE)

    /* sync all metadata, including numrecs, if changed in independent mode.
     * also ensure exiting define mode always entering collective data mode
     */
    if (NC_indep(ncp)) /* exit independent mode, if in independent mode */
        ncmpii_end_indep_data(ncp);

    if (NC_doFsync(ncp)) { /* re-read the header from file */
        err = ncmpii_read_NC(ncp);
        if (err != NC_NOERR) return err;
    }

    /* duplicate a header to be used in enddef() for checking if header grows */
    ncp->old = ncmpii_dup_NC(ncp);
    if (ncp->old == NULL) DEBUG_RETURN_ERROR(NC_ENOMEM)

    /* we are now entering define mode */
    fSet(ncp->flags, NC_INDEF);

    return NC_NOERR;
}

/*----< ncmpiio_sync() >-----------------------------------------------------*/
/* This function must be called collectively, no matter if it is in collective
 * or independent data mode.
 */
int
ncmpiio_sync(ncio *nciop) {
#ifndef DISABLE_FILE_SYNC
    int mpireturn;

    if (nciop->independent_fh != MPI_FILE_NULL) {
        TRACE_IO(MPI_File_sync)(nciop->independent_fh);
        if (mpireturn != MPI_SUCCESS)
            return ncmpii_handle_error(mpireturn, "MPI_File_sync");
    }

    /* nciop->collective_fh is never MPI_FILE_NULL as collective mode is
     * default in PnetCDF */
    TRACE_IO(MPI_File_sync)(nciop->collective_fh);
    if (mpireturn != MPI_SUCCESS)
        return ncmpii_handle_error(mpireturn, "MPI_File_sync");

    TRACE_COMM(MPI_Barrier)(nciop->comm);
#endif
    return NC_NOERR;
}

/*----< ncmpii_begin_indep_data() >------------------------------------------*/
/* This is a collective subroutine. */
int
ncmpii_begin_indep_data(void *ncdp)
{
    int err=NC_NOERR;
    NC *ncp = (NC*)ncdp;

    if (NC_indef(ncp))  /* must not be in define mode */
        DEBUG_RETURN_ERROR(NC_EINDEFINE)

    if (NC_indep(ncp))  /* already in indep data mode */
        return NC_NOERR;

    /* we need no MPI_File_sync() here. If users want a stronger data
     * consistency, they can call ncmpi_sync()
     */
#if 0 && !defined(DISABLE_FILE_SYNC)
    if (!NC_readonly(ncp) && NC_collectiveFhOpened(ncp->nciop)) {
        /* calling file sync for those already open the file */
        int err, mpireturn;
        /* MPI_File_sync() is collective */
        TRACE_IO(MPI_File_sync)(ncp->nciop->collective_fh);
        if (mpireturn != MPI_SUCCESS) {
            err = ncmpii_handle_error(mpireturn, "MPI_File_sync");
            if (err == NC_NOERR) return err;
        }
        TRACE_COMM(MPI_Barrier)(ncp->nciop->comm);
    }
#endif

    /* raise independent flag */
    fSet(ncp->flags, NC_INDEP);

    /* open MPI_COMM_SELF file handler, if not yet opened */
    err = ncmpii_check_mpifh(ncp, 0);

    return err;
}

/*----< ncmpii_end_indep_data() >--------------------------------------------*/
/* This is a collective subroutine.
 * It can be called from:
 *    1. ncmpi_end_indep_data()
 *    2. ncmpi_redef() from independent data mode entering to define more
 *    3. ncmpii_close() when closing the file
 */
int
ncmpii_end_indep_data(void *ncdp)
{
    int status=NC_NOERR;
    NC *ncp = (NC*)ncdp;

    if (!NC_indep(ncp)) /* must be in independent data mode */
        DEBUG_RETURN_ERROR(NC_ENOTINDEP)

    if (!NC_readonly(ncp)) {
        if (ncp->vars.num_rec_vars > 0) {
            /* numrecs dirty bit may not be the same across all processes.
             * force sync in memory no matter if dirty or not.
             */
            set_NC_ndirty(ncp);
            status = ncmpiio_sync_numrecs(ncp, ncp->numrecs);
            /* the only possible dirty part of the header is numrecs */
        }

#ifndef DISABLE_FILE_SYNC
        /* calling file sync for those already open the file */
        if (NC_doFsync(ncp) && ncp->nciop->independent_fh != MPI_FILE_NULL) {
            int mpireturn;
            /* MPI_File_sync() is collective */
            TRACE_IO(MPI_File_sync)(ncp->nciop->independent_fh);
            if (mpireturn != MPI_SUCCESS) {
                int err = ncmpii_handle_error(mpireturn, "MPI_File_sync");
                if (status == NC_NOERR) status = err;
            }
            TRACE_COMM(MPI_Barrier)(ncp->nciop->comm);
            if (mpireturn != MPI_SUCCESS)
                return ncmpii_handle_error(mpireturn, "MPI_Barrier");
        }
#endif
    }

    fClr(ncp->flags, NC_INDEP);

    return status;
}

#define NC_NUMRECS_OFFSET 4

/*----< ncmpii_write_numrecs() >---------------------------------------------*/
/* root process writes the new record number into file.
 * This function is called by:
 * 1. ncmpiio_sync_numrecs
 * 2. collective nonblocking wait API, if the new number of records is bigger
 */
int
ncmpii_write_numrecs(NC         *ncp,
                     MPI_Offset  new_numrecs)
{
    int rank, mpireturn, err;
    MPI_File fh;

    /* root process writes numrecs in file */
    MPI_Comm_rank(ncp->nciop->comm, &rank);
    if (rank > 0) return NC_NOERR;

    /* return now if there is no record variabled defined */
    if (ncp->vars.num_rec_vars == 0) return NC_NOERR;

    fh = ncp->nciop->collective_fh;
    if (NC_indep(ncp))
        fh = ncp->nciop->independent_fh;

    if (new_numrecs > ncp->numrecs || NC_ndirty(ncp)) {
        int len;
        char pos[8], *buf=pos;
        MPI_Offset max_numrecs;
        MPI_Status mpistatus;

        max_numrecs = MAX(new_numrecs, ncp->numrecs);

        if (ncp->format < 5) {
            if (max_numrecs != (int)max_numrecs)
                DEBUG_RETURN_ERROR(NC_EINTOVERFLOW)
            len = X_SIZEOF_SIZE_T;
            err = ncmpix_put_uint32((void**)&buf, (uint)max_numrecs);
            if (err != NC_NOERR) DEBUG_RETURN_ERROR(err)
        }
        else {
            len = X_SIZEOF_INT64;
            err = ncmpix_put_uint64((void**)&buf, (uint64)max_numrecs);
            if (err != NC_NOERR) DEBUG_RETURN_ERROR(err)
        }
        /* ncmpix_put_xxx advances the 1st argument with size len */

        /* root's file view always includes the entire file header */
        TRACE_IO(MPI_File_write_at)(fh, NC_NUMRECS_OFFSET, (void*)pos, len,
                                    MPI_BYTE, &mpistatus);
        if (mpireturn != MPI_SUCCESS) {
            err = ncmpii_handle_error(mpireturn, "MPI_File_write_at");
            if (err == NC_EFILE) DEBUG_RETURN_ERROR(NC_EWRITE)
        }
        else {
            ncp->nciop->put_size += len;
        }
    }
    return NC_NOERR;
}

/*----< ncmpiio_sync_numrecs() >----------------------------------------------*/
/* Synchronize the number of records in memory and write numrecs to file.
 * This function is called by:
 * 1. ncmpi_sync_numrecs(): by the user
 * 2. ncmpi_sync(): by the user
 * 3. ncmpii_end_indep_data(): exit from independent data mode
 * 4. all blocking collective put APIs when writing record variables
 * 5. ncmpii_close(): file close and currently in independent data mode
 *
 * This function is collective.
 */
int
ncmpiio_sync_numrecs(NC         *ncp,
                     MPI_Offset  new_numrecs)
{
    int status=NC_NOERR, mpireturn;
    MPI_Offset max_numrecs;

    assert(!NC_readonly(ncp));
    assert(!NC_indef(ncp)); /* can only be called by APIs in data mode */

    /* return now if there is no record variabled defined */
    if (ncp->vars.num_rec_vars == 0) return NC_NOERR;

    /* find the max new_numrecs among all processes
     * Note new_numrecs may be smaller than ncp->numrecs
     */
    TRACE_COMM(MPI_Allreduce)(&new_numrecs, &max_numrecs, 1, MPI_OFFSET,
                              MPI_MAX, ncp->nciop->comm);
    if (mpireturn != MPI_SUCCESS)
        return ncmpii_handle_error(mpireturn, "MPI_Allreduce");

    /* root process writes max_numrecs to file */
    status = ncmpii_write_numrecs(ncp, max_numrecs);

    if (ncp->safe_mode == 1) {
        /* broadcast root's status, because only root writes to the file */
        int root_status = status;
        TRACE_COMM(MPI_Bcast)(&root_status, 1, MPI_INT, 0, ncp->nciop->comm);
        if (mpireturn != MPI_SUCCESS)
            return ncmpii_handle_error(mpireturn, "MPI_Bcast");
        /* root's write has failed, which is serious */
        if (root_status == NC_EWRITE) DEBUG_ASSIGN_ERROR(status, NC_EWRITE)
    }

    /* update numrecs in all processes's memory only if the new one is larger.
     * Note new_numrecs may be smaller than ncp->numrecs
     */
    if (max_numrecs > ncp->numrecs) ncp->numrecs = max_numrecs;

    /* clear numrecs dirty bit */
    fClr(ncp->flags, NC_NDIRTY);

    return status;
}

/*----< ncmpii_sync_numrecs() >-----------------------------------------------*/
/* this API is collective, but can be called in independent data mode.
 * Note numrecs is always sync-ed in memory and update in file in collective
 * data mode.
 */
int
ncmpii_sync_numrecs(void *ncdp)
{
    int err=NC_NOERR;
    NC *ncp=(NC*)ncdp;

    /* cannot be in define mode */
    if (NC_indef(ncp)) DEBUG_RETURN_ERROR(NC_EINDEFINE)

    /* check if we have defined record variables */
    if (ncp->vars.num_rec_vars == 0) return NC_NOERR;

    if (!NC_indep(ncp)) /* in collective data mode, numrecs is always sync-ed */
        return NC_NOERR;
    else /* if called in independent mode, we force sync in memory */
        set_NC_ndirty(ncp);

    /* sync numrecs in memory and file */
    err = ncmpiio_sync_numrecs(ncp, ncp->numrecs);

#ifndef DISABLE_FILE_SYNC
    if (NC_doFsync(ncp)) { /* NC_SHARE is set */
        int mpierr, mpireturn;
        if (NC_indep(ncp)) {
            TRACE_IO(MPI_File_sync)(ncp->nciop->independent_fh);
        }
        else {
            TRACE_IO(MPI_File_sync)(ncp->nciop->collective_fh);
        }
        if (mpireturn != MPI_SUCCESS) {
            mpierr = ncmpii_handle_error(mpireturn, "MPI_File_sync");
            if (err == NC_NOERR) err = mpierr;
        }
        TRACE_COMM(MPI_Barrier)(ncp->nciop->comm);
        if (mpireturn != MPI_SUCCESS)
            return ncmpii_handle_error(mpireturn, "MPI_Barrier");
    }
#endif
    return err;
}

/*----< ncmpii_sync() >------------------------------------------------------*/
/* This API is a collective subroutine, and must be called in data mode, no
 * matter if it is in collective or independent data mode.
 */
int
ncmpii_sync(void *ncdp)
{
    int err;
    NC *ncp = (NC*)ncdp;

    /* cannot be in define mode */
    if (NC_indef(ncp)) DEBUG_RETURN_ERROR(NC_EINDEFINE)

    if (NC_readonly(ncp))
        /* calling sync for file opened for read only means re-read header */
        return ncmpii_read_NC(ncp);

    /* the only part of header that can be dirty is numrecs (caused only by
     * independent APIs) */
    if (ncp->vars.num_rec_vars > 0 && NC_indep(ncp)) {
        /* sync numrecs in memory among processes and in file */
        set_NC_ndirty(ncp);
        err = ncmpiio_sync_numrecs(ncp, ncp->numrecs);
        if (err != NC_NOERR) return err;
    }

    /* calling MPI_File_sync() on both collective and independent handlers */
    return ncmpiio_sync(ncp->nciop);
}

/*----< ncmpii_abort() >-----------------------------------------------------*/
/* This API is a collective subroutine */
int
ncmpii_abort(void *ncdp)
{
   /*
    * In data mode, same as ncmpiio_close.
    * In define mode, descard new definition.
    * If file is just created, remove the file.
    */
    int status=NC_NOERR, err, doUnlink = 0;
    NC *ncp = (NC*)ncdp;

    /* delete the file if it is newly created by ncmpi_create() */
    doUnlink = NC_IsNew(ncp);

    if (ncp->old != NULL) {
        /* a plain redef, not a create */
        assert(!NC_IsNew(ncp));
        assert(fIsSet(ncp->flags, NC_INDEF));
        ncmpii_free_NC(ncp->old);
        ncp->old = NULL;
        fClr(ncp->flags, NC_INDEF);
    }

    if (!doUnlink) {
        if (!NC_readonly(ncp) &&  /* file is open for write */
             NC_indep(ncp)) {     /* in independent data mode */
            /* exit independent mode, if in independent mode */
            status = ncmpii_end_indep_data(ncp); /* will sync header */
        }

        if (NC_doFsync(ncp)) {
            err = ncmpiio_sync(ncp->nciop); /* calling MPI_File_sync() */
            if (status == NC_NOERR ) status = err;
        }
    }

    /* close the file */
    err = ncmpiio_close(ncp->nciop, doUnlink);
    if (status == NC_NOERR ) status = err;

    ncp->nciop = NULL;

    /* free up space occupied by the header metadata */
    ncmpii_free_NC(ncp);

    return status;
}

/*----< ncmpii_inq() >-------------------------------------------------------*/
int
ncmpii_inq(void *ncdp,
           int  *ndimsp,
           int  *nvarsp,
           int  *nattsp,
           int  *xtendimp)
{
    NC *ncp = (NC*)ncdp;

    if (ndimsp   != NULL) *ndimsp   = ncp->dims.ndefined;
    if (nvarsp   != NULL) *nvarsp   = ncp->vars.ndefined;
    if (nattsp   != NULL) *nattsp   = ncp->attrs.ndefined;
    if (xtendimp != NULL) *xtendimp = ncp->dims.unlimited_id;

    return NC_NOERR;
}

/*----< ncmpii_inq_misc() >--------------------------------------------------*/
/* This is an independent subroutine. */
int
ncmpii_inq_misc(void       *ncdp,
                int        *pathlen,
                char       *path,
                int        *num_fix_varsp,
                int        *num_rec_varsp,
                int        *striping_size,
                int        *striping_count,
                MPI_Offset *header_size,
                MPI_Offset *header_extent,
                MPI_Offset *recsize,
                MPI_Offset *put_size,
                MPI_Offset *get_size,
                MPI_Info   *info_used,
                int        *nreqs,
                MPI_Offset *usage,
                MPI_Offset *buf_size)
{
    int i, flag, mpireturn;
    char value[MPI_MAX_INFO_VAL];
    NC *ncp=(NC*)ncdp;

    /* Get the file pathname which was used to open/create the ncid's file.
     * path must already be allocated. Ignored if NULL */
    if (ncp->nciop->path == NULL) {
        if (pathlen != NULL) *pathlen = 0;
        if (path    != NULL) *path = '\0';
    } else {
        if (pathlen != NULL) *pathlen = (int)strlen(ncp->nciop->path);
        if (path    != NULL) strcpy(path, ncp->nciop->path);
    }

    /* obtain the number of fixed-size variables */
    if (num_fix_varsp != NULL) {
        if (NC_indef(ncp)) {
            /* if in define mode, recalculate the number of record variables */
            *num_fix_varsp = 0;
            for (i=0; i<ncp->vars.ndefined; i++)
                *num_fix_varsp += IS_RECVAR(ncp->vars.value[i]);
        }
        else
            *num_fix_varsp = ncp->vars.num_rec_vars;

        /* no. fixed-size == ndefined - no. record variables */
        *num_fix_varsp = ncp->vars.ndefined - *num_fix_varsp;
    }

    /* obtain the number of record variables */
    if (num_rec_varsp != NULL) {
        if (NC_indef(ncp)) {
            /* if in define mode, recalculate the number of record variables */
            *num_rec_varsp = 0;
            for (i=0; i<ncp->vars.ndefined; i++)
                *num_rec_varsp += IS_RECVAR(ncp->vars.value[i]);
        }
        else
            *num_rec_varsp = ncp->vars.num_rec_vars;
    }

    /* obtain file (system) striping settings, striping size and count, if they
     * are available from MPI-IO hint. Otherwise, 0s are returned.
     */
    if (striping_size != NULL) {
        MPI_Info_get(ncp->nciop->mpiinfo, "striping_unit", MPI_MAX_INFO_VAL-1,
                     value, &flag);
        *striping_size = 0;
        if (flag) {
            errno = 0;  /* errno must set to zero before calling strtoll */
            *striping_size = (int)strtol(value,NULL,10);
            if (errno != 0) *striping_size = 0;
        }
    }

    if (striping_count != NULL) {
        MPI_Info_get(ncp->nciop->mpiinfo, "striping_factor", MPI_MAX_INFO_VAL-1,
                     value, &flag);
        *striping_count = 0;
        if (flag) {
            errno = 0;  /* errno must set to zero before calling strtoll */
            *striping_count = (int)strtol(value,NULL,10);
            if (errno != 0) *striping_count = 0;
        }
    }

    /* the amount of writes, in bytes, committed to file system so far */
    if (put_size != NULL) *put_size = ncp->nciop->put_size;

    /* the amount of reads, in bytes, obtained from file system so far */
    if (get_size != NULL) *get_size = ncp->nciop->get_size;

    if (recsize != NULL) *recsize = ncp->recsize;

    if (header_size != NULL) *header_size = ncp->xsz;

    if (header_extent != NULL) *header_extent = ncp->begin_var;

    if (info_used != NULL) {
        mpireturn = MPI_Info_dup(ncp->nciop->mpiinfo, info_used);
        if (mpireturn != MPI_SUCCESS)
            return ncmpii_handle_error(mpireturn, "MPI_Info_dup");

        sprintf(value, "%lld", ncp->h_align);
        MPI_Info_set(*info_used, "nc_header_align_size", value);

        sprintf(value, "%lld", ncp->v_align);
        MPI_Info_set(*info_used, "nc_var_align_size", value);

        sprintf(value, "%lld", ncp->r_align);
        MPI_Info_set(*info_used, "nc_record_align_size", value);

        sprintf(value, "%lld", ncp->chunk);
        MPI_Info_set(*info_used, "nc_header_read_chunk_size", value);

#ifdef ENABLE_SUBFILING
        sprintf(value, "%d", ncp->subfile_mode);
        MPI_Info_set(*info_used, "pnetcdf_subfiling", value);
        sprintf(value, "%d", ncp->num_subfiles);
        MPI_Info_set(*info_used, "nc_num_subfiles", value);
#endif
    }

    if (nreqs != NULL) {
        /* cannot just use *nreqs = ncp->numGetReqs + ncp->numPutReqs;
         * because some request IDs are repeated, such as record variables and
         * varn requests
         */
        *nreqs = 0;
        for (i=0; i<ncp->numGetReqs; i++) {
            if (i > 0 && ncp->get_list[i].id == ncp->get_list[i-1].id)
                continue;
            (*nreqs)++;
        }
        for (i=0; i<ncp->numPutReqs; i++) {
            if (i > 0 && ncp->put_list[i].id == ncp->put_list[i-1].id)
                continue;
            (*nreqs)++;
        }
    }

    if (usage != NULL) {
        /* check if the buffer has been previously attached */
        if (ncp->abuf == NULL) DEBUG_RETURN_ERROR(NC_ENULLABUF)
        /* return the current usage in bytes */
        *usage = ncp->abuf->size_used;
    }

    if (buf_size != NULL) {
        /* check if the buffer has been previously attached */
        if (ncp->abuf == NULL) DEBUG_RETURN_ERROR(NC_ENULLABUF)
        /* return the current usage in bytes */
        *buf_size = ncp->abuf->size_allocated;
    }

    return NC_NOERR;
}

/*----< ncmpi_delete() >-----------------------------------------------------*/
/* doesn't do anything to release resources. Users are advised to call
 * ncmpi_close() before calling this function.
 *
 * filename: the name of the file we will remove.
 * info: MPI info object, in case underlying file system needs hints.
 */
int
ncmpi_delete(const char *filename,
             MPI_Info    info)
{
    int err=NC_NOERR, mpireturn;

    TRACE_IO(MPI_File_delete)((char*)filename, info);
    if (mpireturn != MPI_SUCCESS)
        err = ncmpii_handle_error(mpireturn, "MPI_File_delete");
    return err;
}

