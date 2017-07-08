/*
 *  Copyright (C) 2003, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* $Id$ */

/*
 * This file implements the corresponding APIs defined in src/dispatchers/file.c
 *
 * ncmpi_close()            : dispatcher->close()
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <errno.h>

#include <mpi.h>

#include <pnc_debug.h>
#include <common.h>
#include "nc.h"
#include "fbits.h"
#ifdef ENABLE_SUBFILING
#include "subfile.h"
#endif

/*----< ncmpiio_close() >----------------------------------------------------*/
int
ncmpiio_close(NC *ncp, int doUnlink) {
    int mpireturn;

    assert(ncp != NULL); /* this should never occur */

    if (ncp->independent_fh != MPI_FILE_NULL) {
        TRACE_IO(MPI_File_close)(&ncp->independent_fh);
        if (mpireturn != MPI_SUCCESS)
            return ncmpii_handle_error(mpireturn, "MPI_File_close");
    }

    if (ncp->collective_fh != MPI_FILE_NULL) {
        TRACE_IO(MPI_File_close)(&ncp->collective_fh);
        if (mpireturn != MPI_SUCCESS)
            return ncmpii_handle_error(mpireturn, "MPI_File_close");
    }

    if (doUnlink) {
        /* called from ncmpi_abort, if the file is being created and is still
         * in define mode, the file is deleted */
        TRACE_IO(MPI_File_delete)((char *)ncp->path, ncp->mpiinfo);
        if (mpireturn != MPI_SUCCESS)
            return ncmpii_handle_error(mpireturn, "MPI_File_delete");
    }

    /* free MPI objects */
    if (ncp->mpiinfo != MPI_INFO_NULL) MPI_Info_free(&(ncp->mpiinfo));
    if (ncp->comm != MPI_COMM_NULL)    MPI_Comm_free(&(ncp->comm));
    NCI_Free(ncp->path);

    return NC_NOERR;
}

/*----< ncmpii_close() >------------------------------------------------------*/
/* This function is collective */
int
ncmpii_close(void *ncdp)
{
    int err=NC_NOERR, status=NC_NOERR;
    NC *ncp = (NC*)ncdp;

    if (NC_indef(ncp)) { /* currently in define mode */
        status = ncmpii__enddef(ncp, 0, 0, 0, 0); /* TODO: defaults */

        if (status != NC_NOERR ) {
            /* To do: Abort new definition, if any */
            if (ncp->old != NULL) {
                ncmpii_free_NC(ncp->old);
                ncp->old = NULL;
                fClr(ncp->flags, NC_INDEF);
            }
        }
    }

    if (!NC_readonly(ncp) &&  /* file is open for write */
         NC_indep(ncp)) {     /* exit independent data mode will sync header */
        err = ncmpii_end_indep_data(ncp);
        if (status == NC_NOERR ) status = err;
    }

    /* if entering this function in  collective data mode, we do not have to
     * update header in file, as file header is always up-to-date */

#ifdef ENABLE_SUBFILING
    /* ncmpii__enddef() will update ncp->num_subfiles */
    /* TODO: should check ncid_sf? */
    /* if the file has subfiles, close them first */
    if (ncp->num_subfiles > 1)
        ncmpii_subfile_close(ncp);
#endif

    /* We can cancel or complete all outstanding nonblocking I/O.
     * For now, cancelling makes more sense. */
#ifdef COMPLETE_NONBLOCKING_IO
    if (ncp->numGetReqs > 0) {
        ncmpii_wait(ncp, NC_GET_REQ_ALL, NULL, NULL, INDEP_IO);
        if (status == NC_NOERR ) status = NC_EPENDING;
    }
    if (ncp->numPutReqs > 0) {
        ncmpii_wait(ncp, NC_PUT_REQ_ALL, NULL, NULL, INDEP_IO);
        if (status == NC_NOERR ) status = NC_EPENDING;
    }
#else
    if (ncp->numGetReqs > 0) {
        int rank;
        MPI_Comm_rank(ncp->comm, &rank);
        printf("PnetCDF warning: %d nonblocking get requests still pending on process %d. Cancelling ...\n",ncp->numGetReqs,rank);
        ncmpii_cancel(ncp, NC_GET_REQ_ALL, NULL, NULL);
        if (status == NC_NOERR ) status = NC_EPENDING;
    }
    if (ncp->numPutReqs > 0) {
        int rank;
        MPI_Comm_rank(ncp->comm, &rank);
        printf("PnetCDF warning: %d nonblocking put requests still pending on process %d. Cancelling ...\n",ncp->numPutReqs,rank);
        ncmpii_cancel(ncp, NC_PUT_REQ_ALL, NULL, NULL);
        if (status == NC_NOERR ) status = NC_EPENDING;
    }
#endif

    /* If the user wants a stronger data consistency by setting NC_SHARE */
    if (fIsSet(ncp->iomode, NC_SHARE))
        ncmpiio_sync(ncp); /* calling MPI_File_sync() */

    /* calling MPI_File_close() */
    ncmpiio_close(ncp, 0);

    /* free up space occupied by the header metadata */
    ncmpii_free_NC(ncp);

    return status;
}

