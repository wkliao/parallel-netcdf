/*
 *  Copyright (C) 2003, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* $Id$ */

/*
 * This file implements the corresponding APIs defined in src/dispatchers/file.c
 *
 * ncmpi_open() : dispatcher->open()
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>  /* strcpy() */
#ifdef HAVE_ACCESS
#include <unistd.h>  /* access() */
#endif

#include <mpi.h>

#include <pnc_debug.h>
#include <common.h>
#include "nc.h"
#include "fbits.h"
#ifdef ENABLE_SUBFILING
#include "subfile.h"
#endif

/*----< ncmpii_open() >------------------------------------------------------*/
int
ncmpii_open(MPI_Comm     comm,
            const char  *path,
            int          omode,
            int          ncid,
            MPI_Info     info, /* user's info and env info combined */
            void       **ncpp)
{
    char *env_str;
    int i, mpiomode, err, mpireturn;
    MPI_File fh;
    MPI_Comm dup_comm;
    MPI_Info info_used;
    NC *ncp=NULL;
#ifndef SEARCH_NAME_LINEARLY
    NC_nametable *nameT;
#endif

    *ncpp = NULL;

    /* Note path's validity and omode consistency have been checked in
     * ncmpi_open() in src/dispatchers/file.c and
     * path consistency will be done in MPI_File_open */

    /* First, check whether omode is valid or supported ---------------------*/
    /* NC_DISKLESS is not supported yet */
    if (omode & NC_DISKLESS) DEBUG_RETURN_ERROR(NC_EINVAL_OMODE)

    /* NC_MMAP is not supported yet */
    if (omode & NC_MMAP) DEBUG_RETURN_ERROR(NC_EINVAL_OMODE)

#if 0 && defined(HAVE_ACCESS)
    if (mpiomode == MPI_MODE_RDONLY) { /* file should already exit */
        int rank, file_exist;
        MPI_Comm_rank(comm, &rank);
        if (rank == 0) {
            if (access(path, F_OK) == 0) file_exist = 1;
            else                         file_exist = 0;
        }
        TRACE_COMM(MPI_Bcast)(&file_exist, 1, MPI_INT, 0, comm);
        if (!file_exist) DEBUG_RETURN_ERROR(NC_ENOENT)
    }
#endif

    /* open file collectively ---------------------------------------------- */
    mpiomode = fIsSet(omode, NC_WRITE) ? MPI_MODE_RDWR : MPI_MODE_RDONLY;

    TRACE_IO(MPI_File_open)(comm, (char *)path, mpiomode, info, &fh);
    if (mpireturn != MPI_SUCCESS)
        return ncmpii_handle_error(mpireturn, "MPI_File_open");

    /* duplicate MPI communicator as user may free it later */
    mpireturn = MPI_Comm_dup(comm, &dup_comm);
    if (mpireturn != MPI_SUCCESS)
        return ncmpii_handle_error(mpireturn, "MPI_Comm_dup");

    /* get the file info used by MPI-IO */
    mpireturn = MPI_File_get_info(fh, &info_used);
    if (mpireturn != MPI_SUCCESS) {
        MPI_Comm_free(&dup_comm);
        return ncmpii_handle_error(mpireturn, "MPI_File_get_info");
    }

    /* Now the file has been successfully opened, allocate/set NC object */

    /* path's validity and omode consistency have been checked in ncmpi_open()
     * in src/dispatchers/file.c */

    /* allocate buffer for header object NC */
    ncp = (NC*) NCI_Calloc(1, sizeof(NC));
    if (ncp == NULL) DEBUG_RETURN_ERROR(NC_ENOMEM)

    /* PnetCDF default fill mode is no fill */
    fSet(ncp->flags, NC_NOFILL);

    ncp->ncid         = ncid;
    ncp->safe_mode    = 0;
    ncp->numGetReqs   = 0;
    ncp->numPutReqs   = 0;
#ifdef ENABLE_SUBFILING
    ncp->subfile_mode = 0;
    ncp->num_subfiles = 0;
    ncp->ncid_sf      = -1; /* subfile ncid; init to -1 */
#endif

    ncp->chunk        = NC_DEFAULT_CHUNKSIZE;
    ncp->h_align      = 0; /* value 0 indicates the hint is not set */
    ncp->v_align      = 0;
    ncp->r_align      = 0;
    ncp->h_minfree    = 0;
    ncp->v_minfree    = 0;

    /* extract I/O hints from user info */
    ncmpio_set_pnetcdf_hints(ncp, info);

    ncp->get_list     = NULL;
    ncp->put_list     = NULL;
    ncp->abuf         = NULL;
    ncp->old          = NULL;

    ncp->iomode         = omode;
    ncp->comm           = dup_comm;
    ncp->mpiinfo        = info_used;
    ncp->mpiomode       = mpiomode;
    ncp->put_size       = 0;
    ncp->get_size       = 0;
    ncp->collective_fh  = fh;
    ncp->independent_fh = MPI_FILE_NULL;
    ncp->path = (char*) NCI_Malloc(strlen(path) + 1);
    strcpy(ncp->path, path);

#ifdef PNETCDF_DEBUG
    /* PNETCDF_DEBUG is set at configure time, which will be overwritten by
     * the run-time environment variable PNETCDF_SAFE_MODE */
    ncp->safe_mode = 1;
#endif
    /* If environment variable PNETCDF_SAFE_MODE is set to 1, then we perform
     * a strict consistent test, i.e. arguments used in def_dim/def_var APIs
     */
    if ((env_str = getenv("PNETCDF_SAFE_MODE")) != NULL) {
        if (*env_str == '0') ncp->safe_mode = 0;
        else                 ncp->safe_mode = 1;
        /* if PNETCDF_SAFE_MODE is set but without a value, *env_str can
         * be '\0' (null character). In this case, safe_mode is enabled */
    }

    /* read header from file into NC object pointed by ncp -------------------*/
    err = ncmpii_hdr_get_NC(ncp);
    if (err != NC_NOERR) { /* fatal error */
        ncmpio_close_files(ncp, 0);
        ncmpii_free_NC(ncp);
        return err;
    }

#ifdef ENABLE_SUBFILING
    if (ncp->subfile_mode) {
        /* check subfiling attribute */
        err = ncmpii_get_att(ncp, NC_GLOBAL, "num_subfiles", &ncp->num_subfiles,
                             NC_INT);
        if (err == NC_NOERR && ncp->num_subfiles > 1) {
            /* ignore error NC_ENOTATT if this attribute is not defined */
            for (i=0; i<ncp->vars.ndefined; i++) {
                err = ncmpii_get_att(ncp, i, "num_subfiles",
                                     &ncp->vars.value[i]->num_subfiles, NC_INT);
                if (err == NC_ENOTATT) continue;
                if (err != NC_NOERR) return err;

                if (ncp->vars.value[i]->num_subfiles > 1) {
                    err = ncmpii_get_att(ncp, i, "ndims_org",
                                         &ncp->vars.value[i]->ndims_org,NC_INT);
                    if (err != NC_NOERR) return err;
                }
            }
            if (ncp->num_subfiles > 1) {
                err = ncmpii_subfile_open(ncp, &ncp->ncid_sf);
                if (err != NC_NOERR) return err;
            }
        }
    }
    else
        ncp->num_subfiles = 0;
#endif

    /* update the total number of record variables --------------------------*/
    ncp->vars.num_rec_vars = 0;
    for (i=0; i<ncp->vars.ndefined; i++)
        ncp->vars.num_rec_vars += IS_RECVAR(ncp->vars.value[i]);

#ifndef SEARCH_NAME_LINEARLY
    /* initialize dim name lookup table -------------------------------------*/
    nameT = ncp->dims.nameT;
    memset(nameT, 0, sizeof(NC_nametable) * HASH_TABLE_SIZE);

    /* populate dim name lookup table */
    for (i=0; i<ncp->dims.ndefined; i++) {
        /* hash the dim name into a key for name lookup */
        int key = HASH_FUNC(ncp->dims.value[i]->name->cp);
        nameT = &ncp->dims.nameT[key];
        if (nameT->num % NC_NAME_TABLE_CHUNK == 0)
            nameT->list = (int*) NCI_Realloc(nameT->list,
                          (size_t)(nameT->num+NC_NAME_TABLE_CHUNK) *SIZEOF_INT);
        nameT->list[nameT->num] = i;
        nameT->num++;
    }

    /* initialize var name lookup table */
    nameT = ncp->vars.nameT;
    memset(nameT, 0, sizeof(NC_nametable) * HASH_TABLE_SIZE);

    /* populate var name lookup table */
    for (i=0; i<ncp->vars.ndefined; i++) {
        /* hash the var name into a key for name lookup */
        int key = HASH_FUNC(ncp->vars.value[i]->name->cp);
        nameT = &ncp->vars.nameT[key];
        if (nameT->num % NC_NAME_TABLE_CHUNK == 0)
            nameT->list = (int*) NCI_Realloc(nameT->list,
                          (size_t)(nameT->num+NC_NAME_TABLE_CHUNK) *SIZEOF_INT);
        nameT->list[nameT->num] = i;
        nameT->num++;
    }
#endif

    *ncpp = (void*)ncp;

    return NC_NOERR;
}

