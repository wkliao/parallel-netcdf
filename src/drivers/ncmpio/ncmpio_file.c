/*
 *  Copyright (C) 2003, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* $Id$ */

/*
 * This file implements the corresponding APIs defined in src/dispatchers/file.c
 *
 * ncmpi_create()           : dispatcher->create()
 * ncmpi_open()             : dispatcher->open()
 * ncmpi_close()            : dispatcher->close()
 * ncmpi_redef()            : dispatcher->redef()
 * ncmpi_sync()             : dispatcher->sync()
 * ncmpi_abort()            : dispatcher->abort()
 * ncmpi_inq()              : dispatcher->inq()
 * ncmpi_begin_indep_data() : dispatcher->begin_indep_data()
 * ncmpi_end_indep_data()   : dispatcher->end_indep_data()
 * ncmpi_inq_xxx()          : dispatcher->inq_misc()
 * ncmpi_sync_numrecs()     : dispatcher->sync_numrecs()
 * ncmpi_buffer_attach()    : dispatcher->buffer_attach()
 * ncmpi_buffer_detach()    : dispatcher->buffer_detach()
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>  /* strtok(), strcpy(), strchr() */
#include <strings.h> /* strcasecmp() */
#ifdef HAVE_ACCESS
#include <unistd.h>  /* access() */
#endif
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
#include "rnd.h"    /* M_RNDUP() */

/*----< ncio_free() >--------------------------------------------------------*/
static void
ncio_free(ncio *nciop) {
    if (nciop != NULL) {
        if (nciop->mpiinfo != MPI_INFO_NULL)
            MPI_Info_free(&(nciop->mpiinfo));
        if (nciop->comm != MPI_COMM_NULL)
            MPI_Comm_free(&(nciop->comm));

        NCI_Free(nciop);
    }
}

/*----< ncio_new() >---------------------------------------------------------*/
static ncio*
ncio_new(const char *path, int ioflags)
{
    size_t sz_ncio = M_RNDUP(sizeof(ncio));
    size_t sz_path = M_RNDUP(strlen(path) + 1);
    ncio *nciop;

    nciop = (ncio *) NCI_Malloc(sz_ncio + sz_path);
    if (nciop == NULL) return NULL;

    nciop->ioflags  = ioflags;
    nciop->comm     = MPI_COMM_NULL;
    nciop->mpiinfo  = MPI_INFO_NULL;
    nciop->put_size = 0;
    nciop->get_size = 0;

    nciop->path = (char *) ((char *)nciop + sz_ncio);
    (void) strcpy((char *)nciop->path, path);

    return nciop;
}

/*----< set_pnetcdf_hints() >------------------------------------------------*/
/* this is where the I/O hints designated to pnetcdf are extracted */
static
void set_pnetcdf_hints(NC *ncp, MPI_Info  info)
{
    char value[MPI_MAX_INFO_VAL];
    int  flag;

    /* value 0 indicates the hint is not set */
    ncp->h_align      = 0;
    ncp->v_align      = 0;
    ncp->r_align      = 0;
    ncp->h_minfree    = 0;
    ncp->v_minfree    = 0;
    ncp->chunk        = 0;
#ifdef ENABLE_SUBFILING
    ncp->subfile_mode = 0;
    ncp->num_subfiles = 0;
#endif

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

/*----< ncmpiio_create() >---------------------------------------------------*/
static int
ncmpiio_create(MPI_Comm     comm,
               const char  *path,
               int          ioflags,
               MPI_Info     info,
               NC          *ncp)
{
    ncio *nciop;
    int rank, mpireturn, err;
    int mpiomode = MPI_MODE_RDWR | MPI_MODE_CREATE;

    /* checking path consistency is expected to be done in MPI-IO */

    MPI_Comm_rank(comm, &rank);

    /* NC_CLOBBER is the default mode, even if it is not used in cmode.
     * Note ioflags has been checked for consistency before entering this API.
     */
    if (fIsSet(ioflags, NC_NOCLOBBER)) {
        /* check if file exists: NetCDF requires NC_EEXIST returned if the file
         * already exists and NC_NOCLOBBER mode is used in create
         */
#ifdef HAVE_ACCESS
        int file_exist;
        /* if access() is available, use it to check whether file already exists
         * rank 0 calls access() and broadcasts file_exist */
        if (rank == 0) {
            /* remove the file system type prefix name if there is any.
             * For example, path=="lustre:/home/foo/testfile.nc",
             * use "/home/foo/testfile.nc" when calling access()
             */
            char *filename = strchr(path, ':');
            if (filename == NULL) /* no prefix */
                filename = (char*)path;
            else
                filename++;

            if (access(filename, F_OK) == 0) file_exist = 1;
            else                             file_exist = 0;
        }
        TRACE_COMM(MPI_Bcast)(&file_exist, 1, MPI_INT, 0, comm);
        if (file_exist) DEBUG_RETURN_ERROR(NC_EEXIST)
#else
        /* use MPI_MODE_EXCL mode in MPI_File_open and check returned error */
        fSet(mpiomode, MPI_MODE_EXCL);
#endif
    }
    else { /* NC_CLOBBER is the default mode in create */
        /* rank 0 deletes the file and ignores error code for file not exist
         * Note calling MPI_File_set_size is expensive as it calls truncate()
         */
        if (rank == 0) {
#ifdef HAVE_UNLINK
            err = unlink(path);
            if (err < 0 && errno != ENOENT) /* ignore ENOENT: file not exist */
                DEBUG_ASSIGN_ERROR(err, NC_EFILE) /* other error */
            else
                err = NC_NOERR;
#else
            err = NC_NOERR;
            TRACE_IO(MPI_File_delete)((char*)path, MPI_INFO_NULL);
            if (mpireturn != MPI_SUCCESS) {
                int errorclass;
                MPI_Error_class(mpireturn, &errorclass);
                if (errorclass != MPI_ERR_NO_SUCH_FILE) /* ignore this error */
                    err = ncmpii_handle_error(mpireturn, "MPI_File_delete");
            }
#endif
        }
        /* all processes must wait here until file deletion is completed */
        TRACE_COMM(MPI_Bcast)(&err, 1, MPI_INT, 0, comm);
        if (err != NC_NOERR) return err;
    }

    /* ignore if NC_NOWRITE set by user */
    fSet(ioflags, NC_WRITE);

    /* allocate ncio object */
    nciop = ncio_new(path, ioflags);
    if (nciop == NULL) DEBUG_RETURN_ERROR(NC_ENOMEM)

    nciop->mpiomode  = MPI_MODE_RDWR;
    nciop->mpioflags = 0;

    /* open file in parallel */
    TRACE_IO(MPI_File_open)(comm, (char *)path, mpiomode, info,
                            &nciop->collective_fh);
    if (mpireturn != MPI_SUCCESS) {
        ncio_free(nciop);
#ifndef HAVE_ACCESS
        if (fIsSet(ioflags, NC_NOCLOBBER)) {
            /* This is the case when NC_NOCLOBBER is used in file creation and
             * function access() is not available. MPI_MODE_EXCL is set in open
             * mode. When MPI_MODE_EXCL is used and the file already exists,
             * MPI-IO should return error class MPI_ERR_FILE_EXISTS. But, some
             * MPI-IO implementations (older ROMIO) do not correctly return
             * this error class. In this case, we can do the followings: check
             * errno to see if it set to EEXIST. Note usually rank 0 makes the
             * file open call and can be the only one having errno set.
             */
            TRACE_COMM(MPI_Bcast)(&errno, 1, MPI_INT, 0, comm);
            if (errno == EEXIST) DEBUG_RETURN_ERROR(NC_EEXIST)
        }
#endif
        return ncmpii_handle_error(mpireturn, "MPI_File_open");
        /* for NC_NOCLOBBER, MPI_MODE_EXCL was added to mpiomode. If the file
         * already exists, MPI-IO should return error class MPI_ERR_FILE_EXISTS
         * which PnetCDF will return error code NC_EEXIST. This is checked
         * inside of ncmpii_handle_error()
         */
    }

    /* collective I/O mode is the default mode */
    set_NC_collectiveFh(nciop);

    /* duplicate communicator as user may free it later */
    mpireturn = MPI_Comm_dup(comm, &(nciop->comm));
    if (mpireturn != MPI_SUCCESS)
        return ncmpii_handle_error(mpireturn, "MPI_Comm_dup");

    /* get the file info actually used by MPI-IO (maybe alter user's info) */
    mpireturn = MPI_File_get_info(nciop->collective_fh, &nciop->mpiinfo);
    if (mpireturn != MPI_SUCCESS)
        return ncmpii_handle_error(mpireturn, "MPI_File_get_info");

    ncp->nciop = nciop;
    return NC_NOERR;
}

/*----< ncmpiio_open() >-----------------------------------------------------*/
static int
ncmpiio_open(MPI_Comm     comm,
             const char  *path,
             int          ioflags,
             MPI_Info     info,
             NC          *ncp)
{
    ncio *nciop;
    int mpireturn;
    int mpiomode = fIsSet(ioflags, NC_WRITE) ? MPI_MODE_RDWR : MPI_MODE_RDONLY;

    /* Note ioflags has been checked for consistency before entering this API.
     */

    assert(ncp != NULL);

    /* checking path consistency is expected done in MPI-IO */

    /* When open an non-existing file for read, we can either call access() to
     * check and return error code NC_ENOENT, or call MPI_File_open and expect
     * error class MPI_ERR_NO_SUCH_FILE. For now, we let MPI-IO to check.
     */
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

    /* allocate ncio object */
    nciop = ncio_new(path, ioflags);
    if (nciop == NULL) DEBUG_RETURN_ERROR(NC_ENOMEM)

    nciop->mpiomode  = mpiomode;
    nciop->mpioflags = 0;

    /* open file in parallel */
    TRACE_IO(MPI_File_open)(comm, (char *)path, mpiomode, info,
                            &nciop->collective_fh);
    if (mpireturn != MPI_SUCCESS) {
        ncio_free(nciop);
        return ncmpii_handle_error(mpireturn, "MPI_File_open");
    }

    /* default mode is collective */
    set_NC_collectiveFh(nciop);

    /* duplicate MPI communicator as user may free it later */
    mpireturn = MPI_Comm_dup(comm, &(nciop->comm));
    if (mpireturn != MPI_SUCCESS)
        return ncmpii_handle_error(mpireturn, "MPI_Comm_dup");

    /* get the file info used by MPI-IO */
    mpireturn = MPI_File_get_info(nciop->collective_fh, &nciop->mpiinfo);
    if (mpireturn != MPI_SUCCESS)
        return ncmpii_handle_error(mpireturn, "MPI_File_get_info");

    ncp->nciop = nciop;
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

    if (NC_independentFhOpened(nciop)) {
        TRACE_IO(MPI_File_sync)(nciop->independent_fh);
        if (mpireturn != MPI_SUCCESS)
            return ncmpii_handle_error(mpireturn, "MPI_File_sync");
    }
    if (NC_collectiveFhOpened(nciop)) {
        TRACE_IO(MPI_File_sync)(nciop->collective_fh);
        if (mpireturn != MPI_SUCCESS)
            return ncmpii_handle_error(mpireturn, "MPI_File_sync");
    }
    TRACE_COMM(MPI_Barrier)(nciop->comm);
#endif
    return NC_NOERR;
}

/*----< ncmpiio_close() >----------------------------------------------------*/
int
ncmpiio_close(ncio *nciop, int doUnlink) {
    int mpireturn;

    if (nciop == NULL) /* this should never occur */
        DEBUG_RETURN_ERROR(NC_EINVAL)

    if (NC_independentFhOpened(nciop)) {
        TRACE_IO(MPI_File_close)(&(nciop->independent_fh));
        if (mpireturn != MPI_SUCCESS)
            return ncmpii_handle_error(mpireturn, "MPI_File_close");
    }

    if (NC_collectiveFhOpened(nciop)) {
        TRACE_IO(MPI_File_close)(&(nciop->collective_fh));
        if (mpireturn != MPI_SUCCESS)
            return ncmpii_handle_error(mpireturn, "MPI_File_close");
    }

    if (doUnlink) {
        TRACE_IO(MPI_File_delete)((char *)nciop->path, nciop->mpiinfo);
        if (mpireturn != MPI_SUCCESS)
            return ncmpii_handle_error(mpireturn, "MPI_File_delete");
    }
    ncio_free(nciop);

    return NC_NOERR;
}

/*----< ncmpii_create() >----------------------------------------------------*/
int
ncmpii_create(MPI_Comm     comm,
              const char  *path,
              int          cmode,
              int          ncid,
              MPI_Info     info, /* user's info and env info combined */
              void       **ncpp)
{
    char *env_str, value[MPI_MAX_INFO_VAL];
    int i, flag, err, safe_mode=0, mpireturn, default_format;
    MPI_Offset chunksize=NC_DEFAULT_CHUNKSIZE;
    NC *ncp=NULL;

#ifdef PNETCDF_DEBUG
    /* PNETCDF_DEBUG is set at configure time, which will be overwritten by
     * the run-time environment variable PNETCDF_SAFE_MODE */
    safe_mode = 1;
#endif
    /* If environment variable PNETCDF_SAFE_MODE is set to 1, then we perform
     * a strict consistent test, i.e. arguments used in def_dim/def_var APIs
     */
    if ((env_str = getenv("PNETCDF_SAFE_MODE")) != NULL) {
        if (*env_str == '0') safe_mode = 0;
        else                 safe_mode = 1;
        /* if PNETCDF_SAFE_MODE is set but without a value, *env_str can
         * be '\0' (null character). In this case, safe_mode is enabled */
    }

    /* path's validity and cmode consistency have been checked in
     * ncmpi_create() in src/dispatchers/file.c */

    /* use default format, if cmode does not include either NC_64BIT_OFFSET
     * or NC_64BIT_DATA */
    ncmpi_inq_default_format(&default_format);

#if SIZEOF_MPI_OFFSET <  8
    /* check cmode */
    if (fIsSet(cmode, NC_64BIT_DATA)     ||
        fIsSet(cmode, NC_64BIT_OFFSET)   ||
        default_format == NC_FORMAT_CDF5 || 
        default_format == NC_FORMAT_CDF2) {
        /* unlike serial netcdf, we will not bother to support
         * NC_64BIT_OFFSET on systems with off_t smaller than 8 bytes.
         * serial netcdf has proven it's possible if datasets are small, but
         * that's a hassle we don't want to worry about */
        DEBUG_RETURN_ERROR(NC_ESMALL)
    }
#endif

    /* NC_DISKLESS is not supported yet */
    if (cmode & NC_DISKLESS) DEBUG_RETURN_ERROR(NC_EINVAL_CMODE)

    /* NC_MMAP is not supported yet */
    if (cmode & NC_MMAP) DEBUG_RETURN_ERROR(NC_EINVAL_CMODE)

    /* It is illegal to have both NC_64BIT_OFFSET & NC_64BIT_DATA */
    if ((cmode & (NC_64BIT_OFFSET|NC_64BIT_DATA)) ==
                 (NC_64BIT_OFFSET|NC_64BIT_DATA)) {
        DEBUG_RETURN_ERROR(NC_EINVAL_CMODE)
    }

    /* allocate buffer for header object NC */
    ncp = (NC*) NCI_Calloc(1, sizeof(NC));
    if (ncp == NULL) DEBUG_RETURN_ERROR(NC_ENOMEM)

    ncp->safe_mode         = safe_mode;
    ncp->ncid              = ncid;
    ncp->abuf              = NULL;
    ncp->old               = NULL;
    ncp->chunk             = NC_DEFAULT_CHUNKSIZE;
    /* initialize arrays storing pending non-blocking requests */
    ncp->numGetReqs        = 0;
    ncp->numPutReqs        = 0;
    ncp->get_list          = NULL;
    ncp->put_list          = NULL;
    /* initialize unlimited_id as no unlimited dimension yet defined */
    ncp->dims.unlimited_id = -1;
    /* find the true header size (not-yet aligned) */
    ncp->xsz               = ncmpii_hdr_len_NC(ncp);
#ifdef ENABLE_SUBFILING
    ncp->ncid_sf           = -1; /* subfile ncid; init to -1 */
#endif

    /* PnetCDF default fill mode is no fill */
    fSet(ncp->flags, NC_NOFILL);
    fSet(ncp->flags, NC_CREAT);

#ifndef SEARCH_NAME_LINEARLY
    for (i=0; i<HASH_TABLE_SIZE; i++) {
        /* initialize dim name lookup table */
        ncp->dims.nameT[i].num = 0;
        ncp->dims.nameT[i].list = NULL;
        /* initialize var name lookup table */
        ncp->vars.nameT[i].num = 0;
        ncp->vars.nameT[i].list = NULL;
    }
#endif

    /* set the file format version based on the create mode, cmode */
    if (fIsSet(cmode, NC_64BIT_DATA)) {
        fSet(ncp->flags, NC_64BIT_DATA);
        ncp->format = 5;
    } else if (fIsSet(cmode, NC_64BIT_OFFSET)) {
        fSet(ncp->flags, NC_64BIT_OFFSET);
        ncp->format = 2;
    } else {
        if (default_format == NC_FORMAT_CDF5) {
            fSet(ncp->flags, NC_64BIT_DATA);
            ncp->format = 5;
        }
        else if (default_format == NC_FORMAT_CDF2) {
            fSet(ncp->flags, NC_64BIT_OFFSET);
            ncp->format = 2;
        }
        else {
            fSet(ncp->flags, NC_32BIT);
            ncp->format = 1;
        }
    }
    /* extract I/O hints from user info */
    set_pnetcdf_hints(ncp, info);

    /* create file collectively */
    err = ncmpiio_create(comm, path, cmode, info, ncp);
    if (err != NC_NOERR) { /* fatal error */
        if (err == NC_EMULTIDEFINE_OMODE) err = NC_EMULTIDEFINE_CMODE;
        ncmpii_free_NC(ncp);
        DEBUG_RETURN_ERROR(err)
    }

    *ncpp = (void*)ncp;

    return NC_NOERR;
}

/*----< ncmpii_open() >------------------------------------------------------*/
int
ncmpii_open(MPI_Comm     comm,
            const char  *path,
            int          omode,
            int          ncid,
            MPI_Info     info, /* user's info and env info combined */
            void       **ncpp)
{
    char *env_str, value[MPI_MAX_INFO_VAL];
    int i, flag, err, status=NC_NOERR, safe_mode=0, mpireturn;
    MPI_Offset chunksize=NC_DEFAULT_CHUNKSIZE;
    NC *ncp=NULL;
#ifndef SEARCH_NAME_LINEARLY
    NC_nametable *nameT;
#endif

#ifdef PNETCDF_DEBUG
    /* PNETCDF_DEBUG is set at configure time, which will be overwritten by
     * the run-time environment variable PNETCDF_SAFE_MODE */
    safe_mode = 1;
#endif
    /* If environment variable PNETCDF_SAFE_MODE is set to 1, then we perform
     * a strict consistent test, i.e. arguments used in def_dim/def_var APIs
     */
    if ((env_str = getenv("PNETCDF_SAFE_MODE")) != NULL) {
        if (*env_str == '0') safe_mode = 0;
        else                 safe_mode = 1;
        /* if PNETCDF_SAFE_MODE is set but without a value, *env_str can
         * be '\0' (null character). In this case, safe_mode is enabled */
    }

    /* path's validity and omode consistency have been checked in ncmpi_open()
     * in src/dispatchers/file.c */

    /* NC_DISKLESS is not supported yet */
    if (omode & NC_DISKLESS) DEBUG_RETURN_ERROR(NC_EINVAL_OMODE)

    /* NC_MMAP is not supported yet */
    if (omode & NC_MMAP) DEBUG_RETURN_ERROR(NC_EINVAL_OMODE)

    /* allocate buffer for header object NC */
    ncp = (NC*) NCI_Calloc(1, sizeof(NC));
    if (ncp == NULL) DEBUG_RETURN_ERROR(NC_ENOMEM)

    ncp->safe_mode         = safe_mode;
    ncp->ncid              = ncid;
    ncp->abuf              = NULL;
    ncp->old               = NULL;
    ncp->chunk             = NC_DEFAULT_CHUNKSIZE;
    /* initialize arrays storing pending non-blocking requests */
    ncp->numGetReqs        = 0;
    ncp->numPutReqs        = 0;
    ncp->get_list          = NULL;
    ncp->put_list          = NULL;

    /* PnetCDF default fill mode is no fill */
    fSet(ncp->flags, NC_NOFILL);

    /* extract I/O hints from user info */
    set_pnetcdf_hints(ncp, info);

    /* open the file in parallel */
    err = ncmpiio_open(comm, path, omode, info, ncp);
    if (err != NC_NOERR) { /* fatal error */
        ncmpii_free_NC(ncp);
        if (status == NC_NOERR) status = err;
        DEBUG_RETURN_ERROR(status) /* return the 1st error encountered */
    }

    /* read header from file into an NC object pointed by ncp */
    err = ncmpii_hdr_get_NC(ncp);
    if (err != NC_NOERR) { /* fatal error */
        ncmpiio_close(ncp->nciop, 0);
        ncmpii_free_NC(ncp);
        if (status == NC_NOERR) status = err;
        DEBUG_RETURN_ERROR(status) /* return the 1st error encountered */
    }

#ifdef ENABLE_SUBFILING
    if (ncp->subfile_mode) {
        /* check subfiling attribute */
        err = ncmpii_get_att(ncp, NC_GLOBAL, "num_subfiles",
                             &ncp->nc_num_subfiles, NC_INT);
        if (err == NC_NOERR && ncp->nc_num_subfiles > 1) {
            /* ignore error NC_ENOTATT if this attribute is not defined */
            int nvars;

            err = ncmpii_inq(ncp, NULL, &nvars, NULL, NULL);
            if (status == NC_NOERR) status = err;

            for (i=0; i<nvars; i++) {
                err = ncmpii_get_att(ncp, i, "num_subfiles",
                                     &ncp->vars.value[i]->num_subfiles, NC_INT);
                if (err == NC_ENOTATT) continue;
                if (err != NC_NOERR && status == NC_NOERR) { /* other error */
                    status = err;
                    continue;
                }

                if (ncp->vars.value[i]->num_subfiles > 1) {
                    err = ncmpii_get_att(ncp, i, "ndims_org",
                                         &ncp->vars.value[i]->ndims_org, NC_INT);
                    if (status == NC_NOERR) status = err;
                }
            }

            if (ncp->nc_num_subfiles > 1) {
                err = ncmpii_subfile_open(ncp, &ncp->ncid_sf);
                if (status == NC_NOERR) status = err;
            }
        }
    }
    else
        ncp->nc_num_subfiles = 0;
#endif

    /* update the total number of record variables */
    ncp->vars.num_rec_vars = 0;
    for (i=0; i<ncp->vars.ndefined; i++)
        ncp->vars.num_rec_vars += IS_RECVAR(ncp->vars.value[i]);

#ifndef SEARCH_NAME_LINEARLY
    /* initialize dim name lookup table */
    nameT = ncp->dims.nameT;
    memset(nameT, 0, sizeof(NC_nametable) * HASH_TABLE_SIZE);

    /* populate dim name lookup table */
    for (i=0; i<ncp->dims.ndefined; i++) {
        /* hash the dim name into a key for name lookup */
        int key = HASH_FUNC(ncp->dims.value[i]->name->cp);
        nameT = &ncp->dims.nameT[key];
        if (nameT->num % NC_NAME_TABLE_CHUNK == 0)
            nameT->list = (int*) NCI_Realloc(nameT->list,
                          (size_t)(nameT->num+NC_NAME_TABLE_CHUNK) * SIZEOF_INT);
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
                          (size_t)(nameT->num+NC_NAME_TABLE_CHUNK) * SIZEOF_INT);
        nameT->list[nameT->num] = i;
        nameT->num++;
    }
#endif

    *ncpp = (void*)ncp;

    return status;
}

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
    if (NC_indep(ncp))
        ncmpiio_end_indep_data(ncp);

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

    fSet(ncp->flags, NC_INDEP);

    err = ncmpii_check_mpifh(ncp, 0);

    return err;
}

/*----< ncmpii_end_indep_data() >--------------------------------------------*/
/* This is a collective subroutine. */
int
ncmpii_end_indep_data(void *ncdp)
{
    NC *ncp = (NC*)ncdp;

    if (!NC_indep(ncp)) /* must be in independent data mode */
        DEBUG_RETURN_ERROR(NC_ENOTINDEP)

    return ncmpiio_end_indep_data(ncp);
}

/*----< ncmpiio_end_indep_data() >-------------------------------------------*/
/* this function is called when:
 * 1. ncmpi_end_indep_data()
 * 2. ncmpi_redef() from independent data mode entering to define more
 * 3. ncmpii_close() when closing the file
 * This function is collective.
 */
int
ncmpiio_end_indep_data(NC *ncp)
{
    int status=NC_NOERR;

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
        if (NC_doFsync(ncp) && NC_independentFhOpened(ncp->nciop)) {
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
            status = ncmpiio_end_indep_data(ncp); /* sync header */
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

/*----< ncmpi_delete() >-----------------------------------------------------*/
/* doesn't do anything to release resources, so call ncmpi_close before calling
 * this function.
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

/* End Of Dataset Functions */

/*----< ncmpii_check_mpifh() >-----------------------------------------------*/
int
ncmpii_check_mpifh(NC  *ncp,
                   int  collective)
{
    int mpireturn;

    if (collective && NC_indep(ncp)) /* collective handle but in indep mode */
        DEBUG_RETURN_ERROR(NC_EINDEP)

    if (!collective && !NC_indep(ncp)) /* indep handle but in collective mode */
        DEBUG_RETURN_ERROR(NC_ENOTINDEP)

    if (collective && !NC_collectiveFhOpened(ncp->nciop)) {
        TRACE_IO(MPI_File_open)(ncp->nciop->comm, (char*)ncp->nciop->path,
                                ncp->nciop->mpiomode, ncp->nciop->mpiinfo,
                                &ncp->nciop->collective_fh);
        if (mpireturn != MPI_SUCCESS)
            return ncmpii_handle_error(mpireturn, "MPI_File_open");

        set_NC_collectiveFh(ncp->nciop);
    }
    else if (!collective && !NC_independentFhOpened(ncp->nciop)) {
        TRACE_IO(MPI_File_open)(MPI_COMM_SELF, (char*)ncp->nciop->path,
                                ncp->nciop->mpiomode, ncp->nciop->mpiinfo,
                                &ncp->nciop->independent_fh);
        if (mpireturn != MPI_SUCCESS)
            return ncmpii_handle_error(mpireturn, "MPI_File_open");

        set_NC_independentFh(ncp->nciop);
    }

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

