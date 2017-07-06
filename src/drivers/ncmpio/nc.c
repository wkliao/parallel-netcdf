/*
 *  Copyright (C) 2003, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* $Id$ */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#ifdef HAVE_STDLIB_H
#include <stdlib.h>
#endif
#include <stdio.h>
#include <string.h>
#include <assert.h>

#include <mpi.h>

#include <pnc_debug.h>
#include <common.h>
#include "nc.h"
#include "rnd.h"
#include "ncx.h"
#ifdef ENABLE_SUBFILING
#include "subfile.h"
#endif

/* These have to do with version numbers. */
#define MAGIC_NUM_LEN 4
#define VER_CLASSIC 1
#define VER_64BIT_OFFSET 2
#define VER_HDF5 3
#define VER_64BIT_DATA 5

#ifdef _CHECK_HEADER_IN_DETAIL
/*----< NC_check_header() >--------------------------------------------------*/
/*
 * Check the consistency of defined header metadata across all processes and
 * overwrite the local header objects with root's if inconsistency is found.
 * This function is collective.
 */
static int
NC_check_header(NC         *ncp,
                void       *buf,
                MPI_Offset  local_xsz) /* size of buf */
{
    int h_size, rank, g_status, status=NC_NOERR, mpireturn;

    /* root's header size has been broadcasted in NC_begin() and saved in
     * ncp->xsz.
     */

    /* TODO: When root process 0 broadcasts its header,
     * currently the header size cannot be larger than 2^31 bytes,
     * due to the 2nd argument, count, of MPI_Bcast being of type int.
     * Possible solution is to broadcast in chunks of 2^31 bytes.
     */
    h_size = (int)ncp->xsz;
    if (ncp->xsz != h_size)
        DEBUG_RETURN_ERROR(NC_EINTOVERFLOW)

    MPI_Comm_rank(ncp->nciop->comm, &rank);

    if (rank == 0) {
        TRACE_COMM(MPI_Bcast)(buf, h_size, MPI_BYTE, 0, ncp->nciop->comm);
    }
    else {
        bufferinfo gbp;
        void *cmpbuf = (void*) NCI_Malloc((size_t)h_size);

        TRACE_COMM(MPI_Bcast)(cmpbuf, h_size, MPI_BYTE, 0, ncp->nciop->comm);

        if (h_size != local_xsz || memcmp(buf, cmpbuf, h_size)) {
            /* now part of this process's header is not consistent with root's
             * check and report the inconsistent part
             */

            /* Note that gbp.nciop and gbp.offset below will not be used in
             * ncmpii_hdr_check_NC() */
            gbp.nciop  = ncp->nciop;
            gbp.offset = 0;
            gbp.size   = h_size;   /* entire header is in the buffer, cmpbuf */
            gbp.index  = 0;
            gbp.pos    = gbp.base = cmpbuf;

            /* find the inconsistent part of the header, report the difference,
             * and overwrite the local header object with root's.
             * ncmpii_hdr_check_NC() should not have any MPI communication
             * calls.
             */
            status = ncmpii_hdr_check_NC(&gbp, ncp);

            /* header consistency is only checked on non-root processes. The
             * returned status can be a fatal error or header inconsistency
             * error, (fatal errors are due to object allocation), but never
             * NC_NOERR.
             */
        }
        NCI_Free(cmpbuf);
    }

    if (ncp->safe_mode) {
        TRACE_COMM(MPI_Allreduce)(&status, &g_status, 1, MPI_INT, MPI_MIN,
                                  ncp->nciop->comm);
        if (mpireturn != MPI_SUCCESS) {
            return ncmpii_handle_error(mpireturn, "MPI_Allreduce"); 
        }
        if (g_status != NC_NOERR) { /* some headers are inconsistent */
            if (status == NC_NOERR) DEBUG_ASSIGN_ERROR(status, NC_EMULTIDEFINE)
        }
    }

    return status;
}
#endif


#if 0
/* 'defined but not used': seems like a useful function though. why did we
 * write it?  should we be using it? */

static int
NC_check_def(MPI_Comm comm, void *buf, MPI_Offset nn) {
  int rank;
  int errcheck;
  MPI_Offset compare = 0;
  void *cmpbuf;
  MPI_Offset max_size;

  MPI_Comm_rank(comm, &rank);

  if (rank == 0)
    max_size = nn;
  MPI_Bcast(&max_size, 1, MPI_OFFSET, 0, comm);

  compare = max_size - nn;

  MPI_Allreduce(&compare, &errcheck, 1, MPI_OFFSET, MPI_LOR, comm);

  if (errcheck)
    DEBUG_RETURN_ERROR(NC_EMULTIDEFINE)

  if (rank == 0)
    cmpbuf = buf;
  else
    cmpbuf = (void *)NCI_Malloc(nn);

  MPI_Bcast(cmpbuf, nn, MPI_BYTE, 0, comm);

  if (rank != 0) {
    compare = memcmp(buf, cmpbuf, nn);
    NCI_Free(cmpbuf);
  }

  MPI_Allreduce(&compare, &errcheck, 1, MPI_OFFSET, MPI_LOR, comm);

  if (errcheck){
    DEBUG_RETURN_ERROR(NC_EMULTIDEFINE)
  }else{
    return NC_NOERR;
  }
}
#endif

/*----< ncmpii_free_NC() >----------------------------------------------------*/
inline void
ncmpii_free_NC(NC *ncp)
{
    if (ncp == NULL) return;

    ncmpii_free_NC_dimarray(&ncp->dims);
    ncmpii_free_NC_attrarray(&ncp->attrs);
    ncmpii_free_NC_vararray(&ncp->vars);

    NCI_Free(ncp);
}

/*----< ncmpii_dup_NC() >----------------------------------------------------*/
NC *
ncmpii_dup_NC(const NC *ref)
{
    NC *ncp;

    ncp = (NC *) NCI_Calloc(1, sizeof(NC));
    if (ncp == NULL) return NULL;

    if (ncmpii_dup_NC_dimarray(&ncp->dims,   &ref->dims)  != NC_NOERR ||
        ncmpii_dup_NC_attrarray(&ncp->attrs, &ref->attrs) != NC_NOERR ||
        ncmpii_dup_NC_vararray(&ncp->vars,   &ref->vars)  != NC_NOERR) {
        ncmpii_free_NC(ncp);
        return NULL;
    }
    ncp->xsz       = ref->xsz;
    ncp->begin_var = ref->begin_var;
    ncp->begin_rec = ref->begin_rec;
    ncp->recsize   = ref->recsize;

    NC_set_numrecs(ncp, NC_get_numrecs(ref));
    return ncp;
}


/*
 *  Verify that this is a user nc_type
 * Formerly
NCcktype()
 * Sense of the return is changed.
 */
inline int
ncmpii_cktype(int     cdf_ver,
              nc_type type)
{
    /* the max data type supported by CDF-5 is NC_UINT64 */
    if (type <= 0 || type > NC_UINT64)
        DEBUG_RETURN_ERROR(NC_EBADTYPE)

    /* For CDF-1 and CDF-2 files, only classic types are allowed. */
    if (cdf_ver < 5 && type > NC_DOUBLE)
        DEBUG_RETURN_ERROR(NC_ESTRICTCDF2)

    return NC_NOERR;
}


/*
 * How many objects of 'type'
 * will fit into xbufsize?
 */
inline MPI_Offset
ncmpix_howmany(nc_type type, MPI_Offset xbufsize)
{
    switch(type){
        case NC_BYTE:
        case NC_UBYTE:
        case NC_CHAR:   return xbufsize;
        case NC_SHORT:  return xbufsize/X_SIZEOF_SHORT;
        case NC_USHORT: return xbufsize/X_SIZEOF_USHORT;
        case NC_INT:    return xbufsize/X_SIZEOF_INT;
        case NC_UINT:   return xbufsize/X_SIZEOF_UINT;
        case NC_FLOAT:  return xbufsize/X_SIZEOF_FLOAT;
        case NC_DOUBLE: return xbufsize/X_SIZEOF_DOUBLE;
        case NC_INT64:  return xbufsize/X_SIZEOF_INT64;
        case NC_UINT64: return xbufsize/X_SIZEOF_UINT64;
        default:
                assert("ncmpix_howmany: Bad type" == 0);
                return(0);
    }
}

#define NC_NUMRECS_OFFSET 4

/*----< ncmpii_write_numrecs() >-----------------------------------------------*/
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

/*
 * Read in the header
 * It is expensive.
 */

inline int
ncmpii_read_NC(NC *ncp) {
  int status = NC_NOERR;

  ncmpii_free_NC_dimarray(&ncp->dims);
  ncmpii_free_NC_attrarray(&ncp->attrs);
  ncmpii_free_NC_vararray(&ncp->vars);

  status = ncmpii_hdr_get_NC(ncp);

  if (status == NC_NOERR)
      fClr(ncp->flags, NC_NDIRTY);

  return status;
}

inline int
ncmpii_dset_has_recvars(NC *ncp)
{
    /* possible further optimization: set a flag on the header data
     * structure when record variable created so we can skip this loop*/
    int i;
    NC_var **vpp;

    vpp = ncp->vars.value;
    for (i=0; i< ncp->vars.ndefined; i++, vpp++) {
        if (IS_RECVAR(*vpp)) return 1;
    }
    return 0;
}

/*
 * Given a valid ncp, check all variables for their sizes against the maximal
 * allowable sizes. Different CDF formation versions have different maximal
 * sizes. This function returns NC_EVARSIZE if any variable has a bad len
 * (product of non-rec dim sizes too large), else return NC_NOERR.
 */
int
ncmpii_NC_check_vlens(NC *ncp)
{
    NC_var **vpp;
    /* maximum permitted variable size (or size of one record's worth
       of a record variable) in bytes.  This is different for format 1
       and format 2. */
    MPI_Offset ii, vlen_max, rec_vars_count;
    MPI_Offset large_fix_vars_count, large_rec_vars_count;
    int last = 0;

    if (ncp->vars.ndefined == 0)
        return NC_NOERR;

    if (ncp->format >= 5) /* CDF-5 */
        return NC_NOERR;

    /* only CDF-1 and CDF-2 need to continue */

    if (ncp->flags & NC_64BIT_OFFSET) /* CDF2 format */
        vlen_max = X_UINT_MAX - 3; /* "- 3" handles rounded-up size */
    else
        vlen_max = X_INT_MAX - 3; /* CDF1 format */

    /* Loop through vars, first pass is for non-record variables */
    large_fix_vars_count = 0;
    rec_vars_count = 0;
    vpp = ncp->vars.value;
    for (ii = 0; ii < ncp->vars.ndefined; ii++, vpp++) {
        if (!IS_RECVAR(*vpp)) {
            last = 0;
            if (ncmpii_NC_check_vlen(*vpp, vlen_max) == 0) {
                /* check this variable's shape product against vlen_max */
                large_fix_vars_count++;
                last = 1;
            }
        } else {
            rec_vars_count++;
        }
    }
    /* OK if last non-record variable size too large, since not used to
       compute an offset */
    if (large_fix_vars_count > 1)  /* only one "too-large" variable allowed */
        DEBUG_RETURN_ERROR(NC_EVARSIZE)

    /* The only "too-large" variable must be the last one defined */
    if (large_fix_vars_count == 1 && last == 0)
        DEBUG_RETURN_ERROR(NC_EVARSIZE)

    if (rec_vars_count == 0) return NC_NOERR;

    /* if there is a "too-large" fixed-size variable, no record variable is
     * allowed */
    if (large_fix_vars_count == 1)
        DEBUG_RETURN_ERROR(NC_EVARSIZE)

    /* Loop through vars, second pass is for record variables.   */
    large_rec_vars_count = 0;
    vpp = ncp->vars.value;
    for (ii = 0; ii < ncp->vars.ndefined; ii++, vpp++) {
        if (IS_RECVAR(*vpp)) {
            last = 0;
            if (ncmpii_NC_check_vlen(*vpp, vlen_max) == 0) {
                /* check this variable's shape product against vlen_max */
                large_rec_vars_count++;
                last = 1;
            }
        }
    }

    /* For CDF-2, no record variable can require more than 2^32 - 4 bytes of
     * storage for each record's worth of data, unless it is the last record
     * variable. See
     * http://www.unidata.ucar.edu/software/netcdf/docs/file_structure_and_performance.html#offset_format_limitations
     */
    if (large_rec_vars_count > 1)  /* only one "too-large" variable allowed */
        DEBUG_RETURN_ERROR(NC_EVARSIZE)

    /* and it has to be the last one */
    if (large_rec_vars_count == 1 && last == 0)
        DEBUG_RETURN_ERROR(NC_EVARSIZE)

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
    /* ncmpii__enddef() will update nc_num_subfiles */
    /* TODO: should check ncid_sf? */
    /* if the file has subfiles, close them first */
    if (ncp->nc_num_subfiles > 1)
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
        MPI_Comm_rank(ncp->nciop->comm, &rank);
        printf("PnetCDF warning: %d nonblocking get requests still pending on process %d. Cancelling ...\n",ncp->numGetReqs,rank);
        ncmpii_cancel(ncp, NC_GET_REQ_ALL, NULL, NULL);
        if (status == NC_NOERR ) status = NC_EPENDING;
    }
    if (ncp->numPutReqs > 0) {
        int rank;
        MPI_Comm_rank(ncp->nciop->comm, &rank);
        printf("PnetCDF warning: %d nonblocking put requests still pending on process %d. Cancelling ...\n",ncp->numPutReqs,rank);
        ncmpii_cancel(ncp, NC_PUT_REQ_ALL, NULL, NULL);
        if (status == NC_NOERR ) status = NC_EPENDING;
    }
#endif

    /* If the user wants a stronger data consistency by setting NC_SHARE */
    if (fIsSet(ncp->nciop->ioflags, NC_SHARE))
        ncmpiio_sync(ncp->nciop); /* calling MPI_File_sync() */

    /* calling MPI_File_close() */
    ncmpiio_close(ncp->nciop, 0);
    ncp->nciop = NULL;

    /* free up space occupied by the header metadata */
    ncmpii_free_NC(ncp);

    return status;
}

/*----< ncmpii_inq() >--------------------------------------------------------*/
int
ncmpii_inq(void *ncdp,
           int  *ndimsp,
           int  *nvarsp,
           int  *nattsp,
           int  *xtendimp)
{
    NC *ncp = (NC*)ncdp;

    if (ndimsp != NULL)
        *ndimsp = (int) ncp->dims.ndefined;
    if (nvarsp != NULL)
        *nvarsp = (int) ncp->vars.ndefined;
    if (nattsp != NULL)
        *nattsp = (int) ncp->attrs.ndefined;
    if (xtendimp != NULL)
        /* *xtendimp = ncmpii_find_NC_Udim(&ncp->dims, NULL); */
        *xtendimp = ncp->dims.unlimited_id;

    return NC_NOERR;
}

