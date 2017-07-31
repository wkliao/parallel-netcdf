/*
 *  Copyright (C) 2017, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* $Id$ */

/*
 * This file implements the following PnetCDF APIs
 *
 * ncmpi_create()           : dispatcher->create()
 * ncmpi_open()             : dispatcher->open()
 * ncmpi_enddef()           : dispatcher->enddef()
 * ncmpi__enddef()          : dispatcher->_enddef()
 * ncmpi_redef()            : dispatcher->redef()
 * ncmpi_begin_indep_data() : dispatcher->begin_indep_data()
 * ncmpi_end_indep_data()   : dispatcher->end_indep_data()
 * ncmpi_abort()            : dispatcher->abort()
 * ncmpi_inq()              : dispatcher->inq()
 * ncmpi_inq_misc()         : dispatcher->inq_misc()
 * ncmpi_wait()             : dispatcher->wait()
 * ncmpi_wait_all()         : dispatcher->wait()
 * ncmpi_cancel()           : dispatcher->cancel()
 *
 * ncmpi_set_fill()         : dispatcher->set_fill()
 * ncmpi_fill_var_rec()     : dispatcher->fill_rec()
 * ncmpi_def_var_fill()     : dispatcher->def_var_fill()
 * ncmpi_inq_var_fill()     : dispatcher->inq()
 *
 * ncmpi_sync()             : dispatcher->sync()
 * ncmpi_sync_numrecs()     : dispatcher->sync_numrecs()
 *
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>

#include <mpi.h>
#include <pnc_debug.h>
#include <common.h>
#include <foo_driver.h>

int
foo_create(MPI_Comm     comm,
           const char  *path,
           int          cmode,
           int          ncid,
           MPI_Info     info,
           void       **ncpp)
{
    return NC_NOERR;
}

int
foo_open(MPI_Comm     comm,
         const char  *path,
         int          omode,
         int          ncid,
         MPI_Info     info,
         void       **ncpp)
{
    return NC_NOERR;
}

int
foo_enddef(void *ncdp)
{
    return NC_NOERR;
}

int
foo__enddef(void       *ncdp,
            MPI_Offset  h_minfree,
            MPI_Offset  v_align,
            MPI_Offset  v_minfree,
            MPI_Offset  r_align)
{
    return NC_NOERR;
}

int
foo_redef(void *ncdp)
{
    return NC_NOERR;
}

int
foo_begin_indep_data(void *ncdp)
{
    return NC_NOERR;
}

int
foo_end_indep_data(void *ncdp)
{
    return NC_NOERR;
}

int
foo_abort(void *ncdp)
{
    return NC_NOERR;
}

int
foo_inq(void *ncdp,
        int  *ndimsp,
        int  *nvarsp,
        int  *nattsp,
        int  *xtendimp)
{
    return NC_NOERR;
}

int
foo_inq_misc(void       *ncdp,
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
    return NC_NOERR;
}

int
foo_cancel(void *ncdp,
           int   num_req,
           int  *req_ids,
           int  *statuses)
{
    return NC_NOERR;
}

int
foo_wait(void *ncdp,
         int   num_reqs,
         int  *req_ids,
         int  *statuses,
         int   reqMode)
{
    return NC_NOERR;
}

int
foo_set_fill(void *ncdp,
             int   fill_mode,
             int  *old_fill_mode)
{
    return NC_NOERR;
}

int
foo_fill_var_rec(void      *ncdp,
                 int        varid,
                 MPI_Offset recno)
{
    return NC_NOERR;
}

int
foo_def_var_fill(void       *ncdp,
                 int         varid,
                 int         no_fill,
                 const void *fill_value)
{
    return NC_NOERR;
}

int
foo_sync_numrecs(void *ncdp)
{
    return NC_NOERR;
}

int
foo_sync(void *ncdp)
{
    return NC_NOERR;
}

