/*
 *  Copyright (C) 2017, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* $Id$ */

/*
 * This file implements the following PnetCDF APIs.
 *
 * ncmpi_inq_attname() : dispatcher->inq_attname()
 * ncmpi_inq_attid()   : dispatcher->inq_attid()
 * ncmpi_inq_att()     : dispatcher->inq_att()
 * ncmpi_rename_att()  : dispatcher->inq_rename_att()
 * ncmpi_copy_att()    : dispatcher->inq_copy_att()
 * ncmpi_del_att()     : dispatcher->inq_del_att()
 * ncmpi_get_att()     : dispatcher->inq_get_att()
 * ncmpi_put_att()     : dispatcher->inq_put_arr()
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
#include <ncfoo_driver.h>

int
ncfoo_inq_attname(void *ncdp,
                  int   varid,
                  int   attid,
                  char *name)
{
    return NC_NOERR;
}

int
ncfoo_inq_attid(void       *ncdp,
                int         varid,
                const char *name,
                int        *attidp)
{
    return NC_NOERR;
}

int
ncfoo_inq_att(void       *ncdp,
              int         varid,
              const char *name,
              nc_type    *datatypep,
              MPI_Offset *lenp)
{
    return NC_NOERR;
}

int
ncfoo_rename_att(void       *ncdp,
                 int         varid,
                 const char *name,
                 const char *newname)
{
    return NC_NOERR;
}


int
ncfoo_copy_att(void       *ncdp_in,
               int         varid_in,
               const char *name,
               void       *ncdp_out,
               int         varid_out)
{
    return NC_NOERR;
}

int
ncfoo_del_att(void       *ncdp,
              int         varid,
              const char *name)
{
    return NC_NOERR;
}

int
ncfoo_get_att(void         *ncdp,
              int           varid,
              const char   *name,
              void         *buf,
              MPI_Datatype  itype)
{
    return NC_NOERR;
}

int
ncfoo_put_att(void         *ncdp,
              int           varid,
              const char   *name,
              nc_type       xtype,
              MPI_Offset    nelems,
              const void   *buf,
              MPI_Datatype  itype)
{
    return NC_NOERR;
}
