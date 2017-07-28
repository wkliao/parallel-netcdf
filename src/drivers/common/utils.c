/*
 *  Copyright (C) 2014, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* $Id$ */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include <mpi.h>

#include <pnetcdf.h>
#include <dispatch.h>
#include <pnc_debug.h>
#include <common.h>

/*----< ncmpii_nc2mpitype() >------------------------------------------------*/
MPI_Datatype
ncmpii_nc2mpitype(nc_type xtype)
{
    switch(xtype){
        case NC_CHAR :   return MPI_CHAR;
        case NC_BYTE :   return MPI_SIGNED_CHAR;
        case NC_SHORT :  return MPI_SHORT;
        case NC_INT :    return MPI_INT;
        case NC_FLOAT :  return MPI_FLOAT;
        case NC_DOUBLE : return MPI_DOUBLE;
        case NC_UBYTE :  return MPI_UNSIGNED_CHAR;
        case NC_USHORT : return MPI_UNSIGNED_SHORT;
        case NC_UINT :   return MPI_UNSIGNED;
        case NC_INT64 :  return MPI_LONG_LONG_INT;
        case NC_UINT64 : return MPI_UNSIGNED_LONG_LONG;
        default:         return MPI_DATATYPE_NULL;
    }
}

