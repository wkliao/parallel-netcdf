dnl Process this m4 file to produce 'C' language file.
dnl
dnl If you see this line, you can ignore the next one.
/* Do not edit this file. It is produced from the corresponding .m4 source */
dnl
/*
 *  Copyright (C) 2017, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* $Id: attribute.m4 2873 2017-02-14 02:58:34Z wkliao $ */

#if HAVE_CONFIG_H
# include <ncconfig.h>
#endif

#include <stdlib.h>

#include <pnetcdf.h>
#include <dispatch.h>
#include <nctypes.h>

include(`foreach.m4')dnl
include(`utils.m4')dnl

dnl
dnl GET_ATT(fntype)
dnl
define(`GET_ATT',dnl
`dnl
/*----< ncmpi_get_att_$1() >-------------------------------------------------*/
/* This is an independent subroutine.
ifelse(`$1',`text',` * This API never returns NC_ERANGE error, as text is not convertible to numerical types')
 */
int
ncmpi_get_att_$1(int             ncid,
                 int             varid,
                 const char     *name,
                 FUNC2ITYPE($1) *buf)
{
    int err;
    PNC *pncp;
ifelse(`$1',`long',`#if SIZEOF_LONG == SIZEOF_INT
    nc_type itype=NC_INT;
#elif SIZEOF_LONG == SIZEOF_LONG_LONG
    nc_type itype=NC_INT64;
#endif',`    nc_type itype=NC_TYPE($1);')

    /* check if ncid is valid */
    err = PNC_check_id(ncid, &pncp);
    if (err != NC_NOERR) return err;

    /* calling the subroutine that implements ncmpi_get_att_$1() */
    err = pncp->dispatch->get_att(pncp->ncp, varid, name, buf, itype);
    if (err != NC_NOERR) return err;

    return NC_NOERR;
}
')dnl

foreach(`iType', (text,schar,uchar,short,ushort,int,uint,long,float,double,longlong,ulonglong),
        `GET_ATT(iType)
')

/*----< ncmpi_get_att() >----------------------------------------------------*/
/* This is an independent subroutine.
 * The user buffer data type matches the external type defined in file.
 */
int
ncmpi_get_att(int         ncid,
              int         varid,
              const char *name,
              void       *buf)
{
    int err;
    PNC *pncp;

    /* check if ncid is valid */
    err = PNC_check_id(ncid, &pncp);
    if (err != NC_NOERR) return err;

    /* calling the subroutine that implements ncmpi_get_att() */
    err = pncp->dispatch->get_att(pncp->ncp, varid, name, buf, NC_NAT);
    if (err != NC_NOERR) return err;

    return NC_NOERR;
}

dnl
dnl PUT_ATT(fntype)
dnl
define(`PUT_ATT',dnl
`dnl
/*----< ncmpi_put_att_$1() >-------------------------------------------------*/
/* This is a collective subroutine, all arguments should be consistent among
 * all processes.
 *
 * Note from netCDF user guide:
 * Attributes are always single values or one-dimensional arrays. This works
 * out well for a string, which is a one-dimensional array of ASCII characters
 *
 * Note ncmpii_put_att_text will never return NC_ERANGE error, as text is not
 * convertible to numerical types.
 */
int
ncmpi_put_att_$1(int         ncid,
                 int         varid,
                 const char *name,     /* attribute name */
                 ifelse(`$1',`text',,`nc_type xtype,')
                 MPI_Offset  nelems,   /* number of elements in buf */
                 const FUNC2ITYPE($1) *buf) /* user write buffer */
{
    int err;
    PNC *pncp;
ifelse(`$1',`long',`#if SIZEOF_LONG == SIZEOF_INT
    nc_type itype=NC_INT;
#elif SIZEOF_LONG == SIZEOF_LONG_LONG
    nc_type itype=NC_INT64;
#endif',`    nc_type itype=NC_TYPE($1);')

    /* check if ncid is valid */
    err = PNC_check_id(ncid, &pncp);
    if (err != NC_NOERR) return err;

    /* calling the subroutine that implements ncmpi_put_att_$1() */
    err = pncp->dispatch->put_att(pncp->ncp, varid, name,
                                  ifelse(`$1',`text',`NC_CHAR,',`xtype,')
                                  nelems, buf, itype);
    if (err != NC_NOERR) return err;

    return NC_NOERR;
}
')dnl

foreach(`iType', (text,schar,uchar,short,ushort,int,uint,long,float,double,longlong,ulonglong),
        `PUT_ATT(iType)
')

/*----< ncmpi_put_att() >----------------------------------------------------*/
/* This is a collective subroutine, all arguments should be consistent among
 * all processes. This API is for when the user buffer data type matches the
 * external type defined in file.
 */
int
ncmpi_put_att(int         ncid,
              int         varid,
              const char *name,
              nc_type     xtype,  /* external data type, i.e. NC_CHAR etc. */
              MPI_Offset  nelems,
              const void *buf)
{
    int err;
    PNC *pncp;

    /* check if ncid is valid */
    err = PNC_check_id(ncid, &pncp);
    if (err != NC_NOERR) return err;

    /* calling the subroutine that implements ncmpi_put_att() */
    err = pncp->dispatch->put_att(pncp->ncp, varid, name, xtype, nelems, buf, xtype);
    if (err != NC_NOERR) return err;

    return NC_NOERR;
}
