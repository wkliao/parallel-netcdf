/* -*- Mode: C; c-basic-offset:4 ; -*- */
/*  
 *  (C) 2001 by Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 *
 * This file is automatically generated by buildiface -infile=../lib/pnetcdf.h -deffile=defs
 * DO NOT EDIT
 */
#include "mpinetcdf_impl.h"


#ifdef F77_NAME_UPPER
#define nfmpi_put_var_real_ NFMPI_PUT_VAR_REAL
#elif defined(F77_NAME_LOWER_2USCORE)
#define nfmpi_put_var_real_ nfmpi_put_var_real__
#elif !defined(F77_NAME_LOWER_USCORE)
#define nfmpi_put_var_real_ nfmpi_put_var_real
/* Else leave name alone */
#endif


/* Prototypes for the Fortran interfaces */
#include "mpifnetcdf.h"
FORTRAN_API int FORT_CALL nfmpi_put_var_real_ ( int *v1, int *v2, float*v3 ){
    int ierr;
    int l2 = *v2 - 1;
    ierr = ncmpi_put_var_float( *v1, l2, v3 );
    return ierr;
}
