/*
 *  Copyright (C) 2015, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* $Id$ */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <stdio.h>
#ifdef HAVE_STDLIB_H
#include <stdlib.h>
#endif
#include <string.h>
#include <strings.h>  /* strcasecmp() */
#include <ctype.h>
#include <assert.h>
#include <mpi.h>

#include <pnc_debug.h>
#include <common.h>
#include "nc.h"
#include "ncx.h"

/*----< ncmpio_sanity_check() >----------------------------------------------*/
/* check the following errors and in that precedence.
 * NC_EBADID, NC_EPERM, NC_EINDEFINE, NC_EINDEP/NC_ENOTINDEP, NC_ENOTVAR,
 * NC_ECHAR, NC_EINVALCOORDS, NC_EEDGE, NC_ESTRIDE, NC_EINVAL.
 */
int ncmpio_sanity_check(NC                *ncp,
                        int               varid,
                        const MPI_Offset  bufcount,
                        MPI_Datatype      buftype,  /* internal datatype */
                        int               reqMode,
                        NC_var          **varp)  /* OUT */
{
    /* all errors detected here are fatal, must return immediately */
    int err=NC_NOERR;

    /* check file write permission if this is write request */
    if (fIsSet(reqMode, NC_REQ_WR) && NC_readonly(ncp))
        DEBUG_RETURN_ERROR(NC_EPERM)

    if (fIsSet(reqMode, NC_REQ_BLK)) {
        /* blocking APIs must be called in data mode */
        if (NC_indef(ncp))
            DEBUG_RETURN_ERROR(NC_EINDEFINE)

        /* for blocking APIs, check if in the right collective or independent
         * mode, nonblocking APIs can be called in either mode */
        if (fIsSet(reqMode, NC_REQ_INDEP) && !NC_indep(ncp))
            DEBUG_RETURN_ERROR(NC_ENOTINDEP)
        else if (fIsSet(reqMode, NC_REQ_COLL) && NC_indep(ncp))
            DEBUG_RETURN_ERROR(NC_EINDEP)
    }

    if (fIsSet(reqMode, NC_REQ_ZERO)) return NC_NOERR;

    /* check if varid is valid (check NC_ENOTVAR) */
    err = ncmpio_NC_lookupvar(ncp, varid, varp);
    if (err != NC_NOERR) return err;

    /* check NC_ECHAR */
    if (fIsSet(reqMode, NC_REQ_FLEX)) {
        /* when buftype == MPI_DATATYPE_NULL, bufcount is ignored and this API
         * assumes argument buf's data type matches the data type of variable
         * defined in the file - no data conversion will be done.
         */
        if (buftype != MPI_DATATYPE_NULL) {
            int isderived, el_size, buftype_is_contig;
            MPI_Datatype ptype;
            MPI_Offset   bnelems=0;

            err = ncmpii_dtype_decode(buftype, &ptype, &el_size, &bnelems,
                                      &isderived, &buftype_is_contig);
            if (err != NC_NOERR) return err;

            err = NCMPII_ECHAR((*varp)->type, ptype);
            if (err != NC_NOERR) return err;
        }
        /* else case: itype matches xtype */

        /* for flexible APIs, bufcount cannot be negative */
        if (bufcount < 0) DEBUG_RETURN_ERROR(NC_EINVAL)
    }
    else { /* called from a high-level API */
        err = NCMPII_ECHAR((*varp)->type, buftype);
        if (err != NC_NOERR) return err;
    }
    return NC_NOERR;
}

/*----< ncmpio_set_pnetcdf_hints() >-----------------------------------------*/
/* this is where the I/O hints designated to pnetcdf are extracted */
void ncmpio_set_pnetcdf_hints(NC *ncp, MPI_Info info)
{
    char value[MPI_MAX_INFO_VAL];
    int  flag;

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

/*----< ncmpio_xlen_nc_type() >----------------------------------------------*/
/* return the length of external NC data type */
int
ncmpio_xlen_nc_type(nc_type type) {
    switch(type) {
        case NC_BYTE:
        case NC_CHAR:
        case NC_UBYTE:  return X_SIZEOF_CHAR;
        case NC_SHORT:  return X_SIZEOF_SHORT;
        case NC_USHORT: return X_SIZEOF_USHORT;
        case NC_INT:    return X_SIZEOF_INT;
        case NC_UINT:   return X_SIZEOF_UINT;
        case NC_FLOAT:  return X_SIZEOF_FLOAT;
        case NC_DOUBLE: return X_SIZEOF_DOUBLE;
        case NC_INT64:  return X_SIZEOF_INT64;
        case NC_UINT64: return X_SIZEOF_UINT64;
        default: DEBUG_RETURN_ERROR(NC_EBADTYPE);
    }
    return NC_NOERR;
}

#include "utf8proc.h"

/* There are 3 levels of UTF8 checking: 1=> (exact)validating 2=>relaxed
   and 3=>very relaxed
*/
/* Use semi-relaxed check */
#define UTF8_CHECK 2

static int
nextUTF8(const char* cp)
{
    /*  The goal here is to recognize the length of each
	multibyte utf8 character sequence and skip it.
        Again, we assume that every non-ascii character is legal.
        We can define three possible tests of decreasing correctness
        (in the sense that the least correct will allow some sequences that
        are technically illegal UTF8).
        As Regular expressions they are as follows:
        1. most correct:
            UTF8   ([\xC2-\xDF][\x80-\xBF])                       \
                 | (\xE0[\xA0-\xBF][\x80-\xBF])                   \
                 | ([\xE1-\xEC][\x80-\xBF][\x80-\xBF])            \
                 | (\xED[\x80-\x9F][\x80-\xBF])                   \
                 | ([\xEE-\xEF][\x80-\xBF][\x80-\xBF])            \
                 | (\xF0[\x90-\xBF][\x80-\xBF][\x80-\xBF])        \
                 | ([\xF1-\xF3][\x80-\xBF][\x80-\xBF][\x80-\xBF]) \
                 | (\xF4[\x80-\x8F][\x80-\xBF][\x80-\xBF])        \

        2. partially relaxed:
            UTF8 ([\xC0-\xDF][\x80-\xBF])
                 |([\xE0-\xEF][\x80-\xBF][\x80-\xBF])
                 |([\xF0-\xF7][\x80-\xBF][\x80-\xBF][\x80-\xBF])

        3. The most relaxed version of UTF8:
            UTF8 ([\xC0-\xD6].)|([\xE0-\xEF]..)|([\xF0-\xF7]...)

        We use #2 here.

	The tests are derived from the table at
	    http://www.w3.org/2005/03/23-lex-U
    */

/* Define a test macro to test against a range */
#define RANGE(c,lo,hi) (((uchar)c) >= lo && ((uchar)c) <= hi)
/* Define a common RANGE */
#define RANGE0(c) RANGE(c,0x80,0xBF)

    int ch0;

    int skip = -1; /* assume failed */

    ch0 = (uchar)*cp;
    if(ch0 <= 0x7f) skip = 1; /* remove ascii case */
    else

#if UTF8_CHECK == 2
    /* Do relaxed validation check */
    if(RANGE(ch0,0xC0,0XDF)) {/* 2-bytes, but check */
        if(cp[1] != 0 && RANGE0(cp[1]))
		skip = 2; /* two bytes */
    } else if(RANGE(ch0,0xE0,0XEF)) {/* 3-bytes, but check */
        if(cp[1] != 0 && RANGE0(cp[1]) && cp[2] != 0 && RANGE0(cp[1]))
		skip = 3; /* three bytes */
    } else if(RANGE(ch0,0xF0,0XF7)) {/* 3-bytes, but check */
        if(cp[1] != 0 && RANGE0(cp[1]) && cp[2] != 0
           && RANGE0(cp[1]) && cp[3] != 0 && RANGE0(cp[1]))
		skip = 4; /* four bytes*/
    }
#elif UTF8_CHECK == 1
    /* Do exact validation check */
    if(RANGE(ch0,0xC2,0xDF)) {/* non-overlong 2-bytes */
	int ch1 = (uchar)cp[1];
	if(ch1 != 0 && RANGE0(ch1)) skip = 2;
    } else if((ch0 == 0xE0)) {/* 3-bytes, not overlong */
	int ch1 = (uchar)cp[1];
	if(ch1 != 0 && RANGE(ch1,0xA0,0xBF)) {
	    int ch2 = (uchar)cp[2];
	    if(ch2 != 0 && RANGE0(ch2)) skip = 3;
	}
    } else if((ch0 == 0xED)) {/* 3-bytes minus surrogates */
	int ch1 = (uchar)cp[1];
	if(ch1 != 0 && RANGE(ch1,0x80,0x9f)) {
	    int ch2 = (uchar)cp[2];
	    if(ch2 != 0 && RANGE0(ch2)) skip = 3;
	}
    } else if(RANGE(ch0,0xE1,0xEC) || ch0 == 0xEE || ch0 == 0xEF) {
	int ch1 = (uchar)cp[1];
	if(ch1 != 0 && RANGE0(ch1)) {
	    int ch2 = (uchar)cp[2];
	    if(ch2 != 0 && RANGE0(ch2)) skip = 3;
	}
    } else if((ch0 == 0xF0)) {/* planes 1-3 */
	int ch1 = (uchar)cp[1];
	if(ch1 != 0 && RANGE(ch1,0x90,0xBF)) {
	    int ch2 = (uchar)cp[2];
	    if(ch2 != 0 && RANGE0(ch2)) {
	        int ch3 = (uchar)cp[3];
	        if(ch3 != 0 && RANGE0(ch3)) skip = 4;
	    }
	}
    } else if((ch0 == 0xF4)) {/* plane 16 */
	int ch1 = (uchar)cp[1];
	if(ch1 != 0 && RANGE0(ch1)) {
	    int ch2 = (uchar)cp[2];
	    if(ch2 != 0 && RANGE0(ch2)) {
	        int ch3 = (uchar)cp[3];
	        if(ch3 != 0 && RANGE0(ch3)) skip = 4;
	    }
	}
    } else if(RANGE(ch0,0xF1,0xF3)) { /* planes 4-15 */
	int ch1 = (uchar)cp[1];
	if(ch1 != 0 && RANGE0(ch1)) {
	    int ch2 = (uchar)cp[2];
	    if(ch2 != 0 && RANGE0(ch2)) {
	        int ch3 = (uchar)cp[3];
	        if(ch3 != 0 && RANGE0(ch3)) skip = 4;
	    }
	}
    }
#else
#error "Must Define UTF8_CHECK as 1 or 2"
#endif
    return skip;
}


#ifdef _CONFORM_NETCDF_3_5_1
/*
 * For CDF-1, Verify that a name string is valid
 * CDL syntax, eg, all the characters are
 * alphanumeric, '-', '_', or '.'.
 * Also permit ':', '@', '(', or ')' in names for chemists currently making
 * use of these characters, but don't document until ncgen and ncdump can
 * also handle these characters in names.
 */
static int
check_name_CDF1(const char *name)
{
    const char *cp = name;
    assert(name != NULL);

    if (*name == 0)
        DEBUG_RETURN_ERROR(NC_EBADNAME) /* empty names disallowed */

    for (; *cp != 0; cp++) {
        int ch = *cp;
        if (!isalnum(ch)) {
            if (ch != '_' && ch != '-' && ch != '+' && ch != '.' &&
                ch != ':' && ch != '@' && ch != '(' && ch != ')')
                DEBUG_RETURN_ERROR(NC_EBADNAME)
        }
    }
    if (cp - name > NC_MAX_NAME)
        DEBUG_RETURN_ERROR(NC_EMAXNAME)

    return NC_NOERR;
}
#endif

/*
 * Verify that a name string is valid syntax.  The allowed name
 * syntax (in RE form) is:
 *
 * ([a-zA-Z_]|{UTF8})([^\x00-\x1F\x7F/]|{UTF8})*
 *
 * where UTF8 represents a multibyte UTF-8 encoding.  Also, no
 * trailing spaces are permitted in names.  This definition
 * must be consistent with the one in ncgen.l.  We do not allow '/'
 * because HDF5 does not permit slashes in names as slash is used as a
 * group separator.  If UTF-8 is supported, then a multi-byte UTF-8
 * character can occur anywhere within an identifier.  We later
 * normalize UTF-8 strings to NFC to facilitate matching and queries.
 */
static int
check_name_CDF2(const char *name)
{
	int skip;
	int ch;
	const char *cp = name;
	ssize_t utf8_stat;

	assert(name != NULL);

	if(*name == 0		/* empty names disallowed */
	   || strchr(cp, '/'))	/* '/' can't be in a name */
		DEBUG_RETURN_ERROR(NC_EBADNAME)

	/* check validity of any UTF-8 */
	utf8_stat = ncmpii_utf8proc_check((const unsigned char *)name);
	if (utf8_stat < 0)
	    DEBUG_RETURN_ERROR(NC_EBADNAME)

	/* First char must be [a-z][A-Z][0-9]_ | UTF8 */
	ch = (uchar)*cp;
	if(ch <= 0x7f) {
	    if(!('A' <= ch && ch <= 'Z')
	       && !('a' <= ch && ch <= 'z')
               && !('0' <= ch && ch <= '9')
	       && ch != '_' )
		DEBUG_RETURN_ERROR(NC_EBADNAME)
	    cp++;
	} else {
	    if((skip = nextUTF8(cp)) < 0)
		DEBUG_RETURN_ERROR(NC_EBADNAME)
	    cp += skip;
	}

	while(*cp != 0) {
	    ch = (uchar)*cp;
	    /* handle simple 0x00-0x7f characters here */
	    if(ch <= 0x7f) {
                if( ch < ' ' || ch > 0x7E) /* control char or DEL */
		  DEBUG_RETURN_ERROR(NC_EBADNAME)
		cp++;
	    } else {
		if((skip = nextUTF8(cp)) < 0) DEBUG_RETURN_ERROR(NC_EBADNAME)
		cp += skip;
	    }
	    if(cp - name > NC_MAX_NAME)
		DEBUG_RETURN_ERROR(NC_EMAXNAME)
	}
	if(ch <= 0x7f && isspace(ch)) /* trailing spaces disallowed */
	    DEBUG_RETURN_ERROR(NC_EBADNAME)
	return NC_NOERR;
}

/*----< ncmpio_NC_check_name() >---------------------------------------------*/
int
ncmpio_NC_check_name(const char *name,
                     int         file_ver) /* CDF version: 1, 2, or 5 */
{
    /* NetCDF4 has made CDF-1 no different from CDF-2 except the size of
     * OFFSET (i.e. 32-bit vs. 64-bit integer. Both formats support extended
     * names now.
     */
#ifdef _CONFORM_NETCDF_3_5_1
    if (file_ver == 1)
        return check_name_CDF1(name);
#endif

    return check_name_CDF2(name);
}

/*----< ncmpio_last_offset() >-----------------------------------------------*/
/* returns the file offset of the last byte accessed of this request
 * If counts is NULL, this is equivalent to the starting offset of this
 * request
 */
int
ncmpio_last_offset(NC               *ncp,
                   NC_var           *varp,
                   const MPI_Offset  starts[],   /* [varp->ndims] */
                   const MPI_Offset  counts[],   /* [varp->ndims] */
                   const MPI_Offset  strides[],  /* [varp->ndims] */
                   const int         reqMode,
                   MPI_Offset       *offset_ptr) /* return file offset */
{
    MPI_Offset offset, *end_off=NULL;
    int status, i, ndims;

    offset = varp->begin; /* beginning file offset of this variable */
    ndims  = varp->ndims; /* number of dimensions of this variable */

    if (counts != NULL) {
        end_off = (MPI_Offset*) NCI_Malloc((size_t)ndims * SIZEOF_MPI_OFFSET);

        if (strides != NULL) {
            for (i=0; i<ndims; i++)
                end_off[i] = starts[i] + (counts[i] - 1) * strides[i];
        }
        else { /* strides == NULL */
            for (i=0; i<ndims; i++)
                end_off[i] = starts[i] + counts[i] - 1;
        }
    }
    else { /* when counts == NULL strides is of no use */
        end_off = (MPI_Offset*) starts;
    }

    /* check whether end_off is valid */
    status = ncmpii_start_count_stride_check(ncp->format, API_VAR1, varp->ndims,
                                             ncp->numrecs, varp->shape,
                                             end_off, NULL, NULL, reqMode);
    if (status != NC_NOERR) {
#ifdef CDEBUG
        printf("%s(): ncmpii_start_count_stride_check() fails\n",__func__);
#endif
        if (end_off != NULL && end_off != starts) NCI_Free(end_off);
        return status;
    }

    if (ndims > 0) {
        if (IS_RECVAR(varp))
            /* no need to check recsize here: if MPI_Offset is only 32 bits we
               will have had problems long before here */
            offset += end_off[0] * ncp->recsize;
        else
            offset += end_off[ndims-1] * varp->xsz;

        if (ndims > 1) {
            if (IS_RECVAR(varp))
                offset += end_off[ndims - 1] * varp->xsz;
            else
                offset += end_off[0] * varp->dsizes[1] * varp->xsz;

            for (i=1; i<ndims-1; i++)
                offset += end_off[i] * varp->dsizes[i+1] * varp->xsz;
        }
    }
    if (counts != NULL && end_off != NULL)
        NCI_Free(end_off);

    *offset_ptr = offset;
    return NC_NOERR;
}

