/*
 *  Copyright (C) 2003, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* $Id$ */

#ifndef _NC_H
#define _NC_H

/*
 * netcdf library 'private' data structures, objects and interfaces
 */

#include <stddef.h>     /* size_t */
#include <sys/types.h>  /* off_t */

#include <dispatch.h>
#include "ncmpio_dispatch.h"
#include "fbits.h"

#define FILE_ALIGNMENT_DEFAULT 512
#define HEADER_ALIGNMENT_LB    4

/* for put request less than 4KB, copy it to a buffer and do byte swap there,
 * so if the user buffer is immutable (assuming smaller than 4KB), it will not
 * cause seg fault. Not a perfect solution, but should be sufficient for most
 * of the cases.
 */
#define NC_BYTE_SWAP_BUFFER_SIZE 4096

/* define MPI_OFFSET if not defined */
#ifndef HAVE_DECL_MPI_OFFSET
    #ifdef HAVE_DECL_MPI_LONG_LONG_INT
        #define MPI_OFFSET MPI_LONG_LONG_INT
    #else
        #define MPI_OFFSET MPI_INT
    #endif
#endif

#define WRITE_REQ 0
#define READ_REQ  1

/* XXX: this seems really low.  do we end up spending a ton of time mallocing?
 * could we reduce that by increasing this to something 21st century? */
#ifndef NC_ARRAY_GROWBY
#define NC_ARRAY_GROWBY 4
#endif

/* ncmpi_create/ncmpi_open set up header to be 'chunksize' big and to grow
 * by 'chunksize' as new items added. This used to be 4k. 256k lets us read
 * in an entire climate header in one go */
#define NC_DEFAULT_CHUNKSIZE 262144

/* when variable's nctype is NC_CHAR, I/O buffer's MPI type must be MPI_CHAR
 * and vice versa */
#define NCMPII_ECHAR(nctype, mpitype) ((((nctype) == NC_CHAR) == ((mpitype) != MPI_CHAR)) ? NC_ECHAR : NC_NOERR)

/*
 * The extern size of an empty
 * netcdf version 1 file.
 * The initial value of ncp->xsz.
 */
#define MIN_NC_XSZ 32

/* netcdf file format:
     netcdf_file  = header  data
     header       = magic  numrecs  dim_list  gatt_list  var_list
     magic        = 'C'  'D'  'F'  VERSION
     VERSION      = \x01 | \x02 | \x05
     numrecs      = NON_NEG | STREAMING
     dim_list     = ABSENT | NC_DIMENSION  nelems  [dim ...]
     gatt_list    = att_list
     att_list     = ABSENT | NC_ATTRIBUTE  nelems  [attr ...]
     var_list     = ABSENT | NC_VARIABLE   nelems  [var ...]
     ABSENT       = ZERO  ZERO                  // Means list is not present
     ZERO         = \x00 \x00 \x00 \x00         // 32-bit zero

  Minimum happens when nothing is defined, i.e.
     magic              -- 4 bytes
     numrecs            -- 4 bytes for CDF-1 and CDF-2, 8 bytes for CDF-5
     dim_list = ABSENT  -- 8 bytes
     gatt_list = ABSENT -- 8 bytes
     var_list = ABSENT  -- 8 bytes
*/

typedef struct NC NC; /* forward reference */

/*
 *  The internal data types
 */
typedef enum {
    NC_UNSPECIFIED =  0,
/*  NC_BITFIELD    =  7, */
/*  NC_STRING      =  8, */
    NC_DIMENSION   = 10,  /* \x00 \x00 \x00 \x0A */
    NC_VARIABLE    = 11,  /* \x00 \x00 \x00 \x0B */
    NC_ATTRIBUTE   = 12   /* \x00 \x00 \x00 \x0C */
} NCtype;


/*
 * Counted string for names and such
 */
typedef struct {
    /* all xdr'd */
    MPI_Offset  nchars;
    char       *cp;     /* [nchars+1] one additional char for '\0' */
} NC_string;

extern NC *
ncmpio_dup_NC(const NC *ref);

/* Begin defined in string.c */
extern void
ncmpio_free_NC_string(NC_string *ncstrp);

extern int
ncmpio_NC_check_name(const char *name, int file_ver);

extern NC_string *
ncmpio_new_NC_string(size_t slen, const char *str);

extern int
ncmpio_set_NC_string(NC_string *ncstrp, const char *str);

/* End defined in string.c */

/*
 * NC dimension structure
 */
typedef struct {
    /* all xdr'd */
    MPI_Offset size;
    NC_string *name;
#ifdef ENABLE_SUBFILING
    MPI_Offset rcount; /* subfile range count */
    int range[2]; /* subfile range {start, end} */
#endif
} NC_dim;

#define NC_NAME_TABLE_CHUNK 16
#define HASH_TABLE_SIZE 256
/*
#define HASH_FUNC(x) ncmpio_jenkins_one_at_a_time_hash(x)
#define HASH_FUNC(x) ncmpio_additive_hash(x)
#define HASH_FUNC(x) ncmpio_rotating_hash(x)
#define HASH_FUNC(x) ncmpio_Pearson_hash(x)
*/
#define HASH_FUNC(x) ncmpio_Bernstein_hash(x)
/* Look like Bernstein's hashing function performs the best */

/* For the initial naive implementation of hashing:
 * #define HASH_FUNC(x) (unsigned char)x[0]
 * if used this simple hash function, HASH_TABLE_SIZE must be 256 which is the
 * number of possible keys can be stored in an unsigned char
 */

typedef struct NC_nametable {
    int  num;  /* number of elements added in the list array. Its value starts
                  with zero and incremented with new name is created. When its
                  value becomes a multiple of NC_NAME_TABLE_CHUNK, list will be
                  reallocated to a space of size (num + NC_NAME_TABLE_CHUNK) */
    int *list; /* dimension or variable IDs */
} NC_nametable;

/* the dimension ID returned from ncmpi_def_dim() is an integer pointer
 * which means the total number of defined dimension allowed in a file
 * is up to 2^31-1. Thus, the member ndefined below should be of type int.
 */
typedef struct NC_dimarray {
    int            nalloc;        /* number allocated >= ndefined */
    int            ndefined;      /* number of defined dimensions */
    int            unlimited_id;  /* -1 for not defined, otherwise >= 0 */
    NC_dim       **value;
    NC_nametable   nameT[HASH_TABLE_SIZE]; /* table for quick name lookup.
                    * indices 0, 1, ... HASH_TABLE_SIZE-1 are the hash keys.
                    * nameT[i].num is the number of hash collisions. The IDs of
                    * dimensions with names producing the same hash key i are
                    * stored in nameT[i].list[*]
                    */
} NC_dimarray;

/* Begin defined in dim.c */

extern void
ncmpio_free_NC_dim(NC_dim *dimp);

extern NC_dim *
ncmpio_new_x_NC_dim(NC_string *name);

extern int
ncmpio_find_NC_Udim(const NC_dimarray *ncap, NC_dim **dimpp);

extern int
incr_NC_dimarray(NC_dimarray *ncap, NC_dim *newdimp);

extern NC_dim*
dup_NC_dim(const NC_dim *dimp);

/* dimarray */

extern void
ncmpio_free_NC_dimarray(NC_dimarray *ncap);

extern int
ncmpio_dup_NC_dimarray(NC_dimarray *ncap, const NC_dimarray *ref);

extern NC_dim *
ncmpio_elem_NC_dimarray(const NC_dimarray *ncap, int elem);

/* End defined in dim.c */

/*
 * NC attribute
 *
 * Number of attributes is limited by 2^31-1 because the argument attnump in
 *  int nc_inq_attid(int ncid, int varid, const char *name, int *attnump);
 * is a signed 4-byte integer.
 */
typedef struct {
    nc_type    type;     /* the discriminant */
    MPI_Offset nelems;   /* number of attribute elements */
    MPI_Offset xsz;      /* amount of space at xvalue (4-byte aligned) */
    NC_string *name;     /* name of the attributes */
    void      *xvalue;   /* the actual data, in external representation */
} NC_attr;

typedef struct NC_attrarray {
    int       nalloc;    /* number allocated >= ndefined */
    int       ndefined;  /* number of defined attributes */
    NC_attr **value;
} NC_attrarray;

/* Begin defined in attr.c */

extern void
ncmpio_free_NC_attr(NC_attr *attrp);

extern NC_attr *
ncmpio_new_x_NC_attr(NC_string *strp, nc_type type, MPI_Offset nelems);

extern int
incr_NC_attrarray(NC_attrarray *ncap, NC_attr *newelemp);

extern NC_attr*
dup_NC_attr(const NC_attr *rattrp);

extern int
ncmpio_NC_findattr(const NC_attrarray *ncap, const char *uname);

/* attrarray */

extern void
ncmpio_free_NC_attrarray(NC_attrarray *ncap);

extern int
ncmpio_dup_NC_attrarray(NC_attrarray *ncap, const NC_attrarray *ref);

extern NC_attr *
ncmpio_elem_NC_attrarray(const NC_attrarray *ncap, MPI_Offset elem);

/* End defined in attr.c */

/*
 * NC variable: description and data
 */
typedef struct {
    int           varid;   /* variable ID */
    int           xsz;     /* byte size of 1 array element */
    int           ndims;   /* number of dimensions */
    nc_type       type;    /* variable's data type */
    int           no_fill; /* whether fill mode is disabled */
#ifdef ENABLE_SUBFILING
    int           num_subfiles;
    int           ndims_org;  /* ndims before subfiling */
    int          *dimids_org; /* dimids before subfiling */
#endif
    int          *dimids; /* array of dimension IDs */
    MPI_Offset   *shape;  /* dim->size of each dim
                             shape[0] == NC_UNLIMITED if record variable */
    MPI_Offset   *dsizes; /* the right to left product of shape */
    NC_string    *name;   /* name of the variable */
    MPI_Offset    len;    /* this is the "vsize" defined in header format, the
                             total size in bytes of the array variable.
                             For record variable, this is the record size */
    MPI_Offset    begin;  /* starting file offset of this variable */
    NC_attrarray  attrs;  /* attribute array */
} NC_var;

/* note: we only allow less than 2^31-1 variables defined in a file */
typedef struct NC_vararray {
    int            nalloc;      /* number allocated >= ndefined */
    int            ndefined;    /* number of defined variables */
    int            num_rec_vars;/* number of defined record variables */
    NC_var       **value;
    NC_nametable   nameT[HASH_TABLE_SIZE]; /* table for quick name lookup.
                    * indices 0, 1, ... HASH_TABLE_SIZE-1 are the hash keys.
                    * nameT[i].num is the number of hash collisions. The IDs of
                    * variables with names producing the same hash key i are
                    * stored in nameT[i].list[*]
                    */
} NC_vararray;

/* Begin defined in var.c */

extern void
ncmpio_free_NC_var(NC_var *varp);

extern NC_var *
ncmpio_new_x_NC_var(NC_string *strp, int ndims);

extern NC_var*
dup_NC_var(const NC_var *rvarp);

extern int
incr_NC_vararray(NC_vararray *ncap, NC_var *newvarp);

/* vararray */

extern void
ncmpio_free_NC_vararray(NC_vararray *ncap);

extern int
ncmpio_dup_NC_vararray(NC_vararray *ncap, const NC_vararray *ref);

extern int
ncmpio_NC_var_shape64(NC_var *varp, const NC_dimarray *dims);

extern int
ncmpio_NC_check_vlen(NC_var *varp, MPI_Offset vlen_max);

extern int
ncmpio_NC_lookupvar(NC *ncp, int varid, NC_var **varp);

/* End defined in var.c */

#define IS_RECVAR(vp) \
        ((vp)->shape != NULL ? (*(vp)->shape == NC_UNLIMITED) : 0 )

/*
 *  The PnetCDF non-blocking I/O request type
 */
typedef struct NC_req {
    int            id;          /* even number for write, odd for read */
    int            buftype_is_contig;
    int            need_swap_back_buf;
    int            abuf_index;  /* index in the abuf occupy_table
                                   -1 means not using attached buffer */
    MPI_Offset     num_recs;    /* number of records requested (1 for
                                   fixed-size variable) */
    void          *buf;         /* the original user buffer */
    void          *xbuf;        /* the buffer used to read/write, may point to
                                   the same address as buf */
    void          *tmpBuf;      /* tmp buffer to be freed, used only by
                                   nonblocking varn when buftype is noncontig */
    void          *userBuf;     /* user buffer to be unpacked from tmpBuf. used
                                   only by by nonblocking varn when buftype is
                                   noncontig */
    NC_var        *varp;
    MPI_Offset    *start;        /* [varp->ndims] */
    MPI_Offset    *count;        /* [varp->ndims] */
    MPI_Offset    *stride;       /* [varp->ndims] */
    MPI_Offset     bnelems;      /* number of elements in user buffer */
    MPI_Offset     offset_start; /* starting of aggregate access region */
    MPI_Offset     offset_end;   /*   ending of aggregate access region */
    MPI_Offset     bufcount;     /* the number of buftype in this request */
    MPI_Datatype   buftype;      /* user defined derived data type */
    MPI_Datatype   ptype;        /* element data type in buftype */
    MPI_Datatype   imaptype;     /* derived data type constructed from imap */
    int           *status;       /* pointer to user's status */
} NC_req;

#define NC_ABUF_DEFAULT_TABLE_SIZE 128

typedef struct NC_buf_status {
    MPI_Aint   buf_addr;
    MPI_Offset req_size;
    int        is_used;
} NC_buf_status;

typedef struct NC_buf {
    MPI_Offset     size_allocated;
    MPI_Offset     size_used;
    int            table_size;
    int            tail;         /* index of last free entry */
    NC_buf_status *occupy_table; /* [table_size] */
    void          *buf;
} NC_buf;

/* chunk size for allocating read/write nonblocking request lists */
#define NC_REQUEST_CHUNK 1024

/* various file modes stored in flags */
#define NC_INDEP  0x10000   /* in independent data mode, cleared by endindep */
#define NC_CREAT  0x20000   /* in create phase, cleared by enddef */
#define NC_INDEF  0x80000   /* in define mode, cleared by enddef */
#define NC_NSYNC  0x100000  /* synchronise numrecs on change */
#define NC_HSYNC  0x200000  /* synchronise whole header on change */
#define NC_NDIRTY 0x400000  /* numrecs has changed */
#define NC_HDIRTY 0x800000  /* header info has changed */
struct NC {
    int           ncid;         /* file ID */
    int           flags;        /* various modes, i.e. define/data, fill,
                                   indep/coll, header dirty, etc */
    int           iomode;       /* cmode or omode used in ncmpi_create/open */
    int           mpiomode;     /* mode used in MPI_File_open, passed from
                                 * collective open to independent open */
    int           format;       /* 1, 2, or 5 corresponding to CDF-1, 2, or 5 */
    int           safe_mode;    /* 0 or 1, for parameter consistency check */
    int           numGetReqs;   /* number of pending nonblocking get requests */
    int           numPutReqs;   /* number of pending nonblocking put requests */
#ifdef ENABLE_SUBFILING
    int           subfile_mode; /* 0 or 1, for disable/enable subfiling */
    int           num_subfiles; /* number of subfiles */
    int           ncid_sf;      /* ncid of subfile */
#endif
    int           striping_unit; /* file stripe size of the file */
    MPI_Offset    chunk;       /* chunk size for reading header */
    MPI_Offset    h_align;     /* file alignment for header */
    MPI_Offset    v_align;     /* file alignment for each fixed variable */
    MPI_Offset    r_align;     /* file alignment for record variable section */
    MPI_Offset    h_minfree;   /* pad at the end of the header section */
    MPI_Offset    v_minfree;   /* pad at the end of the data section for fixed-size variables */
    MPI_Offset    xsz;       /* external size of this header, <= var[0].begin */
    MPI_Offset    begin_var; /* file offset of the first (non-record) var */
    MPI_Offset    begin_rec; /* file offset of the first 'record' */

    MPI_Offset    recsize;   /* length of 'record': sum of single record sizes
                                of all the record variables */
    MPI_Offset    numrecs;   /* number of 'records' allocated */
    MPI_Offset    put_size;  /* amount of writes committed so far in bytes */
    MPI_Offset    get_size;  /* amount of reads  committed so far in bytes */

    MPI_Comm      comm;           /* MPI communicator */
    MPI_Info      mpiinfo;        /* used MPI info object */
    MPI_File      collective_fh;  /* file handle for collective mode */
    MPI_File      independent_fh; /* file handle for independent mode */

    NC_dimarray   dims;     /* dimensions defined */
    NC_attrarray  attrs;    /* global attributes defined */
    NC_vararray   vars;     /* variables defined */

    NC_req       *get_list; /* list of nonblocking read requests */
    NC_req       *put_list; /* list of nonblocking write requests */
    NC_buf       *abuf;     /* attached buffer, used by bput APIs */

    char         *path;     /* file name */
    struct NC    *old;      /* contains the previous NC during redef. */
};

#define NC_readonly(ncp) \
        (!fIsSet((ncp)->iomode, NC_WRITE))

#define NC_IsNew(ncp) \
        fIsSet((ncp)->flags, NC_CREAT)

#define NC_indep(ncp) \
        fIsSet((ncp)->flags, NC_INDEP)

#define NC_indef(ncp) \
        (NC_IsNew(ncp) || fIsSet((ncp)->flags, NC_INDEF))

#define set_NC_ndirty(ncp) \
        fSet((ncp)->flags, NC_NDIRTY)

#define NC_ndirty(ncp) \
        fIsSet((ncp)->flags, NC_NDIRTY)

#define set_NC_hdirty(ncp) \
        fSet((ncp)->flags, NC_HDIRTY)

#define NC_hdirty(ncp) \
        fIsSet((ncp)->flags, NC_HDIRTY)

#define NC_dofill(ncp) \
        (!fIsSet((ncp)->flags, NC_NOFILL))

#define NC_doFsync(ncp) \
        fIsSet((ncp)->iomode, NC_SHARE)

#define NC_doHsync(ncp) \
        fIsSet((ncp)->flags, NC_HSYNC)

#define NC_doNsync(ncp) \
        fIsSet((ncp)->flags, NC_NSYNC)

#define NC_get_numrecs(ncp) \
        ((ncp)->numrecs)

#define NC_increase_numrecs(ncp, nrecs) \
        {if((nrecs) > (ncp)->numrecs) ((ncp)->numrecs = (nrecs));}

#define ErrIsHeaderDiff(err) \
        (NC_EMULTIDEFINE_FIRST >= (err) && (err) >= NC_EMULTIDEFINE_LAST)

#define IsPrimityMPIType(buftype) (buftype == MPI_FLOAT          || \
                                   buftype == MPI_DOUBLE         || \
                                   buftype == MPI_INT            || \
                                   buftype == MPI_CHAR           || \
                                   buftype == MPI_SIGNED_CHAR    || \
                                   buftype == MPI_UNSIGNED_CHAR  || \
                                   buftype == MPI_SHORT          || \
                                   buftype == MPI_UNSIGNED_SHORT || \
                                   buftype == MPI_UNSIGNED       || \
                                   buftype == MPI_LONG           || \
                                   buftype == MPI_LONG_LONG_INT  || \
                                   buftype == MPI_UNSIGNED_LONG_LONG)

/* Begin defined in nc.c */

extern int
ncmpio_cktype(int cdf_ver, nc_type datatype);

extern MPI_Offset
ncmpix_howmany(nc_type type, MPI_Offset xbufsize);

extern int
ncmpio_dset_has_recvars(NC *ncp);

extern int
ncmpio_write_header(NC *ncp);

extern void
ncmpio_free_NC(NC *ncp);

extern int
ncmpio_read_NC(NC *ncp);

/* End defined in nc.c */

#if 0
/* Begin defined in v1hpg.c */

extern size_t
ncx_len_NC(const NC *ncp, MPI_Offset sizeof_off_t);

extern int
ncx_put_NC(const NC *ncp, void **xpp, MPI_Offset offset, MPI_Offset extent);

extern int
nc_get_NC( NC *ncp);

/* End defined in v1hpg.c */

/* Begin defined in putget.c */

extern int
ncmpio_fill_NC_var(NC *ncp, const NC_var *varp, MPI_Offset recno);

extern int
ncmpio_inq_rec(int ncid, MPI_Offset *nrecvars, MPI_Offset *recvarids, MPI_Offset *recsizes);

extern int
ncmpio_get_rec(int ncid, MPI_Offset recnum, void **datap);

extern int
ncmpio_put_rec(int ncid, MPI_Offset recnum, void *const *datap);
#endif

/* End defined in putget.c */

/* Begin defined in header.c */
typedef struct bufferinfo {
    MPI_Comm    comm;
    MPI_File    collective_fh;
    MPI_Offset  get_size; /* amount of reads  committed so far in bytes */
    MPI_Offset  offset;   /* current read/write offset in the file */
    MPI_Offset  size;     /* size of the buffer */
    int         version;  /* 1, 2, and 5 for CDF-1, 2, and 5 respectively */
    int         safe_mode;/* 0: disabled, 1: enabled */
    void       *base;     /* beginning of read/write buffer */
    void       *pos;      /* current position in buffer */
} bufferinfo;

extern int
ncmpix_len_nctype(nc_type type);

#if 0
extern int
hdr_put_NC_attrarray(bufferinfo *pbp, const NC_attrarray *ncap);
#endif

extern MPI_Offset
ncmpio_hdr_len_NC(const NC *ncp);

extern int
ncmpio_hdr_get_NC(NC *ncp);

extern int
ncmpio_hdr_put_NC(NC *ncp, void *buf);

extern int
ncmpio_NC_computeshapes(NC *ncp);

extern int
ncmpio_hdr_check_NC(bufferinfo *getbuf, NC *ncp);
/* end defined in header.c */

/* begin defined in mpincio.c */
extern int
ncmpio_file_sync(NC *ncp);

extern int
ncmpiio_get_hint(NC *ncp, char *key, char *value, int *flag);

extern int
NC_computeshapes(NC *ncp);

/* end defined in mpincio.h */

/* begin defined in error.c */

int ncmpio_handle_error(int mpi_errorcode, char *msg);

/* end defined in error.c */
/*
 * These functions are used to support
 * interface version 2 backward compatibility.
 * N.B. these are tested in ../nc_test even though they are
 * not public. So, be careful to change the declarations in
 * ../nc_test/tests.h if you change these.
 */

int ncmpio_x_putn_NC_CHAR  (void *xbuf, const void *buf, MPI_Offset nelems,
                           MPI_Datatype datatype);
int ncmpio_x_putn_NC_BYTE  (int cdf_ver,
                           void *xbuf, const void *buf, MPI_Offset nelems,
                           MPI_Datatype datatype, void *fillp);
int ncmpio_x_putn_NC_UBYTE (void *xbuf, const void *buf, MPI_Offset nelems,
                           MPI_Datatype datatype, void *fillp);
int ncmpio_x_putn_NC_SHORT (void *xbuf, const void *buf, MPI_Offset nelems,
                           MPI_Datatype datatype, void *fillp);
int ncmpio_x_putn_NC_USHORT(void *xbuf, const void *buf, MPI_Offset nelems,
                           MPI_Datatype datatype, void *fillp);
int ncmpio_x_putn_NC_INT   (void *xbuf, const void *buf, MPI_Offset nelems,
                           MPI_Datatype datatype, void *fillp);
int ncmpio_x_putn_NC_UINT  (void *xbuf, const void *buf, MPI_Offset nelems,
                           MPI_Datatype datatype, void *fillp);
int ncmpio_x_putn_NC_FLOAT (void *xbuf, const void *buf, MPI_Offset nelems,
                           MPI_Datatype datatype, void *fillp);
int ncmpio_x_putn_NC_DOUBLE(void *xbuf, const void *buf, MPI_Offset nelems,
                           MPI_Datatype datatype, void *fillp);
int ncmpio_x_putn_NC_INT64 (void *xbuf, const void *buf, MPI_Offset nelems,
                           MPI_Datatype datatype, void *fillp);
int ncmpio_x_putn_NC_UINT64(void *xbuf, const void *buf, MPI_Offset nelems,
                           MPI_Datatype datatype, void *fillp);

int ncmpio_x_getn_NC_CHAR  (const void *xbuf, void *buf, MPI_Offset nelems,
                           MPI_Datatype datatype);
int ncmpio_x_getn_NC_BYTE  (int cdf_ver,
                           const void *xbuf, void *buf, MPI_Offset nelems,
                           MPI_Datatype datatype);
int ncmpio_x_getn_NC_UBYTE (const void *xbuf, void *buf, MPI_Offset nelems,
                           MPI_Datatype datatype);
int ncmpio_x_getn_NC_SHORT (const void *xbuf, void *buf, MPI_Offset nelems,
                           MPI_Datatype datatype);
int ncmpio_x_getn_NC_USHORT(const void *xbuf, void *buf, MPI_Offset nelems,
                           MPI_Datatype datatype);
int ncmpio_x_getn_NC_INT   (const void *xbuf, void *buf, MPI_Offset nelems,
                           MPI_Datatype datatype);
int ncmpio_x_getn_NC_UINT  (const void *xbuf, void *buf, MPI_Offset nelems,
                           MPI_Datatype datatype);
int ncmpio_x_getn_NC_FLOAT (const void *xbuf, void *buf, MPI_Offset nelems,
                           MPI_Datatype datatype);
int ncmpio_x_getn_NC_DOUBLE(const void *xbuf, void *buf, MPI_Offset nelems,
                           MPI_Datatype datatype);
int ncmpio_x_getn_NC_INT64 (const void *xbuf, void *buf, MPI_Offset nelems,
                           MPI_Datatype datatype);
int ncmpio_x_getn_NC_UINT64(const void *xbuf, void *buf, MPI_Offset nelems,
                           MPI_Datatype datatype);

int NC_start_count_stride_ck(const NC *ncp, const NC_var *varp,
                const MPI_Offset *start, const MPI_Offset *count,
                const MPI_Offset *stride, const int rw_flag);

int ncmpio_need_convert(int format, nc_type nctype, MPI_Datatype mpitype);

int ncmpio_need_swap(nc_type nctype,MPI_Datatype mpitype);

void ncmpio_swapn(void *dest_p, const void* src_p, MPI_Offset nelems, int esize);

void ncmpio_in_swapn(void *buf, MPI_Offset nelems, int esize);

int ncmpio_is_request_contiguous(NC *ncp, NC_var *varp,
                const MPI_Offset starts[], const MPI_Offset  counts[]);

int ncmpio_get_offset(NC *ncp, NC_var *varp, const MPI_Offset starts[],
                const MPI_Offset counts[], const MPI_Offset strides[],
                const int rw_flag, MPI_Offset *offset_ptr);

int ncmpio_check_mpifh(NC* ncp, int collective);

int ncmpio_write_numrecs(NC *ncp, MPI_Offset new_numrecs);

int ncmpiio_sync_numrecs(NC *ncp);

int ncmpio_vars_create_filetype(NC* ncp, NC_var* varp,
                const MPI_Offset start[], const MPI_Offset count[],
                const MPI_Offset stride[], int rw_flag, int *blocklen,
                MPI_Offset *offset, MPI_Datatype *filetype,
                int *is_filetype_contig);

extern int
ncmpio_igetput_varm(NC *ncp, NC_var *varp, const MPI_Offset *start,
                const MPI_Offset *stride, const MPI_Offset *imap,
                const MPI_Offset *count, void *buf, MPI_Offset bufcount,
                MPI_Datatype datatype, int *reqid, int rw_flag, int use_abuf,
                int isSameGroup);

extern int
ncmpio_inq_files_opened(int *num, int *ncids);

extern MPI_Datatype
ncmpio_nc2mpitype(nc_type type);

extern nc_type
ncmpio_mpi2nctype(MPI_Datatype itype);

extern int                
ncmpio_file_set_view(NC *ncp, MPI_File fh, MPI_Offset *offset, MPI_Datatype filetype);

extern int
ncmpio_create_imaptype(NC_var *varp, const MPI_Offset *count,
                       const MPI_Offset *imap, const MPI_Offset  bnelems,
                       const int el_size, MPI_Datatype ptype,
                       MPI_Datatype *imaptype);

extern int
ncmpio_calc_datatype_elems(NC_var *varp, const MPI_Offset *count,
                           MPI_Datatype buftype,
                           MPI_Datatype *ptype, MPI_Offset *bufcount,
                           MPI_Offset *bnelems, MPI_Offset *nbytes,
                           int *el_size, int *buftype_is_contig);

extern int
ncmpio_fill_vars(NC *ncp);

extern int
ncmpio_sanity_check(NC *ncp, int varid, const MPI_Offset *start,
                    const MPI_Offset *count, const MPI_Offset *stride,
                    MPI_Offset bufcount, MPI_Datatype buftype,
                    api_kind api, int isFlexibleAPI, int mustInDataMode,
                    int rw_flag, int io_method, NC_var **varp);

extern char*
ncmpio_err_code_name(int err);

extern int
ncmpio_jenkins_one_at_a_time_hash(const char *str_name);

extern int
ncmpio_additive_hash(const char *str_name);

extern int
ncmpio_rotating_hash(const char *str_name);

extern int
ncmpio_Bernstein_hash(const char *str_name);

extern int
ncmpio_Pearson_hash(const char *str_name);

extern int
ncmpio_update_name_lookup_table(NC_nametable *nameT, const int id,
                                const char *oldname, const char *newname);

extern int
ncmpio_inq_var_fill(NC_var *varp, void *fill_value);

extern int
ncmpio_inq_default_fill_value(int type, void *fill_value);

extern int
ncmpio_getput_zero_req(NC *ncp, int rw_flag);

extern int
ncmpio_NC_check_vlens(NC *ncp);

extern void
ncmpio_set_pnetcdf_hints(NC *ncp, MPI_Info info);

extern int
ncmpio_close_files(NC *ncp, int doUnlink);

#endif /* _NC_H */
