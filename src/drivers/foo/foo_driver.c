/*
 *  Copyright (C) 2017, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* $Id$ */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <dispatch.h>
#include <foo_driver.h>

static PNC_driver foo_driver = {
    /* FILE APIs */
    foo_create,
    foo_open,
    foo_close,
    foo_enddef,
    foo__enddef,
    foo_redef,
    foo_sync,
    foo_abort,
    foo_set_fill,
    foo_inq,
    foo_inq_misc,
    foo_sync_numrecs,
    foo_begin_indep_data,
    foo_end_indep_data,

    /* DIMENSION APIs */
    foo_def_dim,
    foo_inq_dimid,
    foo_inq_dim,
    foo_rename_dim,

    /* ATTRIBUTE APIs */
    foo_inq_att,
    foo_inq_attid,
    foo_inq_attname,
    foo_copy_att,
    foo_rename_att,
    foo_del_att,
    foo_get_att,
    foo_put_att,

    /* VARIABLE APIs */
    foo_def_var,
    foo_def_var_fill,
    foo_fill_var_rec,
    foo_inq_var,
    foo_inq_varid,
    foo_rename_var,
    foo_get_var,
    foo_put_var,
    foo_get_varn,
    foo_put_varn,
    foo_get_vard,
    foo_put_vard,
    foo_iget_var,
    foo_iput_var,
    foo_bput_var,
    foo_iget_varn,
    foo_iput_varn,
    foo_bput_varn,

    foo_buffer_attach,
    foo_buffer_detach,
    foo_wait,
    foo_cancel
};

PNC_driver* foo_inq_driver(void) {
    return &foo_driver;
}

