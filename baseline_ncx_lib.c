#include <stdio.h>
#include <stdlib.h>
#include <string.h> /* strcpy(), strncpy() */
#include <unistd.h> /* getopt() */
#include <time.h>   /* time() localtime(), asctime() */
#include <assert.h>
#include <mpi.h>
#include <pnetcdf.h>
#include "baseline_ncx_lib.h" 
#include "mem_tracker.h"

#ifdef MEM_TRACKING
#define malloc(size) tracked_malloc(size)
#define free(ptr)    tracked_free(ptr)
#endif

  /* ---------------------------------- Serializaition ----------------------------------------*/

int
meta_xlen_nc_type(nc_type xtype, int *size)
{
    switch(xtype) {
        case NC_BYTE:
        case NC_CHAR:
        case NC_UBYTE:  *size = 1; return 0;
        case NC_SHORT:
        case NC_USHORT: *size = 2; return 0;
        case NC_INT:
        case NC_UINT:
        case NC_FLOAT:  *size = 4; return 0;
        case NC_DOUBLE:
        case NC_INT64:
        case NC_UINT64: *size = 8; return 0;
    }
    return 0;
}


void meta_free_hdr_dim(hdr_dim *dim) {
    if (dim != NULL) {
        free(dim->name);
        free(dim);
    }
}

void meta_free_hdr_dimarray(hdr_dimarray *dims) {
    if (dims != NULL) {
        for (int i = 0; i < dims->ndefined; i++) {
            meta_free_hdr_dim(dims->value[i]);
        }
        free(dims->value);
        //free(dims);
    }
}

void meta_free_hdr_attr(hdr_attr *attr) {
    if (attr != NULL) {
        free(attr->name);
        free(attr->xvalue);
        // free(attr);
    }
}

void meta_free_hdr_attrarray(hdr_attrarray *attrs) {
    if (attrs != NULL) {
        if (attrs->value != NULL) {
            for (int i = 0; i < attrs->ndefined; i++) {
                meta_free_hdr_attr(attrs->value[i]);
            }
            free(attrs->value);
            attrs->value = NULL;
            // free(attrs);
        }
    }
}

void meta_free_hdr_var(hdr_var *var) {
    if (var != NULL) {
        free(var->name);
        free(var->dimids);
        meta_free_hdr_attrarray(&(var->attrs));
        free(var);
    }
}

void meta_free_hdr_vararray(hdr_vararray *vars) {
    if (vars != NULL) {
        for (int i = 0; i < vars->ndefined; i++) {
            meta_free_hdr_var(vars->value[i]);
        }
        free(vars->value);
        // free(vars);
    }
}

void meta_free_hdr(struct hdr *header) {
    if (header != NULL) {
        meta_free_hdr_dimarray(&(header->dims));
        // free_hdr_attrarray_meta(&(header->attrs));
        meta_free_hdr_vararray(&(header->vars));
        // free(header);
    }
}

