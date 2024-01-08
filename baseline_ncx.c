#include <stdio.h>
#include <stdlib.h>
#include <string.h> /* strcpy(), strncpy() */
#include <unistd.h> /* getopt() */
#include <time.h>   /* time() localtime(), asctime() */
#include <assert.h>
#include <mpi.h>
#include <pnetcdf.h>
#include "baseline_ncx.h" 


int
xlen_nc_type(nc_type xtype, int *size)
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

static int 
putn_text(void **xpp, MPI_Offset nelems, const char *tp)
{
	(void) memcpy(*xpp, tp, (size_t)nelems);
	*xpp = (void *)((char *)(*xpp) + nelems);
	return 0;
}

static int
put_uint32(void **xpp, unsigned int ip)
{
    memcpy(*xpp, &ip, 4);
    /* advace *xpp 4 bytes */
    *xpp  = (void *)((char *)(*xpp) + 4);
    return 0;
}

static int
serialize_dim(bufferinfo   *pbp,
               const hdr_dim *dimp)
{
    /* copy name */
    serialize_name(pbp, dimp->name);
    put_uint32((void**)(&pbp->pos), (uint32_t)dimp->size);
    return 0;
}

static int
serialize_name(bufferinfo *pbp,
                const char *name)
{
    size_t nchars = strlen(name);

    put_uint32((void**)(&pbp->pos), (uint32_t)nchars);

    return putn_text((void **)(&pbp->pos), (MPI_Offset)nchars, name);
}


static int
serialize_dimarray(bufferinfo        *pbp,
                    const hdr_dimarray *ncap)
{
    int i, status;
    assert(pbp != NULL);

    if (ncap == NULL || ncap->ndefined == 0) { /* ABSENT */
        status = put_uint32((void**)(&pbp->pos), NC_UNSPECIFIED);
        if (status != NC_NOERR) return status;
        status = put_uint32((void**)(&pbp->pos), 0);
        if (status != NC_NOERR) return status;
    }
    else {
        /* copy NC_DIMENSION */
        status = put_uint32((void**)(&pbp->pos), NC_DIMENSION);
        if (status != NC_NOERR) return status;
        status = put_uint32((void**)(&pbp->pos), (uint32_t)ncap->ndefined);
        if (status != NC_NOERR) return status;
        /* copy [dim ...] */
        for (i=0; i<ncap->ndefined; i++) {
            status = serialize_dim(pbp, ncap->value[i]);
            if (status != NC_NOERR) return status;
        }
    }
    return 0;
}

static int
serialize_attrV(bufferinfo    *pbp,
                 const hdr_attr *attrp)
{

    int xsz;
    MPI_Offset sz;

    /* xlen_nc_type() returns the element size (unaligned) of
     * attrp->xtype attrp->xsz is the aligned total size of attribute values
     */
    xlen_nc_type(attrp->xtype, &xsz);
    sz = attrp->nelems * xsz;
    memcpy(pbp->pos, attrp->xvalue, (size_t)sz);
    pbp->pos = (void *)((char *)pbp->pos + sz);
    return 0;
}

/*----< serialize_NC_attr() >--------------------------------------------------*/
static int
serialize_attr(bufferinfo    *pbp,
                const hdr_attr *attrp)
{
    int status;

    /* copy name */
    status = serialize_name(pbp, attrp->name);
    if (status != NC_NOERR) return status;

    /* copy nc_type */
    status = put_uint32((void**)(&pbp->pos), (uint32_t)attrp->xtype);
    if (status != NC_NOERR) return status;

    /* copy nelems */
    status = put_uint32((void**)(&pbp->pos), (uint32_t)attrp->nelems);
    if (status != NC_NOERR) return status;

    /* copy [values ...] */
    status = serialize_attrV(pbp, attrp);
    if (status != NC_NOERR) return status;

    return NC_NOERR;
}

static int
serialize_attrarray(bufferinfo         *pbp,
                     const hdr_attrarray *ncap)
{

    int i, status;
    assert(pbp != NULL);
    if (ncap == NULL || ncap->ndefined == 0) { /* ABSENT */
        status = put_uint32((void**)(&pbp->pos), NC_UNSPECIFIED);
        if (status != NC_NOERR) return status;
        status = put_uint32((void**)(&pbp->pos), 0);
        if (status != NC_NOERR) return status;
    }
    else {
        /* copy NC_ATTRIBUTE */
        status = put_uint32((void**)(&pbp->pos), NC_ATTRIBUTE);
        if (status != NC_NOERR) return status;
        status = put_uint32((void**)(&pbp->pos), (uint32_t)ncap->ndefined);
        if (status != NC_NOERR) return status;
        /* copy [attr ...] */
        for (i=0; i<ncap->ndefined; i++) {
            status = serialize_attr(pbp, ncap->value[i]);
            if (status != NC_NOERR) return status;
        }
    }

    return NC_NOERR;
}

static int
serialize_var(bufferinfo   *pbp,
               const hdr_var *varp)
{

    int i, status;

    /* copy name */
    status = serialize_name(pbp, varp->name);
    if (status != NC_NOERR) return status;

    /* copy nelems */

    status = put_uint32((void**)(&pbp->pos), (uint32_t)varp->ndims);

    if (status != NC_NOERR) return status;

    /* copy [dim_index ...] i*/
    for (i=0; i<varp->ndims; i++) {

        status = put_uint32((void**)(&pbp->pos), (uint32_t)varp->dimids[i]);
        if (status != NC_NOERR) return status;
    }

    /* copy vatt_list */
    status = serialize_attrarray(pbp, &varp->attrs);
    if (status != NC_NOERR) return status;

    /* copy nc_type */
    status = put_uint32((void**)(&pbp->pos), (uint32_t)varp->xtype);
    if (status != NC_NOERR) return status;

    return NC_NOERR;
}


/*----< serialize_vararray() >----------------------------------------------*/
static int
serialize_vararray(bufferinfo        *pbp,
                    const hdr_vararray *ncap)
{
    int i, status;
    assert(pbp != NULL);
    if (ncap == NULL || ncap->ndefined == 0) { /* ABSENT */
        status = put_uint32((void**)(&pbp->pos), NC_UNSPECIFIED);
        if (status != NC_NOERR) return status;
        status = put_uint32((void**)(&pbp->pos), 0);
        if (status != NC_NOERR) return status;
    }
    else {
        /* copy NC_VARIABLE */
        status = put_uint32((void**)(&pbp->pos), NC_VARIABLE);
        if (status != NC_NOERR) return status;

        /* copy nelems */
        status = put_uint32((void**)(&pbp->pos), (uint32_t)ncap->ndefined);
        if (status != NC_NOERR) return status;

        /* copy [var ...] */
        for (i=0; i<ncap->ndefined; i++) {
            status =serialize_var(pbp, ncap->value[i]);
            if (status != NC_NOERR) return status;
        }
    }
    return NC_NOERR;
}



/*----< serialize_hdr() >----------------------------------------------*/
int
serialize_hdr(struct hdr *ncp, void *buf)
{
    int status;
    bufferinfo putbuf;

    putbuf.pos           = buf;
    putbuf.base          = buf;
    putbuf.size          = ncp->xsz;

    /* copy dim_list */
    status = serialize_dimarray(&putbuf, &ncp->dims);
    if (status != NC_NOERR) return status;

    /* copy gatt_list */
    status = serialize_attrarray(&putbuf, &ncp->attrs);
    if (status != NC_NOERR) return status;

    /* copy var_list */
    status = serialize_vararray(&putbuf, &ncp->vars);
    if (status != NC_NOERR) return status;

    return NC_NOERR;
}

