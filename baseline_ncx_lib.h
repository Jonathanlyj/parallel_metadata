#include <stddef.h> 
#include <mpi.h>
#include <pnetcdf.h>


#define NC_NOERR        0  


typedef enum {
    NC_UNSPECIFIED =  0,  /* ABSENT */
    NC_DIMENSION   = 10,  /* \x00 \x00 \x00 \x0A */
    NC_VARIABLE    = 11,  /* \x00 \x00 \x00 \x0B */
    NC_ATTRIBUTE   = 12   /* \x00 \x00 \x00 \x0C */
} NC_tag;

typedef struct {
    MPI_Offset  size;
    size_t      name_len; /* strlen(name), for faster string compare */
    char       *name;
} hdr_dim;


typedef struct hdr_dimarray {
    int            ndefined;      /* number of defined dimensions */
    // int            unlimited_id;  /* -1 for not defined, otherwise >= 0 */
    hdr_dim       **value;
} hdr_dimarray;

typedef struct {
    MPI_Offset nelems;   /* number of attribute elements */
    nc_type    xtype;    /* external NC data type of the attribute */
    size_t     name_len; /* strlen(name) for faster string compare */
    char      *name;     /* name of the attributes */
    void      *xvalue;   /* the actual data, in external representation */
} hdr_attr;


typedef struct hdr_attrarray {
    int            ndefined;  /* number of defined attributes */
    hdr_attr      **value;
} hdr_attrarray;


typedef struct {
    nc_type       xtype;   /* variable's external NC data type */
    size_t        name_len;/* strlen(name) for faster string compare */
    char         *name;    /* name of the variable */
    int           ndims;   /* number of dimensions */
    int          *dimids;  /* [ndims] array of dimension IDs */
    hdr_attrarray  attrs;   /* attribute array */
} hdr_var;


typedef struct hdr_vararray {
    int            ndefined;    /* number of defined variables */
    hdr_var       **value;
} hdr_vararray;


/* various file modes stored in flags */
struct hdr {
    MPI_Offset    xsz;      /* size occupied on the buffer */
    hdr_dimarray   dims;     /* dimensions defined */
    hdr_attrarray  attrs;    /* global attributes defined */
    hdr_vararray   vars;     /* variables defined */
};



typedef struct metabuffer {
    // MPI_Comm    comm;
    int         size;     /* allocated size of the buffer */
    char       *base;     /* beginning of read/write buffer */
    char       *pos;      /* current position in buffer */
    char       *end;      /* end position of buffer */
} metabuffer;


// /* Function prototypes */
int meta_xlen_nc_type(nc_type xtype, int *size);
void meta_free_hdr_dim(hdr_dim *dim);
void meta_free_hdr_dimarray(hdr_dimarray *dims);
void meta_free_hdr_attr(hdr_attr *attr);
void meta_free_hdr_attrarray(hdr_attrarray *attrs);
void meta_free_hdr_var(hdr_var *var);
void meta_free_hdr_vararray(hdr_vararray *vars);
void meta_free_hdr(struct hdr *header);
