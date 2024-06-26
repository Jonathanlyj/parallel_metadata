/*********************************************************************
 *
 *  Copyright (C) 2023, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 *
 *********************************************************************/
/* $Id$ */


#include <stdio.h>
#include <stdlib.h>
#include <string.h> /* strcpy(), strncpy() */
#include <unistd.h> /* getopt() */
#include <time.h>   /* time() localtime(), asctime() */
#include <mpi.h>
#include <pnetcdf.h>
#include "baseline_ncx_app.h" 



static int verbose;


#define ERR {if(err!=NC_NOERR){printf("Error at %s:%d : %s\n", __FILE__,__LINE__, ncmpi_strerror(err));nerrs++;}}

#define FILE_NAME "/files2/scratch/yll6162/parallel_metadata/nue_slice_panoptic_hdf_merged.nc"

#define OUTPUT_NAME "out.nc"
// #define FILE_NAME "testfile.nc"


/* ---------------------------------- Read Metadata ----------------------------------------*/


void read_metdata(int rank, int size, struct hdr *file_info) {

    int ncid, num_vars, num_dims, tot_num_dims, elem_sz, v_attrV_xsz, status;


    file_info->dims.ndefined = 0;
    file_info->attrs.ndefined = 0;
    file_info->vars.ndefined = 0;

    tot_num_dims = 0;
    // Open the NetCDF file
    if (ncmpi_open(MPI_COMM_WORLD, FILE_NAME, NC_NOWRITE, MPI_INFO_NULL, &ncid) != NC_NOERR) {
        fprintf(stderr, "Error opening NetCDF file.\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    // Get the number of variables
    if (ncmpi_inq_nvars(ncid, &num_vars) != NC_NOERR) {
        fprintf(stderr, "Error getting the number of variables.\n");
        ncmpi_close(ncid);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }


    file_info->vars.ndefined = num_vars;
    file_info->vars.value = (hdr_var **)malloc(file_info->vars.ndefined * sizeof(hdr_var *));
    // Each process reads its all variables
    for (int i = 0; i <  num_vars; ++i) {
        // Get variable information
        hdr_var *variable_info = (hdr_var *)malloc(sizeof(hdr_var));
        variable_info->ndims = 0;  // Initialize the number of dimensions
        variable_info->attrs.ndefined = 0;  // Initialize the number of attributes

        // Get variable information
        char var_name[NC_MAX_NAME + 1];
        ncmpi_inq_varname(ncid, i, var_name);
        variable_info->name_len = strlen(var_name);
        variable_info->name = (char *)malloc((variable_info->name_len + 1) * sizeof(char));
        strcpy(variable_info->name, var_name);

        // Get number of dimensions
        ncmpi_inq_varndims(ncid, i, &(variable_info->ndims));

        num_dims = variable_info->ndims;

        // Allocate memory for dimension IDs
        variable_info->dimids = (int *)malloc(variable_info->ndims * sizeof(int));

        // Get dimension IDs
        ncmpi_inq_vardimid(ncid, i, variable_info->dimids);

        // Read attributes
        nc_type xtype;
        int num_attrs;
        ncmpi_inq_var(ncid, i, NULL, &xtype, NULL, NULL, &num_attrs);
        variable_info->attrs.ndefined = num_attrs;
        variable_info->xtype = xtype;



        //Read and store dimension information 
        file_info->dims.ndefined += num_dims;
        if (tot_num_dims == 0) {
            file_info->dims.value = (hdr_dim **)malloc(file_info->dims.ndefined * sizeof(hdr_dim *));
        }else{
            file_info->dims.value = (hdr_dim **)realloc(file_info->dims.value, file_info->dims.ndefined * sizeof(hdr_dim *));
        }
       
        for (int k = 0; k < num_dims; ++k) {
            hdr_dim *dimension_info = (hdr_dim *)malloc(sizeof(hdr_dim));
            int dimid = variable_info->dimids[k];
            // Get dimension name
            char dim_name[NC_MAX_NAME + 1];
            ncmpi_inq_dimname(ncid, dimid, dim_name);
            dimension_info->name_len = strlen(dim_name);
            dimension_info->name = (char *)malloc((dimension_info->name_len + 1) * sizeof(char));
            strcpy(dimension_info->name, dim_name);

            // Get dimension size
            ncmpi_inq_dimlen(ncid, dimid, &(dimension_info->size));

            file_info->dims.value[k + tot_num_dims] = dimension_info;
            variable_info->dimids[k] = k + tot_num_dims; //overwriting previous global dim id to local dim id
            file_info->xsz += sizeof(uint32_t) + sizeof(char) * dimension_info -> name_len; // dim name
            file_info->xsz += sizeof(uint32_t); //size

        }
        tot_num_dims += num_dims;

        // Add the variable information to the file_info structure

        file_info->vars.value[i] = variable_info;
    }

    // Close the NetCDF file
    ncmpi_close(ncid);
}

/* ---------------------------------- Decode Metadata ----------------------------------------*/
int define_hdr(struct hdr *hdr_data, int ncid, int rank){
    //define dimensions
    int ndims= hdr_data->dims.ndefined;
    int *dimid = (int *)malloc(ndims * sizeof(int));
    int i,j,k,nerrs=0;
    int err;

    for (i=0; i<ndims; i++){
        // err = ncmpi_inq_dimid(ncid, hdr_data->dims.value[i]->name, &dimid[i]);
        // if (err != NC_NOERR) {
        // if (rank == 0) printf("\n %s", hdr_data->dims.value[i]->name);
        err = ncmpi_def_dim(ncid, hdr_data->dims.value[i]->name,  hdr_data->dims.value[i]->size, &dimid[i]); ERR
        // }
    }

    //define variables
    // int nvars = hdr_data->vars.ndefined;
    // int *varid = (int *)malloc(nvars * sizeof(int));
    // int v_ndims, v_namelen, xtype, n_att;
    // int *v_dimids;
    // int att_namelen, att_xtype, att_nelems;

    // for (i=0; i<nvars; i++){
        
    //     v_namelen =  hdr_data->vars.value[i]->name_len;
    //     xtype = hdr_data->vars.value[i]->xtype;

    //     v_ndims = hdr_data->vars.value[i]->ndims;
    //     v_dimids = (int *)malloc(v_ndims * sizeof(int));
    //     for(j=0; j<v_ndims; j++) v_dimids[j] = dimid[hdr_data->vars.value[i]->dimids[j]];
    //     err = ncmpi_def_var(ncid, hdr_data->vars.value[i]->name, xtype, v_ndims,  v_dimids, &varid[i]); ERR
    //     n_att = hdr_data->vars.value[i]->attrs.ndefined;
    //     // printf("\nn_att: %d\n", n_att);
    //     for(k=0; k<n_att; k++){
    //         att_namelen = hdr_data->vars.value[i]->attrs.value[k]->name_len;
    //         att_xtype = hdr_data->vars.value[i]->attrs.value[k]->xtype;
    //         att_nelems = hdr_data->vars.value[i]->attrs.value[k]->nelems;
    //         err = ncmpi_put_att(ncid, varid[i],  hdr_data->vars.value[i]->attrs.value[k]->name, att_xtype, 
    //         att_nelems, &hdr_data->vars.value[i]->attrs.value[k]->xvalue[0]); ERR
    //     }
        
    // }
    return nerrs;
}



int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);
    int rank, size, status, err, nerrs=0;
    double start_time1, end_time1, start_time2, end_time2;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    struct hdr all_hdr;

    read_metdata(rank, size, &all_hdr);




    int ncid, cmode;
    char filename[256];
    cmode = NC_64BIT_DATA | NC_CLOBBER;


    err = ncmpi_create(MPI_COMM_WORLD, OUTPUT_NAME, cmode, MPI_INFO_NULL, &ncid); ERR
    define_hdr(&all_hdr, ncid, rank);
    free_hdr(&all_hdr);
    err = ncmpi_enddef(ncid); ERR


    err = ncmpi_close(ncid); ERR


    MPI_Finalize();
    return 0;



}