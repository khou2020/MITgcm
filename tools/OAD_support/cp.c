

#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <linux/limits.h>
#include <fcntl.h>
#ifdef ALLOW_USE_MPI
#include <mpi.h>
#endif
#include "zfp.h"

static int cp_file_num = 0;
//int fd;
FILE *fh;

void cp_wr_open_(int *num){
    int rank;
    char fname[PATH_MAX];

    if (*num > 0){
        cp_file_num = *num;
    }

#ifdef ALLOW_USE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#else
    rank = 0;
#endif

    sprintf(fname, "oad_cp.%03d.%05d", rank, cp_file_num);
    fh = fopen (fname, "wb");
    //fd = open(fname, O_CREAT | O_WRONLY, 0644);

    if (num == NULL){
        cp_file_num++;
    }
}

void cp_rd_open_(int *num){
    int rank;
    char fname[PATH_MAX];

    if (*num > 0){
        cp_file_num = *num;
    }
    else{
        cp_file_num--;
    }

#ifdef ALLOW_USE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#else
    rank = 0;
#endif

    sprintf(fname, "oad_cp.%03d.%05d", rank, cp_file_num);

    fh = fopen (fname, "rb");
    //fd = open(fname, O_CREAT | O_RDONLY, 0644);
}

void cpc_close_(){
    //close(fd);
    fclose(fh);
}

void compresswr_double(void *data, size_t size, int dim, int *shape){
    zfp_stream *zfp;
    zfp_type type = zfp_type_double;                          
    zfp_field *field;
    bitstream *stream;
    size_t bufsize;  
    size_t zfpsize;
    void *buffer;

    // array metadata
    if (dim == 1){
        field = zfp_field_1d(data, type, (unsigned int)(shape[0]) );
    }
    else if (dim == 2){
        field = zfp_field_2d(data, type, (unsigned int)(shape[0]), (unsigned int)(shape[1]));
    }
    else{
        int i;
        unsigned int totalsize = 1;

        for(i = 2; i < dim; i++){
            totalsize *= shape[i];
        }
        field = zfp_field_3d(data, type, (unsigned int)(shape[0]), (unsigned int)(shape[1]), totalsize);
    }

    // compressed stream and parameters
    zfp = zfp_stream_open(NULL);   

    // set tolerance for fixed-accuracy mode           
    zfp_stream_set_accuracy(zfp, 0.1);                  
    //precision = zfp_stream_set_precision(zfp, 16);           
    //zfp_stream_set_rate(zfp, rate, type, 3, 0);       

    // allocate buffer for compressed data
    bufsize = zfp_stream_maximum_size(zfp, field);    
    buffer = malloc(bufsize);                        

    // associate bit stream with allocated buffer
    stream = stream_open(buffer, bufsize);      
    zfp_stream_set_bit_stream(zfp, stream);                  
    zfp_stream_rewind(zfp);                  

    // compress array
    zfpsize = zfp_compress(zfp, field);               

    fwrite(&zfpsize, sizeof(zfpsize), 1, fh);
    fwrite(buffer, zfpsize, 1, fh);
    //write(fd, &zfpsize, sizeof(zfpsize));
    //write(fd, buffer, zfpsize);

    printf("Write %zu -> %zu\n", size, zfpsize);
}

void compressrd_double(void *data, size_t size, int dim, int *shape){
    zfp_stream *zfp;
    int ret;
    zfp_type type = zfp_type_double;                          
    zfp_field *field;
    bitstream *stream;
    size_t bufsize;  
    size_t zfpsize;
    void *buffer;

    // array metadata
    if (dim == 1){
        field = zfp_field_1d(data, type, (unsigned int)(shape[0]) );
    }
    else if (dim == 2){
        field = zfp_field_2d(data, type, (unsigned int)(shape[0]), (unsigned int)(shape[1]));
    }
    else{
        int i;
        unsigned int totalsize = 1;

        for(i = 2; i < dim; i++){
            totalsize *= shape[i];
        }
        field = zfp_field_3d(data, type, (unsigned int)(shape[0]), (unsigned int)(shape[1]), totalsize);
    }

    // compressed stream and parameters
    zfp = zfp_stream_open(NULL);   

    // set tolerance for fixed-accuracy mode           
    zfp_stream_set_accuracy(zfp, 0.1);                  
    //precision = zfp_stream_set_precision(zfp, 16);           
    //zfp_stream_set_rate(zfp, rate, type, 3, 0);       

    // allocate buffer for compressed data                     
    fread(&bufsize, sizeof(bufsize), 1, fh);
    //read(fd, &bufsize, sizeof(bufsize));
    buffer = malloc(bufsize);   
    fread(buffer, bufsize, 1, fh);
    //read(fd, buffer, bufsize);

    // associate bit stream with allocated buffer
    stream = stream_open(buffer, bufsize);      
    zfp_stream_set_bit_stream(zfp, stream);                  
    zfp_stream_rewind(zfp);                  

    ret = zfp_decompress(zfp, field);

    if (ret > 0) {
        printf("Read %zu -> %zu\n", bufsize, size);
    }
    else {
        printf("Decompress fail: addr: %llx, size: %zu\n", data, size);
    }
}

void compresswr_int(void *data, size_t size, int dim, int *shape){
    zfp_stream *zfp;
    zfp_type type = zfp_type_int32;                          
    zfp_field *field;
    bitstream *stream;
    size_t bufsize;  
    size_t zfpsize;
    void *buffer;

    // array metadata
    if (dim == 1){
        field = zfp_field_1d(data, type, (unsigned int)(shape[0]) );
    }
    else if (dim == 2){
        field = zfp_field_2d(data, type, (unsigned int)(shape[0]), (unsigned int)(shape[1]));
    }
    else{
        int i;
        unsigned int totalsize = 1;

        for(i = 2; i < dim; i++){
            totalsize *= shape[i];
        }
        field = zfp_field_3d(data, type, (unsigned int)(shape[0]), (unsigned int)(shape[1]), totalsize);
    }

    // compressed stream and parameters
    zfp = zfp_stream_open(NULL);   

    // set tolerance for fixed-accuracy mode           
    zfp_stream_set_accuracy(zfp, 0);                  
    //precision = zfp_stream_set_precision(zfp, 16);           
    //zfp_stream_set_rate(zfp, rate, type, 3, 0);       

    // allocate buffer for compressed data
    bufsize = zfp_stream_maximum_size(zfp, field);    
    buffer = malloc(bufsize);                        

    // associate bit stream with allocated buffer
    stream = stream_open(buffer, bufsize);      
    zfp_stream_set_bit_stream(zfp, stream);                  
    zfp_stream_rewind(zfp);                  

    // compress array
    zfpsize = zfp_compress(zfp, field);               

    fwrite(&zfpsize, sizeof(zfpsize), 1, fh);
    fwrite(buffer, zfpsize, 1, fh);
    //write(fd, &zfpsize, sizeof(zfpsize));
    //write(fd, buffer, zfpsize);

    printf("Write %zu -> %zu\n", size, zfpsize);
}

void compressrd_int(void *data, size_t size, int dim, int *shape){
    zfp_stream *zfp;
    int ret;
    zfp_type type = zfp_type_int32;                          
    zfp_field *field;
    bitstream *stream;
    size_t bufsize;  
    size_t zfpsize;
    void *buffer;

    // array metadata
    if (dim == 1){
        field = zfp_field_1d(data, type, (unsigned int)(shape[0]) );
    }
    else if (dim == 2){
        field = zfp_field_2d(data, type, (unsigned int)(shape[0]), (unsigned int)(shape[1]));
    }
    else{
        int i;
        unsigned int totalsize = 1;

        for(i = 2; i < dim; i++){
            totalsize *= shape[i];
        }
        field = zfp_field_3d(data, type, (unsigned int)(shape[0]), (unsigned int)(shape[1]), totalsize);
    }

    // compressed stream and parameters
    zfp = zfp_stream_open(NULL);   

    // set tolerance for fixed-accuracy mode           
    zfp_stream_set_accuracy(zfp, 0);                  
    //precision = zfp_stream_set_precision(zfp, 16);           
    //zfp_stream_set_rate(zfp, rate, type, 3, 0);       

    // allocate buffer for compressed data                     
    fread(&bufsize, sizeof(bufsize), 1, fh);
    //read(fd, &bufsize, sizeof(bufsize));
    buffer = malloc(bufsize);   
    fread(buffer, bufsize, 1, fh);
    //read(fd, buffer, bufsize);

    // associate bit stream with allocated buffer
    stream = stream_open(buffer, bufsize);      
    zfp_stream_set_bit_stream(zfp, stream);                  
    zfp_stream_rewind(zfp);                  

    ret = zfp_decompress(zfp, field);

    if (ret > 0) {
        printf("Read %zu -> %zu\n", bufsize, size);
    }
    else {
        printf("Decompress fail: addr: %llx, size: %zu\n", data, size);
    }
}

void compresswr_real_(double *R, int *size, int *dim, int *shape ) {
    //printf("Write %d bytes from %llx\n", *size, R);
    if (*dim > 0){
        compresswr_double((void*)R, (size_t)(*size), *dim, shape);
    }
    else{
        fwrite(R, (size_t)(*size), 1, fh);
        //write(fd, R, (size_t)(*size));
    }
}

void compressrd_real_(double *D, int *size, int *dim, int *shape  ) {
    //printf("Read %d bytes to %llx\n", *size, D);
    if (*dim > 0){
        compressrd_double((void*)D, (size_t)(*size), *dim, shape);
    }
    else{
        fread(D, (size_t)(*size), 1, fh);
        //read(fd, D, (size_t)(*size));
    }
}


void compresswr_integer_(int *R, int *size, int *dim, int *shape  ) {
    //printf("Write %d bytes from %llx\n", *size, R);
    //write(fd, R, (size_t)(*size));
    //compresswr((void*)R, (size_t)(*size));
    if (*dim > 0){
        compresswr_int((void*)R, (size_t)(*size), *dim, shape);
    }
    else{
        fwrite(R, (size_t)(*size), 1, fh);
        //write(fd, R, (size_t)(*size));
    }
}

void compressrd_integer_(int *D, int *size, int *dim, int *shape  ) {
    //printf("Read %d bytes to %llx\n", *size, D);
    //read(fd, D, (size_t)(*size));
    if (*dim > 0){
        compressrd_int((void*)D, (size_t)(*size), *dim, shape);
    }
    else{
        fread(D, (size_t)(*size), 1, fh);
        //read(fd, D, (size_t)(*size));
    }
}


void compresswr_bool_(int *R, int *size, int *dim, int *shape  ) {
    //printf("Write %d bytes from %llx\n", *size, R);
    //write(fd, R, (size_t)(*size));
    //compresswr((void*)R, (size_t)(*size));
    if (*dim > 0){
        compresswr_int((void*)R, (size_t)(*size), *dim, shape);
    }
    else{
        fwrite(R, (size_t)(*size), 1, fh);
        //write(fd, R, (size_t)(*size));
    }
}

void compressrd_bool_(int *D, int *size, int *dim, int *shape  ) {
    //printf("Read %d bytes to %llx\n", *size, D);
    //read(fd, D, (size_t)(*size));
    if (*dim > 0){
        compressrd_int((void*)D, (size_t)(*size), *dim, shape);
    }
    else{
        fread(D, (size_t)(*size), 1, fh);
        //read(fd, D, (size_t)(*size));
    }
}


void compresswr_string_(char *R, int *size , long l ) {
    //printf("Write %d bytes from %llx\n", *size, R);
    fwrite(R, (size_t)(*size), 1, fh);
    //write(fd, R, (size_t)(*size));
}

void compressrd_string_(char *D, int *size , long l ) {
    //printf("Read %d bytes to %llx\n", *size, D);
    fread(D, (size_t)(*size), 1, fh);
    //read(fd, D, (size_t)(*size));
}

