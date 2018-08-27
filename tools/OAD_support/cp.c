

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
int fd;

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

    fd = open(fname, O_CREAT | O_WRONLY, 0644);

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

    fd = open(fname, O_CREAT | O_RDONLY, 0644);
}

void cpc_close_(){
    close(fd);
}

void compresswr(void *data, size_t size){
    uint precision;
    // input: (void* array, int nx, int ny, int nz, double tolerance)

    // initialize metadata for the 3D array a[nz][ny][nx]
    zfp_type type = zfp_type_double;                          // array scalar type
    zfp_field* field = zfp_field_1d(data, type, (unsigned int)size/8); // array metadata
    //zfp_field* field = zfp_field_3d(array, type, nx, ny, nz); // array metadata

    // initialize metadata for a compressed stream
    zfp_stream* zfp = zfp_stream_open(NULL);                  // compressed stream and parameters
    zfp_stream_set_accuracy(zfp, 0.1);                  // set tolerance for fixed-accuracy mode
    //precision = zfp_stream_set_precision(zfp, 16);             // alternative: fixed-precision mode
    //zfp_stream_set_rate(zfp, rate, type, 3, 0);           // alternative: fixed-rate mode
    //printf("Write %zu, precision %u\n", precision);

    // allocate buffer for compressed data
    size_t bufsize = zfp_stream_maximum_size(zfp, field);     // capacity of compressed buffer (conservative)
    void* buffer = malloc(bufsize);                           // storage for compressed stream
    printf("Write %zu -> (estimate) %zu\n", size, bufsize);

    // associate bit stream with allocated buffer
    bitstream* stream = stream_open(buffer, bufsize);         // bit stream to compress to
    zfp_stream_set_bit_stream(zfp, stream);                   // associate with compressed stream
    zfp_stream_rewind(zfp);                                   // rewind stream to beginning

    // compress array
    size_t zfpsize = zfp_compress(zfp, field);                // return value is byte size of compressed stream

    printf("Write %zu -> %zu\n", size, zfpsize);
    write(fd, &zfpsize, sizeof(zfpsize));
    write(fd, buffer, zfpsize);
}

void compressrd(void *data, size_t size){
    // input: (void* array, int nx, int ny, int nz, double tolerance)

    // initialize metadata for the 3D array a[nz][ny][nx]
    zfp_type type = zfp_type_double;                          // array scalar type
    zfp_field* field = zfp_field_1d(data, type, (unsigned int)size/8); // array metadata
    //zfp_field* field = zfp_field_3d(array, type, nx, ny, nz); // array metadata

    // initialize metadata for a compressed stream
    zfp_stream* zfp = zfp_stream_open(NULL);                  // compressed stream and parameters
    zfp_stream_set_accuracy(zfp, 0.1);                  // set tolerance for fixed-accuracy mode
    //zfp_stream_set_precision(zfp, 16);             // alternative: fixed-precision mode
    //zfp_stream_set_rate(zfp, rate, type, 3, 0);           // alternative: fixed-rate mode

    // allocate buffer for compressed data
    size_t bufsize;
    read(fd, &bufsize, sizeof(bufsize));
    void* buffer = malloc(bufsize);                           // storage for compressed stream
    read(fd, buffer, bufsize);

    // associate bit stream with allocated buffer
    bitstream* stream = stream_open(buffer, bufsize);         // bit stream to compress to
    zfp_stream_set_bit_stream(zfp, stream);                   // associate with compressed stream
    zfp_stream_rewind(zfp);                                   // rewind stream to beginning

    int ret = zfp_decompress(zfp, field);

    if (ret > 0) {
        printf("Read %zu -> %zu\n", bufsize, size);
    }
    else {
        printf("Error: decompress fail\n");
    }
}

void compresswr_real_(double *R, int *size  ) {
    //printf("Write %d bytes from %llx\n", *size, R);
    compresswr((void*)R, (size_t)(*size));
}

void compressrd_real_(double *D, int *size  ) {
    //printf("Read %d bytes to %llx\n", *size, D);
    compressrd((void*)D, (size_t)(*size));
}


void compresswr_integer_(int *R, int *size  ) {
    printf("Write %d bytes from %llx\n", *size, R);
    write(fd, R, (size_t)(*size));
    //compresswr((void*)R, (size_t)(*size));
}

void compressrd_integer_(int *D, int *size  ) {
    printf("Read %d bytes to %llx\n", *size, D);
    read(fd, D, (size_t)(*size));
}


void compresswr_bool_(int *R, int *size  ) {
    //printf("Write %d bytes from %llx\n", *size, R);
    write(fd, R, (size_t)(*size));
}

void compressrd_bool_(int *D, int *size  ) {
    //printf("Read %d bytes to %llx\n", *size, D);
    read(fd, D, (size_t)(*size));
}


void compresswr_string_(char *R, int *size , long l ) {
    //printf("Write %d bytes from %llx\n", *size, R);
    write(fd, R, (size_t)(*size));
}

void compressrd_string_(char *D, int *size , long l ) {
    //printf("Read %d bytes to %llx\n", *size, D);
    read(fd, D, (size_t)(*size));
}

