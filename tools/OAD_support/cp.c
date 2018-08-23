

#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <linux/limits.h>
#include <fcntl.h>
#ifdef ALLOW_USE_MPI
#include <mpi.h>
#endif
#include "zlib.h"

#define BSIZE 1048576

static int cp_file_num = 0;
int fd;
int dsize;
size_t bsize;
char *buffer;

void buffer_init(){
    bsize = BSIZE;
    buffer = (char*)malloc(bsize);
}

void buffer_free(){
    bsize = 0;
    free(buffer);
}

void buffer_resize(size_t size){
    if (size > bsize){
        while(bsize < size){
            bsize <<= 1;
        }
        buffer = (char*)realloc(buffer, bsize);
    }
}   

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

    buffer_init();
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

    buffer_init();

    read(fd, &dsize, sizeof(dsize));
}

void cpc_close_(){
    close(fd);

    buffer_free();
}


inline void compresswr(double *R, int *size) {
    int err;

    buffer_resize((size_t)*size);

    // zlib struct
    z_stream defstream;
    defstream.zalloc = Z_NULL;
    defstream.zfree = Z_NULL;
    defstream.opaque = Z_NULL;
    defstream.avail_in = (uInt)(*size); // input size
    defstream.next_in = (Bytef *)R; // input
    defstream.avail_out = (uInt)bsize; // output buffer size
    defstream.next_out = (Bytef *)(buffer + sizeof(dsize)); // output buffer

    // the actual compression work.
    err = deflateInit(&defstream, Z_BEST_COMPRESSION);
    if (err < 0){
        printf("deflateInit fail: %d: %s\n", err, defstream.msg);
    }
    err = deflate(&defstream, Z_FINISH);
    if (err < 0){
        printf("deflate fail: %d: %s\n", err, defstream.msg);
    }
    err = deflateEnd(&defstream);
    if (err < 0){
        printf("deflateEnd fail: %d: %s\n", err, defstream.msg);
    }

    // Size of variable
    dsize = (int)defstream.total_out;
    *((int*)buffer) = dsize;

    // Write to file
    write(fd, buffer, defstream.total_out + sizeof(dsize));
 
    printf("Write: %d -> %d\n", *size, dsize);
}

inline void compressrd(double *D, int *size) {
    int err;

    buffer_resize((size_t)*size);

    // Read from file
    read(fd, buffer, dsize + sizeof(dsize));

    // zlib struct
    z_stream infstream;
    infstream.zalloc = Z_NULL;
    infstream.zfree = Z_NULL;
    infstream.opaque = Z_NULL;
    infstream.avail_in = (unsigned long) dsize; // input size
    infstream.next_in = (Bytef *)buffer; // input
    infstream.avail_out = (uInt)(*size); // output buffer
    infstream.next_out = (Bytef *)D; // buffer size
     
    // the actual DE-compression work.
    err = inflateInit(&infstream);
    if (err < 0){
        printf("inflateInit fail: %d: %s\n", err, infstream.msg);
    }
    err = inflate(&infstream, Z_NO_FLUSH);
    if (err < 0){
        printf("inflate fail: %d: %s\n", err, infstream.msg);
    }
    err = inflateEnd(&infstream);
    if (err < 0){
        printf("inflateEnd fail: %d: %s\n", err, infstream.msg);
    }
    if (*size != infstream.total_out){
        printf("Size mismatch: origin: %d, decompress: %lu\n", *size, infstream.total_out);
    }

    printf("Read: %lu -> %lu\n", dsize, infstream.total_out);
    
    // Size of next variable
    dsize = *((int*)(buffer + dsize));
}







void compresswr_real_(double *R, size_t *size  ) {
    compresswr(R, size);
}

void compressrd_real_(double *D, size_t *size  ) {
    compressrd(D, size);
}


void compresswr_integer_(int *R, size_t *size  ) {
    compresswr(R, size);
}

void compressrd_integer_(int *D, size_t *size  ) {
    compressrd(D, size);
}


void compresswr_bool_(int *R, size_t *size  ) {
    compresswr(R, size);
}

void compressrd_bool_(int *D, size_t *size  ) {
    compressrd(D, size);
}


void compresswr_string_(char *R, size_t *size , long l ) {
    compresswr(R, size);
}

void compressrd_string_(char *D, size_t *size , long l ) {
    compressrd(D, size);
}


