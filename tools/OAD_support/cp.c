
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
size_t bsize;
void *buffer;

void buffer_init(){
    bsize = BSIZE;
    buffer = malloc(bsize);
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
        buffer = realloc(buffer, bsize);
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
}

void cpc_close_(){
    close(fd);

    buffer_free();
}


inline void compresswr(double *R, int *size) {
    int err;
    int csize;

    buffer_resize((size_t)*size);

    // zlib struct
    z_stream defstream;
    defstream.zalloc = Z_NULL;
    defstream.zfree = Z_NULL;
    defstream.opaque = Z_NULL;
    // setup "a" as the input and "b" as the compressed output
    defstream.avail_in = (uInt)(*size); // size of input, string + terminator
    defstream.next_in = (Bytef *)R; // input char array
    defstream.avail_out = (uInt)bsize; // size of output
    defstream.next_out = (Bytef *)buffer; // output char array
    
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

    // This is one way of getting the size of the output
    csize = (int)defstream.total_out;
    //printf("Write: %d -> %d\n", *size, csize);
    write(fd, &(csize), sizeof(csize));
    write(fd, buffer, defstream.total_out);
}

inline void compressrd(double *D, int *size) {
    int err;
    int dsize;

    buffer_resize((size_t)*size);

    // zlib struct
    z_stream infstream;
    infstream.zalloc = Z_NULL;
    infstream.zfree = Z_NULL;
    infstream.opaque = Z_NULL;
    // setup "b" as the input and "c" as the compressed output
    read(fd, &(dsize), sizeof(dsize));
    infstream.avail_in = (unsigned long) dsize;
    read(fd, buffer, (size_t)infstream.avail_in);
    infstream.next_in = (Bytef *)buffer; // input char array
    infstream.avail_out = (uInt)(*size); // size of output
    infstream.next_out = (Bytef *)D; // output char array
     
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
    
    //printf("Read: %lu -> %lu\n", dsize, infstream.total_out);
}





void compresswr_real_(double *R, size_t size) {
    compresswr(R, size);
}

void compressrd_real_(double *D, size_t size) {
    compressrd(D, size);
}


void compresswr_int_(int *R, size_t size) {
    compresswr(R, size);
}

void compressrd_int_(int *D, size_t size) {
    compressrd(D, size);
}


void compresswr_bool_(int *R, size_t size) {
    compresswr(R, size);
}

void compressrd_bool_(int *D, size_t size) {
    compressrd(D, size);
}


void compresswr_string_(char *R, size_t size) {
    compresswr(R, size);
}

void compressrd_string_(char *D, size_t size) {
    compressrd(D, size);
}

