
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


void compresswr_real_(double *R, int *size) {
    // STEP 1.
    // deflate a into b. (that is, compress a into b)
    
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
    deflateInit(&defstream, Z_BEST_COMPRESSION);
    deflate(&defstream, Z_FINISH);
    deflateEnd(&defstream);
     
    // This is one way of getting the size of the output
    printf("Write: %d -> %lu\n", *size, defstream.total_out);
    write(fd, &(defstream.total_out), sizeof(defstream.total_out));
    write(fd, buffer, defstream.total_out);
}

void compressrd_real_(double *D, int *size) {
    // STEP 2.
    // inflate b into c
    // zlib struct
    z_stream infstream;
    infstream.zalloc = Z_NULL;
    infstream.zfree = Z_NULL;
    infstream.opaque = Z_NULL;
    // setup "b" as the input and "c" as the compressed output
    //infstream.avail_in = (uInt)((char*)defstream.next_out - b); // size of input
    read(fd, &(infstream.avail_in), sizeof(infstream.avail_in));
    read(fd, buffer, (size_t)infstream.avail_in);
    infstream.next_in = (Bytef *)buffer; // input char array
    infstream.avail_out = (uInt)(*size); // size of output
    infstream.next_out = (Bytef *)D; // output char array
     
    // the actual DE-compression work.
    inflateInit(&infstream);
    inflate(&infstream, Z_NO_FLUSH);
    inflateEnd(&infstream);
    
    printf("Read: %lu -> %d\n",  infstream.total_out, *size);
}

void compresswr_int_(int *R, int *size) {
    //printf("Write %d bytes from %llx\n", *size, R);
    //write(fd, R, *size);
}

void compressrd_int_(int *D, int *size) {
    //printf("Read %d bytes to %llx\n", *size, D);
    //read(fd, D, *size);
}

void compresswr_bool_(int *R, int *size) {
    write(fd, R, *size);
}

void compressrd_bool_(int *D, int *size) {
    //printf("Read %d bytes\n", *size);
    read(fd, D, *size);
}

void compresswr_string_(int *R, int *size) {
    write(fd, R, *size);
}

void compressrd_string_(int *D, int *size) {
    //printf("Read %d bytes\n", *size);
    read(fd, D, *size);
}