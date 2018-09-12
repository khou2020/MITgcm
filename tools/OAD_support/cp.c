

#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <linux/limits.h>
#include <fcntl.h>
#ifdef ALLOW_USE_MPI
#include <mpi.h>
#endif
#include "zlib.h"
#include <time.h>

#define BSIZE 1048576

static int cp_file_num = 0;
int fd;
int wr;
int dsize;
size_t bsize;
char *buffer;

static double compress_time, compress_time_old;
static double decompress_time, decompress_time_old;
static double wr_time, wr_time_old;
static double rd_time, rd_time_old;

clock_t topen;

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

    buffer_init();

    sprintf(fname, "oad_cp.%03d.%05d", rank, cp_file_num);

    if (num == NULL){
        cp_file_num++;
    }
    wr = 1;

    topen = clock();

    fd = open(fname, O_CREAT | O_WRONLY, 0644);
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

    buffer_init();
    
    sprintf(fname, "oad_cp.%03d.%05d", rank, cp_file_num);

    wr = 0;

    topen = clock();

    fd = open(fname, O_CREAT | O_RDONLY, 0644);

    read(fd, &dsize, sizeof(dsize));

    wr = 0;
}

void cpc_close_(){
    int rank;
    char fname[PATH_MAX];
    struct stat st;
    clock_t tclose;

    close(fd);

    tclose = clock();

    buffer_free();
    
    if (wr){
#ifdef ALLOW_USE_MPI
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#else
        rank = 0;
#endif
        sprintf(fname, "oad_cp.%03d.%05d", rank, cp_file_num);
        stat(fname, &st);
        printf("#%%$: CP_Size_%d: %lld\n", cp_file_num, (long long)st.st_size);

        printf("#%%$: CP_Com_Time_%d: %lf\n", cp_file_num, compress_time - compress_time_old);
        printf("#%%$: CP_Wr_Time_%d: %lf\n", cp_file_num, wr_time - wr_time_old); 

        printf("#%%$: CP_Store_Time_%d: %lf\n", cp_file_num, (double)(tclose - topen) / CLOCKS_PER_SEC); 

        compress_time_old = compress_time;
        wr_time_old = wr_time;
    }
    else{
        printf("#%%$: CP_Decom_Time_%d: %lf\n", cp_file_num, decompress_time - decompress_time_old);
        printf("#%%$: CP_Rd_Time_%d: %lf\n", cp_file_num, rd_time - rd_time_old); 

        printf("#%%$: CP_Restore_Time_%d: %lf\n", cp_file_num, (double)(tclose - topen) / CLOCKS_PER_SEC); 

        decompress_time_old = decompress_time;
        rd_time_old = rd_time;
    }

}

void cpc_profile_(){
    printf("#%%$: CP_Com_Time_All: %lf\n", compress_time);
    printf("#%%$: CP_Wr_Time_All: %lf\n", wr_time); 
    printf("#%%$: CP_Decom_Time_All: %lf\n", decompress_time);
    printf("#%%$: CP_Rd_Time_All: %lf\n", rd_time); 
}

inline void compresswr(void *R, int *size) {
    int err;
    clock_t t1, t2, t3;

    t1 = clock();

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

    t2 = clock();

    // Size of variable
    dsize = (int)defstream.total_out;
    *((int*)buffer) = dsize;

    // Write to file
    write(fd, buffer, defstream.total_out + sizeof(dsize));

    t3 = clock();
    
    printf("Write: %d -> %d\n", *size, dsize);

    compress_time += (double)(t2 - t1) / CLOCKS_PER_SEC;
    wr_time += (double)(t3 - t2) / CLOCKS_PER_SEC;
}

inline void compressrd(void *D, int *size) {
    int err;
    clock_t t1, t2, t3;

    t1 = clock();

    buffer_resize((size_t)*size);

    // Read from file
    read(fd, buffer, dsize + sizeof(dsize));

    t2 = clock();

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

    t3 = clock();

    printf("Read: %lu -> %lu\n", dsize, infstream.total_out);

    // Size of next variable
    dsize = *((int*)(buffer + dsize));

    decompress_time += (double)(t3 - t2) / CLOCKS_PER_SEC;
    rd_time += (double)(t2 - t1) / CLOCKS_PER_SEC;
}







void compresswr_real_(double *R, int *size  ) {
    compresswr(R, size);
}

void compressrd_real_(double *D, int *size  ) {
    compressrd(D, size);
}


void compresswr_integer_(int *R, int *size  ) {
    compresswr(R, size);
}

void compressrd_integer_(int *D, int *size  ) {
    compressrd(D, size);
}


void compresswr_bool_(int *R, int *size  ) {
    compresswr(R, size);
}

void compressrd_bool_(int *D, int *size  ) {
    compressrd(D, size);
}


void compresswr_string_(char *R, int *size , long l ) {
    compresswr(R, size);
}

void compressrd_string_(char *D, int *size , long l ) {
    compressrd(D, size);
}


