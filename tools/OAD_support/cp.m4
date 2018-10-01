dnl Process this m4 file to produce ]C] language file.
dnl
dnl If you see this line, you can ignore the next one.
dnl /* Do not edit this file. It is produced from the corresponding .m4 source */
dnl
dnl

#define _GNU_SOURCE

#include <stdlib.h>
#include <stdio.h>
#include <linux/limits.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#ifdef ALLOW_USE_MPI
#include <mpi.h>
#endif
#include "sz.h"
#include <string.h>
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#define BSIZE (100 * 1024 * 1024)
#define FSIZE (100 * 1024 * 1024)
#define MAXITR 1024
#define THRESHOLD 1024

typedef struct cp_fd{
    int fd;
    char *buf, *abuf, *cbuf;
} cp_fd;

static int cp_file_num = 0;
cp_fd *fd;

int cur_num;
int wr;

size_t bsize, bused;
void *buffer, *abuffer;

double compress_time, decompress_time, wr_time, rd_time;

double topen, tclose;
double times[MAXITR][10];
unsigned long long fsize[MAXITR];
static int max_itr;

double getwalltime(){
    struct timespec tp;
    clock_gettime(CLOCK_MONOTONIC , &tp);
    return (double)tp.tv_sec + (double)tp.tv_nsec / 1e+9f;
}

void *iobuf;
size_t iobsize;

void buffer_init(){
    int pagesize;

    pagesize = getpagesize();
    buffer = malloc(BSIZE + pagesize);
    abuffer=(void*)((((unsigned long long)buffer+(unsigned long long)pagesize-1)/(unsigned long long)pagesize)*(unsigned long long)pagesize);

    bsize = BSIZE;
    buffer = (float*)malloc(bsize);
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
        buffer = (float*)realloc(buffer, bsize);
    }
} 

cp_fd* cp_wr_open(char* fname, size_t fsize){
    int pagesize;
    cp_fd *fd;

    fd = malloc(sizeof(cp_fd));

    pagesize = getpagesize();
    fd->buf  = (char*)malloc(FSIZE + pagesize);
    fd->abuf = fd->cbuf = (char*)((((size_t)fd->buf + (size_t)pagesize - 1) / (size_t)pagesize) * (size_t)pagesize);

    topen = getwalltime();

    fd->fd = open(fname, O_CREAT | O_WRONLY | O_TRUNC, 0644);

    return fd;
}

int cp_write(cp_fd *fd, void *data, size_t size){
    memcpy(fd->cbuf, data, size);
    fd->cbuf += size;
}

int cp_wr_close(cp_fd *fd){
    int ret = 0;
    off_t wsize;
    ssize_t ioret;
    double t1, t2, t3;

#ifdef ALLOW_USE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif

    t1 = getwalltime();
    fsync(fd->fd);
    t2 = getwalltime();
    
#ifdef ALLOW_USE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif

    t3 = getwalltime();

    wsize = 0;
    while(fd->abuf < fd->cbuf){
        ioret = write(fd->fd, fd->abuf, fd->cbuf - fd->abuf);
        if (ioret <= 0){
            ret = -1;
        }
        fd->abuf += ioret;
    }
    fsync(fd->fd);

    wr_time = getwalltime() - t3 - t2 + t1;

    close(fd->fd);

    tclose = getwalltime();
    
    free(fd->buf);
    free(fd);

    return ret;
}

cp_fd* cp_rd_open(char* fname){
    int pagesize;
    cp_fd *fd;
    struct stat st;
    double t1;

    fd = malloc(sizeof(cp_fd));

    stat(fname, &st);

    pagesize = getpagesize();
    fd->buf  = (char*)malloc(FSIZE + pagesize);
    fd->abuf = fd->cbuf = (char*)((((size_t)fd->buf + (size_t)pagesize - 1) / (size_t)pagesize) * (size_t)pagesize);

    topen = getwalltime();

    fd->fd = open(fname, O_RDONLY, 0644);

#ifdef ALLOW_USE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif

    t1 = getwalltime();

    read(fd->fd, fd->abuf, st.st_size);

    rd_time = getwalltime() - t1;

    return fd;
}

int cp_read(cp_fd *fd, void *data, size_t size){
    memcpy(data, fd->cbuf, size);
    fd->cbuf += size;
}

int cp_rd_close(cp_fd *fd){   
    close(fd->fd);

    tclose = getwalltime();
    
    free(fd->buf);
    free(fd);

    return 0;
}

void cp_wr_open_(int *num){
    int rank;
    char fname[PATH_MAX]; 
    char *envfname;
    sz_params sz;

    if (*num > 0){
        cp_file_num = *num;
    }

#ifdef ALLOW_USE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#else
    rank = 0;
#endif

    buffer_init();

    memset(&sz, 0, sizeof(sz_params));
    sz.sol_ID = SZ;
    sz.sampleDistance = 100;
    sz.quantization_intervals = 0;
    sz.max_quant_intervals = 65536;
    sz.predThreshold = 0.98;
    sz.szMode = SZ_BEST_SPEED;
    sz.losslessCompressor = ZSTD_COMPRESSOR;
    sz.gzipMode = 1;
    sz.errorBoundMode = ABS;
    sz.absErrBound = 1E-3;
    sz.relBoundRatio = 1E-5;

    SZ_Init_Params(&sz);

    envfname = getenv("MITGCM_OAD_CP_PREFIX");
    if (envfname == NULL){
        envfname = "oad_cp";
    }
    sprintf(fname, "%s.%03d.%05d", envfname, rank, cp_file_num);
    cur_num = cp_file_num;

    if (*num <= 0){
        cp_file_num++;
    }
    wr = 1;

    compress_time = 0;
    decompress_time = 0;
    wr_time = 0;
    rd_time = 0;

    fd = cp_wr_open(fname, FSIZE);
}

void cp_rd_open_(int *num){
    int rank;
    char fname[PATH_MAX];
    char *envfname;
    sz_params sz;

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

    memset(&sz, 0, sizeof(sz_params));
    sz.sol_ID = SZ;
    sz.sampleDistance = 100;
    sz.quantization_intervals = 0;
    sz.max_quant_intervals = 65536;
    sz.predThreshold = 0.98;
    sz.szMode = SZ_BEST_SPEED;
    sz.losslessCompressor = ZSTD_COMPRESSOR;
    sz.gzipMode = 1;
    sz.errorBoundMode = ABS;
    sz.absErrBound = 1E-3;
    sz.relBoundRatio = 1E-5;

    SZ_Init_Params(&sz);
    
    envfname = getenv("MITGCM_OAD_CP_PREFIX");
    if (envfname == NULL){
        envfname = "oad_cp";
    }
    sprintf(fname, "%s.%03d.%05d", envfname, rank, cp_file_num);
    cur_num = cp_file_num;

    wr = 0;
    compress_time = 0;
    decompress_time = 0;
    wr_time = 0;
    rd_time = 0;

    fd = cp_rd_open(fname);
}

void cpc_close_(){
    int rank, np;
    char fname[PATH_MAX];
    struct stat st;
    double tio;
    unsigned long long size;

    //buffer_free();
    
    if (wr){
        cp_wr_close(fd);

#ifdef ALLOW_USE_MPI
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &np);
#else
        rank = 0;
        np = 1;
#endif
        sprintf(fname, "oad_cp.%03d.%05d", rank, cur_num);
        stat(fname, &st);
        size = (unsigned long long)st.st_size;

        tio = tclose - topen;

        MPI_Reduce(&size, &(fsize[cur_num]), 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(&tio, &(times[cur_num][0]), 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
        MPI_Reduce(&compress_time, &(times[cur_num][1]), 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(&wr_time, &(times[cur_num][2]), 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        times[cur_num][1] /= np;
        times[cur_num][2] /= np;
    }
    else{
        cp_rd_close(fd);

        tio = tclose - topen;

#ifdef ALLOW_USE_MPI
        MPI_Comm_size(MPI_COMM_WORLD, &np);
#else
        np = 1;
#endif

        MPI_Reduce(&tio, &(times[cur_num][3]), 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
        MPI_Reduce(&decompress_time, &(times[cur_num][4]), 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(&rd_time, &(times[cur_num][5]), 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        times[cur_num][4] /= np;
        times[cur_num][5] /= np;
    }

    if (max_itr < cur_num){
        max_itr = cur_num;
    }

    buffer_free();

    SZ_Finalize();
}

void cpc_profile_(){
    int rank;
    int i, j;
    FILE* f;
    char *fname;

#ifdef ALLOW_USE_MPI
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#else
        rank = 0;
#endif

    if (rank == 0){
        printf("Itr,\tstore_time,\tcom_time,\twr_time,\trestore_time,\tdecom_time,\trd_time,\tfsize\n");
        for(i = 0; i <= max_itr; i++){
            printf("%d,\t", i);
            for (j = 0; j < 6; j++){
                printf("%lf,\t", times[i][j]);
            }
            printf("%llu\n", fsize[i]);
        }
        fname = getenv("MITGCM_PROFILE_PATH");
        if (fname != NULL){
            f = fopen(fname, "w");
        }
        else{
            f = fopen("profile_origin.csv", "w");
        }
        fprintf(f, "Itr,\tstore_time,\tcom_time,\twr_time,\trestore_time,\tdecom_time,\trd_time,\tfsize\n");
        for(i = 0; i <= max_itr; i++){
            fprintf(f, "%d,\t", i);
            for (j = 0; j < 6; j++){
                fprintf(f, "%lf,\t", times[i][j]);
            }
            fprintf(f, "%llu\n", fsize[i]);
        }
        fclose(f);
    }
}






void compresswr_real(void *data, size_t size, int dim, int *shape){
    int i;
    double t1, t2, t3;
    unsigned char *buf;
    size_t outsize;
    size_t r[4];

    t1 = getwalltime();

    for(i = 0; i < 4; i++){
        if (i < dim){
            r[i] = shape[i];
        }
        else{
            r[i] = 0;
        }
    }
    for(i = 4; i < dim; i++){
        r[3] *= shape[i];
    }
    
    for(i = 3; i > -1; i--){
        if (r[i] == 1){
            r[i] = 0;
        }
        else {
            break;
        }
    }

    buf = SZ_compress(SZ_DOUBLE, data, &outsize, 0, r[3], r[2], r[1], r[0]);

    t2 = getwalltime();

    cp_write(fd, &outsize, sizeof(outsize));
    cp_write(fd, (void*)buf, outsize);

    t3 = getwalltime();

    compress_time += t2 - t1;
    wr_time += t3 - t2;

    free(buf);

    //printf("Write %zu -> %zu\n", size, zfpsize);
}

void compressrd_real(void *data, size_t size, int dim, int *shape){
    int i;
    int ret = 0;
    double t1, t2, t3;
    unsigned char *out;
    size_t outsize, bufsize;
    size_t r[4];
    
    t1 = getwalltime();

    // allocate buffer for compressed data                     
    cp_read(fd, &bufsize, sizeof(bufsize));
    buffer_resize((size_t)bufsize);   
    cp_read(fd, (void*)buffer, bufsize);

    t2 = getwalltime();

    for(i = 0; i < 4; i++){
        if (i < dim){
            r[i] = shape[i];
        }
        else{
            r[i] = 0;
        }
    }
    for(i = 4; i < dim; i++){
        r[3] *= shape[i];
    }
    for(i = 3; i > -1; i--){
        if (r[i] == 1){
            r[i] = 0;
        }
        else {
            break;
        }
    }

    out = SZ_decompress(SZ_DOUBLE, buffer, bufsize, 0, r[3], r[2], r[1], r[0]);
    if (out == NULL){
        ret = -1;
    }

    t3 = getwalltime();

    memcpy(data, out, size);
    free(out);

    decompress_time += t3 - t2;
    rd_time += t2 - t1;
    
    if (ret < 0) {
        printf("Decompress fail: addr: %llx, size: %zu\n", data, size);   
    }
    else {
        //printf("Read %zu -> %zu\n", bufsize, size);
    }
}

include(`foreach.m4')dnl
include(`forloop.m4')dnl

define(`CONCAT', `$1$2')dnl

define(`CMP_REAL',changequote(`[', `]')dnl
[dnl

void compresswr_$1_($2 *R, int *size, int *dim, int *shape ) {
    //printf("Write %d bytes from %llx\n", *size, R);
    cp_write(fd, R, (size_t)(*size));
}

void compressrd_$1_($2 *D, int *size, int *dim, int *shape) {
    cp_read(fd, D, (size_t)(*size));
}

]changequote([`], [']))dnl
dnl

foreach(`i', (`(`integer', `int', `int')',dnl
                `(`bool', `int', `int')'), `CMP_REAL(translit(i, `()'))')

void compresswr_real_(double *R, int *size, int *dim, int *shape ) {
    if (*size > THRESHOLD){
        compresswr_real((void*)R, (size_t)(*size), *dim, shape);
    }
    else{
        cp_write(fd, R, (size_t)(*size));
    }
}

void compressrd_real_(double *D, int *size, int *dim, int *shape  ) {
    if (*size > THRESHOLD){
        compressrd_real((void*)D, (size_t)(*size), *dim, shape);
    }
    else{
        cp_read(fd, D, (size_t)(*size));
    }
}

void compresswr_string_(char *R, int *size , long l ) {
    cp_write(fd, R, (size_t)(*size));
}

void compressrd_string_(char *D, int *size , long l ) {
    cp_read(fd, D, (size_t)(*size));
}
