

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
#include "zfp.h"
#include <string.h>
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#define BSIZE (1 * 1024 * 1024)
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

    fd->fd = open(fname, O_CREAT | O_WRONLY | O_TRUNC | O_DIRECT | O_SYNC, 0644);

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

    wsize = 0;
    while(fd->abuf < fd->cbuf){
        ioret = write(fd->fd, fd->abuf, fd->cbuf - fd->abuf);
        if (ioret <= 0){
            ret = -1;
        }
        fd->abuf += ioret;
    }

    wr_time = getwalltime() - t1;

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

    fd->fd = open(fname, O_RDONLY | O_DIRECT | O_SYNC, 0644);

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

    if (*num > 0){
        cp_file_num = *num;
    }

#ifdef ALLOW_USE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#else
    rank = 0;
#endif

    buffer_init();

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
    zfp_stream *zfp;
    zfp_type type = zfp_type_double;                          
    zfp_field *field;
    bitstream *stream;
    size_t bufsize;  
    size_t zfpsize;
    double t1, t2, t3;

    t1 = getwalltime();

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
    zfp_stream_set_accuracy(zfp, 0.0001);                     

    // allocate buffer for compressed data
    bufsize = zfp_stream_maximum_size(zfp, field);    
    buffer_resize((size_t)bufsize);              

    // associate bit stream with allocated buffer
    stream = stream_open((void*)buffer, bufsize);      
    zfp_stream_set_bit_stream(zfp, stream);                  
    zfp_stream_rewind(zfp);                  

    // compress array
    zfpsize = zfp_compress(zfp, field);               

    t2 = getwalltime();

    cp_write(fd, &zfpsize, sizeof(zfpsize));
    cp_write(fd, (void*)buffer, zfpsize);

    t3 = getwalltime();

    compress_time += t2 - t1;

    //printf("Write %zu -> %zu\n", size, zfpsize);
}

void compressrd_real(void *data, size_t size, int dim, int *shape){
    zfp_stream *zfp;
    int ret;
    zfp_type type = zfp_type_double;                          
    zfp_field *field;
    bitstream *stream;
    size_t bufsize;  
    size_t zfpsize;
    double t1, t2, t3;
    
    t1 = getwalltime();

    // allocate buffer for compressed data                     
    cp_read(fd, &bufsize, sizeof(bufsize));
    buffer_resize((size_t)bufsize);   
    cp_read(fd, (void*)buffer, bufsize);

    t2 = getwalltime();

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
    zfp_stream_set_accuracy(zfp, 0.0001);                      


    // associate bit stream with allocated buffer
    stream = stream_open((void*)buffer, bufsize);      
    zfp_stream_set_bit_stream(zfp, stream);                  
    zfp_stream_rewind(zfp);                  

    ret = zfp_decompress(zfp, field);

    t3 = getwalltime();

    decompress_time += t3 - t2;
    
    if (ret < 0) {
        printf("Decompress fail: addr: %llx, size: %zu\n", data, size);   
    }
    else {
        //printf("Read %zu -> %zu\n", bufsize, size);
    }
}

/*
void compresswr_int(void *data, size_t size, int dim, int *shape){
    zfp_stream *zfp;
    zfp_type type = zfp_type_int32;                          
    zfp_field *field;
    bitstream *stream;
    size_t bufsize;  
    size_t zfpsize;
    double t1, t2, t3;

    t1 = getwalltime();

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

    // allocate buffer for compressed data
    bufsize = zfp_stream_maximum_size(zfp, field);    
    buffer_resize((size_t)bufsize);              

    // associate bit stream with allocated buffer
    stream = stream_open((void*)buffer, bufsize);      
    zfp_stream_set_bit_stream(zfp, stream);                  
    zfp_stream_rewind(zfp);                  

    // compress array
    zfpsize = zfp_compress(zfp, field);               

    t2 = getwalltime();

    cp_write(fd, &zfpsize, sizeof(zfpsize));
    cp_write(fd, (void*)buffer, zfpsize);

    t3 = getwalltime();

    compress_time += t2 - t1;

    //printf("Write %zu -> %zu\n", size, zfpsize);
}

void compressrd_int(void *data, size_t size, int dim, int *shape){
    zfp_stream *zfp;
    int ret;
    zfp_type type = zfp_type_int32;                          
    zfp_field *field;
    bitstream *stream;
    size_t bufsize;  
    size_t zfpsize;
    double t1, t2, t3;
    
    t1 = getwalltime();

    // allocate buffer for compressed data                     
    cp_read(fd, &bufsize, sizeof(bufsize));
    buffer_resize((size_t)bufsize);   
    cp_read(fd, (void*)buffer, bufsize);

    t2 = getwalltime();

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


    // associate bit stream with allocated buffer
    stream = stream_open((void*)buffer, bufsize);      
    zfp_stream_set_bit_stream(zfp, stream);                  
    zfp_stream_rewind(zfp);                  

    ret = zfp_decompress(zfp, field);

    t3 = getwalltime();

    decompress_time += t3 - t2;
    
    if (ret < 0) {
        printf("Decompress fail: addr: %llx, size: %zu\n", data, size);   
    }
    else {
        //printf("Read %zu -> %zu\n", bufsize, size);
    }
}
*/







void compresswr_integer_(int *R, int *size, int *dim, int *shape ) {
    /*
    if (*size > THRESHOLD){
        compressrd_int((void*)D, (size_t)(*size), *dim, shape);
    }
    else{
    */
    double t1, t2;
    cp_write(fd, R, (size_t)(*size));
    //}
}

void compressrd_integer_(int *D, int *size, int *dim, int *shape  ) {
    /*
    if (*size > THRESHOLD){
        compressrd_int((void*)D, (size_t)(*size), *dim, shape);
    }
    else{
    */
    double t1, t2;
    cp_read(fd, D, (size_t)(*size));
    //}
}


void compresswr_bool_(int *R, int *size, int *dim, int *shape ) {
    /*
    if (*size > THRESHOLD){
        compressrd_int((void*)D, (size_t)(*size), *dim, shape);
    }
    else{
    */
    double t1, t2;
    cp_write(fd, R, (size_t)(*size));
    //}
}

void compressrd_bool_(int *D, int *size, int *dim, int *shape  ) {
    /*
    if (*size > THRESHOLD){
        compressrd_int((void*)D, (size_t)(*size), *dim, shape);
    }
    else{
    */
    double t1, t2;
    cp_read(fd, D, (size_t)(*size));
    //}
}



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
