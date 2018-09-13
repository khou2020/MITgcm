


#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <linux/limits.h>
#include <fcntl.h>
#ifdef ALLOW_USE_MPI
#include <mpi.h>
#endif
#include "zfp.h"
#include <time.h>
#define BSIZE 1048576
#define THRESHOLD 1024

static int cp_file_num = 0;
int cur_num;
int fd;
int wr;

size_t bsize;
char *buffer;

static double compress_time_all, decompress_time_all, wr_time_all, rd_time_all, store_time_all, restore_time_all;
double compress_time, decompress_time, wr_time, rd_time;

double topen;

double getwalltime(){
    struct timespec tp;
    clock_gettime(CLOCK_MONOTONIC , &tp);
    return (double)tp.tv_sec + (double)tp.tv_nsec / 1e+9f;
}

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
    cur_num = cp_file_num;

    if (*num <= 0){
        cp_file_num++;
    }
    wr = 1;

    compress_time = 0;
    decompress_time = 0;
    wr_time = 0;
    rd_time = 0;

    topen = getwalltime();

    fd = open(fname, O_CREAT | O_WRONLY | O_TRUNC, 0644);
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
    cur_num = cp_file_num;

    wr = 0;
    compress_time = 0;
    decompress_time = 0;
    wr_time = 0;
    rd_time = 0;

    topen = getwalltime();

    fd = open(fname, O_CREAT | O_RDONLY, 0644);
}

void cpc_close_(){
    int rank;
    char fname[PATH_MAX];
    struct stat st;
    double tclose;

    //buffer_free();
    
    if (wr){
        fsync(fd);
        close(fd);
        
        tclose = getwalltime();

#ifdef ALLOW_USE_MPI
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#else
        rank = 0;
#endif
        sprintf(fname, "oad_cp.%03d.%05d", rank, cur_num);
        stat(fname, &st);
        printf("#%%$: CP_Size_%d: %lld\n", cur_num, (long long)st.st_size);

        printf("#%%$: CP_Com_Time_%d: %lf\n", cur_num, compress_time);
        printf("#%%$: CP_Wr_Time_%d: %lf\n", cur_num, wr_time); 

        printf("#%%$: CP_Store_Time_%d: %lf\n", cur_num, tclose - topen); 

        compress_time_all += compress_time;
        wr_time_all += wr_time;
        store_time_all += tclose - topen;
    }
    else{
        close(fd);

        tclose = getwalltime();

        printf("#%%$: CP_Decom_Time_%d: %lf\n", cur_num, decompress_time);
        printf("#%%$: CP_Rd_Time_%d: %lf\n", cur_num, rd_time); 

        printf("#%%$: CP_Restore_Time_%d: %lf\n", cur_num, tclose - topen); 

        decompress_time_all += decompress_time;
        rd_time_all += rd_time;
        restore_time_all += tclose - topen;
    }

    buffer_free();
}

void cpc_profile_(){
    printf("#%%$: CP_Com_Time_All: %lf\n", compress_time_all);
    printf("#%%$: CP_Wr_Time_All: %lf\n", wr_time_all); 
    printf("#%%$: CP_Store_Time_All: %lf\n", store_time_all); 
    printf("#%%$: CP_Decom_Time_All: %lf\n", decompress_time_all);
    printf("#%%$: CP_Rd_Time_All: %lf\n", rd_time_all); 
    printf("#%%$: CP_Restore_Time_All: %lf\n", restore_time_all); 
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
    zfp_stream_set_accuracy(zfp, 0.000001);                     

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

    write(fd, &zfpsize, sizeof(zfpsize));
    write(fd, (void*)buffer, zfpsize);

    t3 = getwalltime();

    compress_time += t2 - t1;
    wr_time += t3 - t2;

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
    read(fd, &bufsize, sizeof(bufsize));
    buffer_resize((size_t)bufsize);   
    read(fd, (void*)buffer, bufsize);

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
    zfp_stream_set_accuracy(zfp, 0.000001);                      


    // associate bit stream with allocated buffer
    stream = stream_open((void*)buffer, bufsize);      
    zfp_stream_set_bit_stream(zfp, stream);                  
    zfp_stream_rewind(zfp);                  

    ret = zfp_decompress(zfp, field);

    t3 = getwalltime();

    decompress_time += t3 - t2;
    rd_time += t2 - t1;
    
    if (ret < 0) {
        printf("Decompress fail: addr: %llx, size: %zu\n", data, size);   
    }
    else {
        //printf("Read %zu -> %zu\n", bufsize, size);
    }
}


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

    write(fd, &zfpsize, sizeof(zfpsize));
    write(fd, (void*)buffer, zfpsize);

    t3 = getwalltime();

    compress_time += t2 - t1;
    wr_time += t3 - t2;

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
    read(fd, &bufsize, sizeof(bufsize));
    buffer_resize((size_t)bufsize);   
    read(fd, (void*)buffer, bufsize);

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
    rd_time += t2 - t1;
    
    if (ret < 0) {
        printf("Decompress fail: addr: %llx, size: %zu\n", data, size);   
    }
    else {
        //printf("Read %zu -> %zu\n", bufsize, size);
    }
}




void compresswr_real_(double *R, int *size, int *dim, int *shape ) {
    if (*size > THRESHOLD){
        compresswr_real((void*)R, (size_t)(*size), *dim, shape);
    }
    else{
        double t1, t2;
        t1 = getwalltime();
        write(fd, R, (size_t)(*size));
        t2 = getwalltime();
        wr_time += t2 - t1;
    }
}

void compressrd_real_(double *D, int *size, int *dim, int *shape  ) {
    if (*size > THRESHOLD){
        compressrd_real((void*)D, (size_t)(*size), *dim, shape);
    }
    else{
        double t1, t2;
        t1 = getwalltime();
        read(fd, D, (size_t)(*size));
        t2 = getwalltime();
        rd_time += t2 - t1;
    }
}


void compresswr_integer_(int *R, int *size, int *dim, int *shape ) {
    if (*size > THRESHOLD){
        compresswr_int((void*)R, (size_t)(*size), *dim, shape);
    }
    else{
        double t1, t2;
        t1 = getwalltime();
        write(fd, R, (size_t)(*size));
        t2 = getwalltime();
        wr_time += t2 - t1;
    }
}

void compressrd_integer_(int *D, int *size, int *dim, int *shape  ) {
    if (*size > THRESHOLD){
        compressrd_int((void*)D, (size_t)(*size), *dim, shape);
    }
    else{
        double t1, t2;
        t1 = getwalltime();
        read(fd, D, (size_t)(*size));
        t2 = getwalltime();
        rd_time += t2 - t1;
    }
}


void compresswr_bool_(int *R, int *size, int *dim, int *shape ) {
    if (*size > THRESHOLD){
        compresswr_int((void*)R, (size_t)(*size), *dim, shape);
    }
    else{
        double t1, t2;
        t1 = getwalltime();
        write(fd, R, (size_t)(*size));
        t2 = getwalltime();
        wr_time += t2 - t1;
    }
}

void compressrd_bool_(int *D, int *size, int *dim, int *shape  ) {
    if (*size > THRESHOLD){
        compressrd_int((void*)D, (size_t)(*size), *dim, shape);
    }
    else{
        double t1, t2;
        t1 = getwalltime();
        read(fd, D, (size_t)(*size));
        t2 = getwalltime();
        rd_time += t2 - t1;
    }
}




void compresswr_string_(char *R, int *size , long l ) {
        double t1, t2;
        t1 = getwalltime();
        write(fd, R, (size_t)(*size));
        t2 = getwalltime();
        wr_time += t2 - t1;
}

void compressrd_string_(char *D, int *size , long l ) {
    double t1, t2;
    t1 = getwalltime();
    read(fd, D, (size_t)(*size));
    t2 = getwalltime();
    rd_time += t2 - t1;
}
