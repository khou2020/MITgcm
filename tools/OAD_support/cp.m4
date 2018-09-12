dnl Process this m4 file to produce ]C] language file.
dnl
dnl If you see this line, you can ignore the next one.
dnl /* Do not edit this file. It is produced from the corresponding .m4 source */
dnl
dnl


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
int fd;
int wr;

size_t bsize;
char *buffer;

static double compress_time, compress_time_old;
static double decompress_time, decompress_time_old;
static double wr_time, wr_time_old;
static double rd_time, rd_time_old;
static double store_time, restore_time;

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

    wr = 0;

    topen = clock();

    fd = open(fname, O_CREAT | O_RDONLY, 0644);
}

void cpc_close_(){
    int rank;
    char fname[PATH_MAX];
    struct stat st;
    clock_t tclose;

    if (wr){
        fsync(fd);
        close(fd);

        tclose = clock();

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
        store_time += (double)(tclose - topen) / CLOCKS_PER_SEC;
    }
    else{
        close(fd);

        tclose = clock();

        printf("#%%$: CP_Decom_Time_%d: %lf\n", cp_file_num, decompress_time - decompress_time_old);
        printf("#%%$: CP_Rd_Time_%d: %lf\n", cp_file_num, rd_time - rd_time_old); 

        printf("#%%$: CP_Restore_Time_%d: %lf\n", cp_file_num, (double)(tclose - topen) / CLOCKS_PER_SEC); 

        decompress_time_old = decompress_time;
        rd_time_old = rd_time;
        restore_time += (double)(tclose - topen) / CLOCKS_PER_SEC;
    }

    buffer_free();
}

void cpc_profile_(){
    printf("#%%$: CP_Com_Time_All: %lf\n", compress_time);
    printf("#%%$: CP_Wr_Time_All: %lf\n", wr_time); 
    printf("#%%$: CP_Store_Time_All: %lf\n", store_time); 
    printf("#%%$: CP_Decom_Time_All: %lf\n", decompress_time);
    printf("#%%$: CP_Rd_Time_All: %lf\n", rd_time); 
    printf("#%%$: CP_Restore_Time_All: %lf\n", restore_time); 
}

include(`foreach.m4')dnl
include(`forloop.m4')dnl

define(`CONCAT', `$1$2')dnl

define(`ZFP',dnl
`dnl

void compresswr_$1(void *data, size_t size, int dim, int *shape){
    zfp_stream *zfp;
    zfp_type type = zfp_type_$2;                          
    zfp_field *field;
    bitstream *stream;
    size_t bufsize;  
    size_t zfpsize;
    clock_t t1, t2, t3;

    t1 = clock();

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
    zfp_stream_set_accuracy(zfp, $3);                     

    // allocate buffer for compressed data
    bufsize = zfp_stream_maximum_size(zfp, field);    
    buffer_resize((size_t)bufsize);              

    // associate bit stream with allocated buffer
    stream = stream_open((void*)buffer, bufsize);      
    zfp_stream_set_bit_stream(zfp, stream);                  
    zfp_stream_rewind(zfp);                  

    // compress array
    zfpsize = zfp_compress(zfp, field);               

    t2 = clock();

    write(fd, &zfpsize, sizeof(zfpsize));
    write(fd, (void*)buffer, zfpsize);

    t3 = clock();

    compress_time += (double)(t2 - t1) / CLOCKS_PER_SEC;
    wr_time += (double)(t3 - t2) / CLOCKS_PER_SEC;

    //printf("Write %zu -> %zu\n", size, zfpsize);
}

void compressrd_$1(void *data, size_t size, int dim, int *shape){
    zfp_stream *zfp;
    int ret;
    zfp_type type = zfp_type_$2;                          
    zfp_field *field;
    bitstream *stream;
    size_t bufsize;  
    size_t zfpsize;
    clock_t t1, t2, t3;
    
    t1 = clock();

    // allocate buffer for compressed data                     
    read(fd, &bufsize, sizeof(bufsize));
    buffer_resize((size_t)bufsize);   
    read(fd, (void*)buffer, bufsize);

    t2 = clock();

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
    zfp_stream_set_accuracy(zfp, $3);                      


    // associate bit stream with allocated buffer
    stream = stream_open((void*)buffer, bufsize);      
    zfp_stream_set_bit_stream(zfp, stream);                  
    zfp_stream_rewind(zfp);                  

    ret = zfp_decompress(zfp, field);

    t3 = clock();

    decompress_time += (double)(t3 - t2) / CLOCKS_PER_SEC;
    rd_time += (double)(t2 - t1) / CLOCKS_PER_SEC;
    
    if (ret < 0) {
        printf("Decompress fail: addr: %llx, size: %zu\n", data, size);   
    }
    else {
        //printf("Read %zu -> %zu\n", bufsize, size);
    }
}

')dnl
dnl

define(`CMP_REAL',changequote(`[', `]')dnl
[dnl

void compresswr_$1_($2 *R, int *size, int *dim, int *shape ) {
    if (*size > THRESHOLD){
        compresswr_$3((void*)R, (size_t)(*size), *dim, shape);
    }
    else{
        clock_t t1, t2;
        t1 = clock();
        write(fd, R, (size_t)(*size));
        t2 = clock();
        wr_time += (double)(t2 - t1) / CLOCKS_PER_SEC;
    }
}

void compressrd_$1_($2 *D, int *size, int *dim, int *shape  ) {
    if (*size > THRESHOLD){
        compressrd_$3((void*)D, (size_t)(*size), *dim, shape);
    }
    else{
        clock_t t1, t2;
        t1 = clock();
        read(fd, D, (size_t)(*size));
        t2 = clock();
        rd_time += (double)(t2 - t1) / CLOCKS_PER_SEC;
    }
}

]changequote([`], [']))dnl
dnl

foreach(`_arg', (`(`real', `double', `0.000001')',dnl
                `(`int', `int32', `0')'), `ZFP(translit(_arg, `()'))')

foreach(`i', (`(`real', `double', `real')',dnl
                `(`integer', `int', `int')',dnl
                `(`bool', `int', `int')'), `CMP_REAL(translit(i, `()'))')


void compresswr_string_(char *R, int *size , long l ) {
        clock_t t1, t2;
        t1 = clock();
        write(fd, R, (size_t)(*size));
        t2 = clock();
        wr_time += (double)(t2 - t1) / CLOCKS_PER_SEC;
}

void compressrd_string_(char *D, int *size , long l ) {
    clock_t t1, t2;
    t1 = clock();
    read(fd, D, (size_t)(*size));
    t2 = clock();
    rd_time += (double)(t2 - t1) / CLOCKS_PER_SEC;
}
