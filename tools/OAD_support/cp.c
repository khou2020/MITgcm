


#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <linux/limits.h>
#include <fcntl.h>
#ifdef ALLOW_USE_MPI
#include <mpi.h>
#endif
#include <time.h>
#define BSIZE 1048576
     
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

    buffer_free();
    
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
}

void cpc_profile_(){
    printf("#%%$: CP_Com_Time_All: %lf\n", compress_time_all);
    printf("#%%$: CP_Wr_Time_All: %lf\n", wr_time_all); 
    printf("#%%$: CP_Store_Time_All: %lf\n", store_time_all); 
    printf("#%%$: CP_Decom_Time_All: %lf\n", decompress_time_all);
    printf("#%%$: CP_Rd_Time_All: %lf\n", rd_time_all); 
    printf("#%%$: CP_Restore_Time_All: %lf\n", restore_time_all); 
}

inline void float2double(double *dst, float *src, int n){
    int i;
    for(i = 0; i < n; i++){
        dst[i] = (double)src[i];
    }
}

inline void double2float(float *dst, double *src, int n){
    int i;
    for(i = 0; i < n; i++){
        dst[i] = (float)src[i];
    }
}

void compresswr_real_(double *R, int* size  ) {
    int i, n;
    double t1, t2, t3;

    //printf("Write %d bytes from %llx\n", *size, R);
    n = (*size) >> 3;
    t1 = getwalltime();
    buffer_resize(n * sizeof(float));
    double2float((float*)buffer, R, n);
    t2 = getwalltime();
    write(fd, buffer, n * sizeof(float));
    t3 = getwalltime();

    compress_time += t2 - t1;
    wr_time += t3 - t2;
}

void compressrd_real_(double *D, int *size  ) {
    int i, n;
    double t1, t2, t3;

    //printf("Read %d bytes to %llx\n", *size, D);
    n = (*size) >> 3;
    t1 = getwalltime();
    // Buffersize is already guaranteed on writing
    read(fd, buffer, n * sizeof(float));
    t2 = getwalltime();
    float2double(D, (float*)buffer, n);
    t3 = getwalltime();

    rd_time += t2 - t1;
    decompress_time += t3 - t2;
}




void compresswr_integer_(int *R, int* size  ) {
    double t1, t2;
    //printf("Write %d bytes from %llx\n", *size, R);
    t1 = getwalltime();
    write(fd, R, (size_t)(*size));
    t2 = getwalltime();
    wr_time += t2 - t1;
}

void compressrd_integer_(int *D, int *size  ) {
    double t1, t2;
    //printf("Read %d bytes to %llx\n", *size, D);
    t1 = getwalltime();
    read(fd, D, (size_t)(*size));
    t2 = getwalltime();
    rd_time += t2 - t1;
}


void compresswr_bool_(int *R, int* size  ) {
    double t1, t2;
    //printf("Write %d bytes from %llx\n", *size, R);
    t1 = getwalltime();
    write(fd, R, (size_t)(*size));
    t2 = getwalltime();
    wr_time += t2 - t1;
}

void compressrd_bool_(int *D, int *size  ) {
    double t1, t2;
    //printf("Read %d bytes to %llx\n", *size, D);
    t1 = getwalltime();
    read(fd, D, (size_t)(*size));
    t2 = getwalltime();
    rd_time += t2 - t1;
}


void compresswr_string_(char *R, int* size , long l ) {
    double t1, t2;
    //printf("Write %d bytes from %llx\n", *size, R);
    t1 = getwalltime();
    write(fd, R, (size_t)(*size));
    t2 = getwalltime();
    wr_time += t2 - t1;
}

void compressrd_string_(char *D, int *size , long l ) {
    double t1, t2;
    //printf("Read %d bytes to %llx\n", *size, D);
    t1 = getwalltime();
    read(fd, D, (size_t)(*size));
    t2 = getwalltime();
    rd_time += t2 - t1;
}

