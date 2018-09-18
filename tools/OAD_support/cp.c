




#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <linux/limits.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#ifdef ALLOW_USE_MPI
#include <mpi.h>
#endif
#include <time.h>

#define MAXITR 1024
#define BSIZE 1048576

static int cp_file_num = 0;
int cur_num;
int fd;
int wr;

size_t bsize;
char *buffer;

double compress_time, decompress_time, wr_time, rd_time;

double times[MAXITR][10];
unsigned long long fsize[MAXITR];
static int max_itr;

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

    MPI_Barrier(MPI_COMM_WORLD);

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

    MPI_Barrier(MPI_COMM_WORLD);
    
    topen = getwalltime();

    fd = open(fname, O_CREAT | O_RDONLY, 0644);
}

void cpc_close_(){
    int rank, np;
    char fname[PATH_MAX];
    struct stat st;
    double tclose, tio;
    unsigned long long size;

    buffer_free();
    
    if (wr){
        fsync(fd);
        close(fd);
        
        tclose = getwalltime();
        tio = tclose - topen;

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

        MPI_Reduce(&size, &(fsize[cur_num]), 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(&tio, &(times[cur_num][0]), 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
        MPI_Reduce(&compress_time, &(times[cur_num][1]), 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(&wr_time, &(times[cur_num][2]), 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        times[cur_num][1] /= np;
        times[cur_num][2] /= np;
    }
    else{
        close(fd);

        tclose = getwalltime();
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

