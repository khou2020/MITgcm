

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
int fd;
int wr;

size_t bsize;
char *buffer;

static double compress_time, compress_time_old;
static double decompress_time, decompress_time_old;
static double wr_time, wr_time_old;
static double rd_time, rd_time_old;

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

    //buffer_init();

    wr = 1;
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

    //buffer_init();

    wr = 0;
}

void cpc_close_(){
    int rank;
    char fname[PATH_MAX];
    struct stat st;

    close(fd);

    //buffer_free();
    
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

        compress_time_old = compress_time;
        wr_time_old = wr_time;
    }
    else{
        printf("#%%$: CP_Decom_Time_%d: %lf\n", cp_file_num, decompress_time - decompress_time_old);
        printf("#%%$: CP_Rd_Time_%d: %lf\n", cp_file_num, rd_time - rd_time_old); 

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




void compresswr_real_(double *R, int* size  ) {
    //printf("Write %d bytes from %llx\n", *size, R);
    clock_t t1, t2;
    t1 = clock();
    write(fd, R, (size_t)(*size));
    t2 = clock();
    wr_time += (double)(t2 - t1) / CLOCKS_PER_SEC;
}

void compressrd_real_(double *D, int *size  ) {
    //printf("Read %d bytes to %llx\n", *size, D);
    clock_t t1, t2;
    t1 = clock();
    read(fd, D, (size_t)(*size));
    t2 = clock();
    rd_time += (double)(t2 - t1) / CLOCKS_PER_SEC;
}


void compresswr_integer_(int *R, int* size  ) {
    //printf("Write %d bytes from %llx\n", *size, R);
    clock_t t1, t2;
    t1 = clock();
    write(fd, R, (size_t)(*size));
    t2 = clock();
    wr_time += (double)(t2 - t1) / CLOCKS_PER_SEC;
}

void compressrd_integer_(int *D, int *size  ) {
    //printf("Read %d bytes to %llx\n", *size, D);
    clock_t t1, t2;
    t1 = clock();
    read(fd, D, (size_t)(*size));
    t2 = clock();
    rd_time += (double)(t2 - t1) / CLOCKS_PER_SEC;
}


void compresswr_bool_(int *R, int* size  ) {
    //printf("Write %d bytes from %llx\n", *size, R);
    clock_t t1, t2;
    t1 = clock();
    write(fd, R, (size_t)(*size));
    t2 = clock();
    wr_time += (double)(t2 - t1) / CLOCKS_PER_SEC;
}

void compressrd_bool_(int *D, int *size  ) {
    //printf("Read %d bytes to %llx\n", *size, D);
    clock_t t1, t2;
    t1 = clock();
    read(fd, D, (size_t)(*size));
    t2 = clock();
    rd_time += (double)(t2 - t1) / CLOCKS_PER_SEC;
}


void compresswr_string_(char *R, int* size , long l ) {
    //printf("Write %d bytes from %llx\n", *size, R);
    clock_t t1, t2;
    t1 = clock();
    write(fd, R, (size_t)(*size));
    t2 = clock();
    wr_time += (double)(t2 - t1) / CLOCKS_PER_SEC;
}

void compressrd_string_(char *D, int *size , long l ) {
    //printf("Read %d bytes to %llx\n", *size, D);
    clock_t t1, t2;
    t1 = clock();
    read(fd, D, (size_t)(*size));
    t2 = clock();
    rd_time += (double)(t2 - t1) / CLOCKS_PER_SEC;
}

