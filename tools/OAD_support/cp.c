

#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <linux/limits.h>
#include <fcntl.h>
#ifdef ALLOW_USE_MPI
#include <mpi.h>
#endif

static int cp_file_num = 0;
int fd;

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
}

void cpc_close_(){
    close(fd);
}





void compresswr_real_(double *R, size_t*size  ) {
    //printf("Write %d bytes from %llx\n", *size, R);
    write(fd, R, *size);
}

void compressrd_real_(double *D, size_t *size  ) {
    //printf("Read %d bytes to %llx\n", *size, D);
    read(fd, D, *size);
}


void compresswr_integer_(int *R, size_t*size  ) {
    //printf("Write %d bytes from %llx\n", *size, R);
    write(fd, R, *size);
}

void compressrd_integer_(int *D, size_t *size  ) {
    //printf("Read %d bytes to %llx\n", *size, D);
    read(fd, D, *size);
}


void compresswr_bool_(int *R, size_t*size  ) {
    //printf("Write %d bytes from %llx\n", *size, R);
    write(fd, R, *size);
}

void compressrd_bool_(int *D, size_t *size  ) {
    //printf("Read %d bytes to %llx\n", *size, D);
    read(fd, D, *size);
}


void compresswr_string_(char *R, size_t*size , long l ) {
    //printf("Write %d bytes from %llx\n", *size, R);
    write(fd, R, *size);
}

void compressrd_string_(char *D, size_t *size , long l ) {
    //printf("Read %d bytes to %llx\n", *size, D);
    read(fd, D, *size);
}

