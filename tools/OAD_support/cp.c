#define _GNU_SOURCE

#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <linux/limits.h>
#include <fcntl.h>
#ifdef ALLOW_USE_MPI
#include <mpi.h>
#endif
#include <string.h>
#include <time.h>

#define BSIZE (1 * 1024 * 1024)
#define FSIZE (100 * 1024 * 1024)

typedef struct cp_fd{
    int fd;
    char *buf, *abuf, *cbuf;
};

static int cp_file_num = 0;
int fd;
int wr;

size_t bsize, bused;
void *buffer, *abuffer;

static double compress_time, compress_time_old;
static double decompress_time, decompress_time_old;
static double wr_time, wr_time_old;
static double rd_time, rd_time_old;
static double store_time, restore_time;

clock_t topen;

void *iobuf;
size_t iobsize;

void buffer_init(){
    int pagesize;

    pagesize = getpagesize();
    buffer = malloc(BSIZE + pagesize);
    abuffer=(void*)((((unsigned long long)realbuff+(unsigned long long)pagesize-1)/(unsigned long long)pagesize)*(unsigned long long)pagesize);

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

cp_fd cp_wr_open(char* fname, size_t fsize){
    int pagesize;
    cp_fd tmp;

    pagesize = getpagesize();
    tmp->buf  = (char*)malloc(BSIZE + pagesize);
    tmp->abuf = tmp->cbuf = (char*)((((size_t)tmp->buf + (size_t)pagesize - 1) / (size_t)pagesize) * (size_t)pagesize);

    tmp->fd = open(fname, O_CREAT | O_WRONLY | O_TRUNC | O_DIRECT | O_SYNC, 0644);

    return tmp;
}

int cp_write(cp_fd fd, void *data, size_t size){
    memcpy(fd->cbuf, data, size);
    fd->cbuf += size;
}

int cp_wr_close(cp_fd fd){
    int ret = 0;
    off_t wsize;
    ssize_t ioret;

    wsize = 0;
    while(fd->abuf < fd->cbuf){
        ioret = write(fd->fd, fd->abuf, fd->cbuf - fd->abuf);
        if (ioret <= 0){
            ret = -1;
        }
        fd->abuf += ioret;
    }
    
    close(tmp->fd);
    
    free(tmp->buf);

    return ret;
}

cp_fd cp_rd_open(char* fname, int flag){
    int pagesize;
    cp_fd tmp;
    struct stat st;

    stat(fname, &st);
    printf("#%%$: CP_Size_%d: %lld\n", cp_file_num, (long long)st.st_size);

    pagesize = getpagesize();
    tmp->buf  = (char*)malloc(BSIZE + pagesize);
    tmp->abuf = tmp->cbuf = (char*)((((size_t)tmp->buf + (size_t)pagesize - 1) / (size_t)pagesize) * (size_t)pagesize);

    tmp->fd = open(fname, O_RDONLY | O_TRUNC | O_DIRECT | O_SYNC, 0644);

    return tmp;
}

int cp_read(cp_fd fd, void *data, size_t size){
    memcpy(data, fd->cbuf, size);
    fd->cbuf += size;
}

int cp_rd_close(cp_fd fd){   
    close(tmp->fd);
    
    free(tmp->buf);

    return 0;
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

    fd = cp_wr_open(fname, FSIZE);
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

    //buffer_init();
    
    sprintf(fname, "oad_cp.%03d.%05d", rank, cp_file_num);

    wr = 0;

    topen = clock();

    fd = cp_rd_open(fname);
}

void cpc_close_(){
    int rank;
    char fname[PATH_MAX];
    struct stat st;

    //buffer_free();
    
    if (wr){
        cp_wr_close(fd);

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
        cp_rd_close(fd);

        printf("#%%$: CP_Decom_Time_%d: %lf\n", cp_file_num, decompress_time - decompress_time_old);
        printf("#%%$: CP_Rd_Time_%d: %lf\n", cp_file_num, rd_time - rd_time_old); 

        printf("#%%$: CP_Restore_Time_%d: %lf\n", cp_file_num, (double)(tclose - topen) / CLOCKS_PER_SEC); 

        decompress_time_old = decompress_time;
        rd_time_old = rd_time;
        restore_time += (double)(tclose - topen) / CLOCKS_PER_SEC;
    }
}

void cpc_profile_(){
    printf("#%%$: CP_Com_Time_All: %lf\n", compress_time);
    printf("#%%$: CP_Wr_Time_All: %lf\n", wr_time); 
    printf("#%%$: CP_Store_Time_All: %lf\n", store_time); 
    printf("#%%$: CP_Decom_Time_All: %lf\n", decompress_time);
    printf("#%%$: CP_Rd_Time_All: %lf\n", rd_time); 
    printf("#%%$: CP_Restore_Time_All: %lf\n", restore_time); 
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

