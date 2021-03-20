#include <mpi.h>
#include <iostream>
#include <complex>
#include <assert.h>
#include "omp.h"
#include <cmath>
#include "time.h"
#include "sys/time.h"
#include <stdio.h>
#include <stdlib.h>

using namespace std;
typedef complex<double> complexd;









complexd *read(char *f, int rank, unsigned long long seg_size) {
    MPI_File file;
    MPI_File_open(MPI_COMM_WORLD, f, MPI_MODE_RDONLY, MPI_INFO_NULL,&file);

    //auto *A = new complexd[seg_size];
    complexd *A;
    A = (complexd*) malloc(sizeof(complexd) * seg_size);



    double d[2];
    MPI_File_seek(file, 2 * seg_size * rank * sizeof(double), MPI_SEEK_SET);
    for (int i = 0; i < seg_size; ++i) {
        MPI_File_read(file, &d, 2, MPI_DOUBLE, MPI_STATUS_IGNORE);
        A[i].real(d[0]);
        A[i].imag(d[1]);
    }
    MPI_File_close(&file);
    return A;
}








void write(char *f, complexd *B, int n, int rank, int size, unsigned long long seg_size) {
    MPI_File file;
    MPI_File_open(MPI_COMM_WORLD, f, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &file);
    double d[2];
    MPI_File_seek(file, 2 * seg_size * rank * sizeof(double), MPI_SEEK_SET);
    for (int i = 0; i < seg_size; ++i) {
        d[0] = B[i].real();
        d[1] = B[i].imag();
        MPI_File_write(file, &d, 2, MPI_DOUBLE, MPI_STATUS_IGNORE);
    }
    MPI_File_close(&file);
}










complexd* generate_condition(unsigned long long seg_size, int rank, int size){ 
    double module = 0;
    unsigned int seed = time(NULL) + rank;
    complexd *V = new complexd[seg_size];
    for (long long unsigned  i = 0; i < seg_size; i++){
        V[i].real(rand_r(&seed)%100 + 1);
        V[i].imag(rand_r(&seed)%100 + 1);
        module += abs(V[i] * V[i]);
    }
    
    int rc;
    double new_m;
    MPI_Status stat;
    if(rank != 0){
        module += 1;
        rc = MPI_Send(&module, 1, MPI_DOUBLE, 0, 999, MPI_COMM_WORLD);
        MPI_Recv(&module, 1, MPI_DOUBLE, MPI_ANY_SOURCE, 3, MPI_COMM_WORLD, &stat);
    }
    else{
        for(int i = 1; i < size; i++){
            MPI_Recv(&new_m, 1, MPI_DOUBLE, MPI_ANY_SOURCE, 999, MPI_COMM_WORLD, &stat);
            module += new_m;
        }
        module = sqrt(module);
        for(int i = 1; i < size; i++){
            rc = MPI_Send(&module, 1, MPI_DOUBLE, i, 3, MPI_COMM_WORLD);
        }
    }
    for (long long unsigned j = 0; j < seg_size; j++) {
        V[j] /= module;
    }
    return V;
}









void OneQubitEvolution(complexd *in,complexd *out,complexd U[2][2],int n, int q, int rank, unsigned long long seg_size){
    int shift = n - q;
    int pow = 1 << (shift);
    int N = 1 << n;
    for (int i = rank * seg_size; i < (rank+1)*seg_size; i++) {
        int i0 = i & ~pow;
        int i1 = i | pow;
        int iq = (i & pow) >> shift;
        out[i-rank * seg_size] = U[iq][0] * in[i0] + U[iq][1] * in[i1];
    }
}










int main(int argc, char **argv) {
    int was_read = 0;
    int test = 0;
    char *input, *output, *test_file;
    unsigned k, n;
    complexd *V;
    complexd *need;
    for (int i = 1; i < argc; i++) { 
        string option(argv[i]);
        if ((option.compare("k") == 0)) {
            k = atoi(argv[++i]);
        }
        if (option.compare("n") == 0) {
            n = atoi(argv[++i]);
        }
        if ((option.compare("file_read") == 0)) {
            input = argv[++i];
            was_read = 1;
        }
        if ((option.compare("file_write") == 0)) {
            output = argv[++i];
        }
        if ((option.compare("test") == 0)) {
            test = 1;
        }
        if ((option.compare("file_test") == 0)) {
            test_file = argv[++i];
        }
    }
    MPI_Init(&argc, &argv);
    int rank;
    int size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    unsigned long long index = 1LLU << n;
    unsigned long long seg_size = index / size;
    need = (complexd*) malloc(sizeof(complexd) * index);
    if (was_read == 0) {
        V = generate_condition(seg_size, rank, size); 
    } else {
        V = read(input, rank, seg_size);
    }
    struct timeval start, stop;
    complexd U[2][2];
    U[0][0] = 1 / sqrt(2);
    U[0][1] = 1 / sqrt(2);
    U[1][0] = 1 / sqrt(2);
    U[1][1] = -1 / sqrt(2);
    MPI_Status stat1;
    int error = 0;
    int rc1;
    if(rank != 0){
        rc1 = MPI_Send(V, seg_size, MPI_DOUBLE_COMPLEX, 0, 4, MPI_COMM_WORLD);
        MPI_Recv(need, index, MPI_DOUBLE_COMPLEX, 0, 5, MPI_COMM_WORLD, &stat1);
    }
    else{
        for(int i = 0; i < seg_size; i++){
            need[i] = V[i];
        }
        for(int i = 1; i < size; i++){
            MPI_Recv(V, seg_size, MPI_DOUBLE_COMPLEX, i, 4, MPI_COMM_WORLD, &stat1);
            for(int j = 0; j < seg_size; j++){
                need[j + seg_size*i] = V[j];
            }
        }
        for(int i = 1; i < size; i++){
            rc1 = MPI_Send(need, index, MPI_DOUBLE_COMPLEX, i, 5, MPI_COMM_WORLD);
        }
    }
    double begin = MPI_Wtime();
    OneQubitEvolution(need, V, U, n, k, rank, seg_size);
    double end = MPI_Wtime();
    if (test == 1) {
        int rc;
        complexd *test_vector = read(test_file, rank, seg_size);
        for (int i = 0; i < seg_size; i++) {
            if (abs(test_vector[i].real() - V[i].real()) > 0.000001 || abs(test_vector[i].imag() - V[i].imag()) > 0.000001) {
                error = 1;
            }
        }
        if(rank > 0){
            rc = MPI_Send(&error, 1, MPI_INT, 0, 1, MPI_COMM_WORLD);
        }
        else{
            MPI_Status stat;
            for(int i= 1; i < size; i++){
                MPI_Recv(&rc, 1, MPI_INT, MPI_ANY_SOURCE, 1, MPI_COMM_WORLD, &stat);
                error += rc;
            }
            if(error != 0){
                cout << "error had been done" << endl;
            }
            else{
                cout << "no errors" << endl;
            }
        }
    } else {
        write(output, V, n, rank, size, seg_size);
    }
    cout << end - begin << endl;
    MPI_Finalize();
    free(V);
    free(need);
}