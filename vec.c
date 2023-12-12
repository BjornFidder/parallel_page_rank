
#include "vec.h"

// long rand_r(unsigned int *seed) {
//     *seed = *seed * 1103515245 + 12345;
//     return (long)((unsigned int)(*seed/(RAND_MAX)) % (RAND_MAX/2+1));
// }   

long rand_long(long N, unsigned int* seed) 
{
    return rand_r(seed) % N;
}

void print_vec(long N, double* x) 
{
    for (long i = 0; i < N; i++) 
        printf("%f ", x[i]);
    printf("\n");
}

void gen_rand_prob_vec(long N, double* u, unsigned int* seed) 
{
    //Generate a random stochastic vector by generating an array
    //of random integers, and then normalize by dividing by the sum

    double sum = 0;
    for (long i = 0; i<N; i++)
    {
        u[i] = (double)rand_long(RAND_MAX/(N*10), seed);
        sum += u[i];
    }
    for (long i = 0; i<N; i++) 
        u[i] /= sum;
}

double norm_vec(long N, double* x) 
{
    //compute the squared norm of a vector x of size N
    double norm = 0;
    for (long i = 0; i < N; i++)
        norm += x[i]*x[i];
    return norm;
}

void add_vec(long N, double* x, double* y) 
{
    //add vector x to vector y
    //giving the output in vector x
    for (long i = 0; i < N; i++)
        x[i] = x[i]+y[i];
}

double* div_vec(long N, double* x, uint8_t* y) 
{
    //divide vector x by vector y
    //giving the output in vector z
    double* z = (double*)malloc(N*sizeof(double));
    for (long i = 0; i < N; i++)
        z[i] = x[i]/y[i];
    return z;
}

void initi(long N, long* x, long a) 
{
    for (long i = 0; i < N; i++) 
    {
        x[i] = a;
    }
}

void initi8(long N, uint8_t* x, uint8_t a) 
{
    for (long i = 0; i < N; i++) 
    {
        x[i] = a;
    }
}

void initb(long N, bool* x, bool a) 
{
    for (long i = 0; i < N; i++) 
    {
        x[i] = a;
    }
}

void initd(long N, double* x, double a) 
{
    for (long i = 0; i < N; i++) 
        {
            x[i] = a;
        }
}