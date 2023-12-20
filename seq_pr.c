#include "gen_graph.h"
#include <stdbool.h>
#include <time.h>

double q = 0.85;

double* mul_mat_vec(long N, long* cols, long* start, double* x) 
{
    //matrix vector multiplication, 
    //where matrix is 0-1 valued, given in CRS format
    //output given in double* y
    double* y = vecallocd(N);
    initd(N, y, 0);
    for (long i = 0; i < N; i++) 
        for (long k = start[i]; k < start[i+1]; k++)
            y[i] += x[cols[k]];
    return y;
}

double* comp_residual(long N, long* cols, long* start, uint8_t* D, double* u) 
{
    //Compute the residual according to r0 = e - (I - qGD^{-1})u0
    double* r = vecallocd(N);

    double* Du = div_vec(N, u, D);
    double* GDu = mul_mat_vec(N, cols, start, Du);
    vecfreed(Du);
    for (long i = 0; i < N; i++)
        r[i] = 1 - u[i] + q*GDu[i];
    vecfreed(GDu);
    return r;
}

void iterate(long N, long* cols, long* start, uint8_t* D, double* u, double* r)
{
    //Perform one iteration of the PageRank algorithm
    add_vec(N, u, r);
    double* Dr = div_vec(N, r, D);
    double* GDr = mul_mat_vec(N, cols, start, Dr);
    vecfreed(Dr);
    for (long i = 0; i < N; i++) r[i] = q*GDr[i];
    vecfreed(GDr);
}

void print_vecs(long N, double* u, double* r)
{
    //Print vectors u and r of size N, and their L1/L2 norm
    //we only print the components if N<=10.

    if (N<=10) {
    printf("u: ");
    print_vec(N, u);
    }
    double sum = 0; 
    for (long i = 0; i < N; i++) sum += u[i];
    printf("sum: %f\n", sum);

    if (N<=10) {
    printf("u normalized: ");
    for (long i = 0; i < N; i++) printf("%f ", u[i]/sum);
    printf("\n");}

    bool prob_vec = true;
    for (long i = 0; i < N; i++) if (u[i] < 0) prob_vec = false;
    if (!prob_vec) printf("Not a probability vector!\n");

    printf("r: ");
    if (N<=10) print_vec(N, r);

    printf("norm: %f\n", sqrt(norm_vec(N,r)));

}

long solve_pr(long N, long* cols, long* start, uint8_t* D, unsigned int* seed)
{
    //solve the PageRank problem (I - pGD^{-1})u = e iteratively
    //N is the system size
    //G is given in CRS format by start and cols
    //D is given as an array containing the diagonal elements

    //initial vectors
    double* u = vecallocd(N);
    gen_rand_prob_vec(N, u, seed);
    
    double* r = comp_residual(N, cols, start, D, u);
    if (N <= 10) print_vecs(N, u, r);

    long count = 0;
    while (norm_vec(N, r) >= pow(10, -12))
    { 
        iterate(N, cols, start, D, u, r);
        count++;
    }
    if (N <= 10) print_vecs(N, u, r);
    vecfreed(u);
    vecfreed(r);

    return count;
}

int main(int argc, char **argv) 
{
    long N;
    printf("Please enter N:\n");
    if (scanf("%ld", &N) != 1) 
        printf("Please input a number\n");
    printf("Running PageRank algorithm with N = %ld\n", N);

    // srand(31415);
    
    clock_t clock0, clock1, clock2, clock3;

    clock0 = clock();

    /////////////////////////
    unsigned int seed = (unsigned int)(RAND_MAX*clock0);
    long* start;

    // long* cols = gen_graph(N, N, start, &seed);
    // clock1 = clock();
    // uint8_t* D  = outlinks(N, N, cols, start);
    
    clock2 = clock();

    if (N <= 10) 
        {print_graph(N, cols, start);
        print_D(N, D);}

    long count = solve_pr(N, cols, start, D, &seed);
    
    vecfreei(cols);
    vecfreei(start);
    vecfreei8(D);

    clock3 = clock();
    double gen_time = (double)(clock2 - clock0) / CLOCKS_PER_SEC;
    double outlink_time = (double)(clock2 - clock1) / CLOCKS_PER_SEC;
    double solve_time = (double)(clock3 - clock2) / CLOCKS_PER_SEC;
    double total_time = (double)(clock3 - clock0) / CLOCKS_PER_SEC;

    printf("Needed %ld iterations\n", count);
    printf("Generation run-time: %f\n", gen_time);
    printf("Finding outlinks: %f\n", outlink_time);
    printf("Solving run-time: %f\n", solve_time);
    printf("Time per iteration: %f\n", solve_time / count);
    printf("Total run-time: %f\n", total_time);
    printf("\n");
}

//test graph
// start = vecalloci(N+1); 
// start[0] = 0;
// start[1] = 2;
// start[2] = 3;
// start[3] = 4;
// start[4] = 6;
// start[5] = 7;
// start[6] = 9;

// long* cols = vecalloci(10);
// cols[0] = 3;
// cols[1] = 5;
// cols[2] = 0;
// cols[3] = 1;
// cols[4] = 1;
// cols[5] = 2;
// cols[6] = 2;
// cols[7] = 0;
// cols[8] = 2;

// uint8_t* D = vecalloci8(N);
// D[0] = 2;
// D[1] = 2;
// D[2] = 3;
// D[3] = 1;
// D[4] = 1;
// D[5] = 1;
