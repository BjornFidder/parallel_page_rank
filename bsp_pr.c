#include "gen_graph.h"

long P;
double q = 0.85;

double vec_total_norm(long n, double* x, double* norms)
{
    //compute total norm of u in parallel
    long p = bsp_nprocs();
    long s = bsp_pid();

    double norm = norm_vec(n, x);
    
    for (long t = 0; t < p; t++)
        bsp_put(t, &norm, norms, s*sizeof(double), sizeof(double));
    bsp_sync();
    	
    norm = 0;
    for (long t = 0; t < p; t++)
        norm += norms[t];

    return norm;
}

double gen_rand_vec(long n, double* u, unsigned int* seed) 
{
    //Generate n random vector components and return the sum
    double sum = 0;
    for (long i = 0; i<n; i++)
    {
        u[i] = (double)rand_long(RAND_MAX/(n*10), seed);
        sum += u[i];
    }
    return sum;
}

void div_vec_pr(long n, double* x, uint8_t* y, double* z) 
{
    //divide vector x by vector y
    //giving the output in vector z
    for (long i = 0; i < n; i++)
        z[i] = x[i]/y[i];
}


void print_graphs(long n, uint8_t* D, long* start, long* cols)
{
    long p = bsp_nprocs();
    long s = bsp_pid();

    for (long i = 0; i < p; i++) 
    {
        if (s == i) 
        {
            printf("processor %ld: \n", s);
            print_graph(n, cols, start);
            print_D(n, D); printf("\n\n");
        }
        bsp_sync();
    }
    
}
void print_vecs_pr(long n, double* u, double* r, bool print) 
{
    long p = bsp_nprocs();
    long s = bsp_pid();

    for (long i = 0; i < p; i++) 
    {
        if (s == i) 
        {
            printf("processor %ld: \n", s);
            if (print)
                {printf("u: "); 
                print_vec(n, u);}
            printf("r: ");
            if (print) print_vec(n, r);
            printf("local squared norm: %f\n", norm_vec(n, r));
        }
        bsp_sync();
    } 
}

uint8_t* outlinks_pr(long N, long n, long* cols, long* start) 
{
    long p = bsp_nprocs();
    long s = bsp_pid();

    //compute local outlinks
    uint8_t* ds = vecalloci8(n*p);
    initi8(n*p,ds,0);
    bsp_push_reg(ds, n*p*sizeof(uint8_t));
    bsp_sync();
    
    uint8_t* d = outlinks(N, n, cols, start); //local outlinks

    // 'put' local outlinks in local ds
    for (long j = 0; j < N; j++)
        ds[j/p+s*n] = d[j];

    //broadcast local outlinks
    long t;
    for (long j = 0; j < N; j++)
    {
        t = j%p;
        if (t != s && d[j] > 0) 
            bsp_put(t, &d[j], ds, (j/p + s*nloc(p, t, N))*sizeof(uint8_t), sizeof(uint8_t));
    }
    free(d);
    bsp_sync();
    
    //compute total outlinks for local rows
    uint8_t* D = vecalloci8(n);
    initi8(n, D, 0);
    for (long i = 0; i < n; i++) 
        for (long t = 0; t < p; t++)
            D[i] += ds[i+n*t];

    outlinks_noZeroes(n, D);
    bsp_pop_reg(ds);
    free(ds);

    return D;
}

double* initial_vector(long N, long n, unsigned int* seed) 
{
    long p = bsp_nprocs();
    long s = bsp_pid();
    
    //  Generate random initial vector u0 //
    //register
    double* sums = vecallocd(p);
    bsp_push_reg(sums, p*sizeof(double));
    bsp_sync();

    //generate and distribute sum to all processors
    double* u = vecallocd(n);
    double sum = gen_rand_vec(n, u, seed);
    for (long t = 0; t < p; t++)
        bsp_put(t, &sum, sums, s*sizeof(double), sizeof(double));
    bsp_sync();

    //compute total sum and normalize local components
    sum = 0;
    for (long t = 0; t < p; t++)
        sum += sums[t];
    bsp_pop_reg(sums);
    vecfreed(sums);
    for (long i = 0; i < n; i++)
        u[i] /= sum;

    return u;
}


void mul_GD(long N, long n, double* u, uint8_t* D, long* start, long* cols, bool* cols_get, double* Du, double* vals, double* v) 
{
    //compute v = G*Dinv*u in parallel
    long p = bsp_nprocs();
    long s = bsp_pid();

    initd(n, v, 0);
    div_vec_pr(n, u, D, Du);

    bsp_sync();

    //Slow method:
    // for (int i = 0; i < n; i++)
    //     for (int k = start[i]; k < start[i+1]; k++)
    //             bsp_get(cols[k] % p, Du, cols[k]/p * sizeof(double), &vals[k], sizeof(double));
    // bsp_sync();
    // for (int i = 0; i < n; i++)
    //     for (int k = start[i]; k < start[i+1]; k++)
    //         v[i] += vals[k];

    // 'Get' values from this processor
    for (int i = 0; i < n; i++)
        vals[s*p+i] = Du[i];

    // //Determine which components of v to get from other processors
    for (int j = 0; j < N; j++)
        if (cols_get[j]) bsp_get(j % p, Du, j/p * sizeof(double), &vals[j], sizeof(double));

    bsp_sync();

    for (int i = 0; i < n; i++)
        for (int k = start[i]; k < start[i+1]; k++)
            {v[i] += vals[cols[k]];}

}

void bsp_pr() 
{
    bsp_begin(P);
    long p = bsp_nprocs();
    long s = bsp_pid();
    long N;

    bsp_push_reg(&N, sizeof(long));

    if (s==0) {

        printf("Please enter N:\n");
        if (scanf("%ld", &N) != 1) 
            printf("Please input one number\n");
        printf("Running PageRank algorithm for N = %ld\n", N);
    }
    bsp_push_reg(&N,sizeof(long));
    bsp_sync();

    bsp_get(0,&N,0,&N,sizeof(long));
    bsp_sync();
    bsp_pop_reg(&N);

    double tg0, tg1, tD0, ts0, ts1;
    if (s == 0) tg0 = bsp_time();
    ///////////////////

    //Distribute the rows cyclically
    long n = nloc(p, s, N);

    //Random seed
    unsigned int seed = (unsigned int)(RAND_MAX*bsp_time()) ^ s;

    //Generation graph
    long* start; 
    start = vecalloci(n+1);
    long* cols = gen_graph(N, n, start, &seed);
 
    //Outlinks
    bsp_sync();
    if (s == 0) tD0 = bsp_time(); 
    uint8_t* D = outlinks_pr(N, n, cols, start);
    bsp_sync();
    if (s==0) {tg1 = bsp_time();}

    if (N <= 10) print_graphs(n, D, start, cols);
    
    //Initial vector
    double* u = initial_vector(N, n, &seed);

    // Allocations
    double* Du = vecallocd(n);
    bsp_push_reg(Du, n*sizeof(double));

    double* vals = vecallocd(N*sizeof(double));

    double* norms = vecallocd(p);
    bsp_push_reg(norms, p*sizeof(double));

    bsp_sync();

    //Determine which columns to get when multiplying with G
    bool* cols_get = vecallocb(N);
    initb(N, cols_get, false);
    for (int i = 0; i < n; i++)
        for (int k = start[i]; k < start[i+1]; k++)
            if (cols[k] % p != s) cols_get[cols[k]] = true;//(cols[k] % p != s) cols_get[cols[k]] = false;
    
    //Find r0
    double* r = vecallocd(n);
    double* GDu = vecallocd(n);
    mul_GD(N, n, u, D, start, cols, cols_get, Du, vals, GDu);
    for (int i = 0; i < n; i++)
        r[i] = 1 - u[i] + q*GDu[i];
    bsp_sync();
    vecfreed(GDu);

    //print_vecs_pr(n, u, r, n <= 10);

    //Iterate
    if (s==0) {ts0 = bsp_time();}
    long count = 0;
    double* GDr = vecallocd(n);

    double eps = pow(10, -12);
    double expectedCount = log(eps / N) / log(q) / 2;
    if (s == 0) printf("Expected %f iterations\n", expectedCount);

    while ((count < 0.95 * expectedCount) | (vec_total_norm(n, r, norms) >= eps))
        {
            add_vec(n, u, r);
            mul_GD(N, n, r, D, start, cols, cols_get, Du, vals, GDr);
            for (int i = 0; i < n; i++) r[i] = q*GDr[i];
            count++;
            bsp_sync();
        }

    //////////////////////////////////
    
    bsp_sync();
    
    if (s==0) ts1 = bsp_time();

    if (s == 0) printf("After %ld iterations:\n", count);
    //print_vecs_pr(n, u, r, n <= 10);   

    if (s == 0) {printf("Generation run-time: %f\n", tg1-tg0);
        printf("Finding outlinks: %f\n", tg1-tD0);
        printf("Solving run-time: %f\n", ts1 - tg1);
        printf("Time per iteration: %f\n", (ts1 - ts0) / count);
        printf("Total run-time: %f\n", (ts1 - tg0));
        printf("\n");}
    
    //Pop registries and free
    bsp_pop_reg(Du); bsp_pop_reg(norms);
    vecfreei(start); vecfreei(cols); vecfreei8(D);
    vecfreed(u); vecfreed(r); vecfreed(Du); vecfreed(norms);
    vecfreed(vals); vecfreed(GDr); vecfreeb(cols_get);
}

int main(int argc, char **argv) 
{
    bsp_init(bsp_pr, argc, argv);

    /* Sequential part */
    printf("How many processors do you want to use?\n");
    fflush(stdout);

    if ((scanf("%ld",&P) == 0) || (P > bsp_nprocs())){
        printf("Sorry, only %u processors available.\n",
                bsp_nprocs());
        fflush(stdout);
        exit(EXIT_FAILURE);
    }
    printf("Using %ld processors\n", P);
    bsp_pr();
    
    exit(EXIT_SUCCESS); 
}