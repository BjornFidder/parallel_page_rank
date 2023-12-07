#include "gen_graph.h"
#include "bsp_gen_graph.h"

long P;
double q = 0.85;

// long phi(long J, long p, long N) 
// {   
//     return J % p;
// }

// long j2J(long j, long s, long p, long N)
// {
//     //convert local index j on processor s to global index J
//     return s + j*p;
// }

// long J2j(long J, long p, long N) 
// {
//     //convert global index J to local index j
//     return J / p;
// }

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

double gen_rand_vec(long n, double* u) 
{
    //Generate n random vector components and return the sum
    double sum = 0;
    for (long i = 0; i<n; i++)
    {
        u[i] = (double)(rand() % (RAND_MAX / n));
        sum += u[i];
    }
    return sum;
}

void div_vec_pr(long n, double* x, long* y, double* z) 
{
    //divide vector x by vector y
    //giving the output in vector z
    for (long i = 0; i < n; i++)
        z[i] = x[i]/y[i];
}


void print_graphs(long n, long* D, long* start, long* cols)
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

void mul_GD(long N, long n, double* u, long* D, long* start, long* cols, double* Du, double* vals, double* v) 
{
    //compute v = G*Dinv*u in parallel
    long p = bsp_nprocs();

    initd(n, v, 0);
    div_vec_pr(n, u, D, Du);

    for (int i = 0; i < n; i++)
        for (int k = start[i]; k < start[i+1]; k++)
            bsp_get(cols[k] % p, Du, cols[k]/p * sizeof(double), &vals[k], sizeof(double));
    bsp_sync();
    
    for (int i = 0; i < n; i++)
        for (int k = start[i]; k < start[i+1]; k++)
            v[i] += vals[k];
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
        while (scanf("%ld", &N) != 1) 
            printf("Please input one number\n");
        printf("Running PageRank algorithm for N = %ld\n", N);
    }
    bsp_push_reg(&N,sizeof(long));
    bsp_sync();

    bsp_get(0,&N,0,&N,sizeof(long));
    bsp_sync();
    bsp_pop_reg(&N);

    double t0, t1, t2;
    if (s == 0) t0 = bsp_time();
    ///////////////////

    //Distribute the rows cyclically
    long n = nloc(p, s, N);

    long* start; 
    start = vecalloci(n+1);
    long* cols = gen_graph(N, n, start);

    ////compute local outlinks
    long* ds = vecalloci(n*p);
    bsp_push_reg(ds, n*p*sizeof(long));
    bsp_sync();
    
    long* d = outlinks(N, n, cols, start); //local outlinks

    //broadcast local outlinks
    for (long J = 0; J < N; J++)
        bsp_put(J%p, &d[J], ds, (J/p + s*n)*sizeof(long), sizeof(long));
       
    vecfreei(d);
    bsp_sync();
    
    //compute total outlinks for local rows
    long* D = vecalloci(n);
    initi(n, D, 0);
    for (long i = 0; i < n; i++) 
        for (long t = 0; t < p; t++)
            D[i] += ds[i + t*n];
    
    outlinks_noZeroes(n, D);
    bsp_pop_reg(ds);
    vecfreei(ds);
    bsp_sync();

    if (s==0) {t1 = bsp_time();}
    
    //  Generate random initial vector u0 //
    //register
    double* sums = vecallocd(p);
    bsp_push_reg(sums, p*sizeof(double));
    bsp_sync();

    //generate and distribute sum to all processors
    double* u = vecallocd(n);
    double sum = gen_rand_vec(n, u);
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


    // Allocations
    double* Du = vecallocd(n);
    bsp_push_reg(Du, n*sizeof(double));

    double* vals = vecallocd(start[n]*sizeof(double));

    double* norms = vecallocd(p);
    bsp_push_reg(norms, p*sizeof(double));

    bsp_sync();

    //Find r0
    double* r = vecallocd(n);
    double* GDu = vecallocd(n);
    mul_GD(N, n, u, D, start, cols, Du, vals, GDu);
    for (int i = 0; i < n; i++)
        r[i] = 1 - u[i] + q*GDu[i];
    bsp_sync();
    vecfreed(GDu);

    //print_vecs_pr(n, u, r, n <= 10);

    //Iterate
    long count = 0;
    double* GDr = vecallocd(n);
    while (vec_total_norm(n, r, norms) >= pow(10, -12))
        {
            add_vec(n, u, r);
            mul_GD(N, n, r, D, start, cols, Du, vals, GDr);
            for (int i = 0; i < n; i++) r[i] = q*GDr[i];
            count++;
            bsp_sync();
        }

    //////////////////////////////////
    
    bsp_sync();
    
    if (s==0) t2 = bsp_time();

    if (s == 0) printf("After %ld iterations:\n", count);
    //print_vecs_pr(n, u, r, n <= 10);   

    if (s == 0) {printf("Generation run-time: %f\n", t1-t0);
        printf("Solving run-time: %f\n", t2 - t1);}

    bsp_pop_reg(Du);
    bsp_pop_reg(norms);
    vecfreei(start);
    vecfreei(cols);
    vecfreei(D);
    vecfreed(u);
    vecfreed(r);
    vecfreed(Du);
    vecfreed(norms);
    vecfreed(vals);
    vecfreed(GDr);
}

int main(int argc, char **argv) 
{
    bsp_init(bsp_pr, argc, argv);

    /* Sequential part */
    printf("How many processors do you want to use?\n");
    fflush(stdout);

    scanf("%ld",&P);
    if (P > bsp_nprocs()){
        printf("Sorry, only %u processors available.\n",
                bsp_nprocs());
        fflush(stdout);
        exit(EXIT_FAILURE);
    }
    printf("Using %ld processors\n", P);
    bsp_pr();
    
    exit(EXIT_SUCCESS); 
}