#include "gen_graph.h"

long min(long a, long b) 
{
    if (a < b) return a;
    else return b;
}

long rand_long(long N) {
    return (long)rand() % N;
}

int cmp_long(const void * a, const void * b) {
   return ( *(int*)a - *(int*)b );
}

void gen_nlinks(long N, long n, long* start)
{
    //for each row, determine a random number between 1 and min(N,10), which
    //is the number of links to that node, which determines the start array
    //The last element gives the total number of links.

    start[0] = 0;
    for (long i = 1; i <= n; i++) 
        start[i] = start[i-1] + rand_long(min(N,10)-1)+1;
}

void gen_cols(long N, long n, long* cols, long* start) 
{
    //for each row i, choose start[i+1] - start[i] columns, 
    //corresponding to the edges that connect to node i

    long nlinks;
    for (long i = 0; i < n; i++) 
    {   
        nlinks = start[i+1]-start[i];
        long* cols_i = vecalloci(nlinks);
        for (long j = 0; j < nlinks; j++)
            cols_i[j] = rand_long(N);
        qsort(cols_i, nlinks, sizeof(long), cmp_long);

        for (long k = start[i]; k < start[i+1]; k++) 
            cols[k] = cols_i[k-start[i]];
    }   
}

long* gen_graph(long N, long n, long* start) 
{
    //Generate a graph with a random number of inlinks between 1 and N,
    //store the (nxN) adjacency matrix in CRS format: 
    //for each row, store the column indices in cols,
    //and let start specify the indexes where the rows start.
    //start should be allocated (n+1) longs

    gen_nlinks(N, n, start);
    long* cols = vecalloci(start[n]);
    gen_cols(N, n, cols, start);

    return cols;
}

long* outlinks(long N, long n, long* cols, long* start) 
{ 
    long* D = vecalloci(N);
    initi(N, D, 0);
    for (long k = 0; k < start[n]; k++) 
        D[cols[k]]++;
    for (long i = 0; i < N; i++) 
        if (D[i] == 0) D[i] = 1;
    return D;
    // long outlinks;
    // for (long i = 0; i < N; i++) 
    // {
    //     outlinks = 0;
    //     for (long k = 0; k < start[n]; k++)
    //         if (cols[k] == i) outlinks++;

    //     if (outlinks == 0) D[i] = 1;
    //     else D[i] = outlinks;
    // }   
    // return D;
}

void print_graph(long N, long* cols, long* start) 
{
    printf("start[i]:\n");
    for (long i = 0; i <= N; i++)
        printf("%ld ", start[i]);
    printf("\n");
    printf("j[k]:\n");
    for (long i = 0; i < N; i++) 
    {
        for (long k = start[i]; k < start[i+1]; k++) 
            printf("%ld ", cols[k]);
        printf("\n");
    }
    printf("\n");
}

void print_D(long N, long* D) 
{
   printf("D[i]:\n");
    for (long i = 0; i < N; i++)
        printf("%ld ", D[i]);
    printf("\n"); 
}
