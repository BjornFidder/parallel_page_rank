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
    long maxlinks = min(N, 10);
    start[0] = 0;
    for (long i = 1; i <= n; i++) 
        start[i] = start[i-1] + rand_long(maxlinks)+1;
}

void gen_cols(long N, long n, long* cols, long* start) 
{
    //choose random values for the columns
    for (long k = 0; k < start[n]; k++) 
        cols[k] = rand_long(N);
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

void outlinks_noZeroes(long n, uint8_t* D) 
{
    for (long i = 0; i < n; i++) 
        if (D[i] == 0) D[i] = 1;
}

uint8_t* outlinks(long N, long n, long* cols, long* start) 
{ 
    uint8_t* D = calloc(N, sizeof(uint8_t));
    //initi(N, D, 0);
    for (long k = 0; k < start[n]; k++) 
        D[cols[k]]++;
    return D;
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

void print_D(long N, uint8_t* D) 
{
   printf("D[i]:\n");
    for (long i = 0; i < N; i++)
        printf("%d ", D[i]);
    printf("\n"); 
}
