#include <stdlib.h>
#include <stdio.h>
#include "bspedupack.h"
#include "vec.h"

long rand_long(long N);
long* gen_graph(long N, long n, long* start);
void print_graph(long N, long* cols, long* start);
void print_D(long N, long* D);
long* outlinks(long N, long n, long* cols, long* start);
void outlinks_noZeroes(long n, long* D);
