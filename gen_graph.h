#include <stdlib.h>
#include <stdio.h>
#include "bspedupack.h"
#include "vec.h"

long* gen_graph(long N, long n, long* start, unsigned int* seed);
void print_graph(long N, long* cols, long* start);
void print_D(long N, uint8_t* D);
uint8_t* outlinks(long N, long n, long* cols, long* start);
void outlinks_noZeroes(long n, uint8_t* D);
