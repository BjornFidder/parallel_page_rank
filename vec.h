#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <stdint.h>

void print_vec(long N, double* x);
void gen_rand_prob_vec(long N, double* u);
double norm_vec(long N, double* x);
void add_vec(long N, double* x, double* y);
double* div_vec(long N, double* x, uint8_t* y);
double* mul_mat_vec(long N, long* cols, long* start, double* x);
void initi(long N, long* x, long a);
void initi8(long N, uint8_t* x, uint8_t a);
void initb(long N, bool* x, bool a);
void initd(long N, double* x, double a);