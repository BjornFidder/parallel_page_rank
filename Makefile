CFLAGS= -std=c99 -Wall -O3 
LFLAGS= -lm

OBJSEQ = seq_pr.c gen_graph.c vec.c bspedupack.c

OBJBSP = bsp_pr.c gen_graph.c vec.c bspedupack.c

OBJBEN = bspbench.c bspedupack.c

seq: $(OBJGEN)
	bspcc $(CFLAGS) -o seq_pr $(OBJSEQ) $(LFLAGS)

bsp: $(OBJGEN)
	bspcc $(CFLAGS) -o bsp_pr $(OBJBSP) $(LFLAGS)

bench: $(OBJBEN)
	bspcc $(CFLAGS) -o bench $(OBJBEN) $(LFLAGS)