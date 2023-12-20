# parallel_page_rank
 Parallel implementation of the iterative PageRank algorithm, using MulticoreBSP for C (http://www.multicorebsp.com/documentation/quickC/) and the BSPedupack by Rob Bisseling.
 This project was an assignment for the master course Parallel Algorithms at Utrecht University, in 2023 lectured by Jonas Thies and Martin van Gijzen from TU Delft.

 The sequential code is given in seq_pr.c, where pr is short for PageRank, while the parallel BSP implementation is given in bsp_pr.c.

 The code is benchmarked on the Snellius supercomputer (https://www.surf.nl/en/services/snellius-the-national-supercomputer). The results from the tests are found in the folders bench_snellius and bsp_snellius, while in the folder plots one can find graphical representations.

 See also the file report.pdf for the report that describes and motivates the assignment.
