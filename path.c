#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <omp.h>
#include <mpi.h>
#include "mt19937p.h"

//ldoc on
/**
 * # The basic recurrence
 *
 * At the heart of the method is the following basic recurrence.
 * If $l_{ij}^s$ represents the length of the shortest path from
 * $i$ to $j$ that can be attained in at most $2^s$ steps, then
 * $$
 *   l_{ij}^{s+1} = \min_k \{ l_{ik}^s + l_{kj}^s \}.
 * $$
 * That is, the shortest path of at most $2^{s+1}$ hops that connects
 * $i$ to $j$ consists of two segments of length at most $2^s$, one
 * from $i$ to $k$ and one from $k$ to $j$.  Compare this with the
 * following formula to compute the entries of the square of a
 * matrix $A$:
 * $$
 *   a_{ij}^2 = \sum_k a_{ik} a_{kj}.
 * $$
 * These two formulas are identical, save for the niggling detail that
 * the latter has addition and multiplication where the former has min
 * and addition.  But the basic pattern is the same, and all the
 * tricks we learned when discussing matrix multiplication apply -- or
 * at least, they apply in principle.  I'm actually going to be lazy
 * in the implementation of `square`, which computes one step of
 * this basic recurrence.  I'm not trying to do any clever blocking.
 * You may choose to be more clever in your assignment, but it is not
 * required.
 *
 * The return value for `square` is true if `l` and `lnew` are
 * identical, and false otherwise.
 */

int square(int n,               // Number of nodes
           int* restrict l,     // Partial distance at step s
           int* restrict lnew)  // Partial distance at step s+1
{
    int done = 1;
    for (int j = 0; j < n; ++j) {
        for (int i = 0; i < n; ++i) {
            int lij = lnew[j*n+i];
            for (int k = 0; k < n; ++k) {
                int lik = l[k*n+i];
                int lkj = l[j*n+k];
                if (lik + lkj < lij) {
                    lij = lik+lkj;
                    done = 0;
                }
            }
            lnew[j*n+i] = lij;
        }
    }
    return done;
}

int square_(int n, int* restrict l, int* restrict la, int* restrict lb)
{
    int done = 1;
    for (int j = 0; j < n; ++j) {
        for (int i = 0; i < n; ++i) {
            int lij = l[j*n+i];
            for (int k = 0; k < n; ++k) {
                int lik = la[k*n+i];
                int lkj = lb[j*n+k];
                if (lik + lkj < lij) {
                    lij = lik+lkj;
                    done = 0;
                }
            }
            l[j*n+i] = lij;
        }
    }
    return done;
}

/**
 *
 * The value $l_{ij}^0$ is almost the same as the $(i,j)$ entry of
 * the adjacency matrix, except for one thing: by convention, the
 * $(i,j)$ entry of the adjacency matrix is zero when there is no
 * edge between $i$ and $j$; but in this case, we want $l_{ij}^0$
 * to be "infinite".  It turns out that it is adequate to make
 * $l_{ij}^0$ longer than the longest possible shortest path; if
 * edges are unweighted, $n+1$ is a fine proxy for "infinite."
 * The functions `infinitize` and `deinfinitize` convert back 
 * and forth between the zero-for-no-edge and $n+1$-for-no-edge
 * conventions.
 */

static inline void infinitize(int n, int* l)
{
    for (int i = 0; i < n*n; ++i)
        if (l[i] == 0)
            l[i] = n+1;
}

static inline void deinfinitize(int n, int* l)
{
    for (int i = 0; i < n*n; ++i)
        if (l[i] == n+1)
            l[i] = 0;
}

/**
 *
 * Of course, any loop-free path in a graph with $n$ nodes can
 * at most pass through every node in the graph.  Therefore,
 * once $2^s \geq n$, the quantity $l_{ij}^s$ is actually
 * the length of the shortest path of any number of hops.  This means
 * we can compute the shortest path lengths for all pairs of nodes
 * in the graph by $\lceil \lg n \rceil$ repeated squaring operations.
 *
 * The `shortest_path` routine attempts to save a little bit of work
 * by only repeatedly squaring until two successive matrices are the
 * same (as indicated by the return value of the `square` routine).
 */

void shortest_paths(int n, int* restrict l)
{
    // Generate l_{ij}^0 from adjacency matrix representation
    infinitize(n, l);
    for (int i = 0; i < n*n; i += n+1)
        l[i] = 0;

    // Repeated squaring until nothing changes
    int* restrict lnew = (int*) calloc(n*n, sizeof(int));
    memcpy(lnew, l, n*n * sizeof(int));
    for (int done = 0; !done; ) {
        done = square(n, l, lnew);
        memcpy(l, lnew, n*n * sizeof(int));
    }
    free(lnew);
    deinfinitize(n, l);
}

/**
 * # The random graph model
 *
 * Of course, we need to run the shortest path algorithm on something!
 * For the sake of keeping things interesting, let's use a simple random graph
 * model to generate the input data.  The $G(n,p)$ model simply includes each
 * possible edge with probability $p$, drops it otherwise -- doesn't get much
 * simpler than that.  We use a thread-safe version of the Mersenne twister
 * random number generator in lieu of coin flips.
 */

int* gen_graph(int n, double p)
{
    int* l = calloc(n*n, sizeof(int));
    struct mt19937p state;
    sgenrand(10302011UL, &state);
    for (int j = 0; j < n; ++j) {
        for (int i = 0; i < n; ++i)
            l[j*n+i] = (genrand(&state) < p);
        l[j*n+j] = 0;
    }
    return l;
}

/**
 * # Result checks
 *
 * Simple tests are always useful when tuning code, so I have included
 * two of them.  Since this computation doesn't involve floating point
 * arithmetic, we should get bitwise identical results from run to
 * run, even if we do optimizations that change the associativity of
 * our computations.  The function `fletcher16` computes a simple
 * [simple checksum][wiki-fletcher] over the output of the
 * `shortest_paths` routine, which we can then use to quickly tell
 * whether something has gone wrong.  The `write_matrix` routine
 * actually writes out a text representation of the matrix, in case we
 * want to load it into MATLAB to compare results.
 *
 * [wiki-fletcher]: http://en.wikipedia.org/wiki/Fletcher's_checksum
 */

int fletcher16(int* data, int count)
{
    int sum1 = 0;
    int sum2 = 0;
    for(int index = 0; index < count; ++index) {
          sum1 = (sum1 + data[index]) % 255;
          sum2 = (sum2 + sum1) % 255;
    }
    return (sum2 << 8) | sum1;
}

void write_matrix(const char* fname, int n, int* a)
{
    FILE* fp = fopen(fname, "w+");
    if (fp == NULL) {
        fprintf(stderr, "Could not open output file: %s\n", fname);
        exit(-1);
    }
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) 
            fprintf(fp, "%d ", a[j*n+i]);
        fprintf(fp, "\n");
    }
    fclose(fp);
}

/**
 * # MPI Comunication Helper Functions & Functions
 */
void copy_out(int* restrict dst, const int* restrict src, int nx, int ny, int stride) {
    for (int iy = 0; iy < ny; ++iy) {
        for (int ix = 0; ix < nx; ++ix) {
            dst[iy*nx+ix] = src[iy*stride+ix];
        }
    }
}

void copy_in(const int* restrict src, int* restrict dst, int nx, int ny, int stride) {
    for (int iy = 0; iy < ny; ++iy) {
        for (int ix = 0; ix < nx; ++ix) {
            dst[iy*stride+ix] = src[iy*nx+ix];
        }
    }
}

void send_initial_graph(int* restrict l, int world_rank, int world_size, 
                        int size_block, int n,
                        const int* restrict offset_x, const int* restrict offset_y) {
    int len = size_block * size_block;
    for (int i = 1; i < world_size; ++i) {
        int* buffer = malloc(len * sizeof(int));
        copy_out(buffer, l + offset_x[i-1] + offset_y[i-1] * n, size_block, size_block, n);
        MPI_Send(buffer, len, MPI_INT, i, 
                0, MPI_COMM_WORLD);
        free(buffer);
        printf("Process %d sent initial graph to process %d\n", world_rank, i);
    }
}

void recv_initial_graph(int* restrict l, int world_rank, int size_block) {
    int len = size_block * size_block;
    MPI_Recv(l, len, MPI_INT, 0, 
            0, MPI_COMM_WORLD, NULL);
    printf("Process %d recv initial graph from process %d\n", world_rank, 0);
}

void send_current_graph(int* restrict l, int world_rank, int world_size, 
                        int size_block, int n,
                        const int* restrict offset_ax, const int* restrict offset_ay, 
                        const int* restrict offset_bx, const int* restrict offset_by) {
    int len = size_block * size_block;
    for (int i = 1; i < world_size; ++i) {
        int* buffer = malloc(len * sizeof(int));
        copy_out(buffer, l + offset_ax[i-1] + offset_ay[i-1] * n, size_block, size_block, n);
        MPI_Send(buffer, len, MPI_INT, i, 
                0, MPI_COMM_WORLD);
        copy_out(buffer, l + offset_bx[i-1] + offset_by[i-1] * n, size_block, size_block, n);
        MPI_Send(buffer, len, MPI_INT, i, 
                0, MPI_COMM_WORLD);
        free(buffer);
        printf("Process %d sent current graph to process %d\n", world_rank, i);
    }
}

void recv_current_graph(int* restrict la, int* restrict lb, int world_rank, int size_block) {
    int len = size_block * size_block;
    MPI_Recv(la, len, MPI_INT, 0, 
            0, MPI_COMM_WORLD, NULL);
    MPI_Recv(lb, len, MPI_INT, 0, 
            0, MPI_COMM_WORLD, NULL);
    printf("Process %d recv current graph from process %d\n", world_rank, 0);
}

void send_updated_graph(int* restrict l, int* done, int world_rank, int size_block) {
    int len = size_block * size_block;
    MPI_Send(done, 1, MPI_INT, 0, 
            0, MPI_COMM_WORLD);
    MPI_Send(l, len, MPI_INT, 0, 
            0, MPI_COMM_WORLD);
    printf("Process %d sent updated graph to process %d\n", world_rank, 0);
}

int recv_updated_graph(int* restrict l, int world_rank, int world_size, 
                        int size_block, int n,
                        const int* restrict offset_x, const int* restrict offset_y) {
    int done = 1;
    for (int i = 1; i < world_size; ++i) {
        int done_;
        int len = size_block * size_block;
        int* buffer = malloc(len * sizeof(int));
        MPI_Recv(&done_, 1, MPI_INT, i, 
                0, MPI_COMM_WORLD, NULL);
        done = done && done_;
        MPI_Recv(buffer, len, MPI_INT, i, 
                0, MPI_COMM_WORLD, NULL);
        copy_in(buffer, l + offset_x[i-1] + offset_y[i-1] * n, size_block, size_block, n);
        free(buffer);
        printf("Process %d recv updated graph from process %d\n", world_rank, i);
    }
    return done;
}

/**
 * # The `main` event
 */

const char* usage =
    "path.x -- Parallel all-pairs shortest path on a random graph\n"
    "Flags:\n"
    "  - n -- number of nodes (200)\n"
    "  - p -- probability of including edges (0.05)\n"
    "  - i -- file name where adjacency matrix should be stored (none)\n"
    "  - o -- file name where output matrix should be stored (none)\n";

int main(int argc, char** argv)
{
    // Initialize MPI environment
    MPI_Init(&argc, &argv);

    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    printf("---- Starting Process %d / %d ----\n", world_rank, world_size);



    int n    = 200;            // Number of nodes
    double p = 0.05;           // Edge probability
    const char* ifname = NULL; // Adjacency matrix file name
    const char* ofname = NULL; // Distance matrix file name

    // Option processing
    extern char* optarg;
    const char* optstring = "hn:d:p:o:i:";
    int c;
    while ((c = getopt(argc, argv, optstring)) != -1) {
        switch (c) {
        case 'h':
            fprintf(stderr, "%s", usage);
            return -1;
        case 'n': n = atoi(optarg); break;
        case 'p': p = atof(optarg); break;
        case 'o': ofname = optarg;  break;
        case 'i': ifname = optarg;  break;
        }
    }



    int num_block;
    switch (world_size-1) {
        case 1: num_block = 1; break;
        case 4: num_block = 2; break;
        case 9: num_block = 3; break;
    }
    int size_block = n / num_block;
    int* offset_x = malloc((world_size-1) * sizeof(int));
    int* offset_y = malloc((world_size-1) * sizeof(int));
    int* offset_ax = malloc((world_size-1) * sizeof(int));
    int* offset_ay = malloc((world_size-1) * sizeof(int));
    int* offset_bx = malloc((world_size-1) * sizeof(int));
    int* offset_by = malloc((world_size-1) * sizeof(int));
    for (int i = 0; i < num_block; ++i) {
        for (int j = 0; j < num_block; ++j) {
            offset_x[j+i*num_block] = n * i / num_block;
            offset_y[j+i*num_block] = n * j / num_block;
            offset_ax[j+i*num_block] = n * i / num_block;
            offset_ay[j+i*num_block] = n * ((i+j)%num_block) / num_block;
            offset_bx[j+i*num_block] = n * ((i+j)%num_block) / num_block;
            offset_by[j+i*num_block] = n * j / num_block;
        }
    }
    if (world_rank == 0) {
        printf("Process %d does output file IO and command line IO\n", world_rank);
        for (int i = 0; i < world_size-1; ++i) {
            printf("Process %d has offset_x[%d]: %d, offset_y[%d]: %d, offset_ax[%d]: %d, offset_ay[%d]: %d, offset_bx[%d]: %d, offset_by[%d]: %d\n",
                i+1, i, offset_x[i], i, offset_y[i], i, offset_ax[i], i, offset_ay[i], i, offset_bx[i], i, offset_by[i]);
        }
        printf("There are %d blocks in each direction\n", num_block);
        printf("Each block has size %d\n", size_block);
    }



    if (world_rank == 0) {
        // Graph generation + output
        int* l = gen_graph(n, p);
        if (ifname)
            write_matrix(ifname,  n, l);

        // Generate l_{ij}^0 from adjacency matrix representation
        infinitize(n, l);
        for (int i = 0; i < n*n; i += n+1)
            l[i] = 0;



        // Time the shortest paths code
        double t0 = omp_get_wtime();

        send_initial_graph(l, world_rank, world_size, size_block, n,
                        offset_x, offset_y);
        
        for (int done = 0; !done; ) {
            send_current_graph(l, world_rank, world_size, size_block, n,
                    offset_ax, offset_ay, offset_bx, offset_by);

            done = recv_updated_graph(l, world_rank, world_size, size_block, n,
                    offset_x, offset_y);

            for (int i = 1; i < world_size; ++i) {
                MPI_Send(&done, 1, MPI_INT, i, 
                        0, MPI_COMM_WORLD);
            }
        }

        double t1 = omp_get_wtime();

        

        deinfinitize(n, l);

        printf("n:     %d\n", n);
        printf("p:     %g\n", p);
        printf("Time:  %g\n", t1-t0);
        printf("Check: %X\n", fletcher16(l, n*n));

        // Generate output file
        if (ofname)
            write_matrix(ofname, n, l);

        // Clean up
        free(l);
    }
    else if (world_rank != 0) {
        int world_rank_0 = world_rank - 1;
        int r_rank_0 = (world_rank_0 / num_block) * num_block + (world_rank_0 + 1) % num_block;
        int l_rank_0 = (world_rank_0 / num_block) * num_block + (world_rank_0 + num_block - 1) % num_block;
        int t_rank_0 = (world_rank_0 + world_size - num_block - 1) % (world_size - 1);
        int b_rank_0 = (world_rank_0 + num_block) % (world_size - 1);
        printf("Process %d: left has world rank %d, right has world rank %d, top has world rank %d, bottom has world rank %d\n", 
            world_rank_0+1, l_rank_0+1, r_rank_0+1, t_rank_0+1, b_rank_0+1);

        int* l = malloc(size_block * size_block * sizeof(int));
        recv_initial_graph(l, world_rank, size_block);

        int* la = malloc(size_block * size_block * sizeof(int));
        int* lb = malloc(size_block * size_block * sizeof(int));

        for (int done = 0; !done; ) {
            recv_current_graph(la, lb, world_rank, size_block);

            int len = size_block * size_block;
            for (int iter = 0; iter < num_block; ++iter) {
                int done_ = 1;
                if (iter != 0) {
                    if (world_rank-1 != l_rank_0) {
                        if ((world_rank-1) % num_block == 0) {
                            int* buffer = malloc(len * sizeof(int));
                            copy_out(buffer, la, size_block, size_block, size_block);
                            MPI_Send(buffer, len, MPI_INT, l_rank_0+1, 
                                    0, MPI_COMM_WORLD);
                            MPI_Recv(la, len, MPI_INT, r_rank_0+1, 
                                    0, MPI_COMM_WORLD, NULL);
                            free(buffer);
                        }
                        else {
                            int* buffer = malloc(len * sizeof(int));
                            copy_out(buffer, la, size_block, size_block, size_block);
                            MPI_Recv(la, len, MPI_INT, r_rank_0+1, 
                                    0, MPI_COMM_WORLD, NULL);
                            MPI_Send(buffer, len, MPI_INT, l_rank_0+1, 
                                    0, MPI_COMM_WORLD);
                            free(buffer);
                        }
                    }

                    if (world_rank-1 != t_rank_0) {
                        if ((world_rank-1) / num_block == 0) {
                            int* buffer = malloc(len * sizeof(int));
                            copy_out(buffer, lb, size_block, size_block, size_block);
                            MPI_Send(buffer, len, MPI_INT, t_rank_0+1, 
                                    0, MPI_COMM_WORLD);
                            MPI_Recv(lb, len, MPI_INT, b_rank_0+1, 
                                    0, MPI_COMM_WORLD, NULL);
                            free(buffer);
                        }
                        else {
                            int* buffer = malloc(len * sizeof(int));
                            copy_out(buffer, lb, size_block, size_block, size_block);
                            MPI_Recv(lb, len, MPI_INT, b_rank_0+1, 
                                    0, MPI_COMM_WORLD, NULL);
                            MPI_Send(buffer, len, MPI_INT, t_rank_0+1, 
                                    0, MPI_COMM_WORLD);
                            free(buffer);
                        }
                    }
                }

                done = done_ && square_(size_block, l, la, lb);
            }

            send_updated_graph(l, &done, world_rank, size_block);

            MPI_Recv(&done, 1, MPI_INT, 0, 
                    0, MPI_COMM_WORLD, NULL);
        }
    }

    // Finalize MPI environment.
    MPI_Finalize();
    
    return 0;
}
