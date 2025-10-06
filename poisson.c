#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>

/**
 * poisson.c
 * Implementation of a Poisson solver.
 *
 * This template handles the basic program launch, argument parsing, and memory
 * allocation required to implement the solver *at its most basic level*. You
 * will likely need to allocate more memory, add threading support, account for
 * cache locality, etc...
 *
 * TODO:
 * 1 - Read through this example, understand what it does and what it gives you
 *     to work with.
 * 2 - Implement the basic algorithm and get a correct output.
 * 3 - Add a timer to track how long your execution takes.
 * 4 - Profile your solution and identify weaknesses.
 * 5 - Improve it!
 * 6 - Remember that this is now *your* code and *you* should modify it however
 *     needed to solve the assignment.
 *
 * See the lab notes for a guide on profiling and an introduction to
 * multithreading (see also threads.c which is reference by the lab notes).
 */

/**
 * @brief Solve Poissons equation for a given cube with mixed boundary conditions.
 *
 * @param n             The edge length of the cube. n^3 number of elements.
 * @param source        Pointer to the source term cube, a.k.a. forcing function.
 * @param iterations    Number of iterations to perform.
 * @param threads       Number of threads to use for solving.
 * @param delta         Grid spacing.
 * @param debug         Set to true when operating in debug mode to enable verbose logging.
 * @return double*      Solution to Poissons equation.  Caller must free.
 */

 
 double* run_poisson(int n, double *source, int iterations, int threads, float delta, bool debug)
{
    if (debug)
    {
        printf("Starting solver with:\n"
               "n = %i\n"
               "iterations = %i\n"
               "threads = %i\n"
               "delta = %f\n",
               n, iterations, threads, delta);
    }

    // Allocate two buffers for double buffering
    double *curr = (double*)calloc(n * n * n, sizeof(double));
    double *next = (double*)calloc(n * n * n, sizeof(double));

    if (curr == NULL || next == NULL)
    {
        fprintf(stderr, "Error: ran out of memory when trying to allocate %i sized cube\n", n);
        exit(EXIT_FAILURE);
    }

    // Initialise Dirichlet boundary conditions
    // Left side (i=0): V = -2.0
    // Right side (i=n-1): V = +1.0
    for (int k = 0; k < n; k++)
    {
        for (int j = 0; j < n; j++)
        {
            curr[(k * n + j) * n + 0] = -2.0;
            curr[(k * n + j) * n + (n-1)] = 1.0;
            next[(k * n + j) * n + 0] = -2.0;
            next[(k * n + j) * n + (n-1)] = 1.0;
        }
    }

    double delta_sq = delta * delta;

    // Jacobi iteration
    for (int iter = 0; iter < iterations; iter++)
    {
        // Update interior points (not on any boundary)
        for (int k = 1; k < n-1; k++)
        {
            for (int j = 1; j < n-1; j++)
            {
                for (int i = 1; i < n-1; i++)
                {
                    int idx = (k * n + j) * n + i;
                    
                    double v_im1 = curr[(k * n + j) * n + (i-1)];
                    double v_ip1 = curr[(k * n + j) * n + (i+1)];
                    double v_jm1 = curr[(k * n + (j-1)) * n + i];
                    double v_jp1 = curr[(k * n + (j+1)) * n + i];
                    double v_km1 = curr[((k-1) * n + j) * n + i];
                    double v_kp1 = curr[((k+1) * n + j) * n + i];
                    
                    next[idx] = (v_ip1 + v_im1 + v_jp1 + v_jm1 + v_kp1 + v_km1 
                                 - delta_sq * source[idx]) / 6.0;
                }
            }
        }

        // Update Neumann boundaries using ghost points
        // For zero gradient: ghost point equals interior point across boundary
        
        // j=0 boundary (excluding Dirichlet boundaries at i=0 and i=n-1)
        for (int k = 1; k < n-1; k++)
        {
            for (int i = 1; i < n-1; i++)
            {
                int idx = (k * n + 0) * n + i;
                
                double v_im1 = curr[(k * n + 0) * n + (i-1)];
                double v_ip1 = curr[(k * n + 0) * n + (i+1)];
                double v_jm1 = curr[(k * n + 1) * n + i];      // Ghost: j=-1 = j=1
                double v_jp1 = curr[(k * n + 1) * n + i];      // j=1
                double v_km1 = curr[((k-1) * n + 0) * n + i];
                double v_kp1 = curr[((k+1) * n + 0) * n + i];
                
                next[idx] = (v_ip1 + v_im1 + v_jp1 + v_jm1 + v_kp1 + v_km1 
                             - delta_sq * source[idx]) / 6.0;
            }
        }
        
        // j=n-1 boundary
        for (int k = 1; k < n-1; k++)
        {
            for (int i = 1; i < n-1; i++)
            {
                int idx = (k * n + (n-1)) * n + i;
                
                double v_im1 = curr[(k * n + (n-1)) * n + (i-1)];
                double v_ip1 = curr[(k * n + (n-1)) * n + (i+1)];
                double v_jm1 = curr[(k * n + (n-2)) * n + i];  // j=n-2
                double v_jp1 = curr[(k * n + (n-2)) * n + i];  // Ghost: j=n = j=n-2
                double v_km1 = curr[((k-1) * n + (n-1)) * n + i];
                double v_kp1 = curr[((k+1) * n + (n-1)) * n + i];
                
                next[idx] = (v_ip1 + v_im1 + v_jp1 + v_jm1 + v_kp1 + v_km1 
                             - delta_sq * source[idx]) / 6.0;
            }
        }
        
        // k=0 boundary
        for (int j = 1; j < n-1; j++)
        {
            for (int i = 1; i < n-1; i++)
            {
                int idx = (0 * n + j) * n + i;
                
                double v_im1 = curr[(0 * n + j) * n + (i-1)];
                double v_ip1 = curr[(0 * n + j) * n + (i+1)];
                double v_jm1 = curr[(0 * n + (j-1)) * n + i];
                double v_jp1 = curr[(0 * n + (j+1)) * n + i];
                double v_km1 = curr[(1 * n + j) * n + i];      // Ghost: k=-1 = k=1
                double v_kp1 = curr[(1 * n + j) * n + i];      // k=1
                
                next[idx] = (v_ip1 + v_im1 + v_jp1 + v_jm1 + v_kp1 + v_km1 
                             - delta_sq * source[idx]) / 6.0;
            }
        }
        
        // k=n-1 boundary
        for (int j = 1; j < n-1; j++)
        {
            for (int i = 1; i < n-1; i++)
            {
                int idx = ((n-1) * n + j) * n + i;
                
                double v_im1 = curr[((n-1) * n + j) * n + (i-1)];
                double v_ip1 = curr[((n-1) * n + j) * n + (i+1)];
                double v_jm1 = curr[((n-1) * n + (j-1)) * n + i];
                double v_jp1 = curr[((n-1) * n + (j+1)) * n + i];
                double v_km1 = curr[((n-2) * n + j) * n + i];  // k=n-2
                double v_kp1 = curr[((n-2) * n + j) * n + i];  // Ghost: k=n = k=n-2
                
                next[idx] = (v_ip1 + v_im1 + v_jp1 + v_jm1 + v_kp1 + v_km1 
                             - delta_sq * source[idx]) / 6.0;
            }
        }
        
        // Handle corner and edge cases where multiple Neumann boundaries meet
        // Edges: j=0, k=0
        for (int i = 1; i < n-1; i++)
        {
            int idx = (0 * n + 0) * n + i;
            double v_im1 = curr[(0 * n + 0) * n + (i-1)];
            double v_ip1 = curr[(0 * n + 0) * n + (i+1)];
            double v_jm1 = curr[(0 * n + 1) * n + i];      // Ghost: j=-1 = j=1
            double v_jp1 = curr[(0 * n + 1) * n + i];
            double v_km1 = curr[(1 * n + 0) * n + i];      // Ghost: k=-1 = k=1
            double v_kp1 = curr[(1 * n + 0) * n + i];
            
            next[idx] = (v_ip1 + v_im1 + v_jp1 + v_jm1 + v_kp1 + v_km1 
                         - delta_sq * source[idx]) / 6.0;
        }
        
        // Edges: j=0, k=n-1
        for (int i = 1; i < n-1; i++)
        {
            int idx = ((n-1) * n + 0) * n + i;
            double v_im1 = curr[((n-1) * n + 0) * n + (i-1)];
            double v_ip1 = curr[((n-1) * n + 0) * n + (i+1)];
            double v_jm1 = curr[((n-1) * n + 1) * n + i];  // Ghost: j=-1 = j=1
            double v_jp1 = curr[((n-1) * n + 1) * n + i];
            double v_km1 = curr[((n-2) * n + 0) * n + i];  // k=n-2
            double v_kp1 = curr[((n-2) * n + 0) * n + i];  // Ghost: k=n = k=n-2
            
            next[idx] = (v_ip1 + v_im1 + v_jp1 + v_jm1 + v_kp1 + v_km1 
                         - delta_sq * source[idx]) / 6.0;
        }
        
        // Edges: j=n-1, k=0
        for (int i = 1; i < n-1; i++)
        {
            int idx = (0 * n + (n-1)) * n + i;
            double v_im1 = curr[(0 * n + (n-1)) * n + (i-1)];
            double v_ip1 = curr[(0 * n + (n-1)) * n + (i+1)];
            double v_jm1 = curr[(0 * n + (n-2)) * n + i];  // j=n-2
            double v_jp1 = curr[(0 * n + (n-2)) * n + i];  // Ghost: j=n = j=n-2
            double v_km1 = curr[(1 * n + (n-1)) * n + i];  // Ghost: k=-1 = k=1
            double v_kp1 = curr[(1 * n + (n-1)) * n + i];
            
            next[idx] = (v_ip1 + v_im1 + v_jp1 + v_jm1 + v_kp1 + v_km1 
                         - delta_sq * source[idx]) / 6.0;
        }
        
        // Edges: j=n-1, k=n-1
        for (int i = 1; i < n-1; i++)
        {
            int idx = ((n-1) * n + (n-1)) * n + i;
            double v_im1 = curr[((n-1) * n + (n-1)) * n + (i-1)];
            double v_ip1 = curr[((n-1) * n + (n-1)) * n + (i+1)];
            double v_jm1 = curr[((n-1) * n + (n-2)) * n + i];  // j=n-2
            double v_jp1 = curr[((n-1) * n + (n-2)) * n + i];  // Ghost: j=n = j=n-2
            double v_km1 = curr[((n-2) * n + (n-1)) * n + i];  // k=n-2
            double v_kp1 = curr[((n-2) * n + (n-1)) * n + i];  // Ghost: k=n = k=n-2
            
            next[idx] = (v_ip1 + v_im1 + v_jp1 + v_jm1 + v_kp1 + v_km1 
                         - delta_sq * source[idx]) / 6.0;
        }
        
        // Dirichlet boundaries remain fixed (already initialised)
        // Copy them to next buffer
        for (int k = 0; k < n; k++)
        {
            for (int j = 0; j < n; j++)
            {
                next[(k * n + j) * n + 0] = -2.0;
                next[(k * n + j) * n + (n-1)] = 1.0;
            }
        }

        // Swap buffers
        double *temp = curr;
        curr = next;
        next = temp;
    }

    free(next);

    if (debug)
    {
        printf("Finished solving.\n");
    }

    return curr;
}
