#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>

#include "poisson.h"

/**
 * CAUTION: The command parsing logic below must function
 * as provided in order to run benchmark testing.
 */ 

// Default location and amplitude of source point
#define DEFAULT_X(n) (n/2)
#define DEFAULT_Y(n) (n/2)
#define DEFAULT_Z(n) (n/2)
#define DEFAULT_AMPLITUDE 1.0

static void populate_source_cube(double* source, FILE* file, int n);

int main (int argc, char **argv)
{
    // Default settings for solver
    int iterations = 10;
    int n = 5;
    int threads = 1;
    float delta = 1;
    FILE* source_coordinates_file = NULL;
    bool debug = false;

    int opt;

    // parse the command line arguments
    while  ( (opt = getopt (argc, argv, "h:n:i:s:t:d:")) != -1)
    {
        switch(opt)
        {
        case 'h':
            printf ("Usage: poisson [-n size] [-s source-coordinates-file] [-i iterations] "
                "[-t threads] [-d] (for debug mode)\n");
            return EXIT_SUCCESS;
        case 'n':
            n = atoi(optarg);
            break;
        case 'i':
            iterations = atoi(optarg);
            break;
        case 's':
            source_coordinates_file = fopen(optarg, "r");
            if (source_coordinates_file == NULL)
            {
                fprintf (stderr, "Could not open source coordinate file or file does not exist.\n");
                exit (EXIT_FAILURE);
            }
            break;
        case 't':
            threads = atoi(optarg);
            break;
        case 'd':
            debug = true;
            break;
        default:
            fprintf (stderr, "Usage: poisson [-n size] [-s-source-coordinates-file] [-i iterations] "
                "[-t threads] [-d] (for debug mode)\n");
            exit (EXIT_FAILURE);
        }
    }

    // Ensure we have an odd sized cube
    if (n % 2 == 0)
    {
        fprintf (stderr, "Error: n should be an odd number!\n");
        exit (EXIT_FAILURE);
    }

    // Create the source cube
    double *source = (double*)calloc (n * n * n, sizeof (double));
    if (source == NULL)
    {
        fprintf (stderr, "Error: failed to allocated source cube (n=%i)\n", n);
        exit (EXIT_FAILURE);
    }
    populate_source_cube(source, source_coordinates_file, n);

    // Calculate the resulting field
    double *result = run_poisson (n, source, iterations, threads, delta, debug);

    // Print out the middle slice of the cube for validation
    for (int y = 0; y < n; ++y)
    {
        for (int x = 0; x < n; ++x)
        {
            printf ("%0.5f ", result[((n / 2) * n + y) * n + x]);
        }
        printf ("\n");
    }

    free (source);
    free (result);

    return EXIT_SUCCESS;
}

static void populate_source_cube(double* source, FILE* file, int n)
{
    int x, y, z;
    double amplitude;

    // If no file was provided, set default source point at center of cube.
    if(file == NULL)
    {
        source[(DEFAULT_Z(n) * n + DEFAULT_Y(n)) * n + DEFAULT_X(n)] = DEFAULT_AMPLITUDE;
        return;
    }

    // Otherwise, populate from list of coordinates in file.
    for(int ii = 1; !feof(file); ii++) // scan till end of file
    {
        int num_matched = fscanf(file, "%d,%d,%d,%lf\n", &x,&y,&z,&amplitude);
        if (num_matched != 4) // Check elements are found in this line
        {
            fprintf (stderr, "Incorrect format for source coordinates (item %d) - please check README.\n", ii);
            exit (EXIT_FAILURE);
        }
        if (x < 0 || x > n-1 || y < 0 || y > n-1 || z < 0 || z > n-1) 
        {
            fprintf (stderr, "One or more coordinates out of range (item %d).\n", ii);
            exit (EXIT_FAILURE);
        }
        source[(z * n + y) * n + x] = amplitude;
    }
}