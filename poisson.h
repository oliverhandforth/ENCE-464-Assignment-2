/**
 * poisson.h
 * Implementation of a Poisson solver.
 *
*/

/**
 * @brief Solve Poissons equation for a given cube.
 *
 * @param n             The edge length of the cube. n^3 number of elements.
 * @param source        Pointer to the source term cube, a.k.a. forcing function.
 * @param iterations    Number of iterations to perform.
 * @param threads       Number of threads to use for solving.
 * @param delta         Grid spacing.
 * @param debug         Set to true when operating in debug mode to enable verbose logging.
 * @return double*      Solution to Poissons equation.  Caller must free.
 */
double* run_poisson (int n, double *source, int iterations, int threads, float delta, bool debug);