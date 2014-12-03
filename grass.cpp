#include <chrono>
#include <cstdio>
#include <cstring>
#include <ctime>
#include <vector>

#include <unistd.h>
using namespace std;

//#define DEBUG

const int maxlen = 1000;
char in_filename[maxlen];

// Define standard values for parameters
bool approx_reconstr_err = false; // Use approximation when computing the reconstruction error
bool do_queries = true;		    // If true, run queries.
int param_error_type = 1; // The type of reconstruction error to use. '1' for L1, '2' for L2, '3' for cut-norm.
int param_k = 18; // Number of supernodes at the end
float param_c = 1.0; // Controls the sample size. At each time t, we sample c*n(t) candidate pairs for merging, where n(t) is the number of supernodes at time t

#include "grass.h"

/*
 * Print usage on stderr
 *
 */
void usage(char *binary_name) {
    fprintf(stderr, "Usage: %s [-a] [-c] [-h] [-k INT] [-t {1|2|3] [EDGE_FILE]\n", binary_name);
    fprintf(stderr, " -a : use approximation when computing the reconstruction error.\n");
    fprintf(stderr, " -c FLOAT : controls the sample size. Must be in (0.0,1.0].\n");
    fprintf(stderr, " -h : print usage and exit.\n");
    fprintf(stderr, " -k INT: number of supernodes of the final summary. Must be positive.\n");
    fprintf(stderr, " -q: don't run queries (optional).\n");
    fprintf(stderr, " -t {1|2|3}: type of reconstruction error to use. L1, L2, or cut norm\n.");
}

/*
 * parse command line arguments
 *
 * -a : use approximation when computing the reconstruction error.
 * -c FLOAT : controls the sample size. Must be in (0.0,1.0]
 * -h : print usage and exit.
 * -k INT: number of supernodes of the final summary. Must be positive.
 * -t {1|2|3}: type of reconstruction error to use. L1, L2, or cut norm.
 *
 */
void parse_cmd_args(int argc, char* argv[]) {
    int opt;
    while ((opt = getopt(argc, argv, "ac:hi:k:t:")) != -1) {
        switch (opt) {
	    case 'a':
            approx_reconstr_err = true;
            break;
	    case 'c':
            param_c = atof(optarg);
            if (param_c <= 0.0 || param_c > 1.0) {
                fprintf(stderr, "Sample size parameter must be in (0.0,1.0]\n");
                usage(argv[0]);
                exit(1);
            }
            break;
	    case 'h':
	    usage(argv[0]);
            exit(0);
            break;
	    case 'k':
            param_k = atoi(optarg);
            if (param_k <= 0) {
                fprintf(stderr, "Number of supernodes must be positive\n");
                usage(argv[0]);
                exit(1);
            }
            break;
	    case 'q':
	    do_queries = false;
	    break;
	    case 't':
            param_error_type = atoi(optarg);
            if (param_error_type <= 0 || param_error_type > 3) {
                fprintf(stderr, "Error type must be one of 1,2,3\n");
                usage(argv[0]);
                exit(1);
            }
            break;
        default:
            fprintf(stderr, "Wrong option.\n");
            usage(argv[0]);
            exit(1);
        }
    }

    if (optind < argc) {
        strcpy(in_filename, argv[optind]);
        printf("Input file = %s\n", in_filename);
    } else {
        printf("Input file = stdin\n");
    }
    printf("param_k = %i, param_c = %f, param_error_type = %i, approx_reconstr_err = %i\n", param_k, param_c, param_error_type, approx_reconstr_err);
}

/*
 * Main function
 */
int main(int argc, char* argv[]) {
    setbuf(stdout, NULL);

    // Parse command line arguments
    parse_cmd_args(argc, argv);
    
    // Initialize random generator
    srand(time(NULL));

    // Read input graph 
    read_graph();

    // Save start time. The steady clock is guaranteed to be monotonic
    auto start = std::chrono::steady_clock::now();

    // Run the actual algorithm to build the summary
    grass_sample_pairs(param_k, param_c, param_error_type, approx_reconstr_err);

    // Save end time
    auto end = std::chrono::steady_clock::now();

    // Compute elapsed time
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    printf("Elapsed time: %lf\n", elapsed.count() / 1000.0);

    // Analyze the summary and collect the statistics (also print some)
    vector<double> stats = analyze_summary(approx_reconstr_err, do_queries);

    // Print on stderr info and stats about the graph and the summary in a CSV format
    fprintf(stderr, "grass, %s, %d, %d, %lf, %lf, %lld, %d, %d", in_filename, n, num_edges,
            avg_deg, avg_dens, triangles, param_k, approx_reconstr_err);
    for (double stat : stats) {
	    fprintf(stderr, ", %lf", stat);
    }
    fprintf(stderr, ", %d, %f", param_error_type, param_c);
    fprintf(stderr, ", %f", elapsed.count() / 1000.0);
    fprintf(stderr, "\n");

    printf("\nDone\n");

    return 0;
}

