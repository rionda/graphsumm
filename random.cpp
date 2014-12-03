#include <algorithm>
#include <cassert>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <map>
#include <numeric>
#include <set>
#include <vector>
#include <unordered_map>

#include <unistd.h>
using namespace std;

//#define DEBUG

const int maxlen = 1000;
char in_filename[maxlen];

// Define standard values for parameters
bool approx_reconstr_err = false; // If true, use approximation when computing the reconstruction error
bool do_queries = true;		    // If true, run queries.
int param_k = 18; // Number of supernodes at the end

#include "grass.h"
#include "graph.h"

/*
 * Print usage on stderr
 *
 */
void usage(char *binary_name) {
    fprintf(stderr, "Usage: %s [-a] [-h] [-k INT] [EDGE_FILE]\n", binary_name);
    fprintf(stderr, " -a : use approximation when computing the reconstruction error.\n");
    fprintf(stderr, " -h : print usage and exit.\n");
    fprintf(stderr, " -k INT: number of supernodes of the final summary. Must be positive.\n");
    fprintf(stderr, " -q: don't run queries (optional).\n");
}

/*
 * Parse command line arguments
 *
 * -a : use approximation when computing the reconstruction error.
 * -h : print usage and exit.
 * -k INT: number of supernodes of the final summary.
 *
 */
void parse_cmd_args(int argc, char* argv[]) {
    int opt;
    extern char *optarg;
    extern int optind;
    while ((opt = getopt(argc, argv, "ahk:")) != -1) {
        switch (opt) {
	    case 'a':
		approx_reconstr_err = true;
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

    printf("param_k = %i, approx_reconstr_err= %i\n", param_k, approx_reconstr_err);
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
    random_summary(param_k);

    // Save end time
    auto end = std::chrono::steady_clock::now();

    // Compute elapsed time
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    printf("Elapsed time: %lf\n", elapsed.count() / 1000.0);

    // Analyze the summary and collect the statistics (also print some)
    vector<double> stats = analyze_summary(approx_reconstr_err, do_queries);

    // Print on stderr info and stats about the graph and the summary in a CSV format
    fprintf(stderr, "random, %s, %d, %d, %lf, %lf, %lld, %d, %d", in_filename, n, num_edges,
            avg_deg, avg_dens, triangles, param_k, approx_reconstr_err);
    for (double stat : stats) {
	    fprintf(stderr, ", %lf", stat);
    }
    fprintf(stderr, ", %lf", elapsed.count() / 1000.0);
    fprintf(stderr, "\n");

    printf("\nDone\n");

    return 0;
}

