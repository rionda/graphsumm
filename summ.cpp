#include <cstdio>
#include <cstring>
#include <ctime>
#include <vector>

#include <unistd.h>
using namespace std;

const int maxlen = 1000;
char in_filename[maxlen];

bool use_random_summary = false;
bool approx_alg = false;            // If true, use approximation algorithm for clustering, not greedy
bool approx_reconstr_err = false;   // If true, use approximation when computing the reconstruction error
bool mini_batch = false;            // If true, use mini-batch algorithm for clustering
bool do_queries = true;		    // If true, run queries.
int param_error_type = 2;           // The type of reconstruction error to use. '1' for L1, '2' for L2, '3' for cut-norm
int param_k = 100;                  // Number of supernodes at the end
int param_dim = 800;                // Used for dimensionality reduction
int max_att_impr = 20000;           // max. attempted improvements
int param_iters = 15;               // number of iterations for k-means

#include "graph.h"
#include "szem_clust.h"

/*
 * Print usage on stderr
 *
 */
void usage(char *binary_name) {
    fprintf(stderr, "Usage: %s [-a] [-b] [d INT] [-h] [-k INT] [-p INT] [-q] [-r] [-s] [-t {1,2,3}] [EDGE_FILE]\n", binary_name);
    fprintf(stderr, " -a : use approximation when computing the reconstruction error.\n");
    fprintf(stderr, " -b : use mini-batch clustering algorithm (faster but potentially worse).\n");
    fprintf(stderr, " -d INT : number of dimensions to reduce the points to using Johnson-Lindenstrauss transform. If negative, don't use dimensionality reduction.\n");
    fprintf(stderr, " -h : print usage and exit.\n");
    fprintf(stderr, " -k INT: number of supernodes of the final summary. Must be positive.\n");
    fprintf(stderr, " -p[INT]: use approximation algorithm with the specified number of attempted improvements (optional).\n");
    fprintf(stderr, " -q: don't run queries (optional).\n");
    fprintf(stderr, " -r : build a random summary. Useful for testing and baseline.\n");
//    fprintf(stderr, " -s : use approximation when squaring the adjacency matrix.\n");
    fprintf(stderr, " -t {1|2|3}: type of reconstruction error to use. L1, L2, or cut norm\n");
    fprintf(stderr, " -m INT: maximum number of k-means iterations\n");
}

/*
 * parse command line arguments
 *
 * -a : use approximation when computing the reconstruction error.
 * -b : use the mini-batch clustering algorithm, not standard greedy.
 * -d INT : number of dimensions to reduce the points to using
 *  Johnson-Lindenstrauss transform. Must be positive.
 * -h : print usage and exit.
 * -k INT: number of supernodes of the final summary. Must be positive.
 * -p : use approximation algorithm, not standard 'greedy' clustering.
 * -r : build a random summary. Useful for testing and baseline.
 * -t {1|2|3}: type of reconstruction error to use. L1, L2, or cut norm.
 *
 */
void parse_cmd_args(int argc, char* argv[]) {
    int opt;
    while ((opt = getopt(argc, argv, "abd:hi:k:p::qrt:m:")) != -1) {
        switch (opt) {
        case 'a':
            approx_reconstr_err = true;
            break;
        case 'b':
            mini_batch = true;
            break;
        case 'd':
            param_dim = atoi(optarg);
            if (param_dim <= 0 && param_dim != -1) {
                fprintf(stderr, "Dimensions must be positive or '-1'\n");
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
        case 'm':
            param_iters = atoi(optarg);
            if (param_k < 0) {
                fprintf(stderr, "Number of k-means iterations must be nonnegative\n");
                usage(argv[0]);
                exit(1);
            }
            break;
        case 'p':
            approx_alg = true;
            if (optarg) {
                max_att_impr = atoi(optarg);
                if (max_att_impr <= 0) {
                    fprintf(stderr, "Number of attempted improvements must be positive\n");
                    usage(argv[0]);
                    exit(1);
                }
            }
            break;
	case 'q':
	    do_queries = false;
	    break;
        case 'r':
            use_random_summary = true;
            break;
        case 't':
            param_error_type = atoi(optarg);
            if (param_error_type <= 0 || param_error_type > 3) {
                fprintf(stderr, "Error type must be one of 1,2,3.\n");
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

    if (mini_batch && approx_alg) {
        fprintf(stderr, "Only one between -b and -p can be specified. Abort.\n");
        exit(1);
    }

    if (optind < argc) {
        strcpy(in_filename, argv[optind]);
        printf("Input file = %s\n", in_filename);
    } else {
        printf("Input file = stdin\n");
    }

    printf("param_k = %i param_dim = %i, param_error_type = %i, mini_batch = %i, approx_alg = %i, max_att_impr = %i, use_random_summary = %i, approx_reconstr_err = %i, dont_square = %i\n",
            param_k, param_dim, param_error_type, mini_batch, approx_alg, max_att_impr, use_random_summary, approx_reconstr_err, 0);
}

int main(int argc, char* argv[]) {
    setbuf(stdout, NULL);           // no buffering for stdout

    // Parse command line arguments
    parse_cmd_args(argc, argv);

    // Initialize random generator
#ifndef DEBUG
    srand(time(NULL));
#endif

    // Read input graph
    read_graph();
    printf("n = %i m = %i avg_deg = %lf avg_dens = %lf\n\n", n, num_edges, avg_deg, avg_dens);

    // Run the actual algorithm to build the summary
    vector<double> elapsed_times =
        summarize_by_clustering(param_k, param_error_type, mini_batch, approx_alg, use_random_summary, param_dim, param_iters);

    // Analyze the summary and collect the statistics (also print some)
    vector<double> stats = analyze_summary(approx_reconstr_err, do_queries);

    fprintf(stderr, "summ, %s, %d, %d, %lf, %lf, %lld, %d, %d", in_filename, n, num_edges,
            avg_deg, avg_dens, triangles, param_k, approx_reconstr_err);
    for (double stat : stats) {
        fprintf(stderr, ", %lf", stat);
    }
    fprintf(stderr, ", %d, %d, %d, %d, %d, %d", param_error_type, param_dim, 0, use_random_summary, mini_batch, approx_alg);
    fprintf(stderr, ", %lf, %lf, %lf", elapsed_times[0] + elapsed_times[1], elapsed_times[0], elapsed_times[1]);
    fprintf(stderr, "\n");

    printf("\nDone\n");

    return 0;
}
