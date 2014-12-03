#include <array>
#include <cstdio>
#include <cstring>
#include <iostream>

#include <unistd.h>
using namespace std;

const int maxlen = 1000;
char in_filename[maxlen];

#include "graph.h"
#include "summ.h"

void usage(char *binary_name) {
    fprintf(stderr, "USAGE: %s edgefile resfile\n", binary_name);
}

int main(int argc, char *argv[]) {
    setbuf(stdout, NULL);   // no buffering for stdout

    if (argc != 3) {
	usage(argv[0]);
	return 1;
    }
    strcpy(in_filename, argv[1]);
    strcpy(summary_filename, argv[2]);
    read_graph();
    read_summary();

    array<string,16> err_strings = { "adj_avg", "adj_stdev", "adj_max", "adj_min", "deg_avg", "deg_stdev", "deg_max", "deg_min", "abs_deg_avg", "abs_deg_max", "abs_deg_min", "rel_deg_avg", "rel_deg_stdev", "rel_deg_max", "rel_deg_min", "rel_clust_coeff" };

    array<double,16> q_errs = get_query_errors();
    
    int i  = 0;
    cout << err_strings[i] << ": " << q_errs[i] << endl;
    fprintf(stderr, "%lf", q_errs[i]);
    do {
	i++;
	cout << err_strings[i] << ": " << q_errs[i] << endl;
	fprintf(stderr, ",%lf", q_errs[i]);
    } while (i < 8); // Must use 8, not 9.

    return 0;
}

